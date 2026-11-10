"""RNF-05: one-way layering between physics/simulation and the generic optimizer.

Static AST check of direct imports (it does not execute the code). Relative
imports are resolved against each file's package, ``from pkg import name`` is
expanded to both ``pkg`` and ``pkg.name`` (``name`` may be a submodule), and
calls to ``importlib.import_module`` / ``__import__`` are flagged as
violations because their target cannot be verified statically.

Files under the physics rule (dynamics/, orbital/, simulation/, visualization/,
src/physics/ and the root orbital_optimizer/__init__.py) must not import
optimization/ or problems/, and must not use the dynamic-import machinery at
all: any import of importlib, runpy or pkgutil and any exec/eval call is a
violation, which also closes aliasing tricks such as
``getattr(importlib, "import_module")``. The problems/ adapter layer is
deliberately outside both rules: it may import physics and optimization.
A static check cannot see every runtime trick (e.g. code loaded from data
files); it guards against ordinary and accidental violations. Import-time
behaviour of the root package is a separate matter: see docs/phase3-api.md.
"""
import ast
import sys
import unittest
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1] / "src"
PACKAGE = SOURCE / "orbital_optimizer"
PHYSICS_DIRS = [PACKAGE / "dynamics", PACKAGE / "orbital", PACKAGE / "simulation",
                PACKAGE / "visualization", PACKAGE / "constraints", SOURCE / "physics"]
ROOT_INIT = PACKAGE / "__init__.py"
# Modules of the physics layer that live directly in the package root.
PHYSICS_MODULES = [PACKAGE / "constants.py", PACKAGE / "spacecraft.py"]
OPTIMIZATION_DIR = PACKAGE / "optimization"
OPTIMIZATION = "orbital_optimizer.optimization"
# Layers above the physics, lowest first. A layer may import the ones before it, never the ones after.
UPPER_LAYERS = ("orbital_optimizer.optimization", "orbital_optimizer.problems", "orbital_optimizer.benchmarking",
                "orbital_optimizer.mission")
FORBIDDEN_FOR_PHYSICS = UPPER_LAYERS
# Interface libraries of the Verification Console: the core must not depend on them.
UI_LIBRARIES = ("streamlit", "plotly", "pandas", "altair")
DYNAMIC_IMPORT_MODULES = ("importlib", "runpy", "pkgutil")
ALLOWED_THIRD_PARTY = {"numpy"}
DYNAMIC = "<dynamic import>"
EXEC = "<exec/eval call>"


def package_of(path: Path) -> str:
    """Dotted package containing a source file under src/."""
    parts = list(path.relative_to(SOURCE).with_suffix("").parts)
    return ".".join(parts[:-1])


def _resolve(package: str, level: int, module: str | None) -> str:
    if not level:
        return module or ""
    base = package.split(".")
    if level - 1 >= len(base):
        raise ValueError(f"relative import beyond top-level package from {package!r}")
    base = base[: len(base) - (level - 1)]
    return ".".join(base + ([module] if module else []))


def imported_names(source: str, package: str) -> list[str]:
    """Absolute dotted names a file may import; DYNAMIC for import_module/__import__ calls."""
    names: list[str] = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            names.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            target = _resolve(package, node.level, node.module)
            names.append(target)
            names.extend(f"{target}.{alias.name}" for alias in node.names if alias.name != "*")
        elif isinstance(node, ast.Call):
            function = node.func
            called = function.id if isinstance(function, ast.Name) else (
                function.attr if isinstance(function, ast.Attribute) else None)
            if called in ("import_module", "__import__"):
                names.append(DYNAMIC)
        # Any reference, not only a call, so that aliases such as ``f = eval`` are caught too.
        referenced = node.id if isinstance(node, ast.Name) else (
            node.attr if isinstance(node, ast.Attribute) else None)
        if referenced in ("exec", "eval"):
            names.append(EXEC)
        elif referenced == "__import__":
            names.append(DYNAMIC)
    return names


def _within(name: str, prefix: str) -> bool:
    return name == prefix or name.startswith(prefix + ".")


def physics_violations(source: str, package: str) -> list[str]:
    """Imports of optimization/problems, dynamic-import machinery or exec/eval in a physics-rule file."""
    return [name for name in imported_names(source, package)
            if name in (DYNAMIC, EXEC)
            or any(_within(name, prefix) for prefix in FORBIDDEN_FOR_PHYSICS + DYNAMIC_IMPORT_MODULES)]


def optimization_violations(source: str, package: str) -> list[str]:
    """Imports other than stdlib, NumPy, __future__ and the optimization subpackage itself."""
    violations = []
    for name in imported_names(source, package):
        root = name.split(".")[0]
        allowed = name != DYNAMIC and (
            _within(name, OPTIMIZATION) or root == "__future__"
            or root in sys.stdlib_module_names or root in ALLOWED_THIRD_PARTY)
        if not allowed:
            violations.append(name)
    return violations


def source_files(directory: Path) -> list[Path]:
    return sorted(path for path in directory.rglob("*.py") if "__pycache__" not in path.parts)


class CheckerSelfTests(unittest.TestCase):
    """The checker must catch known-bad snippets and accept known-good ones."""

    def test_physics_checker_flags_known_bad_imports(self) -> None:
        cases = [
            ("from orbital_optimizer import optimization", "orbital_optimizer.orbital"),
            ("from orbital_optimizer import optimization as o", "orbital_optimizer.orbital"),
            ("from orbital_optimizer import problems", "orbital_optimizer.orbital"),
            ("import orbital_optimizer.optimization.problem as p", "orbital_optimizer.orbital"),
            ("from orbital_optimizer.optimization.problem import Bounds", "orbital_optimizer.dynamics"),
            ("from .. import optimization", "orbital_optimizer.simulation"),
            ("from ..optimization.problem import Bounds", "orbital_optimizer.simulation"),
            ("from ..problems import x", "orbital_optimizer.orbital"),
            ("from ... import optimization", "orbital_optimizer.simulation.nested"),
            ("import importlib\nm = importlib.import_module('orbital_optimizer.optimization')", "orbital_optimizer.orbital"),
            ("m = __import__('orbital_optimizer')", "orbital_optimizer.orbital"),
            ("from importlib import import_module\nimport_module('x')", "physics"),
            # Bypass attempts that only the stricter physics rule catches.
            ("from importlib import import_module as load\nload('orbital_optimizer.optimization')", "orbital_optimizer.orbital"),
            ("import importlib as il\ngetattr(il, 'import_' + 'module')('orbital_optimizer.problems')", "orbital_optimizer.orbital"),
            ("import importlib.util", "orbital_optimizer.simulation"),
            ("exec('import orbital_optimizer.optimization')", "orbital_optimizer.dynamics"),
            ("f = eval\nf('1')", "orbital_optimizer.orbital"),
            ("load = __import__\nload('orbital_optimizer.optimization')", "orbital_optimizer.orbital"),
            ("import runpy\nrunpy.run_module('orbital_optimizer.optimization')", "orbital_optimizer.visualization"),
            ("import pkgutil", "orbital_optimizer"),
            ("import builtins\nbuiltins.exec('x = 1')", "orbital_optimizer.orbital"),
            ("from .optimization import DifferentialEvolution", "orbital_optimizer"),
            ("from . import problems", "orbital_optimizer"),
        ]
        for source, package in cases:
            with self.subTest(source=source, package=package):
                self.assertTrue(physics_violations(source, package))

    def test_physics_checker_accepts_known_good_imports(self) -> None:
        cases = [
            ("import numpy as np", "orbital_optimizer.orbital"),
            ("from ..dynamics.two_body import EARTH_MU", "orbital_optimizer.orbital"),
            ("from .mechanics import specific_energy", "orbital_optimizer.orbital"),
            ("from orbital_optimizer import EARTH_MU", "physics"),
            ("from ... import dynamics", "orbital_optimizer.simulation.nested"),
            ("from scipy.integrate import solve_ivp", "orbital_optimizer.simulation"),
            ("from .dynamics.two_body import EARTH_MU", "orbital_optimizer"),
            ("from ..simulation.propagator import PropagationResult", "orbital_optimizer.visualization"),
        ]
        for source, package in cases:
            with self.subTest(source=source, package=package):
                self.assertEqual(physics_violations(source, package), [])

    def test_optimization_checker_flags_known_bad_imports(self) -> None:
        cases = [
            ("import scipy", OPTIMIZATION),
            ("import matplotlib.pyplot as plt", OPTIMIZATION),
            ("import orbital_optimizer", OPTIMIZATION),
            ("from orbital_optimizer import dynamics", OPTIMIZATION),
            ("from orbital_optimizer import EARTH_MU", OPTIMIZATION),
            ("from .. import orbital", OPTIMIZATION),
            ("from ..simulation.propagator import propagate", OPTIMIZATION),
            ("from ... import simulation", OPTIMIZATION + ".nested"),
            ("import importlib\nimportlib.import_module('numpy')", OPTIMIZATION),
        ]
        for source, package in cases:
            with self.subTest(source=source, package=package):
                self.assertTrue(optimization_violations(source, package))

    def test_optimization_checker_accepts_known_good_imports(self) -> None:
        cases = [
            ("from __future__ import annotations", OPTIMIZATION),
            ("import math\nfrom collections.abc import Mapping", OPTIMIZATION),
            ("import numpy as np\nfrom numpy.typing import NDArray", OPTIMIZATION),
            ("from .problem import Bounds\nfrom . import result", OPTIMIZATION),
            ("from orbital_optimizer.optimization.problem import Bounds", OPTIMIZATION),
            ("from .. import problem", OPTIMIZATION + ".nested"),
        ]
        for source, package in cases:
            with self.subTest(source=source, package=package):
                self.assertEqual(optimization_violations(source, package), [])

    def test_relative_import_beyond_top_level_is_an_error(self) -> None:
        with self.assertRaises(ValueError):
            imported_names("from ... import x", "orbital_optimizer.orbital")


class LayeringTests(unittest.TestCase):
    def test_source_tree_is_present(self) -> None:
        self.assertTrue(source_files(OPTIMIZATION_DIR))
        for directory in PHYSICS_DIRS:
            self.assertTrue(source_files(directory), directory)

    def test_physics_and_simulation_do_not_import_optimization(self) -> None:
        files = [ROOT_INIT] + PHYSICS_MODULES + [path for directory in PHYSICS_DIRS for path in source_files(directory)]
        for path in files:
            with self.subTest(file=str(path.relative_to(SOURCE))):
                self.assertTrue(path.is_file())
                self.assertEqual(physics_violations(path.read_text(encoding="utf-8"), package_of(path)), [])

    def test_constants_module_is_a_leaf(self) -> None:
        # RD-01: one place for the constants, importable by every layer, importing none.
        source = (PACKAGE / "constants.py").read_text(encoding="utf-8")
        self.assertEqual(imported_names(source, "orbital_optimizer"), [])

    def test_upper_layers_do_not_import_the_layers_above_them(self) -> None:
        for index, layer in enumerate(UPPER_LAYERS):
            directory = PACKAGE / layer.rsplit(".", 1)[1]
            if not directory.is_dir():          # a layer of a later phase that does not exist yet
                continue
            above = UPPER_LAYERS[index + 1:]
            for path in source_files(directory):
                names = imported_names(path.read_text(encoding="utf-8"), package_of(path))
                with self.subTest(file=str(path.relative_to(SOURCE))):
                    self.assertEqual([name for name in names if any(_within(name, prefix) for prefix in above)], [])

    def test_upper_layer_checker_flags_known_bad_imports(self) -> None:
        above_problems = UPPER_LAYERS[2:]
        for source in ("from ..benchmarking import run_benchmark", "import orbital_optimizer.mission",
                       "from orbital_optimizer.benchmarking.runner import run_once"):
            names = imported_names(source, "orbital_optimizer.problems")
            with self.subTest(source=source):
                self.assertTrue([name for name in names if any(_within(name, prefix) for prefix in above_problems)])

    def test_core_does_not_import_interface_libraries(self) -> None:
        # The Verification Console depends on the core; the core never depends on its libraries.
        for path in source_files(SOURCE):
            names = imported_names(path.read_text(encoding="utf-8"), package_of(path))
            with self.subTest(file=str(path.relative_to(SOURCE))):
                self.assertEqual([name for name in names if name.split(".")[0] in UI_LIBRARIES], [])

    def test_root_package_does_not_reexport_optimization_or_problems(self) -> None:
        # Q5: the optimizer is imported from orbital_optimizer.optimization, not from the root.
        import orbital_optimizer

        exported = set(orbital_optimizer.__all__)
        for name in ("DifferentialEvolution", "DifferentialEvolutionSettings", "Bounds", "TerminationCriteria",
                     "OptimizationResult", "CoplanarTwoImpulseProblem", "TwoImpulsePlan",
                     "verify_two_impulse_plan", "optimization", "problems"):
            with self.subTest(name=name):
                self.assertNotIn(name, exported)
        self.assertFalse(hasattr(orbital_optimizer, "DifferentialEvolution"))
        self.assertFalse(hasattr(orbital_optimizer, "CoplanarTwoImpulseProblem"))

    def test_optimization_imports_only_stdlib_numpy_and_itself(self) -> None:
        for path in source_files(OPTIMIZATION_DIR):
            with self.subTest(file=str(path.relative_to(SOURCE))):
                self.assertEqual(optimization_violations(path.read_text(encoding="utf-8"), package_of(path)), [])


if __name__ == "__main__":
    unittest.main()
