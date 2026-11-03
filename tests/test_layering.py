"""RNF-05: one-way layering between physics/simulation and the generic optimizer.

Static AST check of direct imports (it does not execute the code). Relative
imports are resolved against each file's package, ``from pkg import name`` is
expanded to both ``pkg`` and ``pkg.name`` (``name`` may be a submodule), and
calls to ``importlib.import_module`` / ``__import__`` are flagged as
violations because their target cannot be verified statically. Import-time
behaviour of the root package is a separate matter: see docs/phase3-api.md.
"""
import ast
import sys
import unittest
from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1] / "src"
PACKAGE = SOURCE / "orbital_optimizer"
PHYSICS_DIRS = [PACKAGE / "dynamics", PACKAGE / "orbital", PACKAGE / "simulation", SOURCE / "physics"]
OPTIMIZATION_DIR = PACKAGE / "optimization"
OPTIMIZATION = "orbital_optimizer.optimization"
FORBIDDEN_FOR_PHYSICS = ("orbital_optimizer.optimization", "orbital_optimizer.problems")
ALLOWED_THIRD_PARTY = {"numpy"}
DYNAMIC = "<dynamic import>"


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
    return names


def _within(name: str, prefix: str) -> bool:
    return name == prefix or name.startswith(prefix + ".")


def physics_violations(source: str, package: str) -> list[str]:
    """Imports of optimization/problems (or dynamic imports) from a physics/simulation file."""
    return [name for name in imported_names(source, package)
            if name == DYNAMIC or any(_within(name, prefix) for prefix in FORBIDDEN_FOR_PHYSICS)]


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
        for directory in PHYSICS_DIRS:
            for path in source_files(directory):
                with self.subTest(file=str(path.relative_to(SOURCE))):
                    self.assertEqual(physics_violations(path.read_text(encoding="utf-8"), package_of(path)), [])

    def test_optimization_imports_only_stdlib_numpy_and_itself(self) -> None:
        for path in source_files(OPTIMIZATION_DIR):
            with self.subTest(file=str(path.relative_to(SOURCE))):
                self.assertEqual(optimization_violations(path.read_text(encoding="utf-8"), package_of(path)), [])


if __name__ == "__main__":
    unittest.main()
