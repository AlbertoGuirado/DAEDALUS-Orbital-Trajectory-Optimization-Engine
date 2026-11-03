"""Generic optimization engine (Phase 3): problems, results and Differential Evolution.

Physics independence is static and limited to direct imports: the modules of
this subpackage import only NumPy, the standard library and each other
(checked by tests/test_layering.py). It is not import-time isolation. Python
executes the parent package first, so ``import orbital_optimizer.optimization``
runs ``orbital_optimizer/__init__.py``, which loads the dynamics, orbital and
simulation modules and, through the propagator, SciPy (``scipy.integrate``).
Consequences: SciPy must be installed and its import cost is paid even for
purely generic optimization; there is no circular dependency, because nothing
in the physics layers imports this subpackage. Lazy loading in the root
package (e.g. a module-level ``__getattr__``, PEP 562) would avoid that cost;
it is a possible future change, not implemented. See docs/phase3-api.md.
Orbital problems must adapt to these interfaces from a separate layer.
"""

from .benchmarks import (
    BenchmarkFunction, rastrigin, rastrigin_benchmark, rosenbrock, rosenbrock_benchmark,
    sphere, sphere_benchmark,
)
from .differential_evolution import ALGORITHM, DifferentialEvolution, DifferentialEvolutionSettings
from .optimizer import EvaluationBudget, Optimizer, environment_info, make_rng, problem_info, validate_seed
from .problem import Bounds, Evaluation, EvaluationStatus, FunctionProblem, Problem
from .result import (
    STATUS_ORDER, TERMINATION_PRECEDENCE, OptimizationHistory, OptimizationResult,
    TerminationCriteria, TerminationReason,
)

__all__ = [
    "Bounds", "Evaluation", "EvaluationStatus", "FunctionProblem", "Problem",
    "TerminationCriteria", "TerminationReason", "TERMINATION_PRECEDENCE", "STATUS_ORDER",
    "OptimizationHistory", "OptimizationResult",
    "Optimizer", "EvaluationBudget", "environment_info", "make_rng", "problem_info", "validate_seed",
    "ALGORITHM", "DifferentialEvolution", "DifferentialEvolutionSettings",
    "BenchmarkFunction", "sphere", "rosenbrock", "rastrigin",
    "sphere_benchmark", "rosenbrock_benchmark", "rastrigin_benchmark",
]
