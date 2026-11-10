"""Generic optimization engine (Phase 3): problems, results and optimizers.

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

Optimizers (all implement ``Optimizer``; identical results for the same seed,
settings, problem and environment): DifferentialEvolution (reference
algorithm), ParticleSwarm, GeneticAlgorithm, RandomSearch (null baseline),
NelderMead (local refinement) and HybridOptimizer (global + final local
refinement). None of them guarantees a global optimum.

Constraints (Phase 4): a problem may report a total normalized violation
with each evaluation. The optimizers rank candidates through a constraint
handling strategy (constraint_handling: FeasibilityRules by default,
DeathPenalty, StaticPenalty) and always report the best feasible point.
"""

from .benchmarks import (
    BenchmarkFunction, rastrigin, rastrigin_benchmark, rosenbrock, rosenbrock_benchmark,
    sphere, sphere_benchmark,
)
from .constraint_handling import (
    ConstraintHandling, DeathPenalty, FeasibilityRules, StaticPenalty, feasibility_key, handling_from_mapping,
)
from .differential_evolution import ALGORITHM, DifferentialEvolution, DifferentialEvolutionSettings
from .genetic_algorithm import GeneticAlgorithm, GeneticAlgorithmSettings
from .hybrid import HybridOptimizer, HybridSettings, merge_results
from .nelder_mead import NelderMead, NelderMeadSettings
from .objective import ObjectiveTerm, ScalarizedProblem, WeightedObjective
from .optimizer import (
    EvaluationBudget, GenerationLog, Optimizer, check_termination, environment_info, finish_run,
    make_rng, problem_info, validate_seed,
)
from .particle_swarm import ParticleSwarm, ParticleSwarmSettings
from .problem import Bounds, Evaluation, EvaluationStatus, FunctionProblem, Problem
from .random_search import RandomSearch, RandomSearchSettings
from .result import (
    STATUS_ORDER, TERMINATION_PRECEDENCE, OptimizationHistory, OptimizationResult,
    TerminationCriteria, TerminationReason,
)

__all__ = [
    "Bounds", "Evaluation", "EvaluationStatus", "FunctionProblem", "Problem",
    "TerminationCriteria", "TerminationReason", "TERMINATION_PRECEDENCE", "STATUS_ORDER",
    "OptimizationHistory", "OptimizationResult",
    "Optimizer", "EvaluationBudget", "GenerationLog", "check_termination", "finish_run",
    "environment_info", "make_rng", "problem_info", "validate_seed",
    "ALGORITHM", "DifferentialEvolution", "DifferentialEvolutionSettings",
    "ParticleSwarm", "ParticleSwarmSettings", "GeneticAlgorithm", "GeneticAlgorithmSettings",
    "RandomSearch", "RandomSearchSettings", "NelderMead", "NelderMeadSettings",
    "HybridOptimizer", "HybridSettings", "merge_results",
    "ObjectiveTerm", "WeightedObjective", "ScalarizedProblem",
    "ConstraintHandling", "FeasibilityRules", "DeathPenalty", "StaticPenalty", "feasibility_key",
    "handling_from_mapping",
    "BenchmarkFunction", "sphere", "rosenbrock", "rastrigin",
    "sphere_benchmark", "rosenbrock_benchmark", "rastrigin_benchmark",
]
