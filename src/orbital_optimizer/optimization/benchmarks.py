"""Standard continuous test functions with known optima (RV-3.02).

Dimensionless; independent of astrodynamics. Default bounds are the usual
literature domains. Because a symmetric box centers the optimum of Sphere and
Rastrigin, callers can pass asymmetric bounds to detect center bias.

    Sphere:      f(x) = sum_j x_j^2,                               x* = 0, f* = 0
    Rosenbrock:  f(x) = sum_j [100 (x_{j+1} - x_j^2)^2 + (1 - x_j)^2],  x* = 1, f* = 0, D >= 2
    Rastrigin:   f(x) = 10 D + sum_j [x_j^2 - 10 cos(2 pi x_j)],   x* = 0, f* = 0
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .problem import Bounds, FunctionProblem, _real_array


def _vector(x: ArrayLike, minimum_size: int = 1) -> NDArray[np.float64]:
    values = _real_array(x, "x")
    if values.ndim != 1 or values.size < minimum_size:
        raise ValueError(f"x must be a one-dimensional array with at least {minimum_size} entries")
    return values


def sphere(x: ArrayLike) -> float:
    """f(x) = sum_j x_j^2; global minimum 0 at x = 0."""
    values = _vector(x)
    return float(np.dot(values, values))


def rosenbrock(x: ArrayLike) -> float:
    """f(x) = sum_j [100 (x_{j+1} - x_j^2)^2 + (1 - x_j)^2], D >= 2; minimum 0 at x = 1."""
    values = _vector(x, 2)
    return float(np.sum(100.0 * (values[1:] - values[:-1] ** 2) ** 2 + (1.0 - values[:-1]) ** 2))


def rastrigin(x: ArrayLike) -> float:
    """f(x) = 10 D + sum_j [x_j^2 - 10 cos(2 pi x_j)]; minimum 0 at x = 0, many local minima."""
    values = _vector(x)
    return float(10.0 * values.size + np.sum(values**2 - 10.0 * np.cos(2.0 * np.pi * values)))


@dataclass(frozen=True, eq=False)
class BenchmarkFunction:
    """Test function with its bounds and known global minimizer and minimum."""

    name: str
    function: Callable[[NDArray[np.float64]], float]
    bounds: Bounds
    minimizer: NDArray[np.float64]
    minimum: float

    def __post_init__(self) -> None:
        minimizer = _vector(self.minimizer)
        if not self.bounds.contains(minimizer):
            raise ValueError("the known minimizer must lie within the bounds")
        minimizer.setflags(write=False)
        object.__setattr__(self, "minimizer", minimizer)
        object.__setattr__(self, "minimum", float(self.minimum))

    def problem(self) -> FunctionProblem:
        """FunctionProblem over these bounds; non-finite values are NON_FINITE."""
        return FunctionProblem(self.function, self.bounds, name=self.name)


def _bounds(dimension: int, lower: float, upper: float) -> Bounds:
    if isinstance(dimension, bool) or not isinstance(dimension, (int, np.integer)) or dimension < 1:
        raise ValueError("dimension must be an integer >= 1")
    return Bounds(np.full(dimension, lower, dtype=float), np.full(dimension, upper, dtype=float))


def sphere_benchmark(dimension: int, lower: float = -5.12, upper: float = 5.12) -> BenchmarkFunction:
    """Sphere on [lower, upper]^D (default [-5.12, 5.12]^D)."""
    return BenchmarkFunction(f"sphere-{dimension}d", sphere, _bounds(dimension, lower, upper),
                             np.zeros(dimension), 0.0)


def rosenbrock_benchmark(dimension: int, lower: float = -5.0, upper: float = 10.0) -> BenchmarkFunction:
    """Rosenbrock on [lower, upper]^D (default [-5, 10]^D), D >= 2."""
    if isinstance(dimension, (int, np.integer)) and not isinstance(dimension, bool) and dimension < 2:
        raise ValueError("Rosenbrock requires dimension >= 2")
    return BenchmarkFunction(f"rosenbrock-{dimension}d", rosenbrock, _bounds(dimension, lower, upper),
                             np.ones(dimension), 0.0)


def rastrigin_benchmark(dimension: int, lower: float = -5.12, upper: float = 5.12) -> BenchmarkFunction:
    """Rastrigin on [lower, upper]^D (default [-5.12, 5.12]^D)."""
    return BenchmarkFunction(f"rastrigin-{dimension}d", rastrigin, _bounds(dimension, lower, upper),
                             np.zeros(dimension), 0.0)
