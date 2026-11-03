# Fase 3: API del motor de optimización

Documento parcial, en construcción, de la rama `feature/phase-3-optimization`.
Por ahora solo recoge la sección sobre capas e importación. El resto de la
API se documentará al cerrar el incremento.

## Capas y dependencias en tiempo de importación (Layering and import-time dependencies)

**Implemented:** el subpaquete `orbital_optimizer.optimization` (`problem.py`,
`result.py`, `optimizer.py`, `differential_evolution.py`, `benchmarks.py`)
importa directamente solo NumPy, la biblioteca estándar y sus propios módulos.
Lo comprueba `tests/test_layering.py`, que analiza el código con `ast` sin
ejecutarlo. Ese mismo test verifica también que `dynamics/`, `orbital/`,
`simulation/` y `src/physics/` no importan `optimization` ni `problems`.

**Alcance de la garantía:** la independencia es **estática y de imports
directos**. No es un aislamiento en tiempo de importación. Python ejecuta
primero el paquete padre: `import orbital_optimizer.optimization` ejecuta
`orbital_optimizer/__init__.py`, que importa la API de las fases 1–2.

Comprobación ejecutada el 2026-11-03 en el entorno del proyecto (Python 3.12.1,
NumPy 2.3.3, SciPy 1.16.3), tras importar solo `orbital_optimizer.optimization`:

- se cargaron `orbital_optimizer.dynamics.two_body`, `orbital.elements`,
  `orbital.kepler`, `orbital.mechanics`, `orbital.transfers`,
  `simulation.impulsive`, `simulation.kepler_propagator` y
  `simulation.propagator`;
- se cargaron `scipy` y `scipy.integrate`, este último a través de `simulation.propagator`;
- no se cargó `matplotlib`.

**Consecuencias:**

- SciPy debe estar instalado, y su coste de importación se paga aunque solo
  se use la optimización genérica.
- No hay dependencia circular, porque ninguna capa física importa este subpaquete.
- La API de optimización no se reexporta desde `orbital_optimizer`; se importa
  desde `orbital_optimizer.optimization`.

**Planned (no implementado):** carga perezosa en el paquete raíz, por ejemplo
con un `__getattr__` de módulo (PEP 562), para que importar el optimizador no
cargue la física ni SciPy. Exigiría cambiar `orbital_optimizer/__init__.py`,
que forma parte de la API pública de las fases 1–2, y por tanto requiere
aprobación previa.
