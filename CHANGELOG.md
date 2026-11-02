# Changelog

## Unreleased — integración de Fase 1 (Kepler) y Fase 2 en `develop`

- Fase 1: ecuación de Kepler (`solve_kepler`) y propagación kepleriana
  analítica (`propagate_kepler`); documentación de cierre de la Fase 1.
- Fase 2: transferencia de Hohmann (`hohmann_transfer`, `HohmannTransfer`),
  propelente por Tsiolkovsky (`propellant_mass`, `STANDARD_GRAVITY`) y
  propagación con impulsos instantáneos (`propagate_with_impulses`).
- 47 tests (32 de Fase 1, 15 de Fase 2).
- Configuración exploratoria LEO similar a la ISS en `configs/experiments/`;
  la referencia validada `configs/earth_orbit.json` no cambia.
- Documentación: deuda técnica DT-10 a DT-16, Fase 6 dividida en 6A
  (benchmark, antiguos RF-6.0n → RF-6A.0n) y 6B (misión, borrador).
- Retirados de `develop` los recursos no científicos `visualization/` y
  `vault/.obsidian/` (conservados en `feature/phase-2-hohmann`).
- Interfaz Streamlit aislada en `feature/streamlit-ui`, fuera del núcleo.

## 0.1.0 — Fase 1

- Dinámica de dos cuerpos con estado cartesiano y unidades SI.
- Conversiones de elementos elípticos y tratamiento de casos degenerados.
- Propagación RK4 y RK45; energía, momento angular, excentricidad y Δv.
- 20 tests de física, geometría, propagación y experimentos.
- Ejemplo terrestre configurable con resultados JSON, CSV y figuras opcionales.
- Arquitectura `src`, documentación científica y compatibilidad `physics`.
- CI de tests y experimento reproducible en GitHub Actions.
