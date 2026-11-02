---
tags:
  - moc
aliases:
  - Home
  - Inicio
---

# DAEDALUS — Vault de conocimiento

Base de conocimiento de **DAEDALUS — Orbital Trajectory Optimization Engine**: los conceptos físicos y algorítmicos en los que se apoya el proyecto, enlazados entre sí y con el código que los implementa.

> [!info] Cómo leer este vault
> Cada nota de concepto lleva en sus propiedades:
> - `estado` — `implementado` (existe en `src/`), `oraculo-de-tests` (solo se usa como referencia en tests), `planificado` (fase futura) o `teoria` (contexto, no se implementará).
> - `fase` — la fase del roadmap a la que pertenece. Ver [[Roadmap por fases]].
> - `codigo` — ruta del módulo cuando el concepto está implementado.

## Mapas de contenido

| Área | Qué contiene |
|---|---|
| [[MOC Física orbital]] | Dos cuerpos, invariantes, elementos orbitales, maniobras, propulsión, perturbaciones |
| [[MOC Métodos numéricos]] | Integración RK4/RK45, error numérico, resolución de ecuaciones |
| [[MOC Optimización]] | Formulación, metaheurísticas, métodos locales, control óptimo, multiobjetivo, restricciones |
| [[MOC Validación y benchmark]] | Invariantes como oráculo, niveles de evidencia, estadística de comparación |
| [[MOC Proyecto DAEDALUS]] | Arquitectura, roadmap, principios, resultados medidos |

## Ruta de lectura recomendada

1. [[Problema de dos cuerpos]] → [[Vector de estado]] → [[Runge-Kutta 4]] / [[Dormand-Prince RK45]]
2. [[Energía específica]] · [[Momento angular específico]] · [[Vector de excentricidad]] → [[Invariantes como oráculo de validación]]
3. [[Elementos orbitales clásicos]] → [[Casos degenerados de los elementos orbitales]]
4. [[Maniobra impulsiva]] → [[Transferencia de Hohmann]] (Implemented, Phase 2) → [[Ecuación de Tsiolkovsky]]
5. [[Formulación del problema de optimización]] → [[Differential Evolution]] → [[Metodología de benchmark]]

## Estado de los conceptos

Conceptos ya implementados en el código:

```query
[estado:implementado]
```

Conceptos planificados para fases futuras:

```query
[estado:planificado]
```

## Corte documental — 2026-11-02

**Implemented:** fases 1 y 2 obligatorias, con 47 tests locales OK.
**Planned:** optimización y fases posteriores; también bielíptica, cambio de plano y Lambert.
**Known Limitations:** [[Validación de la Fase 2 e integración]]. Los resultados históricos
conservan sus fechas y no describen por sí solos el checkout actual.
