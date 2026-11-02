---
tags:
  - proyecto
aliases:
  - Fases
  - Roadmap
---

# Roadmap por fases

Cada fase debe estar **validada antes de empezar la siguiente**: un error en la dinámica contamina en silencio cualquier optimización construida encima (ver [[Principios de ingeniería]]).

| Fase | Nombre | Contenido | Estado (2026-11-02) |
|---|---|---|---|
| 1 | Orbital Dynamics Engine | [[Problema de dos cuerpos]], [[Elementos orbitales clásicos]], [[Runge-Kutta 4]], [[Dormand-Prince RK45]], invariantes, coste Δv | Implemented: baseline y Kepler integrados |
| 2 | Baseline analítico | [[Transferencia de Hohmann]], [[Ecuación de Tsiolkovsky]]; opcional [[Transferencia bielíptica]], [[Cambio de plano]], [[Problema de Lambert]] | Implemented: Hohmann, Tsiolkovsky e impulsos; alternativas opcionales Planned |
| 3 | Optimization Engine | [[Formulación del problema de optimización]], [[Differential Evolution]], interfaz común de optimizadores, multiobjetivo | Planned |
| 4 | Constraints | [[Manejo de restricciones]]: masa, empuje, Δv, tiempo, altitud mínima, propelente | Planned |
| 5 | Perturbations | [[Encendido de duración finita]], masa variable, [[Perturbación J2]], [[Arrastre atmosférico]] | Planned |
| 6A | Benchmark científico | Comparación reproducible de algoritmos | Planned |
| 6B | Mission Optimization Engine | Requisito de misión → estrategia de transferencia | Planned |

## Criterios de cierre

> [!abstract] Fase 1 (cumplido)
> Paquete importable, tests obligatorios en verde, experimento ejecutable que informa la deriva de invariantes, tolerancias justificadas antes de ejecutar. Ver [[Resultados de validación de la Fase 1]].

Reproducción y auditoría documental: [[Cierre documental de la Fase 1]]
(2026-10-16). Se conservan 20 tests correctos y las cinco comprobaciones del
experimento. RV-1.10 se amplió después a 201 muestras comunes RK4/RK45 para una
circular y una elipse inclinada, y pasó.

La ampliación Kepler integrada implementa solver de Kepler y
propagación analítica, todavía fuera de `v0.1.0`. Su ejecución final pasó 32 tests tras justificar el término de
acondicionamiento del estado inicial en el modelo de error de los tests
analíticos. El historial del fallo inicial de e=0.95 se conserva en el
registro técnico. El tag anterior conserva solo el baseline de 20 tests.

> [!abstract] Fase 2 (requisitos obligatorios cumplidos en los casos probados)
> El Δv de Hohmann se verifica **aplicando los impulsos al propagador** y comprobando que la órbita final coincide con la objetivo (RV-2.01). El wrapper impulsivo ya está implementado; `propagate` conserva vuelo libre. Ver [[Maniobra impulsiva]] y [[Validación de la Fase 2 e integración]].

> [!todo] Fase 3
> Differential Evolution reproduce el óptimo de Hohmann en el caso en que este es el óptimo conocido (RV-3.01).

> [!todo] Fase 5
> Con todas las perturbaciones desactivadas se recupera exactamente el modelo de dos cuerpos (RV-5.01).
