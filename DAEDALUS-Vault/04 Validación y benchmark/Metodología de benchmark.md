---
tags:
  - validacion
  - metodologia
  - algoritmia/optimizacion
aliases:
  - Benchmark
  - Comparación de algoritmos
estado: planificado
fase: 3
---

# Metodología de benchmark

> [!abstract] Principio rector
> Una única ejecución de cada algoritmo **no demuestra nada**: la varianza entre ejecuciones de un mismo algoritmo estocástico puede superar la diferencia entre algoritmos.

## Diseño experimental
| Elemento | Requisito |
|---|---|
| Semillas | Múltiples, independientes, **registradas**; número fijado a priori |
| Presupuesto | **Idéntico nº de evaluaciones de función** para todos (RF-6.03) |
| Escenarios | **Fijados antes** de ejecutar (RF-6.04) |
| Criterio de parada | Presupuesto agotado, no "convergencia" subjetiva |
| Hiperparámetros | Mismo esfuerzo de ajuste para todos, o valores por defecto de la literatura |

> [!important] Evaluaciones, no tiempo de pared
> El coste está dominado por la propagación, común a todos. El tiempo de pared mide además la calidad de la implementación y el estado de la máquina; las evaluaciones aíslan la **eficiencia algorítmica**. El tiempo se registra, pero nunca es criterio de parada.

## Métricas
| Métrica | Informa de |
|---|---|
| Mejor objetivo | Calidad de la mejor solución |
| **Mediana e IQR** | Comportamiento típico y dispersión (preferidas a media y desviación: distribuciones asimétricas) |
| Tasa de éxito | Fracción de ejecuciones factibles bajo umbral declarado |
| Evaluaciones hasta objetivo | Velocidad de convergencia |
| Curva de convergencia | Mejor valor vs evaluaciones, con banda intercuartílica |
| Violaciones | Magnitud y tipo |

## Gráficas (ejes siempre con unidades)
Trayectoria con impulsos marcados · curvas de convergencia · Δv frente a tiempo de cómputo · diagramas de caja por algoritmo · [[Dominancia de Pareto|frente de Pareto]] con Hohmann como referencia.

## Amenazas a la validez
| Amenaza | Mitigación |
|---|---|
| Ajuste desigual de hiperparámetros | Mismo esfuerzo para todos |
| Escenarios favorables | Fijados antes |
| Semillas insuficientes | Número fijado a priori, no ampliado hasta "obtener" el resultado |
| Generalización excesiva | Declarar la clase de problemas ([[Teorema No Free Lunch]]) |
| Rivales mal configurados | Implementaciones de referencia con configuración documentada |

Relacionado: [[Tests estadísticos no paramétricos]], [[Formulación del problema de optimización]].
