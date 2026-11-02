---
tags:
  - algoritmia/numerico
  - algoritmia/integracion
  - decision-de-diseno
estado: implementado
fase: 1
---

# Comparativa RK4 frente a RK45

| Criterio | [[Runge-Kutta 4]] | [[Dormand-Prince RK45]] |
|---|---|---|
| Control de error | No | Sí, **local** |
| Coste por paso | 4 evaluaciones | 6 efectivas (FSAL) |
| Malla temporal | Fija, predecible | Adaptativa |
| Órbitas casi circulares | Adecuado | Adecuado |
| Órbitas muy excéntricas | Ineficiente o impreciso | Claramente preferible (concentra pasos en el periapsis) |
| Estudios de convergencia | **Ideal** (paso prescrito) | No aplicable |
| Reproducibilidad exacta | Total | Alta, depende de la implementación |
| Muestreo | Altera la malla | Interpolación densa |
| Simpléctico | No | No |

> [!success] Decisión del proyecto
> **RK45 por defecto** para simulación de misión. **RK4 reservado** a estudios de orden de convergencia.

## Validación cruzada
RV-1.10: RK4 con paso fino y RK45 con tolerancias estrictas deben coincidir (< 1 m y < 0.001 m/s). Son implementaciones independientes (una a mano, otra de SciPy), así que su acuerdo es evidencia de que ninguna tiene un error grosero. Test: `test_agreement_between_integrators`.

La ampliación local de Fase 1 compara directamente 201 muestras comunes de
una órbita circular y una elipse inclinada. También contrasta ambos métodos
con `propagate_kepler`; ver [[Cierre documental de la Fase 1]].

Relacionado: [[Integradores simplécticos]], [[Orden de convergencia]].
