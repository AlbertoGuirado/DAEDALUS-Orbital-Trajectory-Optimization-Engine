---
tags:
  - validacion
  - metodologia
aliases:
  - Tolerancias de test
  - RNF-15
---

# Tolerancias justificadas

> [!danger] Regla (RNF-15)
> **Ajustar una tolerancia hasta que un test pase invalida el test.** Convierte una comprobación de corrección en una descripción del comportamiento actual, incluido su posible error.

## Cómo se justifica una tolerancia
**Antes** de ejecutar, a partir de:
- el orden del método,
- el número de pasos,
- la precisión de la aritmética (épsilon, ver [[Umbral de degeneración numérica]]),
- la escala física del problema.

## Ejemplos en DAEDALUS
| Tolerancia | Justificación |
|---|---|
| < 1 m y < 0.001 m/s en un periodo | Resolución submétrica para escala terrestre; mayor que `atol` para permitir acumulación |
| Orden RK4 ∈ [3.7, 4.5] | Entorno asintótico de 4 admitiendo términos de orden superior |
| `rtol=1e-12` en ida y vuelta de elementos | Margen sobre el redondeo de float64 |
| Invariantes < 10⁻⁷ en 10 periodos | Coherente con `rtol=1e-10` acumulado sobre miles de pasos |

## Si un test falla
Se explica y se justifica cualquier cambio de tolerancia. **Nunca** se relaja, se elimina la aserción, se salta el test o se captura la excepción solo para que pase.

Relacionado: [[Tolerancias rtol y atol]], [[Batería de validación de la Fase 1]], [[Análisis de error numérico]].
