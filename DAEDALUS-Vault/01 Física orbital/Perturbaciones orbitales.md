---
tags:
  - fisica/perturbaciones
aliases:
  - Método de Cowell
  - Método de Encke
estado: planificado
fase: 5
---

# Perturbaciones orbitales

Aceleraciones que se añaden al término kepleriano del [[Problema de dos cuerpos]]:

$$\ddot{\mathbf{r}} = -\frac{\mu\mathbf{r}}{r^{3}} + \mathbf{a}_{J_2} + \mathbf{a}_{\text{arrastre}} + \mathbf{a}_{\text{empuje}} + \ldots$$

## Cowell frente a Encke
| Método | Idea | Veredicto |
|---|---|---|
| **Cowell** | Integrar la ecuación completa directamente | ✅ Adoptado: simple y compatible con el propagador actual |
| **Encke** | Integrar solo la desviación respecto a una cónica de referencia | Solo si se demuestra una limitación de precisión concreta |

## Régimen de relevancia
| Perturbación | Cuándo importa | ¿En DAEDALUS? |
|---|---|---|
| [[Perturbación J2]] | Dominante en LEO; cae como $r^{-4}$ | Fase 5 |
| [[Arrastre atmosférico]] | < ~600 km; despreciable > ~1000 km | Fase 5 |
| Tercer cuerpo (Luna, Sol) | Crece con la altitud; GEO y superiores | Fuera de alcance |
| Presión de radiación solar | Relación área/masa alta | Fuera de alcance |
| Armónicos superiores | Correcciones sobre J2 | Fuera de alcance |

## Requisitos clave
- RF-5.01: cada perturbación se activa/desactiva de forma **independiente**.
- RF-5.06: la configuración de perturbaciones activas se registra con cada resultado.
- RV-5.01: con todo desactivado, se recupera **exactamente** el modelo de dos cuerpos.

> [!note] Efecto sobre la validación
> El [[Vector de excentricidad]] deja de conservarse en cuanto hay cualquier perturbación. Los oráculos pasan a ser fórmulas seculares (J2) o comportamiento cualitativo (arrastre).
