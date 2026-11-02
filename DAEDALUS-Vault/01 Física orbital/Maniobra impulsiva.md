---
tags:
  - fisica/maniobras
aliases:
  - Δv
  - Delta-v
  - Impulso instantáneo
  - delta_v
  - total_delta_v
estado: implementado
fase: [1, 2]
codigo: src/orbital_optimizer/orbital/mechanics.py
---

# Maniobra impulsiva

> [!abstract] En una frase
> El encendido dura un tiempo despreciable frente al periodo: la posición no cambia y la velocidad salta.

$$\mathbf{r}^{+} = \mathbf{r}^{-}, \qquad \mathbf{v}^{+} = \mathbf{v}^{-} + \Delta\mathbf{v}$$

## Coste
$$\Delta v_{\text{total}} = \sum_{k=1}^{N} \lVert\Delta\mathbf{v}_k\rVert \qquad [\text{m/s}]$$

> [!danger] Suma de módulos, no de vectores
> Dos impulsos opuestos **no se cancelan**: cuestan el doble. Confundirlo produce una función objetivo que premia maniobras absurdas. Test: `test_impulse_costs_do_not_cancel`.

## Qué existe y qué no
| | Estado |
|---|---|
| Calcular el **coste** de un impulso (`delta_v`, `total_delta_v`) | ✅ Implementado |
| **Aplicar** un impulso durante la propagación | Implemented: `simulation/impulsive.py`, `propagate_with_impulses` |

> [!warning] No confundir
> Calcular un Δv **no es** aplicar un burn al propagador. Hohmann (RV-2.01) se comprueba con arcos de vuelo libre encadenados con saltos de velocidad.

## Validez del modelo
La aproximación requiere una duración de encendido despreciable para la dinámica orbital considerada. El código no evalúa esa condición. Los encendidos prolongados requieren otro modelo → [[Encendido de duración finita]] en la Fase 5.

Relacionado: [[Transferencia de Hohmann]], [[Ecuación de Tsiolkovsky]], [[Función objetivo]].

## Contrato actual y límites

**Implemented:** tiempos estrictamente crecientes en `[0,duration]`, en s;
vectores `(3,)` en m/s en el marco inercial del estado. Permite impulsos en
ambos extremos. La duración es positiva. Usa RK4/RK45 entre impulsos.
**Known Limitations:** solo devuelve el estado final; no orienta vectores,
integra masa ni admite tiempos coincidentes. Cada arco exige dominio elíptico;
el salto final no recibe otra validación orbital. Los invariantes de vuelo
libre pueden cambiar en cada salto. **Planned:** empuje finito y restricciones.
Trazabilidad: [[Validación de la Fase 2 e integración]].
