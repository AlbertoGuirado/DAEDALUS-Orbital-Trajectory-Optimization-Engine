---
tags:
  - algoritmia/optimizacion
  - algoritmia/restricciones
aliases:
  - Restricciones
  - Penalización
  - Constraints
estado: planificado
fase: 4
---

# Manejo de restricciones

Las metaheurísticas operan de forma natural **sin restricciones**. Incorporarlas es una decisión de diseño explícita que afecta materialmente al resultado.

## Estrategias
| Estrategia | Idea | Problema |
|---|---|---|
| **Penalización estática** | $J + \sum_j \rho_j \max(0, g_j)^2$ | $\rho$ pequeño → acepta infactibles; grande → paisaje abrupto |
| **Penalización dinámica** | $\rho$ crece con las generaciones | Un parámetro más que ajustar |
| **Pena de muerte** | Descartar toda infactible | Desperdicia información; falla si la región factible es pequeña |
| **[[Reglas de factibilidad de Deb]]** | Comparar por factibilidad primero | — ✅ **Recomendada** |
| **Reparación** | Proyectar sobre la región factible | Solo si la proyección es barata y bien definida (no es el caso general) |

## Restricciones del proyecto (Fase 4)
| Restricción | Expresión | Naturaleza |
|---|---|---|
| Δv disponible | $\Delta v_{\text{total}} \le \Delta v_{\max}$ | Punto final |
| Propelente | $m_{\text{prop}} \le m_{\text{disp}}$ | Punto final |
| Tiempo de misión | $t_f - t_0 \le t_{\max}$ | Punto final |
| Empuje máximo | $F \le F_{\max}$ | Solo con encendidos finitos |
| **Altitud mínima** | $\min_t r(t) \ge R_\oplus + h_{\min}$ | **De trayectoria** |
| Masa máxima | $m_0 \le m_{\max}$ | Configuración |

> [!danger] La altitud mínima es distinta
> Es una restricción **de trayectoria**: evaluarla solo en los extremos es un error grave y silencioso — la solución parece válida y es físicamente imposible. Además, al ser $r$ un **radio**, el umbral es $R_\oplus + h_{\min}$, no $h_{\min}$ ([[Convención de unidades]]).

## Requisitos
- Toda violación se penaliza o descarta **explícitamente** (RF-4.07) y se registra su **magnitud** (RF-4.08).
- Comparar al menos dos estrategias sobre el mismo problema (RF-4.10).
- Test con solución deliberadamente infactible (RV-4.01) y con periapsis intermedio bajo el umbral (RV-4.02).

Relacionado: [[Formulación del problema de optimización]].
