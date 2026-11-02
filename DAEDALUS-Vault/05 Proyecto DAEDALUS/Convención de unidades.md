---
tags:
  - proyecto
  - fisica
aliases:
  - Unidades
  - SI
estado: implementado
fase: 1
---

# Convención de unidades

**Todo el interior del sistema trabaja en SI y radianes.** Las conversiones (km, grados, días) ocurren solo en las fronteras de entrada y salida.

| Magnitud | Unidad interna |
|---|---|
| Longitud / radio | m |
| Tiempo | s |
| Velocidad, Δv | m/s |
| Ángulo | rad |
| Parámetro gravitacional μ | m³/s² |
| Energía específica | J/kg = m²/s² |
| Momento angular específico | m²/s |
| Masa | kg |
| Empuje | N |
| Impulso específico | s |

> [!warning] Radio ≠ altitud
> Las posiciones son **radios desde el centro del cuerpo**, no altitudes sobre la superficie. Una órbita "a 622 km" tiene radio $r \approx 7\,000\,000$ m.

> [!quote] Por qué tanta rigidez
> El Mars Climate Orbiter se perdió en 1999 por mezclar unidades imperiales y métricas en el software de navegación.

Consecuencia práctica: los errores de **posición y velocidad se reportan por separado**, porque tienen unidades distintas y no tiene sentido combinarlos en una sola norma. Ver [[Tolerancias rtol y atol]].

Relacionado: [[Vector de estado]], [[Constantes físicas]].
