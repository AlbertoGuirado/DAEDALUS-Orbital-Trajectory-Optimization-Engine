---
tags:
  - fisica/invariantes
aliases:
  - ε
  - Energía orbital
  - specific_energy
estado: implementado
fase: 1
codigo: src/orbital_optimizer/orbital/mechanics.py
---

# Energía específica

$$\varepsilon = \frac{v^{2}}{2} - \frac{\mu}{r} = -\frac{\mu}{2a} \qquad [\text{J/kg}]$$

Suma de energía cinética y potencial por unidad de masa. Se **conserva** en el [[Problema de dos cuerpos]].

## Clasificación de la órbita por su signo
| $\varepsilon$ | Órbita |
|---|---|
| $< 0$ | Elipse (cerrada, ligada) — único dominio de la Fase 1 |
| $= 0$ | Parábola |
| $> 0$ | Hipérbola |

De $\varepsilon = -\mu/2a$ se despeja el semieje mayor $a = -\mu/(2\varepsilon)$, usado en `state_to_elements`, y la [[Ecuación de vis-viva]].

## Como herramienta de validación
Detecta especialmente bien la **deriva secular del integrador** y errores en la magnitud de la aceleración. Criterio Fase 1: deriva relativa < 10⁻⁷ en 10 periodos de una elipse inclinada con $e=0.2$. Ver [[Invariantes como oráculo de validación]].

> [!warning]
> Conservar bien la energía **no basta**: un integrador puede hacerlo y a la vez hacer precesar la línea de apsides. Por eso se comprueban también [[Momento angular específico]] y [[Vector de excentricidad]].
