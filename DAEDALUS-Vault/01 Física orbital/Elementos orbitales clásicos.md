---
tags:
  - fisica/geometria
aliases:
  - Elementos keplerianos
  - OrbitalElements
  - COE
  - state_to_elements
  - elements_to_state
estado: implementado
fase: 1
codigo: src/orbital_optimizer/orbital/elements.py
---

# Elementos orbitales clásicos

> [!abstract] En una frase
> Seis números que separan la geometría de la órbita (constante en dos cuerpos) de la posición de la nave sobre ella (lo único que cambia con el tiempo).

| Elemento | Describe | ¿Constante en dos cuerpos? |
|---|---|---|
| $a$ — semieje mayor | Tamaño | Sí |
| $e$ — excentricidad | Forma | Sí |
| $i$ — inclinación | Inclinación del plano | Sí |
| $\Omega$ — RAAN | Orientación del plano (nodo ascendente) | Sí |
| $\omega$ — argumento del periapsis | Orientación de la órbita dentro del plano | Sí |
| $\nu$ — anomalía verdadera | Posición de la nave | **No** |

Por eso las misiones se especifican en elementos: la órbita inicial y la objetivo se declaran así, no como vectores de estado.

## Estado → elementos (`state_to_elements`)
1. $\mathbf{h} = \mathbf{r}\times\mathbf{v}$, nodo $\mathbf{n} = \hat{\mathbf{z}}\times\mathbf{h}$
2. [[Vector de excentricidad]] $\mathbf{e}$ y $e = \lVert\mathbf{e}\rVert$
3. $a = -\mu/(2\varepsilon)$ a partir de la [[Energía específica]]
4. $i = \arccos(h_z/h)$
5. $\Omega$, $\omega$, $\nu$ con resolución de cuadrante

> [!tip] Decisión de implementación
> La teoría clásica usa `arccos` + comprobación de signo para resolver el cuadrante. DAEDALUS usa **`atan2` con productos orientados respecto a $\mathbf{h}$** (función `_angle`), que resuelve el cuadrante sin ambigüedad.

## Elementos → estado (`elements_to_state`)
Se construye el estado en el [[Marco perifocal]] y se rota al inercial con $R_z(\Omega)\,R_x(i)\,R_z(\omega)$.

## Validación
Test de ida y vuelta (RV-1.06): estado → elementos → estado recupera el original con `rtol=1e-12`, `atol=1e-7 m` y `1e-10 m/s`. Cerca de singularidades se compara el **estado reconstruido**, no los ángulos sueltos (mal condicionados).

Ver también: [[Casos degenerados de los elementos orbitales]], [[Anomalías orbitales]].

> [!note] Alternativa no adoptada
> Los **elementos equinocciales modificados** no tienen singularidades, a cambio de perder interpretación geométrica. Se reconsiderarán solo si la optimización sufre por discontinuidades del espacio de búsqueda.
