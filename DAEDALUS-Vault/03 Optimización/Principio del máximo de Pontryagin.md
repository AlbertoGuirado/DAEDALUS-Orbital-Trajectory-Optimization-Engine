---
tags:
  - algoritmia/optimizacion
  - algoritmia/control-optimo
aliases:
  - Pontryagin
  - PMP
  - Vector primer
  - Primer vector
  - Lawden
  - Métodos indirectos
estado: planificado
fase: 3
---

# Principio del máximo de Pontryagin

## Hamiltoniano y condiciones necesarias
$$H = L + \boldsymbol{\lambda}^{\mathsf{T}}\mathbf{f}$$

$$\dot{\mathbf{x}} = \frac{\partial H}{\partial\boldsymbol{\lambda}}, \qquad \dot{\boldsymbol{\lambda}} = -\frac{\partial H}{\partial\mathbf{x}}, \qquad \mathbf{u}^{*} = \arg\min_{\mathbf{u}\in\mathcal{U}} H$$

$\boldsymbol{\lambda}$ son los **costados** (multiplicadores de Lagrange dinámicos). El resultado es un problema de contorno de dos puntos en $(\mathbf{x}, \boldsymbol{\lambda})$.

## Teoría del vector primer (Lawden, 1963)
Para mínimo combustible, el **vector primer** $\mathbf{p} = -\boldsymbol{\lambda}_v$ (costado de la velocidad) cumple en una trayectoria impulsiva óptima:

1. $\lVert\mathbf{p}\rVert \le 1$ en todo el intervalo.
2. $\lVert\mathbf{p}\rVert = 1$ exactamente en los instantes de impulso.
3. $\frac{d}{dt}\lVert\mathbf{p}\rVert = 0$ en los impulsos (interiores).
4. Cada impulso apunta en la dirección de $\mathbf{p}$.

> [!important] El uso más valioso en DAEDALUS: verificar
> Estas condiciones permiten **certificar** una solución producida por una metaheurística. Y diagnostican: si $\lVert\mathbf{p}\rVert > 1$ en algún punto, la solución **no es óptima** y añadir un impulso ahí mejoraría el coste. Ninguna metaheurística da ese diagnóstico por sí sola.

## Por qué no como solver principal
- Hay que derivar analíticamente las ecuaciones de costado.
- Radio de convergencia muy estrecho.
- Necesita estimación inicial de costados, magnitudes sin intuición física.

Relacionado: [[Control óptimo]], [[Métodos directos de transcripción]], [[Maniobra impulsiva]].
