---
tags:
  - algoritmia/optimizacion
  - algoritmia/metaheuristica
aliases:
  - DE
  - DE/rand/1
  - Evolución diferencial
estado: planificado
fase: 3
codigo: optimization/differential_evolution.py (vacío)
---

# Differential Evolution

> [!abstract] En una frase
> Metaheurística poblacional (Storn y Price, 1997) cuya mutación usa **diferencias entre miembros de la población**, de modo que el tamaño de los saltos se adapta solo a la dispersión actual.

Primer algoritmo del proyecto (RF-3.09) por simplicidad, robustez y pocos parámetros.

## Algoritmo (DE/rand/1/bin)
Para cada individuo $\mathbf{x}_i$ de la generación $g$:

**1. Mutación**
$$\mathbf{v}_i = \mathbf{x}_{r_1} + F\,(\mathbf{x}_{r_2} - \mathbf{x}_{r_3}), \qquad r_1 \ne r_2 \ne r_3 \ne i$$

**2. Cruce binomial**
$$u_{i,j} = \begin{cases} v_{i,j} & \text{si } \text{rand}_j \le CR \;\text{ o }\; j = j_{\text{rand}}\\ x_{i,j} & \text{en otro caso}\end{cases}$$
$j_{\text{rand}}$ garantiza que al menos una componente venga del mutante.

**3. Selección voraz**
$$\mathbf{x}_i^{(g+1)} = \begin{cases}\mathbf{u}_i & \text{si } J(\mathbf{u}_i) \le J(\mathbf{x}_i)\\ \mathbf{x}_i & \text{en otro caso}\end{cases}$$
Elitista por construcción: el mejor nunca se pierde.

> [!tip] Por qué funciona la autoadaptación
> Al principio la población está dispersa → diferencias grandes → exploración. Al converger, se concentra → diferencias pequeñas → explotación fina. Sin programar ningún calendario.

## Parámetros
| Parámetro | Rango habitual | Efecto |
|---|---|---|
| $NP$ población | 5–10 × dimensión | Más exploración, más evaluaciones por generación |
| $F$ factor diferencial | 0.4–0.9 | Alto → exploración; bajo → convergencia |
| $CR$ cruce | 0.1–0.9 | Alto → problemas no separables |

## Variantes
- **DE/best/1**: más rápida, más riesgo de estancamiento.
- **DE/rand/2**: más exploradora.
- Autoadaptativas (ajustan $F$ y $CR$ en ejecución).

Implementar **DE/rand/1 primero** (referencia canónica) y comparar variantes por [[Metodología de benchmark|benchmark]] antes de afirmar cuál conviene.

## Validación prevista
1. Funciones de prueba estándar con óptimo conocido (RV-3.02).
2. Comparar con [[Transferencia de Hohmann]] en el problema donde sirve de óptimo de referencia (RV-3.01); investigar formulación y convergencia ante un incumplimiento.
3. Reproducibilidad por semilla (RV-3.03).

Relacionado: [[Algoritmos genéticos]], [[Optimización por enjambre de partículas]], [[Paisaje de optimización orbital]].
