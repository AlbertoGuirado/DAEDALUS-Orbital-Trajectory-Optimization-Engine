---
tags:
  - moc
  - algoritmia/numerico
---

# MOC Métodos numéricos

Volver a [[00 Inicio]].

## Representación
- [[Vector de estado]] — `y = [x, y, z, vx, vy, vz]`, reducción a primer orden

## Integración de EDOs
- [[Runge-Kutta 4]] — paso fijo, orden 4, implementado a mano
- [[Dormand-Prince RK45]] — adaptativo, vía `scipy.integrate.solve_ivp`
- [[Tolerancias rtol y atol]] — por qué la tolerancia absoluta es vectorial
- [[Comparativa RK4 frente a RK45]]
- [[Integradores simplécticos]] — alternativa no adoptada

## Error numérico
- [[Análisis de error numérico]] — truncamiento, redondeo, modelado
- [[Orden de convergencia]] — el test del factor 16
- [[Umbral de degeneración numérica]] — `32·eps`

## Resolución de ecuaciones no lineales
- [[Newton-Raphson]]
- [[Método de bisección]] — usado en el oráculo de Kepler de los tests
