---
tags:
  - validacion
  - metodologia
aliases:
  - Evidence levels
  - Jerarquía de validación
---

# Niveles de evidencia

Vocabulario estricto para hablar del estado de cualquier componente. **No mezclar.**

| Nivel | Significa | No significa |
|---|---|---|
| **Planificado** | No implementado | Que exista ningún resultado |
| **Implementado** | El código existe | Que funcione |
| **Testado** | Un test/comando nombrado se ejecutó y se observó su resultado | Que sea correcto en general |
| **Validado para el caso X** | Se cumplieron criterios cuantitativos para una configuración y referencia definidas | Que valga para otros casos |
| **Benchmarkeado** | Métodos rivales corridos en condiciones comparables documentadas | Que uno sea "el mejor" en abstracto ([[Teorema No Free Lunch]]) |

## Jerarquía de validación (de abajo arriba)
1. Tests de entrada/dominio: dimensiones, finitud, dominio físico, excepciones explícitas.
2. Tests unitarios de fórmulas.
3. Casos de referencia analíticos.
4. Implementación de referencia independiente para resultados no triviales.
5. Invariantes y conservación, con unidades y normalización declaradas.
6. Estudios de convergencia y sensibilidad.
7. Experimento end-to-end reproducible con configuración y resultados guardados.
8. Benchmark de rendimiento — **solo después** de que la corrección esté probada.

> [!danger] Prohibido
> - Inferir que algo funciona leyendo el código.
> - Reportar tests como pasados sin haber observado el código de salida.
> - Presentar números de un caso como cotas universales.

Relacionado: [[Invariantes como oráculo de validación]], [[Principios de ingeniería]].
