---
tags:
  - proyecto
  - metodologia
aliases:
  - Principios
---

# Principios de ingeniería

Vinculantes para todo el proyecto; no se negocian sin justificación explícita.

1. **No inventar resultados.** Ningún número en un informe sin una ejecución real y reproducible.
2. **No afirmar superioridad sin benchmark.** Con múltiples semillas y escenarios. Ver [[Metodología de benchmark]], [[Teorema No Free Lunch]].
3. **Validar antes de construir encima.** Contra solución analítica o invariante. Ver [[Invariantes como oráculo de validación]].
4. **Separar física, simulación y optimización.** Ver [[Arquitectura del código]].
5. **Evitar abstracciones innecesarias.**
6. **Documentar las ecuaciones** que implementa cada función, con su fuente.
7. **Tests unitarios para toda la matemática.**
8. **Reproducibilidad por semilla y configuración.**
9. **Unidades explícitas y consistentes.** Ver [[Convención de unidades]].
10. **Explicar las alternativas de diseño** antes de implementar.

> [!danger] Regla de oro de los tests
> Un test **nunca** compara contra valores producidos antes por el propio código: eso solo verifica que el código no ha cambiado, no que sea correcto. Y **ajustar una tolerancia hasta que el test pase invalida el test** — ver [[Tolerancias justificadas]].

## Riesgos que motivan estos principios
- Deriva hacia un proyecto de IA a través del paquete `agent/` → el motor es determinista; la IA, si llega, aislada y justificada.
- Error silencioso en la dinámica que contamina todo lo demás.
- Comparaciones injustas entre optimizadores (presupuestos distintos, pocas semillas, escenarios elegidos a posteriori).
