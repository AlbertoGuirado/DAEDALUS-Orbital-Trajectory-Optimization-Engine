# Documentación de DAEDALUS

Estado revisado el **2 de noviembre de 2026**, sobre `develop`, commit `27b4116`.

- **Implemented:** Phase 1 (dos cuerpos, RK4/RK45, Kepler) y requisitos
  obligatorios de Phase 2 (Hohmann, Tsiolkovsky e impulsos). Suite local: 47 tests OK.
- **Planned:** optimización, restricciones de misión, perturbaciones y planificador;
  bielíptica, cambio de plano y Lambert siguen fuera de la API.
- **Known Limitations:** elipses no radiales, sin superficie ni colisiones,
  masa no integrada; impulsos instantáneos y solo estado final en su wrapper.

El tag `v0.1.0` identifica el baseline histórico, no toda la integración actual.
Resultados, entorno y alcance: [auditoría de integración](phase1-phase2-audit.md).

| Documento | Uso |
| --- | --- |
| [Mapa de carpetas y ficheros](mapa-del-proyecto.md) | Papel de cada parte del repositorio en párrafos breves |
| [Características](01-caracteristicas-del-proyecto.md) | Alcance del proyecto y estado de la arquitectura |
| [Requerimientos](02-requerimientos.md) | Requisitos identificados y criterios de cierre |
| [Fundamentos](03-fundamentos-cientificos-y-algoritmicos.md) | Teoría y propuestas para todas las fases; no es un inventario de capacidades |
| [Contrato de Fase 1](phase1.md) | Ecuaciones, convenciones y umbrales de aceptación |
| [API y datos](phase1-api.md) | Interfaces reales, configuración, resultados y errores |
| [Validación ejecutada](phase1-validation.md) | Ejecución original y reproducción del cierre |
| [Auditoría y cierre](phase1-closeout.md) | Registro histórico de Phase 1 |
| [API de Phase 2](phase2-api.md) | Hohmann, propelente, impulsos y límites |
| [Auditoría de integración](phase1-phase2-audit.md) | Estado actual, evidencia reproducida y discrepancias |
| [Flujo de datos](data-flow-mission-planning.md) | Conexión real entre APIs y orquestación propuesta |
| [Informe visual](informe-visual.html) | Mapa interactivo de capacidades y planes |

La base científica complementaria está en
[DAEDALUS-Vault](../DAEDALUS-Vault/00%20Inicio.md), incluida como copia documental en el repositorio.
Las notas desarrollan conceptos; `docs/` registra el contrato de implementación.
Ese enlace requiere conservar la disposición local de ambas carpetas.
