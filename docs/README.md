# Documentación de DAEDALUS

Estado revisado el **16 de octubre de 2026**: baseline `fd9d785` (`v0.1.0`)
y ampliación local aún sin commit en `feature/phase-1-kepler`.
La suite ampliada supera 32 tests: ver el registro de cierre y el cambio
justificado del modelo de error de los tests analíticos.
La Fase 1 implementa dinámica de dos cuerpos; Hohmann y optimización pertenecen
a las siguientes fases.

| Documento | Uso |
| --- | --- |
| [Mapa de carpetas y ficheros](mapa-del-proyecto.md) | Papel de cada parte del repositorio en párrafos breves |
| [Características](01-caracteristicas-del-proyecto.md) | Alcance del proyecto y estado de la arquitectura |
| [Requerimientos](02-requerimientos.md) | Requisitos identificados y criterios de cierre |
| [Fundamentos](03-fundamentos-cientificos-y-algoritmicos.md) | Teoría y propuestas para todas las fases; no es un inventario de capacidades |
| [Contrato de Fase 1](phase1.md) | Ecuaciones, convenciones y umbrales de aceptación |
| [API y datos](phase1-api.md) | Interfaces reales, configuración, resultados y errores |
| [Validación ejecutada](phase1-validation.md) | Ejecución original y reproducción del cierre |
| [Auditoría y cierre](phase1-closeout.md) | Inventario, trazabilidad, discrepancias resueltas y pendientes |

La base científica complementaria está en
[DAEDALUS-Vault](../../DAEDALUS-Vault/00%20Inicio.md), junto al repositorio.
Las notas desarrollan conceptos; `docs/` registra el contrato de implementación.
Ese enlace requiere conservar la disposición local de ambas carpetas.
