# Contribuir a DAEDALUS

## Ramas y versiones

- `main`: versión estable y validada.
- `develop`: integración del trabajo de la siguiente entrega.
- `feature/<descripcion>`: una funcionalidad concreta, creada desde `develop`.
- `fix/<descripcion>`: corrección con un test de regresión cuando proceda.
- Tags `vMAJOR.MINOR.PATCH`: hitos reproducibles. `v0.1.0` corresponde a la Fase 1.

La rama inicial `feature/phase-2-hohmann` parte del mismo código validado que
`develop`; no contiene aún una implementación de Hohmann. Las ramas posteriores
se crearán cuando comience su trabajo, para evitar ramas vacías que envejezcan.

Flujo: pull request de `feature/*` a `develop`, revisión de resultados, y pull
request de `develop` a `main` al completar una fase. Se recomienda squash para
funcionalidades y conservar los hitos mediante tags. No reescribir `main` ni
`develop` con force-push. Estas son convenciones del proyecto; no implican que
se hayan activado reglas de protección en la configuración de GitHub.

```powershell
git switch develop
git pull --ff-only
git switch -c feature/nombre-del-cambio
```

## Verificación local

Seguir la instalación del README y ejecutar desde la raíz del repositorio:

```powershell
python -m unittest discover -s tests -v
python -m experiments.validate_two_body --config configs/earth_orbit.json
```

Usar el Python del entorno virtual. El núcleo solo requiere NumPy y SciPy;
las figuras son opcionales. GitHub Actions ejecuta tests en Linux y Windows
con Python 3.11–3.13 y publica los resultados del experimento como artefactos.
Una configuración de CI no equivale a una ejecución ya superada: consultar
el estado de cada run en Actions.

## Organización científica

Mantener el código implementado en `src/orbital_optimizer` con dependencias
unidireccionales: física, geometría orbital y simulación. `src/physics` es una
capa de compatibilidad. Los archivos vacíos originales de fases posteriores
se conservan como estructura histórica; las nuevas implementaciones se añadirán
al paquete `src/orbital_optimizer` cuando llegue su fase.

Documentar ecuaciones, unidades SI, dominio y tolerancias. Incluir referencias
analíticas y tests antes de ampliar el modelo. Fijar semillas solo cuando haya
aleatoriedad; registrar configuración y versiones en los experimentos.
No publicar cifras de benchmarks sin ejecución. Entornos, cachés, builds y
resultados regenerables pertenecen a rutas ignoradas por Git.
