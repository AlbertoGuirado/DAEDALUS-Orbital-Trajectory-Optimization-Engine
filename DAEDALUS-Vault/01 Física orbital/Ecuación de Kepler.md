---
tags:
  - fisica/geometria
  - algoritmia/numerico
aliases:
  - Kepler
estado: implementado
fase: 1
codigo: src/orbital_optimizer/orbital/kepler.py
---

# Ecuación de Kepler

$$M = E - e\sin E$$

Relaciona la anomalía media $M$ (lineal en el tiempo) con la excéntrica $E$. Es **trascendente en $E$**: no tiene solución cerrada. Históricamente, buena parte del análisis numérico nació intentando resolverla.

## Cómo se resuelve
- **[[Newton-Raphson]]** (lo habitual): $E_{k+1} = E_k - \dfrac{E_k - e\sin E_k - M}{1 - e\cos E_k}$, DAEDALUS parte del objetivo reducido y salvaguarda cada paso con un intervalo de raíz.
- **[[Método de bisección]]** (lo que hace DAEDALUS en tests): más lento pero garantizado.

## Uso en DAEDALUS: producción y oráculo independiente

La ampliación Kepler integrada en `develop` incorpora `solve_kepler(M,e)`
y `propagate_kepler`. La primera conserva vueltas de M y usa Newton con
salvaguarda por bisección. La segunda usa coeficientes cartesianos f/g.
Su implementación y su validación son estados distintos; consultar
[[Cierre documental de la Fase 1]] para el resultado de las pruebas.

`kepler_reference` en `tests/test_propagator.py` genera una elipse plana desde el periapsis resolviendo Kepler por **60 bisecciones**. Es deliberadamente **independiente** de las conversiones de producción y del integrador, para que el test `test_eccentric_motion_against_independent_kepler_solution` compare contra una verdad que no comparte código con lo que se prueba.

## Propagación analítica frente a numérica
| | Kepleriana | Numérica |
|---|---|---|
| Exactitud | Sin truncamiento de integración; errores de raíz y redondeo | Acumula truncamiento |
| Coste | Depende de muestras e iteraciones de raíz | Depende de pasos, tolerancias y muestras |
| Extensible a perturbaciones | No | Sí |
| Papel en el proyecto | API analítica y referencia de validación | Motor numérico extensible |

El motor es numérico porque las fases 4 y 5 necesitan perturbaciones y encendidos finitos.

Relacionado: [[Anomalías orbitales]], [[Invariantes como oráculo de validación]].

## Derivación usada por la propagación

Definir $c=e\cos E_0=1-r_0/a$, $s=e\sin E_0=(\mathbf r_0\cdot\mathbf v_0)/\sqrt{\mu a}$,
$n=\sqrt{\mu/a^3}$ y $x=E-E_0$. Restar la ecuación de Kepler inicial de la
final y expandir $\sin(E_0+x)$ da $nt=x-c\sin x+s(1-\cos x)$.
Su derivada $1-c\cos x+s\sin x=r/a$ es positiva para $e<1$.
De $|x-nt|\le2e$ se obtiene el intervalo que salvaguarda Newton.

La evaluación cartesiana usa $\mathbf r=f\mathbf r_0+g\mathbf v_0$,
$\mathbf v=\dot f\mathbf r_0+\dot g\mathbf v_0$, con
$f=1-a(1-\cos x)/r_0$, $g=t-(x-\sin x)/n$,
$\dot f=-\sqrt{\mu a}\sin x/(rr_0)$ y $\dot g=1-a(1-\cos x)/r$.
Aquí $f,\dot g$ son adimensionales; $g$ está en s y $\dot f$ en s⁻¹.
Se eliminan vueltas completas al evaluar tiempo y anomalía en el código.

**Implemented:** raíz elíptica, máximo interno de 100 iteraciones,
terminación a escala float64 y propagación f/g. **Known Limitations:**
acondicionamiento $dE/dM=1/(1-e\cos E)$ cerca de $e=1$, redondeo de estado
inicial y reducción de grandes anomalías; no hay cota universal de precisión.
**Planned:** ampliar fuerzas o dominio requiere otros métodos.
Esta derivación algebraica corresponde al código, no a una nueva validación
independiente. Evidencia: [[Validación de la Fase 2 e integración]].
