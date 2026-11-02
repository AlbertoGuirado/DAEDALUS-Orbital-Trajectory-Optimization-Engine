---
tags:
  - fisica/propulsion
aliases:
  - Tsiolkovsky
  - Ecuación del cohete
  - Impulso específico
  - Isp
estado: implementado
fase: 2
codigo: src/orbital_optimizer/orbital/mechanics.py
---

# Ecuación de Tsiolkovsky

$$\Delta v = v_e \ln\frac{m_0}{m_f} = I_{sp}\,g_0 \ln\frac{m_0}{m_f}$$

$$m_{\text{prop}} = m_0\left[1 - \exp\!\left(-\frac{\Delta v}{I_{sp}g_0}\right)\right]$$

Se deriva de la conservación del momento para un cuerpo de masa variable que expulsa propelente a velocidad efectiva constante. $g_0$ es la constante de definición de [[Constantes físicas]], no la gravedad local.

> [!important] La razón de ser del proyecto
> La razón de masas depende exponencialmente de Δv. El ahorro de propelente depende de la masa inicial y del Isp; no implica por sí solo un aumento determinado de carga útil.

## Matiz sobre la función objetivo
Como $\exp$ es monótona, minimizar Δv equivale a minimizar propelente **para masa inicial fija y un mismo $I_{sp}$ en todas las maniobras**, sin etapas ni otros cambios de masa. Con varios sistemas de propulsión, el objetivo debe expresarse directamente en masa. Ver [[Función objetivo]].

## En el roadmap
- RF-2.04: masa de propelente de la [[Transferencia de Hohmann]].
- RV-5.04: con masa variable integrada, $m_f/m_0$ debe coincidir con Tsiolkovsky. Ver [[Encendido de duración finita]].

## Derivación y contrato

Con $v_e=I_{sp}g_0$ constante y sin pérdidas externas durante el impulso,
el balance de momento da $dv=-v_e\,dm/m$. Integrar de $m_0$ a $m_f$ da
$\Delta v=v_e\ln(m_0/m_f)$; despejar $m_f$ produce la expresión anterior.
$m_0,m_f,m_p$ en kg, $I_{sp}$ en s, $g_0=9.80665$ m/s², Δv en m/s.

**Implemented:** `propellant_mass(m0, delta_v, specific_impulse)` calcula
`-m0*expm1(-delta_v/(specific_impulse*STANDARD_GRAVITY))`.
Exige escalares finitos, $m_0>0$, $I_{sp}>0$, $\Delta v\ge0$.
**Known Limitations:** cálculo algebraico separado de la dinámica, sin masa
seca, reservas, pérdidas ni empuje finito. La masa final es teóricamente
positiva, aunque float64 puede redondear el consumo a $m_0$ para razones extremas.
**Planned:** masa variable integrada en Phase 5.

Fuente: [NASA Glenn, Ideal Rocket Equation](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/),
consultada el 2026-11-02. Evidencia: [[Validación de la Fase 2 e integración]].
