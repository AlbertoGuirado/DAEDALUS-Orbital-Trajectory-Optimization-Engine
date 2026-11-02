---
tags:
  - fisica/propulsion
aliases:
  - Finite burn
  - Masa variable
  - Pérdida por gravedad
estado: planificado
fase: 5
---

# Encendido de duración finita

El empuje actúa durante un tiempo no nulo y se suma a la gravedad:

$$\ddot{\mathbf{r}} = -\frac{\mu\mathbf{r}}{r^{3}} + \frac{F}{m(t)}\,\hat{\mathbf{u}}(t), \qquad \dot{m} = -\frac{F}{I_{sp}g_0}$$

## Consecuencias para el modelo
- El estado pasa de 6 a **7 dimensiones** (la masa entra en el [[Vector de estado]]).
- El sistema **deja de ser autónomo**: el empuje depende de $t$ a través del programa de encendidos.
- Deja de ser conservativo: los invariantes ya no sirven como oráculo durante el encendido.

## Pérdida por gravedad
Un encendido finito siempre es menos eficiente que el impulso equivalente: parte del empuje se gasta contrarrestando la gravedad durante el encendido. Crece con la duración.

## Validación
RV-5.03: en el límite de **empuje muy alto y duración muy corta** debe converger a la [[Maniobra impulsiva]]. Así se verifica que el nuevo modelo *generaliza* al anterior en lugar de contradecirlo.

Relacionado: [[Ecuación de Tsiolkovsky]], [[Perturbaciones orbitales]], [[Control óptimo]].
