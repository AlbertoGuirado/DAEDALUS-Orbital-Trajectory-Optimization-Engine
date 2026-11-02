import numpy as np
import streamlit as st

import app_support as sup

st.caption(
    "Experimento oficial `experiments.validate_two_body`: órbita circular, comparación con la "
    "solución analítica y umbrales de aceptación declarados antes de ejecutar."
)

with st.form("fase1_form"):
    st.subheader("Datos del problema", anchor=False)
    with st.container(horizontal=True):
        radius_km = st.number_input("Radio orbital (km, desde el centro)", min_value=6_400.0,
                                    max_value=500_000.0, value=7_000.0, step=100.0)
        inclination_deg = st.number_input("Inclinación (°)", min_value=0.0, max_value=180.0,
                                          value=28.5, step=0.5)
        periods = st.number_input("Periodos", min_value=1, max_value=50, value=1, step=1)
        samples = st.number_input("Muestras guardadas", min_value=2, max_value=20_000,
                                  value=1001, step=100)
    st.subheader("Integrador", anchor=False)
    with st.container(horizontal=True, vertical_alignment="bottom"):
        method = st.segmented_control("Método", ["RK45", "RK4"], default="RK45", required=True)
        max_step = st.number_input("Paso máximo (s)", min_value=0.1, max_value=3_600.0,
                                   value=60.0, step=10.0)
        rtol = st.number_input("rtol (solo RK45)", min_value=1e-13, max_value=1e-3,
                               value=1e-10, format="%.1e", step=1e-11)
        atol_r = st.number_input("atol posición (m)", min_value=1e-9, max_value=1e3,
                                 value=1e-3, format="%.1e", step=1e-4)
        atol_v = st.number_input("atol velocidad (m/s)", min_value=1e-12, max_value=1.0,
                                 value=1e-6, format="%.1e", step=1e-7)
    st.subheader("Umbrales de aceptación", anchor=False)
    with st.container(horizontal=True):
        acc_pos = st.number_input("Error de posición (m)", min_value=1e-9, value=1.0, format="%.2e")
        acc_vel = st.number_input("Error de velocidad (m/s)", min_value=1e-12, value=1e-3, format="%.2e")
        acc_energy = st.number_input("Deriva de energía", min_value=1e-16, value=1e-7, format="%.1e")
        acc_h = st.number_input("Deriva de momento angular", min_value=1e-16, value=1e-7, format="%.1e")
        acc_e = st.number_input("Deriva de excentricidad", min_value=1e-16, value=1e-7, format="%.1e")
    submitted = st.form_submit_button("Ejecutar experimento", type="primary", icon=":material/play_arrow:")

if submitted:
    config = {
        "radius_m": float(radius_km) * 1000,
        "inclination_rad": float(np.deg2rad(inclination_deg)),
        "mu_m3_s2": sup.EARTH_MU,
        "periods": int(periods),
        "samples": int(samples),
        "method": method,
        "max_step_s": float(max_step),
        "rtol": float(rtol),
        "atol": [float(atol_r)] * 3 + [float(atol_v)] * 3,
        "acceptance": {
            "max_position_error_m": float(acc_pos),
            "max_velocity_error_m_s": float(acc_vel),
            "max_relative_energy_drift": float(acc_energy),
            "max_relative_angular_momentum_drift": float(acc_h),
            "max_absolute_eccentricity_drift": float(acc_e),
        },
    }
    try:
        with st.spinner("Propagando..."):
            st.session_state.fase1 = sup.run_phase1(config)
    except (ValueError, RuntimeError) as exc:
        st.session_state.pop("fase1", None)
        st.error(str(exc), title="No se pudo ejecutar el experimento", icon=":material/error:")

data = st.session_state.get("fase1")
if data is None:
    st.info("Introduce los datos y pulsa **Ejecutar experimento**.", icon=":material/info:")
    st.stop()

report = data["report"]
if report["validation_passed"]:
    st.success("Las cinco comprobaciones cumplen sus umbrales.", title="Validación superada",
               icon=":material/check_circle:")
else:
    failed = int((~data["checks"]["Cumple"]).sum())
    st.error(f"{failed} comprobación(es) superan su umbral.", title="Validación no superada",
             icon=":material/cancel:")

with st.container(horizontal=True):
    st.metric("Periodo orbital", f"{report['period_s'] / 60:.2f} min", border=True)
    st.metric("Error máx. de posición", f"{report['metrics']['max_position_error_m']:.3e} m", border=True)
    st.metric("Deriva de energía", f"{report['metrics']['max_relative_energy_drift']:.2e}", border=True)
    st.metric("Evaluaciones de la derivada", f"{report['function_evaluations']:,}", border=True)
    st.metric("Tiempo de cálculo", f"{report['computation_time_s'] * 1000:.0f} ms", border=True)

st.subheader("Comprobaciones", anchor=False)
st.dataframe(
    data["checks"], hide_index=True,
    column_config={
        "Valor": st.column_config.NumberColumn(format="%.3e"),
        "Límite": st.column_config.NumberColumn(format="%.1e"),
        "Cumple": st.column_config.CheckboxColumn(),
    },
    alt="Tabla de métricas de validación con su valor, límite y resultado",
)

view = st.segmented_control(
    "Gráfico", ["Órbita 3D", "Errores", "Invariantes"], default="Órbita 3D", required=True,
    label_visibility="collapsed",
)
series = data["series"]
if view == "Órbita 3D":
    st.plotly_chart(sup.orbit_3d(data["orbit_km"], data["reference_km"], report["configuration"]["method"]),
                    alt="Órbita propagada en 3D frente a la solución analítica circular")
elif view == "Errores":
    with st.container(horizontal=True):
        st.altair_chart(sup.line_chart(series, "time_min", "position_error_m", "Tiempo (min)",
                                       "Error de posición (m)"),
                        alt="Error de posición frente al tiempo")
        st.altair_chart(sup.line_chart(series, "time_min", "velocity_error_m_s", "Tiempo (min)",
                                       "Error de velocidad (m/s)", color=sup.SERIES[1]),
                        alt="Error de velocidad frente al tiempo")
else:
    st.altair_chart(sup.drift_chart(series),
                    alt="Deriva de energía, momento angular y excentricidad frente al tiempo, escala logarítmica")

st.subheader("Salidas", anchor=False)
with st.container(horizontal=True):
    st.download_button("report.json", sup.to_json(report), file_name="report.json",
                       mime="application/json", icon=":material/download:")
    st.download_button("trajectory.csv", data["trajectory"].to_csv(index=False), file_name="trajectory.csv",
                       mime="text/csv", icon=":material/download:")
    st.download_button("series.csv", series.to_csv(index=False), file_name="series.csv",
                       mime="text/csv", icon=":material/download:")

show_report = st.toggle("Ver report.json")
if show_report:
    st.json(report, expanded=1)
