import streamlit as st

import app_support as sup

st.caption(
    "Transferencia de Hohmann entre órbitas circulares coplanarias (ecuatoriales) y presupuesto de "
    "propelente con Tsiolkovsky. La órbita final se verifica aplicando los impulsos al propagador (RV-2.01)."
)

with st.form("fase2_form"):
    st.subheader("Datos del problema", anchor=False)
    with st.container(horizontal=True):
        alt1 = st.number_input("Altitud inicial (km)", min_value=100.0, max_value=400_000.0,
                               value=300.0, step=50.0)
        alt2 = st.number_input("Altitud objetivo (km)", min_value=100.0, max_value=400_000.0,
                               value=35_786.0, step=50.0)
        mass = st.number_input("Masa inicial (kg)", min_value=1.0, max_value=1e6, value=2_000.0, step=100.0)
        isp = st.number_input("Impulso específico (s)", min_value=1.0, max_value=10_000.0,
                              value=320.0, step=10.0)
    st.subheader("Propagador para la verificación", anchor=False)
    with st.container(horizontal=True, vertical_alignment="bottom"):
        method = st.segmented_control("Método", ["RK45", "RK4"], default="RK45", required=True)
        max_step = st.number_input("Paso máximo (s)", min_value=0.1, max_value=3_600.0, value=60.0, step=10.0)
        rtol = st.number_input("rtol (solo RK45)", min_value=1e-13, max_value=1e-3, value=1e-10,
                               format="%.1e", step=1e-11)
    submitted = st.form_submit_button("Calcular transferencia", type="primary", icon=":material/play_arrow:")
    st.caption(f"Altitudes medidas sobre el radio ecuatorial terrestre, {sup.EARTH_RADIUS_M / 1000:,.3f} km. "
               "35 786 km corresponde a GEO.")

if submitted:
    try:
        with st.spinner("Calculando y propagando..."):
            st.session_state.fase2 = sup.run_phase2(alt1, alt2, mass, isp, method, max_step, rtol)
    except (ValueError, RuntimeError) as exc:
        st.session_state.pop("fase2", None)
        st.error(str(exc), title="No se pudo calcular la transferencia", icon=":material/error:")

data = st.session_state.get("fase2")
if data is None:
    st.info("Introduce los datos y pulsa **Calcular transferencia**.", icon=":material/info:")
    st.stop()

result = data["result"]
hohmann, prop, check = result["hohmann"], result["propellant"], result["verification"]

if check["passed"]:
    st.success("Aplicando los dos impulsos al propagador se alcanza la órbita circular objetivo.",
               title="RV-2.01 superada", icon=":material/check_circle:")
else:
    st.error("La órbita propagada no alcanza la objetivo dentro de los límites.",
             title="RV-2.01 no superada", icon=":material/cancel:")

st.subheader("Maniobra", anchor=False)
with st.container(horizontal=True):
    st.metric("Δv₁", f"{hohmann['first_delta_v_m_s']:,.1f} m/s", border=True)
    st.metric("Δv₂", f"{hohmann['second_delta_v_m_s']:,.1f} m/s", border=True)
    st.metric("Δv total", f"{hohmann['total_delta_v_m_s']:,.1f} m/s", border=True)
    st.metric("Tiempo de transferencia", f"{hohmann['transfer_time_s'] / 3600:,.2f} h", border=True)

st.subheader("Propelente (Tsiolkovsky)", anchor=False)
with st.container(horizontal=True):
    st.metric("Propelente total", f"{prop['total_kg']:,.1f} kg", border=True)
    st.metric("Masa final", f"{prop['final_mass_kg']:,.1f} kg", border=True)
    st.metric("Primer impulso", f"{prop['first_burn_kg']:,.1f} kg", border=True)
    st.metric("Segundo impulso", f"{prop['second_burn_kg']:,.1f} kg", border=True)

st.subheader("Verificación por propagación", anchor=False)
st.dataframe(
    data["checks"], hide_index=True,
    column_config={
        "Valor": st.column_config.NumberColumn(format="%.3e"),
        "Límite": st.column_config.NumberColumn(format="%.1e"),
        "Cumple": st.column_config.CheckboxColumn(),
    },
    alt="Tabla de comprobaciones de la órbita final propagada",
)

view = st.segmented_control(
    "Gráfico", ["Trayectoria", "Perfil", "Δv frente a relación de radios"], default="Trayectoria",
    required=True, label_visibility="collapsed",
)
if view == "Trayectoria":
    st.plotly_chart(sup.transfer_plane(data),
                    alt="Órbita inicial, elipse de transferencia y órbita objetivo en el plano orbital")
elif view == "Perfil":
    profile = data["profile"]
    with st.container(horizontal=True):
        st.altair_chart(sup.line_chart(profile, "time_min", "altitude_km", "Tiempo (min)", "Altitud (km)"),
                        alt="Altitud a lo largo de la elipse de transferencia")
        st.altair_chart(sup.line_chart(profile, "time_min", "speed_km_s", "Tiempo (min)", "Velocidad (km/s)",
                                       color=sup.SERIES[1]),
                        alt="Velocidad a lo largo de la elipse de transferencia")
else:
    st.altair_chart(sup.sweep_chart(data["sweep"], data["ratio"], hohmann["total_delta_v_m_s"] / 1000),
                    alt="Delta-v total de Hohmann frente a la relación de radios, con el caso actual marcado")

st.subheader("Salidas", anchor=False)
st.download_button("fase2_result.json", sup.to_json(result), file_name="fase2_result.json",
                   mime="application/json", icon=":material/download:")
show_json = st.toggle("Ver resultado completo")
if show_json:
    st.json(result, expanded=1)
