import streamlit as st

st.markdown(
    "Interfaz local del motor DAEDALUS. Cada página recibe los datos del problema, "
    "ejecuta el código validado de `src/orbital_optimizer` y muestra resultados, "
    "comprobaciones de aceptación, gráficos y ficheros descargables."
)

with st.container(horizontal=True):
    with st.container(border=True):
        st.subheader(":material/public: Fase 1", anchor=False)
        st.markdown(
            "Propaga una órbita circular con RK45 o RK4 y la compara con la solución "
            "analítica. Comprueba errores de posición y velocidad y la deriva de los invariantes."
        )
        st.page_link("app_pages/fase1.py", label="Abrir Fase 1", icon=":material/arrow_forward:")
    with st.container(border=True):
        st.subheader(":material/rocket_launch: Fase 2", anchor=False)
        st.markdown(
            "Calcula la transferencia de Hohmann entre dos órbitas circulares y el propelente "
            "con Tsiolkovsky. Verifica la órbita final aplicando los impulsos al propagador (RV-2.01)."
        )
        st.page_link("app_pages/fase2.py", label="Abrir Fase 2", icon=":material/arrow_forward:")
    with st.container(border=True):
        st.subheader(":material/checklist: Tests", anchor=False)
        st.markdown(
            "Ejecuta la batería `unittest` completa, la de la Fase 1 o la de la Fase 2, "
            "y muestra el resultado de cada test."
        )
        st.page_link("app_pages/pruebas.py", label="Abrir tests", icon=":material/arrow_forward:")

st.caption(
    "Modelo: dos cuerpos newtoniano, sin perturbaciones. Unidades internas SI; "
    "km, minutos y grados se usan solo para mostrar datos."
)
