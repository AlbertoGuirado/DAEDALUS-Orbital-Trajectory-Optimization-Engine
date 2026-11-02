"""DAEDALUS: interfaz para ejecutar la Fase 1 y la Fase 2 y ver sus resultados.

Ejecutar desde la raíz del repositorio:
    .venv/Scripts/python.exe -m streamlit run streamlit_app.py
"""

import streamlit as st

st.set_page_config(page_title="DAEDALUS", page_icon=":material/satellite_alt:", layout="wide")

page = st.navigation(
    [
        st.Page("app_pages/inicio.py", title="Inicio", icon=":material/home:", default=True),
        st.Page("app_pages/fase1.py", title="Fase 1 · Dinámica orbital", icon=":material/public:"),
        st.Page("app_pages/fase2.py", title="Fase 2 · Hohmann", icon=":material/rocket_launch:"),
        st.Page("app_pages/pruebas.py", title="Batería de tests", icon=":material/checklist:"),
    ],
    position="top",
)

st.title(page.title, icon=page.icon)
page.run()
