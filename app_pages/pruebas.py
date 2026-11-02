import streamlit as st

import app_support as sup

st.caption("Ejecuta `unittest` en un proceso aparte, con el mismo intérprete que la interfaz.")

with st.container(horizontal=True, vertical_alignment="bottom"):
    suite = st.segmented_control("Batería", list(sup.TEST_SUITES), default="Todas", required=True)
    run = st.button("Ejecutar tests", type="primary", icon=":material/play_arrow:")

if run:
    with st.spinner(f"Ejecutando batería: {suite}..."):
        st.session_state.pruebas = {"suite": suite, **sup.run_test_suite(sup.TEST_SUITES[suite])}

data = st.session_state.get("pruebas")
if data is None:
    st.info("Elige una batería y pulsa **Ejecutar tests**.", icon=":material/info:")
    st.stop()

table = data["table"]
correct = int((table["Resultado"] == "Correcto").sum())
if data["passed"]:
    st.success(f"{data['ran']} tests correctos en {data['seconds']:.2f} s.",
               title=f"Batería «{data['suite']}» superada", icon=":material/check_circle:")
else:
    st.error(f"{data['ran'] - correct} de {data['ran']} tests no pasan. Revisa el registro.",
             title=f"Batería «{data['suite']}» con fallos", icon=":material/cancel:")

with st.container(horizontal=True):
    st.metric("Ejecutados", data["ran"], border=True)
    st.metric("Correctos", correct, border=True)
    st.metric("Fallos o errores", data["ran"] - correct, border=True)
    st.metric("Duración", f"{data['seconds']:.2f} s", border=True)

st.dataframe(table, hide_index=True, alt="Resultado de cada test por módulo y clase")

show_log = st.toggle("Ver registro completo")
if show_log:
    st.code(data["log"], language="text")
