import streamlit as st
import pandas as pd
import plotly.express as px

# Título y Descripción del Dashboard
st.title("Dashboard sobre el Material Genético de la Rata")
st.write("""
Este dashboard proporciona información sobre el material genético de la rata, incluyendo la estructura de su ADN, sus cromosomas y otros aspectos relevantes.
""")


# Información básica sobre el material genético de la rata
st.header("¿Qué es el material genético de la rata?")
st.write("""
El material genético de la rata está compuesto por un conjunto de cromosomas que contienen los genes responsables del desarrollo y funcionamiento de su organismo. La rata (Rattus norvegicus) es un modelo biológico comúnmente utilizado en investigación genética, debido a su similitud con los humanos en varios aspectos biológicos.
""")

# Imagen representativa del genoma de la rata (si tienes una)
st.image('assets/rat_genome.png', caption='Esquema del genoma de la rata', use_column_width=True)

# Información sobre los cromosomas
st.header("Cromosomas de la Rata")
st.write("""
La rata tiene un total de 42 cromosomas, distribuidos en 21 pares. Este conjunto incluye cromosomas autosómicos y sexuales. La estructura del ADN en estos cromosomas es clave para comprender la genética de la rata y su utilidad en la investigación biomédica.
""")

# Datos de ejemplo: Genómica de la Rata (p.ej. genes importantes)
data = {
    'Gen': ['BRCA1', 'TP53', 'AKT1', 'EGFR', 'MYC'],
    'Función': ['Reparación del ADN', 'Supresor de tumores', 'Señalización celular', 'Crecimiento celular', 'Regulación de la transcripción'],
    'Enfermedades asociadas': ['Cáncer', 'Cáncer', 'Cáncer, diabetes', 'Cáncer', 'Cáncer']
}

# Convertir a DataFrame
df = pd.DataFrame(data)

# Mostrar tabla de datos de genes
st.subheader("Ejemplos de genes en la rata y sus funciones")
st.write(df)

# Gráfico interactivo sobre la distribución de genes en los cromosomas de la rata
fig = px.bar(df, x='Gen', y='Enfermedades asociadas', color='Función', title="Distribución de Funciones Genéticas en la Rata")
st.plotly_chart(fig)

# Crear una sección sobre aplicaciones del estudio genético de la rata
st.header("Aplicaciones del Estudio Genético de la Rata")
st.write("""
El estudio del material genético de la rata ha llevado a avances importantes en la biomedicina, como la comprensión de enfermedades humanas, el desarrollo de terapias génicas, y la mejora de tratamientos para enfermedades como el cáncer, la diabetes y enfermedades neurodegenerativas.
""")

# Sección para exploración interactiva
st.header("Explora los datos genéticos de la rata")
option = st.selectbox("Selecciona un gen para ver más información", df['Gen'])

if option:
    selected_data = df[df['Gen'] == option]
    st.write(f"**Gen Seleccionado:** {selected_data['Gen'].values[0]}")
    st.write(f"**Función:** {selected_data['Función'].values[0]}")
    st.write(f"**Enfermedades asociadas:** {selected_data['Enfermedades asociadas'].values[0]}")

