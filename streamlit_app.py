import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Entrez, SeqIO
import numpy as np
from collections import Counter
from Bio.SeqUtils import molecular_weight, IsoelectricPoint

st.title("**Dashboard Genético**")

# Función para buscar proteínas en GenBank
def buscar_proteinas(query):
    try:
        Entrez.email = "tu_email@example.com"  
        handle = Entrez.esearch(db="protein", term=query, retmax=10)
        record = Entrez.read(handle)
        handle.close()
        if record["Count"] == "0":
            return None
        return record["IdList"]
    except Exception as e:
        return f"Error al buscar en GenBank: {e}"

# Función para calcular propiedades biofísicas
def calcular_propiedades(sequence):
    # Reemplazamos 'X' por 'A' para evitar problemas
    sequence = sequence.replace("X", "A")

    # Peso molecular de la proteína
    mw = molecular_weight(sequence, seq_type='protein')

    # Punto isoeléctrico de la proteína
    try:
        ip = IsoelectricPoint.IsoelectricPoint(sequence)
        pI = ip.pi()  # Método correcto es pi() en lugar de isoelectric_point()
    except ValueError as e:
        st.error(f"Error de valor: {e}")
        pI = None 
    except Exception as e:
        st.error(f"Error calculando el punto isoeléctrico: {e}")
        pI = None 

    return mw, pI

st.title("Búsqueda en GenBank")
nombre = st.text_input("Introduce un nombre o término de búsqueda:", "")

if nombre:
    st.write(f"Buscando en GenBank para: {nombre}...")
    resultado = buscar_proteinas(nombre)
    if resultado:
        st.write("IDs encontrados en GenBank:")
        id_seleccionado = st.selectbox("Selecciona un ID de GenBank", resultado)
        
        if st.button("Obtener Info"):
            id_full = Entrez.efetch(db="protein", id=id_seleccionado, rettype="genbank", retmode="text")
            id_seq = SeqIO.read(id_full, "genbank")
            sequence = id_seq.seq

            # Información general
            st.subheader("Información General:")
            st.write(f"*Acceso*: {id_seleccionado}")
            st.write(f"*Organismo de origen*: {id_seq.annotations.get('organism', 'No disponible')}")
            st.write(f"*Longitud de la secuencia*: {len(sequence)} pares de bases")
            st.write("*Primeros 200 nucleótidos de la secuencia:*")
            st.write(sequence[:200])

            # Propiedades biofísicas
            st.subheader("Propiedades Biofísicas:")
            mw, pI = calcular_propiedades(sequence)
            st.write(f"**Peso Molecular**: {mw:.2f} Da")
            st.write(f"**Punto Isoeléctrico (pI)**: {pI:.2f}" if pI is not None else "**Punto Isoeléctrico (pI)**: No calculable")

            # Composición de aminoácidos
            st.subheader("Composición de Aminoácidos:")
            aminoacidos = Counter(sequence)
            aminoacidos_list = list(aminoacidos.items())
            aminoacidos_list.sort()
            aminoacidos_names, counts = zip(*aminoacidos_list)
            plt.figure(figsize=(10, 6), facecolor='#0E1117')
            sns.set_theme(style="darkgrid")
            sns.barplot(x=list(aminoacidos_names), y=list(counts), palette="Set2")
            plt.title("Composición de Aminoácidos", fontsize=16, color="White")
            plt.xlabel("Aminoácidos", fontsize=12, color="White")
            plt.ylabel("Frecuencia", fontsize=12, color="White")
            plt.xticks(rotation=45, color="white")
            plt.yticks(color="white")
            plt.gca().set_facecolor('#0E1117')
            st.pyplot(plt)

            # Porcentajes de GC
            st.subheader("Porcentajes de GC:")
            count_c = sequence.count("C")
            count_g = sequence.count("G")
            total = len(sequence)
            gc_content = (count_c + count_g) / total * 100
            plt.figure(figsize=(6, 6), facecolor='#0E1117')  
            sns.set_theme(style="darkgrid") 
            wedges, texts, autotexts = plt.pie([gc_content, 100 - gc_content], 
                                                labels=["GC", "Resto"], 
                                                autopct="%1.1f%%", 
                                                colors=sns.color_palette("Set1", 2),
                                                textprops={'color': 'none'},
                                                shadow=False,
                                                wedgeprops={'edgecolor': 'black'})
            for autotext in autotexts:
                        autotext.set_color('white')
                        autotext.set_fontsize(14)
            for text in texts:
                        text.set_color('white')
                        text.set_fontsize(14)
            plt.pie([gc_content, 100 - gc_content], labels=["GC", "Resto"], autopct="%1.1f%%", colors=sns.color_palette("Set1", 2))
            plt.title("Contenido GC de la Proteína", fontsize=16, color='white') 
            plt.gca().set_facecolor('#0E1117')  
            st.pyplot(plt)   


   
    else:
        st.write("No se encontraron resultados para el ID de GenBank proporcionado.")
        st.write("¡Asegúrate de utilizar el ID correcto o el nombre científico para buscar!")
else:
    st.write("Por favor, introduce un ID de GenBank válido para continuar.")


