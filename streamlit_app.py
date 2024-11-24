import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Entrez, SeqIO
import numpy as np
import requests
from collections import Counter
from Bio.SeqUtils import gc_fraction
from bs4 import BeautifulSoup
from collections import Counter
st.title("**Dashboard Genético**")
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
def get_scientific_name(common_name):
    # Reemplazamos los espacios por guiones y hacemos la búsqueda en Wikipedia
    query = common_name.replace(" ", "_")
    url = f"https://es.wikipedia.org/wiki/{query}"
    try:
        # Hacemos una solicitud HTTP para obtener el contenido de la página
        response = requests.get(url)
        response.raise_for_status()  # Lanza un error si la respuesta no es 200 OK
        
        # Usamos BeautifulSoup para analizar el contenido HTML de la página
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Buscamos el nombre científico en la página (generalmente está en la infobox)
        infobox = soup.find('table', {'class': 'infobox'})
        if infobox:
            rows = infobox.find_all('tr')
            for row in rows:
                th = row.find('th')
                td = row.find('td')
                if th and td and 'nombre científico' in th.get_text().lower():
                    scientific_name = td.get_text(strip=True)
                    return f"El nombre científico de {common_name} es: {scientific_name}"
        return "No se encontró el nombre científico en la página de Wikipedia."
    except requests.RequestException as e:
        return f"Error al buscar la información: {e}"
# Función principal que solicita al usuario el nombre común de la especie
def main():
    print("Bienvenido al sistema de búsqueda de nombres científicos en Wikipedia.")
    common_name = input("Introduce el nombre común de la especie: ")
    result = get_scientific_name(common_name)
    print(result)
if __name__ == "__main__":
    main()
st.title("Búsqueda en GenBank")
nombre = st.text_input("Introduce el nombre cíentifico de la especie para hacer la búsqueda:", "")
if nombre:
    st.write(f"Buscando en GenBank para: {nombre}...")
    resultado = buscar_proteinas(nombre)
    if resultado:
        st.write("IDs encontrados en GenBank:")
        id_seleccionado = st.selectbox("Selecciona un ID de GenBank", resultado)
        if st.button("Obtener Info"):
            id_full=Entrez.efetch(db="protein", id=id_seleccionado, rettype="genbank", retmode="text")
            id_seq=SeqIO.read(id_full, "genbank")
            sequence=id_seq.seq
            #AQUÍ ES LA INFORMACIÓN GENERAL#
            st.subheader("Información General:")
            st.write(f"*Acceso*: {id(id_full)}")
            st.write(f"*Organismo de origen*: {id_seq.annotations.get('organism', 'No disponible')}")
            st.write(f"*Longitud de la secuencia*: {len(sequence)} pares de bases")
            st.write("*Primeros 200 nucleótidos de la secuencia:*")
            st.write(sequence[:200])
            #AQUÍ LA COMPOSICIÓN DE AMINOÁCIDOS#
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
            #PORCENTAJES DE CG"
            st.subheader("Porcentajes de CG:")
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
            #DE APARTIR DE AQUÍ PUEDEN EMPEZAR A ESCRIBIR#
            #3 TABS DE SANGRÍA, USAR st.write Y ESAS MIELDAS#
            #ECHENLE GANAS CABRONES#
            
    else:
        st.write("No se encontraron resultados.")
        st.write("¡Asegurate de utilizar el nombre cientifico para buscar!")
else:
    st.write("¡Asegurate de utilizar el nombre cientifico para buscar!")
