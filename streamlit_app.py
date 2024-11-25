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
from Bio.SeqUtils import molecular_weight, IsoelectricPoint
import wikipediaapi
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
import os
from ete3 import Tree, TreeStyle
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
            # Propiedades biofísicas
            st.subheader("Propiedades Biofísicas:")
            mw, pI = calcular_propiedades(sequence)
            st.write(f"**Peso Molecular**: {mw:.2f} Da")
            st.write(f"**Punto Isoeléctrico (pI)**: {pI:.2f}" if pI is not None else "**Punto Isoeléctrico (pI)**: No calculable")
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
            sns.set_theme(style="dark") 
            wedges, texts, autotexts = plt.pie([gc_content, 100 - gc_content], 
                                            labels=["GC", "Resto"], 
                                            autopct="%1.1f%%", 
                                            colors=sns.color_palette("Set1", 1),
                                            shadow=False,
                                            wedgeprops={'edgecolor': 'white'})
            for autotext in autotexts:
                autotext.set_color('white')

            plt.pie([gc_content, 100 - gc_content], labels=["GC", "Resto"], autopct="%1.1f%%", colors=sns.color_palette("Set1", 2))
            plt.title("Contenido GC de la Proteína", fontsize=16, color='white') 
            plt.gca().set_facecolor('#FFFFFF')  
            st.pyplot(plt)   
            #DE APARTIR DE AQUÍ PUEDEN EMPEZAR A ESCRIBIR#
            #3 TABS DE SANGRÍA, USAR st.write Y ESAS MIELDAS#
            #ECHENLE GANAS CABRONES#
            def obtener_secuencias_genbank(accession_numbers):
    """Obtiene las secuencias de GenBank a partir de los números de acceso"""
    secuencias = {}
    
    Entrez.email = "tu_email@example.com"  # Asegúrate de poner tu correo electrónico
    
    for accession in accession_numbers:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = handle.read()
        handle.close()
        
        # Extraemos la secuencia
        start = record.find("ORIGIN") + len("ORIGIN") + 1
        end = record.find("//", start)
        secuencia = record[start:end].replace("\n", "").replace(" ", "")
        
        secuencias[accession] = secuencia
    
    return secuencias

def alinear_secuencias(secuencias):
    """Alinea las secuencias utilizando ClustalW"""
    with open("secuencias.fasta", "w") as f:
        for accession, secuencia in secuencias.items():
            f.write(f">{accession}\n{secuencia}\n")
    
    clustalw_cline = ClustalwCommandline("clustalw2", infile="secuencias.fasta")
    stdout, stderr = clustalw_cline()
    
    # Leer el archivo alineado
    alineamiento = AlignIO.read("secuencias.aln", "clustal")
    return alineamiento

def construir_arbol(alineamiento):
    """Construye el árbol filogenético a partir del alineamiento"""
    from Bio.Phylo.Applications import FasttreeCommandline
    
    # Guardar el alineamiento en formato PHYLIP
    AlignIO.write(alineamiento, "alineamiento.phy", "phylip")
    
    # Usamos FastTree para generar el árbol
    fasttree_cline = FasttreeCommandline(input="alineamiento.phy")
    stdout, stderr = fasttree_cline()
    
    # Leer el árbol generado
    arbol = Phylo.read(stdout, "newick")
    return arbol

def visualizar_arbol(arbol):
    """Visualiza el árbol filogenético utilizando ETE3"""
    newick_tree = arbol.format("newick")
    tree = Tree(newick_tree)
    
    # Configuración de estilo de árbol
    ts = TreeStyle()
    ts.showLeafName = True
    ts.showBranchSupport = True
    
    return tree, ts

def main():
    st.title("Visualización de Árbol Filogenético desde GenBank")
    
    # Solicitar los números de acceso de GenBank
    accession_numbers = st.text_area("Introduce los números de acceso de GenBank (separados por comas):")
    
    if accession_numbers:
        accession_numbers = [x.strip() for x in accession_numbers.split(",")]
        
        with st.spinner("Obteniendo secuencias de GenBank..."):
            secuencias = obtener_secuencias_genbank(accession_numbers)
        
        with st.spinner("Alineando secuencias..."):
            alineamiento = alinear_secuencias(secuencias)
        
        with st.spinner("Construyendo el árbol filogenético..."):
            arbol = construir_arbol(alineamiento)
        
        with st.spinner("Visualizando el árbol..."):
            tree, ts = visualizar_arbol(arbol)
            st.pyplot(tree.render("%%inline", tree_style=ts))

if __name__ == "__main__":
    main()
            
    else:
        st.write("No se encontraron resultados.")
        st.write("¡Asegurate de utilizar el nombre cientifico para buscar!")
else:
    st.write("¡Asegurate de utilizar el nombre cientifico para buscar!")

