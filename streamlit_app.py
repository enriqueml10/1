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
def get_scientific_name(common_name):
    # Reemplazamos los espacios por guiones bajos y preparamos la búsqueda en Wikipedia
    query = common_name.replace(" ", "_")
    url = f"https://es.wikipedia.org/wiki/{query}"

    try:
        # Hacemos una solicitud HTTP para obtener el contenido de la página
        response = requests.get(url)
        response.raise_for_status()  # Lanza un error si la respuesta no es 200 OK
        
        # Usamos BeautifulSoup para analizar el contenido HTML de la página
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Buscamos la infobox, donde generalmente se encuentra el nombre científico
        infobox = soup.find('table', {'class': 'infobox'})
        if infobox:
            rows = infobox.find_all('tr')
            for row in rows:
                th = row.find('th')
                td = row.find('td')
                if th and td and 'nombre científico' in th.get_text().lower():
                    scientific_name = td.get_text(strip=True)
                    return f"El nombre científico de {common_name} es: {scientific_name}"

        # Si no se encuentra el nombre científico
        return f"No se encontró el nombre científico para {common_name}."
    except requests.RequestException as e:
        return f"Error al intentar acceder a la página: {e}"

# Función principal que solicita al usuario el nombre común de la especie
def main():
    print("Bienvenido al sistema de búsqueda de nombres científicos.")
    common_name = input("Introduce el nombre común del animal: ")
    result = get_scientific_name(common_name)
    print(result)

if __name__ == "__main__":
    main()
