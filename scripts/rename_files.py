# -*- coding: utf-8 -*-
"""

Created on Thu July 18 13:00:00 2024

@author: Marly


"""

#----------- LIBRARIES -----------#
import pandas as pd
import os

#----------- PATHs -----------#
counts_barcodes_original = "/Users/u_erazo/Documents/LABbook/2024/Namibia_sequencing/emu-combined-abundance-tax_id-counts.tsv"
archivo_referencia = "/Users/u_erazo/Documents/LABbook/2024/Namibia_sequencing/metadata.xlsx" 
# Verificar si el archivo existe
if not os.path.exists(counts_barcodes_original):
    raise FileNotFoundError(f"No se encontró el archivo: {counts_barcodes_original}")

#----------- Cargar el archivo de referencia -----------#
df_referencia = pd.read_excel(archivo_referencia)

# Crear un diccionario para mapear los nombres de las columnas
mapeo_nombres = pd.Series(df_referencia.sample_names.values, index=df_referencia.barcode).to_dict()

# Cargar el archivo de datos
df_datos = pd.read_csv(counts_barcodes_original, sep='\t')

# Renombrar las columnas en df_datos basado en el diccionario de mapeo
nuevas_columnas = ['tax_id'] + [mapeo_nombres.get(col, col) for col in df_datos.columns[1:]]
df_datos.columns = nuevas_columnas

# Guardar el archivo modificado
df_datos.to_csv("/Users/u_erazo/Documents/LABbook/2024/Namibia_sequencing/OTU_table.tsv", sep='\t', index=False)