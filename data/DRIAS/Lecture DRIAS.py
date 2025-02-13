import os
import pandas as pd

def process_climate_file(input_path, output_path):
    try:
        # Ouvrir et lire le fichier ligne par ligne
        with open(input_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        
        # Trouver la ligne où commencent les données (la première ligne qui débute par un chiffre)
        start_line = None
        for i, line in enumerate(lines):
            if line.strip() and line[0].isdigit():
                start_line = i
                break
        
        if start_line is None:
            raise ValueError("Aucune donnée trouvée dans le fichier.")
        
        # Lire les données à partir de la ligne détectée
        column_names = [
            "Date", "Latitude", "Longitude", "tasminAdjust", "tasmaxAdjust", "tasAdjust",
            "prtotAdjust", "prsnAdjust", "hussAdjust", "sfcWindAdjust", "evspsblpotAdjust"
        ]
        df = pd.read_csv(input_path, sep=r'\s+', skiprows=start_line, header=None, names=column_names)
        
        # Vérifier que les colonnes nécessaires existent
        if 'prtotAdjust' not in df.columns or 'evspsblpotAdjust' not in df.columns:
            raise ValueError("Les colonnes 'prtotAdjust' ou 'evspsblpotAdjust' sont manquantes.")
        
        # Convertir le format de la colonne "Date" de YYYYMMDD à YYYY-MM-DD
        df['Date'] = pd.to_datetime(df['Date'], format='%Y%m%d', errors='coerce')
        
        # Identifier et signaler les dates invalides
        invalid_dates = df[df['Date'].isna()]
        if not invalid_dates.empty:
            print(f"Lignes avec des dates invalides détectées dans {input_path}:")
            print(invalid_dates)
        
        # Supprimer les lignes avec des dates invalides
        df = df.dropna(subset=['Date'])
        
        # Filtrer les colonnes d'intérêt
        df_filtered = df[["Date", "prtotAdjust", "evspsblpotAdjust"]].copy()
        df_filtered.columns = ["Date", "P", "ETP_hargreave"]
        
        # Convertir les données de kg/m²/s à mm/j
        df_filtered["P"] = df_filtered["P"] * 86400
        df_filtered["ETP_hargreave"] = df_filtered["ETP_hargreave"] * 86400
        
        # Format décimal
        df_filtered['P'] = df_filtered['P'].apply(lambda x: f"{x:.2f}")
        df_filtered['ETP_hargreave'] = df_filtered['ETP_hargreave'].apply(lambda x: f"{x:.2f}")
        
        # Sauvegarder le fichier nettoyé
        df_filtered.to_csv(output_path, index=False, sep="\t")
        print(f"Fichier traité et sauvegardé : {output_path}")
    
    except Exception as e:
        print(f"Erreur lors du traitement de {input_path}: {e}")

# Dossier contenant les fichiers
input_directory = r"C:\Users\User\Documents\Travail\Modèle Sescousse\Projections climatiques\DRIAS source"
output_directory = r"C:\Users\User\Documents\Travail\Modèle Sescousse\Projections climatiques\DRIAS"
os.makedirs(output_directory, exist_ok=True)  # Crée le dossier de sortie s'il n'existe pas

# Liste des fichiers à traiter
file_list = [
    'France_CNRM-CERFACS-CNRM-CM5_KNMI-RACMO22E_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19500101-20051231.txt',
    'France_ICHEC-EC-EARTH_KNMI-RACMO22E_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19500101-20051231.txt',
    'France_ICHEC-EC-EARTH_KNMI-RACMO22E_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-21001231.txt',
    'France_ICHEC-EC-EARTH_SMHI-RCA4_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19700101-20051231.txt',
    'France_ICHEC-EC-EARTH_SMHI-RCA4_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-21001231.txt',
    'France_IPSL-IPSL-CM5A-MR_IPSL-WRF381P_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19510101-20051231.txt',
    'France_IPSL-IPSL-CM5A-MR_IPSL-WRF381P_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-21001231.txt',
    'France_IPSL-IPSL-CM5A-MR_SMHI-RCA4_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19700101-20051231.txt',
    'France_IPSL-IPSL-CM5A-MR_SMHI-RCA4_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-21001231.txt',
    'France_MOHC-HadGEM2-ES_CLMcom-CCLM4-8-17_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19500101-20051231.txt',
    'France_MOHC-HadGEM2-ES_CLMcom-CCLM4-8-17_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-20991231.txt',
    'France_MOHC-HadGEM2-ES_ICTP-RegCM4-6_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19710101-20051231.txt',
    'France_MOHC-HadGEM2-ES_ICTP-RegCM4-6_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-20991231.txt',
    'France_MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19500101-20051231.txt',
    'France_MPI-M-MPI-ESM-LR_CLMcom-CCLM4-8-17_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-21001231.txt',
    'France_MPI-M-MPI-ESM-LR_MPI-CSC-REMO2009_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19500101-20051231.txt',
    'France_MPI-M-MPI-ESM-LR_MPI-CSC-REMO2009_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-21001231.txt',
    'France_NCC-NorESM1-M_DMI-HIRHAM5_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19510101-20051231.txt',
    'France_NCC-NorESM1-M_DMI-HIRHAM5_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-21001231.txt',
    'France_NCC-NorESM1-M_GERICS-REMO2015_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19500101-20051231.txt',
    'France_NCC-NorESM1-M_GERICS-REMO2015_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-21001231.txt',
]

# Parcourir les fichiers et appliquer le traitement
for filename in file_list:
    input_path = os.path.join(input_directory, filename)
    output_path = os.path.join(output_directory, f"processed_{filename}")
    
    if os.path.exists(input_path):
        process_climate_file(input_path, output_path)
    else:
        print(f"Fichier non trouvé : {input_path}")

#%% Lecture RACMO22E CNRM rcp8.5

def process_climate_file_debug(input_path, output_path):
    try:
        # Ouvrir et lire le fichier ligne par ligne
        with open(input_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
        
        # Afficher les premières lignes pour inspection
        print("Aperçu des premières lignes du fichier :")
        print("\n".join(lines[:20]))  # Augmenter l'aperçu pour mieux comprendre la structure
        
        # Trouver la ligne où commencent les données (en supposant que les en-têtes sont sur une ligne avant)
        start_line = None
        for i, line in enumerate(lines):
            if line.strip() and line[0].isdigit():  # La première ligne de données commence par un chiffre
                start_line = i
                break
        
        if start_line is None:
            raise ValueError("Aucune donnée trouvée dans le fichier.")
        
        print(f"Les données commencent à la ligne {start_line + 1}.")

        # Ajuster les colonnes manuellement en cas de décalage des en-têtes
        column_names = [
            "Date", "Latitude", "Longitude", "tasminAdjust", "tasmaxAdjust", "tasAdjust",
            "prtotAdjust", "prsnAdjust", "hussAdjust", "sfcWindAdjust", "evspsblpotAdjust"
        ]
        
        # Lire les données en ignorant les lignes avant la ligne de début
        df = pd.read_csv(input_path, sep=r'\s+', skiprows=start_line, header=None, names=column_names, engine='python')

        # Afficher les colonnes détectées
        print("Colonnes détectées :", df.columns.tolist())
        
        # Vérifier la colonne 'Date'
        if 'Date' not in df.columns:
            raise ValueError("La colonne 'Date' n'a pas été correctement détectée.")
        
        # Convertir le format de la colonne "Date"
        df['Date'] = pd.to_datetime(df['Date'], format='%Y%m%d', errors='coerce')
        
        # Identifier les dates invalides
        invalid_dates = df[df['Date'].isna()]
        if not invalid_dates.empty:
            print("Lignes avec des dates invalides détectées :")
            print(invalid_dates)

        # Supprimer les dates invalides
        df = df.dropna(subset=['Date'])
        
        # Filtrer les colonnes d'intérêt
        df_filtered = df[["Date", "prtotAdjust", "evspsblpotAdjust"]].copy()
        df_filtered.columns = ["Date", "P", "ETP_hargreave"]
        
        # Convertir les unités (de secondes à mm/jour pour les précipitations et ETP)
        df_filtered["P"] = df_filtered["P"] * 86400  # Conversion des précipitations
        df_filtered["ETP_hargreave"] = df_filtered["ETP_hargreave"] * 86400  # Conversion ETP
        
        # Format décimal avec 2 décimales
        df_filtered['P'] = df_filtered['P'].apply(lambda x: f"{x:.2f}")
        df_filtered['ETP_hargreave'] = df_filtered['ETP_hargreave'].apply(lambda x: f"{x:.2f}")
        
        # Sauvegarder le fichier nettoyé
        df_filtered.to_csv(output_path, index=False, sep="\t")
        print(f"Fichier traité et sauvegardé : {output_path}")
    
    except Exception as e:
        print(f"Erreur lors du traitement de {input_path}: {e}")

# Appliquer spécifiquement au fichier problématique
input_path = r"C:\\Users\\User\\Documents\\Travail\\Modèle Sescousse\\Projections climatiques\\DRIAS source\\France_CNRM-CERFACS-CNRM-CM5_KNMI-RACMO22E_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19500101-20051231.txt"
output_path = r"C:\\Users\\User\\Documents\\Travail\\Modèle Sescousse\\Projections climatiques\\DRIAS\\processed_France_CNRM-CERFACS-CNRM-CM5_KNMI-RACMO22E_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-21001231.txt"

process_climate_file_debug(input_path, output_path)

# %% ALADIN 

def process_file(input_path, output_path):
    try:
        # Lire le fichier tout en ignorant les lignes d'entête
        with open(input_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        # Rechercher la première ligne contenant les données (par exemple celle qui commence par des valeurs numériques)
        start_line = None
        for i, line in enumerate(lines):
            if line.strip() and line[0].isdigit():  # Lignes de données
                start_line = i
                break

        if start_line is None:
            raise ValueError("Aucune ligne de données trouvée.")

        # Charger les données à partir de la ligne détectée
        df = pd.read_csv(input_path, sep=',', skiprows=start_line, header=None)

        # Définir les noms de colonnes pour les données
        df.columns = [
            "LambertX", "LambertY", "Date", "tasminAdjust", "tasmaxAdjust", "tasAdjust",
            "prtotAdjust", "prsnAdjust", "hussAdjust", "rsdsAdjust", "rldsAdjust",
            "sfcWindAdjust", "evspsblpotAdjust"
        ]

        # Sélectionner les colonnes pertinentes pour le traitement
        df_filtered = df[["Date", "prtotAdjust", "evspsblpotAdjust"]]

        # Renommer les colonnes pour plus de clarté
        df_filtered.columns = ["Date", "P", "ETP_hargreave"]
        
        # Conversion de P_tot de kg/m²/s à mm/j
        df_filtered["P"] = df_filtered["P"] * 86400
        # Format décimal pour P
        df_filtered['P'] = df_filtered['P'].apply(lambda x: f"{x:.2f}")
        
        # Conversion de P_tot de kg/m²/s à mm/j
        df_filtered["ETP_hargreave"] = df_filtered["ETP_hargreave"] * 86400
        # Format décimal pour P
        df_filtered['ETP_hargreave'] = df_filtered['ETP_hargreave'].apply(lambda x: f"{x:.2f}")
        
        # Convertir le format de la colonne "Date" de YYYYMMDD à YYYY-MM-DD
        df_filtered['Date'] = pd.to_datetime(df_filtered['Date'], format='%Y%m%d')

        # Sauvegarder le fichier filtré dans un nouveau fichier texte
        df_filtered.to_csv(output_path, index=False, sep="\t")  # Utilisation de tabulations pour le fichier texte

        print(f"Fichier traité et sauvegardé : {output_path}")
    except Exception as e:
        print(f"Erreur lors du traitement de {input_path}: {e}")



# # Modèle ALADIN historique
process_file(r"C:\\Users\\User\\Documents\\Travail\\Modèle Sescousse\\Projections climatiques\\DRIAS source\\France_CNRM-CERFACS-CNRM-CM5_CNRM-ALADIN63_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19510101-20051231.txt", r"C:\\Users\\User\\Documents\\Travail\\Modèle Sescousse\\Projections climatiques\\DRIAS\France_CNRM-CERFACS-CNRM-CM5_CNRM-ALADIN63_Historical_METEO-FRANCE_ADAMONT-France_SAFRAN_day_19510101-20051231.txt")
# Modèle ALADIN rcp8.5
process_file(r"C:\\Users\\User\\Documents\\Travail\\Modèle Sescousse\\Projections climatiques\\DRIAS source\\France_CNRM-CERFACS-CNRM-CM5_CNRM-ALADIN63_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-21001231.txt", r"C:\\Users\\User\\Documents\\Travail\\Modèle Sescousse\\Projections climatiques\\DRIAS\France_CNRM-CERFACS-CNRM-CM5_CNRM-ALADIN63_rcp8.5_METEO-FRANCE_ADAMONT-France_SAFRAN_day_20060101-21001231.txt")