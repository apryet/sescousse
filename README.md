
# Modèle d'écoulement

## Construction du modèle 

- `ml_setup.py` : Construction du modèle en régime permanent selon fichier de configuration `sescousse.yml` vers le dossier `ml`
- `ml_transient.py`: Construction du modèle en régime transitoire vers le dossier `ml_transient`
- `pproc.py`: post-traitement d'une simulation en régime transitoire (figures dans le dossier `fig`)
- `pest_setup.py`: préparation des fichiers PEST(++) pour l'estimation des paramètres vers le dossier ``ml_tpl``

## Estimation des paramètres

- Téléchargement du dossier ``ml_tpl`` vers le serveur (avec ``send2hydrolix.sh``)
- Depuis le serveur, lancement de PEST(++) en parallèle avec `start_pest.py`, qui peut être appelé avec la commande shell `python start_pest.py --algo="glm" --unzip="True"`. L'algorithme IES peut aussi être utilisé. Si des modifications sont effectuées sur le dossier `ml_tpl`sur le serveur, il faut appeler la commande avec `--unzip="False` pour éviter que l'archive soit de nouveau extraite. 
- Depuis le serveur, compression du dossier `master_glm`
- Récupération de l'archive `master_glm.tar.gz` avec `getfromhydrolix.sh`
- Mise à jour du modèle avec le "meilleur" jeu de paramètre avec `parrep.sh cal.par cal.pst caleval.pst`
- Calcul de la "meilleure solution", avec NOTPMAX=0 ou NOTPMAX=-1 et `++uncertainty(True)` dans `caleval.pst`, puis `pestpp-glm caleval.pst`
- (ou) Calcul de la "meilleure solution" avec incertitudes avec NOTPMAX=-1 et `++uncertainty(True)` dans `caleval.pst`, puis `pestpp-glm caleval.pst`
-  Post-traitement du calage avec `pproc_glm.py`/`pproc_ies.py`

## Post-traitement et l'analyse des résultats avec :
-  `pproc.py` : post-traitement 
-  `analysis.py` : analyse détaillée 

## Comparaison de simulations 
- `gen_sim.py` : génération des simulations alternatives
- `comp_indic.py` : comparaison des résultats selon les indicateurs

# Simulation des forçages (FS4, F ADES) 
Dans le dossier `sim`, avec `sim_levels.py`
