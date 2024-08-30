# sescousse
Groundwater model for forest drainage

Principales étapes pour la construction du modèle : 

- `ml_setup.py` : Construction du modèle en régime permanent selon fichier de configuration `sescousse.yml` vers le dossier `ml`
- `ml_transient.py`: Construction du modèle en régime transitoire vers le dossier `ml_transient`
- `pproc.py`: post-traitement d'une simulation en régime transitoire (figures dans le dossier `fig`)
- `pest_setup.py`: préparation des fichiers PEST(++) pour l'estimation des paramètres vers le dossier ``ml_tpl``

L'estimation des paramètres se fait sur un serveur de calcul selon les étapes suivantes : 
- Téléchargement du dossier ``ml_tpl`` vers le serveur (avec ``send2hydrolix.sh``)
- Depuis le serveur, lancement de PEST(++) en parallèle avec `start_pest.py`, qui peut être appelé avec la commande shell `python start_pest.py --algo="glm" --unzip="True"`. L'algorithme IES peut aussi être utilisé. Si des modifications sont effectuées sur le dossier `ml_tpl`sur le serveur, il faut appeler la commande avec `--unzip="False` pour éviter que l'archive soit de nouveau extraite. 
- Depuis le serveur, compression du dossier `master_glm`
- Récupération de l'archive `master_glm.tar.gz` avec `getfromhydrolix.sh`
- Post-traitement du calage avec `pproc_glm.py`/`pproc_ies.py`
- Mise à jour du modèle avec le "meilleur" jeu de paramètre avec `parrep.sh cal.par cal.pst caleval.pst`
- Appel de PEST avec NOTPMAX=0 avec `pestpp-glm caleval.pst`

Le post-traitement et l'analyse des résultats avec :
- `analysis.py` 
