# Outils pour les données IPRAL

## installation

Ces scripts sont faits pour être utilisés sur le mesocentre IPSL. Avant de les utiliser, il faut réaliser ces étapes. Les chemins sont spécifiques à climserv pour une utilisation sur ciclad, il faut utilise le dossier `data` au lieu de `homedata`.

### Téléchargement des scripts

```bash
git clone https://gitlab.in2p3.fr/ipsl/sirta/ipral/ipral-tools.git
```

### Installation de l'environnement python

- Charger une version récente de python

    ```bash
    module load python/3.6-anaconda50
    ```

- Créer un dossier pour stocker les environnements

    ```bash
    mkdir /homedata/your_user/python_envs
    ```

- Créer l'environnement

    ```bash
    cd ipral-tools
    conda create -p /homedata/your_user/python_envs/ipral_tools --file environment.yml
    ```

## Activer l'environnment python

```bash
module load python/3.6-anaconda50
source activate /homedata/your_user/python_envs/ipral_tools
```


## ipral_file_infos.py

Ce script à partir de deux dates va checrcher tous les fichiers L1 IPRAL disponibles et déterminer pour chaque les paramètres suivants

- nom du fichier
- date de la première mesure
- date de la dernière mesure
- nombre de profils

```bash
python ipral_file_infos.py yyyy-mm-dd1 yyyy-mm-dd2
```

- `yyyy-mm-dd1`: première date à analyser
- `yyyy-mm-dd1`: dernière date à analyser

## ipral_1a_bck_corrected

Ce script soustrait le fond de ciel des profils des données IPRAL 1a.

```bash
python ipral_1a_bck_corrected.py yyyy-mm-dd dossier_sortie
```

- `yyyy-mm-dd`: la date du fichier à traiter
- `dossier_sortie` le dossier où enregistrer le fichier créé.
