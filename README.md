# Epi'Tortue

Outil de conception et d'optimisation de la coque d'un voilier/bateau
autonome électrique pour le [Microtransat](https://www.microtransat.org/)
(longueur hors-tout < 2,4 m). La coque visée est à déplacement,
**auto-redressante par la forme** (façon aXatlantic), sans quille lestée.

L'objectif est de trouver, parmi toutes les formes de coque possibles, celle
qui **minimise l'énergie de propulsion** tout en respectant les contraintes
du règlement et un budget de masse embarquée fixé.
## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy scipy trimesh cma
```

## Utilisation

Générer et contrôler un maillage avec les paramètres par défaut :

```bash
python3 mesh.py
```

Calculer la puissance de propulsion nécessaire à différentes vitesses pour
la coque de référence :

```bash
python3 prop.py
```

Vérifier que les valeurs par défaut respectent bien les bornes déclarées

## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy scipy trimesh cma
```

## Utilisation

Générer et contrôler un maillage avec les paramètres par défaut :

```bash
python3 mesh.py
```

Calculer la puissance de propulsion nécessaire à différentes vitesses pour
la coque de référence :

```bash
python3 prop.py
```

Vérifier que les valeurs par défaut respectent bien les bornes déclarées
dans `params.py` :

```bash
python3 params.py
```

Lancer l'optimisation CMA-ES multi-départs (peut être interrompue et
reprise, l'état est sauvegardé après chaque génération) :

```bash
python3 resolve.py
```
