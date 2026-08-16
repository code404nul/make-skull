"""
================================================================================
 Epi'Tortue — Paramètres de conception de la coque
================================================================================
Bateau autonome électrique < 2,4 m (Microtransat), coque à DÉPLACEMENT,
AUTO-REDRESSANTE PAR LA FORME (façon aXatlantic, sans quille lestée).

Objectif : minimiser l'énergie consommée, PUIS vérifier l'auto-redressement.

>>> VERSION SYNCHRONISÉE AVEC mesh.py <<<
Les noms de la catégorie 2 correspondent EXACTEMENT aux arguments de
`generate_mesh()`. On peut donc faire directement :

    from params import DEFAUTS_MESH
    from mesh import generate_mesh
    coque = generate_mesh(**DEFAUTS_MESH)

ou, pour l'optimiseur :

    coque = generate_mesh(**candidat)   # candidat = dict tiré de VARIABLES_LIBRES

TROIS CATÉGORIES
----------------
1) INTOUCHABLES     -> valeur fixe. Le programme ne les change JAMAIS.
2) TOUCHABLES       -> tuple (min, max). L'optimiseur teste dans cet intervalle.
                       Ex : L_COQUE = (1.0, 2.4) -> essaie de 1,0 m à 2,4 m.
3) RÉSULTATS CALCULÉS -> NE sont PAS choisis. On les MESURE sur le maillage
                       obtenu (Cp, Cb, LCB, L_WL, tirant d'eau...). On leur
                       donne juste une plage "saine" de vérification.

Unités SI : mètres (m), kilogrammes (kg), secondes (s), degrés (°) si précisé.
"""

# ##############################################################################
#  1) PARAMÈTRES INTOUCHABLES  (valeurs fixes)
# ##############################################################################

# ---- 1.1 Physique (lois de la nature) ---------------------------------------
RHO_EAU        = 1025.0     # masse volumique de l'eau de mer            [kg/m3]
RHO_AIR        = 1.225      # masse volumique de l'air                   [kg/m3]
G              = 9.81       # accélération de la pesanteur               [m/s2]
NU_EAU         = 1.19e-6    # viscosité cinématique de l'eau (~15 °C)    [m2/s]
MU_EAU         = 1.22e-3    # viscosité dynamique de l'eau               [Pa.s]

# ---- 1.2 Règles du Microtransat (non négociables) ---------------------------
LOA_MAX        = 2.4        # longueur hors-tout maximale autorisée      [m]

# ---- 1.3 Densités des matériaux (structure imposée par S1) ------------------
RHO_MOUSSE     = 45.0       # âme mousse rigide PIR/PU                   [kg/m3]
RHO_STRATIFIE  = 1600.0     # stratifié verre/époxy                      [kg/m3]
EP_STRATIFIE   = 0.003      # épaisseur du stratifié de coque            [m]

# ---- 1.4 Budget de masse embarquée (matériel qu'on DOIT transporter) --------
M_BATTERIES    = 8.0        # pack batteries                             [kg]
M_MOTEUR       = 1.5        # moteur brushless + pod                     [kg]
M_ELECTRONIQUE = 2.0        # calculateur, capteurs, satellite, AIS      [kg]
M_SOLAIRE      = 1.0        # panneaux + régulateur                      [kg]
M_DIVERS       = 1.5        # câblage, fixations, étanchéité, marge      [kg]
# La masse de STRUCTURE (coque) est CALCULÉE par le programme (voir cat. 3).

# ---- 1.5 Cibles de mission (cahier des charges) -----------------------------
MASSE_TOTALE_CIBLE = 25.0   # déplacement visé, à faire flotter          [kg]
MARGE_MASSE        = 0.10   # tolérance sur l'équilibre des masses       [-]
FRANC_BORD_MINI    = 0.08   # franc-bord minimal exigé (réserve de       [m]
                            # flottabilité). Le franc-bord est CALCULÉ ;
                            # s'il tombe sous ce seuil -> candidat rejeté.

# ---- 1.6 Discrétisation du maillage (précision, PAS une variable de forme) --
#      Ces valeurs pilotent la finesse du maillage produit par mesh.py.
#      Elles ne changent PAS la coque : elles changent seulement la précision
#      avec laquelle on la décrit. On les fige pour que deux candidats soient
#      comparables (même erreur de discrétisation sur les volumes).
N_SECTIONS     = 61         # nb de stations le long de la coque         [-]
M_FOND         = 20         # nb de points sur le fond                   [-]
M_BOUCHAIN     = 80         # nb de points sur l'arrondi du bouchain     [-]
M_BORDE        = 20         # nb de points sur le bordé                  [-]
M_PONT         = 9          # nb de points sur la largeur du pont        [-]


# ##############################################################################
#  2) PARAMÈTRES TOUCHABLES  (intervalles (min, max) à explorer)
#     /!\ Les NOMS sont ceux des arguments de generate_mesh() (mesh.py).
# ##############################################################################

# ---- 2.1 Dimensions principales de la coque ---------------------------------
L_COQUE        = (1.00, 2.40)   # longueur HORS-TOUT de coque            [m]
                                # (bornée par LOA_MAX ; la longueur à la
                                #  FLOTTAISON L_WL est un RÉSULTAT, cat. 3)
B_MAX          = (0.15, 0.60)   # largeur maxi, mesurée au pont          [m]
CREUX          = (0.15, 0.55)   # hauteur de coque quille -> pont        [m]
                                # -> TIRANT D'EAU et FRANC-BORD en
                                #    DÉCOULENT (flottaison), voir cat. 3.

# ---- 2.2 Forme de la SECTION transversale (Bézier rationnelle) --------------
#      Ces manettes DESSINENT la demi-section de référence (demi_section_gen).
#      Les coefficients de forme (cat. 3) en découleront automatiquement.
DEADRISE       = (0.0, 45.0)    # angle de V du fond (0 = fond plat)     [°]
FLARE          = (-10.0, 30.0)  # évasement du bordé (0 = flanc vertical)[°]
                                # négatif = rentrant (tumblehome) : aide
                                # l'auto-redressement, à tester.
F_BOUCHAIN     = (0.05, 0.95)   # étendue de l'arrondi du bouchain       [-]
                                # (0 = coin vif ponctuel, 1 = arrondi sur
                                #  toute la section). Bornes strictement
                                #  intérieures : 0 et 1 dégénèrent la Bézier.
W_BOUCHAIN     = (0.20, 3.00)   # poids NURBS du bouchain                [-]
                                # < 1 = arrondi, > 1 = anguleux.

# ---- 2.3 Évolution LONGITUDINALE (avant <-> arrière) ------------------------
#      Ce sont CES manettes qui déterminent Cp et LCB (calculés en cat. 3).
X_MAITRE       = (0.40, 0.70)   # position du maître-bau, fraction de
                                # L_COQUE depuis l'étrave (0 = étrave)   [-]
REMPL_AV       = (0.30, 0.80)   # remplissage AVANT (petit = plein,
                                # entrée d'eau fine)                     [-]
REMPL_AR       = (0.30, 0.90)   # remplissage ARRIÈRE (grand = plein,
                                # sortie d'eau porteuse)                 [-]
B_ETRAVE       = (0.02, 0.30)   # largeur résiduelle à l'étrave,
                                # en fraction de B_MAX                   [-]
B_TABLEAU      = (0.30, 1.00)   # largeur au tableau arrière,
                                # en fraction de B_MAX                   [-]

# ---- 2.4 Rocker (courbure longitudinale de la quille) -----------------------
#      Relèvement de la quille aux extrémités : pilote fortement la surface
#      mouillée, l'assiette et le comportement dans la houle.
ROCKER_AV      = (0.000, 0.150) # relèvement de la quille à l'étrave     [m]
ROCKER_AR      = (0.000, 0.100) # relèvement de la quille au tableau     [m]

# ---- 2.5 Pont bombé (élément CLÉ du redressement sans lest) -----------------
HAUTEUR_BOMBE  = (0.000, 0) # flèche du pont bombé au centre         [m]


# ##############################################################################
#  3) RÉSULTATS CALCULÉS  (NE PAS choisir : on les MESURE sur le maillage)
#     Le programme les calcule à partir de la géométrie (cat. 2), puis
#     vérifie qu'ils restent dans une plage "saine". Ce sont des CONTRÔLES,
#     pas des variables d'entrée.
# ##############################################################################

# ---- 3.1 Coefficients de forme (SORTIES, à calculer) ------------------------
#      Plage saine typique d'une coque à déplacement lente. Si le calcul sort
#      de la plage, la coque est jugée non conforme (candidat rejeté).
PLAGE_CP_PRISMATIQUE = (0.52, 0.68)   # coefficient prismatique  [-]
PLAGE_CB_BLOC        = (0.30, 0.55)   # coefficient de bloc      [-]
PLAGE_CW_FLOTTAISON  = (0.65, 0.90)   # coeff. surface flottaison[-]
PLAGE_LCB            = (-0.08, 0.02)  # centre de carène, fraction de L
                                      # depuis le milieu (- vers l'arrière)[-]

# ---- 3.2 Validité du maillage (préalable à TOUTE mesure) --------------------
#      mesh.controle_maillage() doit donner :
#        is_watertight = True, is_volume = True, euler_number = 2,
#        arêtes ouvertes = 0.  Sinon : candidat rejeté sans calcul.

# ---- 3.3 Autres grandeurs calculées (issues de la simulation) ---------------
#      Elles n'ont pas de valeur ici : le programme les REMPLIT.
#        - L_WL (longueur à la flottaison)  [m]  <- MESURÉE, <= L_COQUE
#        - tirant_d_eau (T)        [m]    (trouvé par équilibre de flottaison
#                                          Archimède : rho*V_immerge = masse)
#        - franc_bord = creux - T  [m]    (doit rester > FRANC_BORD_MINI)
#        - volume_immerge          [m3]   (doit satisfaire Archimède)
#        - surface_mouillee        [m2]   (pilote le frottement -> énergie)
#        - masse_structure         [kg]   (mousse + stratifié)
#        - resistance_totale       [N]
#        - puissance_propulsion    [W]    <- objectif à MINIMISER
#        - courbe_GZ(angle 0..180) [m]    <- doit rester favorable
#                                            (contrainte auto-redressement)


# ##############################################################################
#  4) REGROUPEMENTS PRATIQUES
# ##############################################################################

# Toutes les VRAIES variables libres {nom: (min, max)} -> pour l'optimiseur.
# Les clés sont directement utilisables comme kwargs de generate_mesh().
VARIABLES_LIBRES = {
    "L_COQUE":       L_COQUE,
    "B_MAX":         B_MAX,
    "CREUX":         CREUX,
    "DEADRISE":      DEADRISE,
    "FLARE":         FLARE,
    "F_BOUCHAIN":    F_BOUCHAIN,
    "W_BOUCHAIN":    W_BOUCHAIN,
    "X_MAITRE":      X_MAITRE,
    "REMPL_AV":      REMPL_AV,
    "REMPL_AR":      REMPL_AR,
    "B_ETRAVE":      B_ETRAVE,
    "B_TABLEAU":     B_TABLEAU,
    "ROCKER_AV":     ROCKER_AV,
    "ROCKER_AR":     ROCKER_AR,
    "HAUTEUR_BOMBE": HAUTEUR_BOMBE,
}

# Valeurs par défaut de mesh.py : point de départ / coque de référence.
DEFAUTS_MESH = {
    "L_COQUE":       2.35,
    "B_MAX":         0.45,
    "CREUX":         0.30,
    "DEADRISE":      22.0,
    "FLARE":         10.0,
    "F_BOUCHAIN":    0.45,
    "W_BOUCHAIN":    0.60,
    "X_MAITRE":      0.58,
    "REMPL_AV":      0.55,
    "REMPL_AR":      0.55,
    "B_ETRAVE":      0.04,
    "B_TABLEAU":     0.70,
    "ROCKER_AV":     0.075,
    "ROCKER_AR":     0.025,
    "HAUTEUR_BOMBE": 0.01,
}

# Réglages de maillage à passer tels quels à generate_mesh().
MAILLAGE_FIXE = {
    "N_SECTIONS":  N_SECTIONS,
    "M_FOND":      M_FOND,
    "M_BOUCHAIN":  M_BOUCHAIN,
    "M_BORDE":     M_BORDE,
    "M_PONT":      M_PONT,
}

# Plages "saines" des résultats calculés {nom: (min, max)} -> pour VALIDER.
CONTROLES_RESULTATS = {
    "CP_PRISMATIQUE":  PLAGE_CP_PRISMATIQUE,
    "CB_BLOC":         PLAGE_CB_BLOC,
    "CW_FLOTTAISON":   PLAGE_CW_FLOTTAISON,
    "LCB":             PLAGE_LCB,
}

# Constantes intouchables (référence).
CONSTANTES_FIXES = {
    "RHO_EAU": RHO_EAU, "RHO_AIR": RHO_AIR, "G": G,
    "NU_EAU": NU_EAU, "MU_EAU": MU_EAU, "LOA_MAX": LOA_MAX,
    "RHO_MOUSSE": RHO_MOUSSE, "RHO_STRATIFIE": RHO_STRATIFIE,
    "EP_STRATIFIE": EP_STRATIFIE,
    "M_BATTERIES": M_BATTERIES, "M_MOTEUR": M_MOTEUR,
    "M_ELECTRONIQUE": M_ELECTRONIQUE, "M_SOLAIRE": M_SOLAIRE,
    "M_DIVERS": M_DIVERS,
    "MASSE_TOTALE_CIBLE": MASSE_TOTALE_CIBLE, "MARGE_MASSE": MARGE_MASSE,
    **MAILLAGE_FIXE,
}


# ##############################################################################
#  5) CONTRAINTE DE VALIDITÉ GÉOMÉTRIQUE (à tester AVANT de mailler)
# ##############################################################################
#  Dans demi_section_gen(), le fond (angle DEADRISE depuis l'horizontale) et
#  le bordé (angle FLARE depuis la verticale) doivent se couper ENTRE la
#  quille K=(0,0) et le livet D=(B_MAX/2, CREUX). Si le fond est trop relevé,
#  il passe au-dessus du livet : la section est dégénérée (maillage aberrant,
#  volume faux) alors que le solveur ne lève AUCUNE erreur.
#
#  Condition nécessaire :   tan(DEADRISE) * (B_MAX/2) < CREUX
#
#  -> tirer un candidat au hasard dans VARIABLES_LIBRES NE SUFFIT PAS :
#     il faut le passer par section_valide() et le rejeter/retirer sinon.

from math import radians, tan, cos


def section_valide(B_MAX, CREUX, DEADRISE, FLARE, marge=0.05):
    """True si le fond et le bordé se croisent bien à l'intérieur de la
    section. `marge` = fraction de CREUX gardée sous le livet."""
    demi_B = B_MAX / 2.0
    if tan(radians(DEADRISE)) * demi_B >= CREUX * (1.0 - marge):
        return False
    if cos(radians(DEADRISE + FLARE)) <= 1e-6:   # fond et bordé parallèles
        return False
    return True


def candidat_valide(p):
    """Vérifie un dict de paramètres complet (bornes + géométrie)."""
    for nom, (bas, haut) in VARIABLES_LIBRES.items():
        if not (bas <= p[nom] <= haut):
            return False, f"{nom} = {p[nom]} hors de [{bas}, {haut}]"
    if p["L_COQUE"] > LOA_MAX:
        return False, "L_COQUE dépasse LOA_MAX (Microtransat)"
    if not section_valide(p["B_MAX"], p["CREUX"], p["DEADRISE"], p["FLARE"]):
        return False, "section dégénérée (DEADRISE trop fort pour B_MAX/CREUX)"
    return True, "OK"


# ##############################################################################
#  6) PARAMÈTRES RETIRÉS / PAS ENCORE IMPLÉMENTÉS DANS mesh.py
# ##############################################################################
#  Conservés en commentaire pour mémoire : les laisser dans VARIABLES_LIBRES
#  ferait planter generate_mesh() (kwargs inconnus) ou, pire, l'optimiseur
#  passerait du temps à faire varier des manettes SANS EFFET sur la forme.
#
#  Correspondances (ancien params.py -> mesh.py) :
#     L_WL            -> L_COQUE     (mesh donne la longueur hors-tout ;
#                                     L_WL devient un RÉSULTAT mesuré)
#     ANGLE_DEADRISE  -> DEADRISE
#     ANGLE_FLARE     -> FLARE
#     RAYON_BOUCHAIN  -> F_BOUCHAIN  (étendue de l'arrondi)
#     POIDS_NURBS_1   -> W_BOUCHAIN  (poids de la Bézier rationnelle)
#     POIDS_NURBS_2   -> supprimé    (une seule Bézier dans la section)
#     POS_MAITRE_BAU  -> X_MAITRE
#     REMPLISSAGE_AV  -> REMPL_AV
#     REMPLISSAGE_AR  -> REMPL_AR
#
#  Nouveaux dans mesh.py, absents de l'ancien params.py :
#     B_ETRAVE, B_TABLEAU, ROCKER_AV, ROCKER_AR
#
#  Pas de contrepartie dans mesh.py (à coder si on veut les optimiser) :
#     POIDS_PONT      = (0.5, 2.0)     le pont est une parabole imposée
#                                      (flèche = HAUTEUR_BOMBE, pas de poids)
#     ASYM_PONT       = (0.0, 0.15)    le pont est symétrique en y
#     EP_AILES        = (0.001, 0.005) ailes latérales non modélisées
#     ENVERGURE_AILES = (0.0, 0.15)
#     VOLUME_AILES    = (0.0, 0.010)


if __name__ == "__main__":
    # Auto-contrôle : les défauts de mesh.py tiennent-ils dans les plages ?
    ok, msg = candidat_valide(DEFAUTS_MESH)
    print(f"Coque de référence (défauts mesh.py) : {msg}")
    print(f"Variables libres : {len(VARIABLES_LIBRES)}")
    for nom, (bas, haut) in VARIABLES_LIBRES.items():
        val = DEFAUTS_MESH[nom]
        pos = (val - bas) / (haut - bas)
        print(f"  {nom:<15} {bas:>7.3f} |{'-' * int(pos * 20):<20}| "
              f"{haut:>7.3f}   défaut = {val}")