#!/usr/bin/env python3
"""
Demo CMA-ES : 15 curseurs flottants, on cherche le score MAXIMUM.

La fonction objectif est deliberement construite pour reproduire ce que tu decris :
  - deterministe (que du calcul mathematique, aucun alea)
  - non-lineaire et multimodale (plein de faux sommets)
  - NON-SEPARABLE : la bonne valeur d'un curseur depend des autres (epistasie)

Le script compare 4 approches sur le meme budget pour montrer pourquoi
CMA-ES gagne, et surtout pourquoi regler les curseurs un par un echoue.

Lancer :  python3 demo_cmaes.py
"""

import time
import numpy as np
import cma

# =============================================================================
# 1. LA FONCTION OBJECTIF  (a remplacer par la tienne)
# =============================================================================

N = 15  # nombre de curseurs

# Le "vrai" reglage optimal, cache. Dans ton cas reel tu ne le connais pas :
# il sert juste ici a verifier que l'algo le retrouve.
rng = np.random.default_rng(42)
X_OPT = rng.uniform(0.15, 0.85, N)

# Matrice de rotation fixe : c'est ELLE qui cree la non-separabilite.
# Apres rotation, chaque axe du probleme melange plusieurs curseurs,
# donc bouger un seul curseur ne peut pas suffire.
Q, _ = np.linalg.qr(rng.standard_normal((N, N)))

# Echelles differentes par axe : certains curseurs comptent beaucoup plus.
SCALES = np.logspace(0, 1.0, N)


def score(x):
    """Score a MAXIMISER. x = vecteur de 15 flottants dans [0, 1].

    Optimum global = 100.0, atteint exactement en x = X_OPT.
    """
    x = np.asarray(x, dtype=float)

    # decalage par rapport a l'optimum, puis rotation -> non-separable
    z = Q @ (x - X_OPT) * SCALES

    # Rastrigin : une cuvette globale + des centaines d'optima locaux
    rastrigin = 10.0 * N + np.sum(z**2 - 10.0 * np.cos(2.0 * np.pi * z))

    # couplage supplementaire entre curseurs voisins :
    # le terme ne s'annule que si les DEUX curseurs sont bien regles ensemble
    couplage = 8.0 * np.sum((x[:-1] - X_OPT[:-1]) * (x[1:] - X_OPT[1:]))

    return 100.0 - rastrigin - abs(couplage)


# --- compteur d'appels, pour comparer les methodes a budget egal -------------
class Compteur:
    def __init__(self, f):
        self.f, self.n = f, 0

    def __call__(self, x):
        self.n += 1
        return self.f(x)

    def reset(self):
        self.n = 0


f = Compteur(score)
BORNES = [0.0, 1.0]
BUDGET = 50000  # nb d'evaluations accordees a chaque methode



# =============================================================================
# 4. CMA-ES : la boucle ask / tell, expliquee
# =============================================================================

def cmaes_un_run(x0, sigma0=0.3, budget=BUDGET, verbose=False):
    es = cma.CMAEvolutionStrategy(
        x0,                       # point de depart
        sigma0,                   # pas initial ~ 1/4 de l'etendue des curseurs
        {
            'bounds': BORNES,
            'maxfevals': budget,
            'popsize': 16,        # candidats par generation (parallelisables)
            'verbose': -9,
            'seed': 1,
        },
    )
    while not es.stop():
        X = es.ask()                        # 16 combinaisons a tester
        F = [-f(xi) for xi in X]            # CMA-ES MINIMISE -> on inverse
        es.tell(X, F)                       # il ajuste moyenne + covariance
        if verbose and es.countiter % 150 == 0:
            print(f"      gen {es.countiter:4d} | evals {f.n:5d} | "
                  f"meilleur {-es.result.fbest:8.3f} | sigma {es.sigma:.4f}")
    return np.array(es.result.xbest), -es.result.fbest


def cmaes_multistart(n_runs=16, budget_total=BUDGET):
    """Plusieurs departs = protection contre les optima locaux."""
    r = np.random.default_rng(7)
    budget_par_run = budget_total // n_runs
    best_x, best_f = None, -np.inf
    for k in range(n_runs):
        x0 = r.uniform(0.2, 0.8, N)
        x, v = cmaes_un_run(x0, 0.3, budget_par_run)
        marque = ""
        if v > best_f:
            best_f, best_x = v, x
            marque = "  <-- meilleur"
        print(f"    run {k+1}/{n_runs} : score = {v:9.4f}{marque}")
    return best_x, best_f, x0


# =============================================================================
# 5. EXECUTION ET COMPARAISON
# =============================================================================

def erreur(x):
    """Distance au vrai optimum (inconnue dans un cas reel)."""
    return np.linalg.norm(np.asarray(x) - X_OPT)


if __name__ == "__main__":
    print("=" * 68)
    print(f"  Optimisation de {N} curseurs flottants dans [0, 1]")
    print(f"  Budget identique pour chaque methode : {BUDGET} evaluations")
    print(f"  Optimum global = 100.000 (cache, sert de verification)")
    print("=" * 68)

    resultats = []


    print("\n[3] CMA-ES, un seul run depuis le centre")
    f.reset(); t = time.time()
    x_1, v_1 = cmaes_un_run(np.full(N, 0.5), 0.3, BUDGET, verbose=True)
    resultats.append(("CMA-ES x1", v_1, f.n, erreur(x_1), time.time() - t))
    print(f"    score = {v_1:.4f}   ({f.n} evaluations)")

    print("\n[4] CMA-ES multi-start (6 departs, meme budget total)")
    f.reset(); t = time.time()
    x_m, v_m, result = cmaes_multistart(16, BUDGET)
    resultats.append(("CMA-ES x6", v_m, f.n, erreur(x_m), time.time() - t))
    print(result)
    # ---- tableau final ----
    print("\n" + "=" * 68)
    print(f"{'Methode':<12}{'Score':>11}{'Evals':>8}{'Dist.opt':>11}{'Temps':>9}")
    print("-" * 68)
    for nom, v, n, e, dt in resultats:
        print(f"{nom:<12}{v:>11.4f}{n:>8}{e:>11.4f}{dt:>8.1f}s")
    print("=" * 68)

    print("\nMeilleurs reglages trouves (CMA-ES multi-start) :")
    for i, (val, opt) in enumerate(zip(x_m, X_OPT)):
        barre = "#" * int(val * 30)
        print(f"  curseur {i+1:2d} : {val:6.4f} |{barre:<30}| vrai = {opt:6.4f}")

    print("""
Lecture des resultats :
  - L'aleatoire couvre large mais n'affine jamais : plafond bas.
  - "Un par un" se bloque et rend la main sans avoir epuise son budget :
    aucun curseur seul n'ameliore plus rien, alors qu'on est loin du sommet.
    C'est exactement le piege que cree la non-separabilite.
  - Un run unique de CMA-ES converge puis s'arrete tout seul (sigma -> 0) :
    le budget restant est perdu, et il reste coince sur un sommet secondaire.
  - Le multi-start recycle ce budget en plusieurs explorations : meilleur
    score final, sans une seule evaluation de plus.
  - Aucune methode n'atteint 100 : sur un paysage aussi accidente, l'objectif
    realiste est "le meilleur sommet trouve", pas l'optimum garanti.

Adaptation a ton cas (3.5 s par evaluation) :
  - remplace score() par ta fonction, garde les curseurs normalises sur [0,1]
    (indispensable : CMA-ES part d'un seul sigma pour les 15 curseurs)
  - dimensionne le budget : 24000 evals = ~23 h en sequentiel,
    ~3 h sur 8 coeurs. Une nuit sequentielle = ~10000 evals.
  - pour paralleliser, remplace la ligne  F = [-f(xi) for xi in X]  par :
        with Pool(8) as pool:
            F = pool.map(objectif_negatif, X)
    (les 16 candidats d'une generation sont independants entre eux)
  - garde un cache dict{tuple(x): score} : ta fonction est deterministe,
    inutile de recalculer deux fois le meme point.
""")