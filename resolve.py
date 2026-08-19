from mesh import CoqueMesh
from prop import Propulsion
from params import *
import time
import os
import pickle
import numpy as np
import cma

CHECKPOINT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "checkpoints")
MULTISTART_STATE = os.path.join(CHECKPOINT_DIR, "multistart_state.pkl")
RUN_STATE = os.path.join(CHECKPOINT_DIR, "run_state.pkl")


def _save(path, obj):
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        pickle.dump(obj, fh)
    os.replace(tmp, path)


def _load(path):
    if os.path.exists(path):
        with open(path, "rb") as fh:
            return pickle.load(fh)
    return None


def _save_es(path, es):
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        fh.write(es.pickle_dumps())
    os.replace(tmp, path)


def _load_es(path):
    if os.path.exists(path):
        with open(path, "rb") as fh:
            return pickle.loads(fh.read())
    return None

def _fraction(key):
    bas, haut = VARIABLES_LIBRES[key]
    if haut == bas:
        return 0.5
    return (DEFAUTS_MESH[key] - bas) / (haut - bas)

x0 = np.array([_fraction(key) for key in VARIABLES_LIBRES])
BUDGET = 50000

def score(cursor):
    values = {key: VARIABLES_LIBRES[key][0] + (VARIABLES_LIBRES[key][1] - VARIABLES_LIBRES[key][0])*cursor[i] for i, key in enumerate(VARIABLES_LIBRES)}


    coque = CoqueMesh(**values)
    mesh = coque.generate()

    prop = Propulsion(
        mesh,
        masse=60.0,
        x_cg=0.58 * values["L_COQUE"],
        deadrise_deg=coque.DEADRISE,
        rho=1025.0, nu=1.19e-6,
        eta_prop=0.50,
        eta_chaine=0.85,
    )

    if prop.carene.resoudre():
        print("ok")
        return prop.score() 
    print("non ok")
    return 0


class Compteur:
    def __init__(self, f):
        self.f, self.n = f, 0

    def __call__(self, x):
        self.n += 1
        return self.f(x)

    def reset(self):
        self.n = 0


f = Compteur(score)


def cmaes_un_run(x0, sigma0=0.3, budget=BUDGET, verbose=False, resume_es=None):
    es = resume_es if resume_es is not None else cma.CMAEvolutionStrategy(
        x0,                       # point de depart
        sigma0,                   # pas initial ~ 1/4 de l'etendue des curseurs
        {
            'bounds': [0.0, 1.0],
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
        _save_es(RUN_STATE, es)             # checkpoint apres chaque generation
        if verbose and es.countiter % 150 == 0:
            print(f"      gen {es.countiter:4d} | evals {f.n:5d} | "
                  f"meilleur {-es.result.fbest:8.3f} | sigma {es.sigma:.4f}")
    return np.array(es.result.xbest), -es.result.fbest

def cmaes_multistart(x0, n_runs=16, budget_total=BUDGET):
    """Plusieurs departs = protection contre les optima locaux.

    Reprise automatique : l'etat (run en cours + meilleur resultat) est
    sauvegarde sur disque a chaque generation, dans CHECKPOINT_DIR. Si le
    script est interrompu (Ctrl+C, crash, coupure...), le relancer reprend
    exactement la ou il s'etait arrete au lieu de tout refaire.
    """
    budget_par_run = budget_total // n_runs
    state = _load(MULTISTART_STATE)
    same_setup = (
        state is not None
        and np.array_equal(state["x0"], x0)
        and state["n_runs"] == n_runs
        and state["budget_total"] == budget_total
    )
    if same_setup:
        k_start, best_x, best_f = state["k"], state["best_x"], state["best_f"]
        resume_es = _load_es(RUN_STATE)
        print(f"    reprise : run {k_start+1}/{n_runs}, meilleur score actuel = {best_f:.4f}")
    else:
        k_start, best_x, best_f, resume_es = 0, None, -np.inf, None

    for k in range(k_start, n_runs):
        try:
            x, v = cmaes_un_run(x0, 0.3, budget_par_run, resume_es=resume_es)
        except KeyboardInterrupt:
            print("\n    interrompu : progression sauvegardee, relancez pour reprendre.")
            raise
        resume_es = None
        if os.path.exists(RUN_STATE):
            os.remove(RUN_STATE)

        marque = ""
        if v > best_f:
            best_f, best_x = v, x
            marque = "  <-- meilleur"
        print(f"    run {k+1}/{n_runs} : score = {v:9.4f}{marque}")
        _save(MULTISTART_STATE, {
            "k": k + 1, "best_x": best_x, "best_f": best_f,
            "x0": x0, "n_runs": n_runs, "budget_total": budget_total,
        })

    if os.path.exists(MULTISTART_STATE):
        os.remove(MULTISTART_STATE)
    return best_x, best_f, x0


print(cmaes_multistart(x0))