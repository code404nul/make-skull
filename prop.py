"""Estimation de la puissance / energie necessaire pour faire avancer une coque.

Chaine de calcul, reduite au strict necessaire :

    masse            -> tirant d'eau        (equilibre vertical seul)
    carene immergee  -> volume, surface mouillee, Lwl, Bwl, Cb
    vitesse V        -> resistance a l'avancement Rt(V)
    Rt, V, rendements-> puissance a l'arbre, energie sur une distance

Rien d'autre n'est calcule : pas de centre de carene, pas de GM, pas de
courbes hydrostatiques, pas de stabilite.

Hypotheses :
  - eau calme, pas de vent, coque nue (pas d'appendices, pas de windage) ;
  - equilibre en enfoncement seulement (assiette statique nulle) ;
  - regime deplacement : ITTC-57 + facteur de forme + terme de vague forfaitaire ;
  - regime planing      : methode de Savitsky (1964) ;
  - entre les deux : interpolation lineaire (zone ou aucun des deux modeles
    n'est valide, valeur indicative).

Le seul terme reellement empirique est la resistance de vague en regime
deplacement (parametres c_vague / r_hump_max) : c'est le point a caler si
une mesure ou une reference existe. Le reste (equilibre, surface mouillee,
frottement ITTC, Savitsky) est deterministe.
"""

import numpy as np
from trimesh import intersections

G = 9.81


# ==========================================================================
#  1. Carene : tirant d'eau d'equilibre et geometrie mouillee
# ==========================================================================
class Carene:
    """Equilibre vertical d'un maillage etanche + geometrie de la partie immergee.

    Repere du maillage : x longitudinal, y transversal, z vertical vers le haut.
    """

    def __init__(self, mesh, masse, rho=1025.0):
        self.mesh = mesh
        self.masse = masse          # [kg] deplacement total (coque + charge)
        self.rho = rho              # [kg/m3] 1025 mer, 998 eau douce
        self._res = None

    # -- coupe horizontale ------------------------------------------------
    def _immersion(self, z):
        """Coupe le maillage au plan z et renvoie (triangles immerges, volume).

        Le volume sous le plan est obtenu par le theoreme de la divergence
        avec le champ F = (0, 0, z - z0) : le couvercle plan (z = z0) donne une
        contribution nulle, il est donc inutile de refermer le maillage.
            V = somme sur les triangles de  aire * n_z * (z_centroide - z0)
        """
        out = intersections.slice_faces_plane(
            self.mesh.vertices, self.mesh.faces,
            plane_normal=np.array([0.0, 0.0, -1.0]),
            plane_origin=np.array([0.0, 0.0, z]))
        v, f = out[0], out[1]      # (sommets, faces[, index]) selon la version
        if len(f) == 0:
            return None, 0.0
        tri = v[f]
        # produit vectoriel = 2 * aire * normale unitaire (normales sortantes)
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        zc = tri[:, :, 2].mean(axis=1)
        vol = float(np.sum(0.5 * n[:, 2] * (zc - z)))
        return tri, abs(vol)

    # -- equilibre vertical ------------------------------------------------
    def resoudre(self, tol=1e-5):
        """Dichotomie sur z tel que rho * volume_immerge(z) = masse."""
        if self._res is not None:
            return self._res

        v_cible = self.masse / self.rho
        z_min, z_max = self.mesh.bounds[0, 2], self.mesh.bounds[1, 2]
        if self._immersion(z_max)[1] < v_cible:
            return None

        lo, hi = z_min, z_max
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            if self._immersion(mid)[1] < v_cible:
                lo = mid
            else:
                hi = mid
        z = 0.5 * (lo + hi)

        tri, vol = self._immersion(z)

        # surface mouillee = aire de tous les triangles immerges (le plan de
        # flottaison n'est pas maille : rien a retrancher)
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        S = float(0.5 * np.linalg.norm(n, axis=1).sum())

        # trace de la flottaison
        pts = tri.reshape(-1, 3)
        flot = pts[pts[:, 2] > z - 1e-6]
        Lwl = float(flot[:, 0].max() - flot[:, 0].min())
        Bwl = float(2.0 * np.abs(flot[:, 1]).max())
        T = float(z - z_min)

        self._res = dict(
            z_flottaison=z,       # [m] cote de la flottaison dans le repere maillage
            T=T,                  # [m] tirant d'eau (point le plus bas de la quille)
            volume=vol,           # [m3] volume de carene
            S=S,                  # [m2] surface mouillee
            Lwl=Lwl,              # [m] longueur a la flottaison
            Bwl=Bwl,              # [m] largeur max a la flottaison
            Cb=vol / (Lwl * Bwl * T),
            x_tableau=float(self.mesh.bounds[1, 0]),   # x du tableau arriere
        )
        return self._res


# ==========================================================================
#  2. Resistance -> puissance -> energie
# ==========================================================================
class Propulsion:
    """Resistance a l'avancement, puissance et energie.

    Parametres
    ----------
    mesh          : trimesh.Trimesh de la coque (etanche)
    masse         : [kg] deplacement total
    x_cg          : [m] abscisse du centre de gravite, repere du maillage
                    (x=0 a l'etrave, x=L au tableau)
    deadrise_deg  : [deg] angle de V du fond au niveau du CG (parametre DEADRISE)
    rho, nu       : eau (1025 kg/m3, 1.19e-6 m2/s en mer a 15 C ;
                          998  kg/m3, 1.14e-6 m2/s en eau douce a 15 C)
    eta_prop      : rendement propulsif global helice+coque (0.40-0.60 en petite taille)
    eta_chaine    : rendement moteur+variateur+reducteur (1.0 = puissance a l'arbre)
    k_forme       : facteur de forme impose ; None = estime depuis la carene
    c_vague       : coefficient du terme de vague en regime deplacement (a caler)
    r_hump_max    : saturation du terme de vague, en fraction du poids (~0.15-0.25)
    b_chine       : [m] largeur au bouchain pour Savitsky ; None = Bwl
    dCf           : supplement de rugosite (0 pour une coque lisse/petite)
    """

    def __init__(self, mesh, masse, x_cg, deadrise_deg,
                 rho=1025.0, nu=1.19e-6,
                 eta_prop=0.50, eta_chaine=1.0,
                 k_forme=None, c_vague=0.008, r_hump_max=0.20,
                 b_chine=None, dCf=0.0):
        self.carene = Carene(mesh, masse, rho)
        self.masse = masse
        self.x_cg = x_cg
        self.beta = float(deadrise_deg)
        self.rho = rho
        self.nu = nu
        self.eta_prop = eta_prop
        self.eta_chaine = eta_chaine
        self.k_forme = k_forme
        self.c_vague = c_vague
        self.r_hump_max = r_hump_max
        self.b_chine = b_chine
        self.dCf = dCf

    # -- frottement --------------------------------------------------------
    def _cf(self, Re):
        """ITTC-57."""
        return 0.075 / (np.log10(max(Re, 1.0e5)) - 2.0) ** 2

    def _b(self):
        return self.b_chine if self.b_chine else self.carene.resoudre()["Bwl"]

    # -- regime deplacement ------------------------------------------------
    def _r_deplacement(self, V):
        h = self.carene.resoudre()
        Re = V * h["Lwl"] / self.nu
        Rf = 0.5 * self.rho * h["S"] * V * V * (self._cf(Re) + self.dCf)

        if self.k_forme is not None:
            k = self.k_forme
        else:   # Watanabe
            k = -0.095 + 25.6 * h["Cb"] / (
                (h["Lwl"] / h["Bwl"]) ** 2 * np.sqrt(h["Bwl"] / h["T"]))
            k = float(np.clip(k, 0.0, 0.5))

        # --- terme de vague : APPROXIMATION FORFAITAIRE, a caler ----------
        # Rw/W croit en Fn^4 puis sature au niveau de la "bosse" (r_hump_max).
        Fnv = V / np.sqrt(G * h["volume"] ** (1.0 / 3.0))
        x = self.c_vague * Fnv ** 4 / self.r_hump_max
        Rw = self.masse * G * self.r_hump_max * x / (1.0 + x)
        return (1.0 + k) * Rf + Rw

    # -- regime planing (Savitsky 1964) ------------------------------------
    def _r_savitsky(self, V):
        h = self.carene.resoudre()
        b = self._b()
        W = self.masse * G
        beta = self.beta
        lcg = h["x_tableau"] - self.x_cg      # distance CG -> tableau arriere
        if lcg <= 0:
            raise ValueError("x_cg doit etre en avant du tableau arriere")

        Cv = V / np.sqrt(G * b)                          # coef. de vitesse
        CLb = W / (0.5 * self.rho * V * V * b * b)       # portance, carene en V

        # CLb = CL0 - 0.0065 * beta * CL0^0.60  ->  CL0 par point fixe
        CL0 = CLb
        for _ in range(100):
            CL0 = CLb + 0.0065 * beta * CL0 ** 0.60

        LAM_MAX = 12.0

        def lam(tau_deg):
            """lambda (longueur mouillee moyenne / b) pour une assiette donnee.

            Renvoie (valeur, None) ou (None, signe) si la solution sort des
            bornes : signe = +1 -> lambda > LAM_MAX (assiette trop faible),
                     signe = -1 -> lambda ~ 0    (assiette trop forte).
            """
            f = lambda L: (tau_deg ** 1.1 *
                           (0.0120 * np.sqrt(L) + 0.0055 * L ** 2.5 / Cv ** 2) - CL0)
            lo, hi = 1e-6, LAM_MAX
            if f(hi) < 0:
                return None, +1.0
            if f(lo) > 0:
                return None, -1.0
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if f(mid) < 0:
                    lo = mid
                else:
                    hi = mid
            return 0.5 * (lo + hi), None

        def residu(tau_deg):
            """Centre de poussee - CG : nul a l'equilibre en tangage."""
            L, borne = lam(tau_deg)
            if L is None:
                return borne * 1e3
            lp = (0.75 - 1.0 / (5.21 * Cv ** 2 / L ** 2 + 2.39)) * L * b
            return lp - lcg

        lo, hi = 0.2, 25.0            # assiette [deg]
        f_lo, f_hi = residu(lo), residu(hi)
        if f_lo * f_hi > 0:
            return None               # pas d'assiette d'equilibre
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            f_mid = residu(mid)
            if f_lo * f_mid <= 0:
                hi, f_hi = mid, f_mid
            else:
                lo, f_lo = mid, f_mid
        tau = 0.5 * (lo + hi)
        L, _borne = lam(tau)
        if L is None:
            return None

        tr, br = np.radians(tau), np.radians(beta)
        # vitesse moyenne sous carene (survitesse due a la portance)
        Vm = V * np.sqrt(max(1e-6,
                             1.0 - (0.0120 * np.sqrt(L) * tau ** 1.1) / (L * np.cos(tr))))
        Cf = self._cf(Vm * L * b / self.nu) + self.dCf
        Df = self.rho * Vm ** 2 * L * b * b * Cf / (2.0 * np.cos(br))
        Rt = W * np.tan(tr) + Df / np.cos(tr)
        # domaine usuel de validite de Savitsky : lambda <= 4, 2 <= tau <= 15 deg
        valide = (L <= 4.0) and (2.0 <= tau <= 15.0)
        return dict(Rt=Rt, tau=tau, lam=L, Cv=Cv, valide=valide,
                    S_mouillee=L * b * b / np.cos(br))

    # -- selection de regime ------------------------------------------------
    def resistance(self, V):
        """Resistance totale [N] a la vitesse V [m/s], + regime utilise."""
        Cv = V / np.sqrt(G * self._b())
        R_dep = self._r_deplacement(V)
        if Cv <= 1.0:
            return R_dep, "deplacement", None
        sav = self._r_savitsky(V)
        if sav is None or sav["lam"] > 8.0:
            return R_dep, "deplacement (Savitsky sans solution)", None
        if Cv >= 2.0:
            note = "planing" if sav["valide"] else "planing (hors domaine strict)"
            return sav["Rt"], note, sav
        w = Cv - 1.0                      # interpolation lineaire sur 1 < Cv < 2
        return (1.0 - w) * R_dep + w * sav["Rt"], "transition (interpolation)", sav

    # -- resultats ----------------------------------------------------------
    def puissance(self, V):
        """Puissance [W] necessaire a la vitesse V [m/s]."""
        Rt, regime, sav = self.resistance(V)
        Pe = Rt * V                                   # puissance effective (remorquage)
        Pa = Pe / self.eta_prop                       # puissance a l'arbre
        Pin = Pa / self.eta_chaine                    # puissance absorbee (batterie)
        out = dict(V=V, Rt=Rt, P_effective=Pe, P_arbre=Pa, P_absorbee=Pin,
                   regime=regime)
        if sav:
            out.update(assiette_deg=sav["tau"], lambda_=sav["lam"])
        return out

    def energie(self, V, distance):
        """Energie [J et Wh] pour parcourir `distance` [m] a la vitesse V."""
        r = self.puissance(V)
        E = r["P_absorbee"] * distance / V
        r.update(distance=distance, duree_s=distance / V, E_J=E, E_Wh=E / 3600.0)
        return r

    def energie_produisable(self):
        mesh = self.carene.mesh
        b_max = float(mesh.bounds[1, 1] - mesh.bounds[0, 1])
        l_coque = float(mesh.bounds[1, 0] - mesh.bounds[0, 0]) - 0.2
        return l_coque * b_max * 615


    def score(self):
        return self.energie_produisable() / self.energie(1.0, 1000)["E_Wh"]

# =======================================================================
#  Exemple d'utilisation
# ==========================================================================
if __name__ == "__main__":
    from mesh import CoqueMesh      # <- adapter au nom de votre module
    from time import time

    start = time()
    L = 2.35
    coque = CoqueMesh(L_COQUE=L)
    mesh = coque.generate()

    prop = Propulsion(
        mesh,
        masse=60.0,               # [kg] deplacement total
        x_cg=0.58 * L,            # [m] centre de gravite depuis l'etrave
        deadrise_deg=coque.DEADRISE,
        rho=1025.0, nu=1.19e-6,
        eta_prop=0.50,
        eta_chaine=0.85,          # moteur + variateur electriques
    )

    h = prop.carene.resoudre()
    print("=== Carene a l'equilibre ==============================")
    print(f"Tirant d'eau        : {h['T']*1000:7.1f} mm")
    print(f"Volume de carene    : {h['volume']*1000:7.2f} L")
    print(f"Surface mouillee    : {h['S']:7.3f} m2")
    print(f"Lwl / Bwl           : {h['Lwl']:.3f} m / {h['Bwl']:.3f} m")
    print(f"Cb                  : {h['Cb']:7.3f}")

    print("\n=== Puissance =========================================")
    print(f"{'V (m/s)':>8} {'V (kn)':>7} {'Rt (N)':>8} {'P arbre':>9} "
          f"{'P batt':>8} {'Wh/km':>7}  regime")
    for V in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]:
        r = prop.energie(V, 1000.0)
        print(f"{V:8.1f} {V*1.94384:7.2f} {r['Rt']:8.1f} "
              f"{r['P_arbre']:8.0f}W {r['P_absorbee']:7.0f}W {r['E_Wh']:7.1f}  "
              f"{r['regime']}")

    print(prop.energie(1.0, 1000.0))

    print(time() - start)