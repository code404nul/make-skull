import numpy as np
import trimesh
import time


class CoqueMesh:
    """Génère et représente le maillage trimesh de la coque d'un bateau.

    Usage :
        coque = CoqueMesh(L_COQUE=2.35, B_MAX=0.45, ...)
        mesh = coque.generate()
        coque.controle_maillage()
    """

    def __init__(
            self,
            # --- dimensions principales ----------------------------------
            L_COQUE=2.35,        # [m] longueur hors-tout de coque
            B_MAX=0.45,          # [m] largeur maxi (au niveau du pont)
            CREUX=0.30,          # [m] hauteur de coque quille -> pont

            # --- forme de la section transversale -------------------------
            DEADRISE=22.0,       # [deg] angle de V du fond (0 = fond plat)
            FLARE=10.0,          # [deg] evasement du borde (0 = flanc vertical)
            F_BOUCHAIN=0.45,     # [0-1] etendue de l'arrondi du bouchain
            W_BOUCHAIN=0.60,     # [>0]  poids NURBS : <1 arrondi, >1 anguleux

            # --- evolution longitudinale ----------------------------------
            X_MAITRE=0.58,       # [0-1] position du maitre-bau (0=etrave)
            REMPL_AV=0.55,       # [0-1] remplissage AVANT  (petit=plein)
            REMPL_AR=0.55,       # [0-1] remplissage ARRIERE (grand=plein)
            B_ETRAVE=0.04,       # [frac de B_MAX] largeur residuelle a l'etrave
            B_TABLEAU=0.70,      # [frac de B_MAX] largeur au tableau arriere
            ROCKER_AV=0.075,     # [m] relevement de la quille a l'etrave
            ROCKER_AR=0.025,     # [m] relevement de la quille au tableau

            # --- pont ------------------------------------------------------
            HAUTEUR_BOMBE=0.01,  # [m] fleche du pont bombe au centre

            # --- discretisation (precision, pas de la forme) ---------------
            N_SECTIONS=61,       # nb de stations le long de la coque
            M_FOND=20,           # nb de points sur le fond
            M_BOUCHAIN=80,       # nb de points sur l'arrondi du bouchain
            M_BORDE=20,          # nb de points sur le borde
            M_PONT=9,            # nb de points sur la largeur du pont
    ):
        self.L_COQUE = L_COQUE
        self.B_MAX = B_MAX
        self.CREUX = CREUX

        self.DEADRISE = DEADRISE
        self.FLARE = FLARE
        self.F_BOUCHAIN = F_BOUCHAIN
        self.W_BOUCHAIN = W_BOUCHAIN

        self.X_MAITRE = X_MAITRE
        self.REMPL_AV = REMPL_AV
        self.REMPL_AR = REMPL_AR
        self.B_ETRAVE = B_ETRAVE
        self.B_TABLEAU = B_TABLEAU
        self.ROCKER_AV = ROCKER_AV
        self.ROCKER_AR = ROCKER_AR

        self.HAUTEUR_BOMBE = HAUTEUR_BOMBE

        self.N_SECTIONS = N_SECTIONS
        self.M_FOND = M_FOND
        self.M_BOUCHAIN = M_BOUCHAIN
        self.M_BORDE = M_BORDE
        self.M_PONT = M_PONT

        self.mesh = None
        self.face_pont = None

    # --- geometrie de base, independante d'une instance -----------------

    @staticmethod
    def bezier_rationnelle(P0, P1, P2, w, n=80):
        """Courbe passant par P0 et P2, tiree vers P1 avec la force w."""
        t = np.linspace(0.0, 1.0, n).reshape(-1, 1)
        b0 = (1 - t) ** 2
        b1 = 2 * t * (1 - t) * w
        b2 = t ** 2
        numerateur = b0 * P0 + b1 * P1 + b2 * P2
        denominateur = b0 + b1 + b2
        return numerateur / denominateur

    @classmethod
    def demi_section_gen(cls, deadrise, flare, f_bouchain, w_bouchain,
                          demi_largeur, creux,
                          m_fond=20, m_bouchain=80, m_borde=20):
        """Demi-section transversale NORMALISEE : y dans [0,1], z dans [0,1].

        K = quille (0,0), D = livet de pont (demi_largeur, creux).
        Le fond monte a l'angle 'deadrise', le borde descend a l'angle 'flare',
        les deux se coupent en C : c'est le bouchain, arrondi par une Bezier
        rationnelle de poids w_bouchain (<1 arrondi, >1 anguleux).
        """
        K = np.array([0.0, 0.0])
        D = np.array([demi_largeur, creux])

        d, f = np.radians(deadrise), np.radians(flare)
        u = np.array([np.cos(d), np.sin(d)])   # direction du fond
        v = np.array([np.sin(f), np.cos(f)])   # direction du borde

        t_s = np.linalg.solve(np.column_stack([u, -v]), D - K)
        C = K + t_s[0] * u                     # coin theorique du bouchain

        B1 = C + f_bouchain * (K - C)          # debut de l'arrondi
        B2 = C + f_bouchain * (D - C)          # fin de l'arrondi

        fond = np.linspace(K, B1, m_fond)
        bouchain = cls.bezier_rationnelle(B1, C, B2, w_bouchain, m_bouchain)
        borde = np.linspace(B2, D, m_borde)

        return np.vstack([fond, bouchain[1:], borde[1:]]) / D

    # --- generation du maillage ------------------------------------------

    def generate(self):
        """Construit et renvoie le maillage trimesh.Trimesh de la coque."""

        demi_B = self.B_MAX / 2.0
        x_maitre = self.X_MAITRE * self.L_COQUE

        # --- 1. profil transversal de reference (normalise) ---------------
        profil = self.demi_section_gen(
            self.DEADRISE, self.FLARE, self.F_BOUCHAIN, self.W_BOUCHAIN,
            demi_B, self.CREUX,
            self.M_FOND, self.M_BOUCHAIN, self.M_BORDE)

        # --- 2. plan de forme : demi-largeur le long de x ------------------
        avant = self.bezier_rationnelle(
            np.array([0.0, self.B_ETRAVE * demi_B]),
            np.array([self.REMPL_AV * x_maitre, demi_B]),
            np.array([x_maitre, demi_B]), 1.0, 200)

        arriere = self.bezier_rationnelle(
            np.array([x_maitre, demi_B]),
            np.array([x_maitre + self.REMPL_AR * (self.L_COQUE - x_maitre), demi_B]),
            np.array([self.L_COQUE, self.B_TABLEAU * demi_B]), 1.0, 200)

        plan = np.vstack([avant, arriere[1:]])

        X = np.linspace(0.0, self.L_COQUE, self.N_SECTIONS)
        B_loc = np.interp(X, plan[:, 0], plan[:, 1])

        # --- 3. rocker : hauteur de la quille le long de x -----------------
        z_quille = np.where(
            X < x_maitre,
            self.ROCKER_AV * ((x_maitre - X) / x_maitre) ** 2,
            self.ROCKER_AR * ((X - x_maitre) / (self.L_COQUE - x_maitre)) ** 2)
        h_loc = self.CREUX - z_quille

        # --- 4. anneaux fermes (section complete + pont) -------------------
        def anneau(i):
            """Contour ferme de la section n.i, dans l'ordre, sans doublon."""
            y = profil[:, 0] * B_loc[i]
            z = z_quille[i] + profil[:, 1] * h_loc[i]
            droite = np.column_stack([y, z])
            gauche = droite[::-1][:-1] * np.array([-1, 1])

            y_pont = np.linspace(B_loc[i], -B_loc[i], self.M_PONT)[1:-1]
            fleche = self.HAUTEUR_BOMBE * (1.0 - (y_pont / max(B_loc[i], 1e-9)) ** 2)
            pont = np.column_stack([y_pont, self.CREUX + fleche])
            yz = np.vstack([droite, pont, gauche])
            return np.column_stack([np.full(len(yz), X[i]), yz])

        anneaux = np.array([anneau(i) for i in range(self.N_SECTIONS)])
        m_anneau = anneaux.shape[1]

        vertices = anneaux.reshape(-1, 3)

        def idx(i, j):
            return i * m_anneau + (j % m_anneau)

        # --- 5. peau : deux triangles par quadrilatere ---------------------
        # dans anneau(i), le contour est [droite | pont | gauche] :
        # les aretes du pont sont celles dont j va du dernier point de
        # "droite" (livet tribord) au premier point de "gauche" (livet babord).
        n_droite = len(profil)
        n_pont = self.M_PONT - 2
        j_pont = set(range(n_droite - 1, n_droite + n_pont))

        faces = []
        face_pont = []
        for i in range(self.N_SECTIONS - 1):
            for j in range(m_anneau):
                a, b = idx(i, j), idx(i, j + 1)
                c, d = idx(i + 1, j), idx(i + 1, j + 1)
                faces.append([a, b, d])
                faces.append([a, d, c])
                est_pont = j in j_pont
                face_pont.append(est_pont)
                face_pont.append(est_pont)

        # --- 6. bouchons avant / arriere (eventail) ------------------------
        vertices_extra = []

        def bouchon(i_station, sens):
            centre = anneaux[i_station].mean(axis=0)
            centre[1] = 0.0    # la section est symetrique : centre sur l'axe
            i_centre = len(vertices) + len(vertices_extra)
            tris = []
            for j in range(m_anneau):
                a, b = idx(i_station, j), idx(i_station, j + 1)
                tris.append([i_centre, a, b] if sens > 0 else [i_centre, b, a])
            return [centre], tris

        v, t = bouchon(0, +1)
        vertices_extra += v
        faces += t
        face_pont += [False] * len(t)      # bouchons = pas du pont
        v, t = bouchon(self.N_SECTIONS - 1, -1)
        vertices_extra += v
        faces += t
        face_pont += [False] * len(t)

        vertices = np.vstack([vertices, np.array(vertices_extra)])
        faces = np.array(faces)
        face_pont = np.array(face_pont, dtype=bool)

        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=True)
        mesh.fix_normals()
        mesh.face_pont = face_pont         # masque bool. : True = face du pont

        self.mesh = mesh
        self.face_pont = face_pont
        return mesh

    # --- exploitation du maillage genere ---------------------------------

    def get_area(self, mask=None):
        """Aire totale (m^2) du maillage, ou d'un sous-ensemble de faces si
        `mask` est fourni (ex : coque.get_area(coque.face_pont))."""
        if self.mesh is None:
            raise RuntimeError("generate() doit etre appele avant get_area()")
        faces = self.mesh.faces if mask is None else self.mesh.faces[mask]
        return trimesh.triangles.area(self.mesh.vertices[faces]).sum()

    def controle_maillage(self):
        if self.mesh is None:
            raise RuntimeError("generate() doit etre appele avant controle_maillage()")
        mesh = self.mesh
        print("\n=== Controle du maillage =============================")
        print(f"Sommets                : {len(mesh.vertices)}")
        print(f"Triangles              : {len(mesh.faces)}")
        print(f"is_watertight (etanche): {mesh.is_watertight}")
        print(f"is_volume (volume sur) : {mesh.is_volume}")
        n_ouvertes = len(trimesh.grouping.group_rows(mesh.edges_sorted,
                                                     require_count=1))
        print(f"Aretes ouvertes        : {n_ouvertes}")
        print(f"Euler number (2=sphere): {mesh.euler_number}")
        if not mesh.is_watertight:
            print("!! Maillage NON etanche : les mesures n'ont aucun sens.")


if __name__ == "__main__":
    start = time.time()
    coque = CoqueMesh()
    coque.generate()
    print(f"généré en : {time.time() - start}")
    coque.controle_maillage()
