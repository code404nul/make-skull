from mesh import CoqueMesh
from params import *

def get_comsuption(coque_mesh):
    return 616  * (coque_mesh.L_COQUE-0.3) * coque_mesh.B_MAX # Moy Wh /j
    # 4 * 0.22 * 0.7 * 1000 https://epitortue.archibarbu.art/solutions/s7-bilan-energie/



    