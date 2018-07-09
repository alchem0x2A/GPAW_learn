from ase.build import molecule
from ase.parallel import paropen, parprint
from gpaw import GPAW, PW

'''
Calculate the potential energy of H2O
by different grid spacing and BZ k-grid
density
'''

def get_energy(k):
    a = 8
    e_cut = 750
    mol = molecule("H2O")
    mol.set_cell((a, a, a))
    mol.set_pbc((True, True, True))
    mol.center()
    calc = GPAW(mode=PW(e_cut),
                kpts=(k, k, k))
    mol.set_calculator(calc)
    return mol.get_potential_energy()

with paropen("E_convergence_PW_k.txt", "w") as f:
    for k in range(1, 15, 2):
        e = get_energy(k)
        parprint("{0},{1:.4f}".format(k, e),
                 file=f)
