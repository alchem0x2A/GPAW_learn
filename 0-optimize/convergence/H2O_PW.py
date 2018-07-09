from ase.build import molecule
from ase.parallel import paropen, parprint
from gpaw import GPAW, PW

'''
Calculate the potential energy of H2O
by different grid spacing and BZ k-grid
density
'''

def get_energy(e_cut):
    a = 8
    mol = molecule("H2O")
    mol.set_cell((a, a, a))
    mol.set_pbc((True, True, True))
    mol.center()
    calc = GPAW(mode=PW(e_cut))
    mol.set_calculator(calc)
    return mol.get_potential_energy()

with paropen("E_convergence_PW.txt", "w") as f:
    for e_cut in range(300, 800, 50):
        e = get_energy(e_cut)
        parprint("{0},{1:.4f}".format(e_cut, e),
                 file=f)
