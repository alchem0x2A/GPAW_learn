from ase.build import molecule
from ase.parallel import paropen, parprint
from gpaw import GPAW

'''
Calculate the potential energy of H2O
by different grid spacing and BZ k-grid
density
'''

def get_energy(gpt):
    a = 8.0
    h = a / gpt
    mol = molecule("H2O")
    mol.set_cell((a, a, a))
    mol.set_pbc((True, True, True))
    mol.center()
    calc = GPAW(h=h)
    mol.set_calculator(calc)
    return mol.get_potential_energy()

with paropen("E_convergence.txt", "w") as f:
    for gpt in range(12, 64, 4):
        e = get_energy(gpt)
        parprint("{0},{1:.4f}".format(gpt, e),
                 file=f)
