from ase.atoms import Atoms
from ase.parallel import paropen
from ase.parallel import parprint
from gpaw import GPAW, PW
from ase.build import molecule

def calc_energy(formula, a=8.0):
    mol = molecule(formula)
    if len(mol) > 1:
        hund = False             # an atom
    else:
        hund = True            # molecule
    mol.set_cell((a, a, a))
    mol.center()
    mol.set_pbc((True, True, True))
    calc = GPAW(mode=PW(),
                xc="PBE",
                kpts={"density": 2.5},
                hund=hund,
                txt="{0}.out".format(mol.get_chemical_formula()))
    mol.set_calculator(calc)
    E = mol.get_potential_energy()
    return E

with paropen("atomization_H2O_molecule.txt", "w") as f:
    E1 = calc_energy("H")
    parprint("H: {:.4f} eV".format(E1), file=f)
    E2 = calc_energy("O")
    parprint("O: {:.4f} eV".format(E2), file=f)
    E3 = calc_energy("H2O")
    parprint("H2O: {:.4f} eV".format(E3), file=f)
    E = 2 * E1 + E2 - E3
    parprint("Atomization energy H2O: {:.4f} eV".format(E), file=f)   
