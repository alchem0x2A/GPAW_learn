from ase import Atom, Atoms
from ase.io import Trajectory
from gpaw import GPAW, PW, FermiDirac
from ase.optimize import QuasiNewton

a = 8
c = a / 2
H = Atoms("H", positions=[(c, c, c)],
          magmoms=[0], cell=(a, a, a),
          pbc=(True, True, True))

O = Atoms("O", positions=[(c, c, c)],
          magmoms=[0], cell=(a, a, a),
          pbc=(True, True, True))

H2O = Atoms("H2O", positions=[(c, c+0.5, c-0.5),
                              (c, c+0.5, c+0.5),
                              (c, c, c)],
            cell=(a, a, a),
            pbc=(True, True, True))

calc = GPAW(mode=PW(),
            xc="PBE",
            hund=True,
            eigensolver="rmm-diis")

H.set_calculator(calc)
H.center()
E1 = H.get_potential_energy()
calc.write("H.gpw")

O.set_calculator(calc)
O.center()
E2 = O.get_potential_energy()
calc.write("O.gpw")

calc.set(hund=False)
H2O.set_calculator(calc)
H2O.center()
opt = QuasiNewton(H2O, trajectory="H2O.traj")
opt.run(fmax=0.001)
E3 = H2O.get_potential_energy()

fd = open('atomization_H2O.txt', 'w')
print('  hydrogen atom energy:     %5.2f eV' % E1, file=fd)
print('  oxygen atom energy: %5.2f eV' % E2, file=fd)
print('  hydrogen molecule energy: %5.2f eV' % E3, file=fd)
print('  atomization energy:       %5.2f eV' % (2 * E1 + E2 - E3), file=fd)
fd.close()
