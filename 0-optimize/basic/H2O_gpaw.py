from ase import Atoms
from gpaw import GPAW
from ase.optimize import QuasiNewton

system = Atoms('H2O', positions=[[0.0, 0.5, 0.5],
                                 [0.0, 0.5, -0.5],
                                 [0.0, 0.0, 0.0]],
               cell=[[5.0, 0.0, 0.0],
                     [0.0, 5.0, 0.0],
                     [0.0, 0.0, 5.0]],
               pbc=(True, True, True))  # In \AA

system.center()                 # Center the structure

calc = GPAW(xc="LDA")

system.set_calculator(calc)

opt = QuasiNewton(system, trajectory='H2O.gpaw.traj')

opt.run(fmax=0.001)
