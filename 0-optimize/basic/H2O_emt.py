from ase import Atoms
from ase.calculators.emt import EMT  # EMT calculator
from ase.optimize import QuasiNewton

system = Atoms('H2O', positions=[[0.0, 0.5, 0.5],
                                 [0.0, 0.5, -0.5],
                                 [0.0, 0.0, 0.0]])  # In \AA
calc = EMT()

system.set_calculator(calc)

opt = QuasiNewton(system, trajectory='H2O.emt.traj')

opt.run(fmax=0.001)
