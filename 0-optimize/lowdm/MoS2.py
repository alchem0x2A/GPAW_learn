from ase.build import mx2
from ase.constraints import UnitCellFilter
import sys
import os.path
# solve path issue
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
print(sys.path)
from optimize import QuasiNewton
from gpaw import GPAW, PW, FermiDirac, PoissonSolver
from ase.parallel import parprint, paropen

mol = mx2("MoS2", kind="2H", a=3, vacuum=7.5)
mol.set_pbc((True, True, False))
mol.center(axis=2)
calc = GPAW(mode=PW(800),
            kpts=dict(density=4,
                      gamma=True),
            occupations=FermiDirac(0.01),
            xc="PBE",
            basis="dzp",
            poissonsolver=dict(dipolelayer="xy"),
            # txt="MoS2_relax.txt"
)
mol.set_calculator(calc)

if __name__ == "__main__":
    opt = QuasiNewton(mol, trajectory="MoS2.traj")
    opt.run(fmax=0.01, smax=2e-4, smask=[1, 1, 0, 0, 0, 0])
    calc.write("MoS2.gpw")
    parprint(mol.cell)
