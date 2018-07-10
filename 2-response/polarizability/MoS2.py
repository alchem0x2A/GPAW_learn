from ase.io import Trajectory
from ase.parallel import paropen, parprint
from gpaw import GPAW, PW, FermiDirac, PoissonSolver
from gpaw.response.df import DielectricFunction

traj = Trajectory(filename="MoS2.traj")
mol = traj[-1]
mol.set_pbc((True, True, False))

calc = GPAW(mode=PW(800),
            kpts=dict(density=4,
                      gamma=True),
            occupations=FermiDirac(0.001),
            xc="PBE",
            basis="dzp",
            poissonsolver=dict(dipolelayer="xy"))

mol.set_calculator(calc)
# Get ground state
mol.get_potential_energy()
calc.diagonalize_full_hamiltonian(nbands=30)  # Full diagonalization
calc.write("MoS2_gs.gpw", mode="all")

# Dielectric Response
df = DielectricFunction(calc="MoS2_gs.gpw",
                        eta=0.05,
                        domega0=0.02,
                        truncation="2D",  # truncation
                        ecut=50)

alpha0x, alphax = df.get_polarizability(q_c=[0, 0, 0],
                                        direction='x',
                                        pbc=[True, True, False],
                                        filename="MoS2_alpha_x.csv")

alpha0z, alphaz = df.get_polarizability(q_c=[0, 0, 0],
                                        direction='z',
                                        pbc=[True, True, False],
                                        filename="MoS2_alpha_z.csv")




