import sys
from ase import Atoms
from ase.io import Trajectory   # Use Trajectory module

def get_angle(mol):
    # H1 -- 0; H2 -- 1; O -- 2
    return mol.get_angle(a1=0, a2=2, a3=1)

if __name__ == "__main__":
    assert len(sys.argv) == 2
    name = sys.argv[1]
    assert ".traj" in name
    traj = Trajectory(name)
    angles = [get_angle(mol) for mol in traj]
    import matplotlib.pyplot as plt
    plt.plot(range(len(traj)), angles)
    plt.xlabel("Steps")
    plt.ylabel("Angle [$^{\\circ}$]")
    plt.show()
    
