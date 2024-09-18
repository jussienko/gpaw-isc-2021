###
### GPAW performance benchmark: Copper Filament
###

from ase.lattice.cubic import FaceCenteredCubic
from ase.units import Ha
from gpaw.mpi import size, rank
from gpaw import GPAW, Mixer, ConvergenceError
from gpaw.occupations import FermiDirac
from gpaw.eigensolvers.rmmdiis import RMMDIIS
import pytest
#from gpaw.test import equal

# no. of replicates in each dimension (increase to scale up the system)
x = 2
y = 2
z = 3
# setup the system
atoms = FaceCenteredCubic(directions=[[1,-1,0], [1,1,-2], [1,1,1]],
        size=(x,y,z), symbol='Cu', pbc=(0,0,1))
atoms.center(vacuum=6.0, axis=0)
atoms.center(vacuum=6.0, axis=1)

# Simulation parameters
h = 0.22
kpts = (1,1,8)
txt = 'output_M_%i.txt' % size
maxiter = 15

# output benchmark parameters
if rank == 0:
    print("#"*60)
    print("GPAW benchmark: Copper Filament")
    print("  dimensions: x=%d, y=%d, z=%d" % (x, y, z))
    print("  grid spacing: h=%f" % h)
    print("  Brillouin-zone sampling: kpts=" + str(kpts))
    print("  MPI tasks: %d" % size)
    print("#"*60)
    print("")

# setup parameters
args = {'mode' : 'fd',
        'h': h,
        'nbands': -20,
        'occupations': FermiDirac(0.2),
        'kpts': kpts,
        'xc': 'PBE',
        'mixer': Mixer(0.1, 5, 100),
        'eigensolver': RMMDIIS(),
        'maxiter': maxiter,
        'txt': txt}

calc = GPAW(**args)
atoms.calc = calc

# execute the run
try:
    atoms.get_potential_energy()
except ConvergenceError:
    pass

e0 = Ha * calc.hamiltonian.e_total_free
if rank == 0:
    print("Free energy: " + str(e0))

# Check the result
assert e0 == pytest.approx(-226.93897028587003, abs=1e-4)

