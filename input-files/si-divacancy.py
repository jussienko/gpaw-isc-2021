###
### GPAW benchmark: Si divacancy
###
### Includes more intensive MPI_alltoallv communication 
### due to plane wave basis set
###


from ase.build import bulk
from gpaw import GPAW, PW, ConvergenceError
from gpaw.mpi import size
from gpaw.eigensolvers.rmmdiis import RMMDIIS

# nx, ny, nz supercell with divacancy
nx, ny, nz = 3, 3, 3
atoms = bulk('Si', cubic=True)
atoms *= (nx, ny, nz)
atoms.pop(0)
atoms.pop(0)

outfile = "out_{}x{}x{}_n{}.txt".format(nx, ny, nz, size)

atoms.calc = GPAW(mode=PW(180),
                  xc='PBE',
                  kpts=(2,2,2),
                  nbands=-30,
                  eigensolver=RMMDIIS(),
                  txt=outfile,
                  maxiter=10,
                 )

try:
    atoms.get_potential_energy()
except ConvergenceError:
    pass

