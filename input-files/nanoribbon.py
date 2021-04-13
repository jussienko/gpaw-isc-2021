###
### GPAW visualization : Graphene nanoribbon with gold adatom
###
### Writes electron localization function into file elf_ribbon.cube
###


from ase import Atoms
from ase.build import graphene_nanoribbon, add_adsorbate
from ase.io import write
from gpaw import GPAW
from gpaw.occupations import FermiDirac
from gpaw.eigensolvers.rmmdiis import RMMDIIS
from gpaw.mpi import size, rank

from gpaw.elf import ELF

# Graphene nanoribbon
ribbon = graphene_nanoribbon(5, 6, type='armchair', saturated=True,
                            vacuum=3.5)

# Gold adsorbate 
pos = (ribbon[35].position + ribbon[60].position) / 2.0
pos[1] += 2.2
adsorbate = Atoms('Au', (pos,))
ribbon += adsorbate
ribbon.center(axis=1, vacuum=3.5)

txt = 'output_p_{}.txt'.format(size)
ribbon.calc = GPAW(h=0.22,
                   xc='PBE', 
                   txt=txt,
                   occupations=FermiDirac(0.2),
                   eigensolver=RMMDIIS(),
                  )

ribbon.get_potential_energy()

# Calculate electron localization function
elf = ELF(ribbon.calc)
elf.update()
elf_g = elf.get_electronic_localization_function(gridrefinement=2)
if rank == 0:
    write('elf_ribbon.cube', ribbon, data=elf_g)
