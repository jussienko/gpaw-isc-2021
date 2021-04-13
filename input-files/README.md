# Input cases

- [copper.py](copper.py): the main performance benchmark.
- [copper-profile.py](copper-profile.py): the main performance benchmark with Python cProfile 
performance analysis. Produces a profile file for each MPI rank that can be investigated with 
Python [pstats](https://docs.python.org/3/library/profile.html) module. 
- [nanoribbon.py](nanoribbon.py): input for visualization, produces a `.cube` file containing 
electron localization function.
- [si-divacancy.py](si-divacancy.py): **MPI_Alltoallv** communication benchmark.


