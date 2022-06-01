# SNECP2_PDE_OpenMPI_GSL
Solution of task C3 of the SNE benchmark CP2. Definition: https://www.sne-journal.org/fileadmin/user_upload_sne/benchmarks/CP/CP2-definition.pdf

* Build: make all

* parallel implementation:
  * usage: mpirun -np (number of processes) ./para
* sequential implemetation:
  * usage: ./seq

* plot:
  * plot_ut3L4_L2(file_3L4, file_L2) -> plot the solution of the PDE, excitation over time at x = 1/2 L and x = 2/4 L
  * plot_uxt5_t8(file_t5, file_t8) -> plot the solution of the PDE, excitation over space at t = 5 and t = 8
