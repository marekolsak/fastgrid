FastGrid is a fork of [AutoGrid](http://autodock.scripps.edu/wiki/AutoGrid), which computes grid maps of interaction energies for various atom types around a macromolecule. These grid maps are then used by [AutoDock](http://autodock.scripps.edu) for docking calculations, which determine the total interaction energy for a ligand with a macromolecule.

FastGrid is meant to be a replacement of [AutoGrid](http://autodock.scripps.edu/wiki/AutoGrid), having the same behavior as [AutoGrid](http://autodock.scripps.edu/wiki/AutoGrid) 4.2.1, and being more than 100 times faster depending on a molecule and the size of the grid map being computed.

The electrostatic potential computation run on the GPU using NVIDIA CUDA. The rest of the computations run on the CPU using more sophisticated data structures and algorithms with lower computational complexity than those of [AutoGrid](http://autodock.scripps.edu/wiki/AutoGrid). Moreover, the CPU implementation is parallelized to fully exploit computational power of machines that are equipped with multiple CPU cores.

The license of this software is available [here](http://code.google.com/p/fastgrid/source/browse/fastgrid_master/fastgrid/COPYING).

The conference paper about this project can be downloaded [here](http://code.google.com/p/fastgrid/downloads/detail?name=memics09.pdf&can=2&q=).

### References ###
  1. [Olšák, Marek - Filipovič, Jiří - Prokop, Martin. FastGrid -- The Accelerated AutoGrid Potential Maps Generation for Molecular Docking.](http://code.google.com/p/fastgrid/downloads/detail?name=memics09.pdf&can=2&q=) In MEMICS 2009: Fifth Doctoral Workshop on Mathematical and Engineering Methods in Computer Science. Brno, Czech Republic : Masaryk University and Technical University of Brno, 2009. od s. 160-167, 8 s. ISBN 978-80-87342-04-6.
  1. http://autodock.scripps.edu/wiki/AutoGrid