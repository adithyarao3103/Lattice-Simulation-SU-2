# Lattice Simulation of pure SU(2)

- Code written by me to simulate SU(2) pure gauge theory in 1+1 dimensions.
- This code is inefficient and prone to memory leaks. It is only to display an understanding of the lattice simulations
- Especially the definition of a lattice. Here I have implemented a more natural 4-dimensional array that mimics a physical lattice, while the most efficient way would be to encode the lattice serially.
- I have not implemented any observables so far. 
- The only part that has been implemented is the lattice action and the updating of the lattice according to the Markov Chain method.
- Have not used any standard libs in C++ for complex numbers and matrices. Implemented the necessary functionalities from scratch
