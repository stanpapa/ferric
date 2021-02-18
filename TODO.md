Necessary:
1. Build basis set framework
2. Initialize basis on atoms
3. Define orbital struct
4. Overlap (S), kinetic energy (T) and nuclear-electron attraction (V) integrals
5. Initial guess by diagonalizing H = T + V
6. ERIs
7. Construct Fock matrix
8. SCF cycles
9. Find a good name for the program

Extra:
- [ ] Dipole
- [ ] Custom matrix implementation

Rust + cargo related:
- Figure out how to properly make a library inside this crate
- Figure out how to include multiple basis sets which are stored in their own file, without having to type them all out in mod.rs.
- Figure out how to access variables and functions from other modules without super long paths, e.g. crate::gto_integrals::nuclear_repulsion::nuclear_repulsion(&atoms)

Finished:
- Nuclear-nuclear repulsion energy
