# Main program

- [ ] Input reader

# SCF program

- [x] Initial guess by diagonalizing H = T + V
- [x] Construct Fock matrix 
- [x] SCF cycles
- [x] SCF energy
- [x] UHF
- [ ] ROHF

# Library
- [ ] HDF5 support instead of TOML

## Base
- [x] Build basis set framework
- [x] Initialize basis on atoms
- [ ] Define orbital struct
- [ ] `.mwfn` support (?)

## Integrals
- [x] Nuclear-nuclear repulsion energy
- [x] Overlap integral (S)
- [x] kinetic energy integral (T)
- [x] nuclear-electron attraction integral (V)
- [x] ERIs
- [ ] dipole

## Math
- [x] Custom Matrix, Vector, MatrixContainer structs
- [x] Use BLAS (Accelerate)
- [ ] BLAS Level 1
- [ ] BLAS Level 2
- [ ] BLAS Level 3

# Properties
- [ ] Dipole

# Theory
- [ ] Write down equations used in documentation

# GUI?

# Symmetry

# Time-Dependent HF
