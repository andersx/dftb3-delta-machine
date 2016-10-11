#!/usr/bin/env python2
import numpy as np
import cPickle

from numpy.linalg import norm
from numpy import exp, sqrt, sum, square

NUCLEAR_CHARGE = dict()
NUCLEAR_CHARGE["H"] = 1.0
NUCLEAR_CHARGE["C"] = 6.0
NUCLEAR_CHARGE["N"] = 7.0
NUCLEAR_CHARGE["O"] = 8.0
NUCLEAR_CHARGE["S"] = 16.0

COULOMB_MATRIX_SIZE = 23

HOF_DFTB3 = dict()
HOF_DFTB3["H"] = -172.3145
HOF_DFTB3["C"] = -906.4342
HOF_DFTB3["N"] = -1327.2991
HOF_DFTB3["O"] = -1936.6161
HOF_DFTB3["S"] = -1453.3907

class Molecule:

    def __init__(self):
        self.natoms = -1
        self.energy = float("nan")
        self.molid = -1
        self.dftb3_energy = float("nan")
        self.dftb3_hof = float("nan")

        self.atomtypes = []
        self.coordinates = []

    def generate_coulomb_matrix(self):
        self.coulomb_matrix = generate_coulomb_matrix(self.atomtypes, self.coordinates)

def get_lines(filename):

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    return lines


def parse_dft3_energy(molid):

    filename = "logfiles/" + str(molid) + ".log"
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    energy = float("nan")
    for line in lines:
        if "Total Energy" in line:
            tokens = line.split()
            energy = float(tokens[2]) * 627.51

    return energy


def parse_molecules(filename):

    lines = get_lines(filename)

    mols = []

    mol = Molecule()

    for line in lines:

        tokens = line.split()

        if len(tokens) == 1:

            if mol.natoms > 0:
                mols.append(mol)

            mol = Molecule()
            mol.natoms = int(tokens[0])

        if len(tokens) == 2:
            mol.molid = int(tokens[0])
            mol.energy = float(tokens[1])
            mol.dftb3_energy = parse_dft3_energy(mol.molid)


        if len(tokens) == 7:

            mol.atomtypes.append(tokens[0])
            x = float(tokens[4])
            y = float(tokens[5])
            z = float(tokens[6])

            mol.coordinates.append(np.array([x, y, z]))

            mol.dftb3_hof = 0.0
            mol.dftb3_hof += mol.dftb3_energy

            for atom in ["H", "C", "N", "O", "S"]:

                n = mol.atomtypes.count(atom)
                mol.dftb3_hof -= n * HOF_DFTB3[atom] 

    # for mol in mols:
    #     print mol.molid, mol.energy, mol.dftb3_hof

    return mols


def generate_coulomb_matrix(atomtypes, coordinates):

    # Generate row norms for sorting

    row_norms = []

    for i, atomtype_i in enumerate(atomtypes):

        row_norm = 0.0

        for j, atomtype_j in enumerate(atomtypes):

            if i == j:
                row_norm += 0.5 * NUCLEAR_CHARGE[atomtype_i] ** 2.4

            else:
                row_norm += NUCLEAR_CHARGE[atomtype_i] * NUCLEAR_CHARGE[atomtype_j] \
                            / np.linalg.norm(coordinates[i] - coordinates[j])

        row_norms.append((row_norm, i))

    # Sort by row norms

    row_norms.sort(reverse=True)

    sorted_atomtypes = []
    sorted_coordinates = []

    for row_norm in row_norms:

        i = row_norm[1]

        sorted_atomtypes.append(atomtypes[i])
        sorted_coordinates.append(coordinates[i])

    # Fill out 

    Mij = np.zeros((COULOMB_MATRIX_SIZE, COULOMB_MATRIX_SIZE))
    for i, atomtype_i in enumerate(sorted_atomtypes):
        for j, atomtype_j in enumerate(sorted_atomtypes):
            if i == j:
                Mij[i, j] = 0.5 * NUCLEAR_CHARGE[atomtype_i] ** 2.4

            elif j > i:
                continue

            else:
                Mij[i, j] = NUCLEAR_CHARGE[atomtype_i] * NUCLEAR_CHARGE[atomtype_j] \
                            / np.linalg.norm(sorted_coordinates[i] - sorted_coordinates[j])

    return Mij.flatten(COULOMB_MATRIX_SIZE**2)


def load_pickle(filename):
    f = open(filename,"rb")
    p = cPickle.load(f)
    f.close()
    return(p)

# def parse_dftb_energies(filename):
