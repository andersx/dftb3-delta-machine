#!/usr/bin/env python2
from numpy.linalg import norm
from numpy import exp
import numpy as np
from fkernels import kgaussian_kernel as klg
from fkernels import klaplace_kernel as kkl

def klaplace_kernel(a, b, sigma_in):


    na =  a.shape[0]
    nb =  b.shape[0]

    k = np.empty((na, nb), order='F')
    a = np.asfortranarray(a).T
    b = np.asfortranarray(b).T

    kkl(a, na, b, nb, k, sigma_in)

    return np.ascontiguousarray(k)


def kgaussian_kernel(a, b, sigma_in):

    na =  a.shape[0]
    nb =  b.shape[0]

    k = np.empty((na, nb), order='F')
    a = np.asfortranarray(a).T
    b = np.asfortranarray(b).T

    kkg(a, na, b, nb, k, sigma_in)

    return np.ascontiguousarray(k)


gauss_factor = -0.5 / (724.0 * 724.0) 
laplace_factor = -1.0 / (10.0**(3.6))

def gaussian_kernel(mol_a, mol_b):
    return exp(gauss_factor * norm(mol_a.coulomb_matrix - mol_b.coulomb_matrix)) 

def laplace_kernel(mol_a, mol_b):
    return exp(laplace_factor * norm(mol_a.coulomb_matrix - mol_b.coulomb_matrix, ord=1)) 
