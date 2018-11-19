#!/usr/bin/env python

"""
Functions used in recreation of calculations from De Domenico & Insolia 2012. 

@author Francesca Capel
@date October 2018
"""

from energy_loss import DH, get_source_threshold_energy

zmax = 300 / DH

def omega_gzk(z, Earr, alpha):
    """
    See Equation 10. (E in EeV)
    """
    D = z * DH
    Ei = get_source_threshold_energy(Earr, D)[0]
    return (Ei / Earr)**(1 - alpha)

def integrand(z, Earr, alpha):
    """
    Integrand for Omega_gzk calculation.
    """
    D = z * DH
    Ei = []
    for d in D:
        Ei.append(get_source_threshold_energy(Earr, d)[0])
    return [ei**(1 - alpha) for ei in Ei]

def Omega_gzk(z, Earr, alpha, norm):
    """
    See Equation 11.
    """
    numerator, err = integrate.fixed_quad(integrand, z, zmax, args = (Earr, alpha))
    denominator = norm
    return numerator / denominator
