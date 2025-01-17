import os, ROOT
from math import hypot, pi

# ========= UTILITIES FROM FRANK =======================


def deltaPhi(phi1, phi2):
    # Catch if being called with two objects
    if type(phi1) != float and type(phi1) != int:
        phi1 = phi1.phi
    if type(phi2) != float and type(phi2) != int:
        phi2 = phi2.phi
    # Otherwise
    dphi = (phi1 - phi2)
    while dphi > pi:
        dphi -= 2 * pi
    while dphi < -pi:
        dphi += 2 * pi
    return dphi


def deltaR(eta1, phi1, eta2=None, phi2=None):
    # catch if called with objects
    if eta2 == None:
        return deltaR(eta1.eta, eta1.phi, phi1.eta, phi1.phi)
    # otherwise
    return hypot(eta1 - eta2, deltaPhi(phi1, phi2))


def closest(obj, collection, presel=lambda x, y: True):
    ret = None
    drMin = 999
    for x in collection:
        if not presel(obj, x):
            continue
        dr = deltaR(obj, x)
        if dr < drMin:
            ret = x
            drMin = dr
    return (ret, drMin)

# ================= MY UTILITIES ======================

#Returns true if phiTest is in the small angle between phiA and phiB
def isBetween(phiA, phiB, phiTest):
    #Handle object inputs if provided
    if type(phiA) != float and type(phiA) != int:
        phiA = phiA.phi()
    if type(phiB) != float and type(phiB) != int:
        phiB = phiB.phi()
    if type(phiTest) != float and type(phiTest) != int:
        phiTest = phiTest.phi()

    #Reduce all angles > 2pi to be < 2pi
    while phiA > (2 * pi):
        phiA -= (2 * pi)
    while phiB > (2 * pi):
        phiB -= (2 * pi)
    while phiTest > (2 * pi):
        phiTest -= (2 * pi)
        
    #Make sure all angles are measured in the same direction
    if phiA < 0:
        phiA += (2 * pi)
    if phiB < 0:
        phiB += (2 * pi)
    if phiTest < 0:
        phiTest += (2 * pi)

    
    if abs(phiA - phiB) == pi: #If A & B are back-to-back, any angle is "between" and "not between"
        return False
    elif abs(phiA - phiB) < pi: 
        return ( min(phiA, phiB) < phiTest and phiTest < max(phiA, phiB) )
    else:
        return ( min(phiA, phiB) > phiTest or  phiTest > max(phiA, phiB) )


    
