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


    
#----------------------------------------------------------------------------------------------------------

def getSFFile(year, pog, typ=""):
    filename = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/"
    filename += pog + "/"

    if year == "2016":
        filename += "2016preVFP_UL/"
    elif year == "2016post":
        filename += "2016postVFP_UL/"
    elif year == "2017":
        filename += "2017_UL/"
    elif year == "2018":
        filename += "2018_UL/"
    elif year == "2022":
        filename += "2022_Summer22/"
    elif year == "2022post":
        filename += "2022_Summer22EE/"
    elif year == "2023":
        filename += "2023_Summer23/"
    elif year == "2023post":
        filename += "2023_Summer23BPix/"
    else:
        print("ERROR: Unrecognized YEAR was provided to getSFFile")
        exit(2)
    
    if pog == "EGM":
        filename += "electron.json.gz"
    elif pog == "JME":
        if typ == "ID":
            if year in ["2016", "2016post", "2017", "2018"]:
                filename += "jmar.json.gz"
            else:
                filename += "jetid.json.gz"
        elif typ == "MET":
            filename += "met.json.gz"
        elif typ == "JERC":
            filename += "fatJet_jerc.json.gz" 
        elif typ == "VETO":
            filename += "jetvetomaps.json.gz"
        else:
            print("ERROR: No type for JME SF was provided")
            exit(2)
    elif pog == "MUO":
        filename += "muon_Z.json.gz"
    elif pog == "LUM":
        filename += "puWeights.json.gz"
    elif pog == "TAU":
        if year in ["2016", "2016post", "2017", "2018"]:
            filename += "tau.json.gz"
        else:
            if year == "2022":
                filename += "tau_DeepTau2018v2p5_2022_preEE.json.gz"
            elif year == "2022post":
                filename += "tau_DeepTau2018v2p5_2022_postEE.json.gz"
            elif year == "2023":
                filename += 'tau_DeepTau2018v2p5_2023_preBPix.json.gz'
            elif year == "2023post":
                filename += 'tau_DeepTau2018v2p5_2023_postBPix.json.gz'
    else:
        print("ERROR: Unrecognized POG was provided to getSFFile")
        exit(2)

    return filename

#--------------------------------------------------------------------------------------

yearToEGMSfYr = {
    "2016" : "2016preVFP",
    "2016post" : "2016postVFP",
    "2017" : "2017",
    "2018" : "2018",
    "2022" : "2022Re-recoBCD",
    "2022post" : "2022Re-recoE+PromptFG",
    "2023" : "2023PromptC",
    "2023post" : "2023PromptD"
}

yearToJetVeto = {
    "2016" : "Summer19UL16_V1",
    "2016post" : "Summer19UL16_V1",
    "2017" : "Summer19UL17_V1",
    "2018" : "Summer19UL18_V1",
    "2022" : "Summer22_23Sep2023_RunCD_V1",
    "2022post" : "Summer22EE_23Sep2023_RunEFG_V1",
    "2023" : "Summer23Prompt23_RunC_V1",
    "2023post" : "Summer23BPixPrompt23_RunD_V1"
}
