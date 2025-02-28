#Boosted, Z-rest frame reclusterinf of PF candidates corresponding to Z->FatJet events selected by ZProducer
#Requires the information added by the bvnano-prod package: https://github.com/cms-btv-pog/btvnano-prod/tree/NanoAODv12_22Sep2023
# 

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from ROOT import TLorentzVector, TVector3

# ----------------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------------- #

#Returns distance dij from https://arxiv.org/pdf/0802.1189
#objs can either be TLorentzVectors or ReClJets
#R is the radius, p the power (not including the "2") which the momentum is raised two
def distance(obj1, obj2, p, R):
    p = 2 * p
    return min(obj1.Pt()**p, obj2.Pt()**p) * (obj1.DeltaR(obj2)**2) / (R**2)

#Recluster the PF candidate four-vectors pfCands after boosting by the boost TLorentzVector
#p=-1 => anti-kt, p=0 => cambridge/aachen, p=1 =>inclusive kt
#R is the jet cone radius i.e. 0.4 = AK4, 0.8 = AK8, etc.
#Returns a list of four vectors representing the reclustered jets
def recluster(pfCands, p, R):

    reClJets = []
    
    while len(pfCands) > 0:
        minDist = 9999999
        minPairIdxs = None
        minBDist = 9999999
        minBDistIdx = -1

        for pfcN1, pfc1 in enumerate(pfCands):
            
            #Find them minimum distance between pfc1 and all other PFCs
            for pfcN2 in range(pfcN1+1, len(pfCands)):
                pfc2 = pfCands[pfcN2]
                dist = distance(pfc1, pfc2, p, R)
                if dist < minDist:
                    minDist = dist
                    minPairIdxs = (pfcN1, pfcN2)
            
            #Find the minimum "beam distance"
            bDist = pfc1.Pt()**(2*p)
            if bDist < minBDist:
                minBDist = bDist
                minBDistIdx = pfcN1
            
        #If two PFCs are the closest pairing, combine them into the first and delete the second from the list
        if minDist < minBDist and minPairIdxs is not None:
            pfCands[minPairIdxs[0]] += pfCands[minPairIdxs[1]] 
            del pfCands[minPairIdxs[1]]
        else: #If beam distance is minimum, call that a reclustered jet
            reClJets.append(pfCands[minBDistIdx])
            del pfCands[minBDistIdx]

    return reClJets


# ----------------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------------- #

#Return the boost vector which makes the two four vectors closest to being back-to-back as possible
def getBoost(p1, p2):

    pTot3Vec = p1.Vect() + p2.Vect()

    #beta = pTot3Vec / pTot3Vec.Mag()

    boost3Vec = TVector3(pTot3Vec.X() / pTot3Vec.Mag(), pTot3Vec.Y() / pTot3Vec.Mag(), pTot3Vec.Z() / pTot3Vec.Mag())
    
    boost = TLorentzVector()
    boost.SetVectM(boost3Vec, 0)

    return boost

# ----------------------------------------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------------------------- #

class ZJetReclusterProducer(Module):

    def __init__(self):
        pass
    # ----------------------------------------------------------------------------------------------------------------------------- #

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)
    
    # ----------------------------------------------------------------------------------------------------------------------------- #

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("ZReClJ_pow", "I") #"The power used in the reclustering step (-1=anti-kt, 0=cambridge-aachen, 1=kt)"
        self.out.branch("ZReClJ_pt", "F") #"pt of the overall reclustered Z->AK8 Jet"
        self.out.branch("ZReClJ_eta", "F") #"eta of the overall reclustered Z->AK8 Jet"
        self.out.branch("ZReClJ_phi", "F") #"phi of the overall reclustered Z->AK8 Jet"
        self.out.branch("ZReClJ_mass", "F") #"m of the overall reclustered Z->AK8 Jet"
        self.out.branch("ZReClJ_nSJs", "I") #"The number of AK4 subjets after reclustering"
        self.out.branch("ZReClJ_sjPt", "F", lenVar="ZReClJ_nSJs") #"pt of the Z subjets (AK4 clustering results)"
        self.out.branch("ZReClJ_sjEta", "F", lenVar="ZReClJ_nSJs") #"eta of the Z subjets (AK4 clustering results)"
        self.out.branch("ZReClJ_sjPhi", "F", lenVar="ZReClJ_nSJs") #"phi of the Z subjets (AK4 clustering results)"
        self.out.branch("ZReClJ_sjMass", "F", lenVar="ZReClJ_nSJs") #"mass of the Z subjets (AK4 clustering results)"

    # ----------------------------------------------------------------------------------------------------------------------------- #

    def analyze(self, event):
        #POWER=-1 => anti-kt, POWER=0 => cambridge/aachen, POWER=1 =>inclusive kt
        power = 0
        pt = -999.99
        eta = -999.99
        phi = -999.99
        mass = -999.99
        nSJs = 0
        sjPt = []
        sjEta = []
        sjPhi = []
        sjMass = []


        if event.Z_jetIdxAK8 >=0 and event.Z_dm == 0 and event.Z_sJIdx1 >=0 and event.Z_sJIdx2 >=0:

            fjPFCands = Collection(event, "FatJetPFCands")
            pfCands = Collection(event, "PFCands")
            subJets = Collection(event, "SubJet")
            
            zSJ1 = subJets[event.Z_sJIdx1].p4()
            zSJ2 = subJets[event.Z_sJIdx2].p4()
            boost = getBoost(zSJ1, zSJ2)

            zJetPFCs = []
            for fjPFCand in fjPFCands:
                if fjPFCand.jetIdx == event.Z_jetIdxAK8:
                    pfCand = pfCands[fjPFCand.pFCandsIdx].p4()
                    zJetPFCs.append(pfCand + boost)
            
            if len(zJetPFCs) < 2:
                print("WARNING: Did not find at least 2 PFCs matching to Z AK8 jet!")
            else:
                reClAK4Jets = recluster(pfCands=zJetPFCs, p=power, R=0.4)

                reClAK8Jet = TLorentzVector()

                for reClAK4Jet in reClAK4Jets:
                    reClAK4Jet -= boost #Boost back to the lab frame

                    nSJs += 1
                    sjPt.append(reClAK4Jet.Pt())
                    sjEta.append(reClAK4Jet.Eta())
                    sjPhi.append(reClAK4Jet.Phi())
                    sjMass.append(reClAK4Jet.M())

                    reClAK8Jet += reClAK4Jet

                pt = reClAK8Jet.Pt()
                eta = reClAK8Jet.Eta()
                phi = reClAK8Jet.Phi()
                mass = reClAK8Jet.M()


        self.out.fillBranch("ZReClJ_pow", power)
        self.out.fillBranch("ZReClJ_pt", pt)
        self.out.fillBranch("ZReClJ_eta", eta)
        self.out.fillBranch("ZReClJ_phi", phi)
        self.out.fillBranch("ZReClJ_mass", mass)
        self.out.fillBranch("ZReClJ_nSJs", nSJs)
        self.out.fillBranch("ZReClJ_sjPt", sjPt)
        self.out.fillBranch("ZReClJ_sjEta", sjEta)
        self.out.fillBranch("ZReClJ_sjPhi", sjPhi)
        self.out.fillBranch("ZReClJ_sjMass", sjMass)

        return True

    # ----------------------------------------------------------------------------------------------------------------------------- #

zJetReclusterProducerConstr = lambda : ZJetReclusterProducer()
