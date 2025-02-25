#Boosted, Z-rest frame reclusterinf of PF candidates corresponding to Z->FatJet events selected by ZProducer
#Requires the information added by the bvnano-prod package: https://github.com/cms-btv-pog/btvnano-prod/tree/NanoAODv12_22Sep2023
# 

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

from ROOT import TLorentzVector, TVector3

# ----------------------------------------------------------------------------------------------------------------------------- #

#Return the boost vector which makes the two particle flow candidates, pfc1 and pfc2 closest to being back-to-back as possible
def getBoost(pfc1, pfc2):
    p1 = pfc1.p4()
    p2 = pfc2.p4()

    pTot3Vec = p1.Vect() + p2.Vect()

    beta = pTot3Vec / pTot3Vec.Mag()

    boost3Vec = TVector3(*beta)
    boost = TLorentzVector()
    boost.SetVectM(boost3Vec, 0)

    return boost

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
        self.out.branch("ZReCl_nPFCs", "I") #"Number of PF candidates corresponding included in the reclustered jet"
        self.out.branch("ZReCl_pfcIdxs", "I", lenVar="ZReCl_nPFCs") #"Indices to PFCands of the utilized PF candidates"
        self.out.branch("ZReCl_pfcBoostPt", "F", lenVar="ZReCl_nPFCs") #"PT of the PF candidates in the Z frame"
        self.out.branch("ZReCl_pfcBoostEta", "F", lenVar="ZReCl_nPFCs") #"eta of the PF candidates in the Z frame"
        self.out.branch("ZReCl_pfcBoostPhi", "F", lenVar="ZReCl_nPFCs") #"phi of the PF candidates in the Z frame"
        self.out.branch("ZReCl_mass", "F", lenVar="ZReCl_nPFCs") #"Mass of the reclusted jet. idx 0 is of leading PF cand, idx 1, leading+sub-leading PF cand, etc"
        self.out.branch("ZReCl_pt", "F", lenVar="ZReCl_nPFCs") #"pt of the reclusted jet. idx 0 is of leading PF cand, idx 1, leading+sub-leading PF cand, etc"
        self.out.branch("ZReCl_eta", "F", lenVar="ZReCl_nPFCs") #"eta of the reclusted jet. idx 0 is of leading PF cand, idx 1, leading+sub-leading PF cand, etc"
        self.out.branch("ZReCl_phi", "F", lenVar="ZReCl_nPFCs") #"phi of the reclusted jet. idx 0 is of leading PF cand, idx 1, leading+sub-leading PF cand, etc"

    # ----------------------------------------------------------------------------------------------------------------------------- #

    def analyze(self, event):
        MAX_PFCs = 6

        ZReCl_nPFCs = 0
        ZReCl_pfcIdxs = []
        ZReCl_pfcBoostPt  = []
        ZReCl_pfcBoostEta = []
        ZReCl_pfcBoostPhi = []
        ZReCl_mass = []
        ZReCl_pt = []
        ZReCl_eta = []
        ZReCl_phi = []
        
        #First, find all the PF candidates which are in the FatJet selected by ZProducer 
        fatJetPFCands = Collection(event, "FatJetPFCands")
        zJetPFCs = [] # list of (pfCandIdx, pfCandPt)
        for fJPFC in fatJetPFCands:
            if fJPFC.jetIdx == event.Z_jetIdxAK8:
                zJetPFCs.append((fJPFC.pFCandsIdx, fJPFC.pt))

        #Sort the located PF cands by descending pt (in lab frame)
        zJetPFCs.sort(key=lambda fJPFC : fJPFC[1], reverse=True) 

        if len(zJetPFCs) >= 2:
            pfCands = Collection(event, "PFCands")

            #Get the Lorentz boost which makes the two leading PF cands as close to back-to-back as possible
            #Apply the boost to the PFCs in the Z jet and resort by descending pT, this time in the Z rest frame
            boost = getBoost(pfc1=pfCands[zJetPFCs[0][0]], pfc2=pfCands[zJetPFCs[1][0]])
            boostedPFCs = []
            for pfcN in range(len(zJetPFCs)):
                pfcIdx = zJetPFCs[pfcN][0]
                pfc = pfCands[pfcIdx]
                boostedPFCs.append(( pfcIdx, pfCands[pfcIdx].p4() + boost))

            boostedPFCs.sort(key= lambda pfc : pfc[1].Pt(), reverse=True)


            reClZJet = TLorentzVector()
            for pfCandN in range( min(len(boostedPFCs), MAX_PFCs) ):
                ZReCl_nPFCs += 1
                ZReCl_pfcIdxs.append(boostedPFCs[pfCandN][0])
            
                pfCand = boostedPFCs[pfCandN][1]

                ZReCl_pfcBoostPt.append(pfCand.Pt())
                ZReCl_pfcBoostEta.append(pfCand.Eta())
                ZReCl_pfcBoostPhi.append(pfCand.Phi())

                if pfCandN == 0: 
                    #The overall "jet" is just the single pf candidate for the first pfc
                    ZReCl_pt.append(pfCand.Pt())
                    ZReCl_eta.append(pfCand.Eta())
                    ZReCl_phi.append(pfCand.Phi())
                    ZReCl_mass.append(pfCand.M())
                    reClZJet.SetPtEtaPhiM(ZReCl_pt[-1], ZReCl_eta[-1], ZReCl_phi[-1], ZReCl_mass[-1])
                
                else: #In all other iterations, add the next highest pT pfc to the sum of previous
                    reClZJet = reClZJet + pfCand
                    ZReCl_mass.append(reClZJet.M())
                    ZReCl_pt.append(reClZJet.Pt())
                    ZReCl_eta.append(reClZJet.Eta())
                    ZReCl_phi.append(reClZJet.Phi())

        else:
            print("WARNING: <2 Z jet PF cands found! Cannot perform boosed reclustering!")

        self.out.fillBranch("ZReCl_nPFCs", ZReCl_nPFCs)
        self.out.fillBranch("ZReCl_pfcIdxs", ZReCl_pfcIdxs)
        self.out.fillBranch("ZReCl_pfcBoostPt", ZReCl_pfcBoostPt)
        self.out.fillBranch("ZReCl_pfcBoostEta", ZReCl_pfcBoostEta)
        self.out.fillBranch("ZReCl_pfcBoostPhi", ZReCl_pfcBoostPhi)
        self.out.fillBranch("ZReCl_mass", ZReCl_mass)
        self.out.fillBranch("ZReCl_pt", ZReCl_pt)
        self.out.fillBranch("ZReCl_eta", ZReCl_eta)
        self.out.fillBranch("ZReCl_phi", ZReCl_phi)

        return True

    # ----------------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------------------------------------- #
zJetReclusterProducerConstr = lambda : ZJetReclusterProducer()
# ----------------------------------------------------------------------------------------------------------------------------- #