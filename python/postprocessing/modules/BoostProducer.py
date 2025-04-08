#Producer to store information of the Z->jet events in the Z rest frame

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

from ROOT import TLorentzVector

# ----------------------------------------------------------------------------------------------------------------------------- #

class BoostProducer(Module):

    def __init__(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Boost_n", "I") #"The number of boosed jets stored (i.e. len of array variables below). Set to 3 currently."
        self.out.branch("Boost_pt", "F", lenVar="Boost_n") #"Boosted jet pt. First idx is parent FatJet, second idx is leading SubJet, third sub-leading SubJet"
        self.out.branch("Boost_eta", "F", lenVar="Boost_n") #"Boosted jet eta. First idx is parent FatJet, second idx is leading SubJet, third sub-leading SubJet"
        self.out.branch("Boost_phi", "F", lenVar="Boost_n") #"Boosted jet phi. First idx is parent FatJet, second idx is leading SubJet, third sub-leading SubJet"
        self.out.branch("Boost_dR", "F") #"DeltaR between the boosted subJets"
        self.out.branch("Boost_dPhi", "F") #"DeltaPhi between the boosted subJets"
    
    def analyze(self, event):
        Boost_n = 3
        Boost_pt = [0, -999, -999]
        Boost_eta = [0, -999, -999]
        Boost_phi = [0, -999, -999]
        Boost_dR = -999
        Boost_dPhi = -999

        if event.Z_dm == 0 and event.Z_sJIdx1 >= 0 and event.Z_sJIdx2 >=0:
            zFatJet = datamodel.Object(event, "FatJet", event.Z_jetIdxAK8).p4()
            zSubJet1 = datamodel.Object(event, "SubJet", event.Z_sJIdx1).p4()
            zSubJet2 = datamodel.Object(event, "SubJet", event.Z_sJIdx2).p4()

            boostVec = zFatJet.BoostVector()

            #zFatJet_boosted = TLorentzVector(zFatJet)
            #zFatJet_boosted.Boost(-boostVec)
            #Boost_pt[0] = zFatJet_boosted.Pt()
            #Boost_eta[0] = zFatJet_boosted.Eta()
            #Boost_phi[0] = zFatJet_boosted.Phi()

            zSubJet1_boosted = TLorentzVector(zSubJet1)
            zSubJet1_boosted.Boost(-boostVec)
            Boost_pt[1] = zSubJet1_boosted.Pt()
            Boost_eta[1] = zSubJet1_boosted.Eta()
            Boost_phi[1] = zSubJet1_boosted.Phi()

            zSubJet2_boosted = TLorentzVector(zSubJet2)
            zSubJet2_boosted.Boost(-boostVec)
            Boost_pt[2] = zSubJet2_boosted.Pt()
            Boost_eta[2] = zSubJet2_boosted.Eta()
            Boost_phi[2] = zSubJet2_boosted.Phi()

            Boost_dPhi = zSubJet1_boosted.DeltaPhi(zSubJet2_boosted)
            Boost_dR = zSubJet1_boosted.DeltaR(zSubJet2_boosted)
            
        self.out.fillBranch("Boost_n", Boost_n)
        self.out.fillBranch("Boost_pt", Boost_pt)
        self.out.fillBranch("Boost_eta", Boost_eta)
        self.out.fillBranch("Boost_phi", Boost_phi)
        self.out.fillBranch("Boost_dR", Boost_dR)
        self.out.fillBranch("Boost_dPhi", Boost_dPhi)

        return True


# -----------------------------------------------------------------------------------------------------------------------------

boostProducerConstr = lambda : BoostProducer()
