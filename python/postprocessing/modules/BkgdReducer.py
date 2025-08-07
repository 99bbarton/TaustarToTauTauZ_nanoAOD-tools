from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection



class BkgdReducer(Module):

    def __init__(self):
        pass
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("BR_nBTags", "I")  #"Number of b-tagged jets in the event"
        self.out.branch("BR_nLMatch", "O") #"True if the number of leptons in the event matches the channel and Z DM"
        self.out.branch("BR_nJetsMatch", "O") #"True if the number of jets in th event matches the channel and Z DM"

    def analyze(self, event):
        nBTags = 0
        nLMatch = False
        nJetsMatch = False

        nJets = 0
        ak8Jets = Collection(event, "FatJet")
        for jet in ak8Jets:
            if abs(jet.eta) < 2.5 and jet.pt > 100 and jet.jetId == 6:
                nJets += 1
                if jet.btagDeepB > 0.7:
                    nBTags += 1
        nExp_jet = 0
        if event.ETau_isCand or event.MuTau_isCand:
            nExp_jet += 1
        if event.TauTau_isCand:
            nExp_jet += 2
        if event.Z_dm == 0:
            nExp_jet += 2
        nJetsMatch = nExp_jet == nJets

        nExp_el = 0
        nEl = 0
        if event.Z_dm == 1:
            nExp_el += 2
        if event.ETau_isCand:
            nExp_el += 1
        electrons = Collection(event, "Electron")
        for el in electrons:
            if el.pt > 24 and (abs(el.eta + el.deltaEtaSC) < 1.444 or (abs(el.eta + el.deltaEtaSC) > 1.566 and abs(el.eta + el.deltaEtaSC) < 2.5)):
                nEl += 1
        nExp_mu = 0
        nMu = 0
        if event.Z_dm == 2:
            nExp_mu += 2
        if event.MuTau_isCand:
            nExp_mu+= 1
        muons = Collection(event, "Muon")
        for mu in muons:
            if mu.pt > 20.0 and abs(mu.eta) < 2.4 and mu.mediumId:
                nMu += 1

        nLMatch == (nEl == nExp_el) and (nMu == nExp_mu)      

        self.out.fillBranch("BR_nBTags", nBTags)
        self.out.fillBranch("BR_nLMatch", nLMatch)
        self.out.fillBranch("BR_nJetsMatch", nJetsMatch)

        return True
    

bkgdReducerConstr = lambda: BkgdReducer()

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.BkgdReducer import bkgdReducerConstr

files = ["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/BackgroundMC/PFNano/2022post/V0/TTto2L2Nu_2022post.root"]
#files = ["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/TauZ/27Feb2025/taustarToTauZ_m3000_2022post.root"]
p = PostProcessor(".", files, cut="", branchsel=None, postfix="", modules=[bkgdReducerConstr()])
p.run()
