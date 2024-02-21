#Identify taus passing basic ID requirements which have pT>20 and abs(eta)<2.3

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

# ----------------------------------------------------------------------------------------------------------------------------

class TauProducer(Module):

    def __init__(self):
        pass
    
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("SelTaus_n", "I") #"The number of taus which passed basic requirements"
        self.out.branch("SelTaus_idxs", "i", lenVar="SelTaus_n") #"List of length SelTaus_n of indices to Tau collection of taus passing basic requirements"

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        taus = Collection(event, "Tau")
        goodTauIdxs = []
        for idx, tau in enumerate(taus):
            tauID = tau.decayMode != 5 and tau.decayMode != 6 and tau.decayMode != 7
            tauID = tauID and (tau.pt >= 20.0) and (abs(tau.eta) < 2.3)
            tauID = tauID and (32 & tau.idDeepTau2017v2p1VSjet) # 32 = Tight
            tauID = tauID and (8 & tau.idDeepTau2017v2p1VSmu) # 8 = Tight WP
            tauID = tauID and (32 & tau.idDeepTau2017v2p1VSe) # 32 = Tight WP
            if tauID:
                goodTauIdxs.append(idx)
            
        self.out.fillBranch("SelTaus_n", len(goodTauIdxs))
        self.out.fillBranch("SelTaus_idxs", goodTauIdxs)

        return True
    
# ----------------------------------------------------------------------------------------------------------------------------
    
tauProducerConstr = lambda: TauProducer()

# ----------------------------------------------------------------------------------------------------------------------------
files = ["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/TauZ/taustarToTauZ_m3000_2018.root"]
p = PostProcessor(".", files, cut="1>0", branchsel=None, modules=[tauProducerConstr()] )
p.run()
