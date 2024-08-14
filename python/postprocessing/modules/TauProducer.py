#Identify taus passing basic ID requirements 

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

# ----------------------------------------------------------------------------------------------------------------------------

class TauProducer(Module):

    def __init__(self, era):
        self.era = era
    
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("SelTaus_n", "I") #"The number of taus which passed basic requirements"
        self.out.branch("SelTaus_idxs", "i", lenVar="SelTaus_n") #"List of length SelTaus_n of indices to Tau collection of taus passing basic requirements"
        self.out.branch("SelTaus_nProngs", "I", lenVar="SelTaus_n") #"List if len SelTaus_n of the number of prongs in the decay of taus passing basic cuts"
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        taus = Collection(event, "Tau")
        goodTauIdxs = []
        nProngs = []
        for idx, tau in enumerate(taus):
            if self.era == 2: #https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
                tauID = tau.decayMode != 5 and tau.decayMode != 6 and tau.decayMode != 7
                tauID = tauID and tau.pt > 20.0 and abs(tau.eta) < 2.3 and tau.dz < 0.2
                tauID = tauID and (16 & tau.idDeepTau2017v2p1VSjet) # 16 = med, 32 = Tight
                tauID = tauID and (4 & tau.idDeepTau2017v2p1VSmu) # 4 = med, 8 = Tight WP
                tauID = tauID and (16 & tau.idDeepTau2017v2p1VSe) # 16= med, 32 = Tight WP
            elif self.era == 3: #https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3
                tauID = tau.decayMode != 5 and tau.decayMode != 6 and tau.decayMode != 7
                tauID = tauID and (tau.pt > 20.0) and (abs(tau.eta) < 2.5) and tau.dz < 0.2
                tauID = tauID and tau.idDeepTau2018v2p5VSjet == 5 #5= medium
                tauID = tauID and tau.idDeepTau2018v2p5VSmu == 3 #3= medium
                tauID = tauID and tau.idDeepTau2018v2p5VSe == 5 #5= medium
            if tauID:
                goodTauIdxs.append(idx)
                nP = -1
                if tau.decayMode <= 2 and tau.decayMode > 0:
                    nP = 1
                elif tau.decayMode >= 10:
                    nP = 3
                nProngs.append(nP)
                
        self.out.fillBranch("SelTaus_n", len(goodTauIdxs))
        self.out.fillBranch("SelTaus_idxs", goodTauIdxs)
        self.out.fillBranch("SelTaus_nProngs", nProngs)
        
        return True
    
# ----------------------------------------------------------------------------------------------------------------------------
    
tauProducerConstr = lambda era: TauProducer(era = era)

# ----------------------------------------------------------------------------------------------------------------------------
files = ["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/TauZ/taustarToTauZ_m3000_2018.root"]
p = PostProcessor(".", files, cut="1>0", branchsel=None, modules=[tauProducerConstr()] )
p.run()
