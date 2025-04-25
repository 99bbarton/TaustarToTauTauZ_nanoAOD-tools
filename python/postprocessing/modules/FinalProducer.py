#Intended as the final stage of processing.
# Checks that the isCand flag is true for one and only one channel
# For background, checks if the identified particles gen match (i.e. Z is real, taus are real, etc)

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

# ----------------------------------------------------------------------------------------------------------------------------

class FinalProducer(Module):
    placeholder = False
    
    def __init__(self):
        self.placeholder = False

    # ----------------------------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------------

    def endJob(self):
        pass

    # ----------------------------------------------------------------------------------------------------------------------------

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Fin_ch", "I") #"0=tautau, 1=etau, 2=mutau. -1 otherwise. Includes check that ==1 CH_isCand is TRUE"
        self.out.branch("Fin_justRel", "O") #"True if there are no extra leptons/tauhs/z"
        #self.out.branch("Fin_zReal", "O") #"For sig, True if Z gen matches to taustar. For bkgd, True if real Z"
        #self.out.branch("Fin_leg1Real", "O") #"For sig, True if tauh (tau1) gen matches to taustar. For bkgd, True if real tauh"
        #self.out.branch("Fin_leg2Real", "O") #"For sig, True if e/mu/tau2 gen matches to taustar. For bkgd, True if real taue/taumu/tauh"

    # ----------------------------------------------------------------------------------------------------------------------------


    def analyze(self, event):
        ch = -1
        justRel = False
        zReal = False
        leg1Real = False
        leg2Real = False


        #First, reject events which do not have a valid tautauZ triplet
        if not event.Z_isCand:
            return False

        if not hasattr(event, "TauTau_visM"):
            print("No TauTau_isCand")
            return False        
        
        if event.TauTau_isCand and not (event.ETau_isCand or event.MuTau_isCand):
            ch = 0
            if event.Z_dm == 0:
                justRel = (event.nElectron == 0 and event.nMuon == 0 and event.Z_nJetCands == 1)
            elif event.Z_dm == 1:
                justRel = (event.nElectron == 2 and event.nMuon == 0 and event.Z_nJetCands == 0)
            elif event.Z_dm == 2:
                justRel = (event.nElectron == 0 and event.nMuon == 2 and event.Z_nJetCands == 0)
        elif event.ETau_isCand and not (event.TauTau_isCand or event.MuTau_isCand):
            ch = 1
            if event.Z_dm == 0:
                justRel = (event.nElectron == 1 and event.nMuon == 0 and event.Z_nJetCands == 1)
            elif event.Z_dm == 1:
                justRel = (event.nElectron == 3 and event.nMuon == 0 and event.Z_nJetCands == 0)
            elif event.Z_dm == 2:
                justRel = (event.nElectron == 1 and event.nMuon == 2 and event.Z_nJetCands == 0)
        elif event.MuTau_isCand and not (event.TauTau_isCand or event.ETau_isCand):
            ch = 2
            if event.Z_dm == 0:
                justRel = (event.nElectron == 0 and event.nMuon == 1 and event.Z_nJetCands == 1)
            elif event.Z_dm == 1:
                justRel = (event.nElectron == 2 and event.nMuon == 1 and event.Z_nJetCands == 0)
            elif event.Z_dm == 2:
                justRel = (event.nElectron == 0 and event.nMuon == 3 and event.Z_nJetCands == 0)

        if ch < 0:
            return False

        #Now, reject events that don't look signal like or have detector issues
        #TODO
        #e.g. passing Flags,etc


        #Check if the reco objects gen match
        #TODO

        self.out.fillBranch("Fin_ch", ch)
        self.out.fillBranch("Fin_justRel", justRel)

        return True
    # ----------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------------

finalProducerConstr = lambda: FinalProducer()

