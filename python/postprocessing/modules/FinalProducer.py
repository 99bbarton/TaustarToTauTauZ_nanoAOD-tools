#Intended as the final stage of processing.
# Checks that the isCand flag is true for one and only one channel
# For background, checks if the identified particles gen match (i.e. Z is real, taus are real, etc)

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

# ----------------------------------------------------------------------------------------------------------------------------

class FinalProducer(Module):

    def __init__(self):
        pass

    # ----------------------------------------------------------------------------------------------------------------------------
    def beginJob(self):
        pass

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

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

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

        if event.TauTau_isCand and not (event.ETau_isCand or event.MuTau_isCand):
            ch = 0
            if event.Z_dm == 0:
                justRel = (event.nElectrons == 0 and event.nMuons == 0 and event.Z_nJetCands == 1)
            elif event.Z_dm == 1:
                justRel = (event.nElectrons == 2 and event.nMuons == 0 and event.Z_nJetCands == 0)
            elif event.Z_dm == 2:
                justRel = (event.nElectrons == 0 and event.nMuons == 2 and event.Z_nJetCands == 0)
        elif event.ETau_isCand and not (event.TauTau_isCand or event.MuTau_isCand):
            ch = 1
            if event.Z_dm == 0:
                justRel = (event.nElectrons == 1 and event.nMuons == 0 and event.Z_nJetCands == 1)
            elif event.Z_dm == 1:
                justRel = (event.nElectrons == 3 and event.nMuons == 0 and event.Z_nJetCands == 0)
            elif event.Z_dm == 2:
                justRel = (event.nElectrons == 1 and event.nMuons == 2 and event.Z_nJetCands == 0)
        elif event.MuTau_isCand and not (event.TauTau_isCand or event.ETau_isCand):
            ch = 2
            if event.Z_dm == 0:
                justRel = (event.nElectrons == 0 and event.nMuons == 1 and event.Z_nJetCands == 1)
            elif event.Z_dm == 1:
                justRel = (event.nElectrons == 2 and event.nMuons == 1 and event.Z_nJetCands == 0)
            elif event.Z_dm == 2:
                justRel = (event.nElectrons == 0 and event.nMuons == 3 and event.Z_nJetCands == 0)

        if ch < 0:
            return False

        #Now, reject events that don't look signal like or have detector issues
        #TODO
        #e.g. passing Flags,etc


        #Check if the reco objects gen match
        #TODO

        self.out.fillBranch("Fin_ch", ch)
        self.out.fillBranch("Fin_", justRel)

        return True
    # ----------------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------------

finalProducerConstr = lambda: FinalProducer()

