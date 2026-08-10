

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

# ----------------------------------------------------------------------------------------------------------------------------

class SDCMProducer(Module):

    def __init__(self):
        pass
    
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("SDCM_min", "F") #"Min collinear mass for whichever channel isCand[1] is true for."
        self.out.branch("SDCM_max", "F") #"Max collinear mass for whichever channel isCand[1] is true for."

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass


    def analyze(self, event):
        minM = -999.99
        maxM = -999.99

        if event.ETau_isCand[1]:
            minM = event.ETau_minCollM[1]
            maxM = event.ETau_maxCollM[1]
        elif event.MuTau_isCand[1]:
            minM = event.MuTau_minCollM[1]
            maxM = event.MuTau_maxCollM[1]
        elif event.TauTau_isCand[1]:
            minM = event.TauTau_minCollM[1]
            maxM = event.TauTau_maxCollM[1]

        self.out.fillBranch("SDCM_min", minM)
        self.out.fillBranch("SDCM_max", maxM)

        return True
    

sDCMProducerConstr = lambda: SDCMProducer()

