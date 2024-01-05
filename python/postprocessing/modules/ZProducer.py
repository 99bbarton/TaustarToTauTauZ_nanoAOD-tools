from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

# -----------------------------------------------------------------------------------------------------------------------------

class ZProducer(Module):

    def __init__(self):
        pass
    
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Z_ee", "O") #"True if a ee pair is found, false otherwise"
        self.out.branch("Z_mumu", "O") #"True if a mumu pair is found, false otherwise"
        self.out.branch("Z_d1Idx", "I") #"Idx to first daughter of Z in either Muons or Electrons collection (depending on if Z_ee or Z_mumu are true). -1 default"
        self.out.branch("Z_d2Idx", "I") #"Idx to second daughter of Z in either Muons or Electrons collection (depending on if Z_ee or Z_mumu are true). -1 default"
        self.out.branch("Z_pairMass", "F") #"Mass of ee or mumu pair"
        self.out.branch("Z_pairPt", "F") #"Pt of ee or mumu pair"

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        Z_ee = False
        Z_mumu = False
        Z_d1Idx = -1
        Z_d2Idx = -1
        Z_pairMass = 0
        Z_pairPt = 0

        electrons = Collection(event, "Electron")
        for e1Idx, e1 in enumerate(electrons):
            for e2Idx, e2 in enumerate(electrons):
                cuts = True and (e1.charge * e2.charge < 0) #Opposite charge
                cuts = cuts and (abs(e1.eta + e1.deltaEtaSC) >= 1.566 or abs( e1.eta + e1.deltaEtaSC) < 1.444) #Fiducial
                cuts = cuts and (abs(e2.eta + e2.deltaEtaSC) >= 1.566 or abs(e2.eta + e2.deltaEtaSC) < 1.444)
                cuts = cuts and (e1.pt >= 20.0 and abs(e1.eta) < 2.5 and e1.mvaFall17V2noIso_WPL) #ID (basic)
                cuts = cuts and (e2.pt >= 20.0 and abs(e2.eta) < 2.5 and e2.mvaFall17V2noIso_WPL) #ID (basic)

                if cuts:
                    tempPt = (e1.p4() + e2.p4()).Pt()
                    if tempPt > Z_pairPt:
                        Z_pairMass = (e1.p4()+e2.p4()).M()
                        Z_pairPt = (e1.p4()+e2.p4()).Pt()
                        if Z_pairMass >= 50.0 and Z_pairMass < 140.0: 
                            Z_ee = True
                            Z_d1Idx = e1Idx if e1.pt >= e2.pt else e2Idx  #daughter 1 is higher pT e
                            Z_d2Idx = e2Idx if e1.pt >= e2.pt else e1Idx 

        muons = Collection(event, "Muon")
        for mu1Idx, mu1 in enumerate(muons):
            for mu2Idx, mu2 in enumerate(muons):
                if (mu1.charge * mu2.charge < 0): #Opposite charge
                    if (mu1.pt >= 15.0 and abs(mu1.eta) < 2.4 and mu1.looseId) and (mu2.pt >= 15.0 and abs(mu2.eta) < 2.4 and mu2.looseId): #ID
                        tempPt = (mu1.p4() + mu2.p4()).Pt()
                        if tempPt > Z_pairPt:
                            Z_pairMass = (mu1.p4()+mu2.p4()).M()
                            Z_pairPt = (mu1.p4()+mu2.p4()).Pt()
                            if Z_pairMass >= 50.0 and Z_pairMass < 140.0: 
                                Z_mumu = True
                                Z_d1Idx = mu1Idx if mu1.pt >= mu2.pt else mu2Idx  #daughter 1 is higher pT e
                                Z_d2Idx = mu2Idx if mu1.pt >= mu2.pt else mu1Idx
        
        self.out.fillBranch("Z_ee", Z_ee)
        self.out.fillBranch("Z_mumu", Z_mumu)
        self.out.fillBranch("Z_d1Idx", Z_d1Idx)
        self.out.fillBranch("Z_d2Idx", Z_d2Idx)
        self.out.fillBranch("Z_pairMass", Z_pairMass)
        self.out.fillBranch("Z_pairPt", Z_pairPt)

        return True
    
    # -----------------------------------------------------------------------------------------------------------------------------

    zProducerConstr = lambda: ZProducer()
    

