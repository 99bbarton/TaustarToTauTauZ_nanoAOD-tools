from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

# -----------------------------------------------------------------------------------------------------------------------------

class ZProducer(Module):

    def __init__(self, era):
        self.era = era
        pass
    
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Z_dm", "I") #"0 = hadronic, 1=electrons, 2=muons, 3=taus. -1 by default"
        self.out.branch("Z_d1Idx", "I") #"Idx to first daughter of Z in either Muons or Electrons collection (if Z_dm == 1 or Z_dm == 2). -1 default"
        self.out.branch("Z_d2Idx", "I") #"Idx to second daughter of Z in either Muons or Electrons collection (depending on if Z_dm == 1 or Z_dm == 2). -1 default"
        self.out.branch("Z_dauDR", "F") #"DeltaR(zDau1, zDau2). 0 default"
        self.out.branch("Z_mass", "F") #"Mass of ee or mumu pair if either Z_dm == 1 or Z_dm == 2 or jet if Z_dm =0. 0 default"
        self.out.branch("Z_pt", "F") #"Pt of ee or mumu pair or jet. 0 default"
        self.out.branch("Z_jetIdxDT", "I") #"Idx to FatJet collection of the most Z-like jet determined using deepTag score. if Z_dm=0"
        self.out.branch("Z_jetIdxPN", "I") #"Idx to FatJet collection of the most Z-like jet determined using particle net score. if Z_dm=0"
        self.out.branch("Z_nEE", "I") #"Number of candidate ee pairs found"
        self.out.branch("Z_nMuMu", "I") #"Number of candidate mumu pairs found"

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        Z_dm = -1
        Z_d1Idx = -1
        Z_d2Idx = -1
        Z_mass = 0
        Z_pt = 0
        Z_dauDR = 0
        Z_jetIdxDT = -1
        Z_jetIdxPN = -1
        Z_nEE = 0
        Z_nMuMu = 0

        electrons = Collection(event, "Electron")
        for e1Idx, e1 in enumerate(electrons):
            for e2Idx, e2 in enumerate(electrons):
                if e1Idx == e2Idx:
                    continue
                if self.era == 2:
                    cuts = (e1.charge * e2.charge) < 0 #Opposite charge
                    cuts = cuts and (abs(e1.eta + e1.deltaEtaSC) >= 1.566 or abs( e1.eta + e1.deltaEtaSC) < 1.444)#Fiducial
                    cuts = cuts and (abs(e2.eta + e2.deltaEtaSC) >= 1.566 or abs(e2.eta + e2.deltaEtaSC) < 1.444)
                    cuts = cuts and (e1.pt >= 20.0 and abs(e1.eta + e1.deltaEtaSC) < 2.5 and e1.mvaFall17V2noIso_WP80) #ID
                    cuts = cuts and (e2.pt >= 20.0 and abs(e2.eta + e2.deltaEtaSC) < 2.5 and e2.mvaFall17V2noIso_WP80)
                elif self.era == 3:
                    cuts = (e1.charge * e2.charge) < 0 #Opposite charge
                    cuts = cuts and (abs(e1.eta + e1.deltaEtaSC) >= 1.566 or abs( e1.eta + e1.deltaEtaSC) < 1.444) #Fiducial
                    cuts = cuts and (abs(e2.eta + e2.deltaEtaSC) >= 1.566 or abs(e2.eta + e2.deltaEtaSC) < 1.444)
                    cuts = cuts and (e1.pt >= 20.0 and abs(e1.eta + e1.deltaEtaSC) < 2.5 and e1.mvaNoIso_WP80 ) #ID (basic)
                    cuts = cuts and (e2.pt >= 20.0 and abs(e2.eta + e2.deltaEtaSC) < 2.5 and e2.mvaNoIso_WP80) #ID (basic)

                if cuts:
                    Z_nEE += 1
                    tempM = (e1.p4() + e2.p4()).M()
                    if (tempM - 91.18) < (Z_mass - 91.18): #If this pair is closer in mass to the nominal Z
                        Z_mass = (e1.p4()+e2.p4()).M()
                        Z_pt = (e1.p4()+e2.p4()).Pt()
                        if Z_mass >= 60.0 and Z_mass < 120.0: 
                            Z_dm = 1
                            Z_d1Idx = e1Idx if e1.pt >= e2.pt else e2Idx  #daughter 1 is higher pT e
                            Z_d2Idx = e2Idx if e1.pt >= e2.pt else e1Idx 
                            Z_dauDR = e1.DeltaR(e2)
        
        
        muons = Collection(event, "Muon")
        for mu1Idx, mu1 in enumerate(muons):
            for mu2Idx, mu2 in enumerate(muons):
                if mu2Idx == mu1Idx:
                    continue
                if (mu1.charge * mu2.charge < 0): #Opposite charge
                    if (mu1.pt >= 15.0 and abs(mu1.eta) < 2.4 and mu1.mediumId) and (mu2.pt >= 15.0 and abs(mu2.eta) < 2.4 and mu2.mediumId): #ID
                        Z_nMuMu += 1
                        tempM = (mu1.p4() + mu2.p4()).M()
                        if (tempM - 91.18) < (Z_mass - 91.18):
                            Z_mass = (mu1.p4()+mu2.p4()).M()
                            Z_pt = (mu1.p4()+mu2.p4()).Pt()
                            if Z_mass >= 60.0 and Z_mass < 120.0: 
                                Z_dm = 2
                                Z_d1Idx = mu1Idx if mu1.pt >= mu2.pt else mu2Idx  #daughter 1 is higher pT mu
                                Z_d2Idx = mu2Idx if mu1.pt >= mu2.pt else mu1Idx
                                Z_dauDR = mu1.DeltaR(mu2)

                        
        #TODO Add Z-tautau

        
        if Z_dm < 0:
            fatJets = Collection(event, "FatJet")
            score_dt = 0 #Best DeepBoostedJet tagger score
            score_pn = 0 #Best ParticleNet tagger score
            for jetIdx, jet in enumerate(fatJets):
                if jet.mass >= 40 and jet.mass <= 150:
                    if self.era == 2:
                        if jet.deepTag_ZvsQCD > 0.7 and jet.deepTag_ZvsQCD > score_dt:
                            score_dt = jet.deepTag_ZvsQCD
                            Z_jetIdxDT = jetIdx
                        if jet.particleNet_ZvsQCD > 0.7 and jet.particleNet_ZvsQCD > score_pn:
                            score_pn = jet.particleNet_ZvsQCD
                            Z_jetIdxPN = jetIdx
                    elif self.era ==3:
                        if jet.particleNetWithMass_ZvsQCD > 0.7 and jet.particleNetWithMass_ZvsQCD > score_pn:
                            score_pn = jet.particleNetWithMass_ZvsQCD
                            Z_jetIdxPN = jetIdx

            if score_dt > score_pn and score_dt > 0.7: #Use the most confident score for now until 
                Z_mass = fatJets[Z_jetIdxDT].mass
                Z_pt = fatJets[Z_jetIdxDT].pt
                Z_dm = 0
            elif score_pn > 0.7:
                Z_mass = fatJets[Z_jetIdxPN].mass
                Z_pt = fatJets[Z_jetIdxPN].pt
                Z_dm = 0
        
        self.out.fillBranch("Z_dm", Z_dm)
        self.out.fillBranch("Z_d1Idx", Z_d1Idx)
        self.out.fillBranch("Z_d2Idx", Z_d2Idx)
        self.out.fillBranch("Z_mass", Z_mass)
        self.out.fillBranch("Z_pt", Z_pt)
        self.out.fillBranch("Z_dauDR", Z_dauDR)
        self.out.fillBranch("Z_jetIdxDT", Z_jetIdxDT)
        self.out.fillBranch("Z_jetIdxPN", Z_jetIdxPN)
        self.out.fillBranch("Z_nEE", Z_nEE)
        self.out.fillBranch("Z_nMuMu", Z_nMuMu)

        return True
    
    # -----------------------------------------------------------------------------------------------------------------------------

zProducerConstr = lambda era: ZProducer(era = era)

#from PhysicsTools.NanoAODTools.postprocessing.modules.GenProducerZTau import genProducerZTauConstr

#files = ["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/TauZ/taustarToTauZ_m500_2018.root"]
#p = PostProcessor(".", files, cut="1>0", branchsel=None, postfix="", modules=[genProducerZTauConstr(), zProducerConstr()] )
#p.run()
