from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

from ROOT import TH1F

# -----------------------------------------------------------------------------------------------------------------------------

class ZProducer(Module):

    def __init__(self, era):
        self.era = era
        self.writeHistFile = True

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)
        self.h_ak4Mass = TH1F('h_ak4Mass', 'Mass of AK4 ;Mass [GeV];# of Jets', 75, 0, 150)
        self.addObject(self.h_ak4Mass)
        self.h_ak4MassCuts = TH1F('h_ak4MassCuts', 'Mass of AK4 Jets Passing ID Reqs.;Mass [GeV];# of Jets', 30, 60, 120)
        self.addObject(self.h_ak4MassCuts)
        self.h_ak8Mass = TH1F('h_ak8Mass', 'Mass of AK8 Jets;Mass [GeV];# of Jets', 75, 0, 150)
        self.addObject(self.h_ak8Mass)
        self.h_ak8MassCuts = TH1F('h_ak8MassCuts', 'Mass of AK8 Jets Passing ID Reqs;Mass [GeV];# of Jets', 30, 60, 120)
        self.addObject(self.h_ak8MassCuts)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Z_dm", "I") #"0 = hadronic, 1=electrons, 2=muons, 3=taus. -1 by default"
        self.out.branch("Z_d1Idx", "I") #"Idx to first daughter of Z in either Muons or Electrons collection (if Z_dm == 1 or Z_dm == 2). -1 default"
        self.out.branch("Z_d2Idx", "I") #"Idx to second daughter of Z in either Muons or Electrons collection (depending on if Z_dm == 1 or Z_dm == 2). -1 default"
        self.out.branch("Z_dauID", "O") #"For Z_dm=1 or 2, True if both Z daughter candidates pass the appropriate ID" 
        self.out.branch("Z_dauDR", "F") #"DeltaR(zDau1, zDau2). 0 default"
        self.out.branch("Z_mass", "F") #"Mass of ee or mumu pair if either Z_dm == 1 or Z_dm == 2 or jet if Z_dm =0. 0 default"
        self.out.branch("Z_pt", "F") #"Pt of ee or mumu pair or jet. 0 default"
        self.out.branch("Z_eta", "F") #"Eta of the overall Z candidate"
        self.out.branch("Z_phi", "F") #"Phi of the overall Z candidate"
        self.out.branch("Z_jetIdxAK8", "I") #"Idx to FatJet collection of the most Z-like AK8 jet (using particle net). if Z_dm=0"
        self.out.branch("Z_jetIdxAK4", "I") #"Idx to Jet collection of the most Z-like AK4 jet (using particle net). if Z_dm=0"
        self.out.branch("Z_jetR", "I") #"Either 4 or 8 if the best matching jet is an AK4/8 jet respectively if Z_Dm=0. -1 otherwise"
        self.out.branch("Z_nEE", "I") #"Number of candidate ee pairs found"
        self.out.branch("Z_nMuMu", "I") #"Number of candidate mumu pairs found"
        self.out.branch("Z_isCand", "O") #"True if Z is a good candidate (DM, mass, etc)"

    def analyze(self, event):
        Z_dm = -1
        Z_d1Idx = -1
        Z_d2Idx = -1
        Z_mass = 0
        Z_pt = 0
        Z_eta = -999.99
        Z_phi = -999.99
        Z_dauDR = 0
        Z_dauID = False
        Z_jetIdxAK8 = -1
        Z_jetIdxAK4 = -1
        Z_jetR = -1
        Z_nEE = 0
        Z_nMuMu = 0
        Z_isCand = False

        electrons = Collection(event, "Electron")
        for e1Idx, e1 in enumerate(electrons):
            for e2Idx in range(e1Idx + 1, event.nElectron):
                e2 = electrons[e2Idx]
                if self.era == 2:
                    cuts = (e1.charge * e2.charge) < 0 #Opposite charge
                    cuts = cuts and (abs(e1.eta + e1.deltaEtaSC) >= 1.566 or abs( e1.eta + e1.deltaEtaSC) < 1.444)#Fiducial
                    cuts = cuts and (abs(e2.eta + e2.deltaEtaSC) >= 1.566 or abs(e2.eta + e2.deltaEtaSC) < 1.444)
                    #cuts = cuts and (e1.pt >= 20.0 and abs(e1.eta + e1.deltaEtaSC) < 2.5 and e1.mvaFall17V2noIso_WP80) #ID
                    #cuts = cuts and (e2.pt >= 20.0 and abs(e2.eta + e2.deltaEtaSC) < 2.5 and e2.mvaFall17V2noIso_WP80)
                elif self.era == 3:
                    cuts = (e1.charge * e2.charge) < 0 #Opposite charge
                    cuts = cuts and (abs(e1.eta + e1.deltaEtaSC) >= 1.566 or abs( e1.eta + e1.deltaEtaSC) < 1.444) #Fiducial
                    cuts = cuts and (abs(e2.eta + e2.deltaEtaSC) >= 1.566 or abs(e2.eta + e2.deltaEtaSC) < 1.444)
                    #cuts = cuts and (e1.pt >= 20.0 and abs(e1.eta + e1.deltaEtaSC) < 2.5 and e1.mvaNoIso_WP80 ) #ID (basic)
                    #cuts = cuts and (e2.pt >= 20.0 and abs(e2.eta + e2.deltaEtaSC) < 2.5 and e2.mvaNoIso_WP80) #ID (basic)

                if cuts:
                    Z_nEE += 1
                    tempM = (e1.p4() + e2.p4()).M()
                    #print("Found a Z->ee candidate with mass = " + str(tempM) + " : current Z_mass = " + str(Z_mass))
                    if abs(tempM - 91.18) < abs(Z_mass - 91.18) and tempM >= 60.0 and tempM <= 120.0: 
                        Z_dm = 1
                        Z_mass = tempM
                        Z_pt = (e1.p4()+e2.p4()).Pt()
                        Z_eta = (e1.p4()+e2.p4()).Eta()
                        Z_phi = (e1.p4()+e2.p4()).Phi()
                        Z_d1Idx = e1Idx if e1.pt >= e2.pt else e2Idx  #daughter 1 is higher pT e
                        Z_d2Idx = e2Idx if e1.pt >= e2.pt else e1Idx 
                        Z_dauDR = e1.DeltaR(e2)
                        if self.era == 2:
                            Z_dauID = (e1.pt >= 20.0 and abs(e1.eta + e1.deltaEtaSC) < 2.5 and e1.mvaFall17V2noIso_WP80)
                            Z_dauID = Z_dauID and (e2.pt >= 20.0 and abs(e2.eta + e2.deltaEtaSC) < 2.5 and e2.mvaFall17V2noIso_WP80)
                        elif self.era == 3:
                            Z_dauID = (e1.pt >= 20.0 and abs(e1.eta + e1.deltaEtaSC) < 2.5 and e1.mvaNoIso_WP80 )
                            Z_dauID = Z_dauID and (e2.pt >= 20.0 and abs(e2.eta + e2.deltaEtaSC) < 2.5 and e2.mvaNoIso_WP80)
        
        muons = Collection(event, "Muon")
        for mu1Idx, mu1 in enumerate(muons):
            for mu2Idx in range(mu1Idx+1, event.nMuon):
                mu2 = muons[mu2Idx]
                
                if (mu1.charge * mu2.charge < 0): #Opposite charge
                    if (mu1.pt >= 15.0 and abs(mu1.eta) < 2.4) and (mu2.pt >= 15.0 and abs(mu2.eta) < 2.4): 
                        Z_nMuMu += 1
                        tempM = (mu1.p4() + mu2.p4()).M()
                        #print("Found a Z->mumu candidate with mass = " + str(tempM) + " : current Z_mass = " + str(Z_mass))
                        if abs(tempM - 91.18) < abs(Z_mass - 91.18) and tempM >= 60.0 and tempM <= 120.0: 
                            Z_dm = 2
                            Z_mass = tempM
                            Z_pt = (mu1.p4()+mu2.p4()).Pt()
                            Z_eta = (mu1.p4()+mu2.p4()).Eta()
                            Z_phi = (mu1.p4()+mu2.p4()).Phi()
                            Z_d1Idx = mu1Idx if mu1.pt >= mu2.pt else mu2Idx  #daughter 1 is higher pT mu
                            Z_d2Idx = mu2Idx if mu1.pt >= mu2.pt else mu1Idx
                            Z_dauDR = mu1.DeltaR(mu2)
                            Z_dauID = mu1.mediumId and mu2.mediumId

        
        if Z_dm < 0:
            ak8Jets = Collection(event, "FatJet")
            ak4Jets = Collection(event, "Jet")
            
        
            for jetIdx, jet in enumerate(ak8Jets):
                if abs(jet.eta) < 2.5 and jet.pt > 100:
                    self.h_ak8Mass.Fill(jet.mass)
                jetID = jet.jetId == 6 #2 = pass tight ID but fail tight lepton veto, 6 = pass both
                jetID = jetID and (jet.mass >= 60.0 and jet.mass <= 120.0)
                jetID = jetID and (abs(jet.eta) < 2.5 and jet.pt > 100)
                if self.era == 2:
                    jetID = jetID and jet.particleNet_ZvsQCD > 0.9
                elif self.era == 3:
                    jetID = jetID and jet.particleNetWithMass_ZvsQCD > 0.9

                if jetID:
                    self.h_ak8MassCuts.Fill(jet.mass)
                    if abs(jet.mass - 91.18) < abs(Z_mass - 91.18):
                        Z_mass = jet.mass 
                        Z_jetIdxAK8 = jetIdx
                        Z_jetR = 8
                        Z_dm = 0
                        Z_pt = jet.pt
                        Z_eta = jet.eta
                        Z_phi = jet.phi
        
            bestAK4Mass = -999.99
            for jetIdx, jet in enumerate(ak4Jets):
                if abs(jet.eta) < 2.5 and jet.pt > 100:
                    self.h_ak4Mass.Fill(jet.mass)
                jetID = jet.jetId == 6 #2 = pass tight ID but fail tight lepton veto, 6 = pass both
                jetID = jetID and (jet.mass >= 60.0 and jet.mass <= 120.0)
                jetID = jetID and (abs(jet.eta) < 2.5 and jet.pt > 250)

                if jetID:
                    self.h_ak4MassCuts.Fill(jet.mass)
                    if abs(jet.mass - 91.18) < abs(bestAK4Mass - 91.18):
                        Z_jetIdxAK4 = jetIdx
                        bestAK4Mass = jet.mass 

            #IF AK4 jet was a better match than AK8, use AK4
            if abs(bestAK4Mass - 91.18) < abs(Z_mass - 91.18):
                theJet = ak4Jets[Z_jetIdxAK4]
                Z_dm = 0
                Z_jetR = 4
                Z_mass = jet.mass
                Z_pt = theJet.pt
                Z_eta = theJet.eta
                Z_phi = theJet.phi

        Z_isCand = Z_dm == 0 or Z_dm == 1 or Z_dm==2 #Z->jets, ee, mumu
        Z_isCand = Z_isCand and (Z_mass > 81 and Z_mass < 101) #Mass range

        self.out.fillBranch("Z_dm", Z_dm)
        self.out.fillBranch("Z_d1Idx", Z_d1Idx)
        self.out.fillBranch("Z_d2Idx", Z_d2Idx)
        self.out.fillBranch("Z_mass", Z_mass)
        self.out.fillBranch("Z_pt", Z_pt)
        self.out.fillBranch("Z_eta", Z_eta)
        self.out.fillBranch("Z_phi", Z_phi)
        self.out.fillBranch("Z_dauDR", Z_dauDR)
        self.out.fillBranch("Z_dauID", Z_dauID)
        self.out.fillBranch("Z_jetIdxAK8", Z_jetIdxAK8)
        self.out.fillBranch("Z_jetIdxAK4", Z_jetIdxAK4)
        self.out.fillBranch("Z_jetR",Z_jetR)
        self.out.fillBranch("Z_nEE", Z_nEE)
        self.out.fillBranch("Z_nMuMu", Z_nMuMu)
        self.out.fillBranch("Z_isCand", Z_isCand)

        return True
    
    # -----------------------------------------------------------------------------------------------------------------------------

zProducerConstr = lambda era: ZProducer(era = era)

#from PhysicsTools.NanoAODTools.postprocessing.modules.GenProducerZTau import genProducerZTauConstr

#files = ["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/TauZ/taustarToTauZ_m500_2018.root"]
#p = PostProcessor(".", files, cut="1>0", branchsel=None, postfix="", modules=[genProducerZTauConstr(), zProducerConstr()] )
#p.run()
