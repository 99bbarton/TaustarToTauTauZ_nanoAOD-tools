from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from ROOT import TLorentzVector
import os

class VisualizationDumper(Module):

    def __init__(self, hasPFInfo=False):
        self.evtNum = -1
        self.hasPFInfo = hasPFInfo
    # ----------------------------------------------------------------------------------------------------------------------------- #

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)
    
    # ----------------------------------------------------------------------------------------------------------------------------- #

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

    def analyze(self, event):
        self.evtNum += 1
        evtNumStr = str(self.evtNum)

        if not (event.Gen_isCand and event.Z_isCand):
            return False
        
        #Print the Z and any of it's decay products
        theZ = TLorentzVector()
        theZ.SetPtEtaPhiM(event.Z_pt, event.Z_eta, event.Z_phi, event.Z_mass)
        zStr = evtNumStr+",Z," + str(theZ.Px()) + "," + str(theZ.Py()) + "," + str(theZ.Pz()) + "," + str(theZ.M()) +  ",0,0,0,"
        if event.Z_dm==0: #Z->FatJet
            if event.Z_sJIdx1 < 0 or event.Z_sJIdx2 < 0:
                return False
            
            zStr += "0.8"
            zStr =zStr.replace("Z", "Z(AK8)")
            print(zStr)

            subJets = Collection(event, "SubJet")
            subJet1 = subJets[event.Z_sJIdx1]
            print(evtNumStr+",Z_AK4SJ," + str(subJet1.Px()) + "," + str(subJet1.Py()) + "," + str(subJet1.Pz()) + "," + str(subJet1.M()) + ",0,0,0,0.4")
            subJet2 = subJets[event.Z_sJIdx2]
            print(evtNumStr+",Z_AK4SJ," + str(subJet2.Px()) + "," + str(subJet2.Py()) + "," + str(subJet2.Pz()) + "," + str(subJet2.M()) + ",0,0,0,0.4")

            if self.hasPFInfo:
                fjPFCands = Collection(event, "FatJetPFCands")
                pfCands = Collection(event, "PFCands")
            
                for fjPFCand in fjPFCands:
                    if fjPFCand.jetIdx == event.Z_jetIdxAK8:
                        pfCand = pfCands[fjPFCand.pFCandsIdx].p4()
                        print(evtNumStr, ",Z_pfc("+str(pfCand.pdgId)+"),"+str(pfCand.Px())+","+str(pfCand.Py())+","+str(pfCand.Pz())+","+str(pfCand.M())+","+str(pfCand.d0)+","+str(pfCand.dz)+",1,0")
        elif event.Z_dm == 1:
            zStr += "1.0"
            zStr =zStr.replace("Z", "Z(ee)")
            print(zStr)

            electrons = Collection(event, "Electron")
            z_d1 = electrons[event.Z_d1Idx]
            z_d2 = electrons[event.Z_d2Idx]
            print(evtNumStr+",e,"+str(z_d1.Px())+","+str(z_d1.Py())+","+str(z_d1.Pz())+",0.000511,"+str(z_d1.dxy)+","+str(z_d1.dz)+",1,0")
            print(evtNumStr+",e,"+str(z_d2.Px())+","+str(z_d2.Py())+","+str(z_d2.Pz())+",0.000511,"+str(z_d1.dxy)+","+str(z_d1.dz)+",1,0")
        elif event.Z_dm == 2:
            zStr += "1.0"
            zStr =zStr.replace("Z", "Z(mumu)")
            print(zStr)

            muons = Collection(event, "Muon")
            z_d1 = muons[event.Z_d1Idx]
            z_d2 = muons[event.Z_d2Idx]
            print(evtNumStr+",mu,"+str(z_d1.Px())+","+str(z_d1.Py())+","+str(z_d1.Pz())+",0.106,"+str(z_d1.dxy)+","+str(z_d1.dz)+",1,0")
            print(evtNumStr+",mu,"+str(z_d2.Px())+","+str(z_d2.Py())+","+str(z_d2.Pz())+",0.106,"+str(z_d1.dxy)+","+str(z_d1.dz)+",1,0")
        

        taus = Collection(event, "Tau")

        #ETau
        if event.ETau_isCand:
            electrons = Collection(event, "Electron")
            theEl = electrons[event.ETau_eIdx]
            print(evtNumStr+",tau(e),"+str(theEl.Px())+","+str(theEl.Py())+","+str(theEl.Pz())+",0.000511,"+str(theEl.dxy)+","+str(theEl.dz)+",1,0")
            
            theTau = taus[event.ETau_tauIdx]
            print(evtNumStr+",tau(h),"+str(theTau.Px())+","+str(theTau.Py())+","+str(theTau.Pz())+",1.78,"+str(theTau.dxy)+","+str(theTau.dz)+",0,0.5")
        #MuTau
        if event.MuTau_isCand:
            muons = Collection(event, "Muon")
            theMu = electrons[event.MuTau_muIdx]
            print(evtNumStr+",tau(mu),"+str(theMu.Px())+","+str(theMu.Py())+","+str(theMu.Pz())+",0.106,"+str(theMu.dxy)+","+str(theMu.dz)+",1,0")
            
            theTau = taus[event.MuTau_tauIdx]
            print(evtNumStr+",tau(h),"+str(theTau.Px())+","+str(theTau.Py())+","+str(theTau.Pz())+",1.78,"+str(theTau.dxy)+","+str(theTau.dz)+",0,0.5")
        #TauTau
        if event.TauTau_isCand:
            tau1 = taus[event.TauTau_tau1Idx]
            print(evtNumStr+",tau(h),"+str(tau1.Px())+","+str(tau1.Py())+","+str(tau1.Pz())+",1.78,"+str(tau1.dxy)+","+str(tau1.dz)+",0,0.5")
            tau2 = taus[event.TauTau_tau2Idx]
            print(evtNumStr+",tau(h),"+str(tau2.Px())+","+str(tau2.Py())+","+str(tau2.Pz())+",1.78,"+str(tau2.dxy)+","+str(tau2.dz)+",0,0.5")

        return False #Don't need an output file
    
visualizationDumperConstr = lambda hasPFInfo: VisualizationDumper(hasPFInfo=hasPFInfo)

from PhysicsTools.NanoAODTools.postprocessing.modules.VisualizationDumper import visualizationDumperConstr

files = [os.environ["SIG_R3"] + "/taustarToTauZ_m3000_2022.root"]
p = PostProcessor(".", files, cut="Gen_isCand&&Z_isCand&&(ETau_isCand||MuTau_isCand||TauTau_isCand)", branchsel=None, postfix="", modules=[visualizationDumperConstr(True)] )
p.run()