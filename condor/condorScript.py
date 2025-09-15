
import os
import sys


from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.GenProducer import genProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TrigProducer import trigProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ZProducer import zProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.BoostProducer import boostProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ZJetReclusterProducer import zJetReclusterProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ETauProducer import eTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.MuTauProducer import muTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TauTauProducer import tauTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ObjCounter import objCounterConstr

year = sys.argv[1]
era = 0
if year in ["2016", "2016post", "2017", "2018"]:
    era = 2
elif year in ["2022", "2022post", "2023", "2023post"]:
    era = 3


inputFiles = []
files = os.listdir("./")
for fN, inpFile in enumerate(files):
    print(inpFile)
    if inpFile.endswith(".root"):
        inputFiles.append(inpFile)

print("\n\ncondorScript found the following input ROOT files:", inputFiles)
        
modules = [genProducerConstr(era), trigProducerConstr(year), zProducerConstr(year), zJetReclusterProducerConstr(), eTauProducerConstr(year), muTauProducerConstr(year), tauTauProducerConstr(year), objCounterConstr()]

if era ==2:
    cut_etau = "(Sum$(TMath::Abs(Electron_eta)<2.5 && Electron_pt>=24. && Electron_mvaFall17V2Iso_WP80 && (TMath::Abs(Electron_eta+El\
ectron_deltaEtaSC)>=1.566||TMath::Abs(Electron_eta+Electron_deltaEtaSC)<1.444))>0 && Sum$(Tau_pt>20. && TMath::Abs(Tau_eta)<2.3 && TM\
ath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7))"
    cut_mutau = "(Sum$(TMath::Abs(Muon_eta)<2.4 && Muon_pt>=20. && Muon_mediumId)>0 && Sum$(Tau_pt>20. && TMath::Abs(Tau_eta)<2.3 && \
TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7))"
    cut_tautau = "(Sum$(Tau_pt>20. && TMath::Abs(Tau_eta)<2.3 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Ta\
u_decayMode!=7))"
    cut_Z = "(nFatJet > 0 || nElectron >= 2 || nMuon >= 2)"
elif era == 3:
    cut_etau   = "(Sum$(TMath::Abs(Electron_eta)<2.5 && Electron_pt>=24. && Electron_mvaIso_WP80 && (TMath::Abs(Electron_eta+Electron_deltaEtaSC)>=1.566||TMath::Abs(Electron_eta+Electron_deltaEtaSC)<1.444))>0 && Sum$(Tau_pt>=30. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7 && (Tau_idDeepTau2018v2p5VSjet>=4) && (Tau_idDeepTau2018v2p5VSmu>=4) && (Tau_idDeepTau2018v2p5VSe>=2))>0)"
    cut_mutau  = "(Sum$(TMath::Abs(Muon_eta)<2.4 && Muon_pt>=20. && Muon_mediumId)>0 && Sum$(Tau_pt>=27. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7 && (Tau_idDeepTau2018v2p5VSjet>=4) && (Tau_idDeepTau2018v2p5VSmu>=4) && (Tau_idDeepTau2018v2p5VSe>=2))>0)"
    cut_tautau = "(Sum$(Tau_pt>=35. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7 && (Tau_idDeepTau2018v2p5VSjet>=4) && (Tau_idDeepTau2018v2p5VSmu>=4) && (Tau_idDeepTau2018v2p5VSe>=2))>=2)"
    cut_Z = "(nFatJet > 0 || nElectron >= 2 || nMuon >= 2)"

preSelection = "(" + cut_Z + "&&(" + cut_tautau + "||" + cut_etau + "||" + cut_mutau + "))"


p = PostProcessor("outputs", inputFiles, cut=preSelection, branchsel="../keep_and_drop.txt", postfix="output", modules=modules, fwkJobReport=True)
p.run()
