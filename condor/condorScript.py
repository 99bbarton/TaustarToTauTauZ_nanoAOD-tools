
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


year = sys.argv[1]
era = 0
if year in ["2022", "2022post", "2023", "2023post"]:
    era = 3

inputFiles = os.listdir("./")
for fN, inpFile in enumerate(inputFiles):
    if not inpFile.endswith(".root"):
        del inputFiles[fN]

print("\n\ncondorScript found the following input ROOT files:", inputFiles)
        
modules = [genProducerConstr(era), trigProducerConstr(year), zProducerConstr(era), zJetReclusterProducerConstr(), eTauProducerConstr(era), muTauProducerConstr(era), tauTauProducerConstr(era)]


if era == 3:
    cut_etau   = "(Sum$(TMath::Abs(Electron_eta)<2.5 && Electron_pt>=24. && Electron_mvaIso_WP80 && (TMath::Abs(Electron_eta+Electron_deltaEtaSC)>=1.566||TMath::Abs(Electron_eta+Electron_\
deltaEtaSC)<1.444))>0 && Sum$(Tau_pt>=30. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7 && (Tau_idDeepTau2018v2p5VSjet>=\
4) && (Tau_idDeepTau2018v2p5VSmu>=4) && (Tau_idDeepTau2018v2p5VSe>=2))>0)"
    cut_mutau  = "(Sum$(TMath::Abs(Muon_eta)<2.4 && Muon_pt>=20. && Muon_mediumId)>0 && Sum$(Tau_pt>=27. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_de\
cayMode!=6 && Tau_decayMode!=7 && (Tau_idDeepTau2018v2p5VSjet>=4) && (Tau_idDeepTau2018v2p5VSmu>=4) && (Tau_idDeepTau2018v2p5VSe>=2))>0)"
    cut_tautau = "(Sum$(Tau_pt>=35. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7 && (Tau_idDeepTau2018v2p5VSjet>=4) &&\
 (Tau_idDeepTau2018v2p5VSmu>=4) && (Tau_idDeepTau2018v2p5VSe>=2))>=2)"
    cut_Z = "(nFatJet > 0 || nElectron >= 2 || nMuon >= 2)"

preSelection = "(" + cut_Z + "&&(" + cut_tautau + "||" + cut_etau + "||" + cut_mutau + "))"


p = PostProcessor("outputs/", inputFiles, cut=preSelection, branchsel="../crab/keep_and_drop.txt", postfix="", modules=modules, histFileName="hists.root", histDirName="Hists")
p.run()
