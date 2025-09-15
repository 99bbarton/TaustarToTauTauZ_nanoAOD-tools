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


testFiles = ["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/SigPFNano/2023/taustarToTauZ_m3000_2023.root"]
#testFiles = ["./data/taustarToTauZ_m3000_2023.root"]

#testFiles = "root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m250_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/2810000/16d5e0fd-d03b-4131-8d1f-3796807217a2.root"
#testFiles = ["root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m3000_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/2820000/be403352-9d04-41ac-8099-b00fc6304bec.root", "root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m3000_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/2820000/242c32f7-a63d-4aaf-99e3-e7ef93be748b.root", "root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m3000_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/120000/fb2e39e2-fc8e-43be-be61-4080b89ffd88.root"]
#testFiles = ["root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m3000_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/120000/fb2e39e2-fc8e-43be-be61-4080b89ffd88.root"]

#RUN2
testFiles = ["root://cmsxrootd.fnal.gov//store/mc/RunIISummer20UL18NanoAODv9/TaustarToTauZ_m3000_TuneCP5_13TeV_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2830000/79B9792D-6D5F-2D4E-BDDC-697739B52CF7.root"]

year = "2018"

if year in ["2016", "2016post", "2017", "2018"]:
    era = 2
elif year in ["2022", "2022post", "2023", "2023post"]:
    era = 3

if era == 2:
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
    cut_tautau = "(Sum$(Tau_pt>=35. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7 && (Tau_idDeepTau2018v2p5VSjet>=4) &&\
 (Tau_idDeepTau2018v2p5VSmu>=4) && (Tau_idDeepTau2018v2p5VSe>=2))>=2)"
    cut_Z = "(nFatJet > 0 || nElectron >= 2 || nMuon >= 2)"

preSelection = "(" + cut_Z + "&&(" + cut_tautau + "||" + cut_etau + "||" + cut_mutau + "))"

modules = [genProducerConstr(era), trigProducerConstr(year), zProducerConstr(year), zJetReclusterProducerConstr(),  eTauProducerConstr(year), muTauProducerConstr(year), tauTauProducerConstr(year), objCounterConstr()]
#modules = [genProducerConstr(era), trigProducerConstr(year), zProducerConstr(year)]

p = PostProcessor("data", testFiles, cut=preSelection, branchsel="keep_and_drop.txt", postfix="", modules=modules)#, histFileName="hists.root", histDirName="Hists")
p.run()
