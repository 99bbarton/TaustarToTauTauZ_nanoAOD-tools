from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.GenProducer import genProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TrigProducer import trigProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ZProducer import zProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.BoostProducer import boostProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ZJetReclusterProducer import zJetReclusterProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ETauProducer import eTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.MuTauProducer import muTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TauTauProducer import tauTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.FinalProducer import finalProducerConstr


baseDir = "root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/"
testFiles = ["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/SigPFNano/2023/taustarToTauZ_m3000_2023.root"]
#testFiles = ["./data/taustarToTauZ_m3000_2023.root"]

#testFile = "root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m250_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/2810000/16d5e0fd-d03b-4131-8d1f-3796807217a2.root"
#testFiles = ["root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m3000_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/2820000/be403352-9d04-41ac-8099-b00fc6304bec.root", "root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m3000_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/2820000/242c32f7-a63d-4aaf-99e3-e7ef93be748b.root", "root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m3000_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/120000/fb2e39e2-fc8e-43be-be61-4080b89ffd88.root"]
#testFiles = ["root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m3000_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/120000/fb2e39e2-fc8e-43be-be61-4080b89ffd88.root"]
#masses = ["250","500","750","1000","1500","2000","2500","3000","3500","4000","4500","5000"]
#masses = ["4000", "4500", "5000"]

masses = ["3000"]
decay = "TauZ"
year = "2023"
era = 3

if era == 3:
    cut_etau   = "(Sum$(TMath::Abs(Electron_eta)<2.5 && Electron_pt>=24. && Electron_mvaIso_WP80 && (TMath::Abs(Electron_eta+Electron_deltaEtaSC)>=1.566||TMath::Abs(Electron_eta+Electron_deltaEtaSC)<1.444))>0 && Sum$(Tau_pt>=30. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7 && (Tau_idDeepTau2018v2p5VSjet>=4) && (Tau_idDeepTau2018v2p5VSmu>=4) && (Tau_idDeepTau2018v2p5VSe>=2))>0)"
    cut_mutau  = "(Sum$(TMath::Abs(Muon_eta)<2.4 && Muon_pt>=20. && Muon_mediumId)>0 && Sum$(Tau_pt>=27. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7 && (Tau_idDeepTau2018v2p5VSjet>=4) && (Tau_idDeepTau2018v2p5VSmu>=4) && (Tau_idDeepTau2018v2p5VSe>=2))>0)"
    cut_tautau = "(Sum$(Tau_pt>=35. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7 && (Tau_idDeepTau2018v2p5VSjet>=4) &&\
 (Tau_idDeepTau2018v2p5VSmu>=4) && (Tau_idDeepTau2018v2p5VSe>=2))>=2)"
    cut_Z = "(nFatJet > 0 || nElectron >= 2 || nMuon >= 2)"

preSelection = "(" + cut_Z + "&&(" + cut_tautau + "||" + cut_etau + "||" + cut_mutau + "))"



for mass in masses:
    #files = [baseDir + decay +"/taustarTo" + decay + "_m" + mass + "_" + year + ".root"]
    #files = [testFile]
    if decay == "TauZ":
        modules = [genProducerConstr(era), trigProducerConstr(year), zProducerConstr(era), zJetReclusterProducerConstr(),  eTauProducerConstr(era), muTauProducerConstr(era), tauTauProducerConstr(era)]#, finalProducerConstr()]
        #modules = [tauProducerConstr(era)]
        #modules = [finalProducerConstr()]
        #modules = [genProducerConstr(era)]
        #modules = [genProducerConstr(era), zProducerConstr(era), finalProducerConstr()]#zJetReclusterProducerConstr()]
    #elif decay == "WNu":
    #    modules = [genProducerWNuConstr(), tauProducerConstr()]
    else: 
        modules = []
    p = PostProcessor("data", testFiles, cut=preSelection, branchsel="crab/keep_and_drop.txt", postfix="", modules=modules, histFileName="hists.root", histDirName="Hists")
    p.run()
