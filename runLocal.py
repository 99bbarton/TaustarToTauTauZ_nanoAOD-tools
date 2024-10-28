from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.GenProducer import genProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TrigProducer import trigProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ZProducer import zProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ETauProducer import eTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.MuTauProducer import muTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TauTauProducer import tauTauProducerConstr



baseDir = "root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/"
#testFile = "root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m250_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/2810000/16d5e0fd-d03b-4131-8d1f-3796807217a2.root"
testFiles = ["root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m3000_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/2820000/be403352-9d04-41ac-8099-b00fc6304bec.root", "root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m3000_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/2820000/242c32f7-a63d-4aaf-99e3-e7ef93be748b.root", "root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m3000_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/120000/fb2e39e2-fc8e-43be-be61-4080b89ffd88.root"]
#testFiles = ["root://cmsxrootd.fnal.gov//store/mc/Run3Summer23NanoAODv12/TaustarToTauZ_m3000_TuneCP5_13p6TeV_pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/120000/fb2e39e2-fc8e-43be-be61-4080b89ffd88.root"]
#masses = ["250","500","750","1000","1500","2000","2500","3000","3500","4000","4500","5000"]
#masses = ["4000", "4500", "5000"]

masses = ["3000"]
decay = "TauZ"
year = "2023"
era = 3

for mass in masses:
    #files = [baseDir + decay +"/taustarTo" + decay + "_m" + mass + "_" + year + ".root"]
    #files = [testFile]
    if decay == "TauZ":
        #modules = [genProducerConstr(era), trigProducerConstr(year), zProducerConstr(era), eTauProducerConstr(era), muTauProducerConstr(era), tauTauProducerConstr(era)]
        #modules = [tauProducerConstr(era)]
        modules = [genProducerConstr(era), zProducerConstr(era)]
    #elif decay == "WNu":
    #    modules = [genProducerWNuConstr(), tauProducerConstr()]
    else: 
        modules = []
    p = PostProcessor("data", testFiles, cut="1>0", branchsel=None, postfix="", modules=modules, histFileName="hists.root", histDirName="Hists")
    p.run()
