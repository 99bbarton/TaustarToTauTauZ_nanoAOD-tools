from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.ZProducer import zProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.GenProducerZTau import genProducerZTauConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.GenProducerWNu import genProducerWNuConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TauProducer import tauProducerConstr


baseDir = "root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/"
#masses = ["250","500","750","1000","1500","2000","2500","3000","3500","4000","4500","5000"]
#masses = ["4000", "4500", "5000"]
masses = ["500"]
decay = "TauZ"
year = "2018"

for mass in masses:
    files = [baseDir + decay +"/taustarTo" + decay + "_m" + mass + "_" + year + ".root"]
    if decay == "TauZ":
        modules = [genProducerZTauConstr(), zProducerConstr(), tauProducerConstr()]
    elif decay == "WNu":
        modules = [genProducerWNuConstr(), tauProducerConstr()]
    else: 
        modules = []
    p = PostProcessor("data", files, cut="1>0", branchsel=None, postfix="", modules=modules )
    p.run()
