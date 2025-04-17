
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

modules = [genProducerConstr(era), trigProducerConstr(year), zProducerConstr(era), zJetReclusterProducerConstr(), eTauProducerConstr(era), muTauProducerConstr(era), tauTauProducerConstr(era)]

preSelection = "1>0"

p = PostProcessor("./", inputFiles, cut=preSelection, branchsel="../crab/keep_and_drop.txt", postfix="", modules=modules, histFileName="hists.root", histDirName="Hists")
p.run()
