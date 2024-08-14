#!/usr/bin/env python
import os
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import *
# this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.utils.crabhelper import inputFiles, runsAndLumis

from PhysicsTools.NanoAODTools.postprocessing.modules.GenProducerZTau import genProducerZTauConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ZProducer import zProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TauProducer import tauProducerConstr

outDir = "./"
preSelection = "1>0"
modules = [genProducerZTauConstr(), zProducerConstr(), tauProducerConstr()]

#files=["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/taustarToTauZ_m250_2018.root"]
p = PostProcessor(outDir, inputFiles(), cut=preSelection, modules=modules)
p.run()

print("DONE")
