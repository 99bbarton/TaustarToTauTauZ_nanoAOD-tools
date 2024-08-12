#!/usr/bin/env python
import os
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import *
# this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles, runsAndLumis

from PhysicsTools.NanoAODTools.postprocessing.modules.GenProducerZTau import genProducerZTauConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ZProducer import zProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TauProducer import tauProducerConstr

outDir = "./"
preSelection = "1>0"
modules = [genProducerZTauConstr(), zProducerConstr(), tauProducerConstr()]

p = PostProcessor(outDir, inputFiles(), preSelection, modules=modules, provenance=True, fwkJobReport=True,)
p.run()

print("DONE")