#!/usr/bin/env python
import os
import sys
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import *
# this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.utils.crabhelper import inputFiles, runsAndLumis

from PhysicsTools.NanoAODTools.postprocessing.modules.GenProducerZTau import genProducerZTauConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ZProducer import zProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TauProducer import tauProducerConstr


#if len(sys.argv) == 3:
#    era = sys.argv[2]
#    era = int(era)
#    if era != 2 and era != 3:
#        print("Era can only be 2 or 3")
#        exit(1)
#else:
#    print("Expected era argument!")
#    exit(1)

era = 3
outDir = "./"
preSelection = "1>0"
modules = [genProducerZTauConstr(), zProducerConstr(era), tauProducerConstr(era)]
############################## TODO add era argument to constructors for z and tau, in .sh and config maker


#files=["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/TauZ/taustarToTauZ_m250_2018.root"]

#NB: Despite not being in most of the supposedly working examples and defaulting to False, both fwkJObReport and provenance must be True for jobs to succeed fully.
p = PostProcessor(outDir, inputFiles(), cut=preSelection, modules=modules, haddFileName="tree.root", fwkJobReport=True, provenance=True)
p.run()

print("DONE")
