#!/usr/bin/env python3
import os
import sys
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import *
# this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.utils.crabhelper import inputFiles, runsAndLumis

from PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule import exampleModuleConstr

from PhysicsTools.NanoAODTools.postprocessing.modules.GenProducer import genProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TrigProducer import trigProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ZProducer import zProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ETauProducer import eTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.MuTauProducer import muTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TauTauProducer import tauTauProducerConstr


if len(sys.argv) == 3:
    year = sys.argv[2]
    if year in ["2022", "2022post", "2023", "2023post"]:
        era = 3
    elif year in ["2016", "2016post", "2017", "2018"]:
        era = 2
    else:
        print("ERROR: Year not recognized! Year was " + year)
        exit(1)
else:
    print("ERROR: Number of arguments was not 3! syst.argv was : " + str(sys.argv))
    exit(1)

outDir = "./"
preSelection = "1>0"
modules = [genProducerConstr(era), trigProducerConstr(year), zProducerConstr(era), eTauProducerConstr(era), muTauProducerConstr(era), tauTauProducerConstr(era)]
#modules = [exampleModuleConstr()]
############################## TODO add era argument to constructors for z and tau, in .sh and config maker


#files=["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/TauZ/taustarToTauZ_m250_2018.root"]

#NB: Despite not being in most of the supposedly working examples and defaulting to False, both fwkJObReport and provenance must be True for jobs to succeed fully.
p = PostProcessor(outDir, inputFiles(), cut=preSelection, modules=modules, haddFileName="tree.root", fwkJobReport=True, provenance=True, histFileName="hists.root", histDirName="Hists" )
p.run()

print("DONE")
