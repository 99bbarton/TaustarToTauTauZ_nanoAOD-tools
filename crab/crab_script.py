#!/usr/bin/env python3
import os
import sys
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import *
# this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.utils.crabhelper import inputFiles, runsAndLumis

#from PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule import exampleModuleConstr

from PhysicsTools.NanoAODTools.postprocessing.modules.GenProducer import genProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TrigProducer import trigProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ZProducer import zProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.BoostProducer import boostProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ZJetReclusterProducer import zJetReclusterProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.ETauProducer import eTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.MuTauProducer import muTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.TauTauProducer import tauTauProducerConstr
from PhysicsTools.NanoAODTools.postprocessing.modules.FinalProducer import finalProducerConstr


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

cut_etau   = "(Sum$(TMath::Abs(Electron_eta)<2.5 && Electron_pt>=20. && (Electron_mvaFall17V2Iso_WP90||Electron_mvaFall17V2noIso_WP90) && (TMath::Abs(Electron_eta+Electron_deltaEtaSC)>=1.566||TMath::Abs(Electron_eta+Electron_deltaEtaSC)<1.444))>0 && Sum$(Tau_pt>=30. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7 && (4&Tau_idDeepTau2017v2p1VSjet) && (4&Tau_idDeepTau2017v2p1VSmu) && (2&Tau_idDeepTau2017v2p1VSe))>0)"
cut_mutau  = "(Sum$(TMath::Abs(Muon_eta)<2.4 && Muon_pt>=20. && Muon_mediumId)>0 && Sum$(Tau_pt>=27. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7 && (4&Tau_idDeepTau2017v2p1VSjet) && (4&Tau_idDeepTau2017v2p1VSmu) && (2&Tau_idDeepTau2017v2p1VSe))>0)"
cut_tautau = "(Sum$(Tau_pt>=35. && TMath::Abs(Tau_eta)<2.1 && TMath::Abs(Tau_dz)<0.2 && Tau_decayMode!=5 && Tau_decayMode!=6 && Tau_decayMode!=7 && (4&Tau_idDeepTau2017v2p1VSjet) && (4&Tau_idDeepTau2017v2p1VSmu) && (2&Tau_idDeepTau2017v2p1VSe))>=2)"
cut_Z = "(nFatJet > 0 || nElectron >= 2 || nMuon >= 2)"

preSelection = "(" + cut_Z + "&&(" + cut_tautau + "||" + cut_etau + "||" + cut_mutau + "))"

modules = [genProducerConstr(era), trigProducerConstr(year), zProducerConstr(era), boostProducerConstr(), zJetReclusterProducerConstr(), eTauProducerConstr(era), muTauProducerConstr(era), tauTauProducerConstr(era), finalProducerConstr()]
#modules = [exampleModuleConstr()]
############################## TODO add era argument to constructors for z and tau, in .sh and config maker


#files=["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/TauZ/taustarToTauZ_m250_2018.root"]

#NB: Despite not being in most of the supposedly working examples and defaulting to False, both fwkJObReport and provenance must be True for jobs to succeed fully.
p = PostProcessor(outDir, inputFiles(), cut=preSelection, modules=modules, haddFileName="tree.root", fwkJobReport=True, provenance=True, histFileName="hists.root", histDirName="Hists", outputbranchsel="keep_and_drop.txt" )
p.run()

print("DONE")
