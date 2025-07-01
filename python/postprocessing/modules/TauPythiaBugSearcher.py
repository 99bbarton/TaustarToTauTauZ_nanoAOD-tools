#Analyzer to produce histograms of the visible energy fraction of DM==0 (tau->pi + nu) in + vs - charged tau decays
#Needed to examine if samples are affected by the bug discussed here:
# https://indico.cern.ch/event/1549765/contributions/6526658/attachments/3084210/5459973/TauPOG-JuneCMSWeek%20(3).pdf 
# https://indico.cern.ch/event/1495529/contributions/6535498/attachments/3084314/5460656/TauModellingIssues_Lucas.pdf


from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.utils.GenTools import getDecayChain

from ROOT import TH1F, TLorentzVector

# ----------------------------------------------------------------------------------------------------------------------------- #

class TauPythiaBugSearcher(Module):

    def __init__(self):
        self.writeHistFile = True

    # ----------------------------------------------------------------------------------------------------------------------------- #

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)
        self.h_visEFrac_plus = TH1F("h_visEFrac_plus", "Visible Energy Fraction in DM=0 Tau Decays;Energy Fraction;Events", 20, 0, 100)
        self.addObject(self.h_visEFrac_plus)
        self.h_visEFrac_neg = TH1F("h_visEFrac_neg", "Visible Energy Fraction in DM=0 Tau Decays;Energy Fraction;Events", 20, 0, 100)
        self.addObject(self.h_visEFrac_neg)   

    # ----------------------------------------------------------------------------------------------------------------------------- #

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

    def analyze(self, event):

        genParts = Collection(event, "GenPart")

        for gpIdx, genPart in enumerate(genParts):
            if not genPart.statusflag('isLastCopy'):
                continue

            if abs(genPart.pdgId) == 15:
                decChain = getDecayChain(gpIdx, genParts)
                if len(decChain) == 2:
                    tauP4 = TLorentzVector()
                    tauP4.SetPtEtaPhiM(genPart.pt, genPart.eta, genPart.phi, 1.776) 
                    
                    if abs(genParts[decChain[0]].pdgId) == 211:
                        pion = genParts[decChain[0]]
                    elif abs(genParts[decChain[1]].pdgId) == 211:
                        pion = genParts[decChain[1]]
                    else:
                        print("Warning: Neither tau decay product was a charged pion!")
                        continue

                    piP4 = TLorentzVector()
                    piP4.SetPtEtaPhiM(pion.pt, pion.eta, pion.phi, 0.1395)
                    
                    if pion.pdgId > 0:
                        self.h_visEFrac_plus.Fill(piP4.E() / tauP4.E())
                    else:
                        self.h_visEFrac_neg.Fill(piP4.E() / tauP4.E())
        
        return False #We only care about the hist file, not the output tree
    
        # ----------------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------------------------------------- #

tauPythiaBugSearcherConstr = lambda : TauPythiaBugSearcher()

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.TauPythiaBugSearcher import tauPythiaBugSearcherConstr
import os



masses = ["250","500","750","1000","1500","2000","2500","3000","3500","4000","4500","5000"]
years = ["2022", "2022post", "2023", "2023post"]
files = []
for year in years:
    for mass in masses:
        files.append(os.environ["SIG_R3"] + "taustarToTauZ_m" + mass + "_" + year + ".root")

files = [os.environ["SIG_R3"] + ""]
p = PostProcessor(".", files, cut="", branchsel=None, postfix="", modules=[tauPythiaBugSearcherConstr()], histFileName="hists.root", histDirName="Hists")
p.run()
