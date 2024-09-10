#Producer which identifies with triggers of interest are satisfied per event
#After this is run, in principle all HLT_* and L1_* branches can be dropped from the trees

#Trigger recommendations taken from:
# Tau: https://twiki.cern.ch/twiki/bin/view/CMS/TauTrigger
# Met: https://twiki.cern.ch/twiki/bin/view/CMS/JetMETPathsRun2#Single_PF_MET_paths_Unprescaled & RUN3 UNAVAILABLE AS OF 10Sep2024
# El : https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary & https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIIISummary 
# Mu : https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT#Details_for_each_year


from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

# -----------------------------------------------------------------------------------------------------------------------------

class TrigProducer(Module):

    def __init__(self, year):
        self.year = year
        pass
    
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Trig_eTau", "O") #"Whether e+tau trigger is satisfied"
        self.out.branch("Trig_muTau", "O") #"Whether mu+tau trigger is satisfied"
        self.out.branch("Trig_tauTau", "O") #"Whether di-tau trigger is satisfied"
        self.out.branch("Trig_tau", "O") #"Whether single tau trigger is satisfied"
        self.out.branch("Trig_tauOR", "O") #"OR of all of the above"
        self.out.branch("Trig_MET", "O") #"Whether the MET trigger is satisfied"
        self.out.branch("Trig_eIso", "O") #"Whether the single iso el trigger is satisfied"
        self.out.branch("Trig_muIso", "O") #"Whether the single iso mu trigger is satisfied"

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        Trig_eTau = False
        Trig_muTau = False
        Trig_tauTau = False
        Trig_tau = False
        Trig_tauOR = False
        Trig_MET = False
        Trig_eIso = False
        Trig_muIso = False

        if self.year == "2016" or self.year == "2016post":
            Trig_eTau = False
            Trig_muTau = False
            Trig_tauTau = False
            Trig_tau = False #No single-tau trig for Run2
            Trig_tauOR = Trig_eTau or Trig_muTau or Trig_tauTau
            Trig_MET = False
            Trig_eIso = event.HLT_Ele27_WPTight_Gsf or event.HLT_Photon175
            Trig_muIso = event.HLT_IsoMu24 or event.HLT_IsoTkMu24
        elif self.year == "2017":
            Trig_eTau = False
            Trig_muTau = False
            Trig_tauTau = False
            Trig_tau = False #No single-tau trig for Run2
            Trig_tauOR = Trig_eTau or Trig_muTau or Trig_tauTau
            Trig_MET = False
            Trig_eIso = event.HLT_Ele35_WPTight_Gsf or event.HLT_Photon200
            Trig_muIso = event.HLT_IsoMu27
        elif self.year == "2018":
            Trig_eTau = event.HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1
            Trig_muTau = event.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1
            Trig_tauTau = event.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg
            Trig_tau = False #No single-tau trig for Run2
            Trig_tauOR = Trig_eTau or Trig_muTau or Trig_tauTau
            Trig_MET = event.HLT_PFMET200_HBHE_BeamHaloCleaned
            Trig_eIso = event.HLT_Ele32_WPTight_Gsf
            Trig_muIso = event.HLT_IsoMu24
        elif self.year == "2022" or self.year == "2022post":
            Trig_eTau = event.HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1
            Trig_muTau = event.HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1
            Trig_tauTau = event.HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1
            Trig_tau = event.HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1
            Trig_tauOR = Trig_eTau or Trig_muTau or Trig_tauTau or Trig_tau
            Trig_MET = event.HLT_PFMET200_BeamHaloCleaned #NOTE: No official trig recommendation from JetMET as of 4Sep2024
            Trig_eIso = event.HLT_Ele30_WPTight_Gsf
            Trig_muIso = event.HLT_IsoMu24
        elif self.year == "2023" or self.year == "2023post":
            Trig_eTau = event.HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1
            Trig_muTau = event.HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1
            Trig_tauTau = event.HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1
            Trig_tau = event.HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1
            Trig_tauOR = Trig_eTau or Trig_muTau or Trig_tauTau or Trig_tau
            Trig_MET = event.HLT_PFMET200_BeamHaloCleaned
            Trig_eIso = event.HLT_Ele30_WPTight_Gsf
            Trig_muIso = event.HLT_IsoMu24

        self.out.fillBranch("Trig_eTau", Trig_eTau)
        self.out.fillBranch("Trig_muTau", Trig_muTau)
        self.out.fillBranch("Trig_tauTau", Trig_tauTau)
        self.out.fillBranch("Trig_tau", Trig_tau)
        self.out.fillBranch("Trig_tauOR", Trig_tauOR)
        self.out.fillBranch("Trig_MET", Trig_MET)
        self.out.fillBranch("Trig_eIso", Trig_eIso)
        self.out.fillBranch("Trig_muIso", Trig_muIso)

        return True
    
    # -----------------------------------------------------------------------------------------------------------------------------

trigProducerConstr = lambda year: TrigProducer(year = year)