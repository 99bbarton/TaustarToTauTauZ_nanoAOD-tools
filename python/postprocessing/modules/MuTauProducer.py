#Identify ETau channel events 

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

from ROOT import TLorentzVector
from math import cos

# ----------------------------------------------------------------------------------------------------------------------------

class MuTauProducer(Module):

    def __init__(self, era):
        self.era = era
    
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("MuTau_muIdx", "I") #"Index to Muons of muon"
        self.out.branch("MuTau_tauIdx", "I") #"Index to Taus of hadronic tau"
        self.out.branch("MuTau_tauProngs", "I") #"Number if prongs of tau (1 or 3)"
        self.out.branch("MuTau_havePair", "O") #"True if have a good mu and tau"
        self.out.branch("MuTau_muTauDR", "F") #"DeltaR betwenn mu and tau"
        self.out.branch("MuTau_muTauDPhi", "F") #"Delta phi between mu and tau"
        self.out.branch("MuTau_visM", "F") #"Visible mass of the mu+tau pair"
        self.out.branch("MuTau_haveTrip", "O") #"True if have a good mu, tau, and Z"
        self.out.branch("MuTau_minCollM", "F") #"The smaller collinear mass of e+nu+Z or tau+nu+z"
        self.out.branch("MuTau_maxCollM", "F") #"The larger collinear mass of either mu+nu+Z or tau+nu+Z"
        self.out.branch("MuTau_isCand", "O") #"True if the event is good mu+tau+Z event"

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        muIdx = -1
        tauIdx = -1
        tauProngs = -1
        muTauDR = -999.99
        muTauDPhi = -999.99
        visM = -999.99
        havePair = False
        haveTrip = False
        maxCollM = -999.99
        minCollM = -999.99
        isCand = False
        
        taus = Collection(event, "Tau")
        muons = Collection(event, "Muon")

        currTauVsJet = 0
        currTauPt = 0
        theTau = None
        for tauI, tau in enumerate(taus):
            if self.era == 2: #TODO
                print("ERROR: run2 tauID not implemented in mutau producer!")
            elif self.era == 3:
                tauID = tau.pt > 27 and abs(tau.eta) < 2.1 and abs(tau.dz) < 0.2 
                #WPs chosen based on existing tau pog SFs
                tauID = tauID and tau.idDeepTau2018v2p5VSjet >= 4 #4= loose
                tauID = tauID and tau.idDeepTau2018v2p5VSmu >= 4 #4= tight
                tauID = tauID and tau.idDeepTau2018v2p5VSe >= 2 #2= VVLoose

            if tauID and tau.idDeepTau2018v2p5VSjet >= currTauVsJet:
                if tau.idDeepTau2018v2p5VSjet == currTauVsJet:
                    if tau.pt < currTauPt:
                        continue
                tauIdx = tauI
                theTau = tau
                currTauPt = tau.pt
                currTauVsJet = tau.idDeepTau2018v2p5VSjet
        
        if theTau != None:
            if theTau.decayMode <= 2 and theTau.decayMode > 0:
                tauProngs = 1
            elif theTau.decayMode >= 10:
                tauProngs = 3
        
        currMuPt = 0
        theMu = None
        for muI, mu in enumerate(muons):
            if event.Z_dm == 2 and (muI == event.Z_d1Idx or muI == event.Z_d2Idx):#Make sure we don't select one of the Z->mumu muons
                continue
            if self.era == 2:#TODO
                print("ERROR: run2 mu ID not implemented in mutau producer!")
            elif self.era == 3:
                muID = mu.pt > 20.0 and abs(mu.eta) < 2.4 and mu.mediumId
            if muID and mu.pt > currMuPt:
                muIdx = muI
                theMu = mu
                currMuPt = mu.pt

        if theTau != None and theMu != None: #If we have an mu+tau pair, calculate relevant quantities
            havePair = True

            muTauDR = theMu.DeltaR(theTau)
            muTauDPhi = theMu.phi - theTau.phi
            muPlusTau = theTau.p4() + theMu.p4()
            visM = muPlusTau.M()

            #If the event also has a good Z candidate, we can calculate collinear mass
            if event.Z_dm >= 0 and event.Z_dm <= 2:
                haveTrip = True

                #collinear approximation 
                nuTau = TLorentzVector()
                nuMu = TLorentzVector()
                cos_nuTau_MET = cos(theTau.phi - event.MET_phi)
                cos_nuMu_MET = cos(theMu.phi - event.MET_phi)
                cos_tau_mu = cos(theTau.phi - theMu.phi)

                if (1.0 - cos_tau_mu) < 0.001: #Avoid divide by zero issues if tau and el have same phi coord
                    cos_tau_mu_temp = 0.999
                else:
                    cos_tau_mu_temp = cos_tau_mu

#TODO check the direction of the decomposed MET vs the visible objects
                nuTau_mag = event.MET_pt * (cos_nuTau_MET - (cos_nuMu_MET * cos_tau_mu)) / (1. - (cos_tau_mu_temp**2))
                nuEl_mag = ((event.MET_pt * cos_nuTau_MET) - nuTau_mag) / cos_tau_mu

                nuTau.SetPtEtaPhiM(nuTau_mag, theTau.eta, theTau.phi, 0.)
                nuMu.SetPtEtaPhiM(nuEl_mag, theMu.eta, theMu.phi, 0.)

                theZ = TLorentzVector()
                theZ.SetPtEtaPhiM(event.Z_pt, event.Z_eta, event.Z_phi, event.Z_mass)

                collM_tauZ = (theTau.p4() + nuTau + theZ).M() 
                collM_muZ = (theMu.p4() + nuMu + theZ).M()
                minCollM = min(collM_tauZ, collM_muZ)
                maxCollM = max(collM_tauZ, collM_muZ)

                isCand = haveTrip #A good triplet
                isCand = isCand and (event.Trig_tau or event.Trig_muTau)  #Appropriate trigger
                isCand = isCand and (abs(theMu.DeltaR(theTau)) > 0.4) #Separation of mu and tau
                isCand = isCand and (cos_tau_mu**2 < 0.95) #DPhi separation of the mu and tau

        self.out.fillBranch("MuTau_muIdx", muIdx) 
        self.out.fillBranch("MuTau_tauIdx", tauIdx)
        self.out.fillBranch("MuTau_tauProngs", tauProngs)
        self.out.fillBranch("MuTau_havePair", havePair)
        self.out.fillBranch("MuTau_muTauDR", muTauDR)
        self.out.fillBranch("MuTau_muTauDPhi", muTauDPhi)
        self.out.fillBranch("MuTau_visM", visM) 
        self.out.fillBranch("MuTau_haveTrip", haveTrip) 
        self.out.fillBranch("MuTau_minCollM", minCollM)
        self.out.fillBranch("MuTau_maxCollM", maxCollM)
        self.out.fillBranch("MuTau_isCand", isCand) 

        return True
    
# ----------------------------------------------------------------------------------------------------------------------------
    
muTauProducerConstr = lambda era: MuTauProducer(era = era)
