#Identify ETau channel events 

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

from ROOT import TLorentzVector
from math import cos

# ----------------------------------------------------------------------------------------------------------------------------

class TauTauProducer(Module):

    def __init__(self, era):
        self.era = era

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("TauTau_tau1Idx", "I") #"Index to Taus of the best hadronic tau (sorted by vsJets, pt)"
        self.out.branch("TauTau_tau2Idx", "I") #"Index to Taus of second best hadronic tau (sorted by vsJets, pt, and DR>0.4 from best tau)"
        self.out.branch("TauTau_tau1Prongs", "I") #"Number if prongs of tau1 (1 or 3)"
        self.out.branch("TauTau_tau2Prongs", "I") #"Number if prongs of tau2 (1 or 3)"
        self.out.branch("TauTau_havePair", "O") #"True if have two good taus"
        self.out.branch("TauTau_tausDR", "F") #"DeltaR betwenn the two taus"
        self.out.branch("TauTau_tausDPhi", "F") #"Delta phi between the two taus"
        self.out.branch("TauTau_visM", "F") #"Visible mass of the tau pair"
        self.out.branch("TauTau_haveTrip", "O") #"True if have two good taus and a Z"
        self.out.branch("TauTau_minCollM", "F") #"The smaller collinear mass of tau1+nu+Z or tau2+nu+Z"
        self.out.branch("TauTau_maxCollM", "F") #"The larger collinear mass of either tau1+nu+Z or tau2+nu+Z"
        self.out.branch("TauTau_isCand", "O") #"True if the event is good tau+tau+Z event"

    def analyze(self, event):
        tau1Idx = -1
        tau2Idx = -1
        tau1Prongs = -1
        tau2Prongs = -1
        tausDR = -999.99
        tausDPhi = -999.99
        visM = -999.99
        havePair = False
        haveTrip = False
        maxCollM = -999.99
        minCollM = -999.99
        isCand = False
        
        taus = Collection(event, "Tau")
        if len(taus) < 2:
            return True
        
        goodTaus = [] #List of (tauIdx, tau.IDvsJets, tau.pt)
        for tauI, tau in enumerate(taus):
            
            if self.era == 2:#TODO
                print("ERROR: run2 tauID not implemented in tautau producer!")
            elif self.era == 3:
                #Thesholds chosen based on trigger acceptances
                tauID = tau.pt > 35 and abs(tau.eta) < 2.1 and abs(tau.dz) < 0.2 
                #WPs chosen based on existing tau pog SFs
                tauID = tauID and tau.idDeepTau2018v2p5VSjet >= 4 #4= loose
                tauID = tauID and tau.idDeepTau2018v2p5VSmu >= 4 #4= tight
                tauID = tauID and tau.idDeepTau2018v2p5VSe >= 2 #2= VVLoose
            if tauID:
                goodTaus.append((tauI, tau.idDeepTau2018v2p5VSjet, tau.pt))
                            
        if len(goodTaus) >= 2:
            havePair = True
            #Sorts by vsJet score, then pt. Then we can find DR separated pairs
            goodTaus.sort(key=lambda goodTau : goodTau[2], reverse=True)
            goodTaus.sort(key=lambda goodTau : goodTau[1], reverse=True)

            tau1Idx = goodTaus[0][0]
            tau1 = taus[tau1Idx]
            
            for i in range(1, len(goodTaus)):
                tau2 = taus[goodTaus[i][0]]
                
                if abs(tau1.DeltaR(tau2)) < 0.4:
                    tau2Idx = goodTaus[i][0]
                    break #Since list is sorted, as soon as we find a DR separated tau, we're done
            if tau2Idx < 0: #No DR separated tau so choose the second best vsJet and pt ranked tau 
                tau2Idx = goodTaus[1][0]
                tau2 = taus[tau2Idx]
            
            if tau1.decayMode <= 2 and tau1.decayMode > 0:
                tau1Prongs = 1
            elif tau1.decayMode >= 10:
                tau1Prongs = 3
            if tau2.decayMode <= 2 and tau2.decayMode > 0:
                tau2Prongs = 1
            elif tau2.decayMode >= 10:
                tau2Prongs = 3


            tausDR = tau1.DeltaR(tau2)
            tausDPhi = tau1.phi - tau2.phi
            tauPlusTau = tau1.p4() + tau2.p4()
            visM = tauPlusTau.M()

            #If the event also has a good Z candidate, we can calculate collinear mass
            if event.Z_dm >= 0 and event.Z_dm <= 2:
                haveTrip = True

                #collinear approximation 
                nuTau1 = TLorentzVector()
                nuTau2 = TLorentzVector()
                cos_nuTau1_MET = cos(tau1.phi - event.MET_phi)
                cos_nuTau2_MET = cos(tau2.phi - event.MET_phi)
                cos_tau1_tau2 = cos(tau1.phi - tau2.phi)

                if (1.0 - cos_tau1_tau2) < 0.001: #Avoid divide by zero issues if taus have same phi coord
                    cos_tau1_tau2_temp = 0.999
                else:
                    cos_tau1_tau2_temp = cos_tau1_tau2

#TODO check the direction of the decomposed MET vs the visible objects
                nuTau1_mag = event.MET_pt * (cos_nuTau1_MET - (cos_nuTau2_MET * cos_tau1_tau2)) / (1. - (cos_tau1_tau2_temp**2))
                nuTau2_mag = ((event.MET_pt * cos_nuTau1_MET) - nuTau1_mag) / cos_tau1_tau2

                nuTau1.SetPtEtaPhiM(nuTau1_mag, tau1.eta, tau1.phi, 0.)
                nuTau2.SetPtEtaPhiM(nuTau2_mag, tau2.eta, tau2.phi, 0.)

                theZ = TLorentzVector()
                theZ.SetPtEtaPhiM(event.Z_pt, event.Z_eta, event.Z_phi, event.Z_mass)

                collM_tau1Z = (tau1.p4() + nuTau1 + theZ).M() 
                collM_tau2Z= (tau2.p4() + nuTau2 + theZ).M()
                minCollM = min(collM_tau1Z, collM_tau2Z)
                maxCollM = max(collM_tau1Z, collM_tau2Z)

                isCand = haveTrip #A good triplet
                isCand = isCand and (event.Trig_tau or event.Trig_tauTau)  #Appropriate trigger
                isCand = isCand and (abs(tau2.DeltaR(tau1)) > 0.4) #Separation of e and tau
                isCand = isCand and (cos_tau1_tau2**2 < 0.95) #DPhi separation of the e and tau


        self.out.fillBranch("TauTau_tau1Idx", tau1Idx)
        self.out.fillBranch("TauTau_tau2Idx", tau2Idx)
        self.out.fillBranch("TauTau_tau1Prongs", tau1Prongs)
        self.out.fillBranch("TauTau_tau2Prongs", tau2Prongs)
        self.out.fillBranch("TauTau_havePair", havePair)
        self.out.fillBranch("TauTau_tausDR", tausDR)
        self.out.fillBranch("TauTau_tausDPhi", tausDPhi)
        self.out.fillBranch("TauTau_visM", visM)
        self.out.fillBranch("TauTau_haveTrip", haveTrip)
        self.out.fillBranch("TauTau_minCollM", minCollM)
        self.out.fillBranch("TauTau_maxCollM", maxCollM)
        self.out.fillBranch("TauTau_isCand", isCand)

        return True
    
# ----------------------------------------------------------------------------------------------------------------------------
    
tauTauProducerConstr = lambda era: TauTauProducer(era = era)