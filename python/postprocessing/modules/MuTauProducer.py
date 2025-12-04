#Identify MuTau channel events 

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel
from PhysicsTools.NanoAODTools.postprocessing.utils.Tools import deltaPhi, deltaR, isBetween, getSFFile
from PhysicsTools.NanoAODTools.postprocessing.utils.GenTools import prodChainContains, getProdChain

from correctionlib import _core as corrLib
import gzip
from ROOT import TLorentzVector
from math import cos

# ----------------------------------------------------------------------------------------------------------------------------

class MuTauProducer(Module):

    def __init__(self, year):
        self.year = year
        if year in ["2016", "2016post", "2017", "2018"]:
            self.era = 2
        elif year in ["2022", "2022post", "2023", "2023post", "2024"]:
            self.era = 3
        else:
            print("ERROR: Unrecognized year passed to MuTauProducer!")  
            exit(1)

        if self.year != "2024":
            sfFileName = getSFFile(year=year, pog="MUO")
            with gzip.open(sfFileName,'rt') as fil:
                unzipped = fil.read().strip()
            self.muSFs = corrLib.CorrectionSet.from_string(unzipped)

            sfFileName = getSFFile(year=year, pog="TAU")
            with gzip.open(sfFileName,'rt') as fil:
                unzipped = fil.read().strip()
            self.tauSFs = corrLib.CorrectionSet.from_string(unzipped)
        

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("MuTau_muIdx", "I") #"Index to Muons of muon"
        self.out.branch("MuTau_tauIdx", "I") #"Index to Taus of hadronic tau"
        self.out.branch("MuTau_tauProngs", "I") #"Number if prongs of tau (1 or 3)"
        self.out.branch("MuTau_havePair", "O") #"True if have a good mu and tau"
        self.out.branch("MuTau_MuTauDR", "F") #"DeltaR betwenn mu and tau"
        self.out.branch("MuTau_MuTauDPhi", "F") #"Delta phi between the muon and tau"
        #self.out.branch("MuTau_MuTauCos2DPhi", "F") #"cos^2(tau.phi-mu.phi)"
        self.out.branch("MuTau_visM", "F") #"Visible mass of the mu+tau pair"
        self.out.branch("MuTau_sign", "I") #"The product of tau.charge and mu.charge"
        self.out.branch("MuTau_haveTrip", "O") #"True if have a good mu, tau, and Z"
        self.out.branch("MuTau_minCollM", "F") #"The smaller collinear mass of e+nu+Z or tau+nu+z. Uses best Z of recl vs reco Z mass"
        self.out.branch("MuTau_maxCollM", "F") #"The larger collinear mass of either mu+nu+Z or tau+nu+Z. Uses best Z of recl vs reco Z mass"
        #self.out.branch("MuTau_highPtGenMatch", "O") #"True if the higher pt tau decay matched to the GEN taustar tau"
        #self.out.branch("MuTau_highPtCollM", "F") #"Either the min or max coll m, whichever was from the higher pt tau decay"
        #Scale factors
        self.out.branch("MuTau_tauESCorr" , "F", 3) #"The energy scale correction applied to the tau [down, nom, up]"
        self.out.branch("MuTau_tauVsESF", "F", 3) #"DeepTau tau vs e SFs [down, nom, up]"
        self.out.branch("MuTau_tauVsMuSF", "F", 3) #"DeepTau tau vs mu SFs [down, nom, up]"
        self.out.branch("MuTau_tauVsJetSF", "F", 3) #"DeepTau tau vs jet SFs [down, nom, up]"
        self.out.branch("MuTau_muIDSF", "F", 3) #"Muon ID SFs [down, nom, up]"
        #Trigger matching
        self.out.branch("MuTau_trigMatchTau", "O") #"True if the tau matches the single-tau trigger obj"
        self.out.branch("MuTau_trigMatchMuTau", "O") #"True if the reco mu and tau fired the mu-tau cross trigger"
        self.out.branch("MuTau_isCand", "O") #"True if the event is good mu+tau+Z event"

    def analyze(self, event):
        muIdx = -1
        tauIdx = -1
        tauProngs = -1
        muTauDR = -999.99
        muTauDPhi = -999.99
        cos_tau_mu = -999.99
        visM = -999.99
        sign = 0
        havePair = False
        haveTrip = False
        maxCollM = -999.99
        minCollM = -999.99
        highPtGenMatch = False
        highPtCollM = -999.99
        isCand = False
        trigMatchTau = False
        trigMatchMuTau = False

        if self.era == 3:
            tauESCorr = [0.97, 1, 1.03]
            tauVsESF = [0.94, 1, 1.06] 
            tauVsJetSF = [0.94, 1, 1.06]
            tauVsMuSF = [0.94, 1, 1.06]
        else:
            tauESCorr = [1, 1, 1]
            tauVsESF = [1, 1, 1] 
            tauVsJetSF = [1, 1, 1]
            tauVsMuSF = [1, 1, 1]
        muIDSF = [1, 1, 1]
        
        taus = Collection(event, "Tau")
        muons = Collection(event, "Muon")

        currTauVsJet = 0
        currTauPt = 0
        theTau = None
        for tauI, tau in enumerate(taus):
            if self.era == 2:
                if tau.decayMode == 5 or tau.decayMode == 6:
                    continue
                esCorr = self.tauSFs["tau_energy_scale"].evaluate(tau.pt, abs(tau.eta), tau.decayMode, tau.genPartFlav, "DeepTau2017v2p1", "nom")
                tauCorrPt = tau.pt * esCorr 
                tauID = tauCorrPt> 20 and abs(tau.eta) < 2.3 and abs(tau.dz) < 0.2 
                #WPs chosen to match run3 choices which were based on existing tau pog SFs
                tauID = tauID and (tau.idDeepTau2017v2p1VSjet & 8) #8= loose
                tauID = tauID and (tau.idDeepTau2017v2p1VSmu & 8) #8= tight
                tauID = tauID and (tau.idDeepTau2017v2p1VSe & 2) #2= VVLoose

                if tauID and tau.idDeepTau2017v2p1VSjet >= currTauVsJet:
                    if tau.idDeepTau2017v2p1VSjet == currTauVsJet:
                        if tauCorrPt < currTauPt:
                            continue
                    tauIdx = tauI
                    theTau = tau
                    theTau.pt = tauCorrPt
                    theTau.mass = theTau.mass * esCorr
                    currTauPt = tauCorrPt
                    currTauVsJet = tau.idDeepTau2017v2p1VSjet
            elif self.era == 3:
                #Tau POG recommendations https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3
                esCorr = 1.00
                #esCorr = self.tauSFs["tau_energy_scale"].evaluate(tau.pt, abs(tau.eta), tau.decayMode, tau.genPartFlav, "Loose", "VVLoose", "nom")
                tauCorrPt = tau.pt * esCorr
                tauID = tauCorrPt > 20 and abs(tau.eta) < 2.5 and abs(tau.dz) < 0.2 
                #WPs chosen based on existing tau pog SFs
                tauID = tauID and tau.idDeepTau2018v2p5VSjet >= 4 #4= loose
                tauID = tauID and tau.idDeepTau2018v2p5VSmu >= 4 #4= tight
                tauID = tauID and tau.idDeepTau2018v2p5VSe >= 2 #2= VVLoose

                #I believe this is already applied but included here anyway for safety
                tauID = tauID and tau.decayMode != 5 and tau.decayMode != 6 

                if tauID and tau.idDeepTau2018v2p5VSjet >= currTauVsJet:
                    if tau.idDeepTau2018v2p5VSjet == currTauVsJet:
                        if tauCorrPt < currTauPt:
                            continue
                    tauIdx = tauI
                    theTau = tau
                    theTau.pt = tauCorrPt
                    theTau.mass = theTau.mass * esCorr
                    currTauPt = tauCorrPt
                    currTauVsJet = tau.idDeepTau2018v2p5VSjet
        
        if theTau != None:
            if theTau.decayMode >= 0 and theTau.decayMode <= 2:
                tauProngs = 1
            elif theTau.decayMode == 10 or theTau.decayMode == 11:
                tauProngs = 3
        
        currMuPt = 0
        theMu = None
        for muI, mu in enumerate(muons):
            if event.Z_dm == 2 and (muI == event.Z_d1Idx or muI == event.Z_d2Idx):#Make sure we don't select one of the Z->mumu muons
                continue
            if self.era == 2:
                muID = mu.pt > 20.0 and abs(mu.eta) < 2.4 and mu.mediumId
            elif self.era == 3:
                muID = mu.pt > 20.0 and abs(mu.eta) < 2.4 and mu.mediumId
            if muID and mu.pt > currMuPt:
                muIdx = muI
                theMu = mu
                currMuPt = mu.pt

        if theTau != None and theMu != None: #If we have an mu+tau pair, calculate relevant quantities
            havePair = True

            muTauDR = theMu.DeltaR(theTau)
            muTauDPhi = deltaPhi(theTau.phi, theMu.phi)
            muPlusTau = theTau.p4() + theMu.p4()
            visM = muPlusTau.M()
            sign = theTau.charge * theMu.charge

            #Pythia bug means we have to use placeholder SFs for run3
            if self.era == 2:
                for i, syst in enumerate(["down", "nom", "up"]):
                    tauESCorr[i] = self.tauSFs["tau_energy_scale"].evaluate(theTau.pt, abs(theTau.eta), theTau.decayMode, theTau.genPartFlav, "DeepTau2017v2p1", syst)
                    tauVsESF[i] = self.tauSFs["DeepTau2017v2p1VSe"].evaluate(abs(theTau.eta), theTau.genPartFlav, "VVLoose", syst)
                    tauVsMuSF[i] = self.tauSFs["DeepTau2017v2p1VSmu"].evaluate(abs(theTau.eta), theTau.genPartFlav, "Tight", syst)
                    tauVsJetSF[i] = self.tauSFs["DeepTau2017v2p1VSjet"].evaluate(theTau.pt, theTau.decayMode, theTau.genPartFlav, "Loose", "VVLoose", syst, "pt")
            for i, syst in enumerate(["systdown", "nominal", "systup"]):
                muIDSF[i] = self.muSFs["NUM_MediumID_DEN_TrackerMuons"].evaluate(abs(theMu.eta), theMu.pt, syst)

            #If the event also has a good Z candidate, we can calculate collinear mass
            if event.Z_dm >= 0 and event.Z_dm <= 2:
                haveTrip = True

                #collinear approximation 
                nuTau = TLorentzVector()
                nuMu = TLorentzVector()
                cos_nuTau_MET = cos(deltaPhi(theTau.phi, event.MET_phi))
                cos_nuMu_MET = cos(deltaPhi(theMu.phi, event.MET_phi))
                cos_tau_mu = cos(deltaPhi(theTau.phi, theMu.phi))
                cos_tau_mu_sqrd = cos_tau_mu * cos_tau_mu
                
                if cos_tau_mu_sqrd > 0.999: #Avoid divide by zero issues if tau and mu have same phi coord
                    cos_tau_mu_sqrd_div0Safe = 0.999
                else:
                    cos_tau_mu_sqrd_div0Safe = cos_tau_mu_sqrd

                nuTau_mag = event.MET_pt * (cos_nuTau_MET - (cos_nuMu_MET * cos_tau_mu)) / (1. - cos_tau_mu_sqrd_div0Safe)
                nuMu_mag = ((event.MET_pt * cos_nuTau_MET) - nuTau_mag) / cos_tau_mu

                nuTau.SetPtEtaPhiM(nuTau_mag, theTau.eta, theTau.phi, 0.)
                nuMu.SetPtEtaPhiM(nuMu_mag, theMu.eta, theMu.phi, 0.)

                fullMuDecay = theMu.p4() + nuMu
                fullTauDecay = theTau.p4() + nuTau

                theZ = TLorentzVector()
                if abs(event.ZReClJ_mass - 91.19) < abs(event.Z_mass - 91.19) and event.Z_dm == 0:
                    theZ.SetPtEtaPhiM(event.ZReClJ_pt, event.ZReClJ_eta, event.ZReClJ_phi, event.ZReClJ_mass)
                else:
                    theZ.SetPtEtaPhiM(event.Z_pt, event.Z_eta, event.Z_phi, event.Z_mass)

                collM_tauZ = (theTau.p4() + nuTau + theZ).M() 
                collM_muZ = (theMu.p4() + nuMu + theZ).M()
                minCollM = min(collM_tauZ, collM_muZ)
                maxCollM = max(collM_tauZ, collM_muZ)

                
                isCand = haveTrip #A good triplet
                if self.era == 2: #Appropriate trigger
                    isCand = isCand and event.Trig_MET
                elif self.era == 3:
                    isCand = isCand and event.Trig_tau  
                isCand = isCand and abs(theMu.DeltaR(theTau)) > 0.5 #Separation of mu and tau
                isCand = isCand and cos_tau_mu**2 < 0.99 #DPhi separation of the mu and tau
                isCand = isCand and abs(deltaR(theZ.Eta(), theZ.Phi(), theTau.eta, theTau.phi)) > 0.5 #Separation of the Z and tau
                isCand = isCand and abs(deltaR(theZ.Eta(), theZ.Phi(), theMu.eta, theMu.phi)) > 0.5 #Separation of the Z and mu
                isCand = isCand and isBetween(theTau.phi, theMu.phi, event.MET_phi) #MET is in small angle between tau & mu
                isCand = isCand and minCollM > visM # Collinear mass should be greater than visible mass
                isCand = isCand and not event.ETau_isCand

        if isCand and self.era == 3:
            trigObjs = Collection(event, "TrigObj")
            tauLeg = False
            muLeg = False
            for trigObj in trigObjs:
                if abs(trigObj.id) == 15: #Per TAU twiki, filter bits are incorrect in nanoAODv12-v13 so only use DR matching of trig objs
                    if deltaR(trigObj, theTau) < 0.5:
                        trigMatchTau = True
                        tauLeg = True
                elif abs(trigObj.id) == 13: 
                    if trigObj.filterBits & (2**2) and trigObj.filterBits & (2**6) and deltaR(trigObj, theMu) < 0.1: 
                        muLeg  = True
            trigMatchMuTau = tauLeg and muLeg
        elif isCand and self.era == 2:
            #NB: For run2 where the MET trigger is used, no trigger matching can actually be performed. The variable is used for easier run2+run3 combined cuts
            trigMatchTau = True 

        isCand = isCand and trigMatchTau

        self.out.fillBranch("MuTau_muIdx", muIdx) 
        self.out.fillBranch("MuTau_tauIdx", tauIdx)
        self.out.fillBranch("MuTau_tauProngs", tauProngs)
        self.out.fillBranch("MuTau_havePair", havePair)
        self.out.fillBranch("MuTau_MuTauDR", muTauDR)
        self.out.fillBranch("MuTau_MuTauDPhi", muTauDPhi)
        self.out.fillBranch("MuTau_visM", visM)
        self.out.fillBranch("MuTau_sign", sign)
        self.out.fillBranch("MuTau_haveTrip", haveTrip) 
        self.out.fillBranch("MuTau_minCollM", minCollM)
        self.out.fillBranch("MuTau_maxCollM", maxCollM)
        #self.out.fillBranch("MuTau_highPtGenMatch", highPtGenMatch)
        #self.out.fillBranch("MuTau_highPtCollM", highPtCollM)
        self.out.fillBranch("MuTau_tauESCorr", tauESCorr)
        self.out.fillBranch("MuTau_tauVsESF",tauVsESF)
        self.out.fillBranch("MuTau_tauVsMuSF", tauVsMuSF)
        self.out.fillBranch("MuTau_tauVsJetSF", tauVsJetSF)
        self.out.fillBranch("MuTau_muIDSF", muIDSF)
        self.out.fillBranch("MuTau_trigMatchTau", trigMatchTau)
        self.out.fillBranch("MuTau_trigMatchMuTau", trigMatchMuTau)
        self.out.fillBranch("MuTau_isCand", isCand) 

        return True
    
# ----------------------------------------------------------------------------------------------------------------------------
    
muTauProducerConstr = lambda year: MuTauProducer(year = year)
