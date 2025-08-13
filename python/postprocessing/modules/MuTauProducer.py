#Identify ETau channel events 

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel
from PhysicsTools.NanoAODTools.postprocessing.utils.Tools import deltaPhi, deltaR, isBetween
from PhysicsTools.NanoAODTools.postprocessing.utils.GenTools import prodChainContains, getProdChain

from ROOT import TLorentzVector
from math import cos

# ----------------------------------------------------------------------------------------------------------------------------

class MuTauProducer(Module):

    def __init__(self, era):
        self.era = era
    
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
        self.out.branch("MuTau_haveTrip", "O") #"True if have a good mu, tau, and Z"
        self.out.branch("MuTau_minCollM", "F") #"The smaller collinear mass of e+nu+Z or tau+nu+z"
        self.out.branch("MuTau_maxCollM", "F") #"The larger collinear mass of either mu+nu+Z or tau+nu+Z"
        self.out.branch("MuTau_highPtGenMatch", "O") #"True if the higher pt tau decay matched to the GEN taustar tau"
        self.out.branch("MuTau_highPtCollM", "F") #"Either the min or max coll m, whichever was from the higher pt tau decay"
        self.out.branch("MuTau_isCand", "O") #"True if the event is good mu+tau+Z event"

    def analyze(self, event):
        muIdx = -1
        tauIdx = -1
        tauProngs = -1
        muTauDR = -999.99
        muTauDPhi = -999.99
        cos_tau_mu = -999.99
        visM = -999.99
        havePair = False
        haveTrip = False
        maxCollM = -999.99
        minCollM = -999.99
        highPtGenMatch = False
        highPtCollM = -999.99
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

                #I believe this is already applied but included here anyway for safety
                tauID = tauID and tau.decayMode != 5 and tau.decayMode != 6 

            if tauID and tau.idDeepTau2018v2p5VSjet >= currTauVsJet:
                if tau.idDeepTau2018v2p5VSjet == currTauVsJet:
                    if tau.pt < currTauPt:
                        continue
                tauIdx = tauI
                theTau = tau
                currTauPt = tau.pt
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
            muTauDPhi = deltaPhi(theTau.phi, theMu.phi)
            muPlusTau = theTau.p4() + theMu.p4()
            visM = muPlusTau.M()

            #If the event also has a good Z candidate, we can calculate collinear mass
            if event.Z_dm >= 0 and event.Z_dm <= 2:
                haveTrip = True

                #collinear approximation 
                nuTau = TLorentzVector()
                nuMu = TLorentzVector()
                cos_nuTau_MET = cos(deltaPhi(theTau.phi, event.MET_phi))
                cos_nuMu_MET = cos(deltaPhi(theMu.phi, event.MET_phi))
                cos_tau_mu = cos(deltaPhi(theTau.phi, theMu.phi))

                if (1.0 - cos_tau_mu) < 0.001: #Avoid divide by zero issues if tau and el have same phi coord
                    cos_Tau_mu_div0Safe = 0.999
                else:
                    cos_Tau_mu_div0Safe = cos_tau_mu

                nuTau_mag = event.MET_pt * (cos_nuTau_MET - (cos_nuMu_MET * cos_Tau_mu_div0Safe)) / (1. - (cos_Tau_mu_div0Safe**2))
                nuMu_mag = ((event.MET_pt * cos_nuTau_MET) - nuTau_mag) / cos_Tau_mu_div0Safe

                nuTau.SetPtEtaPhiM(nuTau_mag, theTau.eta, theTau.phi, 0.)
                nuMu.SetPtEtaPhiM(nuMu_mag, theMu.eta, theMu.phi, 0.)

                fullMuDecay = theMu.p4() + nuMu
                fullTauDecay = theTau.p4() + nuTau

                theZ = TLorentzVector()
                if event.ZReClJ_mass > 0:
                    theZ.SetPtEtaPhiM(event.ZReClJ_pt, event.ZReClJ_eta, event.ZReClJ_phi, event.ZReClJ_mass)
                else:
                    theZ.SetPtEtaPhiM(event.Z_pt, event.Z_eta, event.Z_phi, event.Z_mass)

                collM_tauZ = (theTau.p4() + nuTau + theZ).M() 
                collM_muZ = (theMu.p4() + nuMu + theZ).M()
                minCollM = min(collM_tauZ, collM_muZ)
                maxCollM = max(collM_tauZ, collM_muZ)

                genParts = Collection(event, "GenPart")
                if fullMuDecay.Pt() > fullTauDecay.Pt():
                    highPtCollM = collM_muZ
                    prodChain = getProdChain(theMu.genPartIdx, genParts)
                    highPtGenMatch = prodChainContains(prodChain, idx=event.Gen_tsTauIdx)
                elif theTau.genPartIdx >= 0:
                    highPtCollM = collM_tauZ
                    prodChain = getProdChain(event.GenVisTau_genPartIdxMother[theTau.genPartIdx], genParts)
                    highPtGenMatch = prodChainContains(prodChain, idx=event.Gen_tsTauIdx)
                else:
                    #print("In MuTau tau.genPartIdx is negative!")
                    pass

                isCand = haveTrip #A good triplet
                isCand = isCand and event.Trig_tau  #Appropriate trigger
                isCand = isCand and abs(theMu.DeltaR(theTau)) > 0.5 #Separation of mu and tau
                isCand = isCand and cos_tau_mu**2 < 0.95 #DPhi separation of the mu and tau
                isCand = isCand and abs(deltaR(theZ.Eta(), theZ.Phi(), theTau.eta, theTau.phi)) > 0.5 #Separation of the Z and tau
                isCand = isCand and abs(deltaR(theZ.Eta(), theZ.Phi(), theMu.eta, theMu.phi)) > 0.5 #Separation of the Z and mu
                isCand = isCand and isBetween(theTau.phi, theMu.phi, event.MET_phi) #MET is in small angle between tau & mu
                isCand = isCand and minCollM > visM # Collinear mass should be greater than visible mass
                isCand = isCand and not event.ETau_isCand

        trigObjs = Collection(event, "TrigObj")
        tauLeg = False
        muLeg = False
        for trigObj in trigObjs:
            if abs(trigObj.id) == 15:
                if trigObj.filterBits & (2**3) and trigObj.filterBits & (2**10) and trigObj.deltaR(theTau) < 0.1:
                    trigMatchTau = True
                elif trigObj.filterBits & (2**3) and trigObj.filterBits & (2**9) and trigObj.deltaR(theTau) < 0.1:
                    tauLeg = True
            elif abs(trigObj.id) == 13: #TODO verify, these are educated guesses for trig bits
                if trigObj.filterBits & (2**2) and trigObj.filterBits & (2**6) and trigObj.deltaR(theMu) < 0.1: 
                    muLeg  = True
        trigMatchMuTau = tauLeg and muLeg
                

        self.out.fillBranch("MuTau_muIdx", muIdx) 
        self.out.fillBranch("MuTau_tauIdx", tauIdx)
        self.out.fillBranch("MuTau_tauProngs", tauProngs)
        self.out.fillBranch("MuTau_havePair", havePair)
        self.out.fillBranch("MuTau_MuTauDR", muTauDR)
        self.out.fillBranch("MuTau_MuTauDPhi", muTauDPhi)
        self.out.fillBranch("MuTau_visM", visM) 
        self.out.fillBranch("MuTau_haveTrip", haveTrip) 
        self.out.fillBranch("MuTau_minCollM", minCollM)
        self.out.fillBranch("MuTau_maxCollM", maxCollM)
        self.out.fillBranch("MuTau_highPtGenMatch", highPtGenMatch)
        self.out.fillBranch("MuTau_highPtCollM", highPtCollM)
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
    
muTauProducerConstr = lambda era: MuTauProducer(era = era)
