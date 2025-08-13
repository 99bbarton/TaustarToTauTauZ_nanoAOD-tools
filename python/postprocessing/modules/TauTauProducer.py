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
        self.out.branch("TauTau_TauTauDR", "F") #"DeltaR between the two taus"
        self.out.branch("TauTau_TauTauDPhi", "F") #"Delta phi between the two taus"
        #self.out.branch("TauTau_TauTauCos2DPhi", "F") #"cos^2(tau1.phi-tau2.phi)"
        self.out.branch("TauTau_visM", "F") #"Visible mass of the tau pair"
        self.out.branch("TauTau_haveTrip", "O") #"True if have two good taus and a Z"
        self.out.branch("TauTau_minCollM", "F") #"The smaller collinear mass of tau1+nu+Z or tau2+nu+Z"
        self.out.branch("TauTau_maxCollM", "F") #"The larger collinear mass of either tau1+nu+Z or tau2+nu+Z"
        self.out.branch("TauTau_highPtGenMatch", "O") #"True if the higher pt tau decay matched to the GEN taustar tau"
        self.out.branch("TauTau_highPtCollM", "F") #"Either the min or max coll m, whichever was from the higher pt tau decay"
        self.out.branch("TauTau_isCand", "O") #"True if the event is good tau+tau+Z event"

        #self.out.branch("TauTau_trigMatchTau", "O") #"True if the event passes the single tau trigger and one tau matches to the trigObj"
        #self.out.branch("TauTau_trigMatchTauTau", "O") #"True if the event passes the d-tau trigger and both taus matche to the trigObj"
        

    def analyze(self, event):
        tau1Idx = -1
        tau2Idx = -1
        tau1Prongs = -1
        tau2Prongs = -1
        tausDR = -999.99
        tausDPhi = -999.99
        cos_tau1_tau2 = -999.99
        visM = -999.99
        havePair = False
        haveTrip = False
        maxCollM = -999.99
        minCollM = -999.99
        highPtGenMatch = False
        highPtCollM = -999.99

        trigMatchTau = False
        trigMatchTauTau = False

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

                #I believe this is already applied but included here anyway for safety
                tauID = tauID and tau.decayMode != 5 and tau.decayMode != 6 
                
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
                
                if abs(tau1.DeltaR(tau2)) < 0.5:
                    tau2Idx = goodTaus[i][0]
                    break #Since list is sorted, as soon as we find a DR separated tau, we're done
            if tau2Idx < 0: #No DR separated tau so choose the second best vsJet and pt ranked tau 
                tau2Idx = goodTaus[1][0]
                tau2 = taus[tau2Idx]
            
            if tau1.decayMode >= 0 and tau1.decayMode <= 2:
                tau1Prongs = 1
            elif tau1.decayMode == 10 or tau1.decayMode == 11:
                tau1Prongs = 3
            if tau2.decayMode >= 0 and tau2.decayMode <= 2:
                tau2Prongs = 1
            elif tau2.decayMode == 10 or tau2.decayMode == 11:
                tau2Prongs = 3


            tausDR = tau1.DeltaR(tau2)
            tausDPhi = deltaPhi(tau1.phi, tau2.phi)
            tauPlusTau = tau1.p4() + tau2.p4()
            visM = tauPlusTau.M()

            #If the event also has a good Z candidate, we can calculate collinear mass
            if event.Z_dm >= 0 and event.Z_dm <= 2:
                haveTrip = True

                #collinear approximation 
                nuTau1 = TLorentzVector()
                nuTau2 = TLorentzVector()
                cos_nuTau1_MET = cos(deltaPhi(tau1.phi, event.MET_phi))
                cos_nuTau2_MET = cos(deltaPhi(tau2.phi, event.MET_phi))
                cos_tau1_tau2 = cos(deltaPhi(tau1.phi, tau2.phi))

                if (1.0 - cos_tau1_tau2) < 0.001: #Avoid divide by zero issues if taus have same phi coord
                    cos_tau1_tau2_div0Safe = 0.999
                else:
                    cos_tau1_tau2_div0Safe = cos_tau1_tau2
                    
                nuTau1_mag = event.MET_pt * (cos_nuTau1_MET - (cos_nuTau2_MET * cos_tau1_tau2_div0Safe)) / (1. - (cos_tau1_tau2_div0Safe**2))
                nuTau2_mag = ((event.MET_pt * cos_nuTau1_MET) - nuTau1_mag) / cos_tau1_tau2_div0Safe

                nuTau1.SetPtEtaPhiM(nuTau1_mag, tau1.eta, tau1.phi, 0.)
                nuTau2.SetPtEtaPhiM(nuTau2_mag, tau2.eta, tau2.phi, 0.)

                fullTau1Decay = tau1.p4() + nuTau1
                fullTau2Decay = tau2.p4() + nuTau2

                theZ = TLorentzVector()
                if event.ZReClJ_mass > 0:
                    theZ.SetPtEtaPhiM(event.ZReClJ_pt, event.ZReClJ_eta, event.ZReClJ_phi, event.ZReClJ_mass)
                else:
                    theZ.SetPtEtaPhiM(event.Z_pt, event.Z_eta, event.Z_phi, event.Z_mass)

                collM_tau1Z = (tau1.p4() + nuTau1 + theZ).M() 
                collM_tau2Z= (tau2.p4() + nuTau2 + theZ).M()
                minCollM = min(collM_tau1Z, collM_tau2Z)
                maxCollM = max(collM_tau1Z, collM_tau2Z)

                genParts = Collection(event, "GenPart")
                if fullTau1Decay.Pt() > fullTau2Decay.Pt() and tau1.genPartIdx >= 0:
                    highPtCollM = collM_tau1Z
                    prodChain = getProdChain(event.GenVisTau_genPartIdxMother[tau1.genPartIdx], genParts)
                    highPtGenMatch = prodChainContains(prodChain, idx=event.Gen_tsTauIdx)
                elif tau2.genPartIdx >= 0:
                    highPtCollM = collM_tau2Z
                    prodChain = getProdChain(event.GenVisTau_genPartIdxMother[tau2.genPartIdx], genParts)
                    highPtGenMatch = prodChainContains(prodChain, idx=event.Gen_tsTauIdx)
                else:
                    pass
                    #print("In TauTau tau1 or tau2 .genPartIdx is negative!")

                isCand = haveTrip #A good triplet
                isCand = isCand and (event.Trig_tau or event.Trig_tauTau)  #Appropriate trigger
                isCand = isCand and abs(tau2.DeltaR(tau1)) > 0.5 #Separation of two taus
                isCand = isCand and cos_tau1_tau2**2 < 0.95 #DPhi separation of the two taus
                isCand = isCand and abs(deltaR(theZ.Eta(), theZ.Phi(), tau1.eta, tau1.phi)) > 0.5 #Separation of the Z and tau1
                isCand = isCand and abs(deltaR(theZ.Eta(), theZ.Phi(), tau2.eta, tau2.phi)) > 0.5 #Separation of the Z and tau2
                isCand = isCand and isBetween(tau1.phi, tau2.phi, event.MET_phi) #MET in small angle between taus
                isCand = isCand and minCollM > visM # Collinear mass should be greater than visible mass
                isCand = isCand and not (event.ETau_isCand or event.MuTau_isCand)

            trigObjs = Collection(event, "TrigObj")
            tau1Match = False
            tau2Match = False
            for trigObj in trigObjs:
                if abs(trigObj.id) == 15:
                    if trigObj.filterBits & (2**3) and trigObj.filterBits & (2**10) and (trigObj.deltaR(tau1) < 0.1 or trigObj.deltaR(tau2) < 0.1):
                        trigMatchTau = True
                    elif trigObj.filterBits & (2**3) and trigObj.filterBits & (2**7) and trigObj.deltaR(tau1) < 0.1:
                        tau1Match = True
                    elif trigObj.filterBits & (2**3) and trigObj.filterBits & (2**7) and trigObj.deltaR(tau2) < 0.1:
                        tau2Match = True
            trigMatchTauTau = tau1Match and tau2Match


        self.out.fillBranch("TauTau_tau1Idx", tau1Idx)
        self.out.fillBranch("TauTau_tau2Idx", tau2Idx)
        self.out.fillBranch("TauTau_tau1Prongs", tau1Prongs)
        self.out.fillBranch("TauTau_tau2Prongs", tau2Prongs)
        self.out.fillBranch("TauTau_havePair", havePair)
        self.out.fillBranch("TauTau_TauTauDR", tausDR)
        self.out.fillBranch("TauTau_TauTauDPhi", tausDPhi)
        self.out.fillBranch("TauTau_visM", visM)
        self.out.fillBranch("TauTau_haveTrip", haveTrip)
        self.out.fillBranch("TauTau_minCollM", minCollM)
        self.out.fillBranch("TauTau_maxCollM", maxCollM)
        self.out.fillBranch("TauTau_highPtGenMatch", highPtGenMatch)
        self.out.fillBranch("TauTau_highPtCollM", highPtCollM)
        self.out.fillBranch("TauTau_tauESCorr", tauESCorr)
        self.out.fillBranch("TauTau_tauVsESF",tauVsESF)
        self.out.fillBranch("TauTau_tauVsMuSF", tauVsMuSF)
        self.out.fillBranch("TauTau_tauVsJetSF", tauVsJetSF)
        self.out.fillBranch("TauTau_trigMatchTau", trigMatchTau)
        self.out.fillBranch("TauTau_trigMatchTauTau", trigMatchTauTau)
        self.out.fillBranch("TauTau_isCand", isCand) 

        return True
    
# ----------------------------------------------------------------------------------------------------------------------------
    
tauTauProducerConstr = lambda era: TauTauProducer(era = era)
