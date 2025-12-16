#Identify ETau channel events 

#from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel
from PhysicsTools.NanoAODTools.postprocessing.utils.GenTools import prodChainContains, getProdChain
from PhysicsTools.NanoAODTools.postprocessing.utils.Tools import deltaPhi, deltaR, isBetween, getSFFile, yearToEGMSfYr

from correctionlib import _core as corrLib
import gzip
from ROOT import TLorentzVector
from math import cos


# ----------------------------------------------------------------------------------------------------------------------------

class ETauProducer(Module):

    def __init__(self, year):
        self.year = year
        if year in ["2016", "2016post", "2017", "2018"]:
            self.era = 2
        elif year in ["2022", "2022post", "2023", "2023post", "2024"]:
            self.era = 3
        else:
            print("ERROR: Unrecognized year passed to ETauProducer!")  
            exit(1)
        
        if self.year != "2024":
            sfFileName = getSFFile(year=year, pog="EGM")
            with gzip.open(sfFileName,'rt') as fil:
                unzipped = fil.read().strip()
            self.egmSFs = corrLib.CorrectionSet.from_string(unzipped)

            sfFileName = getSFFile(year=year, pog="TAU")
            with gzip.open(sfFileName,'rt') as fil:
                unzipped = fil.read().strip()
            self.tauSFs = corrLib.CorrectionSet.from_string(unzipped)
        

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("ETau_eIdx", "I") #"Index to Electrons of electron"
        self.out.branch("ETau_tauIdx", "I", 3) #"Index to Taus of hadronic tau. Entries correspond to [down, nom, up] tau ES scale"
        self.out.branch("ETau_tauProngs", "I", 3) #"Number of prongs of tau (1 or 3). Entries correspond to [down, nom, up] tau ES scale "
        
        #e+tau 
        self.out.branch("ETau_havePair", "O", 3) #"True if have a good e and tau. Entries correspond to [down, nom, up] tau ES scale" 
        self.out.branch("ETau_visM", "F", 3) #"Visible mass of the e+tau pair. Entries correspond to [down, nom, up] tau ES scale"
        self.out.branch("ETau_sign", "I", 3) #"The product of tau.charge and el.charge"
        self.out.branch("ETau_ETauDR", "F", 3) #"Delta R between electron and tau"
        self.out.branch("ETau_ETauDPhi", "F", 3) #"Delta phi between the electron and tau"
        self.out.branch("ETau_tau1ZDPhi", "F", 3) #"Delta phi between the tau and the Z. Using 'tau1' here for convenient multichannel plotting"
        self.out.branch("ETau_tau2ZDPhi", "F") #"Delta phi between the el and the Z. Using 'tau2' here for convenient multichannel plotting"
        #self.out.branch("ETau_highPtGenMatch", "O") #"True if the higher pt tau decay matched to the GEN taustar tau"
        #self.out.branch("ETau_highPtCollM", "F") #"Either the min or max coll m, whichever was from the higher pt tau decay"
        #self.out.branch("ETau_ETauCos2DPhi", "F") #"cos^2(tau.phi-e.phi)"

        #e+tau+Z
        self.out.branch("ETau_haveTrip", "O", 3) #"True if have a good e, tau, and Z. Entries correspond to [down, nom, up] tau ES scale"
        self.out.branch("ETau_minCollM", "F", 3) #"The smaller collinear mass of either e+nu+Z or tau+nu+Z. Uses best Z of recl vs reco Z mass. Entries correspond to [down, nom, up] tau ES scale"
        self.out.branch("ETau_maxCollM", "F", 3) #"The larger collinear mass of either e+nu+Z or tau+nu+Z. Uses best Z of recl vs reco Z mass. Entries correspond to [down, nom, up] tau ES scale"

        #Scale factors
        self.out.branch("ETau_tauESCorr" , "F", 3) #"The energy scale correction applied to the tau [down, nom, up]"
        self.out.branch("ETau_tauVsESF", "F", 3) #"DeepTau tau vs e SFs [down, nom, up]"
        self.out.branch("ETau_tauVsMuSF", "F", 3) #"DeepTau tau vs mu SFs [down, nom, up]"
        self.out.branch("ETau_tauVsJetSF", "F", 3) #"DeepTau tau vs jet SFs [down, nom, up]"
        self.out.branch("ETau_eIDSF", "F", 3) #"Electron ID SFs [down, nom, up]"

        self.out.branch("ETau_trigMatchTau", "O", 3) #"True if the tau matches the single-tau trigger obj"
        #self.out.branch("ETau_trigMatchETau", "O", 3) #"True if the reco el and tau fired the e-tau cross trigger"
        self.out.branch("ETau_isCand", "O", 3) #"True if the event is good e+tau+Z event. Entries correspond to [down, nom, up] tau ES scale"
        
    def analyze(self, event):
        eIdx = -1
        tauIdx = [-1, -1, -1]
        tauProngs = [-1, -1, -1]
        eTauDR = [-999.99, -999.99, -999.99]
        eTauDPhi = [-999.99, -999.99, -999.99]
        tauZDPhi = [-999.99, -999.99, -999.99]
        eZDPhi = -999.99
        cos_tau_el = [-999.99, -999.99, -999.99]
        visM = [-999.99, -999.99, -999.99]#All 3 entry lists correspond to [down, nom, up] of tau ES
        havePair = [False, False, False]
        haveTrip = [False, False, False]
        maxCollM = [-999.99, -999.99, -999.99]
        minCollM = [-999.99, -999.99, -999.99]
        isCand = [False, False, False]
        sign = [0, 0, 0]

        trigMatchTau = [False, False, False]
        trigMatchETau = False

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
        eIDSF = [1, 1, 1]
        
        taus = Collection(event, "Tau")
        electrons = Collection(event, "Electron")

        currTauVsJet = 0
        currTauPt = 0
        theTau = None
        
        if self.era == 3:
            esCorr = [0.97, 1.00, 1.03]
        else:
            esCorr = [1, 1, 1]
        esScaleVars = ["down", "nom", "up"]
        passTauPtCut = [False, False, False]
        for i in range(3):
            for tauI, tau in enumerate(taus):

                if self.era == 2:
                    if tau.decayMode == 5 or tau.decayMode == 6:
                        continue
                    
                    esCorr[i] = self.tauSFs["tau_energy_scale"].evaluate(tau.pt, abs(tau.eta), tau.decayMode, tau.genPartFlav, "DeepTau2017v2p1", esScaleVars[i])
                    tauCorrPt = tau.pt * esCorr[i]
                    # Use the nominal for most ID checks but keep pT cut info for use later in determining candidacy 
                    passTauPtCut[i] = tauCorrPt > 20.0
                    tauID = passTauPtCut[i] and abs(tau.eta) < 2.3 and abs(tau.dz) < 0.2 
                    
                    #WPs chosen to match run3 choices which were based on existing tau pog SFs
                    tauID = tauID and (tau.idDeepTau2017v2p1VSjet & 8) #8= loose
                    tauID = tauID and (tau.idDeepTau2017v2p1VSmu & 8) #8= tight
                    tauID = tauID and (tau.idDeepTau2017v2p1VSe & 2) #2= VVLoose

                    if tauID and tau.idDeepTau2017v2p1VSjet >= currTauVsJet:
                        if tau.idDeepTau2017v2p1VSjet == currTauVsJet:
                            if tauCorrPt < currTauPt:
                                continue
                        tauIdx[i] = tauI
                        theTau = tau
                        theTau.pt = tauCorrPt
                        theTau.mass = theTau.mass * esCorr[1]
                        currTauPt = tauCorrPt
                        currTauVsJet = tau.idDeepTau2017v2p1VSjet

                elif self.era == 3:
                    #Tau POG recommendations https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendationForRun3
                    
                    #esCorr = self.tauSFs["tau_energy_scale"].evaluate(tau.pt, abs(tau.eta), tau.decayMode, tau.genPartFlav, "Loose", "VVLoose", "nom")

                    tauCorrPt = tau.pt * esCorr[i] # Use the nominal for most ID checks but keep pT cut info for use later in determining candidacy 
                    passTauPtCut[i] = tauCorrPt > 20.0 
                    tauID = passTauPtCut[i] and abs(tau.eta) < 2.5 and abs(tau.dz) < 0.2 

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
                        tauIdx[i] = tauI
                        theTau = tau
                        theTau.pt = tauCorrPt
                        theTau.mass = theTau.mass * esCorr[1]
                        currTauPt = tauCorrPt
                        currTauVsJet = tau.idDeepTau2018v2p5VSjet
        

            if tauIdx[i] >= 0:
                theTau = taus[tauIdx[i]]
                if theTau.decayMode >= 0 and theTau.decayMode <= 2:
                    tauProngs[i] = 1
                elif theTau.decayMode == 10 or theTau.decayMode == 11:
                    tauProngs[i] = 3

                
        #Choose the best electron, if multiple good ID's electrons, chooses the hightes pt candidate
        currElPt = 0
        theEl = None
        for elI, el, in enumerate(electrons):
            if event.Z_dm == 1 and (elI == event.Z_d1Idx or elI == event.Z_d2Idx): #Make sure we don't select one of the Z->ee els
                continue 
            if self.era == 2:
                elID = el.pt > 24 and (abs(el.eta + el.deltaEtaSC) < 1.444 or (abs(el.eta + el.deltaEtaSC) > 1.566 and abs(el.eta + el.deltaEtaSC) < 2.5))
                elID = elID and el.mvaFall17V2Iso_WP90
            elif self.era == 3:
                elID = el.pt > 24 and (abs(el.eta + el.deltaEtaSC) < 1.444 or (abs(el.eta + el.deltaEtaSC) > 1.566 and abs(el.eta + el.deltaEtaSC) < 2.5))
                elID = elID and el.mvaIso_WP90
            if elID and el.pt > currElPt:
                eIdx = elI
                theEl = el
                currElPt = el.pt

        if theEl != None and (tauIdx[0]>=0 or tauIdx[1]>=0 or tauIdx[2]>=0): #If we have an e+tau pair, calculate relevant quantities
        
            #Pythia bug means we have to use placeholder SFs for run3
            if self.era == 2:
                for i, syst in enumerate(["down", "nom", "up"]):
                    tauESCorr[i] = self.tauSFs["tau_energy_scale"].evaluate(theTau.pt, abs(theTau.eta), theTau.decayMode, theTau.genPartFlav, "DeepTau2017v2p1", syst)
                    tauVsESF[i] = self.tauSFs["DeepTau2017v2p1VSe"].evaluate(abs(theTau.eta), theTau.genPartFlav, "VVLoose", syst)
                    tauVsMuSF[i] = self.tauSFs["DeepTau2017v2p1VSmu"].evaluate(abs(theTau.eta), theTau.genPartFlav, "Tight", syst)
                    tauVsJetSF[i] = self.tauSFs["DeepTau2017v2p1VSjet"].evaluate(theTau.pt, theTau.decayMode, theTau.genPartFlav, "Loose", "VVLoose", syst, "pt")
            
            if self.year != "2024":
                for i, syst in enumerate(["sfdown", "sf", "sfup"]): #2024 not available yet
                    if self.year == "2023" or self.year == "2023post": #2023 has phi-dependent SFs
                        eIDSF[i] = self.egmSFs["Electron-ID-SF"].evaluate(yearToEGMSfYr[self.year], syst, "wp90iso", theEl.eta + theEl.deltaEtaSC, theEl.pt, theEl.phi)
                    elif self.year == "2022" or self.year == "2022post":
                        eIDSF[i] = self.egmSFs["Electron-ID-SF"].evaluate(yearToEGMSfYr[self.year], syst, "wp90iso", theEl.eta + theEl.deltaEtaSC, theEl.pt)
                    else:
                        eIDSF[i] = self.egmSFs["UL-Electron-ID-SF"].evaluate(yearToEGMSfYr[self.year], syst, "wp90iso", theEl.eta + theEl.deltaEtaSC, theEl.pt)
    

            havePair = [passTauPtCut[0], passTauPtCut[1], passTauPtCut[2]] #We only considered nom tauES scale above

            for i in range(3):
                if not havePair[i]:
                    continue
                
                theTau = taus[tauIdx[i]]
                theTau.mass = theTau.mass * esCorr[i]
                theTau.pt = theTau.pt * esCorr[i]

                eTauDR[i] = theEl.DeltaR(theTau)
                eTauDPhi[i] = deltaPhi(theTau.phi, theEl.phi)
                sign[i] = theTau.charge * theEl.charge

                ePlusTau = theTau.p4() + theEl.p4()
                visM[i] = ePlusTau.M()
                
                #If the event also has a good Z candidate, we can calculate collinear mass
                if event.Z_dm >= 0 and event.Z_dm <= 2:
                    haveTrip[i] = True

                    #collinear approximation 
                    nuTau = TLorentzVector()
                    nuEl = TLorentzVector()
                    theZ = TLorentzVector()
                    
                    cos_nuTau_MET = cos(deltaPhi(theTau.phi, event.MET_phi))
                    cos_nuEl_MET = cos(deltaPhi(theEl.phi, event.MET_phi))
                    cos_tau_el = cos(deltaPhi(theTau.phi, theEl.phi))
                    cos_tau_el_sqrd = cos_tau_el * cos_tau_el
                    
                    if cos_tau_el_sqrd > 0.999: #Avoid divide by zero issues if tau and el have same phi coord
                        cos_tau_el_sqrd_div0Safe = 0.999
                    else:
                        cos_tau_el_sqrd_div0Safe = cos_tau_el_sqrd

                    nuTau_mag = event.MET_pt * (cos_nuTau_MET - (cos_nuEl_MET * cos_tau_el)) / (1. - cos_tau_el_sqrd_div0Safe)
                    nuEl_mag = ((event.MET_pt * cos_nuTau_MET) - nuTau_mag) / cos_tau_el

                    nuTau.SetPtEtaPhiM(nuTau_mag, theTau.eta, theTau.phi, 0.)
                    nuEl.SetPtEtaPhiM(nuEl_mag, theEl.eta, theEl.phi, 0.)

                    if abs(event.ZReClJ_mass - 91.19) < abs(event.Z_mass - 91.19) and event.Z_dm == 0: #By default, choose closest mass (reclustered vs reco) to nominal
                        theZ.SetPtEtaPhiM(event.ZReClJ_pt, event.ZReClJ_eta, event.ZReClJ_phi, event.ZReClJ_mass)
                    else:
                        theZ.SetPtEtaPhiM(event.Z_pt, event.Z_eta, event.Z_phi, event.Z_mass)

                    tauZDPhi[i] = deltaPhi(theTau.phi, theZ.Phi())
                    eZDPhi = deltaPhi(theEl.phi, theZ.Phi())

                    fullElDecay = theEl.p4() + nuEl
                    fullTauDecay = theTau.p4() + nuTau

                    collM_tauZ = (theTau.p4() + nuTau + theZ).M() 
                    collM_elZ = (theEl.p4() + nuEl + theZ).M()
                    minCollM[i] = min(collM_tauZ, collM_elZ)
                    maxCollM[i] = max(collM_tauZ, collM_elZ)

                    if self.era == 3:
                        trigObjs = Collection(event, "TrigObj")
                        tauLeg = False
                        eLeg = False
                        for trigObj in trigObjs: #Per TAU twiki, filter bits are incorrect in nanoAODv12-v13 so only use DR matching of trig objs
                            if abs(trigObj.id) == 15:
                                if deltaR(trigObj, theTau) < 0.5:
                                    trigMatchTau[i] = True
                                    tauLeg = True
                            elif abs(trigObj.id) == 11:
                                if trigObj.filterBits & (2**3) and trigObj.filterBits & (2**6) and deltaR(trigObj, theEl) < 0.1:
                                    eLeg  = True
                    else:
                        #NB: For run2 where the MET trigger is used, no trigger matching can actually be performed. The variable is used for easier run2+run3 combined cuts
                        trigMatchTau[i] = True

                    isCand[i] = haveTrip[i] #A good triplet
                    if self.era == 3:
                        isCand[i] = isCand[i] and event.Trig_tau  #Appropriate trigger
                    else:
                        isCand[i] = isCand[i] and event.Trig_MET

                    isCand[i] = isCand[i] and abs(theEl.DeltaR(theTau)) > 0.5 #Separation of e and tau
                    isCand[i] = isCand[i] and cos_tau_el**2 < 0.99 #dphi separation of the e and tau
                    isCand[i] = isCand[i] and abs(deltaR(theZ.Eta(), theZ.Phi(), theTau.eta, theTau.phi)) > 0.5 #Separation of the Z and tau
                    isCand[i] = isCand[i] and abs(deltaR(theZ.Eta(), theZ.Phi(), theEl.eta, theEl.phi)) > 0.5 #Separation of the Z and el
                    isCand[i] = isCand[i] and isBetween(theTau.phi, theEl.phi, event.MET_phi) #MET is in small angle between tau & el
                    isCand[i] = isCand[i] and minCollM[i] > visM[i] # Collinear mass should be greater than visible mass
                    isCand[i] = isCand[i] and trigMatchTau[i]
            
        self.out.fillBranch("ETau_eIdx", eIdx) 
        self.out.fillBranch("ETau_tauIdx", tauIdx)
        self.out.fillBranch("ETau_tauProngs", tauProngs)
        self.out.fillBranch("ETau_havePair", havePair)
        self.out.fillBranch("ETau_ETauDR", eTauDR)
        self.out.fillBranch("ETau_ETauDPhi", eTauDPhi)
        self.out.fillBranch("ETau_tau1ZDPhi", tauZDPhi)
        self.out.fillBranch("ETau_tau2ZDPhi", eZDPhi)
        self.out.fillBranch("ETau_visM", visM)
        self.out.fillBranch("ETau_sign", sign)
        self.out.fillBranch("ETau_haveTrip", haveTrip) 
        self.out.fillBranch("ETau_minCollM", minCollM)
        self.out.fillBranch("ETau_maxCollM", maxCollM)
        #self.out.fillBranch("ETau_highPtGenMatch", highPtGenMatch)
        #self.out.fillBranch("ETau_highPtCollM", highPtCollM)
        self.out.fillBranch("ETau_tauESCorr", tauESCorr)
        self.out.fillBranch("ETau_tauVsESF",tauVsESF)
        self.out.fillBranch("ETau_tauVsMuSF", tauVsMuSF)
        self.out.fillBranch("ETau_tauVsJetSF", tauVsJetSF)
        self.out.fillBranch("ETau_eIDSF", eIDSF)
        self.out.fillBranch("ETau_trigMatchTau", trigMatchTau)
        #self.out.fillBranch("ETau_trigMatchETau", trigMatchETau)
        self.out.fillBranch("ETau_isCand", isCand) 

        return True
    
# ----------------------------------------------------------------------------------------------------------------------------
    
eTauProducerConstr = lambda year: ETauProducer(year = year)
