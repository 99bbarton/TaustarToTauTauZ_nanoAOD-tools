#Identify ETau channel events 

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel

from ROOT import TLorentzVector
from math import cos


# ----------------------------------------------------------------------------------------------------------------------------

class ETauProducer(Module):

    def __init__(self, era):
        self.era = era

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("ETau_eIdx", "I") #"Index to Electrons of electron"
        self.out.branch("ETau_tauIdx", "I") #"Index to Taus of hadronic tau"
        self.out.branch("ETau_tauProngs", "I") #"Number of prongs of tau (1 or 3) "
        
        #e+tau 
        self.out.branch("ETau_havePair", "O") #"True if have a good e and tau"
        self.out.branch("ETau_eTauDR", "F") #"Delta R between electron and tau"
        self.out.branch("ETau_eTauDPhi", "F") #"Delta Phi between electron and tau"
        self.out.branch("ETau_visM", "F") #"Visible mass of the e+tau pair"
        
        #e+tau+Z
        self.out.branch("ETau_haveTrip", "O") #"True if have a good e, tau, and Z"
        self.out.branch("ETau_minCollM", "F") #"The smaller collinear mass of either e+nu+Z or tau+nu+Z"
        self.out.branch("ETau_maxCollM", "F") #"The larger collinear mass of either e+nu+Z or tau+nu+Z"
        self.out.branch("ETau_isCand", "O") #"True if the event is good e+tau+Z event"

    def analyze(self, event):
        eIdx = -1
        tauIdx = -1
        tauProngs = -1
        eTauDR = -999.99
        eTauDPhi = -999.99
        visM = -999.99
        havePair = False
        haveTrip = False
        maxCollM = -999.99
        minCollM = -999.99
        isCand = False

        taus = Collection(event, "Tau")
        electrons = Collection(event, "Electron")

        currTauVsJet = 0
        currTauPt = 0
        theTau = None
        for tauI, tau in enumerate(taus):

            if self.era == 2:#TODO
                print("ERROR: run2 tauID not implemented in etau producer!")
            elif self.era == 3:
                tauID = tau.pt > 30 and abs(tau.eta) < 2.1 and abs(tau.dz) < 0.2 
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
        #Choose the best electron, if multiple good ID's electrons, chooses the hightes pt candidate
        currElPt = 0
        theEl = None
        for elI, el, in enumerate(electrons):
            if event.Z_dm == 1 and (elI == event.Z_d1Idx or elI == event.Z_d2Idx): #Make sure we don't select one of the Z->ee els
                continue 
            if self.era == 2:#TODO
                print("ERROR: run2 elID not implemented in etau producer!")
            elif self.era == 3:
                elID = el.pt > 24 and (abs(el.eta) < 1.444 or (abs(el.eta) > 1.566 and abs(el.eta) < 2.5))
                elID = elID and el.mvaIso_WP80
            if elID and el.pt > currElPt:
                eIdx = elI
                theEl = el
                currElPt = el.pt

        if theTau != None and theEl != None: #If we have an e+tau pair, calculate relevant quantities
            havePair = True

            eTauDR = theEl.DeltaR(theTau)
            eTauDPhi = theEl.phi - theTau.phi
            ePlusTau = theTau.p4() + theEl.p4()
            visM = ePlusTau.M()

            #If the event also has a good Z candidate, we can calculate collinear mass
            if event.Z_dm >= 0 and event.Z_dm <= 2:
                haveTrip = True

                #collinear approximation 
                nuTau = TLorentzVector()
                nuEl = TLorentzVector()
                cos_nuTau_MET = cos(theTau.phi - event.MET_phi)
                cos_nuEl_MET = cos(theEl.phi - event.MET_phi)
                cos_tau_el = cos(theTau.phi - theEl.phi)
                
                if abs(1.0 - cos_tau_el) < 0.001: #Avoid divide by zero issues if tau and el have same phi coord
                    cos_tau_el_temp = 0.999
                else:
                    cos_tau_el_temp = cos_tau_el

#TODO check the direction of the decomposed MET vs the visible objects
                nuTau_mag = event.MET_pt * (cos_nuTau_MET - (cos_nuEl_MET * cos_tau_el)) / (1. - (cos_tau_el_temp**2))
                nuEl_mag = ((event.MET_pt * cos_nuTau_MET) - nuTau_mag) / cos_tau_el

                nuTau.SetPtEtaPhiM(nuTau_mag, theTau.eta, theTau.phi, 0.)
                nuEl.SetPtEtaPhiM(nuEl_mag, theEl.eta, theEl.phi, 0.)

                theZ = TLorentzVector()
                theZ.SetPtEtaPhiM(event.Z_pt, event.Z_eta, event.Z_phi, event.Z_mass)

                collM_tauZ = (theTau.p4() + nuTau + theZ).M() 
                collM_elZ = (theEl.p4() + nuEl + theZ).M()
                minCollM = min(collM_tauZ, collM_elZ)
                maxCollM = max(collM_tauZ, collM_elZ)

                isCand = haveTrip #A good triplet
                isCand = isCand and (event.Trig_tau or event.Trig_eTau)  #Appropriate trigger
                isCand = isCand and (abs(theEl.DeltaR(theTau)) > 0.4) #Separation of e and tau
                isCand = isCand and (cos_tau_el**2 < 0.95) #DPhi separation of the e and tau

        self.out.fillBranch("ETau_eIdx", eIdx) 
        self.out.fillBranch("ETau_tauIdx", tauIdx)
        self.out.fillBranch("ETau_tauProngs", tauProngs)
        self.out.fillBranch("ETau_havePair", havePair)
        self.out.fillBranch("ETau_eTauDR", eTauDR)
        self.out.fillBranch("ETau_eTauDPhi", eTauDPhi)
        self.out.fillBranch("ETau_visM", visM) 
        self.out.fillBranch("ETau_haveTrip", haveTrip) 
        self.out.fillBranch("ETau_minCollM", minCollM)
        self.out.fillBranch("ETau_maxCollM", maxCollM)
        self.out.fillBranch("ETau_isCand", isCand) 

        return True
    
# ----------------------------------------------------------------------------------------------------------------------------
    
eTauProducerConstr = lambda era: ETauProducer(era = era)