#Producer module to identify and store indexes to the most interesting gen particles (taustars + decay products, etc)
#This producer handles taustar -> WNu channel (i.e. tauWNu final state)


from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel
from PhysicsTools.NanoAODTools.postprocessing.framework.GenTools import getDecayChain, getProdChain, prodChainContains


class GenProducerWNu(Module):

    def __init__(self):
        pass
    
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("GenW_tsIdx", "I") #"Idx to last copy of taustar in GenPart"
        self.out.branch("GenW_tauIdx", "I") #"Idx to the last copy of the tau produced alongside the taustar in GenPart"
        self.out.branch("GenW_wIdx", "I") #"Idx to the last copy of the W from the taustar decay in GenPart"
        self.out.branch("GenW_wDM", "I") #"Decay mode of W. 0 = hadronic, 1=electron, 2=muon, 3=tau. -1 default"
        self.out.branch("GenW_wDau1Idx", "I") #"Idx to GenPart of the higher pT daughter of the W if wDM==0 or the charged lepton if wDM>0. -1 default"
        self.out.branch("GenW_wDau2Idx", "I") #"Idx to GenPart of the lower pT daughter of the W if wDM!=0 or the nu if wDM>0. -1 default"
        self.out.branch("GenW_wGenAK8Idx", "I") #"Idx to GenJetAK8 collection jet matching the W from taustar if wDM == 0"
        self.out.branch("GenW_dr_tauW", "F") #"DeltaR between tau and W"

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        
        #Indices default to -1
        tsIdx = -1
        tauIdx = -1
        wIdx = -1
        wDau1Idx = -1
        wDau2Idx = -1
        wDM = -1
        dr_tauW = -999
        wGenAK8Idx = -1

        genParts = Collection(event, "GenPart")

        for idx, genPart in enumerate(genParts):
            #We only care about the last copy of relevant particles
            if not genPart.statusflag('isLastCopy'):
                continue

            if abs(genPart.pdgId) == 4000015: #Found the taustar
                tsIdx = idx
            elif abs(genPart.pdgId) == 15: #Found a tau
                prodChain = getProdChain(idx, genParts)
                if prodChainContains(prodChain, idx = 0): # If this was the tau produced in the CI with the taustar
                    tauIdx = idx
            elif abs(genPart.pdgId) == 24: #Found a W
                prodChain = getProdChain(idx, genParts)
                if prodChainContains(prodChain, pdgID = 4000015): #If this W is the taustar decay product
                    wIdx = idx
                    decayChain = getDecayChain(wIdx, genParts)

                    if len(decayChain) != 2:
                        print("WARNING: Number of W decay chain particles is != 2. Will write -1 indices for this event")
                    else:
                        if abs(event.GenPart_pdgId[wDau1Idx]) <= 6: #Hadronic
                            wDM = 0
                            if event.GenPart_pt[decayChain[0]] > event.GenPart_pt[decayChain[1]]: #Hightest pt daughter gets idx 1 label
                                wDau1Idx = decayChain[0]
                                wDau2Idx = decayChain[1]
                            else:
                                wDau1Idx = decayChain[1]
                                wDau2Idx = decayChain[0]
                            
                            theW = genParts[wIdx]
                            genJetAK8s = Collection(event, "GenJetAK8")
                            for jetIdx, jet in enumerate(genJetAK8s):
                                if jet.DeltaR(theW) < 0.1 and abs(jet.pt - theW.pt) <= 0.1 * theW.pt:
                                    wGenAK8Idx = jetIdx
                                    break
                        elif abs(event.GenPart_pdgId[wDau1Idx]) == 11: #el
                            wDM = 1
                            wDau1Idx = decayChain[1]
                            wDau2Idx = decayChain[2]
                        elif abs(event.GenPart_pdgId[wDau1Idx]) == 12: 
                            wDM = 1
                            wDau1Idx = decayChain[2]
                            wDau2Idx = decayChain[1]
                        elif abs(event.GenPart_pdgId[wDau1Idx]) == 13: #mu
                            wDM = 2
                            wDau1Idx = decayChain[1]
                            wDau2Idx = decayChain[2]
                        elif abs(event.GenPart_pdgId[wDau1Idx]) == 14: 
                            wDM = 2
                            wDau1Idx = decayChain[2]
                            wDau2Idx = decayChain[1]
                        elif abs(event.GenPart_pdgId[wDau1Idx]) == 15: #tau
                            wDM = 3
                            wDau1Idx = decayChain[1]
                            wDau2Idx = decayChain[2]
                        elif abs(event.GenPart_pdgId[wDau1Idx]) == 16: 
                            wDM = 3
                            wDau1Idx = decayChain[2]
                            wDau2Idx = decayChain[1]
            #End W

        #DeltaR calculations
        tau = genParts[tauIdx]
        w = genParts[wIdx]
        dr_tauW = tau.DeltaR(w)

        self.out.fillBranch("GenW_tsIdx", tsIdx)
        self.out.fillBranch("GenW_tauIdx", tauIdx)
        self.out.fillBranch("GenW_wIdx", wIdx)
        self.out.fillBranch("GenW_wDM", wDM)
        self.out.fillBranch("GenW_wDau1Idx", wDau1Idx)
        self.out.fillBranch("GenW_wDau2Idx", wDau2Idx)
        self.out.fillBranch("GenW_wGenAK8Idx", wGenAK8Idx)
        self.out.fillBranch("GenW_dr_tauW", dr_tauW)


# -----------------------------------------------------------------------------------------------------------------------------            

genProducerWNuConstr = lambda: GenProducerWNu()
