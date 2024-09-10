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
        self.out.branch("Gen_tsIdx", "I") #"Idx to last copy of taustar in GenPart"
        self.out.branch("Gen_tauIdx", "I") #"Idx to the last copy of the tau produced alongside the taustar in GenPart"
        self.out.branch("Gen_tauDM", "I") #"Decay mode of the spectator tau:  0=had, 1=e, 2=muon"
        self.out.branch("Gen_wIdx", "I") #"Idx to the last copy of the W from the taustar decay in GenPart"
        self.out.branch("Gen_wDM", "I") #"Decay mode of W. 0 = hadronic, 1=electron, 2=muon, 3=tau. -1 default"
        self.out.branch("Gen_wDau1Idx", "I") #"Idx to GenPart of the higher pT daughter of the W if wDM==0 or the charged lepton if wDM>0. -1 default"
        self.out.branch("Gen_wDau2Idx", "I") #"Idx to GenPart of the lower pT daughter of the W if wDM!=0 or the nu if wDM>0. -1 default"
        self.out.branch("Gen_wGenAK8Idx", "I") #"Idx to GenJetAK8 collection jet matching the W from taustar if wDM == 0"
        self.out.branch("Gen_nuIdx", "I") #"Idx to the last instance of the Nu in GenPart"
        self.out.branch("Gen_tauNuMET_pt", "F") #"MET from spectator tau and nu decays"
        self.out.branch("Gen_tauNuMET_eta", "F") #"MET from spectator tau and nu decays"
        self.out.branch("Gen_tauNuMET_phi", "F") #"MET from spectator tau and nu decays"
        self.out.branch("Gen_wMET_pt", "F") #"MET from W decay"
        self.out.branch("Gen_wMET_eta", "F") #"MET from W decay"
        self.out.branch("Gen_wMET_phi", "F") #"MET from W decay"
        self.out.branch("Gen_totMET_pt", "F") #"MET from W decay, spectator tau decay, and nu"
        self.out.branch("Gen_totMET_eta", "F") #"MET from W decay, spectator tau decay, and nu"
        self.out.branch("Gen_totMET_phi", "F") #"MET from W decay, spectator tau decay, and nu"
        self.out.branch("Gen_dr_tauW", "F") #"DeltaR between tau and W"
        self.out.branch("Gen_dr_tauTotMET", "F") #"DeltaR between the spectator tau and the total MET"
        self.out.branch("Gen_dr_wTotMET", "F") #"DeltaR between the the W and the total MET"
        self.out.branch("Gen_dr_tauNuMETTotMET", "F") #"DeltaR between the MET from the spectator tau + nu and the total MET"


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        
        #Indices default to -1
        tsIdx = -1
        tauIdx = -1
        tauDM = -1
        wIdx = -1
        nuIdx = -1
        wDau1Idx = -1
        wDau2Idx = -1
        wDM = -1
        wGenAK8Idx = -1
        tauNuMET_pt = -999
        tauNuMET_eta = -999
        tauNuMET_phi = -999
        wMET_pt = -999
        wMET_eta = -999
        wMET_phi = -999
        totMET_pt = -999
        totMET_eta = -999
        totMET_phi = -999
        dr_tauW = -999
        dr_tauTotMET = -999
        dr_wTotMET = -999
        dr_tauNuMETTotMET = -999

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
                    decayChain = getDecayChain(tauIdx, genParts)
                    for partIdx in decayChain: #Check for electron/muon or their neutrinos in decay chain
                        if abs(genParts[partIdx].pdgId) == 11 or abs(genParts[partIdx].pdgId) == 12:
                            tauDM = 1
                            break
                        elif abs(genParts[partIdx].pdgId) == 13 or abs(genParts[partIdx].pdgId) == 14:
                            tauDM = 2
                            break
                    if tauDM < 0: #If we didnt find an e/mu or nu_e/nu_mu must be a hadronic decay
                        tauDM = 0
                
            elif abs(genPart.pdgId) == 12 or abs(genPart.pdgId) == 14 or abs(genPart.pdgId) == 16: #Neutrino
                prodChain = getProdChain(idx, genParts)
                if prodChainContains(prodChain, idx=tsIdx): # If this was the neutrino from the taustar decay
                    nuIdx = idx
            elif abs(genPart.pdgId) == 24: #Found a W
                prodChain = getProdChain(idx, genParts)
                if prodChainContains(prodChain, pdgID = 4000015): #If this W is the taustar decay product
                    wIdx = idx
                    decayChain = getDecayChain(wIdx, genParts)

                    if len(decayChain) != 2:
                        print("WARNING: Number of W decay chain particles is != 2. Will write -1 indices for this event")
                    else:
                        if abs(event.GenPart_pdgId[decayChain[0]]) <= 6: #Hadronic
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
                        elif abs(event.GenPart_pdgId[decayChain[0]]) == 11: #el
                            wDM = 1
                            wDau1Idx = decayChain[0]
                            wDau2Idx = decayChain[1]
                        elif abs(event.GenPart_pdgId[decayChain[0]]) == 12: 
                            wDM = 1
                            wDau1Idx = decayChain[1]
                            wDau2Idx = decayChain[0]
                        elif abs(event.GenPart_pdgId[decayChain[0]]) == 13: #mu
                            wDM = 2
                            wDau1Idx = decayChain[0]
                            wDau2Idx = decayChain[1]
                        elif abs(event.GenPart_pdgId[decayChain[0]]) == 14: 
                            wDM = 2
                            wDau1Idx = decayChain[1]
                            wDau2Idx = decayChain[0]
                        elif abs(event.GenPart_pdgId[decayChain[0]]) == 15: #tau
                            wDM = 3
                            wDau1Idx = decayChain[0]
                            wDau2Idx = decayChain[1]
                        elif abs(event.GenPart_pdgId[decayChain[0]]) == 16: 
                            wDM = 3
                            wDau1Idx = decayChain[1]
                            wDau2Idx = decayChain[0]
            #End W
            if tsIdx >= 0 and tauIdx >=0 and wIdx>= 0 and nuIdx >= 0: #Found the interesting particles, do calculations now

                #MET calculations
                decayChain_tau = getDecayChain(tauIdx, genParts)
                invisbleParticles = []
                #First calculate total spectator tau + neutrino MET
                invisbleParticles.append(genParts[nuIdx])
                for partIdx in decayChain_tau: 
                    if abs(genParts[partIdx].pdgId) == 12 or abs(genParts[partIdx].pdgId) == 14 or abs(genParts[partIdx].pdgId) == 16:
                        invisbleParticles.append(genParts[partIdx]) 

                tauNuMET = invisbleParticles[0].p4()
                for inPrtIdx, invisPart in enumerate(invisbleParticles):
                    if inPrtIdx == 0:
                        continue
                    tauNuMET = tauNuMET + invisPart.p4()
                tauNuMET_eta = tauNuMET.Eta()
                tauNuMET_phi = tauNuMET.Phi()
                tauNuMET_pt = tauNuMET.Pt()

                #Now calculate MET from W
                if wDM > 0:
                    invisbleParticles = []
                    decayChain_w = getDecayChain(wIdx, genParts)
                    for i, partIdx in enumerate(decayChain_w):
                        if abs(genParts[partIdx].pdgId) == 12 or abs(genParts[partIdx].pdgId) == 14 or abs(genParts[partIdx].pdgId) == 16:
                            invisbleParticles.append(genParts[partIdx])
                        elif abs(genParts[partIdx].pdgId) == 15: #If the W decayed to tauNu then we also need to search for invisible particles from tau decay
                            decayChain_wTau = getDecayChain(partIdx, genParts)
                            for tauPartIdx in decayChain_wTau:
                                if abs(genParts[tauPartIdx].pdgId) == 12 or abs(genParts[tauPartIdx].pdgId) == 14 or abs(genParts[tauPartIdx].pdgId) == 16:
                                    invisbleParticles.append(genParts[tauPartIdx])

                    wMET = invisbleParticles[0].p4()
                    for inPrtIdx, invisPart in enumerate(invisbleParticles):
                        if inPrtIdx == 0:
                            continue
                        wMET = wMET + invisPart.p4()

                    wMET_eta = wMET.Eta()
                    wMET_phi = wMET.Phi()
                    wMET_pt = wMET.Pt()

                    #Total MET from interesting particles is tau + nu + W MET
                    totMET = tauNuMET + wMET
                    totMET_pt = totMET.Pt()
                    totMET_eta = totMET.Eta()
                    totMET_phi = totMET.Phi()
                else:
                    totMET = tauNuMET
                    totMET_pt = tauNuMET.Pt()
                    totMET_eta = tauNuMET.Eta()
                    totMET_phi = tauNuMET.Phi()

                #DeltaR calculations
                tau = genParts[tauIdx]
                w = genParts[wIdx]
                dr_tauW = tau.DeltaR(w)
                dr_tauTotMET = tau.DeltaR(totMET)
                dr_wTotMET = w.DeltaR(totMET)
                dr_tauNuMETTotMET = tauNuMET.DeltaR(totMET)

                break #Since we are done finding/calculating
        #End GenPart looping

        self.out.fillBranch("Gen_tsIdx", tsIdx)
        self.out.fillBranch("Gen_tauIdx", tauIdx)
        self.out.fillBranch("Gen_tauDM", tauDM)
        self.out.fillBranch("Gen_wIdx", wIdx)
        self.out.fillBranch("Gen_wDM", wDM)
        self.out.fillBranch("Gen_wDau1Idx", wDau1Idx)
        self.out.fillBranch("Gen_wDau2Idx", wDau2Idx)
        self.out.fillBranch("Gen_wGenAK8Idx", wGenAK8Idx)
        self.out.fillBranch("Gen_tauNuMET_pt", tauNuMET_pt)
        self.out.fillBranch("Gen_tauNuMET_eta", tauNuMET_eta)
        self.out.fillBranch("Gen_tauNuMET_phi", tauNuMET_phi)
        self.out.fillBranch("Gen_wMET_pt", wMET_pt)
        self.out.fillBranch("Gen_wMET_eta", wMET_eta)
        self.out.fillBranch("Gen_wMET_phi", wMET_phi)
        self.out.fillBranch("Gen_totMET_pt", totMET_pt)
        self.out.fillBranch("Gen_totMET_eta", totMET_eta)
        self.out.fillBranch("Gen_totMET_phi", totMET_phi)
        self.out.fillBranch("Gen_dr_tauW", dr_tauW)
        self.out.fillBranch("Gen_dr_tauTotMET", dr_tauTotMET)
        self.out.fillBranch("Gen_dr_wTotMET", dr_wTotMET)
        self.out.fillBranch("Gen_dr_tauNuMETTotMET", dr_tauNuMETTotMET)

        return True

# -----------------------------------------------------------------------------------------------------------------------------            

genProducerWNuConstr = lambda: GenProducerWNu()


files = ["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/WNu/taustarToWNu_m1000_2018.root"]
p = PostProcessor(".", files, cut="1>0", branchsel=None, modules=[genProducerWNuConstr()] )
p.run()
