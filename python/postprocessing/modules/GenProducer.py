#Producer module to identify and store indexes to the most interesting gen particles in tau + taustar -> tautauZ decays
#This producer is only intented to be run over signal MC for that decay channel


from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel
from PhysicsTools.NanoAODTools.postprocessing.utils.GenTools import getDecayChain, getProdChain, prodChainContains
from PhysicsTools.NanoAODTools.postprocessing.utils.Tools import deltaR, isBetween

from ROOT import TH1F, TLorentzVector


# -----------------------------------------------------------------------------------------------------------------------------

class GenProducer(Module):

    def __init__(self, era):
        self.era = era
        self.writeHistFile = True
    
    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)
        self.h_genAk4Mass = TH1F('h_genAk4Mass', 'Mass of GEN AK4 ;Mass [GeV];# of Jets', 75, 0, 150)
        self.addObject(self.h_genAk4Mass)
        self.h_genAk8Mass = TH1F('h_genAk8Mass', 'Mass of GEN AK8 Jets;Mass [GeV];# of Jets', 75, 0, 150)
        self.addObject(self.h_genAk8Mass)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Gen_tsIdx", "I") #"Idx to last copy of taustar in GenPart"
        self.out.branch("Gen_tsTauIdx", "I") #"Idx to the last copy of the tau decay product of the taustar in GenPart"
        self.out.branch("Gen_tsTauDM", "I") #"Decay mode of the ts tau:  0=had, 1=e, 2=muon"
        self.out.branch("Gen_tauIdx", "I") #"Idx to the last copy of the spectator tau produced alongside the taustar in GenPart"
        self.out.branch("Gen_tauDM", "I") #"Decay mode of the spectator tau:  0=had, 1=e, 2=muon
        self.out.branch("Gen_tsTau_visInvDR", "F") #"DeltaR between the visible & invisible components of the tau from the taustar"
        self.out.branch("Gen_tau_visInvDR", "F") #"DeltaR between the visible & invisible components of the spectator tau"
        self.out.branch("Gen_tsTauFid", "O") #"True if the tsTau passes fiducial cuts"
        self.out.branch("Gen_tauFid", "O") #"True if the tau passes fiducial cuts"
        self.out.branch("Gen_zIdx", "I") #"Idx to the last copy of the Z from the taustar decay in GenPart"
        self.out.branch("Gen_zDau1Idx", "I") #"Idx to the first daughter of the Z at zIdx in GenPart"
        self.out.branch("Gen_zDau2Idx", "I") #"Idx to the second daughter of the Z at zIdx in GenPart"
        self.out.branch("Gen_zDauFid", "O") #"True if the (gen) Z daughters pass the appropriate fiducial cuts"
        self.out.branch("Gen_zD1RecIdx", "I") #"The idx of the matching Electron/Muon matching to the zDau1Idx GenPart if zDM==1/2"
        self.out.branch("Gen_zD2RecIdx", "I") #"The idx of the matching Electron/Muon matching to the zDau2Idx GenPart if zDM==1/2"
        self.out.branch("Gen_zDM","I") #"0 = hadronic, 1=electrons, 2=muons, 3=taus, 4=invisible"
        self.out.branch("Gen_zGenAK8Idx", "I") #"Idx to GenJetAK8 collection of a possible jet matching the Z from taustar" 
        self.out.branch("Gen_zGenAK4Idx", "I") #"Idx to GenJet collection of a possible jet matching the Z from the taustar"
        self.out.branch("Gen_zRecAK8Idx", "I") #"Idx to FatJet collection of reco AK8 matching the best GenAK8 jet"
        self.out.branch("Gen_zRecAK4Idx", "I") #"Idx to Jet collection of reco AK4 matching the best GenAK4 jet"
        self.out.branch("Gen_zSubJetIdx1", "I") #"Idx to higher pt SubGenJetAK8 collection of subJet matching the zGenAK8Idx jet"
        self.out.branch("Gen_zSubJetIdx2", "I") #"Idx to lower pt SubGenJetAK8 collection of subJet matching the zGenAK8Idx jet"
        self.out.branch("Gen_isCand", "O") #"True if event could be reco'd i.e. z,tau,tsTau DMs good, all fiducial, etc" 
        self.out.branch("Gen_ch", "I") #"0=tautau, 1=etau, 2=mutau. -1 otherwise"
        
        self.out.branch("Gen_tausMET_pt", "F") #"pT of MET from the spectator tau and taustar tau"
        self.out.branch("Gen_tausMET_eta", "F") #"eta of MET from the spectator tau and taustar tau"
        self.out.branch("Gen_tausMET_phi", "F") #"phi of MET from the spectator tau and taustar tau"
        self.out.branch("Gen_totMET_pt", "F") #"pT of MET from the spectator tau and taustar tau, and Z decay chains"
        self.out.branch("Gen_totMET_eta", "F") #"eta of MET from the spectator tau, taustar tau, and Z decay chains"
        self.out.branch("Gen_totMET_phi", "F") #"phi of MET from the spectator tau, taustar tau, and Z decay chains"
        self.out.branch("Gen_dr_tsTauTau", "F") #"DeltaR between tsTau and Tau"
        self.out.branch("Gen_dr_tsTauZ", "F") #"DeltaR between tsTau and Z"
        self.out.branch("Gen_dr_tauZ", "F") #"DeltaR between tau and Z"
        self.out.branch("Gen_dr_zDaus", "F") #"DeltaR between the two Z daughters"
        self.out.branch("Gen_dr_tsTauTausMET", "F") #"DeltaR between tsTau and MET from the 2 taus" 
        self.out.branch("Gen_dr_tauTausMET", "F") #"DeltaR between tau and MET from the 2 taus" 
        self.out.branch("Gen_dr_zTausMET", "F") #"DeltaR between the Z and MET from the 2 taus" 
        self.out.branch("Gen_dr_tsTauTotMET", "F") #"DeltaR between tsTau and total MET of interesting particles" 
        self.out.branch("Gen_dr_tauTotMET", "F") #"DeltaR between tau and total MET of interesting particles" 
        self.out.branch("Gen_dr_zTotMET", "F") #"DeltaR between the Z and total MET of interesting particles" 
        self.out.branch("Gen_dr_tausMETTotMET", "F") #"DeltaR between MET from the 2 taus and total MET of interesting particles" 

    def analyze(self, event):
        
        #Indices default to -1, floats to -999
        tsIdx = -1
        tsTauIdx = -1
        tsTauDM = -1
        tauIdx = -1
        tauDM = -1
        #tsTau_vis_pt = -999
        #tsTau_vis_eta = -999
        #tsTau_vis_phi = -999
        #tsTau_vis_mass = -999
        #tsTau_inv_pt = -999 
        #tsTau_inv_eta = -999
        #tsTau_inv_phi = -999
        tsTau_visInvDR = -999
        #tau_vis_pt = -999
        #tau_vis_eta = -999
        #tau_vis_phi = -999
        #tau_vis_mass = -999
        #tau_inv_pt = -999 
        #tau_inv_eta = -999
        #tau_inv_phi = -999
        tau_visInvDR = -999
        
        tsTauFid = False
        tauFid = False
        zIdx = -1
        zDau1Idx = -1
        zDau2Idx = -1
        zD1RecIdx = -1
        zD2RecIdx = -1
        zDauFid = False
        zDauKin = False
        zDM = -1
        zGenAK8Idx = -1
        zGenAK4Idx = -1
        recAK4JetIdx = -1
        recAK8JetIdx = -1
        subGenJet1Idx = -1
        subGenJet2Idx = -1
        isCand = False
        dr_tsTauTau = -999
        dr_tsTauZ = -999
        dr_tauZ = -999
        dr_zDaus = -999
        tausMET_pt = -999
        tausMET_eta = -999
        tausMET_phi = -999
        totMET_pt = -999
        totMET_eta = -999
        totMET_phi = -999
        dr_tsTauTausMET = -999
        dr_tauTausMET = -999
        dr_zTausMET = -999
        dr_tsTauTotMET = -999
        dr_tauTotMET = -999
        dr_zTotMET = -999
        dr_tausMETTotMET = -999
        ch = -1


        genParts = Collection(event, "GenPart")

        foundAll = False
        for idx, genPart in enumerate(genParts):
            #We only care about the last copy of relevant particles
            if not genPart.statusflag('isLastCopy'):
                continue
            
            if abs(genPart.pdgId) == 4000015: #Found the taustar
                tsIdx = idx
            elif abs(genPart.pdgId) == 15: #Found a tau
                prodChain = getProdChain(idx, genParts)
                if prodChainContains(prodChain, pdgID = 4000015): #If this tau is a taustar decay product
                    tsTauIdx = idx
                    if self.era == 2:
                        tsTauFid = abs(genPart.eta) < 2.3
                    elif self.era == 3:
                        tsTauFid = abs(genPart.eta) < 2.5
                    decayChain = getDecayChain(tsTauIdx, genParts)
                    for partIdx in decayChain: #Check for electron/muon or their neutrinos in decay chain
                        if abs(genParts[partIdx].pdgId) == 11:
                            tsTauDM = 1
                        elif abs(genParts[partIdx].pdgId) == 12:
                            tsTauDM = 1
                        
                        elif abs(genParts[partIdx].pdgId) == 13 or abs(genParts[partIdx].pdgId) == 14:
                            tsTauDM = 2
                            break
                    if tsTauDM < 0: #If we didnt find an e/mu or nu_e/nu_mu must be a hadronic decay
                        tsTauDM = 0
                elif prodChainContains(prodChain, idx = 0): # If this was the tau produced in the CI with the taustar (always idx 0 in GenParts)
                    tauIdx = idx
                    if self.era == 2:
                        tauFid = abs(genPart.eta) < 2.3
                    elif self.era == 3:
                        tauFid = abs(genPart.eta) < 2.5
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

            elif abs(genPart.pdgId) == 23: #Found a Z
                prodChain = getProdChain(idx, genParts)
                if prodChainContains(prodChain, pdgID = 4000015): #If this Z is the taustar decay product
                    zIdx = idx
                    decayChain = getDecayChain(zIdx, genParts)
                    if len(decayChain) != 2:
                        print("WARNING: Number of Z decay chain particles is != 2. Will write -1 indices for this event")
                    else:
                        if event.GenPart_pt[decayChain[0]] > event.GenPart_pt[decayChain[1]]: #Hightest pt daughter gets idx 1 label
                            zDau1Idx = decayChain[0]
                            zDau2Idx = decayChain[1]
                        else:
                            zDau1Idx = decayChain[1]
                            zDau2Idx = decayChain[0]

                        if abs(event.GenPart_pdgId[zDau1Idx]) <= 6: #Quarks -> hadronic decay mode
                            zDM = 0
                        elif abs(event.GenPart_pdgId[zDau1Idx]) == 11: #Z->ee
                            zDM = 1
                            zDau1 = genParts[zDau1Idx]
                            zDau2 = genParts[zDau2Idx]
                            zDauFid = abs(zDau1.eta) < 2.5 and (abs(zDau1.eta) > 1.566 or abs(zDau1.eta) < 1.444)
                            zDauFid = zDauFid and abs(zDau2.eta) < 2.5 and (abs(zDau2.eta) > 1.566 or abs(zDau2.eta) < 1.444)

                            zDauKin = zDau1.pt > 20.0 and zDau2.pt > 20.0
                            
                            #Find the matching reco electron:
                            electrons = Collection(event, "Electron")
                            for elIdx, el in enumerate(electrons):
                                if el.genPartIdx == zDau1Idx:
                                    zD1RecIdx = elIdx
                                elif el.genPartIdx == zDau2Idx:
                                    zD2RecIdx = elIdx
                                
                        elif abs(event.GenPart_pdgId[zDau1Idx]) == 13: #Z->mumu
                            zDM = 2
                            zDau1 = genParts[zDau1Idx]
                            zDau2 = genParts[zDau2Idx]
                            zDauFid = zDau1.eta < 2.4 and zDau2.eta < 2.4    
                            zDauKin = zDau1.pt > 15.0 and zDau2.pt > 15.0
                        
                            #Find the matching reco electron:
                            muons = Collection(event, "Muon")
                            for muIdx, mu in enumerate(muons):
                                if mu.genPartIdx == zDau1Idx:
                                    zD1RecIdx = muIdx
                                elif mu.genPartIdx == zDau2Idx:
                                    zD2RecIdx = muIdx
                                
                        elif abs(event.GenPart_pdgId[zDau1Idx]) == 15: #Z->tautau
                            zDM = 3
                        if abs(event.GenPart_pdgId[zDau1Idx]) == 12 or abs(event.GenPart_pdgId[zDau1Idx]) == 14 or abs(event.GenPart_pdgId[zDau1Idx]) == 16:
                            zDM = 4 #Z->invisible

                        #Check if the Z corresponds to a Gen AK8 jet 
                        #print("Gen_zDM = " + str(zDM))
                        theZ = genParts[zIdx]
                        genJetAK8s = Collection(event, "GenJetAK8")
                        nMatches_AK8 = 0
                        for jetIdx, jet in enumerate(genJetAK8s):
                            if abs(jet.eta) < 2.5 and jet.pt > 100:
                                self.h_genAk8Mass.Fill(jet.mass)
                            if jet.DeltaR(theZ) < 0.1 and abs(jet.pt - theZ.pt) <= 0.1*theZ.pt:# and abs(jet.mass - 91.18) < 9:
                                nMatches_AK8 += 1
                                if zGenAK8Idx >= 0:
                                    if jet.DeltaR(theZ) < genJetAK8s[zGenAK8Idx].DeltaR(theZ):
                                        zGenAK8Idx = jetIdx
                                else:
                                    zGenAK8Idx = jetIdx
                                
                        #print("Number of AK8 Jet Matches = " + str(nMatches_AK8))

                        genJetAK4s = Collection(event, "GenJet")
                        zGenAK4Idx = -1
                        nMatches_AK4 = 0
                        for jetIdx, jet in enumerate(genJetAK4s):
                            if abs(jet.eta) < 2.5 and jet.pt > 100:
                                self.h_genAk4Mass.Fill(jet.mass)
                            if jet.DeltaR(theZ) < 0.1 and abs(jet.pt - theZ.pt) <= 0.1*theZ.pt:# and abs(jet.mass - 91.18) < 9:
                                nMatches_AK4 += 1
                                if zGenAK4Idx >= 0:
                                    if jet.DeltaR(theZ) < genJetAK4s[zGenAK4Idx].DeltaR(theZ):
                                        zGenAK4Idx = jetIdx
                                else:
                                    zGenAK4Idx = jetIdx
                        #print("Number of AK4 Jet Matches = " + str(nMatches_AK4))

                        #Find the reco jets corresponding to the GEN jets 
                        if zGenAK8Idx >= 0:
                            ak8Jets = Collection(event, "FatJet")
                            for jetIdx, jet in enumerate(ak8Jets):
                                if jet.genJetAK8Idx == zGenAK8Idx:
                                    recAK8JetIdx = jetIdx
                                    break
                        if zGenAK4Idx >= 0:
                            ak4Jets = Collection(event, "Jet")
                            for jetIdx, jet in enumerate(ak4Jets):
                                if jet.genJetIdx == zGenAK4Idx:
                                    recAK4JetIdx = jetIdx
                                    break
                            
                        if zDM == 0 and zGenAK8Idx >=0:
                            subGenJets = Collection(event, "SubGenJetAK8")
                            zGenAK8 = genJetAK8s[zGenAK8Idx]
                            for sJIdx, subJet in enumerate(subGenJets):
                                if deltaR(subJet.eta, subJet.phi, zGenAK8.eta, zGenAK8.phi) < 0.5:
                                    if subGenJet1Idx < 0:
                                        subGenJet1Idx = sJIdx
                                    elif subGenJet2Idx < 0:
                                        subGenJet2Idx = sJIdx
                                    else:
                                        print('WARNING: Found more than two "matching" subGenJets to zGenAK8!"')
                            if subGenJet1Idx < 0 or subGenJet2Idx < 0:
                                print("WARNING: Did not find two subGenJets matching to zGenAK8")
                            else:
                                #Make first idx higher pt for consistency with other multi-idx convention
                                if subGenJets[subGenJet1Idx].pt < subGenJets[subGenJet2Idx].pt:
                                    temp = subGenJet1Idx
                                    subGenJet1Idx = subGenJet2Idx
                                    subGenJet2Idx = temp


            foundAll = (tsTauIdx >= 0) and (tauIdx >= 0) and (zIdx >= 0) and (zDau1Idx >=0) and (zDau2Idx >= 0)
            if foundAll:

                isCand = zDM >= 0 and zDM <= 2 #Z->had/ee/mumu
                isCand = isCand and (tsTauDM == 0 or tauDM == 0) #etau,mutau,tautau
                isCand = isCand and (tsTauFid and tauFid) #Taus within fiducial region
                if zDM == 1 or zDM == 2:
                    isCand = isCand and zDauFid
                    isCand = isCand and (zD1RecIdx >= 0 and zD2RecIdx >= 0) #Matching reco particles found
                    isCand = isCand and zDauKin
                #One more requirement is below after MET calculation. Namely, that MET is between the two taus

                if tauDM == 0 and tsTauDM == 0:
                    ch = 0
                elif (tauDM == 0 and tsTauDM == 1) or (tauDM == 1 and tsTauDM == 0):
                    ch = 1
                elif (tauDM == 0 and tsTauDM == 2) or (tauDM == 2 and tsTauDM == 0):
                    ch = 2
                    
                #Calculate the total MET coming from the decay of the two taus
                decayChain_tsTau = getDecayChain(tsTauIdx, genParts)
                decayChain_tau = getDecayChain(tauIdx, genParts)
                invisbleParticles = []
                tau_inv = TLorentzVector()
                tsTau_inv = TLorentzVector()
                for partIdx in decayChain_tsTau:
                    if abs(genParts[partIdx].pdgId) == 12 or abs(genParts[partIdx].pdgId) == 14 or abs(genParts[partIdx].pdgId) == 16:
                        invisbleParticles.append(genParts[partIdx])
                        tsTau_inv = tsTau_inv + genParts[partIdx].p4()
                for partIdx in decayChain_tau:
                    if abs(genParts[partIdx].pdgId) == 12 or abs(genParts[partIdx].pdgId) == 14 or abs(genParts[partIdx].pdgId) == 16:
                        invisbleParticles.append(genParts[partIdx]) 
                        tau_inv = tau_inv + genParts[partIdx].p4()
                if len(invisbleParticles) < 2: #Must have at least two neutrinos from the two taus
                    print("ERROR: In GenProducerZTau, could not find at least two neutrinos from tau decays")
                else:
                    met = invisbleParticles[0].p4()
                    for inPrtIdx, invisPart in enumerate(invisbleParticles):
                        if inPrtIdx == 0:
                            continue
                        met = met + invisPart.p4()
                    tausMET_eta = met.Eta()
                    tausMET_phi = met.Phi()
                    tausMET_pt = met.Pt()
                
                tau_vis = genParts[tauIdx].p4() - tau_inv
                tau_visInvDR = tau_vis.DeltaR(tau_inv)
                tsTau_vis = genParts[tsTauIdx].p4() - tsTau_inv
                tsTau_visInvDR = tsTau_vis.DeltaR(tsTau_inv)

                #Z->tautau and Z->invisible also contribute to MET from interesting particles in the event
                totMet = met
                if zDM == 3: #Z->tautau we need the daughters of the two taus
                    tauDaus = getDecayChain(zDau1Idx, genParts)
                    tauDaus.extend(getDecayChain(zDau2Idx, genParts))
                    for zTauDauIdx in tauDaus:
                        if abs(genParts[zTauDauIdx].pdgId) == 12 or abs(genParts[zTauDauIdx].pdgId) == 14 or abs(genParts[zTauDauIdx].pdgId) == 16:
                            totMet = totMet + genParts[zTauDauIdx].p4()
                elif zDM == 4:#Z->invisible we just need zDaughters 
                    totMet = totMet + genParts[zDau1Idx].p4()
                    totMet = totMet + genParts[zDau2Idx].p4()
                
                totMET_eta = totMet.Eta()
                totMET_phi = totMet.Phi()
                totMET_pt = totMet.Pt()


                isCand = isCand and isBetween(genParts[tsTauIdx].phi, genParts[tauIdx].phi, totMET_phi)
                
                
                #Now do DeltaR calculations
                tsTau = genParts[tsTauIdx]
                tau = genParts[tauIdx]
                z = genParts[zIdx]
                zDau1 = genParts[zDau1Idx]
                zDau2 = genParts[zDau2Idx]

                dr_tsTauTau = tsTau.DeltaR(tau)
                dr_tsTauZ = tsTau.DeltaR(z)
                dr_tauZ = tau.DeltaR(z)
                dr_zDaus = zDau1.DeltaR(zDau2)
                dr_tsTauTotMET = tsTau.DeltaR(totMet)
                dr_tauTotMET = tau.DeltaR(totMet)
                dr_zTotMET = tau.DeltaR(totMet)
                dr_tsTauTausMET = tsTau.DeltaR(met)
                dr_tauTausMET = tau.DeltaR(met)
                dr_zTausMET = z.DeltaR(met)
                dr_tausMETTotMET = met.DeltaR(totMet)

                break #Since we've found all interesting particles and done all calculations
        
        if not foundAll:
            print("WARNING: In GenProducer: All interesting gen particles were not found!")

        #Write output to tree
        self.out.fillBranch("Gen_tsIdx", tsIdx)
        self.out.fillBranch("Gen_tsTauIdx", tsTauIdx)
        self.out.fillBranch("Gen_tsTauDM", tsTauDM)
        self.out.fillBranch("Gen_tauIdx", tauIdx)
        self.out.fillBranch("Gen_tauDM", tauDM)
        self.out.fillBranch("Gen_tau_visInvDR", tau_visInvDR)
        self.out.fillBranch("Gen_tsTau_visInvDR", tsTau_visInvDR)
        self.out.fillBranch("Gen_tsTauFid", tsTauFid)
        self.out.fillBranch("Gen_tauFid", tauFid)
        self.out.fillBranch("Gen_zIdx", zIdx)
        self.out.fillBranch("Gen_zDau1Idx",zDau1Idx)
        self.out.fillBranch("Gen_zDau2Idx",zDau2Idx)
        self.out.fillBranch("Gen_zDauFid", zDauFid)
        self.out.fillBranch("Gen_zD1RecIdx", zD1RecIdx)
        self.out.fillBranch("Gen_zD2RecIdx", zD2RecIdx)
        self.out.fillBranch("Gen_zDM", zDM)
        self.out.fillBranch("Gen_zGenAK8Idx", zGenAK8Idx)
        self.out.fillBranch("Gen_zGenAK4Idx", zGenAK4Idx)
        self.out.fillBranch("Gen_zRecAK8Idx", recAK8JetIdx)
        self.out.fillBranch("Gen_zRecAK4Idx", recAK4JetIdx)
        self.out.fillBranch("Gen_zSubJetIdx1", subGenJet1Idx)
        self.out.fillBranch("Gen_zSubJetIdx2", subGenJet2Idx)
        self.out.fillBranch("Gen_isCand", isCand)
        self.out.fillBranch("Gen_ch", ch)
        self.out.fillBranch("Gen_tausMET_pt", tausMET_pt)
        self.out.fillBranch("Gen_tausMET_eta", tausMET_eta)
        self.out.fillBranch("Gen_tausMET_phi", tausMET_phi)
        self.out.fillBranch("Gen_totMET_pt", totMET_pt)
        self.out.fillBranch("Gen_totMET_eta", totMET_eta)
        self.out.fillBranch("Gen_totMET_phi", totMET_phi)
        self.out.fillBranch("Gen_dr_tsTauTau", dr_tsTauTau)
        self.out.fillBranch("Gen_dr_tsTauZ", dr_tsTauZ)
        self.out.fillBranch("Gen_dr_tauZ", dr_tauZ)
        self.out.fillBranch("Gen_dr_zDaus", dr_zDaus)
        self.out.fillBranch("Gen_dr_tsTauTausMET", dr_tsTauTausMET)
        self.out.fillBranch("Gen_dr_tauTausMET", dr_tauTausMET)
        self.out.fillBranch("Gen_dr_zTausMET", dr_zTausMET)
        self.out.fillBranch("Gen_dr_tsTauTotMET", dr_tsTauTotMET)
        self.out.fillBranch("Gen_dr_tauTotMET", dr_tauTotMET)
        self.out.fillBranch("Gen_dr_zTotMET", dr_zTotMET)
        self.out.fillBranch("Gen_dr_tausMETTotMET", dr_tausMETTotMET)

        return True

# -----------------------------------------------------------------------------------------------------------------------------            

genProducerConstr = lambda era: GenProducer(era)

# -----------------------------------------------------------------------------------------------------------------------------

#files = ["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/TauZ/taustarToTauZ_m3000_2018.root"]
#p = PostProcessor(".", files, cut="1>0", branchsel=None, modules=[genProducerZTauConstr()] )
#p.run()
