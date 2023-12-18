#Producer module to identify and store indexes to the most interesting gen particles (taustars + decay products, etc)


from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import PhysicsTools.NanoAODTools.postprocessing.framework.datamodel as datamodel


# -----------------------------------------------------------------------------------------------------------------------------

#Build a list of tuples of the form (idx, pdgId) recording the production chain in genParts
# idx : the idx to GenPart of the particle to find the chain of (will be last entry in returned list
# genParts : the GenPart collection object
def getProdChain(idx, genParts):

    part = genParts[idx]
    chain = [(idx, part.pdgId)] 
    mothIdx = part.genPartIdxMother
    while mothIdx >= 0:
        moth = genParts[mothIdx]
        chain.insert(0, (mothIdx, moth.pdgId))
        mothIdx = moth.genPartIdxMother

    return chain

# -----------------------------------------------------------------------------------------------------------------------------

def getDecayChain(idx, genParts):
    chain = []
    for dauIdx in range(idx + 1, genParts._len):
        dauPart = genParts[dauIdx]
        if dauPart.genPartIdxMother == idx:
            chain.append(dauIdx)
    return chain


# -----------------------------------------------------------------------------------------------------------------------------

# Check whether the production chain contains a particle matching idx and/xor pdgID
# prodChain : The production chain list as returned by getProdChain()
# idx : the idx to GenPart to try to match. If not specified, will not try to match indices
# pdgID : the pdgID to try to match. If not specified, will not try to match
# Returns a boolean specifying if a match was found
def prodChainContains(prodChain, idx = None, pdgID = None ):

    if idx == None and pdgID == None:
        print("ERROR: idx and pdgID are both None when trying to check production chain for matches")
        exit(1)

    found = False

    for step in prodChain:
        if idx != None and pdgID == None:
            found = found or (step[0] == idx)
        elif idx == None and pdgID != None:
            found = found or (abs(step[1]) == abs(pdgID))
        else:
            found = found or ( (step[0] == idx) and (abs(step[1]) == abs(pdgID)) )

        if found:
            return True

    return found

# -----------------------------------------------------------------------------------------------------------------------------

class GenProducer(Module):

    def __init__(self):
        pass
    
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Gen_tsIdx", "I") #Idx to last copy of taustar in GenPart
        self.out.branch("Gen_tsTauIdx", "I") #Idx to the last copy of the tau decay product of the taustar in GenPart
        self.out.branch("Gen_tauIdx", "I") #Idx to the last copy of the tau produced alongside the taustar in GenPart
        self.out.branch("Gen_zIdx", "I") #Idx to the last copy of the Z from the taustar decay in GenPart
        self.out.branch("Gen_zDau1Idx", "I") #Idx to the first daughter of the Z at zIdx in GenPart
        self.out.branch("Gen_zDau2Idx", "I") #Idx to the second daughter of the Z at zIdx in GenPart
        self.out.branch("Gen_zDM","I") #0 = hadronic, 1=electrons, 2=muons, 3=taus
        self.out.branch("Gen_dr_tsTauTau", "F") #DeltaR between tsTau and Tau
        self.out.branch("Gen_dr_tsTauZ", "F") # DeltaR between tsTau and Z
        self.out.branch("Gen_dr_tauZ", "F") #DeltaR between tau and Z
        self.out.branch("Gen_dr_zDaus", "F") #DeltaR between the two Z daughters

        #self.out.branch("Gen_ch","I") #1 = ETau, 2=MuTau, 3=TauTau

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass


    def analyze(self, event):
        
        #Indices default to -1
        tsIdx = -1
        tsTauIdx = -1
        tauIdx = -1
        zIdx = -1
        zDau1Idx = -1
        zDau2Idx = -1
        zDM = -1
        dr_tsTauTau = -999
        dr_tsTauZ = -999
        dr_tauZ = -999
        dr_zDaus = -999

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
                elif prodChainContains(prodChain, idx = 0): # If this was the tau produced in the CI with the taustar
                    tauIdx = idx
            elif abs(genPart.pdgId) == 23: #Found a Z
                prodChain = getProdChain(idx, genParts)
                if prodChainContains(prodChain, pdgID = 4000015): #If this Z is the taustar decay product
                    zIdx = idx
                    decayChain = getDecayChain(zIdx, genParts)
                    if len(decayChain) != 2:
                        print("WARNING: Number of Z decay chain particles is != 2. Will write -1 indices for this event")
                    else:
                        zDau1Idx = decayChain[0]
                        zDau2Idx = decayChain[1]
                        if abs(event.GenPart_pdgid[zDau1Idx]) <= 6:
                            zDM = 0
                        elif abs(event.GenPart_pdgid[zDau1Idx]) == 11:
                            zDM = 1
                        elif abs(event.GenPart_pdgid[zDau1Idx]) == 13:
                            zDM = 2
                        elif abs(event.GenPart_pdgid[zDau1Idx]) == 15:
                            zDM = 3
            foundAll = (tsTauIdx >= 0) and (tauIdx >= 0) and (zIdx >= 0) and (zDau1Idx >=0) and (zDau2Idx >= 0)
            if foundAll:
                break
        #End idx finding

        #Calculate DeltaR between the interesting particles
        if foundAll:
            tsTau = genParts[tsTauIdx]
            tau = genParts[tauIdx]
            z = genParts[zIdx]
            zDau1 = genParts[zDau1Idx]
            zDau2 = genParts[zDau2Idx]

            dr_tsTauTau = tsTau.DeltaR(tau)
            dr_tsTauZ = tsTau.DeltaR(z)
            dr_tauZ = tau.DeltaR(z)
            dr_zDaus = zDau1.DeltaR(zDau2)
        else:
            print("WARNING: In GenProducer: All interesting gen particles were not found. Will not calculate DeltaR's")

        #Write output to tree
        self.out.fillBranch("Gen_tsIdx", tsIdx)
        self.out.fillBranch("Gen_tsTauIdx", tsTauIdx)
        self.out.fillBranch("Gen_tauIdx", tauIdx)
        self.out.fillBranch("Gen_zIdx", zIdx)
        self.out.fillBranch("Gen_zDau1Idx",zDau1Idx)
        self.out.fillBranch("Gen_zDau2Idx",zDau2Idx)
        self.out.fillBranch("Gen_zDM", zDM)
        self.out.fillBranch("Gen_dr_tsTauTau", dr_tsTauTau)
        self.out.fillBranch("Gen_dr_tsTauZ", dr_tsTauZ)
        self.out.fillBranch("Gen_dr_tauZ", dr_tauZ)
        self.out.fillBranch("Gen_dr_zDaus", dr_zDaus)

        return True
                    
# -----------------------------------------------------------------------------------------------------------------------------            
            
genProducerConstr = lambda: GenProducer()


files = ["root://cmsxrootd.fnal.gov//store/user/bbarton/TaustarToTauTauZ/SignalMC/taustarToTauZ_m3000_2018.root"]
p = PostProcessor(".", files, cut="1>0", branchsel=None, modules=[genProducerConstr()] )
p.run()
