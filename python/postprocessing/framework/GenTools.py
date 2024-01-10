#Tools for working with GenPart particles

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