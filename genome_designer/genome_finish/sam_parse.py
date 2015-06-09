"""
SAM parsing helper functions taken from lumpy's extractSplitReads
authored by Ira Hall under the MIT license
"""

import re


class SAM(object):
    """
    __very__ basic class for SAM input.
    """
    def __init__(self, samList = []):
        if len(samList) > 0:
            self.query    = samList[0]
            self.flag     = int(samList[1])
            self.ref      = samList[2]
            self.pos      = int(samList[3])
            self.mapq     = int(samList[4])
            self.cigar    = samList[5]
            self.matRef   = samList[6]
            self.matePos  = int(samList[7])
            self.iSize    = int(samList[8])
            self.seq      = samList[9]
            self.qual     = samList[10]
            self.tags     = samList[11:]#tags is a list of each tag:vtype:value sets
            self.valid    = 1
        else:
            self.valid = 0
            self.query = 'null'

    def extractTagValue(self, tagID):
        for tag in self.tags:
            tagParts = tag.split(':', 2);
            if (tagParts[0] == tagID):
                if (tagParts[1] == 'i'):
                    return int(tagParts[2]);
                elif (tagParts[1] == 'H'):
                    return int(tagParts[2],16);
                return tagParts[2];
        return None;
    
#-----------------------------------------------
cigarPattern = '([0-9]+[MIDNSHP])'
cigarSearch = re.compile(cigarPattern)
atomicCigarPattern = '([0-9]+)([MIDNSHP])'
atomicCigarSearch = re.compile(atomicCigarPattern)

def extractCigarOps(cigar,flag):
    if (cigar == "*"):
        cigarOps = []
    elif (flag & 0x0010):
        cigarOpStrings = cigarSearch.findall(cigar)
        cigarOps = []
        for opString in cigarOpStrings:
            cigarOpList = atomicCigarSearch.findall(opString)
            # "struct" for the op and it's length
            cigar = cigarOp(cigarOpList[0][0], cigarOpList[0][1])
            # add to the list of cigarOps
            cigarOps.append(cigar)
            cigarOps = cigarOps
        cigarOps.reverse()
        ##do in reverse order because negative strand##
    else:
        cigarOpStrings = cigarSearch.findall(cigar)
        cigarOps = []
        for opString in cigarOpStrings:
            cigarOpList = atomicCigarSearch.findall(opString)
            # "struct" for the op and it's length
            cigar = cigarOp(cigarOpList[0][0], cigarOpList[0][1])
            # add to the list of cigarOps
            cigarOps.append(cigar)
    return(cigarOps)

def calcQueryPosFromCigar(cigarOps):
    qsPos = 0
    qePos = 0
    qLen  = 0
    # if first op is a H, need to shift start position 
    # the opPosition counter sees if the for loop is looking at the first index of the cigar object    
    opPosition = 0  
    for cigar in cigarOps:
        if opPosition == 0 and (cigar.op == 'H' or cigar.op == 'S'):
            qsPos += cigar.length
            qePos += cigar.length
            qLen  += cigar.length
        elif opPosition > 0 and (cigar.op == 'H' or cigar.op == 'S'):
            qLen  += cigar.length
        elif cigar.op == 'M' or cigar.op == 'I':
            qePos += cigar.length
            qLen  += cigar.length
            opPosition += 1
    d = queryPos(qsPos, qePos, qLen);
    return d

class cigarOp(object):
    """
    sturct to store a discrete CIGAR operations
    """
    def __init__(self, opLength, op):
        self.length = int(opLength)
        self.op     = op

class queryPos(object):
    """
    struct to store the start and end positions of query CIGAR operations
    """
    def __init__(self, qsPos, qePos, qLen):
        self.qsPos = int(qsPos)
        self.qePos = int(qePos)
        self.qLen  = int(qLen)


def calcQueryOverlap(s1,e1,s2,e2):
    o = 1 + min(e1, e2) - max(s1, s2)
    return max(0, o)
