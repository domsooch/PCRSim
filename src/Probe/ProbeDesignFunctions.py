print 'Import ProbeDesignFunctions.py'
import os, sys 

parent_dir = os.path.dirname(__file__)
parent_dir_mod = parent_dir.replace(os.path.split(parent_dir)[1], '') # subtract the source file name
sys.path.append(parent_dir_mod)

import sys, os, copy
from math import log
from math import floor

import probe as PROBE

parent_dir = os.path.dirname(__file__)
parent_dir_mod = parent_dir.replace(os.path.split(parent_dir)[1], '') # subtract the source file name
sys.path.append(parent_dir_mod)

import utils.IO
import utils.FASTA as FA
import utils.ListManipulator as LM
import MolBio as MB


location = {}
DEBUG = None
#Types of FASTA Headers to Parse:
# (1) AY768562 (A/chicken/Korea/SNU0073/00  NA (6)  H9N2  yr:2000) L:262
# (2) >gi|9626945|ref|NC_001498.1| Measles virus, complete genome
# (3) AY768562_(A/chicken/Korea/SNU0073/00  NA (6)  H9N2  yr:2000) L:262
# (-1) gnl|UG|Gga#S19183375 Gallus gallus breast cancer 2, early onset(BRCA2),
#      mRNA /cds=(1,10194) /gb=NM_204276 /gi=45383585 /ug=Gga.1 /len=10310



DefaultHybScoreMatrix = {'HP':1000,
               'PD':10,
               'GC':10000,
               'rep':10000,
               'Tm':100,
               'CyclePenalty':1000}

DefaultPCRScoreMatrix = {'HP':1000,
               'PD':10,
               'GC':10000,
               'rep':10000,
               'Tm':100,
               'CyclePenalty':1000,
                'PCRGC':10000,
                'self-Anneal':1}
        
DefaultProbeFilterParameters = {'strict':None,	##irrelevant
                  'TmMin':68,		##<-- relevant
                  'MinLen':20,	##<-- relevant
                  'MaxLen':35,	##<-- relevant, max is 40
                  'MaxGC':65,		##irrelevant
                  'MinGC':35,		##irrelevant
                  'MaxHPTm':55,		##irrelevant: Hairpin max Tm
                  'MaxPalindrome':70,	##irrelevant: Palindrom Max Tm 
                  'MaxRepeats':6,		##irrelevant: Max number of repeats 3.1 is TTT repeats 3.2 is TATATA
                  'MaxSynthCycles':150,    ##relevant for synthesis
                  }

PCRPrimerTilerparams = {'type':1, 		##irrelevant
                  'strict':None,	##irrelevant
                  'parseLabelType':100,	##irrelevant Can be lambdaFunction
                  'ProbeType':'PCR',	##irrelevant
                  'TmMin':60,		##<-- relevant
                  'MaxSynthCycles':150, ##<---relevant for synthesis
                  'MinLen':19,	##<-- relevant
                  'MaxLen':28,	##<-- relevant, max is 40
                  'MaxGC':65,		##irrelevant
                  'MinGC':35,		##irrelevant
                  'half':None,	##Exclusive: Tiling parameter start half-way through previous probe
                  'overlap':None,	##Exclusive: Tiling parameter start x bases before end of previous probe
                  'startIncrement':1,	##Exclusive: Tiling parameter start 11 bases after start of previous probe
                  'MaxHPTm':55,		##irrelevant: Hairpin max Tm
                  'MaxPalindrome':55,	##irrelevant: Palindrom Max Tm 
                  'MaxRepeats':6,		##irrelevant: Max number of repeats 3.1 is TTT repeats 3.2 is TATATA
                  'CaptureEnd': 1,		##See Note below
              'ProbeStartNumber' : 0	##For labelling purposes only
            }

HybProbeTilerparams ={'type':1, 		##irrelevant
                  'strict':None,	##irrelevant
                  'parseLabelType':100,	##irrelevant Can be lambdaFunction
                  'ProbeType':'Hyb',	##irrelevant
                  'TmMin':72,		##<-- relevant
                  'MaxSynthCycles':150, ##<---relevant for synthesis
                  'MinLen':20,	##<-- relevant
                  'MaxLen':35,	##<-- relevant, max is 40
                  'MaxGC':65,		##irrelevant
                  'MinGC':35,		##irrelevant
                  'half':1,	##Exclusive: Tiling parameter start half-way through previous probe
                  'overlap':None,	##Exclusive: Tiling parameter start x bases before end of previous probe
                  'startIncrement':None,	##Exclusive: Tiling parameter start 11 bases after start of previous probe
                  'MaxHPTm':55,		##irrelevant: Hairpin max Tm
                  'MaxPalindrome':75,	##irrelevant: Palindrom Max Tm 
                  'MaxRepeats':6,		##irrelevant: Max number of repeats 3.1 is TTT repeats 3.2 is TATATA
                  'CaptureEnd': None,		##See Note below
              'ProbeStartNumber' : 0}	##For labelling purposes only






def RenderSeqs(SeqData, numProbesToTake):
    """RenderSeqs
    |(SeqData, numProbesToTake)
    |Output: list of fastaObjects
    |Returns the numProbesToTake best probes in alist of FilteredTiled seqs run with
    returnProbeQuality:ON
    |dependent on MB.FilterProbes()|"""
    ##LabelList = ['ProbeGood','HP','PD','repeatNum','GC','HPstruct','PDstruct','failreason']
    LabelLine = SeqData.pop(0)
    LabelDict = {}
    i = 0
    for obj in LabelLine:
        LabelDict[obj] = i
        i+=1
    i = 0;nexti = i + 1000;ProbesToProcess = len(SeqData)
    outProbeData = [] #[labelList]
    sortList = []
    for obj in SeqData:
        goodProbe = int(obj[LabelDict['ProbeGood']])
        repeatNum = int(obj[LabelDict['repeatNum']])
        SSTm = max(int(obj[LabelDict['HP']]),int(obj[LabelDict['PD']]))
        Tm = int(obj[LabelDict['Tm']])
        GCDiff = abs(50 - obj[LabelDict['GC']])
        if GCDiff > 15:
            GCDiff = 100
        SSscore = int((Tm-SSTm)/10)
        if SSscore <0:
            SSscore = 5
        if repeatNum < 3:
            RepeatScore = 0
        else:
            RepeatScore = (repeatNum -3)*1
        score = goodProbe*10000 - RepeatScore*1000 - SSscore*1000 - GCDiff*100
        sortList.append([i,score])
        i +=1
        if i > nexti:
            nexti = i + 1000
            print 'RenderSeqs has processed %i of %i probes in object SeqData.' %(i,ProbesToProcess)
        
    sortList.sort(lambda x,y:cmp(x[1],y[1]))
    numProbesToTake = numProbesToTake*-1
    chooseIndices = sortList[numProbesToTake:]
    outDB = []
    for obj in chooseIndices:
        i = obj[0]
        outProbeData.append(SeqData[i])
        outDB.append([SeqData[i][LabelDict['label']],SeqData[i][LabelDict['sequence']],SeqData[i][LabelDict['length']]])
    return outDB, outProbeData


def ParseFASTAHeader(sLine, FHtype):
    """ParseFASTAHeader
    |(sLine, FHtype)
    |Output:'AY768562'
    |Returns a simplified, truncated label given a FASTA label. This function
    will expect a type of label and return an output depending on what FHtype
    is specified.
    |Input: (FASTA label, 1-5)|"""
    if str(type(FHtype)).find('str')<>-1:
        return FHtype
    if str(type(FHtype)).find('function')<>-1:
        return FHtype(sLine)
    else:
        FHtype = int(FHtype)
    if FHtype == 1:
        if ' ' in sLine:
            label = sLine[0:sLine.find(' ')]
        else:
            label = sLine
    elif FHtype == 2:
        label = sLine[LM.FindNth(sLine,'|',3)+1:LM.FindNth(sLine,'|',4)]
    elif FHtype == 3:
        label = sLine[0:sLine.find('_')]
    elif FHtype == 4:
        label = sLine
    elif FHtype == 5:
        label = sLine[0:sLine.find(' ')]
    elif FHtype == -1:
        label = sLine[sLine.find('#S'):sLine.find(' ')]
    elif FHtype == 7:
        label = sLine.split('|')[3]
    else:
        label = sLine[0:FHtype].strip().replace(' ','')
    return label


def CyclesToSynth(seq, word = 'acgt', Verbose = 0):
    """CyclesToSynth
    |(seq,word = 'acgt')
    |Num_cycles (integer)
    |This function calculates the number of cycles required to synth oligo
    |
    |
    """
    seq = seq.lower()
    seq = list(seq)
    cycle = 0
    b = 2
    while len(seq):
        basemade = ''
        cycle +=1
        b +=1;b = b%4
        base = word[b]
        if not(seq[-1] in word):
            return None
        if seq[-1] == base:
            basemade = seq.pop(-1)
        if Verbose:
            print 'Cycle %i base %s %i basemade %s' %(cycle, base, b, basemade)
    return cycle

def CyclesCutOff(seq,label,max_Cycles,Word = 'acgt'):
    """CyclesCutOff
    |(seq,label,max_Cycles,word = 'acgt')
    |[label,seq]
    |This function does not allow a probe to be made that is > max_Cycles
    | The label returns the cut off and the missing bases
    |
    """
    #print Word

    if CyclesToSynth(seq,word = Word) <= max_Cycles:
        return [label,seq]
    else:
        Word = Word.upper()
        seq = seq.upper()
        seqLst = list(seq)
        seqLen = len(seqLst)
        b = 0
        cycle = 0
        returnSeq = ''
        returnSeqLen = 0
        while not(cycle > max_Cycles -1) and seqLen > 0:
            synthBase = seqLst[-1]
            cycle +=1
            testBase = Word[b]
            b += 1
            if b == 4:
                b = 0
            if testBase == synthBase:
                s = seqLst.pop()
                returnSeq = s + returnSeq
                seqLen = len(seqLst)
                returnSeqLen = len(returnSeq)
            if not(synthBase in Word) and synthBase <> '':
                return 0
            #print 'lenseqlst %i cycle %i seqLen %i' %(len(seqLst),cycle,seqLen)
        
        seq = "".join(seqLst)
        label = label + 'trunctoCy' + str(max_Cycles) + '_' + seq + '|'
        return [label,returnSeq]

def ProbeStats(inDB,UserDir,OFN='ns',seq_index = 1,return_AnnotList = None, Verbose = 1):
    """ProbeStats
    |(inDB,UserDir,OFN = 'ns',seq_index = 1,return_AnnotList = None, Verbose = 1)
    |histogramFileName,(optional: annotatedprobeList)
    |This function calculates length, Tm and %GC, and makes a histogram
    |input does not have to be fastaObj, but it helps
    |
    """
    report_interval = 1000
    header = ['Tm','Len','%GC','Cycles']
    TmHisto =  [0 for i in range(300)]
    LenHisto = [0 for i in range(300)]
    perGCHisto = [0 for i in range(300)]
    cycleHisto = [0 for i in range(300)] 
    sumTm = 0
    sumLen = 0
    sumGC = 0
    sumCycles = 0
    num_recs = len(inDB)
    report = 0
    for i in range(num_recs):
        if Verbose:
            report += 1
            if report == 0 or report == report_interval:
                print 'Processing record %i of %i' %(i,num_recs)
                report = 1
        seq = inDB[i][seq_index]    
        Tm = int(round(MB.Tm(seq)))
        length = len(seq)
        perGC = int(round(MB.perGC(seq)))
        cycles = CyclesToSynth(seq)
        if Tm> 100:
            print seq + '  ' + str(Tm)
        TmHisto[Tm] += 1
        LenHisto[length] += 1
        perGCHisto[perGC]  += 1
        cycleHisto[cycles] +=1
        sumTm += Tm
        sumLen += length
        sumGC += perGC
        sumCycles += cycles
        if return_AnnotList:
            outPiece = [Tm,length,perGC,cycles]
            inDB[i] = inDB[i] + outPiece
    dataList = [TmHisto, LenHisto, perGCHisto,cycleHisto]
    outList = [header[:]]
    for i in range(len(TmHisto)):
        piece = []
        for j in range(len(dataList)):
            piece.append(dataList[j][i])
        outList.append(piece)
    OFN = OFN + '_Histo.txt'
    OFN = LM.OutArray(outList,OFN,UserDir)
    #Compute avgs
    avgTm = sumTm/num_recs
    avgLen = sumLen/num_recs
    avgGC = sumGC/num_recs
    avgCycle = sumCycles/num_recs
    if Verbose:
        print '%i Probes analyzed AverageTm: %.1f AvgLength: %.1f AvgGC: %.1f AvgCycles: %i' %(num_recs,avgTm, avgLen,avgGC,avgCycle)
    if return_AnnotList:
        return inDB
    else:
        return OFN

def GenerateIncludedTargets(inDB, OFN, userDir, headerDict = None, Verbose = 1):
    """GenerateIncludedTargets
    |(inDB, OFN, userDir, headerDict = None,
    Verbose = 1)
    |Output: File|Exports a group of targets into a file of
    CombiProbe Designer format.
    |Input: (List of Sequences, Name of File to Generate, Path to create File
    in, Dictionary to define the file format.)
    |headDict={'index':i, 'seq_accession':None, 'seq_custom':1,'targetname':0,
    'probes':'1', 'replicates':'1', 'notes':2}"""
    
    i = 0
    headDict={'index':i, 'seq_accession':None, 'seq_custom':1, 'targetname':0,'probes':'1', 'replicates':'1', 'notes':2}
    if headerDict:
        headDict.update(headerDict)
    header = ['index', 'seq_accession', 'seq_custom', 'targetname', 'probes', 'replicates', 'notes']
    outDB = [header[:]]
    for obj in inDB:
        i += 1
        outObj = [i]
        for index in header[1:]:
            Index = headDict[index]
            if Index == None:
                outObj.append('')
            elif type(Index) == type(0):
                outObj.append(obj[Index])
            elif type(Index) == type(''):
                outObj.append(Index)
        outDB.append(outObj)
    OFN = OFN[:OFN.find('.')] + '_IncludedTargets.csv'
    CSVWriter(outDB, OFN, userDir)
    if Verbose == 1:
        print 'outDATA: Wrote ' + str(len(inDB)) + ' lines to file: ' + OFN
    return OFN

def GenerateIncludedProbes(inDB, OFN, userDir, headerDict = None, Verbose = 1):
    """GenerateIncludedProbes
    |(inDB, OFN, userDir, headerDict, Verbose)
    |Output: File
    |Exports a group of probes into a file of CombiProbe Designer format.
    |Input: (List of Sequences, File Name, Path to create File
    in, Dictionary to define the file format.)
    |headDict = {'index':i, 'captureprobe':0, 'name':1, 'include':'TRUE',
    'replicates':'1', 'sourcetarget':4, 'notes':2}"""
    
    #index   captureprobe   name   include   replicates   sourcetarget   notes
    i = 0
    headDict = {'index':None,'captureprobe':0, 'name':1, 'include':'TRUE', 'replicates':'1', 'sourcetarget':4,'notes':2}
    if headerDict:
        headDict.update(headerDict)
    header = ['index', 'captureprobe', 'name', 'include', 'replicates', 'sourcetarget', 'notes']
    outDB = [header[:]]
    for obj in inDB:
        i += 1
        outObj = [i]
        for index in header[1:]:
            Index = headDict[index]
            if Index == None:
                outObj.append('')
            elif type(Index) == type(0):
                outObj.append(obj[Index])
            elif type(Index) == type(''):
                outObj.append(Index)
        outDB.append(outObj)
    OFN = OFN[:OFN.find('.')] + '_IncludedProbes.csv'
    CSVWriter(outDB, OFN, userDir)
    if Verbose == 1:
        print 'outDATA: Wrote ' + str(len(outDB)) + ' lines to file: ' + OFN
    return OFN

def CSVWriter(inDB, OFN, userDir):
    """CSVWriter
    |(inDB, OFN, userDir)
    |Output: CSV File.
    |CSVWriter creates a CSV File from the output from either
    GenerateIncludedProbes or GenerateIncludedTargets.
    |Input: List of FASTA Objects, Name of Output file, Path.|"""
    outfile = open(userDir + OFN ,'w')
    for sLine in inDB:
        if type(sLine) == type([]):
            for element in sLine:
                outfile.write(str(element) + ',')
        else:
            outfile.write(str(sLine))
        outfile.write('\n')
    outfile.close()

def Reverse(stre):
    """Reverse
    |(stre)
    |Output: string
    |Reverses a string.||"""
    lcomp = list(stre)
    lcomp.reverse()
    return ''.join(lcomp)

def FindNextNonN(seq,pos):
        foundEnd  = None
        while pos < len(seq):
            if seq[pos] in 'GATCgatc':
                foundEnd = 1
                break
            pos +=1
        if not(foundEnd):
            return 'EOF'
        else:
            return pos

def FindNextN(seq,pos = 0):
        while pos < len(seq):
            if seq[pos] in 'Nn':
                break
            pos +=1
        return pos


     
class TilerObject:
    """TilerObject"""
    def __init__(self,FASTALst,INparams = None, Verbose = 0, FullAuto = None):
        self.FullAuto = FullAuto
        self.FASTALst = FASTALst
        self.n = len(FASTALst)
        self.TilerParams = {}
        self.WriteParameters(INparams)
        self.Verbose = Verbose
        self.OutTile = []
        self.ProbeObjectList = []
        self.LabelLine = ['label','sequence','length','SynthCycles',
                          'Tm','Start','Finish','LabelRoot',
                          'OriginalFASTALabel','ProbeGood',
                          'HP','PD','repeatNum','GC','HPstruct',
                          'PDstruct','failreason', 'score', 'duplicate']
        self.MarkedDuplicates = None
        self.NumProbeObjects = 0
        self.NumOutTile = 0
        self.GenesDesigned = 0
        self.NumDuplicateProbes = None
    def ProcessFASTALabelList(self):
        print 'ProcessFASTALabelList'
        for i in range(len(self.FASTALst)):
            label = self.FASTALst[i][0]
            newLabel = self.processFASTALabel(self, label) + '|'
            self.FASTALst[i][0] = newLabel
    def processFASTALabel(self, label):
        ParseLab = self.TilerParams['parseLabelType']
        return ParseFASTAHeader(label, ParseLab)
    
    def WriteParameters(self,InParams = None):
        TilerParams ={'type':1, 		##irrelevant
                  'strict':None,	##irrelevant
                  'parseLabelType':100,	##irrelevant Can be lambdaFunction
                  'ProbeType':'Hyb',	##irrelevant
                  'TmMin':68,		##<-- relevant
                  'MaxSynthCycles':150,    ##relevant for synthesis
                  'MinLen':20,	##<-- relevant
                  'MaxLen':35,	##<-- relevant, max is 40
                  'MaxGC':65,		##irrelevant
                  'MinGC':35,		##irrelevant
                  'half':1,	##Exclusive: Tiling parameter start half-way through previous probe
                  'overlap':None,	##Exclusive: Tiling parameter start x bases before end of previous probe
                  'startIncrement':None,	##Exclusive: Tiling parameter start 11 bases after start of previous probe
                  'MaxHPTm':55,		##irrelevant: Hairpin max Tm
                  'MaxPalindrome':70,	##irrelevant: Palindrom Max Tm 
                  'MaxRepeats':6,		##irrelevant: Max number of repeats 3.1 is TTT repeats 3.2 is TATATA
                  'CaptureEnd': None,		##See Note below
              'ProbeStartNumber' : 0,	##For labelling purposes only
                'TileFrom':None,       ##Default setting is None to tile whole sequence
                      'TileTo':None}    ##Default setting is None to tile whole sequence

        ##CaptureEnd means that if a new probe extends beyond the end of the sequence, then the probe bulds itself backwards for the last probe.
        TilerParams.update(InParams)
        if TilerParams['ProbeType'] == 'PCR':
            TilerParams.update({'TmSaltConc' : 0.05, 'TmMgConc': 1.0, #
                            'TmProbeConc':200.0e-9})
        elif TilerParams['ProbeType'] == 'Hyb':
            TilerParams.update({'TmSaltConc' : 0.33, 'TmMgConc': 0.0, #
                            'TmProbeConc':1.0e-9})
#         if not(InParams) and not(self.FullAuto):
#             TilerParams = US.AskParams(TilerParams, Title = '_Tiler Parameters_')
        else:
            TilerParams.update(InParams)
        self.TilerParams = TilerParams
    
    def TakeProbe(self,LabelRoot,LabelStart,LabelEnd,ProbeSeq, FASTALabel,ProbeStartNumber = 0):
        params = self.TilerParams
        MinLen = int(params['MinLen'])
        strict = params['strict']
        MaxLen = int(params['MaxLen'])
        strict = params['strict']
        TakeProbe = None
        probeLen = len(ProbeSeq)
        probeTm = MB.Tm(ProbeSeq, float(params['TmSaltConc']), int(params['TmMgConc']), float(params['TmProbeConc']))
        if probeLen >= MinLen:
            TakeProbe = 1
        elif not(strict)and probeLen > 15:
            TakeProbe = 1
        if TakeProbe:
            if LabelStart == None:
                newLabel = LabelRoot
                LabelStart,LabelEnd = 0,0
            else:
                newLabel = LabelRoot + 'sl-' + str(LabelStart+ProbeStartNumber) + '-' + str(LabelEnd+ProbeStartNumber)
            probeDescriptorLst = [newLabel,ProbeSeq, probeLen, probeTm, LabelStart, LabelEnd, LabelRoot, FASTALabel]
            self.OutTile.append(probeDescriptorLst)
            #self.ProbeDB(TileList = [probeDescriptorLst]) ##Moved to TBalanced: self.ProbeDB(self.OutTile)
        self.NumOutTile = len(self.OutTile)

    def INPUTProbes(self,FALst):
        params = self.TilerParams
        for obj in FALst:
            newLabel = obj[0]
            ProbeSeq = obj[1]
            probeLen = len(ProbeSeq)
            LabelStart,LabelEnd,LabelRoot,FASTALabel = None,None,newLabel,newLabel
            probeTm = MB.Tm(ProbeSeq, float(params['TmSaltConc']), int(params['TmMgConc']), float(params['TmProbeConc']))
            self.TakeProbe(LabelRoot,LabelStart,LabelEnd,ProbeSeq, FASTALabel)
            #self.OutTile.append([newLabel,ProbeSeq,probeLen,probeTm, LabelStart,LabelEnd,LabelRoot])
        self.ProbeDB()
    def TBalanced(self):
        print 'TBalanced:'
        params = self.TilerParams
        if self.Verbose:
            print 'Tiler is running with the following specs:'
            for y in params.iteritems():
                print y
        #Pull variables from Dictionary
        if params['ProbeStartNumber'] == 'FromFASTAObject':
            ProbeStartNumber = 'FromFASTAObject'
        else:
            ProbeStartNumber = int(params['ProbeStartNumber'])
        runType = int(params['type'])
        ParseLab = params['parseLabelType']
        sHalf = params['half']
        iStartIncr = params['startIncrement']##This takes precedence over Overlap, it is converted to int below
        iOverlap = params['overlap']##This only fires if there is no iStartIncr, it is converted to int below
        MinLen = int(params['MinLen'])
        strict = params['strict']
        MaxLen = int(params['MaxLen'])
        sProbe = ''
        NaConc = float(params['TmSaltConc'])
        MgConc = int(params['TmMgConc'])
        PrbConc = float(params['TmProbeConc'])
        TmMin = float(params['TmMin'])
        CaptureEnd = params['CaptureEnd'] ##if you cannot make a good probe near the end take it anyway?
        TileFrom = params['TileFrom']
        TileTo = params['TileTo']
        i = 0;nexti = i + 1000;seqsToTile = len(self.FASTALst)
        for fastaSeq in self.FASTALst:
            i +=1
            if self.Verbose:
                print 'Tiling : ' + fastaSeq[0]
            if ProbeStartNumber == 'FromFASTAObject':
                ProbeStartNumber = int(fastaSeq[3])  ##If it is different each time, then get here's a way to send it to a whole list of sequences
            FASTALabel = fastaSeq[0]
            LabelRoot = ParseFASTAHeader(FASTALabel, ParseLab)
            if LabelRoot[-1] <>'|':
                LabelRoot = LabelRoot + '|'
            sSeq = fastaSeq[1]
            ##Allows control on where tiling takes place It is a little primitive though
            ##All tiled sequences will be tiled the same way you cannot control where tiling takes place
            ##The only feature that alloows that is the ProbeStartNumber where you give the Input the number
            ##And you've done the slicing ahead of time.
            seqLen = len(sSeq)
            TileFrom, TileTo, SeqSpaceTiled = self.DetermineTilingRange(seqLen)
            
            if self.Verbose:
                print 'Tiling %s from %i to %i' %(LabelRoot,TileFrom, TileTo)
                #print 'Label:'+ LabelRoot + '  '+ 'seqlen '+ str(len(sSeq))
            
            start = TileFrom
            end = start + MinLen
            probeTm = 0.0
            probeLen = 0
            Looking = 1
            TakeProbe = None
            reINIT = None
            report = 0
            while Looking:
                ProbeSeq = sSeq[start:end]
                #print '%s\t%i\t%i' %(ProbeSeq, start, end)
                LabelStart = start
                LabelEnd = end
                probeTm = MB.Tm(ProbeSeq, NaConc, MgConc, PrbConc)
                probeLen = len(ProbeSeq)
                if probeTm == 0:#if 0 then non-acgt in probe
                    i = 0
                    pos = FindNextN(ProbeSeq,pos = 0)
                    end = start + pos
                    ProbeSeq = sSeq[start:end]
                    LabelEnd = end
                    ##self.TakeProbe(LabelRoot,LabelStart,LabelEnd,ProbeSeq, ProbeStartNumber)
                    start = FindNextNonN(sSeq,end +1)
                    if start <> 'EOF':
                        end = start + MinLen
                    else:
                        Looking = None
                    continue ##Essentially you already reinitialized the loop
                 ###This allows capture up to end of given sequence
                if end >= TileTo:
                    if CaptureEnd:
                        end = TileTo
                        start = end - MinLen
                        LabelStart = start
                        LabelEnd = end
                        probeTm = MB.Tm(ProbeSeq, NaConc, MgConc, PrbConc)
                        probeLen = len(ProbeSeq)
                        while probeTm < TmMin and probeLen <MaxLen and start >TileFrom:
                            start = start -1
                            LabelStart = start
                            ProbeSeq = sSeq[start:end]
                            probeLen = len(ProbeSeq)
                            probeTm = MB.Tm(ProbeSeq, NaConc, MgConc, PrbConc)
                            #print 'start %i end %i sProbe: %s Tm: %s' %(start, end, ProbeSeq, str(probeTm))
                    reINIT = None
                    Looking = None
                #print 'probeTm %i >= %i TmMin or probeLen %i >= %i MaxLen' %(int(probeTm), TmMin, probeLen, MaxLen)
                if probeTm >= TmMin or probeLen >= MaxLen or (end >= TileTo and probeLen >= MinLen):  
                    self.TakeProbe(LabelRoot,LabelStart,LabelEnd,ProbeSeq, FASTALabel, ProbeStartNumber)
                    reINIT = 1
                else:
                    end += 1
                if reINIT: #This is done more for the sake of logic than for neccessity
                    reINIT = None
                    #print 'end is at %i  StartIncr: %s' %(end, str(iStartIncr)),
                    if sHalf:
                        iOverlap = len(ProbeSeq)/2
                        start = start + iOverlap
                        end = start + MinLen
                    else:
                        if iStartIncr:
                            start = start + int(iStartIncr)
                        else:
                            start = end - int(iOverlap)   #new start
                        end = start + MinLen   #arbitrary restart
                        #print 'restarted at position %i iOverlap %s' %(end, str(iOverlap))
                    ProbeSeq = ''
                    probeLen = 0
            if i>nexti:
                nexti = i + 1000
                print 'Tbalanced processed %i of %i' %(i,seqsToTile)
        self.ProbeDB()
    def DetermineTilingRange(self, seqLen):
        ##Usage: TileFrom, TileTo, SeqSpaceTiled = self.DetermineTilingRange(seqLen)
        params = self.TilerParams
        TileFrom = params['TileFrom']
        TileTo = params['TileTo']
        if TileFrom == None:
            TileFrom = 0 ##> seqLen:
        if TileFrom < 0:
            TileFrom = seqLen + TileFrom
        if TileFrom < 0: TileFrom = 0
        if TileTo == None:
            TileTo =seqLen
        if TileTo < 0:
            TileTo = seqLen + TileTo
        if (TileTo > seqLen) or (TileTo <0):
            TileTo = seqLen
        SeqSpaceTiled = TileTo - TileFrom
        return TileFrom, TileTo, SeqSpaceTiled
    
    def ProbeDB(self,TileList = None, inparams = {}):
        params = self.TilerParams
        params.update(inparams)
        if not(TileList):##If nothing is brought in by TileList, then probeDB works on OutTile alone this is legacy
            TileList = self.OutTile
            self.ProbeObjectList = []
        
        for pobj in TileList:
            #Each pobj: [newLabel,ProbeSeq, probeLen, probeTm, LabelStart, LabelEnd, LabelRoot, FASTALabel]
            self.ProbeObjectList.append(PROBE.probeObject(pobj, params))
        if len(TileList)> 1:
            print 'Processed %i probes from %i FASTA sequences:' %(len(self.ProbeObjectList),len(self.FASTALst))
        self.NumProbeObjects = len(self.ProbeObjectList)
    def generateProbeDB(self):
        outDB = []
        for pobj in self.ProbeObjectList:
            outDB.append(pobj.ProbeDBOutput())
        return outDB
    def DeleteBadProbes(self):
        outDB = []
        prevSize = len(self.ProbeObjectList)
        self.OutTile = [] #Otherwise a call to probeDB() will restore all the bad probes
        for probeObj in self.ProbeObjectList:
            if probeObj.ProbeGood:
                outDB.append(probeObj)
        self.ProbeObjectList = outDB
        newSize = len(self.ProbeObjectList)
        print 'DeleteBadProbes: went from %i to %i probes in object' %(prevSize, newSize)
    def MarkRedundantProbes(self):
        self.MarkedDuplicates = 1
        self.NumDuplicateProbes = 0
        sortingFunction = lambda x: x.seq.upper()
        probeDB_BU = LM.LambdaListBreakup(self.ProbeObjectList, sortingFunction)
        for BU_obj in probeDB_BU:
            for probeObj in BU_obj:
                numduplicates = len(BU_obj)
                probeObj.duplicate = numduplicates
                if numduplicates > 1:
                    self.NumDuplicateProbes +=1
    def DeleteDuplicateProbes(self, KeepOneDuplicateCopy = 1):
        print 'DeleteDuplicateProbes:'
        if self.MarkedDuplicates == None:
            self.MarkRedundantProbes()
        outDB = []
        prevSize = len(self.ProbeObjectList)
        duplicateProbeObjLst = []
        duplicateSeqLst = []
        self.OutTile = [] #Otherwise a call to probeDB() will restore all the bad probes
        for probeObj in self.ProbeObjectList:
            if probeObj.duplicate < 2:
                outDB.append(probeObj)
            else:
                if not(probeObj.seq in duplicateSeqLst):
                    duplicateProbeObjLst.append(probeObj)
                    duplicateSeqLst.append(probeObj.seq)##Fixed this
        if KeepOneDuplicateCopy:
            outDB.extend(duplicateProbeObjLst)
            print 'DeleteDuplicateProbes: Keeping one copy of each duplicate:'
        self.ProbeObjectList = outDB
        newSize = len(self.ProbeObjectList)
        print 'DeleteBadProbes: went from %i to %i probes in object' %(prevSize, newSize)
    def RenderGoodSeqs(self):
        outDB = []
        for probeObj in self.ProbeObjectList:
            if probeObj.ProbeGood:
                outDB.append(probeObj.probeLst)
        return outDB
    def RenderGoodProbeObjects(self):
        outDB = []
        for probeObj in self.ProbeObjectList:
            if probeObj.ProbeGood:
                outDB.append(probeObj)
        return outDB
    def RenderAllProbeObjects(self):
        return self.ProbeObjectList
    def RenderAllSeqs(self,returnObject = None):
        outDB = [self.ProbeObjectList[0].ProbeDBOutput(LL=1)]
        for probeObj in self.ProbeObjectList:
            if returnObject:
                outDB.append(probeObj)
            else:
                outDB.append(probeObj.ProbeDBOutput())
        return outDB
    def RenderBestNSeqsPerIntervalLength(self,IntervalLength, scoreMatrix = None, returnObject = None, DefaultMinimun = 2):
        FASTADict = FA.BetterFASTALstObj(inList = self.FASTALst)
        outDB = []
        y = lambda x: x.OrigFASTALabel
        probeBU = LM.FastListBreakup(self.ProbeObjectList,y)
        if scoreMatrix:
            for probeobj in self.ProbeObjectList:
                probeobj.Score(scoreMatrix)
                
        bu_Lst = LM.FastListBreakup(self.ProbeObjectList, y)
        for fastProbeLst in bu_Lst:
            targetLabel = fastProbeLst[0].OrigFASTALabel
            targetLength = float(FASTADict[targetLabel].n)
            TileFrom, TileTo, SeqSpaceTiled = self.DetermineTilingRange(targetLength)
            probesToChoose = int(targetLength/IntervalLength) + DefaultMinimun ##at least two are always chosen
            
            fastProbeLst.sort(lambda x,y:cmp(x.score,y.score))
            choice = fastProbeLst[:probesToChoose]
            for probeObj in choice:
                if returnObject:
                    outlst = probeObj
                else:
                    outlst = probeObj.ProbeDBOutput()
                outDB.append(outlst)
        return outDB
    def RenderBestNSeqs(self,N, scoreMatrix = None, returnObject = None):
        outDB = []
        y = lambda x: x.OrigFASTALabel
        probeBU = LM.FastListBreakup(self.ProbeObjectList,y)
        if scoreMatrix:
            for probeobj in self.ProbeObjectList:
                probeobj.Score(scoreMatrix)
                
        bu_Lst = LM.FastListBreakup(self.ProbeObjectList, y)
        for fastProbeLst in bu_Lst:
            fastProbeLst.sort(lambda x,y:cmp(x.score,y.score))
            choice = fastProbeLst[:N]
            for probeObj in choice:
                if returnObject:
                    outlst = probeObj
                else:
                    outlst = probeObj.ProbeDBOutput()
                outDB.append(outlst)
        return outDB
    def RenderBestPercentSeqs(self,perc, scoreMatrix = None, Min = 1):
        outDB = []
        y = lambda x: x.LabelRoot
        if scoreMatrix:
            for probeobj in self.ProbeObjectList:
                probeobj.Score(scoreMatrix)
        bu_Lst = LM.FastListBreakup(self.ProbeObjectList, y)
        for fastProbeLst in bu_Lst:
            N = int(perc * len(fastProbeLst)) + Min
            fastProbeLst.sort(lambda x,y:cmp(x.score,y.score))
            choice = fastProbeLst[:N]
            for probeObj in choice:
                outlst = probeObj.ProbeDBOutput()
                outDB.append(outlst)
        #print "RenderBestPercentSeqs: N",N, 'len buLst', len(bu_Lst), outDB[0] 
        return outDB
    def RenderBestPercentSeqs_myfilter(self,perc, scoreMatrix = None, Min = 1):
        #filter out seqs that are <100bp
        outDB = []
        y = lambda x: x.LabelRoot
        if scoreMatrix:
            for probeobj in self.ProbeObjectList:
                probeobj.Score(scoreMatrix)
        bu_Lst = LM.FastListBreakup(self.ProbeObjectList, y)
        for fastProbeLst in bu_Lst:
            N = int(perc * len(fastProbeLst)) + Min
            fastProbeLst.sort(lambda x,y:cmp(x.score,y.score))
            choice = fastProbeLst[:N]
            for probeObj in choice:
                #myfilter
                if probeObj.len<100:continue
                outlst = probeObj.ProbeDBOutput()
                outDB.append(outlst)
        print "RenderBestPercentSeqs: N",N, 'len buLst', len(bu_Lst), outDB[0] 
        return outDB
    def BlindTile(self):
        params = self.TilerParams
        startIncr = params['startIncrement']
        overlap = params['overlap']##This only fires if there is no iStartIncr, it is converted to int below
        TargLen = int(params['MinLen'])
        strict = params['strict']
        parseLabelType = params['parseLabelType']
        avoidN = 1
        Strict = 1
       
        #BlindTile(FASTAList, startIncr, TargLen, parseLabelType, LabelStartNumber = 0, overlap = None, AvoidN = 1, strict = 1, Verbose = 0)
        self.OutTile = BlindTile(self.FASTALst, startIncr, TargLen, parseLabelType, overlap = overlap, AvoidN = avoidN, strict = Strict, Verbose = self.Verbose)
        #print 'BlindTile: TargLen = int(params[MinLen])', int(params['MinLen']), self.OutTile[0]
        for obj in self.OutTile:
            self.ProbeObjectList.append(PROBE.probeObject(obj))
    def TallyGeneProbes(self, Rigorous = 1):
        if self.NumProbeObjects <> self.NumOutTile:
            self.ProbeDB()##This means that the probes have not been properly porcessed yet
        if Rigorous:
            y = lambda x: x.LabelRoot
            bu_Lst = LM.FastListBreakup(self.ProbeObjectList, y)
            designedProbesByLabelRoot = len(bu_Lst)
        else:
            designedProbesByLabelRoot = 999999
        y = lambda x: x.OrigFASTALabel
        bu_Lst = LM.FastListBreakup(self.ProbeObjectList, y)
        designedProbesByFASTALabel = len(bu_Lst)
        print 'designedProbesByLabelRoot: %i  designedProbesByFASTALabel: %i' %(designedProbesByLabelRoot,designedProbesByFASTALabel)
        self.GenesDesigned = len(bu_Lst)
        ##Duplicates
        if self.NumDuplicateProbes ==None:
            self.MarkRedundantProbes()
        print 'In the dataset of %i probes, there are %i duplictes' %(self.NumProbeObjects, self.NumDuplicateProbes)
        ## Find genes that have no Probes designed at all and those that have no good probes designed
        designList = []
        goodDesignList = []
        for obj in bu_Lst:
            OrigFASTALabel = obj[0].OrigFASTALabel
            designList.append(OrigFASTALabel)
            numGoodProbes = 0
            for pobj in obj:
                if pobj.ProbeGood:
                    numGoodProbes +=1
                    break
            if numGoodProbes > 0: ##if you have limits get rid of the break !! "How many genes have > 5 good probes"
                goodDesignList.append(pobj.OrigFASTALabel)
        print 'Genes with probes designed %i genes with at least one good probe designed %i' %(len(designList), len(goodDesignList))
        seqsNotDesigned = []
        seqsNotGoodDesigned = []
        for faObj in self.FASTALst:
            label = self.processFASTALabel(faObj[0])
            if not(label in designList):
                seqsNotDesigned.append(faObj)
                continue
            if not(label in goodDesignList):
                seqsNotGoodDesigned.append(faObj)
        NumGenes = self.n
        print '%i genes Input %i genes Designed %i genes Not Designed %i genes Not Good ProbeDesigned' %(NumGenes, self.GenesDesigned, len(seqsNotDesigned), len(seqsNotGoodDesigned))
        
        return seqsNotDesigned, seqsNotGoodDesigned
    def RenderByInterval(self,probesToReturn,seqInterval):
        pass
    def ExportTAB(self,OFN,UD):
        outPATH = os.path.join(UD,OFN)
        outDB = [self.LabelLine] + self.generateProbeDB()
        outPATH.outDB(outDB)
        return outPATH
    def ExportFASTA(self,OFN,UD):
        outPATH = FA.path([UD,OFN],fileType = 'fasta')
        outPATH.outDB(self.OutTile)
    
        
            
        
        


def BlindTile(FASTAList, startIncr, TargLen, parseLabelType, LabelStartNumber = 0, overlap = None, AvoidN = 1, strict = 1, Verbose = 0):
    """BlindTile
    |(fastaList, startIncrement, length, parseLabelType, overlap = None)
    |Output: List of [FASTA labels, Seqs, Len(seq)]
    |Used to Blind Tile a given list of sequences.
    |Input: List of FASTA objects, Number of bases to skip per tile, size of
    sequence to create, Label type, number of bases to overlap per tile.|"""
    outTile = []
    for FASTASeq in FASTAList:
        label = ParseFASTAHeader(FASTASeq[0], parseLabelType)
        print label, FASTASeq[2]
        if label[-1] <>'|':
            label = label + '|'
        if Verbose:
            print label
        seq = FASTASeq[1]
        seqLength = len(seq)
        if overlap:
            Incr = TargLen - overlap
        else:
            Incr = startIncr
        for start in range(0, seqLength, Incr):
            end = start + TargLen
            #if end > seqLength : break
            probe = seq[start:end]
            probeL = len(probe)
            Lstart = start + LabelStartNumber
            newLabel = label + str(Lstart) + '-' + str(Lstart+probeL)+'|'
            if strict and len(probe)!=TargLen:continue
            outTile.append([newLabel, probe, len(probe)])
    return outTile

def SNPBlindTile(fastaList, startIncr, TargLen, sParseLab, overlap = None, HybSeq = None, ProbeAdjNumber = 0, BaseList = ['a', 'c', 'g', 't']):
    """SNPBlindTile
    |(fastaList, startIncr, TargLen, sParseLab, overlap = None)
    |Output: List of [FASTA labels, Seqs, Len(seq)]
    |Used to Blind Tile a given list of sequences and create four sequences
    per tile: one original and 3 mutants with a SNP inserted into the
    middle of the sequence.
    |Input: List of FASTA objects, Number of bases to skip
    per tile, size of sequence to create, Label type, number of bases to
    overlap per tile.|"""
    outTile = []
    p = ProbeAdjNumber
    if overlap:
        Incr = TargLen - overlap
    else:
        Incr = startIncr
        
    for fastaSeq in fastaList:
        label = ParseFASTAHeader(fastaSeq[0], sParseLab)
        if label[-1] <>'|':
            label = label + '|'
        #print label
        seq = fastaSeq[1].upper()
        seqLength = int(fastaSeq[2])
        for start in range(0, seqLength, startIncr):
            end = start + TargLen
            if end > seqLength :
                break
            probe = seq[start:end]
            #Placement of point mutation (iSpot)
            if HybSeq:
                iSpot = 0
                mutType = 'H&S'
            else:
                probeLen = float(len(probe))
                if probeLen % 2 == 1.0:
                    iSpot = int((probeLen/2) - 0.5)
                else:
                    iSpot = int(probeLen/2)
                mutType = 'CMG'
            #Create Probes
            origBase = probe[iSpot].lower()
            for b in BaseList:
                if origBase == b:
                    base = b.upper()
                    mu = mutType + ':_'
                else:
                    base = b
                    mu = mutType + ':mu'
                probe = probe[:iSpot] + base + probe[iSpot+1:]
                newLabel = label + 'sl:%i:%i|Len:%i|%s-%i-%s' %(start+p, end, TargLen, mu, start+iSpot+p, base)
                outTile.append([newLabel,probe,len(probe)])
    return outTile


def FinddGMidPoint(seq):
    dgBalLst = [99,99,99]
    for i in range(3,len(seq)-2):
        left = seq[0:i]
        right = seq[i+1:]
        center = seq[i]
        ratio = abs(1-MB.dG(left)/MB.dG(right))
        dgBalLst.append(ratio)
        if DEBUG: print '%i\t%s [%s] %s %s' %(i, left, center, right, str(ratio))
    return LM.MinArr(dgBalLst, retIndex = 1)


def FindProbe(seq, Tm, TmFunction, start = 0):
    probeTm = 0.0
    end = start + 5
    SeqLen = len(seq)
    while probeTm < Tm:
        end +=1
        probeseq = seq[start:end]
        probeTm = TmFunction(probeseq)
        if end >= SeqLen: break
    return probeseq, start, end, probeTm


def CharacterizeProbe(probeObj, FilterParams = None, returnLabel = None):
    """CharacterizeProbe
    |(probeObj, FilterParams = None, returnLabel = None, Verbose = 1)
    |probe object needs to be of the form [label, seq, len, anything else, , , ]
    |
    |"""
    
    if returnLabel:
        return ['ProbeGood','HP','PD','repeatNum','GC','HPstruct','PDstruct','failreason']
    Limits = {'MaxHPTm':55,
            'MaxPalindrome':60,
            'MaxRepeats':6,
            'MaxGC':65.0,
            'MinGC':35.0}
    if FilterParams:
        Limits.update(FilterParams)
    
    
    MaxHPTm = float(Limits['MaxHPTm'])
    MaxPalindrome = float(Limits['MaxPalindrome'])
    MaxRepeats = int(Limits['MaxRepeats'])
    MaxGC = float(Limits['MaxGC'])
    MinGC = float(Limits['MinGC'])
    outLst = []

    failreason = 'good'
    ProbeGood = 1

    probeSeq = probeObj[1]
    for base in probeSeq:
        if not(base in 'GATCgatc'):
            ProbeGood = 0
            failreason = 'badBase'
            break
    if failreason == 'badBase':
        probeObj.extend([0,0,0,0,0,0,0,failreason])
        return probeObj
    
    probeSeqSSLst = MB.SecStruct(probeSeq, MinAnneal = 4)
    HPLst = []
    PDLst = []
    for obj in probeSeqSSLst:
        if obj[2]>0:
            HPLst.append(obj)
        else:
            PDLst.append(obj)
    PDobj = MB.StabSecStruct(PDLst,Pal = 1,OutType = 1)
    PD = PDobj[0]
    PDstruct = PDobj[1]
    HPobj = MB.StabSecStruct(HPLst, OutType = 1)
    HP = HPobj[0]
    HPstruct = HPobj[1]
    if PD > MaxPalindrome:
        ProbeGood = 0
        failreason = 'PD > ' + str(MaxPalindrome)
    if HP > MaxHPTm:
        ProbeGood = 0
        failreason = 'HP > ' + str(MaxHPTm)
    GC = MB.perGC(probeSeq)
    if GC > MaxGC:
        ProbeGood = 0
        failreason = '%GC > ' + str(MaxGC)
    if GC < MinGC:
        ProbeGood = 0
        failreason = '%GC < ' + str(MinGC)
    repeatNum = MB.repeats(probeSeq)
    
    repeatType = int((repeatNum - floor(repeatNum))*10)
    if  repeatNum > MaxRepeats:
        ProbeGood = 0
        failreason = 'repeats > ' + str(MaxRepeats) + '-' + str(repeatType)
    probeObj = probeObj[:7] ##prevents retarded Concatenations
    probeObj.extend([ProbeGood,HP,PD,repeatNum,GC,HPstruct,PDstruct,failreason])
    
    return probeObj
    


if __name__ == '__main__':
    
    Verbose = 1
    seq = 'gatcgatgctagctagtcgatcgtagctagtcgatcgatgctagctgatcgatcgtagctagctagtcg'
    tiledDBObj = TilerObject([['testSeq',seq,25]], PCRPrimerTilerparams, Verbose = 0)
    tiledDBObj.TBalanced()
    probeObjList = tiledDBObj.RenderBestNSeqs(10,scoreMatrix = DefaultPCRScoreMatrix, returnObject = 1)
    for obj in probeObjList:
        print '%s\t%s\t%s\t%s' %(obj.label, obj.seq, obj.score, obj.failreason)
    pobj = PROBE.probeObject(['testSeq','gatcgatgctagctgatcgatgctagtcgatgctagtcgatcg',25])
    pobj = PROBE.probeObject(pobj.retFASTA())
    pobj.displayProbe()
