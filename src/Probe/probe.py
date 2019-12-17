import sys, os, math, random
from importlib import reload

print ('Import probe.py')

parent_dir = os.path.dirname(__file__)
parent_dir_mod = parent_dir.replace(os.path.split(parent_dir)[1], '') # subtract the source file name
sys.path.append(parent_dir_mod)

import Probe.MolBio as MB#; reload(MB)
reload(MB)
Verbose = None


DefaultProbeFilterParameters = {'strict':None,	##irrelevant
                  'TmMin':68,		##<-- relevant
                  'MaxGC':65,		##irrelevant
                  'MinGC':35,		##irrelevant
                  'MaxHPTm':55,		##irrelevant: Hairpin max Tm
                  'MaxPalindrome':70,	##irrelevant: Palindrom Max Tm 
                  'MaxRepeats':6,		##irrelevant: Max number of repeats 3.1 is TTT repeats 3.2 is TATATA
                  'MaxSynthCycles':150,    ##relevant for synthesis
                  }


def makeRandomProbes(probesToMake, probeLength, FilterParams= {}):
    #Usage: outDB = makeRandomProbes(probesToMake, probeLength, FilterParams= {})
    outDB = []
    c = 0
    while len(outDB) < probesToMake:
        c +=1
        OLST = [random.choice(['a','t','c','g']) for b in range(probeLength)]
        seq = ''.join(OLST)
        faObj = ['randseq%i' %c, seq, len(seq)]
        pobj = probeObject(faObj, FiltParams = FilterParams)
        if pobj.ProbeGood:
            outDB.append(pobj.export())
    print ('makeRandomProbes: Made %i probes, chose %i probes' %(c, len(outDB)))
    outDB = [pobj.export(LL=1)] + outDB
    return outDB
    

class probeObject:
    """ProbeObject"""
    def __init__(self, probeLst, FiltParams = {}, TmRange = 15.0):
        self.TmRange = TmRange
        self.Params = DefaultProbeFilterParameters
        self.Params.update(FiltParams)
        ##Input for probeObject is: probeLst = [newLabel,ProbeSeq, probeLen, probeTm, LabelStart, LabelEnd, LabelRoot, FASTALabel]
        ##At a minimum it has to be [newLabel,ProbeSeq, probeLen]
        self.probeLst = probeLst
        self.SynthCycles = 0
        self.ProbeGood = 0
        self.HP = 0.0
        self.PD = 0.0
        self.len = len(probeLst[1])
        self.repeatNum = 0
        self.GC = 0.0
        self.HPstruct = ''
        self.PDstruct = ''
        self.selfAnnealScore = 0.0
        ##AllVariables
        self.InitializeProbeObject()
    def Antisense(self):
        self.label = self.label+'|as'
        self.seq = MB.Antisense(self.seq)
    def InitializeProbeObject(self):
        if Verbose: print ('InitializeProbeObject():')
        if 'instance' in str(self.probeLst).lower():
            self.InitializeProbeObjectFromProbeObject(self.probeLst)
        else:
            if len(self.probeLst) < 8:
                self.InitializeProbeObjectFromFAObj()
            (label,sequence,probeLen,probeTm, Start,End,LabelRoot,OrigFASTALabel) = self.probeLst[:8]
            self.label = label
            self.seq = sequence
            self.length = len(self.seq)
            self.Tm = float(probeTm)
            self.slStart = Start
            self.slEnd = End
            self.LabelRoot = LabelRoot
            self.OrigFASTALabel = OrigFASTALabel
            self.score = 0
            self.duplicate = 0 #0 means it has not been checked 1 means that it is unique 2is that there are two copies, etc
            self.failreason = ''
            self.ProbeGood = 0
            self.characterize()
    def InitializeProbeObjectFromFAObj(self):
        if Verbose: print (' InitializeProbeObjectFromFAObj():')
        self.probeLst = [self.probeLst[0],
                          self.probeLst[1],
                          self.probeLst[2],
                          MB.Tm(self.probeLst[1]),
                          0,
                          0,
                          'No_LabelRoot',
                          'No_OrigFASTALabel']
    def InitializeProbeObjectFromProbeObject(self, pObj):
        if Verbose: print ('InitializeProbeObjectFromProbeObject:')
        attributes = dir(pObj)
        for attribute in attributes:
            attrType=[]
            cmd= 'attrType = str(type(pObj.%s))'%attribute
            exec(cmd)
            if ('instancemethod' in attrType) or ('__' in attribute):
                continue
            else:
                ##print  'self.%s = pObj.%s' %(attribute, attribute)
                cmd='self.%s = pObj.%s' %(attribute, attribute)
                exec(cmd)
    def characterize(self, filtparams = {}):
        #[SynthCycles, ProbeGood,HP,PD,repeatNum,GC,HPstruct,PDstruct,failreason]
        ret = CharacterizeProbe(self.seq, filtparams = self.Params)
        #print 'ret: ' + str(ret)
        self.SynthCycles, self.ProbeGood, self.HP, self.PD, self.repeatNum, self.GC, self.HPstruct, self.PDstruct, self.failreason =ret[:9]

        ProbeFloor = self.Tm - self.TmRange
        self.BaseRepeat = self.repeatNum - int(self.repeatNum)
        
        self.HPScore = (self.HP - ProbeFloor)/self.TmRange
        if self.HPScore <0: self.HPScore = 0
        
        self.PDScore = (self.PD - ProbeFloor)/self.TmRange
        if self.PDScore <0: self.PDScore = 0
        
        self.GCScore = abs(self.GC - 50.0)/15.0
        if self.GCScore <1.0:
            self.GCScore = 0.0
        self.repScore = 0
        if self.BaseRepeat == 0.1 and self.repeatNum >5:
            self.repScore = self.repeatNum
        if self.BaseRepeat == 0.2 and self.repeatNum >6:
            self.repScore = self.repeatNum
        if not(self.repScore):
            if self.repeatNum >6:
                self.repScore = self.repeatNum
        self.score = self.Score()
    
    def Score(self,scoreMatrix = {'HP':1000,'PD':10,'GC':10000,'rep':10000,'Tm':100,'CyclePenalty':1000}):
        ###score is from best to worst highest score is worst probe
        score = 0
        if 'Tm' in scoreMatrix:
            score = score + abs(self.Params['TmMin'] - self.Tm)*scoreMatrix['Tm']
        if 'HP' in scoreMatrix:
            score = score + self.HPScore * scoreMatrix['HP']
        if 'PD' in scoreMatrix:
            score = score + self.PDScore * scoreMatrix['PD']
        if 'GC' in scoreMatrix:
            score = score + self.GCScore * scoreMatrix['GC']
        if 'rep' in scoreMatrix:
            score = score + self.repScore * scoreMatrix['rep']
        if 'PCRGC' in scoreMatrix:
            gcpcr = MB.perGC(self.seq[-6:])/100*6
            gcpcr = int(gcpcr)
            gcpcrScore = 0
            if gcpcr == 0:
                gcpcrScore = 3
            if gcpcr == 1 or gcpcr == 6:
                gcpcrScore = 2
            if gcpcr == 2 or gcpcr == 5:
                gcpcrScore = 1
            score = score + gcpcrScore * scoreMatrix['PCRGC']
            ##print self.seq + 'gcpcr is ' + str(gcpcr) + ' score = ' + str(score)
#         if 'self-Anneal' in scoreMatrix:
#             val = self.ComputeSelfAnneal()
#             val = max(0, val -30000)
#             score = score + val*scoreMatrix['self-Anneal']
        if 'CyclePenalty' in scoreMatrix:
            cycleIndex = max(0, (self.SynthCycles - self.Params['MaxSynthCycles']))
            score = score + cycleIndex * scoreMatrix['CyclePenalty']
        self.score = score
        return score
#     def ComputeSelfAnneal(self):
#         pass
#         val, q,s = PCR.ComputeSelfPrimerDimer(self.seq)
#         self.selfAnnealScore = val
#         return val
    def retFASTA(self):
        return [self.label,self.seq,self.length]
    def returnLL(self):
        return ['label','seq','length','SynthCycles','Tm',
                'slStart','slEnd','LabelRoot','OriginalFASTALabel',
                'ProbeGood', 'HP', 'PD',
                'repeatNum','GC','HPstruct',
                'PDstruct','failreason', 'score','duplicate']
    def ProbeDBOutput(self, LL = None):
        if LL:
            return self.returnLL()
        return [self.label,self.seq,self.length,self.SynthCycles, self.Tm,
                self.slStart,self.slEnd,self.LabelRoot,self.OrigFASTALabel,
                self.ProbeGood, self.HP, self.PD,
                self.repeatNum,self.GC,self.HPstruct,
                self.PDstruct,self.failreason, self.score,self.duplicate]
    def displayProbe(self):
        LL = self.returnLL()
        probeLst = self.ProbeDBOutput()
        for i in range(len(LL)):
            print ('%s\t\t%s' %(LL[i], str(probeLst[i])))
    def outDB(self, LL = None):
        if LL:
            return self.returnLL()
        else:
            return self.ProbeDBOutput()
    def export(self, LL = None):
        ##TRhis is so it conforms to nameing convention
        if LL:
            return self.returnLL()
        else:
            return self.ProbeDBOutput()

        
        
        


class GTProbeObject(probeObject):
    def __init__(self, probeObj = None, pdbLst = None, accPos = None):
        ## This is just to allocate the memory
        self.BaseRepeat = ''
        self.GC = ''
        self.GCScore = ''
        self.HP = ''
        self.HPScore = ''
        self.HPstruct = ''
        self.LabelRoot = ''
        self.OrigFASTALabel = ''
        self.PD = ''
        self.PDScore = ''
        self.PDstruct = ''
        self.Params = ''
        self.ProbeGood = ''
        self.SynthCycles = ''
        self.Tm = ''
        self.TmRange = ''
        self.duplicate = ''
        self.failreason = ''
        self.label = ''
        self.length = ''
        self.probeLst = ''
        self.repScore = ''
        self.repeatNum = ''
        self.score = ''
        self.seq = ''
        self.slEnd = ''
        self.slStart = ''
        
        if probeObj and ('probeObject' in str(probeObj)):
            self.importProbeObject(probeObj)
        else:
            probeObj = probeObject(pdbLst)
            self.importProbeObject(probeObj)
        
        if accPos:
            self.accNum = self.label.split('|')[accPos]
        else:
            self.accNum = self.LabelRoot
        ##GT probeObject Specific Variables
        self.subjHitLst = []
        self.percHitLst = None
        self.numHitLst = None
        self.ThermoProfile = None
        self.HitSum = 0
        self.ScoreSum = 0
        self.PerfectHits = None
        self.NumFamilyHits = 0
        self.dbSequence = ''##This is used in cases in which you want to walk alignemnts like when you do inf sequnec reconstructions
        self.bestHitTm = 0.0
    def importProbeObject(self, probeObj):
        self.InitializeProbeObjectFromProbeObject(probeObj)

    def MakeHitLists(self,BinObj):
        self.percHitLst = BinObj.returnBlankList()
        self.numHitLst = BinObj.returnBlankList()
        self.PerfectHits = BinObj.returnBlankList()
        self.ThermoProfile = [100 for obj in range(len(self.PerfectHits))]
        
    def scoreHit(self, HitObj, binObj, TmRange = 15):
        HitProbeTmRange = abs(self.Tm - HitObj.xTm)
        subjectLabel = HitObj.subject
        #binName = subjectLabel[:subjectLabel.find('|')]##3-08difference
        #binIndex, numMembers = binObj[binName]
        binIndex, numMembers = binObj.SubjectLabelToBinLabelConverter(subjectLabel)
        ##UpdateThermoProfile
        oldxTm = self.ThermoProfile[binIndex]
        if oldxTm > HitProbeTmRange:
            self.ThermoProfile[binIndex] = HitProbeTmRange
       ##Is Hit relevant ?
        if HitProbeTmRange < TmRange:
             ##Update BinUsage
            binObj.incrementBinUsageLst(binIndex)
             ##Update subjHitLst
            if not(subjectLabel in self.subjHitLst):
                self.subjHitLst.append(subjectLabel)
            ##Update others
            self.HitSum += 1
            self.ScoreSum += HitObj.score
            self.numHitLst[binIndex] +=1
            self.percHitLst[binIndex] = float(self.numHitLst[binIndex])/float(numMembers)
            self.NumFamilyHits = self.CountNonZeros(self.numHitLst)
            if HitProbeTmRange < 2:##This is for PCR
                self.countPerfectHit(binIndex)
##            if HitObj.xTm > self.bestHitTm:
##                self.bestHitTm = HitObj.xTm
##                self.dbSequence = HitObj.dbSequence #This only allows best dbSequence in
    def countPerfectHit(self,binIndex):
        self.PerfectHits[binIndex] +=1

    def ExportLst(self, binObj, AccPos = None, LL = None):
        if AccPos:
            self.accNum = self.label.split('|')[AccPos]
        else:
            self.accNum = self.LabelRoot
        if LL:
            #export labelLine
            binLst = binObj.binLst
            percHitbinList = binObj.generateBinLabelList(suffix='_perc')## [bin+'_perc' for bin in binLst]
            NumHitsbinList = binObj.generateBinLabelList(suffix='_#Hits')## [bin+'_#Hits' for bin in binLst]
            PerfectHitsbinList = binObj.generateBinLabelList(suffix='_perf')## [bin+'_perf' for bin in binLst]
            
            labelLine = ['index', 'label', 'seq','length','Start','Tm','accNum', 'chooseMe']
            labelLine = labelLine + percHitbinList + ['ScoreSum'] + NumHitsbinList + ['NumFamilyHits']
            return labelLine
        else:
            outLine = [0,self.label, self.seq, self.length,self.Start,self.Tm, self.accNum, 0]
            outLine = outLine + binObj.maskOutput(self.percHitLst) + [self.ScoreSum] + binObj.maskOutput(self.numHitLst) + [self.NumFamilyHits]
            return outLine
    def ThermoExportLst(self, binObj = None, AccPos = None, LL = None):
        if AccPos:
            self.accNum = self.label.split('|')[AccPos]
        else:
            self.accNum = self.LabelRoot
        if LL:
            #export labelLine
            binLst = binObj.binLst
            NumHitsbinList = binObj.generateBinLabelList(suffix='_#Hits') ## for bin in binLst]
            ThermoHitLst = binObj.generateBinLabelList(suffix='_dxTm') ##[bin+'_dxTm' for bin in binLst]
            binList = binObj.binLst
            labelLine = ['label', 'seq','length','Start','Tm','accNum']
            labelLine = labelLine + NumHitsbinList + ['ScoreSum'] + ThermoHitLst + ['NumFamilyHits']
            return labelLine
        else:
            outLine = [self.label, self.seq, self.length,self.slStart,self.Tm, self.accNum]
            outLine = outLine + binObj.maskOutput(self.numHitLst) + [self.ScoreSum] + binObj.maskOutput(self.ThermoProfile) + [self.NumFamilyHits]
            return outLine
    def ShortenedExportLst(self, binObj = None, AccPos = None, LL=None):
        if LL:
            labelLine = ['label', 'seq','length','Start','Tm','accNum']
            labelLine = labelLine + [ 'ScoreSum', 'NumFamilyHits'] 
            return labelLine
        if AccPos!=None:
            self.accNum = self.label.split('|')[AccPos]
        else:
            self.accNum = self.LabelRoot
        if binObj:
            #export labelLine
            binLst = binObj.binLst
            percHitbinList = []#[bin+'_perc' for bin in binLst]
            NumHitsbinList = []#[bin+'_#Hits' for bin in binLst]
            PerfectHitsbinList = []#[bin+'_perf' for bin in binLst]
            binList = binObj.binLst
            labelLine = ['label', 'seq','length','Start','Tm','accNum']
            labelLine = labelLine + percHitbinList + ['ScoreSum'] + NumHitsbinList + ['NumFamilyHits']
            return labelLine
        else:
            outLine = [self.label, self.seq, self.length,self.slStart,self.Tm, self.accNum]
            outLine = outLine + [self.ScoreSum] + [self.NumFamilyHits]
            return outLine
    def CountNonZeros(self, Lst):
        i = 0
        for j in Lst:
            if j != 0:
                i +=1
        return i











class geneObj:
    """geneObject"""
    def __init__(self,line,LineNum):
        line = line.replace('\n','')
        inLst = line.split('\t')
        self.inLst = inLst
        if len(inLst) == 29:
            if inLst[1] !='accession':
                self.LineNum = LineNum
                self.line = line + '\n'
                self.rec= inLst[0]
                self.AccNum= inLst[1]
                self.organismLN= inLst[3]
                self.organismPhyl= inLst[4]
                self.defline= inLst[5]
                self.seqLen= inLst[6]
                if inLst[12] == '':
                    self.CDSgene=inLst[15].strip().replace(' ','_')
                else:
                    self.CDSgene= inLst[12].strip().replace(' ','_')
                self.CDSgeneID= inLst[13]
                self.CDSproduct= inLst[15]
                self.CDSstart= inLst[17]
                self.CDSend= inLst[18]
                self.CDSbasespan = inLst[20].replace('[','').replace(']','')
                self.CDSuniqueLabel= inLst[21]
                self.CDSlocus= inLst[22]
                self.CDSnote= inLst[23]
                self.CDSproteinSeq= inLst[25]
                self.sequence= inLst[28]
                self.badRec = None
            else:
                self.badRec = 1
        else:
            self.badRec = 1
    def fastaObjExport(self):
        outLst= [self.AccNum , self.CDSgene , self.CDSlocus , self.CDSbasespan,self.sequence,len(self.sequence),self.LineNum]
        return outLst
    



DefaultProbeFilterParameters = {'strict':None,	##irrelevant
                  'TmMin':68,		##<-- relevant
                  'MaxGC':65,		##irrelevant
                  'MinGC':35,		##irrelevant
                  'MaxHPTm':55,		##irrelevant: Hairpin max Tm
                  'MaxPalindrome':70,	##irrelevant: Palindrom Max Tm 
                  'MaxRepeats':6,		##irrelevant: Max number of repeats 3.1 is TTT repeats 3.2 is TATATA
                  'MaxSynthCycles':150,    ##relevant for synthesis
                  }

def CharacterizeProbe(probeSeq, filtparams = {}):
    """CharacterizeProbe
    |(probeObj, FilterParams = None, returnLabel = None, Verbose = 1)
    |probe object needs to be of the form [label, seq, len, anything else, , , ]
    |
    |"""
    Limits = DefaultProbeFilterParameters
    Limits.update(filtparams) ##this is if you have a special set of conditions

    #['ProbeGood','HP','PD','repeatNum','GC','HPstruct','PDstruct','failreason']
    MaxHPTm = float(Limits['MaxHPTm'])
    MaxPalindrome = float(Limits['MaxPalindrome'])
    MaxRepeats = int(Limits['MaxRepeats'])
    MaxGC = float(Limits['MaxGC'])
    MinGC = float(Limits['MinGC'])
    if 'MaxSynthCycles' in Limits:
        MaxSynthCycles = int(Limits['MaxSynthCycles'])
    else:
        MaxSynthCycles = 160

    failreason = 'good'
    ProbeGood = 1
   
    for base in probeSeq:
        if not(base in 'GATCgatc'):
            ProbeGood = 0
            failreason = 'badBase'
            break
    if failreason == 'badBase':
        return [0,0,0,0,0,0,0,0,0]
    
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
    
    repeatType = int((repeatNum - math.floor(repeatNum))*10)
    if  repeatNum > MaxRepeats:
        ProbeGood = 0
        failreason = 'repeats > ' + str(MaxRepeats) + '-' + str(repeatType)
    SynthCycles = CyclesToSynth(probeSeq)
    if SynthCycles > MaxSynthCycles:
        ProbeGood = 0
        failreason = 'SynthCycles %i>%i(max)' %(SynthCycles, MaxSynthCycles)
    
    probeObj = [SynthCycles, ProbeGood,HP,PD,repeatNum,GC,HPstruct,PDstruct,failreason]
    #print str(probeObj)
    return probeObj



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
            print ('Cycle %i base %s %i basemade %s' %(cycle, base, b, basemade))
    return cycle

if __name__ == '__main__':
    Verbose = 1
    pobj = probeObject(['testSeq','gatcgatgctagctgatcgatgctagtcgatgctagtcgatcg',25])
    pobj = probeObject(pobj.retFASTA())
    pobj.displayProbe()

    
    
