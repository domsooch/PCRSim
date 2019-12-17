import os, sys, copy
sys.path.append(os.path.realpath('.'))
sys.path.append(os.path.realpath('../'))

from utils.IO import *

def ParsePrimer3(PrimerOUT_IFN, UD, InGetVals = [], returnDictLst=False):
    OFN = os.path.join(UD, PrimerOUT_IFN)
    inBuff = open(PrimerOUT_IFN, 'r').read()
    RunLst = inBuff.split('\n=\n')[:-1]
    AllVals = ['SEQUENCE_ID ',
                'SEQUENCE_TEMPLATE',
                'PRIMER_TM_SANTALUCIA',
                'PRIMER_SALT_CORRECTIONS',
                'PRIMER_OPT_SIZE',
                'PRIMER_MIN_TM',
                'PRIMER_WARNING',
                'PRIMER_MAX_SIZE',
                'PRIMER_OPT_TM',
                'PRIMER_GC_CLAMP',
                'PRIMER_MIN_SIZE',
                'PRIMER_PRODUCT_TM_OLIGO_TM_DIFF',
                'PRIMER_PRODUCT_SIZE_RANGE',
                'PRIMER_MAX_TM',
                'PRIMER_LEFT_NUM_RETURNED',
                'PRIMER_RIGHT_NUM_RETURNED',
                'PRIMER_INTERNAL_NUM_RETURNED',
                'PRIMER_PAIR_NUM_RETURNED',
                'PRIMER_PAIR_0_PENALTY',
                'PRIMER_LEFT_0_PENALTY',
                'PRIMER_RIGHT_0_PENALTY',
                'PRIMER_LEFT_0_SEQUENCE',
                'PRIMER_RIGHT_0_SEQUENCE',
                'PRIMER_LEFT_0',
                'PRIMER_RIGHT_0',
                'PRIMER_LEFT_0_TM',
                'PRIMER_RIGHT_0_TM',
                'PRIMER_LEFT_0_GC_PERCENT',
                'PRIMER_RIGHT_0_GC_PERCENT',
                'PRIMER_LEFT_0_SELF_ANY_TH',
                'PRIMER_RIGHT_0_SELF_ANY_TH',
                'PRIMER_LEFT_0_SELF_END_TH',
                'PRIMER_RIGHT_0_SELF_END_TH',
                'PRIMER_LEFT_0_HAIRPIN_TH',
                'PRIMER_RIGHT_0_HAIRPIN_TH',
                'PRIMER_LEFT_0_END_STABILITY',
                'PRIMER_RIGHT_0_END_STABILITY',
                'PRIMER_PAIR_0_COMPL_ANY_TH',
                'PRIMER_PAIR_0_COMPL_END_TH',
                'PRIMER_PAIR_0_PRODUCT_SIZE',
                ]
    GetVals = ['SEQUENCE_ID',
                'SEQUENCE_TEMPLATE',
                'PRIMER_SALT_CORRECTIONS',
                'PRIMER_WARNING',
                'PRIMER_LEFT_0_SEQUENCE',
                'PRIMER_RIGHT_0_SEQUENCE',
                'PRIMER_LEFT_0',
                'PRIMER_RIGHT_0',
                'PRIMER_LEFT_0_TM',
                'PRIMER_RIGHT_0_TM',
                'PRIMER_LEFT_0_GC_PERCENT',
                'PRIMER_RIGHT_0_GC_PERCENT',
                'PRIMER_PAIR_0_PRODUCT_SIZE',
                ]
    for Retreive_Val in InGetVals:
        if Retreive_Val in AllVals:
            GetVals.append(Retreive_Val)
    LL = ['runIndex'] + copy.deepcopy(GetVals)+['pdLabel', 'pdseq', 'pcr_sz']
    outDB = []
    runIndex = 0
    dict_Lst = []
    for robj in RunLst:
        d = {}
        runObj = robj.split('\n')
        runDict = {}
        runIndex +=1
        FailedRun = None
        for outputline in runObj:
            lineSplit = outputline.split('=')
            if len(lineSplit) ==2:
                q, val  = lineSplit
                q = q.strip()
                runDict[q] = val.strip()
        outObj = [runIndex]
        for query in GetVals:
            if query in runDict:
                outObj.append(runDict[query])
            else:
                outObj.append('missing_%s'%query)
                FailedRun = 1
        if FailedRun ==None:
            seq = runDict['SEQUENCE_TEMPLATE']
            Fprimer = runDict['PRIMER_LEFT_0_SEQUENCE']
            f_loc = runDict['PRIMER_LEFT_0']
            r_loc = runDict['PRIMER_RIGHT_0']

            s_start = int(f_loc.split(',')[0])
            s_end = int(r_loc.split(',')[0])+1
            pdseq = seq[s_start:s_end]
            label = runDict['SEQUENCE_ID']
            if not(label[-1] =='|'):
                label = label + '|'
            pdLabel = label + '%spd-%i-%i'%(label, s_start, s_end)
            outObj = outObj + [pdLabel, pdseq, len(pdseq)]
        else:
            outObj.append('FAILEDRUN')
            outObj[0] = 'fail'
        if FailedRun and runDict:
            print ('This run Failed %s'%runDict['SEQUENCE_ID'])
            d['SEQUENCE_TEMPLATE']=runDict['SEQUENCE_TEMPLATE']
            d['PRIMER_WARNING']=runDict['PRIMER_WARNING']
            d['SEQUENCE_ID']=runDict['SEQUENCE_ID']
            d[LL[0]] = outObj[0]
        else:
            for idx in range(len(LL)):
                d[LL[idx]] = outObj[idx]
        outDB.append(outObj)
        dict_Lst.append(d)
    if returnDictLst:
        return dict_Lst
    outDB = [LL] + outDB
    WriteLstToFile(outDB, OFN, delim = '\t', buff_sz = 1000)
    return OFN



def RunPrimer3(inDB, UD, INSettingsDict = {}, OFN_PATHorpath = None, root = 'generic', Verbose=None):
    """RunPrimer3: Runs Primer 3
    |IFN_PATHorpath, OFN_PATHorpath
    |Output: OFN the output is primer3-type output further processing is neccessary
    |
    ||"""
    Primer3_IFN = os.path.join(UD, 'primer3_%s_infile.txt' %root)
    OUTpath = os.path.join(UD, 'primer3_%s_outfile.txt' %root)
    ofp = open(Primer3_IFN, 'w')
    SettingsDict = {
        'PRIMER_PRODUCT_SIZE_RANGE':'180-250',
        'PRIMER_OPT_SIZE':'19',
        'PRIMER_MIN_SIZE':'17',
        'PRIMER_MAX_SIZE':'25',
        'PRIMER_WARNING':'s (*)',
        'PRIMER_MIN_TM':'52',
        'PRIMER_OPT_TM':'55',
        'PRIMER_MAX_TM':'56',
        'PRIMER_TM_SANTALUCIA':'1',
        'PRIMER_SALT_CORRECTIONS':'1',
        'PRIMER_PRODUCT_TM_OLIGO_TM_DIFF':'1.0',
        'PRIMER_GC_CLAMP':'2',}
    SettingsDict.update(INSettingsDict)
    Settings = [key + '=' + SettingsDict[key] for key in SettingsDict.keys()]
    settingsStr = '\n'.join(Settings)
    for seqObj in inDB:
        if len(seqObj)<2:continue
        primerseqid = seqObj[0]
        seq = seqObj[1]
        outStr = 'SEQUENCE_ID = %s\nSEQUENCE_TEMPLATE=%s\n%s\n=\n' %(primerseqid, seq,settingsStr)
        outStr = outStr + ''
        ofp.write(outStr)
    ofp.close()
    cmd = 'primer3_core <%s >%s' %(Primer3_IFN , OUTpath)
    if Verbose:print (cmd)
    os.system(cmd)
    return OUTpath


if __name__ == '__main__':
    IFNorLst = [['seq_a', 'CAACGCCATGCATCGCCATCTTCATGCAGCGTACTGCGCCTTCGCCAGACGGAGCAACCATGTCTGCACCATCAGAGGTTGCGCCGTAGCCAACGATTTC', 100]]
    UD = './'
    PrimerOUT_IFN = RunPrimer3(IFNorLst, UD, INSettingsDict = {'PRIMER_PRODUCT_SIZE_RANGE':'70-100','PRIMER_MIN_TM':'56', 'PRIMER_OPT_TM':'58', 'PRIMER_MAX_TM':'60'}, OFN_PATHorpath = None, root = 'test')
    ofn_or_dict = ParsePrimer3(PrimerOUT_IFN, UD, InGetVals = [], returnDictLst=True)
    print (ofn_or_dict)

