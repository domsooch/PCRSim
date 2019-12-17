import os, sys, inspect, random, glob
import pandas as pd

sys.path.append('../')
import Probe.primer3 as primer3
import utils.IO as IO


"""
This module was written March 2019. It is designed to take different probe output formats and combine thme into one class.
I experimented with the Map(dict) object. This has weird behavior when you reload the module from where it originates.
Both dot and [] access work well enough.

"""
#Global Vriables reset these for the project you want

group_count_dict = {'A2':563,
                    'A3':28,
                    'B':260,
                    'C1':712,
                    'C2':113,
                    'C3':39,
                    'D':366,
                    'E':621,
                    'F':299,
                    'G1':38,
                    'G2':239,
                    'H':441,
                    'I':242,
                    'K1':329,
                    'K2':132,
                    'K3':7,
                    'L':100,
                    'M':37,
                    'N':323,
                    'O':943,
                    'P':846,
                    'Q':183,
                    'S':2683,
                    'T':440,
                    'U':272,
                    'V':303,
                    'W':238,
                    'X':26,
                    'Y':97,}
groupLst = list(group_count_dict.keys())
groupLst.sort()


class Map(dict):
    """
    https://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary
    Example:
    m = Map({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
    """
    def __init__(self, inobj):
        #super(Map, self).__init__(*args, **kwargs)
        super(Map, self).__init__()
        #args = self.makeIterableDict(inobj)#This does not work from here
        for k, v in inobj:
            print('kv', k,v)
            self[k] = v
#         if kwargs:
#             for k, v in kwargs.iteritems():
#                 self[k] = v
    def set_groupLst(self, in_groupLst=[]):
        global groupLst
        if not(in_groupLst):
            self.groupLst = groupLst
        else:
            self.groupLst = in_groupLst
    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(Map, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(Map, self).__delitem__(key)
        del self.__dict__[key]
    def makeIterableDict(self, obj):
        if type(obj) is dict:
            #self.row_dict = obj
            return obj
        if  type(obj) is pd.core.series.Series:
            #self.row_pd = obj
            return obj
        if (type(obj) is list or type(obj) is tuple) and len(obj)==2:
            print('dict')
            keys, values = obj
            return dict(zip(keys, values))
        if inspect.isclass(obj.__class__):
            return obj.__dict__
        raise "Bad input for MakeIterable %s"%obj
    def process_specificity_profiles(self, pd_row):
        self.percHitLL=[]
        percHit=False
        self.MASimLL =[]
        MAsim=False
        d={}
        for k in list(pd_row.keys()):
            if '_perchit' in k:
                percHit=True
            elif 'MAsim' in k:
                MAsim=True
            else:
                d[k]=pd_row[k]
        if MAsim:
            self.MASimLL = ['%s_MAsim'%g for g in self.groupLst]
            self.group_MASimLst=[float(pd_row['%s_MAsim'%g]) for g in self.groupLst]
        if percHit:
            self.percHitLL = ['%s_perchit'%g for g in groupLst]
            self.group_percentHitLst=[float(pd_row['%s_perchit'%g]) for g in self.groupLst]
    def export_dict(self, LL=[]):
        if not(LL):
            LL=self.LL
        return dict( (k,self[k]) for k in LL)
    def export_pdrow(self, LL=[]):
        return pd.Series(self.export_dict(LL=LL), index=LL)

class Map_nd():
    """
    https://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary
    Example:
    m = Map({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
    """
    def __init__(self, inobj):
        for k, v in inobj:
            print('kv', k,v)
            self.__dict__[k] = v
    def set(self, k,v):
        self.__dict__[k] = v
    def __setattr__(self, key, value):
        self.__dict__[key] =  value
    def set_groupLst(self, in_groupLst=[]):
        global groupLst
        if not(in_groupLst):
            self.groupLst = groupLst
        else:
            self.groupLst = in_groupLst
    def makeIterableDict(self, obj):
        if type(obj) is dict:
            #self.row_dict = obj
            return obj
        if  type(obj) is pd.core.series.Series:
            #self.row_pd = obj
            return obj
        if (type(obj) is list or type(obj) is tuple) and len(obj)==2:
            print('dict')
            keys, values = obj
            return dict(zip(keys, values))
        if inspect.isclass(obj.__class__):
            return obj.__dict__
        raise "Bad input for MakeIterable %s"%obj
    def process_specificity_profiles(self, pd_row):
        self.percHitLL=[]
        percHit=False
        self.MASimLL =[]
        MAsim=False
        d={}
        for k in list(pd_row.keys()):
            if '_perchit' in k:
                percHit=True
            elif 'MAsim' in k:
                MAsim=True
            else:
                d[k]=pd_row[k]
        if MAsim:
            self.MASimLL = ['%s_MAsim'%g for g in self.groupLst]
            self.group_MASimLst=[float(pd_row['%s_MAsim'%g]) for g in self.groupLst]
        if percHit:
            self.percHitLL = ['%s_perchit'%g for g in groupLst]
            self.group_percentHitLst=[float(pd_row['%s_perchit'%g]) for g in self.groupLst]
    def export_dict(self, LL=[]):
        if not(LL):
            LL=self.LL
        return dict( (k,self[k]) for k in LL)
    def export_pdrow(self, LL=[]):
        return pd.Series(self.export_dict(LL=LL), index=LL)

class ProbeMap(Map_nd):
    def __init__(self, inobj, in_groupLst=[]):
        """
        Use this as a template
        Universal genotyper probe importer.
        possible inputs:
            pandas series row
            dict
            (keyLst, valLst)
        """
        self.set_groupLst(in_groupLst=in_groupLst)
        self.MASimLL=[]
        self.percHitLL=[]
        indict = self.makeIterableDict(inobj)
        self.process_specificity_profiles(indict)
        super(Map, self).__init__(indict)
        self.LL = ['probe_label','seq','group', 'num_hits','num_groups_hit','self_group_hits'] + self.percHitLL + self.MASimLL

def count_unique(inlst, in_set):
    for rec in inlst:in_set.add(rec)
    print('set: [n=%i] dedupe: %i'%(len(inlst), len(in_set)))
    return in_set

def collect_probes(glob_wildcardORLst=None, pathLst =[], out_path=''):
    """
    Imports either from probe list or from glob wildcard file or list
    :param glob_wildcardORLst:
    :param pathLst:
    :param out_path:
    :return:
    """
    if not(glob_wildcardORLst is None):
        if type(glob_wildcardORLst) is str:
            glob_wildcardORLst = [glob_wildcardORLst]
        for glob_wildcard in glob_wildcardORLst:
            print('glob_wildcard', glob_wildcard)
            pathLst.extend(glob.glob(glob_wildcard))
    print("collect_probes: [n=%i]"%(len(pathLst)))
    LL_Lst = []
    probeLst = []
    seqLst = []
    s=0
    for csv_path in pathLst:
        plst=IO.simple_inCSV(csv_path)
        LL=plst.pop(0)
        seq_idx = LL.index('seq')
        if LL_Lst:
            assert len(LL_Lst[-1]) == len(LL)
        LL_Lst.append(LL)
        probe_lst=[]
        for line in plst:
            if len(line)<len(LL):break
            seq = line[seq_idx]
            if not(seq in seqLst):
                probe_lst.append(line)
                seqLst.append(seq)
                s+=1
        source = "%s [n=%i]"%(os.path.basename(csv_path), len(probe_lst))
        for p in probe_lst:
            p.append(source)
            probeLst.append(p)
    LL.append('source')
    print('collect_probes: %i probes %i probes deduped'%(s, len(probeLst)))
    probe_df = pd.DataFrame.from_records(probeLst, columns=LL)
    probe_df.to_csv(out_path)
    return pd.read_csv(out_path, header=0)

def MakeProbeLst(probe_pathLst=[], probe_pd=None, Unique=True, ReturnDF=True, ProbeClass=ProbeMap, Verbosity=0.002):
    """
    Takes probelsit after probe design and imports it for probe choice
    Also tales in a dataframe as probe_pd
    :param probe_pathLst:
    :param Unique:
    :param ReturnDF:
    :param ProbeClass:
    :param Verbosity:
    :return:
    """
    if probe_pathLst and probe_pd is None:
        dfLst = [];in_set = set()
        for pth in probe_pathLst:
            df = pd.read_csv(pth, header=0)
            dfLst.append(df)
            unique_set = count_unique(df['seq'].values, in_set)
            print("MakeProbeLst %s  probes: %i unique: %i"%(os.path.basename(pth), len(df.index), len(unique_set)))
        probe_pd= pd.concat(dfLst, ignore_index=True)
        sz = len(probe_pd.index)
        print("MakeProbeLst: num_files: %i  total: %i"%(len(probe_pathLst), sz))
        if Unique:
            probe_pd = probe_pd.drop_duplicates(['seq'], keep='last')
            print("after dedupe on seq: %i"%len(probe_pd.index))
            sz = len(probe_pd.index)
        if ReturnDF: return probe_pd
    sz = len(probe_pd.index)
    pLst = []
    print("making probeDB n=%i"%sz)
    i = 0
    for row in probe_pd.iterrows():
        i+=1
        pLst.append(ProbeClass(row[1]))
        if random.random() < Verbosity:
            print('probedb [%i of %i]'%(i, sz), end='\r')
    return pLst

def ExportProbeObjLstTodf(probeDB, opath=None, LL=[]):
    if not(LL):LL = probeDB[0].LL
    probe_df = pd.DataFrame([x.export_dict(LL=LL) for x in probeDB], columns=LL)
    if opath:
        probe_df.to_csv(opath)
    return probe_df



#Custom Probe Genotyper Objects

class ProbeChoiceProbeMap(Map):
    """
    This class is used for chooseing probesets
    """
    def __init__(self, inobj, in_groupLst=[], group_hit_threshold = 0.95, cross_hyb_threshold = 0.05, KeepLL=True):
        """
        Use this as a template
        Universal genotyper probe importer.
        possible inputs:
            pandas series row
            dict
            (keyLst, valLst)
        """
        self.set_groupLst(in_groupLst=in_groupLst)
        self.MASimLL=[]
        self.percHitLL=[]
        self.bit_field_count=-1
        self.group_hit_threshold = group_hit_threshold
        self.cross_hyb_threshold = cross_hyb_threshold

        #Choice Matrics
        self.advantage_at_acceptance = 0
        self.entropy = 0
        self.advantage = 0
        self.advantage_crosshyb_sorter = 0
        self.cross_hyb = 0.0
        self.IsSpecific = False

        #Load Data
        indict = self.makeIterableDict(inobj)
        self.process_specificity_profiles(indict)
        self.n = len(self.group_percentHitLst)
        super(Map, self).__init__(indict)
        self.LL = ['probe_label','seq','group', 'num_hits','num_groups_hit','self_group_hits'] + self.percHitLL + self.MASimLL
        if KeepLL:
            for label in indict.keys():
                if not label in self.LL:
                    self.LL.append(label)
        self.process()
    def process(self):
        self.offTargetGroupHits()
        self.ComputeSelfGroupHitpercent()
        self.ComputeIsSpecific(self.cross_hyb_threshold)
    def Set(self, k,v):
        self[k] = v
    def offTargetGroupHits(self, max_nonhit_sim=0.05, max_offtarget_masim=0.10, hit_thresh_sim=0.95):
        self.offtargetgroups_hit=0;self.offtargetgenomes_hit = 0
        b=[];h=[]
        for g in range(len(self.groupLst)):
            group = self.groupLst[g]
            percHit = float(self['%s_perchit'%group])
            masim = float(self['%s_MAsim'%group])
            if percHit >= hit_thresh_sim:
                b.append('1')
                h.append(g)
            elif max_nonhit_sim < percHit < hit_thresh_sim:
                if percHit > max_offtarget_masim:
                    self.offtargetgroups_hit+=1
                    self.offtargetgenomes_hit+=group_count_dict[group]*percHit
                b.append('0')
            elif percHit <= max_nonhit_sim:
                b.append('0')
        self.bit_field = 'B'+''.join(b)
        self.hits = h
        self.LL.extend(['offtargetgroups_hit', 'offtargetgenomes_hit', 'bitfield', 'bit_field'])
    def ComputeSelfGroupHitpercent(self):
        """
        probe_label: cons_C3-K1|1747326-1747426
        This means that teh sequence came from consensus sequence cons_C3
        That was blasted against the K1 consensus sequence
        """
        group_group = self.probe_label.split('|')[0].replace('cons_','')
        self.group, self.anti_group = group_group.split('-')
        self.self_perchit = float(self['%s_perchit'%self.group])
        self.anti_perchit = float(self['%s_perchit'%self.anti_group])
        self.design_spec_score = (self.self_perchit-self.anti_perchit)/self.self_perchit if self.self_perchit>0 else 0.0
        self.LL.extend(['group', 'anti_group','self_perchit', 'anti_perchit', 'design_spec_score'])
    def ComputeSpecificityScore(self, cross_hyb_threshold):
        self.specificity_score = 0.0
        for s in self.group_percentHitLst:
            if s>self.group_hit_threshold:
                self.specificity_score += abs(1.0-s)
            if s<= cross_hyb_threshold:
                self.specificity_score +=s
    def ComputeIsSpecific(self, cross_hyb_threshold):
        """Determines Specificity v-a-v hits to groups that fall under the group_hit_threshold"""
        self.IsSpecific = True
        self.cross_hyb = 0.0
        for perc_hit in self.group_percentHitLst:
            if  perc_hit < self.group_hit_threshold:
                self.cross_hyb+= perc_hit
                if perc_hit > cross_hyb_threshold:
                    self.IsSpecific = False
        self.ComputeSpecificityScore(cross_hyb_threshold)
        return self.IsSpecific
    def set_bit_dedup_count(self, c):
        self.bit_field_count = c
        self.LL.append('bit_field_count')
    def emit_bitfield(self, lam=lambda x:x.better_bf):
        bfstr = lam(self)
        inlst= [int(b) for b in bfstr.split()[1:]]
        assert len(inlst)<64, "make_bitfield: Undefined behavior for larger bitfields"
        bf = 1<<63
        for i in range(len(inlst)):
            if inlst[i]==1: bf |= (1 << i)
        return bf
    def hits_groupToHitLst(self, groupToHit_Lst):
        num_hits = 0
        for h in self.hits:
            if groupToHit_Lst[h]:num_hits+=1
        return num_hits
    def set_groupToHitLst(self, groupToHit_Lst):
        for h in self.hits:
            if h<len(groupToHit_Lst):
                groupToHit_Lst[h] = 0
    def hit(self, idx):
        if idx in self.hits:
            return True
        return False
    def Compute_advantage_crosshyb_sorter(self):
        #Adds cross_hyb to computation of specificity
        #penalizes advantage score with crosshyb, effect is limited to 5
        self.advantage_crosshyb_sorter = self.advantage + self.cross_hyb/(float(self.n)*5)
    def DesignPCR(self, primerud, PCR_settings={}):
        inlst = [[self.probe_label, self.seq]]
        PCRSettingsDict = {'PRIMER_PRODUCT_SIZE_RANGE': '70-100', 'PRIMER_MIN_TM': '56',
                          'PRIMER_OPT_TM': '58', 'PRIMER_MAX_TM': '60'}
        PCRSettingsDict.update(PCR_settings)
        PrimerOUT_IFN = primer3.RunPrimer3(inlst, primerud,
                                           INSettingsDict=PCRSettingsDict,
                                           OFN_PATHorpath=None, root='test')
        pcr_dict = primer3.ParsePrimer3(PrimerOUT_IFN, primerud, InGetVals=[], returnDictLst=True)
        PCR_DESIGN_SUCCESS=True
        pcrLabels = ['PRIMER_LEFT_0_SEQUENCE', 'PRIMER_RIGHT_0_SEQUENCE', 'pdseq', 'SEQUENCE_TEMPLATE']
        if pcr_dict[0]['runIndex']=='fail':
            #print (self.probe_label, "FAIL PCR DEZ")
            for pl in pcrLabels:
                pcr_dict[0][pl] = 'fail'
            PCR_DESIGN_SUCCESS=False
        else:
            self.seq = pcr_dict[0]['pdseq']
        for pl in pcrLabels:
            self[pl] = pcr_dict[0][pl]
        self.LL.extend(pcrLabels)
        return PCR_DESIGN_SUCCESS
