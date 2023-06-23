#!/usr/bin/env python

from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.AtomPairs import Pairs,Torsions
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import EState
import sys,os,string,time
from random import randint
from stat import *
import cPickle as pickle
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
import math, time
from rdkit.Chem.SaltRemover import SaltRemover
from sets import Set


#########################
## Set the environment ##
#########################
R_path = "/pkg/R/2.15.3/EL6/x86_64/bin"
vsplus_local_path = "/home/ssciabol/bin/vsplus-1.1.1" ##fragollino.biogen.com
adme_model_path = "/UserUCDD/InSilico_ADME/Models" ##Global Path for ADME models
sdf_file = (sys.argv[1])[:-4]
activity_tag = sys.argv[2]
pH = 7.4
descType = sys.argv[4]
rfNtree = 100
nBatches = 24
nCpU = 20
runMode = sys.argv[3]
group_set = sys.argv[6]
#----------------------------------------------------------------------------------------------------------------------------------

remover = SaltRemover(defnData="[Cl,Br]")

###############
# Get Actvity #
###############
sdFile=Chem.SDMolSupplier("%s.sdf" % sdf_file)
act = {}
group = {} ##contains information about TRAIN and TEST sets
projects = {} ##contains project information
for mol in sdFile:
    if "NEW" in sys.argv[5]:
        if mol is None:
            print "mol not found"
        else:
            try:
                molName = mol.GetProp('BIO Number')
            except:
                molName = mol.GetProp('_Name')
            activity = "0.0000"
            act[molName] = "%.2f" % float(activity)
            pgroup = "TEST"
            group[molName] = pgroup
    else:
        if mol is None:
            print "mol not found"
        else:
            try:
                SmIlEs = Chem.MolToSmiles(mol)
            except:
                SmIlEs = "Smiles Not Found"
            try:
                pname = mol.GetProp('Concat Distinct;Project')
                tmp = string.strip(pname)
                tmp2 = tmp.replace("\n",";")
                pname = tmp2.replace(" ","_")
            except:
                pname = 'Project Not Found'
            try:
                molName = mol.GetProp('BIO Number')
            except:
                molName = mol.GetProp('_Name')
            try:
                activity = mol.GetProp('%s' % activity_tag)
            except KeyError:
                pass
            
            if "." in activity:
                ##regression, change the float format to two decimal digits
                act[molName] = "%.2f" % float(activity)
            else:
                ##classification, do not change the format
                act[molName] = activity 
            pgroup = mol.GetProp('%s' % group_set)
            group[molName] = pgroup
            projects[molName] = pname,SmIlEs


######################
#Get VS+ descriptors #
######################
vsfp = {}
if "VSPLUS" in descType:
    sdFile=Chem.SDMolSupplier("%s.sdf" % sdf_file)
    if runMode == "lsf":
        pass
    elif runMode == "pbs":
        pass
    elif runMode == "mc":
        import multiprocessing
        from multiprocessing import Pool         
        def chunks(l, n):
            """Yield successive n-sized chunks from l."""
            for i in xrange(0, len(l), n):
                yield l[i:i+n]
        
        def run_vsplus_cli(folderName):
            name = string.strip(folderName)
            parentDir = os.getcwd()
            dirName = "%s/%s" % (parentDir,name)
            os.chdir(dirName)
            os.system("%s/vsplus-cli vs_commandfile" % vsplus_local_path)
            os.chdir(parentDir)
        ##Prepare folders
        batches = list(chunks(range(len(sdFile)), int(math.ceil(len(sdFile)/float(nBatches)))))
        batchID = 1
        fIDs = open("folderIDs.lst", "w")
        for batch in batches:
            pwd = os.getcwd()
            newFolder = pwd+"/batch_%d" % batchID
            fIDs.write("batch_%d\n" % batchID)
            os.mkdir(newFolder)
            os.chdir(newFolder) 
            sdfout = Chem.SDWriter("batch_%d.sdf" % batchID)
            for molID in batch:
                mol = sdFile[molID]
                if mol is None:
                    continue
                else:
                    try:
                        molName = mol.GetProp('BIO Number')
                    except:
                        molName = mol.GetProp('_Name')                      
                    mol.SetProp("_Name", molName)
                    sdfout.write(mol)
            sdfout.close()       
            ##Write VSPLUS input file
            cfile = open("vs_commandfile", "w")
            cfile.write('set GRUB_PATH = "%s/grub.dat";\n' % vsplus_local_path)
            cfile.write('set PROTONATION = NORMALIZE %f;\n' % pH)
            #cfile.write('set PROTONATION = AS_IS;\n')
            cfile.write('import SDF "batch_%d.sdf" name by MOL_NAME;\n' % batchID)
            #cfile.write('save as "vsDescriptors.vs+";\n')
            cfile.write('export datamatrix as CSV "vsDescriptors_batch%d.csv" using separator ",";\n' % batchID)
            cfile.close()
            os.chdir(pwd)
            batchID+=1
        fIDs.close()
        ##Run vsplus using 
        batch_list = open("folderIDs.lst","r").readlines()
        p = Pool(processes=int(nCpU))
        p.map(run_vsplus_cli, batch_list)
        print "done running vsplus-cli on all the batches"
        p.close()    
        ##Getting descriptors
        os.system("cat batch_*/vsDescriptors_batch*csv > data0")
        os.system("grep -v Ob data0 > data")
        os.system("cat batch_*/vsDescriptors_batch*csv > header0")
        os.system("grep Ob header0 > header1")
        os.system("head -1 header1 > header")
        os.system("cat header > vsDescriptors.csv")
        os.system("cat data >> vsDescriptors.csv")
        vsIN = open("vsDescriptors.csv","r")
        vsheader = string.strip(vsIN.readline()).split(",")[1:]
        for line in vsIN.readlines():
            tokens = (string.strip(line)).split(",")
            vdes = tokens[1:]
            if "n/a" in vdes:
                pass
            else:
                vsfp[tokens[0]] = vdes
        os.system("rm -rf data0 data header0 header1 header batch_* completed n_comp folderIDs.lst vsDescriptors.csv")  
        vsIN.close()
    


#######################
#Get RDKit descriptors#
#######################
rdk = {}
maccs = {}
fcfp4_bit = {}
fcfp6_bit = {}
ecfp4_bit = {}
ecfp6_bit = {}
rdMD = {}
if "ECFP" in descType or "FCFP" in descType or "MACCS" in descType or "rdMolDes" in descType or "APFP" in descType or "RDK" in descType:
    #sdFile=Chem.SDMolSupplier("%s.sdf" % sdf_file)
    for Rmol in sdFile:
        if Rmol is None:
            print "mol not found"
        else:
            mol = remover.StripMol(Rmol)
            #print Chem.MolToSmiles(mol), Chem.MolToSmiles(Rmol)
            #props = list(mol.GetPropNames())
            try:
                molName = mol.GetProp('BIO Number')
            except:
                molName = mol.GetProp('_Name')
            ##rdMD RDKIT Descriptors
            MDlist = []
            MDlist.append(rdMolDescriptors.CalcTPSA(mol))
            MDlist.append(rdMolDescriptors.CalcFractionCSP3(mol))
            MDlist.append(rdMolDescriptors.CalcNumAliphaticCarbocycles(mol))
            MDlist.append(rdMolDescriptors.CalcNumAliphaticHeterocycles(mol))
            MDlist.append(rdMolDescriptors.CalcNumAliphaticRings(mol))
            MDlist.append(rdMolDescriptors.CalcNumAmideBonds(mol))
            MDlist.append(rdMolDescriptors.CalcNumAromaticCarbocycles(mol))
            MDlist.append(rdMolDescriptors.CalcNumAromaticHeterocycles(mol))
            MDlist.append(rdMolDescriptors.CalcNumAromaticRings(mol))
            MDlist.append(rdMolDescriptors.CalcNumHBA(mol))
            MDlist.append(rdMolDescriptors.CalcNumHBD(mol))
            MDlist.append(rdMolDescriptors.CalcNumLipinskiHBA(mol))
            MDlist.append(rdMolDescriptors.CalcNumLipinskiHBD(mol))            
            MDlist.append(rdMolDescriptors.CalcNumHeteroatoms(mol))
            MDlist.append(rdMolDescriptors.CalcNumRings(mol))
            MDlist.append(rdMolDescriptors.CalcNumRotatableBonds(mol))
            MDlist.append(rdMolDescriptors.CalcNumSaturatedCarbocycles(mol))
            MDlist.append(rdMolDescriptors.CalcNumSaturatedHeterocycles(mol))
            MDlist.append(rdMolDescriptors.CalcNumSaturatedRings(mol))
            MDlist.append(rdMolDescriptors.CalcHallKierAlpha(mol))
            MDlist.append(rdMolDescriptors.CalcKappa1(mol))
            MDlist.append(rdMolDescriptors.CalcKappa2(mol))
            MDlist.append(rdMolDescriptors.CalcKappa3(mol))
            MDlist.append(rdMolDescriptors.CalcChi0n(mol))
            MDlist.append(rdMolDescriptors.CalcChi0v(mol))
            MDlist.append(rdMolDescriptors.CalcChi1n(mol))
            MDlist.append(rdMolDescriptors.CalcChi1v(mol))
            MDlist.append(rdMolDescriptors.CalcChi2n(mol))
            MDlist.append(rdMolDescriptors.CalcChi2v(mol))
            MDlist.append(rdMolDescriptors.CalcChi3n(mol))
            MDlist.append(rdMolDescriptors.CalcChi3v(mol))
            MDlist.append(rdMolDescriptors.CalcChi4n(mol))
            MDlist.append(rdMolDescriptors.CalcChi4v(mol))
            MDlist.append(rdMolDescriptors.CalcAsphericity(mol))
            MDlist.append(rdMolDescriptors.CalcEccentricity(mol))   
            MDlist.append(rdMolDescriptors.CalcInertialShapeFactor(mol))
            MDlist.append(rdMolDescriptors.CalcExactMolWt(mol))  
            MDlist.append(rdMolDescriptors.CalcPBF(mol)) #Returns the PBF (plane of best fit) descriptor (http://dx.doi.org/10.1021/ci300293f) 
            MDlist.append(rdMolDescriptors.CalcPMI1(mol))
            MDlist.append(rdMolDescriptors.CalcPMI2(mol))
            MDlist.append(rdMolDescriptors.CalcPMI3(mol))
            MDlist.append(rdMolDescriptors.CalcRadiusOfGyration(mol))
            MDlist.append(rdMolDescriptors.CalcSpherocityIndex(mol))
            #MDlist.append(rdMolDescriptors.CalcNumBridgeheadAtoms(mol))
            #MDlist.append(rdMolDescriptors.CalcNumAtomStereoCenters(mol))
            #MDlist.append(rdMolDescriptors.CalcNumHeterocycles(mol))
            #MDlist.append(rdMolDescriptors.CalcNumSpiroAtoms(mol))
            #MDlist.append(rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol))
            MDlist.append(rdMolDescriptors.CalcLabuteASA(mol))
            MDlist.append(rdMolDescriptors.CalcNPR1(mol))
            MDlist.append(rdMolDescriptors.CalcNPR2(mol))
            #for d in rdMolDescriptors.CalcGETAWAY(mol): #197 descr (http://www.vcclab.org/lab/indexhlp/getades.html)
            #    MDlist.append(d)
            for d in rdMolDescriptors.PEOE_VSA_(mol): #14 descr
                MDlist.append(d)
            for d in rdMolDescriptors.SMR_VSA_(mol): #10 descr
                MDlist.append(d)
            for d in rdMolDescriptors.SlogP_VSA_(mol): #12 descr
                MDlist.append(d)
            for d in rdMolDescriptors.MQNs_(mol): #42 descr
                MDlist.append(d)
            for d in rdMolDescriptors.CalcCrippenDescriptors(mol): #2 descr
                MDlist.append(d)
            for d in rdMolDescriptors.CalcAUTOCORR2D(mol):  #192 descr
                MDlist.append(d)
            #for d in rdMolDescriptors.CalcAUTOCORR3D(mol):  #80 descr
            #    MDlist.append(d)
            #for d in rdMolDescriptors.CalcMORSE(mol):  #224 descr
            #    MDlist.append(d)
            #for d in rdMolDescriptors.CalcRDF(mol): #210 descr
                #MDlist.append(d)
            #for d in rdMolDescriptors.GetUSR(mol):   #12 descr
                #MDlist.append(d)
            #for d in rdMolDescriptors.GetUSRCAT(mol): #60 descr
                #MDlist.append(d)
            #for d in rdMolDescriptors.CalcWHIM(mol): #114 descr
                #MDlist.append(d)

            rdMD[molName] = MDlist
            
            ##Daylight like
            rdkfp = FingerprintMols.FingerprintMol(mol)
            rdk[molName] = rdkfp.ToBitString()        
            
            ##MACCS Keys
            maccs_fp = MACCSkeys.GenMACCSKeys(mol)
            maccs[molName] = maccs_fp.ToBitString()
            
            ##Morgan (Circular) Fingerprints (FCFP4) BitVector
            fcfp4_bit_fp = AllChem.GetMorganFingerprintAsBitVect(mol,2,useFeatures=True,nBits=1024)
            fcfp4_bit[molName] = fcfp4_bit_fp.ToBitString()
            
            ##Morgan (Circular) Fingerprints (FCFP6) BitVector
            fcfp6_bit_fp = AllChem.GetMorganFingerprintAsBitVect(mol,3,useFeatures=True,nBits=2048)
            fcfp6_bit[molName] = fcfp6_bit_fp.ToBitString()      
            
            ##Morgan (Circular) Fingerprints (ECFP4) BitVector
            ecfp4_bit_fp = AllChem.GetMorganFingerprintAsBitVect(mol,2,useFeatures=False,nBits=1024)
            ecfp4_bit[molName] = ecfp4_bit_fp.ToBitString()
            
            ##Morgan (Circular) Fingerprints (ECFP6) BitVector
            ecfp6_bit_fp = AllChem.GetMorganFingerprintAsBitVect(mol,3,useFeatures=False,nBits=2048)
            ecfp6_bit[molName] = ecfp6_bit_fp.ToBitString()
        
        

########################
# Get IMHB descriptors #
########################
imhb_db = {"O=[C&X3;!R]~[C&X3;!R]-[NX3&H1,nH1,NX3&H2]":"aC3a_1",
           "O=[C&X3;!R]~[#7;!R]-[NX3&H1,nH1,NX3&H2]":"aNa",
           "O=[C&X3;!R]~[C&X4;!R]-[NX3&H1,nH1,NX3&H2]":"aC4a_1",
           "O=[c,C&X3;R]~[c,C&X3;R]~[C&X3;!R]-[NX3&H1,nH1,NX3&H2]":"cC3aC3a_1",
           "O=[C&X3;!R]~[c,C&X3;R]~[c,C&X3;R]-!@[NX3&H1,nH1,NX3&H2]":"aC3cC3a_1",
           "O=[c,C&X3;R]~[#7;R]~[C&X3;!R]-[NX3&H1,nH1,NX3&H2]":"cNaC3a",
           "O=[C&X3;!R]~[#7;!R]~[c,C&X3;R]~[NX3&H1,nH1,NX3&H2]":"aNaC3c",
           "O=[C&X3;!R]~[#7;!R]~[C&X3;!R]-[NX3&H1,nH1,NX3&H2]":"aNaC3a_1",
           "O=[C&X3;!R]~[C&X4;!R]~[C&X3;!R]-[NX3&H1,nH1,NX3&H2]":"aC4aC3a_1",
           "O=[C&X3;!R]~[C&X4;!R]~[C&X4;!R]-[NX3&H1,nH1,NX3&H2]":"aC4aC4a_1",
           "O=[C&X3;!R]~[C&X4;!R]~[c,C&X3;R]~[c,C&X3;R]-[NX3&H1,nH1,NX3&H2]":"aC4aC3cC3a",
           "O=[C&X3;!R]~[#7;!R]~[c,C&X3;R]~[c,C&X3;R]-[NX3&H1,nH1,NX3&H2]":"aNaC3cC3a",
           "O=[C&X3;!R]~[#7;1R]~[C&X4;!R]~[C&X4;!R]-[NX3&H1,nH1,NX3&H2]":"aNaC4aC4a",
           "O=[C&X3;!R]~[#7;!R]~[C&X4;!R]~[C&X3;!R]-[NX3&H1,nH1,NX3&H2]":"aNaC4aC3a",
           "O=[C&X3;!R]~[#7;!R]~[C&X4;!R]~[C&X4;!R]~[C&X4;!R]-[NX3&H1,nH1,NX3&H2]":"aNaC4aC4aC4a",
           "O=[C&X3;!R]~[C&X4;!R]-[OH1]":"aC4a_2",
           "O=[C&X3;!R]~[c,C&X3;R]~[c,C&X3;R]-[OH1]":"aC3cC3a_2",
           "O=[C&X3;!R]~[C&X4;R]~[C&X4;R]-[OH1]":"aC4cC4a",
           "O=[C&X3;R]~[C&X4;R]~[C&X4;!R]-[OH1]":"cC4aC4a",
           "O=[C&X3;!R]~[C&X4;!R]~[C&X4;!R]-[OH1]":"aC4aC4a_2",
           "O=[C&X3;!R]~[#7;!R]~[C&X4;R]~[C&X4;R]-[OH1]":"aNaC4cC4a",
           "[nX2,NX2;R]~[c,C;R]~[C&X3;!R]~[NX3&H1,nH1,NX3&H2]":"aC3a_2",
           "[nX2,NX2;R]~[c,C;R]~[C&X4;!R]~[NX3&H1,nH1,NX3&H2]":"aC4a_3",
           "[nX2,NX2;R]~[c,C;R]~[#7;!R]~[C&X3;!R]~[NX3&H1,nH1,NX3&H2]":"aNaC3a_2",
           "[nX2,NX2;R]~[c,C;R]~[C&X3;!R]~[C&X3;!R]~[NX3&H1,nH1,NX3&H2]":"aC3aC3a",
           "[nX2,NX2;R]~[c,C;R]~[C&X4;!R]~[C&X4;!R]~[NX3&H1,nH1,NX3&H2]":"aC4aC4a_3",
           "[nX2,NX2;R]~[c,C;R]~[cX3;R]~[C&X3;!R]~[NX3&H1,nH1,NX3&H2]":"cC3aC3a_2",
           "[nX2,NX2;R]~[c,C;R]~[C&X4;!R]~[OH1]":"aC4a_4",
           "[nX2,NX2;R]~[c,C;R]~[C&X4;!R]~[C&X4;!R]~[OH1]":"aC4aC4a_4",
           "[O&D2;!R]~[C&X4;!R]~[C&X3;!R]~[NX3&H1,nH1,NX3&H2]":"aC4aC3a_2",
           "[O&D2;!R]~[c,C&X3;R]~[c,C&X3;R]~[NX3&H1,nH1,NX3&H2]":"aC3cC3a_3",
           "[O&D2;!R]~[C&X4;!R]~[C&X4;!R]~[NX3&H1,nH1,NX3&H2]":"aC4aC4a_5",
           "[O&D2;!R]~[c,C&X3;R]~[c,C&X3;R]~[C&X3;!R]~[NX3&H1,nH1,NX3&H2]":"aC3cC3aC3a",
           "[O&D2;!R]~[c,C&X3;R]~[c,C&X3;R]-[c,C&X3;R]~[NX3&H1,nH1,NX3&H2]":"aC3cC3aC3c",
           "[O&D2;!R]~[C&X3;!R]~[c,C&X3;R]~[c,C&X3;R]-!@[NX3&H1,nH1,NX3&H2]":"aC3aC3cC3a"}
imhb_vars = imhb_db.values()
imhb_fp = {}
if "IMHB" in descType:
    for mol in sdFile:
        if mol != None:
            #props = list(mol.GetPropNames())
            try:
                molName = mol.GetProp('BIO Number')
            except:
                molName = mol.GetProp('_Name')
            ##Fingerprinting and SSS
            hbfp = [0]*len(imhb_db.keys())
            for SmArTs in imhb_db.keys():
                query = Chem.MolFromSmarts('%s' % SmArTs)
                if mol.HasSubstructMatch(query):
                    sss_name = imhb_db[SmArTs]
                    pf_idx = imhb_vars.index(sss_name)
                    hbfp[pf_idx] = hbfp[pf_idx] + 1
            imhb_fp[molName] = hbfp
              
                

                  
            
                      
##########################
## Get MoKa descriptors ##
##########################
moka_fp = {}



####################
##Merge descriptors#
####################
dlist = descType.split("_")
combinedheader = []
dtable = {}
vsTest = 1
fcfp4Test = 1
ecfp4Test = 1
fcfp6Test = 1
ecfp6Test = 1
maccsTest = 1
rdkTest = 1
imhbTest = 1
rdMDTest = 1

################################################################
##Take the common set of keys among all the descriptors blocks #
################################################################
rdkSet = Set(rdk.keys())
maccsSet = Set(maccs.keys())
fcfp4Set = Set(fcfp4_bit.keys())
fcfp6Set = Set(fcfp6_bit.keys())
ecfp4Set = Set(ecfp4_bit.keys())
ecfp6Set = Set(ecfp6_bit.keys())
rdMDSet = Set(rdMD.keys())
vsSet = Set(vsfp.keys())
actSet = Set(act.keys())
if "VSPLUS" in dlist:
    commonKeys = vsSet & actSet
else:
    commonKeys = actSet
    
    
for key in commonKeys:
    name = key
    if act[key] != "":
        tmpTable = []
        activity = act[key]
        
        if "VSPLUS" in dlist:
            vD = vsfp[key]
            k = 0
            for i in vD:
                tmpTable.append(i)
                if vsTest:
                    combinedheader.append(vsheader[k])
                    k+=1
            vsTest = 0
        
                  
        if "FCFP4" in dlist:
            fcfp4D = fcfp4_bit[key]
            z = fcfp4D.replace('0','0,')
            o = z.replace('1','1,')
            f = o[:-1]   
            fcfp4D = f.split(",") 
            k = 1
            for i in fcfp4D:
                tmpTable.append(i)
                if fcfp4Test:
                    varname = "fcfp4_%d" % k
                    combinedheader.append(varname)
                    k+=1
            fcfp4Test = 0

        if "ECFP4" in dlist:
            ecfp4D = ecfp4_bit[key]
            z = ecfp4D.replace('0','0,')
            o = z.replace('1','1,')
            f = o[:-1]   
            ecfp4D = f.split(",") 
            k = 1
            for i in ecfp4D:
                tmpTable.append(i)
                if ecfp4Test:
                    varname = "ecfp4_%d" % k
                    combinedheader.append(varname)
                    k+=1
            ecfp4Test = 0

        if "FCFP6" in dlist:
            fcfp6D = fcfp6_bit[key]
            z = fcfp6D.replace('0','0,')
            o = z.replace('1','1,')
            f = o[:-1]   
            fcfp6D = f.split(",") 
            k = 1
            for i in fcfp6D:
                tmpTable.append(i)
                if fcfp6Test:
                    varname = "fcfp6_%d" % k
                    combinedheader.append(varname)
                    k+=1
            fcfp6Test = 0

        if "ECFP6" in dlist:
            ecfp6D = ecfp6_bit[key]
            z = ecfp6D.replace('0','0,')
            o = z.replace('1','1,')
            f = o[:-1]   
            ecfp6D = f.split(",") 
            k = 1
            for i in ecfp6D:
                tmpTable.append(i)
                if ecfp6Test:
                    varname = "ecfp6_%d" % k
                    combinedheader.append(varname)
                    k+=1
            ecfp6Test = 0

        if "MACCS" in dlist:
            maccsD = maccs[key]
            z = maccsD.replace('0','0,')
            o = z.replace('1','1,')
            f = o[:-1]   
            maccsD = f.split(",") 
            k = 1
            for i in maccsD:
                tmpTable.append(i)
                if maccsTest:
                    varname = "maccs_%d" % k
                    combinedheader.append(varname)
                    k+=1
            maccsTest = 0
            
        if "RDK" in dlist:
            rdkD = rdk[key]
            z = rdkD.replace('0','0,')
            o = z.replace('1','1,')
            f = o[:-1]   
            rdkD = f.split(",") 
            k = 1
            for i in rdkD:
                tmpTable.append(i)
                if rdkTest:
                    varname = "rdk_%d" % k
                    combinedheader.append(varname)
                    k+=1
            rdkTest = 0

        if "IMHB" in dlist:
            imhbD = imhb_fp[key]
            k = 1
            for i in imhbD:
                tmpTable.append(str(i))
                if imhbTest:
                    varname = "imhb_%d" % k
                    combinedheader.append(varname)
                    k+=1
            imhbTest = 0
            
        
        if "rdMolDes" in dlist:
            rdMD_des = rdMD[key]
            k = 1
            for i in rdMD_des:
                tmpTable.append(str(i))
                if rdMDTest:
                    varname = "rdMD_%d" % k
                    combinedheader.append(varname)
                    k+=1
            rdMDTest = 0
            
        tmpTable.append(activity)
        dtable[key] = tmpTable
combinedheader.append("activity")
#Save out temporary file
rawData = open("rawData.csv","w")
for h in combinedheader[:-1]:
    rawData.write("%s," % h)
rawData.write("%s\n" % combinedheader[-1])
for cmpd in dtable.keys():
    comboD = dtable[cmpd]
    rawData.write("%s," % cmpd)
    for d in comboD[:-1]:
        rawData.write("%s," % d)
    rawData.write("%s\n" % comboD[-1])
rawData.close()



##########################
## Build/Validate model ##
##########################
if "TEST" in sys.argv[5]:
    #--> Split in TRAIN and TEST
    sdFilteredData = open("rawData.csv","r")
    sdHead = sdFilteredData.readline()
    trainSet = open("trainSet.csv","w")
    testSet = open("testSet.csv","w")
    trainSet.write(sdHead)
    testSet.write(sdHead)
    for line in sdFilteredData.readlines():
        tokens = string.strip(line).split(",")
        if group[tokens[0]] == "TRAIN" or group[tokens[0]] == "train" or group[tokens[0]] == "Train" or group[tokens[0]] == "training":
            trainSet.write("%s," % tokens[0])
            for d1 in tokens[1:-1]:
                trainSet.write("%s," % d1)
            trainSet.write("%s\n" % tokens[-1])
        elif group[tokens[0]] == "TEST" or group[tokens[0]] == "test" or group[tokens[0]] == "Test":
            testSet.write("%s," % tokens[0])
            for d2 in tokens[1:-1]:
                testSet.write("%s," % d2)
            testSet.write("%s\n" % tokens[-1])
    trainSet.close()
    testSet.close()
    
    #--> Save Project Information File
    prjInfo = open("ProjectInfo.csv","w")
    for mname in projects.keys():
        prjInfo.write("%s,%s,%s\n" % (mname,projects[mname][1],projects[mname][0]))
    prjInfo.close()    

    if "RF" in sys.argv[5]:
        #--> 1. Build RF model using TRAIN data
        #--> 2. Predict PROPERTY for data in TEST
        rscript2 = open("rscript2.R","w")
        rscript2.write("library(randomForest)\n")
        rscript2.write("trainset = read.table('trainSet.csv',head=TRUE,sep=',')\n")
        rscript2.write("testset = read.table('testSet.csv',head=TRUE,sep=',')\n")
        rscript2.write("RFmodel = randomForest(activity ~ .,data=trainset,importance=TRUE,ntree=%d)\n" % rfNtree)
        rscript2.write("summary(RFmodel)\n")
        rscript2.write("RFmodel\n")
        rscript2.write("round(importance(RFmodel),3)\n")
        rscript2.write("RFpred = predict(RFmodel,testset[,1:(length(testset)-1)])\n")
        rscript2.write("cor(RFpred,testset[,length(testset)])\n") ##External testset correlation
        rscript2.write("cor(RFmodel$predicted,trainset[,length(trainset)])\n") ##OOB correlation
        #rscript2.write("write.table(RFpred,'RFpred.csv',quote=FALSE,sep=',')\n")
        rscript2.write("write.table(cbind(RFpred,testset[,length(testset)]),'RFpred.csv',quote=FALSE,sep=',')\n") #predicted_test,experimental_test
        #rscript2.write("write.table(RFmodel$predicted,'RFmodel_predicted_OOB.csv',quote=FALSE,sep=',')\n")
        rscript2.write("write.table(cbind(RFmodel$predicted,trainset[,length(trainset)]),'RFmodel_predicted_OOB.csv',quote=FALSE,sep=',')\n") #predicted_train_OOB,experimental_train
        rscript2.write("saveRDS(RFmodel, 'RFmodel.rds')")
        rscript2.close()
        os.system("%s/R CMD BATCH rscript2.R" % R_path)
        

elif "BUILD" in sys.argv[5]:
    #sdFilter()
    #os.system("cp sdFilteredData.csv trainSet.csv")
    os.system("cp rawData.csv trainSet.csv")
    if "RF" in sys.argv[5]:
        #--> 1. Build RF model using ALL data
        rscript2 = open("rscript2.R","w")
        rscript2.write("library(randomForest)\n")
        rscript2.write("trainset = read.table('trainSet.csv',head=TRUE,sep=',')\n")
        rscript2.write("RFmodel = randomForest(activity ~ .,data=trainset,importance=TRUE,ntree=%d)\n" % rfNtree)
        rscript2.write("summary(RFmodel)\n")
        rscript2.write("RFmodel\n")
        rscript2.write("round(importance(RFmodel),3)\n")
        rscript2.write("RFfit = predict(RFmodel,trainset[,1:(length(trainset)-1)])\n")
        rscript2.write("cor(RFfit,trainset[,length(trainset)])\n")
        rscript2.write("write.table(RFfit,'RFfit.csv',quote=FALSE,sep=',')\n")
        rscript2.write("write.table(RFmodel$predicted,'RFmodel_predicted_OOB.csv',quote=FALSE,sep=',')\n")
        rscript2.write("saveRDS(RFmodel, 'RFmodel.rds')")
        rscript2.close()
        os.system("%s/R CMD BATCH rscript2.R" % R_path)
        
    
elif "NEW" in sys.argv[5]:
    if "RF" in sys.argv[5]:
        print "Start running ADME predictions"
        #HLMmodel
        ##assume that a XXmodel.rds object has been created previously and stored in the working directory
        rscript2 = open("rscript2.R","w")
        rscript2.write("library(randomForest)\n")
        rscript2.write("predSet = read.table('rawData.csv',head=TRUE,sep=',')\n")
        rscript2.write("RFmodel = readRDS('%s/HLM_model.rds')\n" % adme_model_path)
        rscript2.write("RFpred = predict(RFmodel,predSet[,1:(length(predSet)-1)])\n")
        rscript2.write("write.table(RFpred,'HLM_pred.csv',quote=FALSE,sep=',')")
        rscript2.close()
        os.system("%s/R CMD BATCH rscript2.R" % R_path)
        os.system("rm -rf rscript1.Rout .RData rscript1.R trainSet.csv testSet.csv sdFilteredData.csv rscript2.R del rscript2.Rout RFmodel_predicted_OOB.csv RFfit.csv")        
        
        #hPPBmodel
        ##assume that a XXmodel.rds object has been created previously and stored in the working directory
        rscript2 = open("rscript2.R","w")
        rscript2.write("library(randomForest)\n")
        #rscript2.write("rawDataSet = read.table('rawData.csv',head=TRUE,sep=',')\n")
        rscript2.write("predSet = read.table('rawData.csv',head=TRUE,sep=',')\n")
        rscript2.write("RFmodel = readRDS('%s/hPPB_model.rds')\n" % adme_model_path)
        rscript2.write("RFpred = predict(RFmodel,predSet[,1:(length(predSet)-1)])\n")
        rscript2.write("write.table(RFpred,'hPPB_pred.csv',quote=FALSE,sep=',')")
        rscript2.close()
        os.system("%s/R CMD BATCH rscript2.R" % R_path)
        os.system("rm -rf rscript1.Rout .RData rscript1.R trainSet.csv testSet.csv sdFilteredData.csv rscript2.R del rscript2.Rout RFmodel_predicted_OOB.csv RFfit.csv")        
        
        #MDR1-MDCKmodel
        ##assume that a XXmodel.rds object has been created previously and stored in the working directory
        rscript2 = open("rscript2.R","w")
        rscript2.write("library(randomForest)\n")
        #rscript2.write("rawDataSet = read.table('rawData.csv',head=TRUE,sep=',')\n")
        rscript2.write("predSet = read.table('rawData.csv',head=TRUE,sep=',')\n")
        rscript2.write("RFmodel = readRDS('%s/MDR1-MDCK_model.rds')\n" % adme_model_path)
        rscript2.write("RFpred = predict(RFmodel,predSet[,1:(length(predSet)-1)])\n")
        rscript2.write("write.table(RFpred,'MDR1-MDCK_pred.csv',quote=FALSE,sep=',')")
        rscript2.close()
        os.system("%s/R CMD BATCH rscript2.R" % R_path)
        os.system("rm -rf rscript1.Rout .RData rscript1.R trainSet.csv testSet.csv sdFilteredData.csv rscript2.R del rscript2.Rout RFmodel_predicted_OOB.csv RFfit.csv")        
        
        #RLMmodel
        ##assume that a XXmodel.rds object has been created previously and stored in the working directory
        rscript2 = open("rscript2.R","w")
        rscript2.write("library(randomForest)\n")
        #rscript2.write("rawDataSet = read.table('rawData.csv',head=TRUE,sep=',')\n")
        rscript2.write("predSet = read.table('rawData.csv',head=TRUE,sep=',')\n")
        rscript2.write("RFmodel = readRDS('%s/RLM_model.rds')\n" % adme_model_path)
        rscript2.write("RFpred = predict(RFmodel,predSet[,1:(length(predSet)-1)])\n")
        rscript2.write("write.table(RFpred,'RLM_pred.csv',quote=FALSE,sep=',')")
        rscript2.close()
        os.system("%s/R CMD BATCH rscript2.R" % R_path)
        os.system("rm -rf rscript1.Rout .RData rscript1.R trainSet.csv testSet.csv sdFilteredData.csv rscript2.R del rscript2.Rout RFmodel_predicted_OOB.csv RFfit.csv")        
        
        #rPPBmodel
        ##assume that a XXmodel.rds object has been created previously and stored in the working directory
        rscript2 = open("rscript2.R","w")
        rscript2.write("library(randomForest)\n")
        #rscript2.write("rawDataSet = read.table('rawData.csv',head=TRUE,sep=',')\n")
        rscript2.write("predSet = read.table('rawData.csv',head=TRUE,sep=',')\n")
        rscript2.write("RFmodel = readRDS('%s/rPPB_model.rds')\n" % adme_model_path)
        rscript2.write("RFpred = predict(RFmodel,predSet[,1:(length(predSet)-1)])\n")
        rscript2.write("write.table(RFpred,'rPPB_pred.csv',quote=FALSE,sep=',')")
        rscript2.close()
        os.system("%s/R CMD BATCH rscript2.R" % R_path)
        os.system("rm -rf rscript1.Rout .RData rscript1.R trainSet.csv testSet.csv sdFilteredData.csv rscript2.R del rscript2.Rout RFmodel_predicted_OOB.csv RFfit.csv")        
        
        #Solubilitymodel
        ##assume that a XXmodel.rds object has been created previously and stored in the working directory
        rscript2 = open("rscript2.R","w")
        rscript2.write("library(randomForest)\n")
        #rscript2.write("rawDataSet = read.table('rawData.csv',head=TRUE,sep=',')\n")
        rscript2.write("predSet = read.table('rawData.csv',head=TRUE,sep=',')\n")
        rscript2.write("RFmodel = readRDS('%s/Solubility_model.rds')\n" % adme_model_path)
        rscript2.write("RFpred = predict(RFmodel,predSet[,1:(length(predSet)-1)])\n")
        rscript2.write("write.table(RFpred,'Solubility_pred.csv',quote=FALSE,sep=',')")
        rscript2.close()
        os.system("%s/R CMD BATCH rscript2.R" % R_path)
        os.system("rm -rf rscript1.Rout .RData rscript1.R trainSet.csv testSet.csv sdFilteredData.csv rscript2.R del rscript2.Rout RFmodel_predicted_OOB.csv RFfit.csv")        
        print "Finished running ADME predictions" 

        #################
        ## Cleaning Up ##
        #################
        os.system("rm -rf rscript1.Rout .RData rscript1.R rawData.csv trainSet.csv testSet.csv sdFilteredData.csv rscript2.R del rscript2.Rout RFmodel_predicted_OOB.csv RFfit.csv")
        
        ########################################
        # Merging ADME endpoints with SDF file #
        ########################################
        def load_preds(fileCSV):
            openFile = open(fileCSV,"r")
            openFile.readline()
            d_adme = {}
            for line in openFile.readlines():
                t = (string.strip(line)).split(",")
                d_adme[t[0]] = t[1] 
            return d_adme
        print "Start merging d360, insilico and adme data"
        Solubility_preds = load_preds("Solubility_pred.csv")
        rPPB_preds = load_preds("rPPB_pred.csv")
        RLM_preds = load_preds("RLM_pred.csv")
        MDR1_MDCK_preds = load_preds("MDR1-MDCK_pred.csv")
        hPPB_preds = load_preds("hPPB_pred.csv")
        HLM_preds = load_preds("HLM_pred.csv")
        
        sdFileOut = Chem.SDWriter("%s_cADME.sdf" % sdf_file)
        for mol in sdFile:
            if mol is None:
                print "mol not found"
            else:
                #props = list(mol.GetPropNames())
                try:
                    molName = mol.GetProp('BIO Number')
                except:
                    molName = mol.GetProp('_Name')
                try:
                    get_HLM = pow(10,float(HLM_preds[molName]))
                    get_hPPB = pow(10,float(hPPB_preds[molName]))
                    get_MDR1_MDCK = pow(10,float(MDR1_MDCK_preds[molName]))
                    get_RLM = pow(10,float(RLM_preds[molName]))
                    get_rPPB = pow(10,float(rPPB_preds[molName]))
                    get_Solubility = float(Solubility_preds[molName])
                    mol.SetProp("cHLM",str(get_HLM))
                    mol.SetProp("cHPPB",str(get_hPPB))
                    mol.SetProp("cMDR1",str(get_MDR1_MDCK))
                    mol.SetProp("cRLM",str(get_RLM))
                    mol.SetProp("cRPPB",str(get_rPPB))
                    mol.SetProp("cSolubility",str(get_Solubility))
                    sdFileOut.write(mol)
                except:
                    sdFileOut.write(mol)
        
        sdFileOut.close()   
        os.system("rm Solubility_pred.csv rPPB_pred.csv RLM_pred.csv MDR1-MDCK_pred.csv hPPB_pred.csv HLM_pred.csv")
        print "Finished merging d360, insilico and adme data"
