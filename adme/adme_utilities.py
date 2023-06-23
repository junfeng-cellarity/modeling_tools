__author__ = 'jfeng1'
import tempfile,os
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oequacpac import *
from openeye.oeszybki import *

class MyTmpFile:
    def __init__(self,mySuffix):
        tuple = tempfile.mkstemp(suffix = mySuffix, text = False)
        if tuple is not None and len(tuple)==2:
            self.file_handle = tuple[0]
            self.full_pathname = tuple[1]
        else:
            self.file_handle = None
            self.full_pathname = None
        return

    def getFullPath(self):
        return self.full_pathname

    def getFileHandle(self):
        return self.file_handle

    def __del__(self):
        if self.full_pathname is not None:
            os.remove(self.full_pathname)

class MolUtilities():
    def __init__(self):
        return

    def convertToMol(self, molString, format):
        ifs = oemolistream()
        ifs.SetFormat(format)
        ifs.openstring(molString)
        mol = OEGraphMol()
        OEReadMolecule(ifs,mol)
        return mol

    def convertToInchiKey(self,mol):
        ofs = oemolostream()
        ofs.SetFormat(OEFormat_INCHIKEY)
        ofs.openstring()
        OEWriteMolecule(ofs,mol)
        return ofs.GetString()

    def molToMDL(self, mol):
        ofs = oemolostream()
        ofs.SetFormat(OEFormat_MDL)
        ofs.openstring()
        OEWriteMolecule(ofs,mol)
        return ofs.GetString()

    def molToSDF(self, mol):
        ofs = oemolostream()
        ofs.SetFormat(OEFormat_SDF)
        ofs.openstring()
        OEWriteMolecule(ofs, mol)
        return ofs.GetString()

    def get3DMol(self, mol):
        mol = OEMol(mol)
        if mol.GetDimension()!=3:
            omega = OEOmega()
            omega.SetStrictStereo(False)
            omega.SetMaxConfs(1)
            if omega(mol):
                return mol
            else:
                return None
        else:
            return mol

    def minimizeMol(self,mol):
        sz = OESzybki()
        option = OESzybkiOptions()
        option.SetOptimizerType(OEOptType_CG)
        option.SetRunType(OERunType_CartesiansOpt)
        option.SetOptCriteria(500,0.01)
        result = OESzybkiResults()
        if not sz(mol,result):
            return None
        else:
            return mol

    def prepareMoleculeForGlide(self,mol):
        if mol is not None:
            OETheFunctionFormerlyKnownAsStripSalts(mol)
            OEGetReasonableProtomer(mol)
            OEAddExplicitHydrogens(mol)
            mol = self.get3DMol(mol)
            self.minimizeMol(mol)
        return mol

    def splitMols(self, sdfString, nTrunks):
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(sdfString)
        mol = OEGraphMol()
        molTrunks = []
        outputList = []
        hasMol = {}
        for i in range(nTrunks):
            hasMol[i] = False
            ofs = oemolostream()
            ofs.SetFormat(OEFormat_SDF)
            ofs.openstring()
            outputList.append(ofs)

        i = 0
        while OEReadMolecule(ifs,mol):
            id = i%nTrunks
            OEWriteMolecule(outputList[id],mol)
            if not hasMol[id]:
                hasMol[id] = True
            i += 1
        ifs.close()

        for id,ofs in enumerate(outputList):
            if hasMol[id]:
                molTrunks.append(ofs.GetString())
            ofs.close()

        return molTrunks
