import tempfile,os
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oequacpac import *
from openeye.oeszybki import *
TAGS = ["mmff94smod_Sheff","mmff94smod_NoEstat"]
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

    def generate_relative_energies_oeomega(self, input_filename, output_filename):
        ifs = oemolistream()
        lowestEnergyDict = {}
        ifs.open(input_filename)
        mol = OEGraphMol()
        while (OEReadMolecule(ifs, mol)):
            title = mol.GetTitle()
            energy = float(OEGetSDData(mol, TAGS[1]))
            if title not in lowestEnergyDict:
                lowestEnergyDict[title] = energy
            else:
                if energy < lowestEnergyDict[title]:
                    lowestEnergyDict[title] = energy
        ifs.close()

        ifs.open(input_filename)
        ofs = oemolostream()
        ofs.open(output_filename)
        mol = OEGraphMol()
        while (OEReadMolecule(ifs, mol)):
            title = mol.GetTitle()
            energy = float(OEGetSDData(mol, TAGS[1]))
            lowest_energy = lowestEnergyDict[title]
            OESetSDData(mol, "deltaEnergy", "%5.2f" % (energy - lowest_energy))
            OEWriteMolecule(ofs, mol)
        ifs.close()
        ofs.close()
        return

    def convertToMolList(self,molString):
        molList = []
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(molString)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            molList.append(OEGraphMol(mol))
        return molList

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
            if mol.GetDimension()<3:
                OETheFunctionFormerlyKnownAsStripSalts(mol)
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


if __name__ == "__main__":
    molUtilities = MolUtilities()
    mol = molUtilities.convertToMol("CC1(N)CCN(CC2=C3N(C=C2)N=CN=C3NC2=CC(Br)=CC(O)=C2)CC1",OEFormat_SMI)
    mol = molUtilities.prepareMoleculeForGlide(mol)
    print (molUtilities.molToMDL(mol))