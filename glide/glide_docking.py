#!/usr/bin/env python
__author__ = 'jfeng1'

import csv
import os,glob,tempfile,shutil, subprocess
import json
import platform
import imp
imp.load_source("typing_extensions","/home/jfeng/schrodinger.ve/lib/python3.8/site-packages/pip/_vendor/typing_extensions.py")
from twisted.web import server
from twisted.web import xmlrpc as twisted_xmlrpc
from twisted.internet import threads
from glide_utilities import MolUtilities
from openeye.oechem import *
import schrodinger.forcefield.minimizer as minimizer
from schrodinger.forcefield.ffld_options import *
from schrodinger import structure
import base64
from parse_mopac import MOPAC
import molvs
from rdkit import Chem
from molvs import Standardizer
from molvs.tautomer import TautomerCanonicalizer
from molvs.fragment import LargestFragmentChooser
from molvs.charge import Uncharger
import re
import xmlrpc
from parse import parse

SCHRODINGER = os.getenv("SCHRODINGER")
OE_DIR = "/home/jfeng/openeye/bin/"
SUBSET_COMMAND = "%s/utilities/structsubset"%SCHRODINGER
CONVERT_COMMAND = "%s/utilities/sdconvert"%SCHRODINGER
OE_OMEGA_INPUT = "ligands.sdf"
OE_DOCKING_INPUT = "ligands_multiconfs.sdf"
OE_DOCKING_GRID = "receptor.oeb.gz"
OE_MAKE_RECEPTOR_COMMAND = "%s -protein %s -bound_ligand %s -receptor %s"
OE_HYBRID_COMMAND = "%s -mpi_np 40 -receptor %s -dbase %s -annotate_scores -dock_resolution High -hitlist_size 5000 -save_component_scores -num_poses %d -docked_molecule_file %s"
OE_POSIT_COMMAND = "%s -mpi_np 40 -receptor %s -dbase %s -hitlist_size 20000 -num_poses %d -allowed_clashes %s -minimum_probability 0.01 -sortby probability -docked_molecule_file %s"
#OE_OMEGA_COMMAND = "%s -mpi_np 40 -in %s -out %s -maxConfRange 200,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600 -rangeIncrement 1 -addtorlib %s"
#OE_OMEGA_COMMAND_PARAM = "%s -mpi_np 40 -in %s -out %s -param %s -addtorlib %s -ewindow %f -rms %f -maxconfs %d"
OE_OMEGA_COMMAND = "%s -mpi_np 40 -in %s -out %s -sdEnergy -maxConfRange 200,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600 -rangeIncrement 1 -addtorlib %s -sampleHydrogens true -ewindow 5"
#OE_OMEGA_COMMAND_PARAM = "%s -mpi_np 40 -in %s -out %s -fromCT -maxConfRange 200,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600 -rangeIncrement 1 -param %s -addtorlib %s -sampleHydrogens true -ewindow %f -rms %f "
OE_OMEGA_COMMAND_PARAM = "%s -mpi_np 40 -in %s -out %s -sdEnergy -fromCT -maxConfRange 200,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600 -rangeIncrement 1 -sampleHydrogens true -ewindow %f -rms %f "

OE_OMEGA_COMMAND_MACROCYCLE = "%s macrocycle -mpi_np 16 -in %s -out %s -maxconfs 1000 -ewindow %f -rms %f "
OE_FREEFORM_COMMAND = "%s -in %s -track %s -prefix %s"
#OE_MINIMIZE_RECEPTOR_COMMAND = "-MMFF94S -optGeometry Honly %s %s"
#OE_MINIMIZE_COMPLEX_COMMAND = "-MMFF94S -p %s -in %s -out %s -protein_elec PB -optGeometry cart"
OE_MINIMIZE_LIGAND_SIMPLE = "%s -MMFF94S %s %s"
OE_MINIMIZE_LIGAND_FIX = "%s -MMFF94S -fix_file %s %s %s"
OE_MINIMIZE_LIGAND_SHEFFIELD = "%s -MMFF94S -sheffield -am1bcc %s %s"
OE_MINIMIZE_LIGAND_PB = "%s -MMFF94S -solventPB -solventCA 0.005 %s %s"

CONF_GEN_PY = "/home/jfeng/bin/generate_confs.py %s %s %s %s" #/home/jfeng/bin/generate_confs.py input.sdf output.oeb ewindow rmsd
OE_ROCS = "%s/rocs -mpi_np 24 -bestHits 0 -query %s -dbase %s -hitsfile %s -outputquery false"
OE_ROCS_EON_INPUT = "%s/rocs -mpi_np 24 -bestHits 0 -outputquery false -query %s -dbase %s -eon_input true -eon_maxconfs 50 -eon_input_size 0 -eon_input_file %s"
OE_EON = "%s/eon -mpi_np 24 -fixpka_query false -fixpka_dbase false -bestHits 0 -outputquery false -query %s -dbase %s -hitsfile %s"

DOCK_GRID_DIR = "/home/jfeng/Models/DockingModels"
TEMPLATE_DIR = "/home/jfeng/Models/LigandModels"
OE_PARAM_DIR = "/home/jfeng/Models/OEModels"
PYMOL_EXE = "/home/jfeng/Programming/pymol/bin/pymol"

MACROMODEL = "/opt/schrodinger/installations/default/macromodel -WAIT %s"
DOCKING_METHOD = "confgen"
PRECISION = "SP"
SHAPE_SCREEN_FLEX = "/opt/schrodinger/installations/default/shape_screen -WAIT -flex -flexSearchMethod thorough -flexMaxConfs %d -flexMaxRelEnergy %f -shape %s -screen %s -JOB %s"
SHAPE_SCREEN_RIGID = "/opt/schrodinger/installations/default/shape_screen -WAIT -distinct -shape %s -screen %s -JOB %s"
NMR_CONVERT = "/home/jfeng/Programming/insilicotools/scripts/python/convert_to_lib.py"
MOKA_COMMAND = "/home/jfeng/MoKa/moka-4.0.9-linux/moka_cli"
ALOGD_CMD = "python /opt/schrodinger/installations/default/mmshare-v6.1/python/common/ld_protocols/ld_alogD.py -ph %f %s"
class MacroModelCmd:
    def __init__(self, keyword, arg1 = 0, arg2 = 0, arg3 = 0, arg4 = 0, arg5 = 0.0, arg6 = 0.0, arg7 = 0.0, arg8 = 0.0):
        self.keyword = keyword
        self.arg1 = arg1
        self.arg2 = arg2
        self.arg3 = arg3
        self.arg4 = arg4
        self.arg5 = arg5
        self.arg6 = arg6
        self.arg7 = arg7
        self.arg8 = arg8

    def get_cmd_line(self):
        return " %4s "%self.keyword+ " %6d"%self.arg1 + " %6d"%self.arg2 + " %6d"%self.arg3 + " %6d"%self.arg4 +" %10.04f"%self.arg5 + " %10.04f"%self.arg6 + " %10.04f"%self.arg7 + " %10.04f\n"%self.arg8

    @staticmethod
    def generate_default_com(com_fname, ewindow, rmsd):
        com_file = open(com_fname,"w")
        ewindow = ewindow*4.17828
        commands = []
        commands.append("input.mae\n")
        commands.append("output.mae\n")
        #commands.append(MacroModelCmd("NPRC",arg1=48,arg2=1000).get_cmd_line())
        commands.append(MacroModelCmd("DEBG",arg1=55).get_cmd_line())
        commands.append(MacroModelCmd("FFLD",arg1=16, arg2=1, arg5=1.0).get_cmd_line())
        commands.append(MacroModelCmd("SOLV",arg1=3,arg2=1).get_cmd_line())
        commands.append(MacroModelCmd("EXNB").get_cmd_line())
        commands.append(MacroModelCmd("BDCO",arg5=59.4427, arg6=99999.00000).get_cmd_line())
        commands.append(MacroModelCmd("READ").get_cmd_line())
        commands.append(MacroModelCmd("CRMS",arg6=rmsd, arg8=2.0).get_cmd_line())
        commands.append(MacroModelCmd("LMCS",arg1=1000,arg7=3.0,arg8=6.0).get_cmd_line())
        commands.append(MacroModelCmd("NANT").get_cmd_line())
        commands.append(MacroModelCmd("MCNV",arg1=1,arg2=5).get_cmd_line())
        commands.append(MacroModelCmd("MCSS",arg1=2,arg4=10.0).get_cmd_line())
        commands.append(MacroModelCmd("MCOP",arg1=1,arg5=0.5).get_cmd_line())
        commands.append(MacroModelCmd("DEMX",arg2=833,arg5=ewindow, arg6=2*ewindow).get_cmd_line())
        commands.append(MacroModelCmd("MSYM").get_cmd_line())
        commands.append(MacroModelCmd("AUOP",arg5=100.0).get_cmd_line())
        commands.append(MacroModelCmd("AUTO",arg2=1,arg3=1,arg4=1,arg6=1.0,arg8=8).get_cmd_line())
        commands.append(MacroModelCmd("CONV",arg1=2,arg5=0.01).get_cmd_line())
        commands.append(MacroModelCmd("MULT").get_cmd_line())
        commands.append(MacroModelCmd("MINI",arg1=1,arg3=2500).get_cmd_line())
        com_file.writelines(commands)
        return
    @staticmethod
    def generate_fast_com(com_fname, ewindow, rmsd):
        com_file = open(com_fname,"w")
        ewindow = ewindow*4.17828
        commands = []
        commands.append("input.mae\n")
        commands.append("output.mae\n")
        #commands.append(MacroModelCmd("NPRC",arg1=48,arg2=1000).get_cmd_line())
        commands.append(MacroModelCmd("DEBG",arg1=55).get_cmd_line())
        commands.append(MacroModelCmd("FFLD",arg1=16, arg2=1, arg5=1.0).get_cmd_line())
        commands.append(MacroModelCmd("SOLV",arg1=3,arg2=1).get_cmd_line())
        commands.append(MacroModelCmd("EXNB").get_cmd_line())
        commands.append(MacroModelCmd("BDCO",arg5=89.4427, arg6=99999.00000).get_cmd_line())
        commands.append(MacroModelCmd("READ").get_cmd_line())
        commands.append(MacroModelCmd("CRMS",arg6=rmsd, arg8=2.0).get_cmd_line())
        commands.append(MacroModelCmd("MCMM",arg1=1000).get_cmd_line())
        commands.append(MacroModelCmd("NANT").get_cmd_line())
        commands.append(MacroModelCmd("MCNV",arg1=1,arg2=5).get_cmd_line())
        commands.append(MacroModelCmd("MCSS",arg1=2,arg4=21.0).get_cmd_line())
        commands.append(MacroModelCmd("MCOP",arg1=1).get_cmd_line())
        commands.append(MacroModelCmd("DEMX",arg2=833,arg5=ewindow, arg6=2*ewindow).get_cmd_line())
        commands.append(MacroModelCmd("MSYM").get_cmd_line())
        commands.append(MacroModelCmd("AUOP",arg5=100.0).get_cmd_line())
        commands.append(MacroModelCmd("AUTO",arg2=1,arg3=1,arg4=1,arg6=1.0,arg8=8).get_cmd_line())
        commands.append(MacroModelCmd("CONV",arg1=2,arg5=0.05).get_cmd_line())
        commands.append(MacroModelCmd("MINI",arg1=1,arg3=500).get_cmd_line())
        com_file.writelines(commands)
        return

class pKa:
    def __init__(self, type,value,idx,std):
        self.type = type
        self.value = float(value)
        self.idx = int(idx)
        self.std = float(std)


    def __str__(self):
        dict = self.getDict()
        return json.dumps(dict)

    def getDict(self):
        dict = {}
        dict['type'] = self.type
        dict['value'] = self.value
        dict['atomid'] = self.idx
        dict['stddev'] = self.std
        return dict

class GlideDockingServer(twisted_xmlrpc.XMLRPC):
    def xmlrpc_generatePymolSession(self, jsonString):
        tmpdir = tempfile.mkdtemp(prefix="pymol_")
        f = open(os.path.join(tmpdir,"pymol.pml"),"w")

        dict = json.loads(jsonString)
        if "receptor" in dict:
            receptorStr = dict['receptor'].encode()
            if receptorStr is not None:
                f.write("receptor_str = %s\n"%receptorStr)
                f.write("cmd.read_pdbstr(receptor_str.decode(),'receptor')\n")

        if "reference" in dict:
            reference = dict['reference'].encode()
            if reference is not None:
                f.write("reference_str = %s\n"%reference)
                f.write("cmd.read_molstr(reference_str.decode(),'reference')\n")
                f.write("cmd.hide('lines','reference')\n")
                f.write("cmd.show('sticks','reference')\n")

        if "ligands" in dict:
            ligands = dict['ligands']
            for idx,ligandStr in enumerate(ligands):
                f.write("ligand_str = %s\n"%ligandStr.encode())
                f.write("cmd.read_molstr(ligand_str.decode(),'mol_%d')\n"%idx)
                f.write("cmd.hide('lines','mol_%d')\n"%(idx))
                f.write("cmd.show('sticks','mol_%d')\n"%(idx))
        f.write("cmd.set('valence',1)\n")
        f.write("cmd.hide('everything','h.')\n")
        f.write("cmd.save(os.path.join('%s','pymol.pse'))\n"%tmpdir)
        f.write("cmd.quit()\n")
        f.close()

        if os.path.exists(os.path.join(tmpdir,"pymol.pml")):
            pml = os.path.join(tmpdir,"pymol.pml")
            p = subprocess.Popen([PYMOL_EXE,"-c",pml])
            p.communicate()

        if os.path.exists(os.path.join(tmpdir,"pymol.pse")):
            pse = open(os.path.join(tmpdir,"pymol.pse"),"rb").read()
            return twisted_xmlrpc.Binary(pse)
        else:
            return None

    def xmlrpc_getAvailableDockingGrids(self):
        result = threads.deferToThread(self.getAvailableDockingGrids)
        return result

    def xmlrpc_getAvailableTemplates(self):
        result = threads.deferToThread(self.getAvailableTemplates)
        return result

    # ffld_version (integer module constant) - Use one of the valid force fields from mmcommon_get_valid_forcefields()
    # struct (schrodinger.structure.Structure) - The structure to be minimized.
    # cleanup (bool) - If True, attempts to automatically clean up the structure will be made. (This uses the C function mmlewis_apply().) Note that this can modify the atom types of the passed in structure.
    # min_method (enum) - The minimizer method. Valid values are MinBFGS, MinAUTO and MinCG. Default is MinAUTO, which uses BFGS if number of atoms is less than 500, CG otherwise.
    # max_steps (int) - The maximum number of steps to run the mimization for. Default is 1500.
    # energy_criterion (float) - The energy convergence criterion. Default is 5.0e-09.
    # gradient_criterion (float) - The gradient convergence criterion. Default is 0.05.
    # verbose (bool) - Printing verbosity. Default is False.
    # energy_no_electrostatics (bool) - Whether to turn off electrostatics. Default is False, i.e. to use electrostatics. NOTE: Can not be modified after the instance is created
    # energy_lj_repulsive_only (bool) - Whether to use only the repulsive portion of the Lennard-Jones term. Default is False. NOTE: Can not be modified after the instance is created
    # nonbonded_cutoff (float) - The cutoff for non-bonded interactions. Default is 14.0.
    #     maintain_planarity (bool) - Enable scaling of forces to artificially maintain planarity of certain sub-structures. This sets associated scale factors to their default values. Default is False. NOTE: Can not be modified after the instance is created
    # dielectric_constant (float) - The strength of the constant dielectric field used in the energy calculation. The default is 1.0 (vacuum).
    def xmlrpc_convert_nmr_lib(self,pdbString):
        result = threads.deferToThread(self.convert_nmr_lib,pdbString)
        return result

    def convert_nmr_lib(self,pdbString):
        tmpdir = tempfile.mkdtemp(prefix="nmr")
        input_file = os.path.join(tmpdir,"input.pdb")
        open(input_file,"w").write(pdbString)
        output_file = os.path.join(tmpdir,"output.lib")
        p = subprocess.Popen([NMR_CONVERT,input_file,output_file])
        p.communicate()
        result = open(output_file,"r").read()
        return result

    def xmlrpc_minimize_structure(self,mol_string_3d):
        result = threads.deferToThread(self.minimize_structure,mol_string_3d)
        return result

    def xmlrpc_minimize_structure_fixed(self, mol_string_3d, idList):
        result = threads.deferToThread(self.minimize_structure_fixed, mol_string_3d, idList)
        return result

    def xmlrpc_szybki(self,mol_string_3d):
        return threads.deferToThread(self.oe_minimize_structure,mol_string_3d)

    def xmlrpc_szybki_fixed(self,mol_string_3d, id_list_str):
        return threads.deferToThread(self.oe_minimize_structure_fix,mol_string_3d,id_list_str)

    def xmlrpc_mopac(self,mol_string_3d, formal_charge):
        return threads.deferToThread(self.mopac,mol_string_3d,formal_charge)

    def xmlrpc_canonicalizer(self,mol_string_2d):
        return threads.deferToThread(self.canonicalizer,mol_string_2d)

    def canonicalizer(self,mol_string_2d):
        standardizer = Standardizer()
        canonicalizer = TautomerCanonicalizer(max_tautomers=100)
        frament_chooser = LargestFragmentChooser()
        uncharger = Uncharger()
        mol = Chem.MolFromMolBlock(mol_string_2d)
        mol = frament_chooser.choose(mol)
        mol = uncharger.uncharge(mol)
        mol = standardizer.standardize(mol)
        mol = canonicalizer.canonicalize(mol)
        molfile = Chem.MolToMolBlock(mol)
        return molfile

    def mopac(self,mol_string_3d, formal_charge):
        tmpdir = tempfile.mkdtemp(prefix="mopac")
        input_file = os.path.join(tmpdir, "input.sdf")
        f = open(input_file,"w")
        f.write(mol_string_3d)
        f.close()
        mopac_input_file = os.path.join(tmpdir, "input.mae")
        p = subprocess.Popen([CONVERT_COMMAND,"-isd",input_file,"-omae",mopac_input_file])
        p.communicate()

        mopac_output = "%s" % os.path.join(tmpdir, "mopac")
        p = subprocess.Popen(["%s/run" % SCHRODINGER, "%s/mmshare-v6.1/bin/Linux-x86_64/semi_emp.py" % SCHRODINGER, "%s" % mopac_input_file, "-WAIT", "-jobname",
                              mopac_output, "-method", "AM1", "-keywords", "charge=%d mullik super vectors allvec" % formal_charge], stdout=subprocess.PIPE)
        p.communicate()
        mopac_result = MOPAC("%s.out"%mopac_output)
        homo = mopac_result.homo
        lumo = mopac_result.lumo

        result = {}
        result['HOMO'] = homo
        result['LUMO'] = lumo
        result['SAR'] = mopac_result.sar_dict
        return json.dumps(result)

    def minimize_structure(self, mol_string_3d):
        mol = next(structure.StructureReader.fromString(mol_string_3d,format=structure.SD))
        min_options = MinimizationOptions(gradient_convergence=0.05)
        minimizer.minimize_structure(mol,options=min_options)
        molstring_min = mol.writeToString(structure.SD)
        return molstring_min

    def minimize_structure_fixed(self, mol_string_3d, idListStr):
        if len(idListStr) == 0:
            return self.minimize_structure(mol_string_3d)
        else:
            atom_list = []
            for idStr in idListStr.split(","):
                atom_list.append(int(idStr)+1)
            mol = next(structure.StructureReader.fromString(mol_string_3d, format=structure.SD))
            atomIndices = mol.getAtomIndices()
            atom_list_to_optimize = []
            for atomId in atomIndices:
                if atomId not in atom_list:
                    atom_list_to_optimize.append(mol.atom[atomId])
            minimizer.minimize_substructure(mol, atom_list_to_optimize)
            return mol.writeToString(structure.SD)

    def oe_minimize_structure(self,mol_string_3d):
        tmpdir = tempfile.mkdtemp(prefix="szybki")
        input_file = os.path.join(tmpdir, "input.sdf")
        f = open(input_file,"w")
        f.write(mol_string_3d)
        f.close()
        output_file = os.path.join(tmpdir, "output.sdf")
        szybki_command = OE_MINIMIZE_LIGAND_SIMPLE%(os.path.join(OE_DIR,"szybki"),input_file,output_file)
        print(szybki_command)
        p = subprocess.Popen(szybki_command.split(),stdout=subprocess.PIPE)
        p.communicate()
        return open(output_file,"r").read()

    def oe_minimize_structure_fix(self, mol_string_3d, idListStr):
        if len(idListStr)==0:
            return self.oe_minimize_structure(mol_string_3d)
        else:
            mol_utilities = MolUtilities()
            mol = mol_utilities.convertToMol(mol_string_3d, OEFormat_SDF)
            tmpdir = tempfile.mkdtemp(prefix="szybki_fix")
            input_file = os.path.join(tmpdir, "input.sdf")
            output_file = os.path.join(tmpdir,"output.sdf")
            f = open(input_file, "w")
            f.write(mol_string_3d)
            f.close()
            fix_file = os.path.join(tmpdir, "fix_file")
            id_file = open(fix_file, "w")
            id_file.write("%s\n"%mol.GetTitle())
            idList = idListStr.split(",")
            if len(idList) > 0:
                for id in idList:
                    if len(id) > 0:
                        id_file.write("%s\n"%id)
            id_file.close()
            szybki_command = OE_MINIMIZE_LIGAND_FIX%(os.path.join(OE_DIR,"szybki"),fix_file, input_file,output_file)
            print(szybki_command)
            p = subprocess.Popen(szybki_command.split(),stdout=subprocess.PIPE)
            p.communicate()
            return open(output_file,"r").read()

    def getAvailableDockingGrids(self):
        grid_dirs = glob.glob("%s/*"%DOCK_GRID_DIR)
        gridNames = []
        for dir in grid_dirs:
            gridNames.append(os.path.basename(dir))
        return gridNames

    def xmlrpc_getLigandAndReceptor(self,gridName):
        result = threads.deferToThread(self.getLigandAndReceptor,gridName)
        return result

    def getLigandAndReceptor(self,gridName):
        try:
            pdbFile = "%s/%s/receptor.pdb"%(DOCK_GRID_DIR,gridName)
            ligandFile = "%s/%s/ligand.sdf"%(DOCK_GRID_DIR,gridName)
            pdb = open(pdbFile,"r").read()
            sdf = open(ligandFile,"r").read()
            return pdb,sdf
        except:
            return "Wrong receptor name!"

    def xmlrpc_getTemplate(self,template):
        result = threads.deferToThread(self.getTemplate,template)
        return result

    def getAvailableTemplates(self):
        template_dirs = glob.glob("%s/*"%TEMPLATE_DIR)
        templateNames = []
        for dir in template_dirs:
            templateNames.append(os.path.basename(dir))
        return templateNames

    def getTemplate(self,templateName):
        try:
            templateFile = "%s/%s"%(TEMPLATE_DIR,templateName)
            sdf = open(templateFile,"r").read()
            return sdf
        except:
            return "Wrong template name!"

    def xmlrpc_docking(self,grid_dir,ligandMolString, poses_per_ligand, docking_method, precision, use_ref = False, refLigandStr = None,
                       core_atoms = None, core_smarts = None):
        result = threads.deferToThread(self.docking, grid_dir,ligandMolString,poses_per_ligand,docking_method, precision, use_ref,
                                       refLigandStr, core_atoms, core_smarts)
        return result

    def xmlrpc_oedocking(self,grid_dir,ligandMolString,docking_method, numOfPoses):
        result = threads.deferToThread(self.oe_docking,grid_dir,ligandMolString,docking_method, numOfPoses)
        return result

    def xmlrpc_oe_rocs(self,queryMolString,ligandMolString,ewindow,rmsd):
        result = threads.deferToThread(self.oe_rocs, queryMolString,ligandMolString,ewindow,rmsd)
        return result

    def xmlrpc_shape_screen_flex(self, queryMolString, ligandMolString, ewindow, rmsd):
        result = threads.deferToThread(self.shape_screen_flex, queryMolString, ligandMolString, ewindow, rmsd)
        return result

    def xmlrpc_shape_screen_rigid(self, queryMolString, ligandMolString):
        result = threads.deferToThread(self.shape_screen_rigid, queryMolString, ligandMolString)
        return result

    def xmlrpc_oe_eon(self,queryMolString,ligandMolString,ewindow,rmsd):
        result = threads.deferToThread(self.oe_eon,queryMolString,ligandMolString,ewindow,rmsd)
        return result

    def xmlrpc_alogd(self, ligandMolString):
        result = threads.deferToThread(self.alogd, ligandMolString)
        return result

#OE_OMEGA_COMMAND = "%s -mpi_np 40 -in %s -out %s -maxConfRange 200,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600 -rangeIncrement 1 -addtorlib %s"
#OE_OMEGA_COMMAND_PARAM = "%s -mpi_np 40 -in %s -out %s -maxConfRange 200,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600 -rangeIncrement 1 -param %s -addtorlib %s -ewindow %f -rms %f "

    def xmlrpc_oeconfgen(self,ligandMolString, ewindow, rms):
        return threads.deferToThread(self.mmod_confgen, ligandMolString, ewindow, rms)

    def xmlrpc_freeform(self, inputSdf):
        return threads.deferToThread(self.freeform_hpc,inputSdf)

    def xmlrpc_Moka(self,inputSdf):
        result = threads.deferToThread(self.moka, inputSdf)
        return result

    @staticmethod
    def parseMokaDescResult(filename):
        ifs = oemolistream()
        ifs.open(filename)
        mol = OEGraphMol()
        result = {}
        while OEReadMolecule(ifs, mol):
            mol_result = {}
            pka_result = []
            pka_data = OEGetSDData(mol, "MoKa")
            if len(pka_data) > 0:
                p = parse("{molName} {covalent_hydration} - {numPka} {result}", pka_data)
                if p is not None:
                    numPka = int(p.named["numPka"])
                    r = p.named["result"]
                    arry = r.strip().split(" ")
                    for i in range(0, numPka):
                        n = i * 4
                        pka = pKa(arry[n], arry[n + 1], arry[n + 2], arry[n + 3])
                        pka_result.append(pka.getDict())
                    mol_result["pKa"] = pka_result
            logd_data = OEGetSDData(mol, "MoKa.LogD")
            if len(logd_data) > 0:
                lines = logd_data.split("\n")
                for line in lines:
                    p = parse("{ph}: {result}", line)
                    ph = p.named['ph']
                    mol_result["MoKa_LogD%s" % ph] = p.named["result"]
            logp_data = OEGetSDData(mol, "MoKa.LogP")
            if len(logp_data) > 0:
                mol_result["MoKa_LogP"] = logp_data
            mol_id = OEGetSDData(mol, "moka_id")
            result[mol_id] = mol_result
        ifs.close()
        return json.dumps(result)

    def moka(self, inputSdf):
        inputSdf = base64.decodebytes(bytes(inputSdf, 'utf-8')).decode('utf-8')
        tmpdir = tempfile.mkdtemp(prefix="moka_")
        moka_input = os.path.join(tmpdir,"moka_input.sdf")
        moka_output = os.path.join(tmpdir,"moka_output.sdf")
        f = open(moka_input,"w")
        f.write(inputSdf)
        f.close()
        p = subprocess.Popen([MOKA_COMMAND,"-s","moka_id","--show-logp","--show-logd=7.4","-o",moka_output,moka_input])
        p.communicate()
        jsonResult = GlideDockingServer.parseMokaDescResult(moka_output)
        return jsonResult


    def freeform_hpc(self,inputSdf):
        import multiprocessing
        manager = multiprocessing.Manager()
        inputSdf = base64.decodestring(bytes(inputSdf,'utf-8')).decode('utf-8')
        n_processors = 4
        trunks = MolUtilities().splitMols(inputSdf, n_processors)
        return_dict = manager.dict()
        jobs = []
        for sdf in trunks:
            p = multiprocessing.Process(target=self.freeform, args=(sdf,return_dict))
            jobs.append(p)
            p.start()
            print ("%s started"%p.name)

        for p in jobs:
            p.join()

        dict = {}
        for key in return_dict.keys():
            dict[key] = return_dict[key]
        return json.dumps(dict)

    def freeform(self, sdf, return_dict=None):
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(sdf)
        mol = OEGraphMol()
        ofs = oemolostream()
        ofs.SetFormat(OEFormat_SDF)
        ofs.openstring()
        resultDict = {}

        while OEReadMolecule(ifs, mol):
            moka_id = OEGetSDData(mol, "moka_id")
            mol.SetTitle(moka_id)
            OEWriteMolecule(ofs, mol)
        new_sdf = ofs.GetString()
        ifs.close()
        ofs.close()

        if len(new_sdf) > 0:
            tmpdir = tempfile.mkdtemp(prefix="freeform_")
            sdf_input = os.path.join(tmpdir, "input.sdf")

            prefix = os.path.join(tmpdir, "freeform")
            f = open(sdf_input, "w")
            f.write(new_sdf.decode('utf-8'))
            f.close()
            oe_freeform_cmd =  OE_FREEFORM_COMMAND%(os.path.join(OE_DIR,"freeform"),sdf_input,sdf_input,prefix)
            p = subprocess.Popen(oe_freeform_cmd.split(),stdout=subprocess.PIPE)
            p.communicate()

            Erel = "Erel"
            deltaG = "deltaG"
            local_strain = "Local Strain"
            global_strain = "Global Strain"

            log_file = os.path.join(tmpdir,"freeform.log")
            f = open(log_file, "r").read().splitlines()
            for line_no, line in enumerate(f):
                if line.strip() == "molname            Idx | Energy  Energy  RMSD  |  Energy  Energy  RMSD  |   Energy   Energy":
                    if len(f) > line_no + 2:
                        args = f[line_no + 2].split()
                        if len(args) == 10:
                            moka_id = args[0]
                            resultDict[moka_id] = {}
                            resultDict[moka_id][Erel] = args[3]
                            resultDict[moka_id][deltaG] = args[2]
                            resultDict[moka_id][local_strain] = args[5]
                            resultDict[moka_id][global_strain] = args[6]

        if return_dict is not None:
            for key in resultDict:
                if not key in return_dict:
                    return_dict[key] = resultDict[key]
        return json.dumps(resultDict)

    def oe_eon(self,refLigandMolString, ligandMolString, ewindow, rms):
        tmpdir = tempfile.mkdtemp(prefix="oerocs")
        query_file = os.path.join(tmpdir,"query.sdf")
        f = open(query_file,"w")
        f.write(refLigandMolString)
        f.close()

        input_file = os.path.join(tmpdir, "input.sdf")
        f = open(input_file,"w")
        f.write(ligandMolString)
        f.close()

        output_file = os.path.join(tmpdir,"output.oeb")
        confgen_command = CONF_GEN_PY%(input_file,output_file,ewindow,rms)
        p = subprocess.Popen(confgen_command.split(),stdout=subprocess.PIPE)
        p.communicate()
        #OE_ROCS_EON_INPUT = "%s/rocs -mpi_np 24 -bestHits 0 -query %s -dbase %s -eon_input -eon_input_size 0 -eon_input_file %s"
        # OE_EON = "%s/eon -mpi_np 24 -bestHits 0 -query %s -dbase %s -hitsfile %s -scoreonly"
        eon_input_file = os.path.join(tmpdir,"eon_input.oeb")
        rocs_command = OE_ROCS_EON_INPUT%(OE_DIR, query_file,output_file,eon_input_file)
        p = subprocess.Popen(rocs_command.split(),stdout=subprocess.PIPE)
        p.communicate()
        hits_file = os.path.join(tmpdir,"hits.sdf")
        eon_input_file = os.path.join(tmpdir,"eon_input_1.oeb")
        eon_command = OE_EON%(OE_DIR,query_file,eon_input_file,hits_file)
        p = subprocess.Popen(eon_command.split(),stdout=subprocess.PIPE)
        p.communicate()
        return open(hits_file,"r").read()

    #SHAPE_SCREEN_FLEX = "/opt/schrodinger/installations/default/shape_screen -WAIT -flex -flexSearchMethod thorough -flexMaxConfs %d -flexMaxRelEnergy %f -shape %s -screen %s -JOB %s"

    def shape_screen_flex(self, refLigandMolString, ligandMolString, ewindow, rmsd):
        tmpdir = tempfile.mkdtemp(prefix="ssflex")
        query_file = os.path.join(tmpdir,"query.sdf")
        f = open(query_file,"w")
        f.write(refLigandMolString)
        f.close()

        target_file = os.path.join(tmpdir,"target.sdf")
        s_writer = structure.StructureWriter(target_file, format=structure.SD)
        for st in structure.StructureReader.fromString(ligandMolString,format=structure.SD):
            #print(st.has3dCoords())
            if not st.has3dCoords():
                st.generate3dConformation()
            #print(st.writeToString(format=structure.SD))
            s_writer.append(st)
        s_writer.close()

        allMolStr = self.mmod_confgen(open(target_file,"r").read(),ewindow,rmsd,fast=True)
        return self.shape_screen_rigid(refLigandMolString,allMolStr)

        # """
        # job_name = os.path.join(tmpdir,"my_result")
        # result_fname = "%s_align.maegz"%job_name
        # sdf_name = "%s.sdf"%job_name
        # flag_fname = "%s_shape.okay"%job_name
        # shape_command = SHAPE_SCREEN_FLEX%(1000, ewindow, query_file, target_file, job_name)
        # p = subprocess.Popen(shape_command.split(), stdout=subprocess.PIPE)
        # p.communicate()
        # if os.path.exists(flag_fname):
        #     with structure.StructureWriter(sdf_name) as writer:
        #         with structure.StructureReader(result_fname) as reader:
        #             for st in reader:
        #                 writer.append(st)
        #         reader.close()
        #     writer.close()
        #     return open(sdf_name,"r").read()
        # return ""
        # """

    #SHAPE_SCREEN_RIGID = "/opt/schrodinger/installations/default/shape_screen -WAIT -distinct -shape %s -screen %s -JOB %s"
    def shape_screen_rigid(self, refLigandMolString, ligandMolString):
        tmpdir = tempfile.mkdtemp(prefix="ss_rigid")
        query_file = os.path.join(tmpdir,"query.sdf")
        f = open(query_file,"w")
        f.write(refLigandMolString)
        f.close()

        target_file = os.path.join(tmpdir,"target.sdf")
        s_writer = structure.StructureWriter(target_file, format=structure.SD)
        for st in structure.StructureReader.fromString(ligandMolString,format=structure.SD):
            #print(st.has3dCoords())
            if not st.has3dCoords():
                st.generate3dConformation()
            #print(st.writeToString(format=structure.SD))
            s_writer.append(st)
        s_writer.close()

        job_name = os.path.join(tmpdir,"my_result")
        result_fname = "%s_align.maegz"%job_name
        flag_fname = "%s_shape.okay"%job_name
        sdf_name = "%s.sdf" % job_name
        shape_command = SHAPE_SCREEN_RIGID%(query_file, target_file, job_name)
        p = subprocess.Popen(shape_command.split(), stdout=subprocess.PIPE)
        p.communicate()
        if os.path.exists(flag_fname):
            with structure.StructureWriter(sdf_name) as writer:
                with structure.StructureReader(result_fname) as reader:
                    for st in reader:
                        writer.append(st)
                reader.close()
            writer.close()
            return open(sdf_name,"r").read()
        else:
            print("Rigid failed.")
            return ""

    def oe_rocs(self,refLigandMolString, ligandMolString, ewindow, rms):
        tmpdir = tempfile.mkdtemp(prefix="oerocs")
        query_file = os.path.join(tmpdir,"query.sdf")
        f = open(query_file,"w")
        f.write(refLigandMolString)
        f.close()

        input_file = os.path.join(tmpdir, "input.sdf")
        f = open(input_file,"w")
        f.write(ligandMolString)
        f.close()

        output_file = os.path.join(tmpdir,"output.oeb")
        confgen_command = CONF_GEN_PY%(input_file,output_file,ewindow,rms)
        p = subprocess.Popen(confgen_command.split(),stdout=subprocess.PIPE)
        p.communicate()
        # OE_ROCS = "%s/rocs -mpi_np 24 -bestHits 0 -query %s -dbase %s -hitsfile %s"
        # OE_EON = "%s/eon -mpi_np 24 -bestHits 0 -query %s -dbase %s -hitsfile %s -scoreonly"
        hits_file = os.path.join(tmpdir,"hits.sdf")
        rocs_command = OE_ROCS%(OE_DIR,query_file,output_file,hits_file)
        p = subprocess.Popen(rocs_command.split(),stdout=subprocess.PIPE)
        p.communicate()
        return open(hits_file,"r").read()


    def mmod_confgen(self, ligandMolString, ewindow, rmsd, fast=False):
        #mol = next(structure.StructureReader.fromString(ligandMolString, format=structure.SD))
        tmpdir = tempfile.mkdtemp(prefix="mmod")
        input_fname = os.path.join(tmpdir,"input.mae")
        with structure.StructureWriter(input_fname,overwrite=True) as writer:
            for mol in structure.StructureReader.fromString(ligandMolString,format=structure.SD):
                writer.append(mol)
        writer.close()
        com_file = os.path.join(tmpdir,"mmod_csearch.com")
        if fast:
            MacroModelCmd.generate_fast_com(com_file,ewindow,rmsd)
        else:
            MacroModelCmd.generate_default_com(com_file, ewindow, rmsd)
        output_fname = os.path.join(tmpdir,"output.mae")
        output_sdf = os.path.join(tmpdir,"output.sdf")
        confgen_command = MACROMODEL%(os.path.basename(com_file))
        p = subprocess.Popen(confgen_command.split(),stdout=subprocess.PIPE, cwd=tmpdir)
        p.communicate()
        with structure.StructureWriter(output_sdf) as sd_writer:
            for st in structure.StructureReader(output_fname):
                sd_writer.append(st)
        sd_writer.close()
        return open(output_sdf,"r").read()

    #ALOGD_CMD = "python /opt/schrodinger/installations/default/mmshare-v6.1/python/common/ld_protocols/ld_alogD.py -ph %f %s"
    def alogd(self, ligandMolString, pH=7.4):
        #print(ligandMolString)
        tmpdir = tempfile.mkdtemp(prefix="logd")
        input_fname = os.path.join(tmpdir, "input.sdf")
        with structure.StructureWriter(input_fname,overwrite=True) as writer:
            for st in structure.StructureReader.fromString(ligandMolString,format=structure.SD):
                writer.append(st)
        alogd_command = ALOGD_CMD%(pH, input_fname)
        p = subprocess.Popen(alogd_command.split(),stdout=subprocess.PIPE, cwd=tmpdir)
        p.communicate()
        result_csv = os.path.join(tmpdir,"results.csv")
        result_dict = {}
        with open(result_csv) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                try:
                    id = row['Corporate ID']
                    logd = float(row['AlogD Custom pH'])
                    result_dict[id] = logd
                except:
                    pass
        print(result_dict)
        return json.dumps(result_dict)

    def oeconfgen(self,ligandMolString, ewindow, rms):
        molUtilities = MolUtilities()
        #param_file = os.path.join(OE_PARAM_DIR,"confgen.param")
        #torsion_file = os.path.join(OE_PARAM_DIR,"my_torlib.txt")
        tmpdir = tempfile.mkdtemp(prefix="oeomega")
        input_file = os.path.join(tmpdir, "input.sdf")
        f = open(input_file,"w")
        f.write(ligandMolString)
        f.close()
        output_file = os.path.join(tmpdir, "output.sdf")
#       omega_command = OE_OMEGA_COMMAND_PARAM%(os.path.join(OE_DIR,"omega2"),input_file,output_file,param_file,torsion_file, ewindow, rms)
        omega_command = OE_OMEGA_COMMAND_PARAM%(os.path.join(OE_DIR,"omega2"),input_file,output_file,ewindow, rms)
        print(omega_command)
        p = subprocess.Popen(omega_command.split(),stdout=subprocess.PIPE)
        p.communicate()
        output_file2 = os.path.join(tmpdir,"output2.sdf")
        molUtilities.generate_relative_energies_oeomega(output_file,output_file2)
        return open(output_file2,"r").read()

    def oe_docking(self, grid_dir, ligandMolString, docking_method, numOfPoses): #docking_method = ["hybrid","posit"]
        # OE_OMEGA_INPUT = "ligands.sdf"
        # OE_DOCKING_INPUT = "ligands_multiconfs.sdf"
        # OE_DOCKING_GRID = "receptor.oeb.gz"
        # OE_MAKE_RECEPTOR_COMMAND = "receptor_setup -protein %s -bound_ligand %s -receptor %s"
        # OE_DOCKING_COMMAND = "%s -receptor %s -dbase %s -docked_molecule_file %s"
        # OE_OMEGA_COMMAND = "omega2 -in %s -out %s -maxConfRange '200,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600' -rangeIncrement 1"

        molUtilities = MolUtilities()
        torsion_file = os.path.join(OE_PARAM_DIR,"my_torlib.txt")
        recep_file = os.path.join(DOCK_GRID_DIR,grid_dir,"receptor.pdb")
        reference_ligand_file = os.path.join(DOCK_GRID_DIR,grid_dir,"ligand.sdf")
        tmpdir = tempfile.mkdtemp(prefix="oe_docking")
        receptor_grid = os.path.join(tmpdir,OE_DOCKING_GRID)
        oe_receptor_grid_command = OE_MAKE_RECEPTOR_COMMAND % (os.path.join(OE_DIR,"receptor_setup"),recep_file, reference_ligand_file, receptor_grid)
        p = subprocess.Popen(oe_receptor_grid_command.split(), stdout=subprocess.PIPE)
        p.communicate()
        molList = molUtilities.convertToMolList(ligandMolString)

        n = len(molList)
        minimum = 1
        maximum = 24

        m = n/5

        if m < minimum:
            m = minimum
        if m > maximum:
            m = maximum

        ofs = oemolostream()
        omega_input = os.path.join(tmpdir, OE_OMEGA_INPUT)
        ofs.open(omega_input)
        for mol in molList:
            OEWriteMolecule(ofs,mol)
        ofs.close()
        oe_docking_input = os.path.join(tmpdir,OE_DOCKING_INPUT)

        omega_command = OE_OMEGA_COMMAND%(os.path.join(OE_DIR,"omega2"),omega_input,oe_docking_input,torsion_file)
        p = subprocess.Popen(omega_command.split(),stdout=subprocess.PIPE)
        p.communicate()

        oedocking_output = os.path.join(tmpdir,"dockedLigands.sdf")
        if docking_method.startswith("posit"):
            docking_method,clashes = docking_method.split("_")
            oedocking_command = OE_POSIT_COMMAND%(os.path.join(OE_DIR,docking_method),receptor_grid,oe_docking_input,numOfPoses,clashes, oedocking_output)
        else: #hybrid
            oedocking_command = OE_HYBRID_COMMAND%(os.path.join(OE_DIR,docking_method),receptor_grid,oe_docking_input,numOfPoses,oedocking_output)
        p = subprocess.Popen(oedocking_command.split(),stdout=subprocess.PIPE)
        p.communicate()

        try:
            ligand = open(os.path.join(tmpdir,"dockedLigands.sdf"),"r").read()
        except:
            ligand = ""

        try:
            receptor = open(recep_file,"r").read()
        except:
            receptor = ""

        try:
            reference = open(reference_ligand_file,"r").read()
        except:
            reference = ""

        rs = {}
        rs['ligand'] = ligand
        rs['receptor'] = receptor
        rs['reference'] = reference
        shutil.rmtree(tmpdir,ignore_errors=True)
        return json.dumps(rs)

    def docking(self,grid_dir,ligandMolString, poses_per_ligand, docking_method, precision, use_ref,
                refLigandStr, core_atoms, core_smarts):

        molUtilities = MolUtilities()
        gridfile = os.path.join(DOCK_GRID_DIR,grid_dir,"grid.zip")
        recep_file = os.path.join(DOCK_GRID_DIR,grid_dir,"receptor.pdb")
        reference_ligand_file = os.path.join(DOCK_GRID_DIR,grid_dir,"ligand.sdf")
        tmpdir = tempfile.mkdtemp(prefix="glide_docking")
        shutil.copy(gridfile,tmpdir)
        shutil.copy(recep_file,tmpdir)
        molList = molUtilities.convertToMolList(ligandMolString)

        n = len(molList)
        minimum = 1
        maximum = 2

        m = n/5

        if m < minimum:
            m = minimum
        if m > maximum:
            m = maximum

        ofs = oemolostream()
        ofs.open(os.path.join(tmpdir,"ligand.sdf"))
        for mol in molList:
            # mol = molUtilities.prepareMoleculeForGlide(mol)
            OEWriteMolecule(ofs,mol)
        ofs.close()


        glide_input = os.path.join(tmpdir,"glide_docking.in")
        inputfile = open(glide_input,"w")
        inputfile.write("GRIDFILE\t%s\n"%os.path.join(tmpdir,"grid.zip"))
        inputfile.write("DOCKING_METHOD\t%s\n"%docking_method)
        inputfile.write("PRECISION\t%s\n"%precision)
        inputfile.write("LIGANDFILE\t%s\n"%(os.path.join(tmpdir,"ligand.sdf")))
        inputfile.write("POSES_PER_LIG\t%s\n"%poses_per_ligand)
        if use_ref:
            if refLigandStr is not None and core_atoms is not None and core_smarts is not None:
                #    USE_REF_LIGAND  True
                #    REF_LIGAND_FILE refLigandH.sdf
                #    CORE_DEFINITION smarts
                #    CORE_ATOMS      1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21
                #    CORE_SMARTS     O=CNc1cccc(c1)C1CCCN(C1)c1ncnc2[nH]ccc12
                refFileName = os.path.join(tmpdir, "reference.sdf")
                refLigandFile = open(refFileName, "w")
                refLigandFile.write(refLigandStr)
                refLigandFile.close()
                inputfile.write("USE_REF_LIGAND\tTrue\n")
                inputfile.write("REF_LIGAND_FILE\t%s\n"%refFileName)
                inputfile.write("CORE_DEFINITION\tsmarts\n")
                inputfile.write("CORE_ATOMS\t%s\n"%core_atoms)
                inputfile.write("CORE_SMARTS\t%s\n"%core_smarts)
                inputfile.write("CORE_RESTRAIN\tTrue\n")
        inputfile.close()

        glide_command = os.path.join(SCHRODINGER,"glide")
        os.chdir(tmpdir) #not thread safe, but time window is very short, so should be ok.
        p = subprocess.Popen([glide_command,"-WAIT","-NJOBS","%d"%m, "-HOST","localhost:%d"%m,glide_input],stdout=subprocess.PIPE)
        p.communicate()
        p = subprocess.Popen([SUBSET_COMMAND,"-n", "2:", os.path.join(tmpdir,"glide_docking_pv.maegz"), os.path.join(tmpdir,"ligand.mae")])
        p.communicate()
        p = subprocess.Popen([CONVERT_COMMAND,"-imae",os.path.join(tmpdir,"ligand.mae"),"-osd",os.path.join(tmpdir,"dockedLigand.sdf")])
        p.communicate()

        try:
            ligand = open(os.path.join(tmpdir,"dockedLigand.sdf"),"r").read()
        except:
            ligand = ""

        try:
            receptor = open(recep_file,"r").read()
        except:
            receptor = ""

        try:
            reference = open(reference_ligand_file,"r").read()
        except:
            reference = ""

        rs = {}
        rs['ligand'] = ligand
        rs['receptor'] = receptor
        rs['reference'] = reference
        return json.dumps(rs)

if __name__ == "__main__":
    from twisted.internet import reactor
    glideServer = GlideDockingServer()
    reactor.listenTCP(9527,server.Site(glideServer)) #production
    #reactor.listenTCP(9999,server.Site(glideServer))  #testing
    reactor.run()

