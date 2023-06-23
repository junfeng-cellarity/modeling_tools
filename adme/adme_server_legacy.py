__author__ = 'jfeng1'
import os,subprocess
import tempfile
import sys
from adme_utilities import *
# sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
# from glide.glide_utilities import *
import json
from parse import *
from twisted.web import xmlrpc, server
from twisted.internet import threads
from openeye.oechem import *
from openeye.oegraphsim import *
import base64
from rpy2 import *
import rpy2.robjects as robjects
import pandas
import numpy
import traceback
from sqlitedict import SqliteDict
import pandas as pd
from check_inventory_evotec import *
import math
from sklearn import svm
from sklearn.metrics import confusion_matrix,cohen_kappa_score
from sklearn.metrics import r2_score
from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel
from sklearn import cross_validation
from sklearn.preprocessing.data import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.metrics.classification import classification_report
from sklearn.metrics.classification import accuracy_score
from sklearn.feature_selection import RFECV
from sklearn.metrics import make_scorer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.ensemble.forest import RandomForestRegressor
import shutil


MOKA_COMMAND = "/home/jfeng1/moldiscovery/moka-2.6.1-x86_64/moka_cli"
BLABBER_COMMAND = "/home/jfeng1/moldiscovery/moka-2.6.1-x86_64/blabber_sd"
hERG_COMMAND = "/UserUCdd/ienyedy/OE/openeye/bin/vsepr_krige_mol"
PPB_COMMAND = "/UserUCdd/ienyedy/OE/openeye/bin/krige_molecule"
VOLSURF_COMMAND = "/home/jfeng1/moldiscovery/vsplus-64/vsplus-cli"
CLOGP_COMMAND = "/pkg/clogp/5.4/Linux/i686/bin/clogp"
KRIG_COMMAND = "/home/jfeng1/Kriging/openeye/bin/krige_mol"
FIXPKA_COMMAND = "fixpka"
KRIG_TRAINING = "/DbUCdd/databases/ADME/Krige_Training_Sets"
FREEFORM_COMMAND = "freeform"
INSILICO_ADME = os.path.join(os.path.dirname(os.path.abspath(__file__)),"physchem_adme_predict.py")
INSILICO_ADME_COMMAND = "%s %s xxxx mc MACCS_FCFP4_rdMolDes RF_NEW xxxx"

# /UserUCdd/ienyedy/OE/openeye/bin/vsepr_krige_mol -in junk.sdf -out output.sdf
# -training_Set /DbUCdd/databases/ADME/hERG/Kriging_model_on_PP_insilico/trainning_set_1884.oeb
# -training_tag IC50_uM -response_tag Predicted_log_IC50 -take_log_of_training_data true
# /UserUCdd/ienyedy/OE/openeye/bin/krige_molecule  -in junk.sdf -out output.sdf
# -predicted_tag Pred_Fu -training_molecules /DbUCdd/databases/ADME/Krige_ADME/Plasma_Protein_Binding/PPB_dataset_2098.sdf
# -response_tag PPB_HUMAN_unbound -take_log_of_response true

command_file = """set FIELD_PARAM = STATIC;
set GRID_SPACING = 0.5;
set PROTONATION = AS_IS;
set GRUB_PATH = "/home/jfeng1/moldiscovery/vsplus-64/grub.dat";
set VOLSURF_LIBRARY = "/DbUCDD/databases/ADME/Clearance/vs_models";
new attribute "id";
import SDF "%s" name by attribute "id" assign ("moka_id" = "id");
export datamatrix as CSV "%s" using separator ",";"""


clearance_command_file = """set FIELD_PARAM = STATIC;
set GRID_SPACING = 0.5;
set PROTONATION = AS_IS;
set GRUB_PATH = "/home/jfeng1/moldiscovery/vsplus-64/grub.dat";
set VOLSURF_LIBRARY = "/DbUCDD/databases/ADME/Clearance/vs_models";
new attribute "id";
import SDF "%s" name by attribute "id" assign ("moka_id" = "id");
project on library model "hCLpathway483_42.PCA2";
export "hCLpathway483_42.PCA2" predicted scores as CSV "%s" using separator ",";"""


VDSS_DESCRIPTOR_NAMES = ["MolName","Total_Anionic_Charge","Total_Cationic_Charge","Species_abundance",
                         "LogD8","LogD8.5","LogD9","LogD9.5",
                         "LogD10","LogD10.5","LogD11",
                         "LogD11.5","LogD12","LogD12.5",
                         "VolSurf_W7", "VolSurf_W8", "VolSurf_WO1", "VolSurf_WO2", "VolSurf_WO3",
                         "VolSurf_WO4", "VolSurf_WN5", "VolSurf_WN6", "VolSurf_CW7", "VolSurf_CW8",
                         "VolSurf_DIFF", "VolSurf_PSA", "VolSurf_%FU4", "VolSurf_%FU5",
                         "VolSurf_%FU6", "VolSurf_%FU7", "VolSurf_%FU8", "VolSurf_%FU9",
                         "VolSurf_%FU10","SmartsString14"]

CLEARANCE_MECHANISM_DESCRIPOR_NAMES = ["MolName","VolSurf_V","VolSurf_S","VolSurf_R","VolSurf_G","VolSurf_IW1","VolSurf_IW2","VolSurf_IW3","VolSurf_IW4","VolSurf_ID1","VolSurf_ID2",
                                       "VolSurf_ID3","VolSurf_ID4","VolSurf_HL1","VolSurf_HL2","VolSurf_A","VolSurf_CP","VolSurf_POL","VolSurf_MW","VolSurf_FLEX",
                                       "VolSurf_FLEX_RB","VolSurf_NCC","VolSurf_LOGP n-Oct","VolSurf_LOGP c-Hex","VolSurf_PSA","VolSurf_HSA","VolSurf_PSAR",
                                       "VolSurf_PHSAR","VolSurf_LgD5","VolSurf_LgD6","VolSurf_LgD7","VolSurf_LgD7.5","VolSurf_LgD8","VolSurf_LgD9","VolSurf_LgD10",
                                       "VolSurf_AUS7.4","VolSurf_%FU4","VolSurf_%FU5","VolSurf_%FU6","VolSurf_%FU7","VolSurf_%FU8","VolSurf_%FU9","VolSurf_%FU10"]
#
#
# CLEARANCE_META_VOLSURF = ["VolSurf_W7","VolSurf_D3","VolSurf_D4","VolSurf_WO1","VolSurf_WO2","VolSurf_WN6","VolSurf_IW4","VolSurf_A","VolSurf_CP","VolSurf_FLEX","VolSurf_LOGP n-Oct",
#                           "VolSurf_LgD5","VolSurf_LgD8","VolSurf_LgD9","VolSurf_LgD10","VolSurf_AUS7.4","VolSurf_%FU4","VolSurf_%FU5","VolSurf_DRDRAC","VolSurf_DRACAC"]
#
# CLEARANCE_META_SMARTS = ["[CH](C)=C","[O]=S","[CH3]c","[CH2](C(C)N)c(c)c","[C](O)(=O)CC",
#                          "[c](Cl)(cc)cc","[OH]C(c)=O","[CH3]n","[CH](CC)=CC","[c](O)(cc)c(O)c",
#                          "[C](=O)(CC)OC","[c](c(O)c)(c(=O)c)c(c)c","[NH]=C","[c]1(cc)c(C)cnc1c",
#                          "[CH2](C(C)C)C(C)C","[c](C(C)C)(cc)cc","[O](c)c","[CH](C)(c)O","[cH](c(N)c)c(O)c","[c](NC)(cc)c(C)c","[C](C)(NC)=C(C)C",
#                          "[CH3]C(=C)N","[CH2](Cc)C(C)N","[NH](C(c)C)C(C)=O","[S](C)(N)(=O)=O","[CH2](CC)n(c)c",
#                          "[c](=O)(nc)n(C)c","[CH](NC)(C(C)O)C(N)=O","[cH]1ccsc1","[nH]1c(C)ccc1c",
#                          "[c](C(c)=N)(cc)c(N)c","[N](C)(C)N","[n](c(N)n)c(c)c","[CH](C=C)=C(C)C",
#                          "[c](-c(c)n)(cc)cc","[CH2]1C(C)CCC1O","[O](CC)C(C)C","[c](C(C)(C)O)(cc)c(C)c","[CH](CC)(CC)N(C)C",
#                          "[c](C(c)=N)(cc)cc","[CH](C)(CC)CC","[c](=O)(c(C)c)n(C)c","[c](N)(nc)c(c)c",
#                          "[NH](C(N)=O)S(c)(=O)=O","[c](C(c)C)(cc)c(C)c","[C](O)(=O)C(C)C","[CH2](N=C)C(N)=O",
#                          "[cH](cc)c(N)n","[CH2](CC)C(C)(C)c","[CH](CC)(c(c)c)c(c)c","[c](OC)(cc)c(-c)c","[CH]1(CCCO1)n(c)c",
#                          "[c]1(C)nccn1C","[N]1(C)C(C)CCC1C","[cH](cc)c(Br)c","[CH2](CO)c(c)c",
#                          "[n](C)(c(=O)c)c(=O)n","[CH3]N(C)S","[C](c)(=C)O","[cH](cc)c(S)c"]
#
# CLEARANCE_RENAL_VOLSURF = ["VolSurf_G","VolSurf_W4","VolSurf_D2","VolSurf_D3","VolSurf_D4","VolSurf_D5","VolSurf_WO1","VolSurf_WO2","VolSurf_IW2","VolSurf_CW4","VolSurf_CW5","VolSurf_ID3","VolSurf_ID4",
#                            "VolSurf_CD1","VolSurf_CD4","VolSurf_CD5","VolSurf_CD6","VolSurf_CD7","VolSurf_CD8","VolSurf_HL1","VolSurf_HL2","VolSurf_A","VolSurf_NCC","VolSurf_PSAR",
#                            "VolSurf_LgD7","VolSurf_LgD7.5","VolSurf_LgD8","VolSurf_LgD9","VolSurf_LgD10","VolSurf_%FU5","VolSurf_%FU6","VolSurf_%FU7","VolSurf_%FU8","VolSurf_DRDRAC"]
#
# CLEARANCE_RENAL_SMARTS = ["[cH](cc)c(S)c","[cH](c)c","[c](C)(c)c","[n](C)(c)c","[O](C(C)C)C(C)=O","[O]=N",
#                           "[cH](cc)c(F)c","[CH2](C)CC","[n](C)(c)n","[CH2](CN)NC","[N](C)=C","[NH](C)S",
#                           "[C](C)(=O)Nc","[c](c(O)c)(c(c)c)n(C)c","[cH](c(N)c)c(c)n",
#                           "[CH](N)(CC)C(O)O","[c](S)(n)s","[N](C)(C)Cc","[c](=O)(nc)c(C)c",
#                           "[C]1(=O)CCN1C","[n](c(N)n)c(c)n"]
#
# META_MODEL_COEFFS = {'[n](C)(c(=O)c)c(=O)n': -0.532, '[C](c)(=C)O': -0.703, '[N](C)(C)N': -0.3142,
#                           'VolSurf_LOGP n-Oct': -0.02285, '[cH](c(N)c)c(O)c': -0.3473, '[C](O)(=O)C(C)C': -0.8262,
#                           '[c](C(c)=N)(cc)cc': -0.2225, '[CH2]1C(C)CCC1O': 0.6233, '[S](C)(N)(=O)=O': 0.3997,
#                           '[CH2](C(C)N)c(c)c': 0.1091, '[CH2](CC)C(C)(C)c': 0.3662, '[c](Cl)(cc)cc': -0.02195,
#                           '[c](=O)(c(C)c)n(C)c': -0.5144, '[CH](CC)(CC)N(C)C': 0.3357, '[c]1(C)nccn1C': -0.4456,
#                           '[CH]1(CCCO1)n(c)c': 0.2999, '[c](O)(cc)c(O)c': 0.1752, '[N]1(C)C(C)CCC1C': 0.2787,
#                           '[OH]C(c)=O': -0.2136, '[C](=O)(CC)OC': 0.232, '[CH3]C(=C)N': 0.0631, '[c](NC)(cc)c(C)c': -0.2779,
#                           '[c](-c(c)n)(cc)cc': 0.5068, '[n](c(N)n)c(c)c': -0.2328, '[CH](CC)=CC': -0.157, '[c](C(c)C)(cc)c(C)c': 0.3581,
#                           '[CH2](C(C)C)C(C)C': -0.2468, '[CH](CC)(c(c)c)c(c)c': -0.7462, 'VolSurf_WO2': 0.0008443, 'VolSurf_WO1': 0.0005021,
#                           '[CH3]c': -0.123, 'VolSurf_%FU4': -0.0007448, 'VolSurf_%FU5': -0.0005609, '[CH3]n': -0.07349, '[cH](cc)c(Br)c': -0.3183,
#                           '[nH]1c(C)ccc1c': -0.5472, '[c]1(cc)c(C)cnc1c': 0.1629, 'VolSurf_LgD9': 0.01551, 'VolSurf_LgD8': 0.008123,
#                           'VolSurf_DRACAC': -0.001847, '[cH]1ccsc1': 0.1819, '[C](C)(NC)=C(C)C': 0.0631, '[c](C(c)=N)(cc)c(N)c': -0.158,
#                           'VolSurf_A': -0.01914, 'VolSurf_DRDRAC': -0.00227, 'VolSurf_LgD5': -0.03726, '[c](OC)(cc)c(-c)c': 0.5185, '[CH](C)=C': -0.02948,
#                           '[cH](cc)c(N)n': -0.625, 'VolSurf_AUS7.4': 0.01186, '[CH](NC)(C(C)O)C(N)=O': -0.2645, '[c](C(C)C)(cc)cc': -0.2397, '[O](CC)C(C)C': -0.7625,
#                           '[NH]=C': 0.2205, '[CH](C)(c)O': -0.5027, '[CH2](CO)c(c)c': -0.3782, '[CH2](Cc)C(C)N': 0.3098, '[O](c)c': -0.6415, '[CH2](N=C)C(N)=O': -0.2876,
#                           'VolSurf_D4': -0.002212, '[CH3]N(C)S': -0.703, '[c](=O)(nc)n(C)c': 0.177, 'VolSurf_D3': -0.001015, '[c](c(O)c)(c(=O)c)c(c)c': 0.25,
#                           '[C](O)(=O)CC': 0.2586, '[c](N)(nc)c(c)c': -0.2328, '[O]=S': -0.117, '[NH](C(N)=O)S(c)(=O)=O': -0.5355, 'VolSurf_LgD10': 0.0151,
#                           'VolSurf_WN6': -0.001014, '[c](C(C)(C)O)(cc)c(C)c': -0.7784, '[cH](cc)c(S)c': -0.04735,
#                           '[CH2](CC)n(c)c': 0.3293, 'VolSurf_W7': -0.0008883, '[CH](C)(CC)CC': -0.5229,
#                           'VolSurf_IW4': -0.1014, '[CH](C=C)=C(C)C': -0.3753, '[NH](C(c)C)C(C)=O': 0.3373,
#                           'VolSurf_CP': 0.003049, 'VolSurf_FLEX': 0.02683}
#
# RENAL_MODEL_COEFFS = {'VolSurf_NCC': -0.03145, 'VolSurf_LgD7.5': 0.007247, 'VolSurf_LgD8': 0.007011, 'VolSurf_ID3': 0.2065,
#                       'VolSurf_ID4': 0.1985, '[cH](c)c': -0.01384, '[N](C)(C)Cc': 0.2756, '[n](C)(c)n': -0.1634,
#                       '[c](S)(n)s': -0.3203, '[n](C)(c)c': 0.1067, '[C]1(=O)CCN1C': 0.2714, '[c](=O)(nc)c(C)c': 0.266,
#                       '[NH](C)S': 0.2556, '[O](C(C)C)C(C)=O': 0.4339, '[CH2](C)CC': -0.2385, 'VolSurf_%FU8': -0.0003735,
#                       'VolSurf_WO2': 0.0001443, 'VolSurf_WO1': 7.4e-05, '[CH2](CN)NC': 0.1266, '[c](c(O)c)(c(c)c)n(C)c': -0.1605,
#                       'VolSurf_%FU5': -0.0009373, 'VolSurf_%FU6': -0.0007606, 'VolSurf_%FU7': -0.0007033, '[C](C)(=O)Nc': -0.1337,
#                       'VolSurf_G': -0.0957, 'VolSurf_HL1': 7.87e-05, 'VolSurf_A': -0.01034, 'VolSurf_DRDRAC': -0.000649,
#                       '[n](c(N)n)c(c)n': 0.1374, 'VolSurf_LgD7': 0.007859, 'VolSurf_PSAR': -0.2119, 'VolSurf_CD1': 0.2061,
#                       'VolSurf_CD6': 1.798, 'VolSurf_CD7': 2.866, 'VolSurf_CD4': 1.121, 'VolSurf_CD5': 1.477, 'VolSurf_CD8': 3.283,
#                       'VolSurf_LgD9': 0.006609, 'VolSurf_D4': -0.00234, 'VolSurf_D5': -0.003013, 'VolSurf_HL2': 0.0001944, '[c](C)(c)c': -0.01987,
#                       'VolSurf_D2': -0.0006458, 'VolSurf_D3': -0.001759, '[cH](cc)c(F)c': -0.4295, '[cH](cc)c(S)c': -0.3639, 'VolSurf_CW5': -0.09199,
#                       'VolSurf_CW4': -0.07129, 'VolSurf_LgD10': 0.005995, '[O]=N': 0.07896, '[cH](c(N)c)c(c)n': 0.1853, 'VolSurf_W4': -4.23e-05,
#                       '[N](C)=C': 0.3579, '[CH](N)(CC)C(O)O': -0.1769, 'VolSurf_IW2': 0.1432}

CLEARANCE_META_VOLSURF = [
    "VolSurf_W7",
    "VolSurf_D3",
    "VolSurf_D4",
    "VolSurf_WO1",
    "VolSurf_WO2",
    "VolSurf_WN6",
    "VolSurf_IW4",
    "VolSurf_A",
    "VolSurf_CP",
    "VolSurf_FLEX",
    "VolSurf_LOGP n-Oct",
    "VolSurf_LgD5",
    "VolSurf_LgD8",
    "VolSurf_LgD9",
    "VolSurf_LgD10",
    "VolSurf_AUS7.4",
    "VolSurf_%FU4",
    "VolSurf_%FU5",
    "VolSurf_DRDRAC",
    "VolSurf_DRACAC"]

CLEARANCE_META_SMARTS = [
    "[CH](C)=C",
    "[O]=S",
    "[c]([CH](C)C)([cH]c)[cH]c",
    "[c]([CH](c)C)([cH]c)c(C)c",
    "[C]([CH3])([NH]C)=C(C)C",
    "[c]([NH]C)([cH]c)c(C)c",
    "[c]([NH2])(nc)c(c)c",
    "[c]([OH])([cH]c)c(O)c",
    "[C]([OH])(=O)[CH](C)C",
    "[C]([OH])(=O)[CH2]C",
    "[C](=O)([CH2]C)OC",
    "[c](=O)(c(C)c)n(C)c",
    "[c](=O)(nc)n(C)c",
    "[c](C(C)(C)O)([cH]c)c(C)c",
    "[c](C(c)=N)([cH]c)[cH]c",
    "[c](C(c)=N)([cH]c)c(N)c",
    "[c](-c(c)n)([cH]c)[cH]c",
    "[c](c(O)c)(c(=O)c)c(c)c",
    "[C](c)(=C)O",
    "[c](Cl)([cH]c)[cH]c",
    "[c](OC)([cH]c)c(-c)c",
    "[c]1(C)nccn1C",
    "[c]1(cc)c(C)cnc1c",
    "[CH]([CH]=C)=C(C)C",
    "[cH]([cH]c)c(Br)c",
    "[cH]([cH]c)c(N)n",
    "[cH]([cH]c)c(S)c",
    "[CH]([CH2]C)([CH2]C)N(C)C",
    "[CH]([CH2]C)(c(c)c)c(c)c",
    "[CH]([CH2]C)=[CH]C",
    "[CH]([CH3])([CH2]C)[CH2]C",
    "[CH]([NH]C)([CH](C)O)C(N)=O",
    "[cH](c(N)c)c(O)c",
    "[CH](C)(c)O",
    "[CH]1(CCCO1)n(c)c",
    "[cH]1ccsc1",
    "[CH2]([CH](C)C)[CH](C)C",
    "[CH2]([CH](C)N)c(c)c",
    "[CH2]([CH2]c)[CH](C)N",
    "[CH2]([CH2]C)C(C)(C)c",
    "[CH2]([CH2]C)n(c)c",
    "[CH2]([CH2]O)c(c)c",
    "[CH2](N=C)C(N)=O",
    "[CH2]1C(C)CCC1O",
    "[CH3]C(=C)N",
    "[CH3]c",
    "[CH3]N(C)S",
    "[CH3]n",
    "[n]([CH3])(c(=O)c)c(=O)n",
    "[n](c(N)n)c(c)c",
    "[N](C)(C)N",
    "[N]1(C)C(C)CCC1C",
    "[NH]([CH](c)C)C(C)=O",
    "[NH](C(N)=O)S(c)(=O)=O",
    "[NH]=C",
    "[nH]1c(C)ccc1c",
    "[O]([CH2]C)[CH](C)C",
    "[O](c)c",
    "[OH]C(c)=O",
    "[S](C)(N)(=O)=O"]

CLEARANCE_RENAL_VOLSURF = [
    "VolSurf_G",
    "VolSurf_W4",
    "VolSurf_D2",
    "VolSurf_D3",
    "VolSurf_D4",
    "VolSurf_D5",
    "VolSurf_WO1",
    "VolSurf_WO2",
    "VolSurf_IW2",
    "VolSurf_CW4",
    "VolSurf_CW5",
    "VolSurf_ID3",
    "VolSurf_ID4",
    "VolSurf_CD1",
    "VolSurf_CD4",
    "VolSurf_CD5",
    "VolSurf_CD6",
    "VolSurf_CD7",
    "VolSurf_CD8",
    "VolSurf_HL1",
    "VolSurf_HL2",
    "VolSurf_A",
    "VolSurf_NCC",
    "VolSurf_PSAR",
    "VolSurf_LgD7",
    "VolSurf_LgD7.5",
    "VolSurf_LgD8",
    "VolSurf_LgD9",
    "VolSurf_LgD10",
    "VolSurf_%FU5",
    "VolSurf_%FU6",
    "VolSurf_%FU7",
    "VolSurf_%FU8",
    "VolSurf_DRDRAC"]

CLEARANCE_RENAL_SMARTS = [
    "[cH]([cH]c)c(S)c",
    "[cH](c)c",
    "[c](C)(c)c",
    "[n](C)(c)c",
    "[O]([CH](C)C)[CH2](C)=O",
    "[O]=N",
    "[cH]([cH]c)c(F)c",
    "[CH2]([CH3])[CH2]C",
    "[n](C)(c)n",
    "[CH2]([CH2]N)[NH]C",
    "[N](C)=C",
    "[NH](C)S",
    "[C]([CH3])(=O)[NH]c",
    "[c](c(O)c)(c(c)c)n(C)c",
    "[cH](c(N)c)c(c)n",
    "[CH]([NH2])([CH2]C)[CH](O)O",
    "[c](S)(n)s",
    "[N]([CH3])([CH3])[CH2]c",
    "[c](=O)(nc)c(C)c",
    "[C]1(=O)CCN1C",
    "[n](c(N)n)c(c)n"]

META_MODEL_COEFFS = {
    '[CH](C)=C':  -0.03225,
    '[O]=S':  -0.09293,
    '[c]([CH](C)C)([cH]c)[cH]c':  -0.281,
    '[c]([CH](c)C)([cH]c)c(C)c':  0.2908,
    '[C]([CH3])([NH]C)=C(C)C':  0.1034,
    '[c]([NH]C)([cH]c)c(C)c':  -0.3599,
    '[c]([NH2])(nc)c(c)c':  -0.2894,
    '[c]([OH])([cH]c)c(O)c':  0.2773,
    '[C]([OH])(=O)[CH](C)C':  -0.7765,
    '[C]([OH])(=O)[CH2]C':  0.395,
    '[C](=O)([CH2]C)OC':  0.4299,
    '[c](=O)(c(C)c)n(C)c':  -0.5968,
    '[c](=O)(nc)n(C)c':  0.3724,
    '[c](C(C)(C)O)([cH]c)c(C)c':  -0.3343,
    '[c](C(c)=N)([cH]c)[cH]c':  -0.1909,
    '[c](C(c)=N)([cH]c)c(N)c':  -0.1815,
    '[c](-c(c)n)([cH]c)[cH]c':  0.5276,
    '[c](c(O)c)(c(=O)c)c(c)c':  0.2926,
    '[C](c)(=C)O':  -0.6562,
    '[c](Cl)([cH]c)[cH]c':  -0.09005,
    '[c](OC)([cH]c)c(-c)c':  0.6225,
    '[c]1(C)nccn1C':  -0.5236,
    '[c]1(cc)c(C)cnc1c':  0.1867,
    '[CH]([CH]=C)=C(C)C':  -0.3358,
    '[cH]([cH]c)c(Br)c':  -0.2688,
    '[cH]([cH]c)c(N)n':  -0.4731,
    '[cH]([cH]c)c(S)c':  -0.09717,
    '[CH]([CH2]C)([CH2]C)N(C)C':  0.3592,
    '[CH]([CH2]C)(c(c)c)c(c)c':  -0.5984,
    '[CH]([CH2]C)=[CH]C':  -0.1652,
    '[CH]([CH3])([CH2]C)[CH2]C':  -0.1655,
    '[CH]([NH]C)([CH](C)O)C(N)=O':  -0.2754,
    '[cH](c(N)c)c(O)c':  -0.1608,
    '[CH](C)(c)O':  -0.2961,
    '[CH]1(CCCO1)n(c)c':  0.4224,
    '[cH]1ccsc1':  0.2887,
    '[CH2]([CH](C)C)[CH](C)C':  -0.2495,
    '[CH2]([CH](C)N)c(c)c':  0.1386,
    '[CH2]([CH2]c)[CH](C)N':  0.3426,
    '[CH2]([CH2]C)C(C)(C)c':  0.2628,
    '[CH2]([CH2]C)n(c)c':  0.3035,
    '[CH2]([CH2]O)c(c)c':  -0.2682,
    '[CH2](N=C)C(N)=O':  -0.3695,
    '[CH2]1C(C)CCC1O':  0.6055,
    '[CH3]C(=C)N':  0.1034,
    '[CH3]c':  -0.1061,
    '[CH3]N(C)S':  -0.6562,
    '[CH3]n':  -0.1342,
    '[n]([CH3])(c(=O)c)c(=O)n':  -0.3285,
    '[n](c(N)n)c(c)c':  -0.2894,
    '[N](C)(C)N':  -0.3267,
    '[N]1(C)C(C)CCC1C':  0.3069,
    '[NH]([CH](c)C)C(C)=O':  0.5035,
    '[NH](C(N)=O)S(c)(=O)=O':  -0.5549,
    '[NH]=C':  0.2925,
    '[nH]1c(C)ccc1c':  -0.3169,
    '[O]([CH2]C)[CH](C)C':  -0.7976,
    '[O](c)c':  -0.3963,
    '[OH]C(c)=O':  -0.3049,
    '[S](C)(N)(=O)=O':  0.4043,
    'VolSurf_%FU4': -0.0007448,
    'VolSurf_%FU5': -0.0005609,
    'VolSurf_A': -0.01914,
    'VolSurf_AUS7.4': 0.01186,
    'VolSurf_CP': 0.003049,
    'VolSurf_D3': -0.001015,
    'VolSurf_D4': -0.002212,
    'VolSurf_DRACAC': -0.001847,
    'VolSurf_DRDRAC': -0.00227,
    'VolSurf_FLEX': 0.02683,
    'VolSurf_IW4': -0.1014,
    'VolSurf_LgD10': 0.0151,
    'VolSurf_LgD5': -0.03726,
    'VolSurf_LgD8': 0.008123,
    'VolSurf_LgD9': 0.01551,
    'VolSurf_LOGP n-Oct': -0.02285,
    'VolSurf_W7': -0.0008883,
    'VolSurf_WN6': -0.001014,
    'VolSurf_WO1': 0.0005021,
    'VolSurf_WO2': 0.0008443}

RENAL_MODEL_COEFFS = {
    '[cH]([cH]c)c(S)c':  -0.3635,
    '[cH](c)c':  -0.01832,
    '[c](C)(c)c':  -0.01557,
    '[n](C)(c)c':  0.06009,
    '[O]([CH](C)C)[CH2](C)=O':  0.4131,
    '[O]=N':  0.1264,
    '[cH]([cH]c)c(F)c':  -0.5028,
    '[CH2]([CH3])[CH2]C':  -0.2435,
    '[n](C)(c)n':  -0.2675,
    '[CH2]([CH2]N)[NH]C':  0.09271,
    '[N](C)=C':  0.308,
    '[NH](C)S':  0.287,
    '[C]([CH3])(=O)[NH]c':  -0.1525,
    '[c](c(O)c)(c(c)c)n(C)c':  -0.1833,
    '[cH](c(N)c)c(c)n':  0.134,
    '[CH]([NH2])([CH2]C)[CH](O)O':  -0.2143,
    '[c](S)(n)s':  -0.3126,
    '[N]([CH3])([CH3])[CH2]c':  0.3654,
    '[c](=O)(nc)c(C)c':  0.2314,
    '[C]1(=O)CCN1C':  0.2277,
    '[n](c(N)n)c(c)n':  0.1178,
    'VolSurf_%FU5': -0.0009373,
    'VolSurf_%FU6': -0.0007606,
    'VolSurf_%FU7': -0.0007033,
    'VolSurf_%FU8': -0.0003735,
    'VolSurf_A': -0.01034,
    'VolSurf_CD1': 0.2061,
    'VolSurf_CD4': 1.121,
    'VolSurf_CD5': 1.477,
    'VolSurf_CD6': 1.798,
    'VolSurf_CD7': 2.866,
    'VolSurf_CD8': 3.283,
    'VolSurf_CW4': -0.07129,
    'VolSurf_CW5': -0.09199,
    'VolSurf_D2': -0.0006458,
    'VolSurf_D3': -0.001759,
    'VolSurf_D4': -0.00234,
    'VolSurf_D5': -0.003013,
    'VolSurf_DRDRAC': -0.000649,
    'VolSurf_G': -0.0957,
    'VolSurf_HL1': 7.87e-05,
    'VolSurf_HL2': 0.0001944,
    'VolSurf_ID3': 0.2065,
    'VolSurf_ID4': 0.1985,
    'VolSurf_IW2': 0.1432,
    'VolSurf_LgD10': 0.005995,
    'VolSurf_LgD7.5': 0.007247,
    'VolSurf_LgD7': 0.007859,
    'VolSurf_LgD8': 0.007011,
    'VolSurf_LgD9': 0.006609,
    'VolSurf_NCC': -0.03145,
    'VolSurf_PSAR': -0.2119,
    'VolSurf_W4': -4.23e-05,
    'VolSurf_WO1': 7.4e-05,
    'VolSurf_WO2': 0.0001443}

VOLSURF_LOGBB = 'VolSurf_LgBB'

class VolSurfDB():
    def __init__(self):
        self.db = SqliteDict("/DbUCdd/databases/ADME/Clearance/VolSurf3.db",autocommit=True)

    def isMolInDb(self,smiles):
        return self.db.has_key(smiles)

    def getDescriptors(self,smiles):
        try:
            if self.isMolInDb(smiles):
                return self.db[smiles]
            else:
                return None
        except:
            return None

    def insertDescriptor(self,smiles,jsonString):
        self.db[smiles] = jsonString

    def __del__(self):
        print >>sys.stderr, "Volsurf Database closed."
        self.db.close()

class AdmeServer(xmlrpc.XMLRPC):
    def __init__(self):
        xmlrpc.XMLRPC.__init__(self)

    def xmlrpc_Test(self,inputSdf):
        tmpdir = tempfile.mkdtemp(prefix="moka_")
        moka_input = os.path.join(tmpdir,"moka_input.sdf")
        f = open(moka_input,"w")
        sdf = base64.decodestring(inputSdf)
        f.write(sdf)
        f.close()
        return inputSdf

    def xmlrpc_MoKaDescriptor(self,inputSdf):
        result = threads.deferToThread(self.moka_descriptors, inputSdf)
        return result

    def xmlrpc_VolSurfDescriptor(self,inputSdf):
        result = threads.deferToThread(self.volsulf_descriptors, inputSdf)
        return result

    def xmlrpc_Moka(self,inputSdf):
        result = threads.deferToThread(self.moka, inputSdf)
        return result

    def xmlrpc_predictBBB(self,inputSdf):
        result = threads.deferToThread(self.predictBBB, inputSdf)
        return result

    def xmlrpc_RLM(self,inputSdf):
        result = threads.deferToThread(self.RLM,inputSdf)
        return result

    def xmlrpc_inSilicoADME(self,inputSdf):
        result = threads.deferToThread(self.inSilicoADME,inputSdf)
        return result

    def xmlrpc_efflux(self,inputSdf):
        result = threads.deferToThread(self.efflux,inputSdf)
        return result

    def xmlrpc_hERG(self,inputSdf):
        result = threads.deferToThread(self.hERG, inputSdf)
        return result

    def xmlrpc_PPB(self, inputSdf):
        result = threads.deferToThread(self.PPB, inputSdf)
        return result

    def xmlrpc_CLogP(self,smiles):
        try:
            result = subprocess.check_output([CLOGP_COMMAND,smiles]).split()[0]
            return result
        except:
            return ""

    def xmlrpc_CLogPBatch(self,smiles):
        try:
            result = threads.deferToThread(self.clogp_batch, smiles)
            return result
        except:
            return ""


    def xmlrpc_getTareWeight(self,barcode):
        targetweightDb = SqliteDict("/DbUCdd/databases/ADME/TareWeightDB/TareWeight.db",autocommit=True)
        result = "Not Found"
        if targetweightDb.has_key(barcode):
            result = targetweightDb[barcode]
        targetweightDb.close()
        return result


    def xmlrpc_VDssDescriptors(self,inputSdf):
        #Just for test
        result = threads.deferToThread(self.VDssDescriptors(inputSdf))
        print result

    def xmlrpc_predictVdss(self,inputSdf):
        result = threads.deferToThread(self.predictVdss, inputSdf)
        return result

    def xmlrpc_predictAMES(self,inputSdf):
        result = threads.deferToThread(self.predictAMESWrapper,inputSdf)
        return result

    def xmlrpc_check_evotec_inventory(self,amount_uL):
        result = threads.deferToThread(self.check_evotec_inventory, amount_uL)
        return result

    def xmlrpc_efflux_svm(self,inputSdf):
        result = threads.deferToThread(self.predictEffluxRatioSVM,inputSdf)
        return result

    def xmlrpc_freeform(self, inputSdf):
        return threads.deferToThread(self.freeform_hpc,inputSdf)

    def freeform_hpc(self,inputSdf):
        import multiprocessing
        manager = multiprocessing.Manager()
        inputSdf = base64.decodestring(inputSdf)
        n_processors = 4
        trunks = MolUtilities().splitMols(inputSdf, n_processors)
        return_dict = manager.dict()
        jobs = []
        for sdf in trunks:
            p = multiprocessing.Process(target=self.freeform, args=(sdf,return_dict))
            jobs.append(p)
            p.start()
            print "%s started"%p.name

        for p in jobs:
            p.join()

        dict = {}
        for key in return_dict.keys():
            dict[key] = return_dict[key]
        return json.dumps(dict)
# >  <cHLM>  (1)
# 14.493323966
#
# >  <cHPPB>  (1)
# 30.2146568756
#
# >  <cMDR1>  (1)
# 1.25192698312
#
# >  <cRLM>  (1)
# 68.4851902673
#
# >  <cRPPB>  (1)
# 33.1734774995
    def inSilicoADME(self,inputSdf):
        inputSdf = base64.decodestring(inputSdf)
        properties = ["cHLM","cHPPB","cMDR1","cRLM","cRPPB","cSolubility"]
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(inputSdf)
        mol = OEGraphMol()
        ofs = oemolostream()
        ofs.SetFormat(OEFormat_SDF)
        ofs.openstring()
        resultDict = {}

        while OEReadMolecule(ifs,mol):
            moka_id = OEGetSDData(mol, "moka_id")
            mol.SetTitle(moka_id)
            OEWriteMolecule(ofs,mol)
        new_sdf = ofs.GetString()
        ifs.close()
        ofs.close()

        if len(new_sdf)>0:
            tmpdir = tempfile.mkdtemp(prefix="inSilicoADME_")

            sdf_input = os.path.join(tmpdir,"input.sdf")
            sdf_output = os.path.join(tmpdir,"input_cADME.sdf")
            f = open(sdf_input,"w")
            f.write(new_sdf)
            f.close()
            command = INSILICO_ADME_COMMAND%(INSILICO_ADME,sdf_input)

            p = subprocess.Popen(command.split(),cwd=tmpdir)
            p.communicate()
            ifs = oemolistream()
            ifs.open(sdf_output)
            mol = OEGraphMol()
            while OEReadMolecule(ifs,mol):
                try:
                    name = mol.GetTitle()
                    compoundDict = {}
                    for property in properties:
                        value = float(OEGetSDData(mol,property))
                        compoundDict[property] = value
                    resultDict[name] = compoundDict
                except:
                    continue
        return json.dumps(resultDict)

    def freeform(self,sdf,return_dict=None):
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(sdf)
        mol = OEGraphMol()
        ofs = oemolostream()
        ofs.SetFormat(OEFormat_SDF)
        ofs.openstring()
        resultDict = {}

        while OEReadMolecule(ifs,mol):
            moka_id = OEGetSDData(mol, "moka_id")
            mol.SetTitle(moka_id)
            OEWriteMolecule(ofs,mol)
        new_sdf = ofs.GetString()
        print new_sdf
        ifs.close()
        ofs.close()

        subResult = {}
        if len(new_sdf)>0:
            tmpdir = tempfile.mkdtemp(prefix="freeform_")
            sdf_input = os.path.join(tmpdir,"input.sdf")
            tracked_free = os.path.join(tmpdir,"freeform.tracked_free.oeb")
            deltaG_tag = "conf_dGSheff"
            strainLocal_tag = "ElocStrainFreeSheff"
            strainGlobal_tag = "GlblStrainSheff"

            tracked_rstr = os.path.join(tmpdir,"freeform.tracked_rstr.oeb")
            prefix = os.path.join(tmpdir,"freeform")
            f = open(sdf_input,"w")
            f.write(new_sdf)
            f.close()
            p = subprocess.Popen([FREEFORM_COMMAND,"-in",sdf_input,"-track",sdf_input,"-prefix",prefix])
            p.communicate()

            ifs = oemolistream()
            ifs.open(tracked_free)
            mol = OEGraphMol()

            while OEReadMolecule(ifs,mol):
                subResult[mol.GetTitle()] = {}
                subResult[mol.GetTitle()][deltaG_tag] = OEGetSDData(mol,deltaG_tag)
            ifs.close()

            ifs.open(tracked_rstr)
            while OEReadMolecule(ifs,mol):
                if subResult.has_key(mol.GetTitle()):
                    subResult[mol.GetTitle()][strainLocal_tag]= OEGetSDData(mol,strainLocal_tag)
                    subResult[mol.GetTitle()][strainGlobal_tag]= OEGetSDData(mol,strainGlobal_tag)
            ifs.close()

            # subResult = parseVolsurfResult(descriptor_file,False)
            # for moka_id in subResult.keys():
            #     smiles = smilesDict[moka_id]
            #     if not volsurfDb.isMolInDb(smiles):
            #         volsurfDb.insertDescriptor(smiles,json.dumps(subResult[moka_id]))
        result = dict(resultDict.items()+subResult.items())

        if return_dict is not None:
            for key in result:
                if not return_dict.has_key(key):
                    return_dict[key] = result[key]
        return json.dumps(result)


    def check_evotec_inventory(self,amount_uL):
        csv = evotec_inventory_check(amount_uL)
        result = {}
        result["csv"] = csv
        result["count"] = len(csv.split("\n"))
        return json.dumps(result)

    def predictAMESWrapper(self,inputSdf):
        inputSdf = base64.decodestring(inputSdf)
        mol_list = []
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(inputSdf)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            mol_list.append(OEGraphMol(mol))
        ifs.close()

        resultStr = self.predictAMES(mol_list,"/DbUCdd/databases/ADME/AMES/AMES_Tree_2016_Jan_15.rfModel")
        return resultStr

    def predictAMES(self, mol_list, modelName):
        dict = {}
        columnNames = []
        for line_no,mol in enumerate(mol_list):
            fp = OEFingerPrint()
            OEMakeFP(fp,mol,OEFPType_Tree)
            if line_no == 0:
                columnNames = ["%d"%i for i in range(0,fp.GetSize())]
                for colName in columnNames:
                    dict[colName] = []

            for colName in columnNames:
                idx = int(colName)
                desc = 0
                if fp.IsBitOn(idx):
                    desc = 1
                dict[colName].append(desc)
        frame = pd.DataFrame(data=dict,columns=columnNames)
        tmpdir = tempfile.mkdtemp(prefix="randomForest")
        input = os.path.join(tmpdir,"input.csv")
        frame.to_csv(input)

        r = robjects.r
        r.library("randomForest")
        rTraining =  """
                    predictRF<-function(){
                    input <- read.table("%s",sep=",",header=TRUE)
                    load("%s")
                    result<-predict(ames.rf,input)
                    result
                    }
                    """%(input,modelName)
        r(rTraining)
        result = r['predictRF']()
        ames_result = numpy.array(result)-1
        result_dict = {}
        AMES = {0:"Negative",1:"Positive"}
        for id,mol in enumerate(mol_list):
            moka_id = OEGetSDData(mol, "moka_id")
            result_dict[moka_id] = AMES[ames_result[id]]
        return json.dumps(result_dict)

    def predictEffluxRatioSVM(self,inputSdf):
        mdck_train_file = "/DbUCdd/databases/ADME/MDCK/MDCK_SVM_TrainingSet.sdf"
        data_tag = "Ratio(B-A/A-B)1uM"
        volsurfRes = json.loads(self.volsurf_descriptors_hpc(inputSdf))
        mdck_test_file = base64.decodestring(inputSdf)
        #54
        #descriptor_names = ['ChemAxon LogD', 'CLogP', 'ChemAxon Acidic pKa', 'hydrogen-bond acceptors', 'hydrogen-bond donors', 'VolSurf_WO1', 'VolSurf_WO2', 'VolSurf_WO3', 'VolSurf_CW7', 'VolSurf_PSA', 'VolSurf_%FU4', 'VolSurf_%FU5', 'VolSurf_%FU7', 'VolSurf_LgBB', 'VolSurf_LOGP c-Hex', 'VolSurf_DD5', 'VolSurf_DD1', 'VolSurf_CD8', 'VolSurf_IW2', 'VolSurf_CD5', 'VolSurf_DRDRDR', 'VolSurf_W6', 'VolSurf_W3', 'VolSurf_W4', 'VolSurf_W2', 'VolSurf_L4LgS', 'VolSurf_D2', 'VolSurf_ID1', 'VolSurf_DRDODO', 'VolSurf_ID2', 'VolSurf_ID4', 'VolSurf_D7', 'VolSurf_LgS10', 'VolSurf_LgS11', 'VolSurf_LgD7', 'VolSurf_LgD5', 'VolSurf_LgD6', 'VolSurf_V', 'VolSurf_LgS9', 'VolSurf_LgS8', 'VolSurf_ACDODO', 'VolSurf_LgS5', 'VolSurf_LgS7.5', 'VolSurf_CW4', 'VolSurf_CW5', 'VolSurf_CW3', 'VolSurf_NCC', 'VolSurf_A', 'VolSurf_R', 'VolSurf_FLEX_RB', 'VolSurf_LgD10', 'VolSurf_SOLY', 'VolSurf_DODODO', 'VolSurf_WN1']

        #59
        descriptor_names = ['2d PSA', 'ChemAxon LogD', 'CLogP', 'ChemAxon Acidic pKa', 'ChemAxon Basic pKa 1', 'ChemAxon Basic pKa 2','hydrogen-bond acceptors', 'hydrogen-bond donors', 'VolSurf_W7', 'VolSurf_WO1', 'VolSurf_WN6', 'VolSurf_CW8', 'VolSurf_%FU5', 'VolSurf_%FU7', 'VolSurf_%FU10', 'VolSurf_LgBB', 'VolSurf_LOGP c-Hex', 'VolSurf_DD4', 'VolSurf_DD2', 'VolSurf_DD3', 'VolSurf_DD1', 'VolSurf_CD7', 'VolSurf_CD8', 'VolSurf_WO6', 'VolSurf_DD6', 'VolSurf_LOGP n-Oct', 'VolSurf_IW4', 'VolSurf_L4LgS', 'VolSurf_POL', 'VolSurf_D5', 'VolSurf_D4', 'VolSurf_D3', 'VolSurf_D2', 'VolSurf_CP', 'VolSurf_D1', 'VolSurf_ID1', 'VolSurf_DRDODO', 'VolSurf_ID2', 'VolSurf_ID3', 'VolSurf_ID4', 'VolSurf_D8', 'VolSurf_LgS10', 'VolSurf_LgD7', 'VolSurf_SKIN', 'VolSurf_AUS7.4', 'VolSurf_LgS7', 'VolSurf_LgS6', 'VolSurf_LgS3', 'VolSurf_ACACAC', 'VolSurf_PHSAR', 'VolSurf_CW5', 'VolSurf_NCC', 'VolSurf_A', 'VolSurf_R', 'VolSurf_DRACDO', 'VolSurf_FLEX_RB', 'VolSurf_LgD10', 'VolSurf_PSAR', 'VolSurf_DODODO', 'VolSurf_WN3']
        columns = list(descriptor_names)
        print ("%d descriptors are being used."% len(columns))
        columns.insert(0,"Activity")
        dict = {}
        for colName in columns:
            dict[colName] = []
        ifs = oemolistream()
        ifs.open(mdck_train_file)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            if OEHasSDData(mol,data_tag):
                try:
                    efflux_ratio = math.log10(float(OEGetSDData(mol,data_tag)))
                    dict["Activity"].append(efflux_ratio)
                    for name in columns[1:]:
                        if OEHasSDData(mol,name):
                            data = float(OEGetSDData(mol,name))
                        else:
                            data = 0.0
                        dict[name].append(data)
                except:
                    traceback.print_exc(file=sys.stdout)
                    continue
            else:
                continue
        ifs.close()

        train_frame = pd.DataFrame(data=dict,columns=columns)

        dict = {}
        for colName in columns:
            dict[colName] = []
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(mdck_test_file)
        mol = OEGraphMol()
        mol_names = []
        while OEReadMolecule(ifs,mol):
            molId = OEGetSDData(mol,"moka_id")
            if volsurfRes.has_key(molId):
                volRes = volsurfRes[molId]
                mol_names.append(molId)
                dict.pop("Activity",None)
                for name in columns[1:]:
                    if OEHasSDData(mol,name):
                        data = float(OEGetSDData(mol,name))
                    elif volRes.has_key(name):
                        data = float(volRes[name])
                    else:
                        data = 0.0
                    dict[name].append(data)
        ifs.close()
        test_frame = pd.DataFrame(data=dict,columns=columns[1:])

        X_train = train_frame[columns[1:]]
        features = train_frame['Activity'].values

        clf = svm.SVR(kernel="rbf", C=20, gamma=0.01, epsilon=0.1)

        scaler = StandardScaler()
        train_values = scaler.fit_transform(X_train)
        cv = KFold(n_splits=10,shuffle=False)
        a = cross_val_score(clf, train_values, features, cv=cv,scoring="r2",n_jobs=-1)
        clf.fit(train_values,features)

        test_values = scaler.transform(test_frame[columns[1:]])
        predicted = clf.predict(test_values)
        result_dict = {}
        for mol_id, predicted in zip(mol_names,predicted):
            result_dict[mol_id] = round(10**predicted,1)
        return json.dumps(result_dict)
	

    def predictEffluxRatioWrapper(self,inputSdf):
        inputSdf = base64.decodestring(inputSdf)
        mol_list = []
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(inputSdf)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            mol_list.append(OEGraphMol(mol))
        ifs.close()

        resultStr = self.predictEffluxRatio(mol_list,"/DbUCdd/databases/ADME/MDCK/MDCK_2017_Jan_11.svmModel")
        return resultStr

    def predictEffluxRatio(self, mol_list, modelName):
        columnNames =["2d PSA","molecular weight","ChemAxon LogD","CLogP","ChemAxon Acidic pKa","ChemAxon Basic pKa","hydrogen-bond acceptors","hydrogen-bond donors","sum of formal charges"]
        dict = {}
        for line_no,mol in enumerate(mol_list):
            if line_no == 0:
                for colName in columnNames:
                    dict[colName] = []

            for colName in columnNames:
                data = 0
                if OEHasSDData(mol,colName):
                    data = float(OEGetSDData(mol,colName))
                dict[colName].append(data)
        frame = pd.DataFrame(data=dict,columns=columnNames)
        tmpdir = tempfile.mkdtemp(prefix="svm")
        input = os.path.join(tmpdir,"input.csv")
        frame.to_csv(input)

        r = robjects.r
        r.library("e1071")
        rTraining =  """
                    predictSVM<-function(){
                    input <- read.table("%s",sep=",",header=TRUE)
                    load("%s")
                    result<-predict(model.svm,input)
                    result
                    }
                    """%(input,modelName)
        r(rTraining)
        result = r['predictSVM']()
        res = numpy.array(result)-1
        result_dict = {}
        EFFLUX_RATIO = {0:"Low Efflux(<=3)",1:"High Efflux(>3)"}
        for id,mol in enumerate(mol_list):
            moka_id = OEGetSDData(mol, "moka_id")
            result_dict[moka_id] = EFFLUX_RATIO[res[id]]
        return json.dumps(result_dict)


    def predictBBB(self,inputSdf):
        volsurfRes = json.loads(self.volsurf_descriptors_hpc(inputSdf))
        resultDict = {}
        inputSdf = base64.decodestring(inputSdf)
        oemol = OEGraphMol()
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(inputSdf)
        while OEReadMolecule(ifs,oemol):
            molId = OEGetSDData(oemol,"moka_id")
            if volsurfRes.has_key(molId):
                volRes = volsurfRes[molId]
                print volRes.keys()
                resultDict[molId] = volRes[VOLSURF_LOGBB]
        return json.dumps(resultDict)


    def predictVdss(self,inputSdf):
        molecule_list,dict = self.VDssDescriptors(inputSdf)
        pd = pandas.DataFrame(data=dict,index=VDSS_DESCRIPTOR_NAMES,columns=molecule_list).transpose()
        tmpdir = tempfile.mkdtemp(prefix="vdss")
        input = os.path.join(tmpdir,"input.csv")
        pd.to_csv(input)
        r = robjects.r
        r.library("randomForest")
        rTraining =  """
			trainRF<-function(){
			VDss_training <- read.csv("/DbUCdd/databases/ADME/Human_PK_Data/VDss_model/RF_model_July6_2015/External_Test_34/Trainingset_1096cmpds_33descriptors_VDss.csv")
			set.seed(171)
            VDss.rf <- randomForest(logVDss ~ .,data=VDss_training[,2:35],ntree=500,nodesize=10,replace=TRUE,importance=TRUE,promixity=TRUE,nPerm=1,na.action=na.omit)
            VDss_predict <- read.csv("%s")
            prediction <- predict(VDss.rf,VDss_predict)
            cbind(prediction)
			}
			"""%input
        r(rTraining)
        rfObj = r['trainRF']()
        resultDict = {}
        for idx,molId in enumerate(molecule_list):
            resultDict[molId] = 10**(float(rfObj[idx]))
        return json.dumps(resultDict)

    def xmlrpc_predictClearanceMode(self,inputSdf):
        result = threads.deferToThread(self.predictClearanceModeVolSurf,inputSdf)
        return result

    def predictClearanceModeVolSurf(self,inputSdf):
        dict = self.volsurf_clearance_hpc(inputSdf)
        for molId in dict.keys():
            pc_1 = dict[molId]["pc_1"]
            pc_2 = dict[molId]["pc_2"]
            grayMet = 4*pc_1 - 2
            grayRen = 4*pc_1 + 4
            grayRen2 = -2*pc_1
            bileMet = -1.33*pc_1 + 1
            bileRen = 0.6*pc_1 + 0.5
            if pc_2 < grayMet:
                dict[molId]['route'] = "metabolic"
            elif pc_2 > grayRen and pc_2 < grayRen2:
                dict[molId]['route'] = "renal"
            else:
                dict[molId]['route'] = "uncertain"
            if pc_2 < bileRen and pc_2 < bileMet and pc_2 < 0:
                dict[molId]['route'] = "possible billary"
        return json.dumps(dict)

    def predictClearanceMode(self,inputSdf):
        molecule_list,dict = self.calculateClearanceModeDescriptors(inputSdf)
        pd = pandas.DataFrame(data=dict,index=CLEARANCE_MECHANISM_DESCRIPOR_NAMES[1:],columns=molecule_list).transpose()
        tmpdir = tempfile.mkdtemp(prefix="clearance")
        input = os.path.join(tmpdir,"input.csv")
        pd.to_csv(input)
        r = robjects.r
        rTraining =  """
                    predictPCA<-function(){
                    pca_training <- read.csv("/DbUCdd/databases/ADME/Clearance/mechanism_prediction/clearance_volsurf_training.txt")
                    pca_model<-prcomp(pca_training,center=TRUE,scale=TRUE)
                    pca_predict <- read.csv("%s")
                    prediction <- predict(pca_model,pca_predict)
                    prediction[,1:2]
                    }
                    """%input
        r(rTraining)
        rResult = r['predictPCA']()
        rArray = numpy.array(rResult)
        resultDict = {}
        for idx,molId in enumerate(molecule_list):
            resultDict[molId] = {}
            if rArray.ndim == 1:
                pc_1 = rArray[0]
                pc_2 = rArray[1]
            else:
                pc_1 = rArray[idx][0]
                pc_2 = rArray[idx][1]

            is_metabolic = pc_2<(3*pc_1)-1
            if is_metabolic:
                resultDict[molId]["route"]="metabolic_clearance"
            else:
                resultDict[molId]["route"] ="renal_clearance"

            resultDict[molId]["pc_1"]= pc_1
            resultDict[molId]["pc_2"]= pc_2
        return json.dumps(resultDict)

    def xmlrpc_predictMetabolicClearance(self,inputSdf):
        result = threads.deferToThread(self.predictClearance,inputSdf,True)
        return result

    def xmlrpc_predictRenalClearance(self,inputSdf):
        result = threads.deferToThread(self.predictClearance,inputSdf,False)
        return result

    def predictClearance(self,inputSdf, isMetabolic,volsurfRes=None):
        if volsurfRes is None:
            volsurfRes = json.loads(self.volsurf_descriptors_hpc(inputSdf))

        resultDict = {}
        inputSdf = base64.decodestring(inputSdf)
        oemol = OEGraphMol()
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(inputSdf)

        smartsList = CLEARANCE_RENAL_SMARTS
        volsurfList = CLEARANCE_RENAL_VOLSURF
        constant = 0.7735 #renal excreted
        coeff_dict = RENAL_MODEL_COEFFS
        if isMetabolic:
            constant = 1.199 #metabolic cleared
            smartsList = CLEARANCE_META_SMARTS
            volsurfList = CLEARANCE_META_VOLSURF
            coeff_dict = META_MODEL_COEFFS

        while OEReadMolecule(ifs,oemol):
            sum = 0.0
            sum += constant
            molId = OEGetSDData(oemol,"moka_id")
            if volsurfRes.has_key(molId):
                volRes = volsurfRes[molId]
                for volsurf_descname in volsurfList:
                    volsurf_coeff = coeff_dict[volsurf_descname]
                    sum += volsurf_coeff * volRes[volsurf_descname]
                    print volsurf_descname, volsurf_coeff, volsurf_coeff * volRes[volsurf_descname], sum

                for smarts in smartsList:
                    subsearch = OESubSearch(smarts)
                    smarts_coeff = coeff_dict[smarts]
                    matchCount = 0
                    OEPrepareSearch(oemol,subsearch)
                    for dummy in subsearch.Match(oemol, True):
                        matchCount += 1
                    sum += matchCount*smarts_coeff
                    if matchCount >0:
                        print smarts,"%d"%matchCount, smarts_coeff, sum
            resultDict[molId] = 10.0**sum
        return json.dumps(resultDict)

    def calculateClearanceModeDescriptors(self, inputSdf, volsurfRes=None):
        if volsurfRes is None:
            volsurfRes = json.loads(self.volsurf_descriptors_hpc(inputSdf))
        inputSdf = base64.decodestring(inputSdf)
        oemol = OEGraphMol()
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(inputSdf)
        molecule_list = []
        descDict = {}

        while OEReadMolecule(ifs,oemol):
            molId = OEGetSDData(oemol,"moka_id")
            mol_descriptor = []
            if volsurfRes.has_key(molId):
                volRes = volsurfRes[molId]
                for descriptor in CLEARANCE_MECHANISM_DESCRIPOR_NAMES[1:]:
                    mol_descriptor.append(volRes[descriptor])
            if len(mol_descriptor) == len(CLEARANCE_MECHANISM_DESCRIPOR_NAMES)-1:
                molecule_list.append(molId)
                descDict[molId] = mol_descriptor
        ifs.close()
        return molecule_list,descDict

    def VDssDescriptors(self,inputSdf, volsurfRes=None):
        if volsurfRes is None:
            volsurfRes = json.loads(self.volsurf_descriptors_hpc(inputSdf))
        mokaRes = json.loads(self.moka_descriptors(inputSdf))
        blabberRes = json.loads(self.moka_species_abundance(inputSdf))
        inputSdf = base64.decodestring(inputSdf)
        oemol = OEGraphMol()
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(inputSdf)
        molecule_list = []
        descDict = {}

        while OEReadMolecule(ifs,oemol):
            molId = OEGetSDData(oemol,"moka_id")
            mol_descriptor = []
            if volsurfRes.has_key(molId) and mokaRes.has_key(molId) and blabberRes.has_key(molId):
                molRes = mokaRes[molId]
                species_abundance = float(blabberRes[molId]['ABUNDANCE'])/100.0
                volRes = volsurfRes[molId]
                total_anion_charge = 0
                total_cation_charge = 0
                if molRes.has_key("pKa"):
                    all_pKa = molRes["pKa"]
                    for pKa in all_pKa:
                        if pKa['type'] == 'b':
                            total_cation_charge += 1.0/(1.0+10**(7.4-pKa['value']))
                        if pKa['type'] == 'a':
                            total_anion_charge += -1.0/(1.0+10**(pKa['value']-7.4))
                mol_descriptor.append(molId)
                mol_descriptor.append(total_anion_charge)
                mol_descriptor.append(total_cation_charge)
                mol_descriptor.append(species_abundance)
                for descriptor in VDSS_DESCRIPTOR_NAMES[4:14]:
                    mol_descriptor.append(molRes[descriptor])
                for descriptor in VDSS_DESCRIPTOR_NAMES[14:33]:
                    mol_descriptor.append(volRes[descriptor])
                hasSulfur = False
                for atom in oemol.GetAtoms():
                    if atom.GetAtomicNum()== 16:
                        hasSulfur = True
                        break
                if hasSulfur:
                    mol_descriptor.append("TRUE")
                else:
                    mol_descriptor.append("FALSE")
            if len(mol_descriptor) == len(VDSS_DESCRIPTOR_NAMES):
                molecule_list.append(molId)
                descDict[molId] = mol_descriptor
        ifs.close()
        return molecule_list,descDict


    def moka_species_abundance(self, inputSdf):
        tmpdir = tempfile.mkdtemp(prefix="blabber")
        moka_input = os.path.join(tmpdir,"moka_input.sdf")
        moka_output = os.path.join(tmpdir,"moka_output.sdf")
        f = open(moka_input,"w")
        sdf = base64.decodestring(inputSdf)
        f.write(sdf)
        f.close()
        p = subprocess.Popen([BLABBER_COMMAND,moka_input,"-o",moka_output])
        p.communicate()
        jsonResult = parseBlabberResult(moka_output)
        return jsonResult

    def clogp_batch(self,inputSmi):
        tmpdir = tempfile.mkdtemp(prefix="clogp_")
        clogp_input = os.path.join(tmpdir,"clogp.smi")
        outputfile = os.path.join(tmpdir, "clogp_out.txt")
        clogp_output = open(outputfile, "w")
        f = open(clogp_input,"w")
        f.write(inputSmi)
        f.close()
        subprocess.call([CLOGP_COMMAND,"-f",clogp_input],stdout=clogp_output)
        jsonResult = parseCLogPResult(outputfile)
        return jsonResult

    def moka_descriptors(self,inputSdf):
        tmpdir = tempfile.mkdtemp(prefix="moka_desc_")
        moka_input = os.path.join(tmpdir,"moka_input.sdf")
        moka_output = os.path.join(tmpdir,"moka_output.sdf")
        f = open(moka_input,"w")
        sdf = base64.decodestring(inputSdf)
        f.write(sdf)
        f.close()
        p = subprocess.Popen([MOKA_COMMAND,"--load-model=/DbUCDD/databases/ADME/MoKa_Current.mkd","-s","moka_id",
                              "--show-logp","--show-logd=7.5-13","-o",moka_output,moka_input])
        p.communicate()
        jsonResult = parseMokaDescResult(moka_output)
        return jsonResult

    def volsurf_descriptors_hpc(self,inputSdf):
        import multiprocessing
        manager = multiprocessing.Manager()
        inputSdf = base64.decodestring(inputSdf)
        n_processors = 8
        trunks = MolUtilities().splitMols(inputSdf, n_processors)
        return_dict = manager.dict()
        jobs = []
        for sdf in trunks:
            p = multiprocessing.Process(target=self.volsulf_descriptors, args=(sdf,return_dict))
            jobs.append(p)
            p.start()
            print "%s started"%p.name

        for p in jobs:
            p.join()

        dict = {}
        for key in return_dict.keys():
            dict[key] = return_dict[key]
        return json.dumps(dict)

    def volsulf_descriptors(self,sdf,return_dict=None):
        volsurfDb = VolSurfDB()
        ifs = oemolistream()
        ifs.SetFormat(OEFormat_SDF)
        ifs.openstring(sdf)
        mol = OEGraphMol()
        ofs = oemolostream()
        ofs.SetFormat(OEFormat_SDF)
        ofs.openstring()
        resultDict = {}
        smilesDict = {}
        while OEReadMolecule(ifs,mol):
            smiles = OEMolToSmiles(mol)
            moka_id = OEGetSDData(mol,"moka_id")
            smilesDict[moka_id] = smiles
            if volsurfDb.isMolInDb(smiles):
                print "descriptor found: "+smiles
                descriptor = json.loads(volsurfDb.getDescriptors(smiles))
                resultDict[moka_id] = descriptor
            else:
                OEWriteMolecule(ofs,mol)

        new_sdf = ofs.GetString()
        ifs.close()
        ofs.close()
        subResult = {}
        if len(new_sdf)>0:
            tmpdir = tempfile.mkdtemp(prefix="volsurf_")
            sdf_input = os.path.join(tmpdir,"volsurf_input.sdf")
            f = open(sdf_input,"w")
            f.write(new_sdf)
            f.close()
            cmd = os.path.join(tmpdir,"volsurf.vsl")
            descriptor_file = os.path.join(tmpdir,"volsurf_descriptors.csv")
            f = open(cmd,"w")
            f.write(command_file%(sdf_input,descriptor_file))
            f.close()
            p = subprocess.Popen([VOLSURF_COMMAND,cmd])
            p.communicate()
            subResult = parseVolsurfResult(descriptor_file,False)
            for moka_id in subResult.keys():
                smiles = smilesDict[moka_id]
                if not volsurfDb.isMolInDb(smiles):
                    volsurfDb.insertDescriptor(smiles,json.dumps(subResult[moka_id]))
        result = dict(resultDict.items()+subResult.items())

        if return_dict is not None:
            for key in result:
                if not return_dict.has_key(key):
                    return_dict[key] = result[key]
        return json.dumps(result)

    def volsurf_clearance_hpc(self,inputSdf):
        import multiprocessing
        manager = multiprocessing.Manager()
        inputSdf = base64.decodestring(inputSdf)
        n_processors = 8
        trunks = MolUtilities().splitMols(inputSdf, n_processors)
        return_dict = manager.dict()
        jobs = []
        for sdf in trunks:
            p = multiprocessing.Process(target=self.volsulf_clearance, args=(sdf,return_dict))
            jobs.append(p)
            p.start()
            print "%s started"%p.name

        for p in jobs:
            p.join()

        dict = {}
        for key in return_dict.keys():
            dict[key] = return_dict[key]
        return dict

    def volsulf_clearance(self,sdf,return_dict=None):
        subResult = {}
        if len(sdf)>0:
            tmpdir = tempfile.mkdtemp(prefix="volsurf_")
            sdf_input = os.path.join(tmpdir,"volsurf_input.sdf")
            f = open(sdf_input,"w")
            f.write(sdf)
            f.close()
            cmd = os.path.join(tmpdir,"volsurf.vsl")
            descriptor_file = os.path.join(tmpdir,"volsurf_clearance.csv")
            f = open(cmd,"w")
            f.write(clearance_command_file%(sdf_input,descriptor_file))
            f.close()
            p = subprocess.Popen([VOLSURF_COMMAND,cmd])
            p.communicate()
            pc1_dict,pc2_dict = parseCleranceResult(descriptor_file)
            for moka_id in pc1_dict.keys():
                subResult[moka_id] = {}
                subResult[moka_id]["pc_1"] = pc1_dict[moka_id]
                subResult[moka_id]["pc_2"] = pc2_dict[moka_id]

        if return_dict is not None:
            for key in subResult:
                if not return_dict.has_key(key):
                    return_dict[key] = subResult[key]
        return subResult


    def moka(self, inputSdf):
        tmpdir = tempfile.mkdtemp(prefix="moka_")
        moka_input = os.path.join(tmpdir,"moka_input.sdf")
        moka_output = os.path.join(tmpdir,"moka_output.sdf")
        f = open(moka_input,"w")
        sdf = base64.decodestring(inputSdf)
        f.write(sdf)
        f.close()
        p = subprocess.Popen([MOKA_COMMAND,"--load-model=/DbUCDD/databases/ADME/MoKa_Current.mkd","-s","moka_id",
                              "--show-logp","--show-logd=7.4","-o",moka_output,moka_input])
        p.communicate()
        jsonResult = parseMokaResult(moka_output)
        return jsonResult

    # krige_mol.py -train /UserUCdd/ienyedy/OE/ADME/TRAINING_SETS/hERG_IC50_combined_solv.sdf
    # -in probe2.sdf -response_tag "HERG_IC50_uM" -krige_log
    # -classification_boundries 1 3 10  -report_significant_figures 1
    # -title "hERG inhibition"  -unmeasured_values 10.0 30.0 50.0 100.0
    # -local_krige 0  -response_name "IC50(uM)"  -prefix  tmp

    def hERG(self,inputSdf):
        tmpdir = tempfile.mkdtemp(prefix="hERG_")
        raw = os.path.join(tmpdir,"raw.sdf")
        input = os.path.join(tmpdir,"input.sdf")
        tag = "IC50(uM)"
        f = open(raw,"w")
        sdf = base64.decodestring(inputSdf)
        f.write(sdf)
        f.close()
        prefix = "tmp"
        output = os.path.join(tmpdir,"%s_molecules.sdf"%prefix)

        p = subprocess.Popen([FIXPKA_COMMAND,"-in",raw,"-out",input])
        p.communicate()

        p = subprocess.Popen([KRIG_COMMAND,"-train",
                              "/DbUCdd/databases/ADME/Krige_Training_Sets/hERG_IC50_combined.sdf",
                              "-in",input,
                              "-response_tag","HERG_IC50_uM",
                              "-krige_log",
                              "-classification_boundries","1","3","10",
                              "-report_significant_figures","1",
                              "-title","hERG inhibition",
                              "-unmeasured_values","10.0", "30.0", "50.0", "100.0",
                              "-local_krige","0",
                              "-response_name","%s"%tag,
                              "-prefix",os.path.join(tmpdir,prefix)
                              ])
        p.communicate()
        return parse_Krig_result(output,tag)

## Predict efflux
# krige_mol.py -train /UserUCdd/ienyedy/OE/ADME/TRAINING_SETS/PgP/MDCK-MDR1_1uM_solv_hits4M2T_A.sdf \
# -in /UserUCdd/ienyedy/NEW_TARGETS/PLD1/KRIGING/RESULTS/EAP01032017_solv_hERG_molecules.sdf \
#      -response_tag "Ratio(B-A/A-B)1uM" \
#                    -universal_tag "Electrostatic_Solv_Energy" \
#                                   -universal_tag "Molecular_Weight" \
#                                                  -krige_log \
#                                                  -classification_boundries 3 10 \
#                                                                              -report_significant_figures 2 \
#                                                                                                          -report ./RESULTS/EAP01032017_solv_Efflux.pdf \
#                                                                                                          -title "Efflux Ratio" \
#                                                                                                                 -local_krige 0 \
#                                                                                                                              -response_name "Ratio(B-A/A-B)1uM" \
#                                                                                                                                             -prefix ./RESULTS/EAP01032017_solv_Efflux

    def efflux(self,inputSdf):
        tmpdir = tempfile.mkdtemp(prefix="efflux_")
        raw = os.path.join(tmpdir,"raw.sdf")
        input = os.path.join(tmpdir,"input.sdf")
        tag = "Ratio(B-A/A-B)1uM"
        f = open(raw,"w")
        sdf = base64.decodestring(inputSdf)
        f.write(sdf)
        f.close()
        prefix = "tmp"
        output = os.path.join(tmpdir,"%s_molecules.sdf"%prefix)

        p = subprocess.Popen([FIXPKA_COMMAND,"-in",raw,"-out",input])
        p.communicate()

        p = subprocess.Popen([KRIG_COMMAND,"-train",
                              "/DbUCdd/databases/ADME/Krige_Training_Sets/MDCK-MDR1_1uM.sdf",
                              "-in",input,
                              "-response_tag","Ratio(B-A/A-B)1uM_mean",
                              "-krige_log",
                              "-classification_boundries","3","10",
                              "-report_significant_figures","2",
                              "-title","Efflux Ratio",
                              "-local_krige","0",
                              "-response_name","%s"%tag,
                              "-prefix",os.path.join(tmpdir,prefix)
                              ])
        p.communicate()
        return parse_Krig_result(output,tag)

# RLM turnover

# krige_mol.py -train /UserUCdd/ienyedy/OE/ADME/TRAINING_SETS/Microsomal_Stability_Qh_Rat_solv.sdf \
# -in /UserUCdd/ienyedy/NEW_TARGETS/PLD1/KRIGING/RESULTS/EAP01032017_solv_Efflux_molecules.sdf \
#      -response_tag "Rat_Qh" \
#                    -krige_log \
#                    -classification_boundries 10 75 \
#                                                 -report_significant_figures 1 \
#                                                                             -report ./RESULTS/EAP01032017_solv_RLM.pdf \
#                                                                             -title "RLM Stability" \
#                                                                                    -unmeasured_values 10.0 30.0 50.0 100.0 \
#                                                                                                                      -local_krige 0 \
#                                                                                                                                   -response_name "%Qh" \
#                                                                                                                                                  -prefix ./RESULTS/EAP01032017_solv_hits_RLM
    def RLM(self,inputSdf):
        tmpdir = tempfile.mkdtemp(prefix="RLM_")
        raw = os.path.join(tmpdir,"raw.sdf")
        input = os.path.join(tmpdir,"input.sdf")
        tag = "%Qh"
        f = open(raw,"w")
        sdf = base64.decodestring(inputSdf)
        f.write(sdf)
        f.close()
        prefix = "tmp"
        output = os.path.join(tmpdir,"%s_molecules.sdf"%prefix)

        p = subprocess.Popen([FIXPKA_COMMAND,"-in",raw,"-out",input])
        p.communicate()

        p = subprocess.Popen([KRIG_COMMAND,"-train",
                              "/DbUCdd/databases/ADME/Krige_Training_Sets/Microsomal_Stability_Qh_Rat.sdf",
                              "-in",input,
                              "-response_tag","Qh_Mean",
                              "-krige_log",
                              "-classification_boundries","10","75",
                              "-report_significant_figures","1",
                              "-title","RLM Stability",
                              "-local_krige","0",
                              "-unmeasured_values", "10.0", "30.0", "50.0", "100.0",
                              "-response_name","%s"%tag,
                              "-prefix",os.path.join(tmpdir,prefix)
                              ])
        p.communicate()
        return parse_Krig_result(output,tag)


    def PPB(self,inputSdf):
        tmpdir = tempfile.mkdtemp(prefix="PPB_")
        input = os.path.join(tmpdir,"input.sdf")
        output = os.path.join(tmpdir,"output.sdf")
        f = open(input,"w")
        sdf = base64.decodestring(inputSdf)
        f.write(sdf)
        f.close()
        p = subprocess.Popen([PPB_COMMAND,"-in",input,"-out",output,
                              "-predicted_tag","Pred_Fu","-training_molecules",
                              "/DbUCdd/databases/ADME/Krige_ADME/Plasma_Protein_Binding/PPB_dataset_2098.sdf",
                              "-response_tag","PPB_HUMAN_unbound","-take_log_of_response","true"])
        p.communicate()
        return parseResultSimple(output,"Pred_Fu", True)



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

def parse_Krig_result(filename,tag):
    data_tag = "Krige(%s)"%tag
    lower_limit_tag = "Krige(%s) 95%% confidence lower"%tag
    upper_limit_tag = "Krige(%s) 95%% confidence upper"%tag
    ifs = oemolistream()
    ifs.open(filename)
    mol = OEGraphMol()
    result = {}
    while OEReadMolecule(ifs,mol):
        try:
            data = float(OEGetSDData(mol, data_tag))
            lower_data = OEGetSDData(mol, lower_limit_tag)
            lower = float(lower_data)
            upper = float(OEGetSDData(mol, upper_limit_tag))
            mol_id = OEGetSDData(mol,"moka_id")
            result[mol_id] = [data,lower,upper]
        except:
            continue
    ifs.close()
    return json.dumps(result)


def parseResultSimple(filename,tag, useLog):
    ifs = oemolistream()
    ifs.open(filename)
    mol = OEGraphMol()
    result = {}
    while OEReadMolecule(ifs,mol):
        data = OEGetSDData(mol, tag)
        mol_id = OEGetSDData(mol,"moka_id")
        if useLog:
            try:
                result[mol_id] = 10.0**float(data)
            except:
                pass
        else:
            result[mol_id] = data
    ifs.close()
    return json.dumps(result)

def parseBlabberResult(filename):
    ifs = oemolistream()
    ifs.open(filename)
    mol = OEGraphMol()
    result = {}
    while OEReadMolecule(ifs,mol):
        dict = {}
        data = OEGetSDData(mol, "ABUNDANCE")
        p = parse("{percent}% at pH 7.4",data)
        if p is not None:
            dict["ABUNDANCE"] = float(p['percent'])
        else:
            dict["ABUNDANCE"] = 0.0
        mol_id = OEGetSDData(mol,"moka_id")
        result[mol_id] = dict
    ifs.close()
    return json.dumps(result)


def parseMokaResult(filename):
    ifs = oemolistream()
    ifs.open(filename)
    mol = OEGraphMol()
    result = {}
    while OEReadMolecule(ifs,mol):
        mol_result = {}
        pka_result = []
        pka_data = OEGetSDData(mol, "MoKa")
        if len(pka_data) > 0:
            p = parse("{molName} {covalent_hydration} - {numPka} {result}", pka_data)
            if p is not None:
                numPka = int(p.named["numPka"])
                r = p.named["result"]
                arry = r.strip().split(" ")
                for i in range(0,numPka):
                    n = i*4
                    pka = pKa(arry[n],arry[n+1],arry[n+2],arry[n+3])
                    pka_result.append(pka.getDict())
                mol_result["pKa"] = pka_result
        logd_data = OEGetSDData(mol,"MoKa.LogD")
        if len(logd_data) > 0:
            p = parse("7.4: {result}",logd_data)
            mol_result["MoKa_LogD7.4"] = p.named["result"]
        logp_data = OEGetSDData(mol,"MoKa.LogP")
        if len(logp_data) > 0:
            mol_result["MoKa_LogP"] = logp_data
        mol_id = OEGetSDData(mol,"moka_id")
        result[mol_id] = mol_result
    ifs.close()
    return json.dumps(result)

def parseCLogPResult(filename):
    dict = {}
    f = open(filename,"r")
    for line in f.read().splitlines():
        args = line.split()
        dict[args[3]] = args[0]
    return json.dumps(dict)

def parseMokaDescResult(filename):
    ifs = oemolistream()
    ifs.open(filename)
    mol = OEGraphMol()
    result = {}
    while OEReadMolecule(ifs,mol):
        mol_result = {}
        pka_result = []
        pka_data = OEGetSDData(mol, "MoKa")
        if len(pka_data) > 0:
            p = parse("{molName} {covalent_hydration} - {numPka} {result}", pka_data)
            if p is not None:
                numPka = int(p.named["numPka"])
                r = p.named["result"]
                arry = r.strip().split(" ")
                for i in range(0,numPka):
                    n = i*4
                    pka = pKa(arry[n],arry[n+1],arry[n+2],arry[n+3])
                    pka_result.append(pka.getDict())
                mol_result["pKa"] = pka_result
        logd_data = OEGetSDData(mol,"MoKa.LogD")
        if len(logd_data) > 0:
            lines = logd_data.split("\n")
            for line in lines:
                p = parse("{ph}: {result}",line)
                ph = p.named['ph']
                mol_result["LogD%s"%ph] = p.named["result"]
        logp_data = OEGetSDData(mol,"MoKa.LogP")
        if len(logp_data) > 0:
            mol_result["MoKa_LogP"] = logp_data
        mol_id = OEGetSDData(mol,"moka_id")
        result[mol_id] = mol_result
    ifs.close()
    return json.dumps(result)

def parseVolsurfResult(filename,convert=True):
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    result = {}
    tags = []
    for line_no,line in enumerate(lines):
        if line_no == 0:
            tags = line.strip().split(",")[1:]
        else:
            data = line.strip().split(",")
            molid = data[0]
            result[molid] = {}
            for i,d in enumerate(data[1:]):
                result[molid]["VolSurf_"+tags[i]] = float(d)
    if convert:
        return json.dumps(result)
    else:
        return result

def parseCleranceResult(filename):
    dataframe = pandas.read_csv(filename)
    dict = dataframe.to_dict()
    id_dict = dict['Unnamed: 0']
    pc1_dict = {}
    for key in dict['PC 1']:
        pc1_dict[id_dict[key]] = dict['PC 1'][key]
    pc2_dict = {}
    for key in dict['PC 2']:
        pc2_dict[id_dict[key]] = dict['PC 2'][key]
    return pc1_dict,pc2_dict


if __name__ == "__main__":
    if __name__ == "__main__":
        from twisted.internet import reactor
        # sdf = open("/Users/jfeng1/tmp_aligned.sdf","r").read()
        admeServer = AdmeServer()
        # print admeServer.freeform_hpc(base64.encodestring(sdf))
        #for production
        #reactor.listenTCP(9528,server.Site(admeServer))
        #for Dev
        reactor.listenTCP(9555,server.Site(admeServer))
        reactor.run()
