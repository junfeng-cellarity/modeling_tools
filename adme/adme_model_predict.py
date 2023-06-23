#!/usr/bin/env python
###############################################################
#          Must have deepchem and rdkit                       #
#          OpenEye is not used                                #
###############################################################
import os,sys
import deepchem as dc
import tempfile
from rdkit import Chem
import math,random
import csv
from deepchem.models.tensorgraph.models.graph_models import GraphConvModel
import tensorflow as tf
import numpy as np
from adme_model import ADMEModel

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage:%s model(hERG or mdck, or msm) input.sdf output.sdf"%sys.argv[0])
    else:
        if os.path.exists(sys.argv[2]):
            model_name = sys.argv[1]
            print("|"+model_name+"|")
            if model_name != "hERG" and model_name != "mdck" and model_name != "msm":
                print("Unknown model: use hERG or mdck or msm")
                sys.exit(1)

            log_transformer = dc.trans.LogTransformer(transform_y=True)
            mols = [mol for mol in Chem.SDMolSupplier(sys.argv[2]) if mol != None]
            metric = dc.metrics.Metric(dc.metrics.pearson_r2_score, np.mean)
            if model_name == "hERG":
                hERG_sdf = "/home/jfeng/apache-tomcat-9.0.13/webapps/data/hERG_project.sdf"
                hERG_model_dir = "/home/jfeng/apache-tomcat-9.0.13/webapps/data/adme/hERG_model"
                data_tag = "hERG IC50 Pharmaron"
                hERG_model = ADMEModel(log_transformer,hERG_model_dir)
                hERG_model.load()
                hERG_model.predict_mols(mols,data_tag)
            elif model_name == "mdck":
                mdck_sdf = "/home/jfeng/apache-tomcat-9.0.13/webapps/data/mdck_project.sdf"
                mdck_model_dir = "/home/jfeng/apache-tomcat-9.0.13/webapps/data/adme/mdck_model"
                mdck_data_tag = "Papp A-B"
                mdck_model = ADMEModel(log_transformer,mdck_model_dir)
                mdck_model.load()
                mdck_model.predict_mols(mols,mdck_data_tag)
            elif model_name == "msm":
                msm_sdf = "/home/jfeng/apache-tomcat-9.0.13/webapps/data/msm_project.sdf"
                msm_model_dir = "/home/jfeng/apache-tomcat-9.0.13/webapps/data/adme/msm_model"
                msm_data_tag = "General Metabolic Stability: T1/2 (min)"
                msm_model = ADMEModel(log_transformer,msm_model_dir)
                msm_model.load()
                msm_model.predict_mols(mols,msm_data_tag)

            sw_writer = Chem.SDWriter(sys.argv[3])
            for id,mol in enumerate(mols):
                sw_writer.write(mol)
            sw_writer.close()
        else:
            print("Input file does not exist.")
            sys.exit(1)
