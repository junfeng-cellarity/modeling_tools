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
    model_names = ["msm","hERG","mdck"]
    rsq = []
    rsq_file = open("rsq.txt","w")
    for i in range(0,2):
        for model_name in model_names:
            log_transformer = dc.trans.LogTransformer(transform_y=True)
            metric = dc.metrics.Metric(dc.metrics.pearson_r2_score, np.mean)
            if model_name == "hERG":
                hERG_sdf = "/home/jfeng/apache-tomcat-9.0.13/webapps/data/hERG_project.sdf"
                hERG_model_dir = "/home/jfeng/apache-tomcat-9.0.13/webapps/data/adme/hERG_model"
                data_tags = ["hERG IC50 Pharmaron"]
                hERG_model = ADMEModel(log_transformer,hERG_model_dir)
                hERG_model.build_dataset(hERG_sdf,data_tags)
                hERG_model.train_test_model(150)
                rsq.append(hERG_model.rsq)
                print(hERG_model.rsq)
            elif model_name == "mdck":
                mdck_sdf = "/home/jfeng/apache-tomcat-9.0.13/webapps/data/mdck_project.sdf"
                mdck_model_dir = "/home/jfeng/apache-tomcat-9.0.13/webapps/data/adme/mdck_model"
                mdck_data_tags = ["Papp A-B"]
                mdck_model = ADMEModel(log_transformer,mdck_model_dir)
                mdck_model.build_dataset(mdck_sdf,mdck_data_tags)
                mdck_model.train_test_model(50)
                rsq.append(mdck_model.rsq)
                print(mdck_model.rsq)
            elif model_name == "msm":
                msm_sdf = "/home/jfeng/apache-tomcat-9.0.13/webapps/data/msm_project.sdf"
                msm_model_dir = "/home/jfeng/apache-tomcat-9.0.13/webapps/data/adme/msm_model"
                msm_data_tags = ["General Metabolic Stability: T1/2 (min)"]
                msm_model = ADMEModel(log_transformer,msm_model_dir)
                msm_model.build_dataset(msm_sdf,msm_data_tags)
                msm_model.train_test_model(50)
                rsq.append(msm_model.rsq)
                print(msm_model.rsq)

    print(rsq,file=rsq_file)
