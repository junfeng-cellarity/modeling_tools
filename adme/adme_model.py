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
import re

class ADMEModel:
    def __init__(self,y_transformer,model_dir):
        self.y_transformer = y_transformer
        self.model_dir = model_dir
        self.featurizer_func = dc.feat.ConvMolFeaturizer()

    # @classmethod
    # def load_from_file(cls, model_dir, data_tag,  y_transformer):
    #     new_model = cls(None,data_tag,y_transformer)
    #     new_model.model = dc.models.TensorGraph.load_from_dir(model_dir)
    #     new_model.model.model_dir = model_dir
    #     new_model.featurizer_func = dc.feat.ConvMolFeaturizer()
    #     return new_model

    def build_dataset(self,sdf_name,data_tags):
        csvfname = tempfile.mktemp(".csv","adme")
        with open(csvfname,"w") as csvfile:
            datawriter = csv.writer(csvfile)
            datawriter.writerow(["smiles","name","activity"])
            training_mols = [mol for mol in Chem.SDMolSupplier(sdf_name) if mol is not None]
            for mol in training_mols:
                data = []
                for data_tag in data_tags:
                    if mol.HasProp(data_tag) and len(mol.GetProp(data_tag)) > 0:
                        #activity = -math.log10(1e-6*float(mol.GetProp(data_tag)))
                        activity = float(re.sub(">|<|=","",mol.GetProp(data_tag)))
                        data.append(activity)
                        # activity = 0
                        # if raw_act>5:
                        #     activity = 1
                if len(data)>=1:
                    datawriter.writerow([Chem.MolToSmiles(mol),mol.GetProp("_Name"),np.mean(activity)])

        loader = dc.data.CSVLoader(tasks=['activity'],smiles_field="smiles",featurizer=self.featurizer_func)
        self.dataset = loader.featurize(csvfname)
        if self.y_transformer is not None:
            self.dataset = self.y_transformer.transform(self.dataset)

        from deepchem.splits.splitters import IndexSplitter
        splitter=IndexSplitter()
        self.train_data,self.test_data=splitter.train_test_split(self.dataset,frac_train=0.8)

    def train_model(self,n_epoch):
        config = tf.ConfigProto()
        config.gpu_options.allow_growth = False
        self.model = GraphConvModel(1,mode='regression',configproto=config,model_dir=self.model_dir)
        self.model.fit(self.dataset,nb_epoch=n_epoch)

    def train_test_model(self,n_epoch):
        import tensorflow as tf
        config = tf.ConfigProto()
        config.gpu_options.allow_growth = False
        self.model = GraphConvModel(1,mode='regression',configproto=config,model_dir=self.model_dir)
        self.model.fit(self.train_data,nb_epoch=n_epoch)

        import numpy as np
        metric = dc.metrics.Metric( dc.metrics.pearson_r2_score, np.mean )
        self.rsq = self.model.evaluate(self.test_data, [metric])['mean-pearson_r2_score']

    def predict_sdf(self,input_sdf, output_sdf, data_tag):
        mols = [ mol for mol in Chem.SDMolSupplier(input_sdf) if mol != None ]
        X=self.featurizer_func.featurize(mols)
        predicted = self.model.predict_on_batch(X).astype('float')
        if self.y_transformer is not None:
            predicted = self.y_transformer.untransform(predicted)
        sw_writer = Chem.SDWriter(output_sdf)
        for id,mol in enumerate(mols):
            mol.SetProp("%s_prediction"%data_tag,"%5.2f"%predicted[id][0])
            sw_writer.write(mol)
        sw_writer.close()
        return

    def predict_mols(self, mols, data_tag):
        if mols is not None and len(mols)>0:
            X=self.featurizer_func.featurize(mols)
            predicted = self.model.predict_on_batch(X).astype('float')
            if self.y_transformer is not None:
                predicted = self.y_transformer.untransform(predicted)
            for id,mol in enumerate(mols):
                mol.SetProp("%s_prediction"%data_tag,"%5.2f"%predicted[id][0])

    def save(self):
        if self.model is not None:
            self.model.save()

    def load(self):
        self.model = dc.models.TensorGraph.load_from_dir(self.model_dir)

