__author__ = 'Jun Feng'
import os,subprocess
import tempfile
import sys
import json
from twisted.web import xmlrpc, server
from twisted.internet import threads
from rdkit import Chem
import base64

ADME_COMMAND = "/home/jfeng/Programming/modeling/adme/adme_model_predict.py"
hERG_tag = "hERG IC50 Pharmaron_prediction"
mdck_tag = "Papp A-B_prediction"
msm_tag = "General Metabolic Stability: T1/2 (min)_prediction"

class AdmeServer(xmlrpc.XMLRPC):
    def __init__(self):
        xmlrpc.XMLRPC.__init__(self)

    def xmlrpc_hERG(self,inputSdf):
        result = threads.deferToThread(self.hERG, inputSdf)
        return result

    def xmlrpc_mdck(self,inputSdf):
        result = threads.deferToThread(self.mdck,inputSdf)
        return result

    def xmlrpc_msm(self,inputSdf):
        return threads.deferToThread(self.msm,inputSdf)

    def msm(self, inputSdf):
        tmpdir = tempfile.mkdtemp(prefix="msm_")
        input = os.path.join(tmpdir,"input.sdf")
        f = open(input,"w")
        f.write(inputSdf)
        f.close()
        output = os.path.join(tmpdir,"msm_prediction.sdf")
        p = subprocess.Popen([ADME_COMMAND,"msm",input,output])
        p.communicate()
        return parse_result(output,msm_tag)

    def mdck(self,inputSdf):
        tmpdir = tempfile.mkdtemp(prefix="mdck_")
        input = os.path.join(tmpdir,"input.sdf")
        f = open(input,"w")
        f.write(inputSdf)
        f.close()
        output = os.path.join(tmpdir,"mdck_prediction.sdf")
        p = subprocess.Popen([ADME_COMMAND,"mdck",input,output])
        p.communicate()
        return parse_result(output,mdck_tag)

    def hERG(self,inputSdf):
        tmpdir = tempfile.mkdtemp(prefix="hERG_")
        input = os.path.join(tmpdir,"input.sdf")
        f = open(input,"w")
        f.write(inputSdf)
        f.close()
        output = os.path.join(tmpdir,"hERG_prediction.sdf")
        p = subprocess.Popen([ADME_COMMAND,"hERG",input,output])
        p.communicate()
        return parse_result(output,hERG_tag)

def parse_result(filename,tag):
    result = {}
    sd_reader = Chem.SDMolSupplier(filename)
    for mol in sd_reader:
        mol_dict = mol.GetPropsAsDict()
        if "moka_id" in mol_dict and tag in mol_dict:
            try:
                result[mol_dict["moka_id"]] = float(mol_dict[tag])
            except:
                pass
    return json.dumps(result)

if __name__ == "__main__":
    if __name__ == "__main__":
        from twisted.internet import reactor
        admeServer = AdmeServer()
        reactor.listenTCP(9555,server.Site(admeServer))
        reactor.run()
