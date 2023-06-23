#!/usr/bin/env python
import requests
import sys
import pandas as pd
import sqlalchemy
from StringIO import StringIO
from pandas.io import sql


def batchId2CorporateId(row):
    return row['Batch Id'][0:11]

def batchId2BioNumber(batchId):
    if batchId is not None:
        return batchId[0:11]
    else:
        return None

def evotec_inventory_check(bionumber_list):
    #login_url = "https://cmct.evotec.com/EvotecCM/loginAction.action"
    #session = requests.Session()
    #response = session.post(login_url,data={"username":"jun.feng@biogen.com","password":"Beard301","buttonName":"Log In"})
    #download_url = "https://cmct.evotec.com/EvotecCM/fullInventoryDownload.action"
    #response = session.get(download_url,stream=True)
    #csv_string = response.content
    #csv_string = open("/DbUCdd/databases/Evotec/evotec_inventory_current.csv","r").read()
    corporateIds = pd.DataFrame(data=bionumber_list,columns=['Corporate Id'])['Corporate Id']
    #corporateIds = df[df['MaxAmount']<20]['Corporate ID']
    print corporateIds
    csv_string = open("/Users/jfeng1/BiogenDB/Evotec/evotec_inventory_current.csv","r").read()
    dataframe = pd.read_csv(StringIO(csv_string))

    dataframe['Corporate ID'] = dataframe.apply(lambda row:batchId2CorporateId(row),axis = 1)

    engine = sqlalchemy.create_engine("sqlite://")
    sql.execute("drop table if exists data", engine)
    dataframe.to_sql("data",engine)

    df = pd.read_sql_query("select *, Max(`Available Amount`) as MaxAmount from data where "
                            "`Conc. Unit`='mM' and `Labware Type`='BGN_Vial_MLR96_(1.4mL)' "
                            "and `Available Amount Unit`='uL' and `Concentration`=10 group by `Corporate Id` order by MaxAmount",engine)

    wanted = df[df['Corporate ID'].isin(corporateIds)]
    return wanted.to_csv()


if __name__ == "__main__":
    if len(sys.argv)!=2:
        print ("Usage:%s bionumber_list.txt"%sys.argv[0])
    else:
        bionumbers = open(sys.argv[1],"r").readlines()
        bionumber_list = []
        for lot_id in bionumbers:
            bionumber_list.append(batchId2BioNumber(lot_id).strip())
        print evotec_inventory_check(bionumber_list)
