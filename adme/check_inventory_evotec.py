import requests
import pandas as pd
import sqlalchemy
from StringIO import StringIO
from pandas.io import sql

def batchId2CorporateId(row):
    return row['Batch Id'][0:11]

def evotec_inventory_check(amount_uL):
    #login_url = "https://cmct.evotec.com/EvotecCM/loginAction.action"
    #session = requests.Session()
    #response = session.post(login_url,data={"username":"jun.feng@biogen.com","password":"Beard301","buttonName":"Log In"})
    #download_url = "https://cmct.evotec.com/EvotecCM/fullInventoryDownload.action"
    #response = session.get(download_url,stream=True)
    #csv_string = response.content
    csv_string = open("/DbUCdd/databases/Evotec/evotec_inventory_current.csv","r").read()
    dataframe = pd.read_csv(StringIO(csv_string))

    dataframe['Corporate ID'] = dataframe.apply(lambda row:batchId2CorporateId(row),axis = 1)

    engine = sqlalchemy.create_engine("sqlite://")
    sql.execute("drop table if exists data", engine)
    dataframe.to_sql("data",engine)

    screen_set_df = pd.read_sql_query("select * from data where "
                           "`Conc. Unit`='mM' and `Labware Type`='BGN_Vial_MLR96_(1.4mL)' "
                           "and `Available Amount Unit`='uL' and `Concentration`=10",engine)
    df = pd.read_sql_query("select *, Max(`Available Amount`) as MaxAmount from data where "
                            "`Conc. Unit`='mM' and `Labware Type`='BGN_Vial_MLR96_(1.4mL)' "
                            "and `Available Amount Unit`='uL' and `Concentration`=10 group by `Corporate Id` order by MaxAmount",engine)

    corporateIds = df[df['MaxAmount']<float(amount_uL)]['Corporate ID']
    wanted = screen_set_df[screen_set_df['Corporate ID'].isin(corporateIds)]
    return wanted.to_csv()


if __name__ == "__main__":
    print evotec_inventory_check(50)
