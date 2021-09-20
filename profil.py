import arcpy
import pandas as pd
import matplotlib.pyplot as plt
import time


def profil(inputfc, insurface, out_table):
    arcpy.ddd.StackProfile(inputfc, insurface, out_table)
    with arcpy.da.SearchCursor(out_table, ["FIRST_DIST", 'FIRST_Z']) as cursor:
        m_list = []
        z_list = []
        
        for row in cursor:
            m, z = row
            m_list.append(m)
            z_list.append(z)
              
        #Etablerer Pandas dataframe for forenkling av vidare databehandling
    df = pd.DataFrame(list(zip(m_list, z_list)), columns =['M', 'Z'])
    return df

#fgdb = arcpy.GetParameterAsText(2)
#arcpy.env.workspace = fgdb

#Inputdata til beregning
inputfc = arcpy.GetParameterAsText(0)
insurface = arcpy.GetParameterAsText(1)
tittel = arcpy.GetParameterAsText(2)
xlabel = arcpy.GetParameterAsText(3)
ylabel = arcpy.GetParameterAsText(4)

out_table = "temp_tabell" + str(time.time())[-4:]

df = profil(inputfc, insurface, out_table)

plt.plot(df['M'], df['Z'])
#plt.plot(df['M'], df['snitt'])
plt.title(tittel)
#plt.xlabel("Lengde (m)")
plt.xlabel(xlabel)
#plt.ylabel("Høyde (m)")
plt.ylabel(ylabel)
plt.grid(True)
plt.show()

