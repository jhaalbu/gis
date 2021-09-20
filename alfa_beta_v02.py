#Utvikla av Jan Helge Aalbu, Asplan Viak AS

#TODO: Innarbeide feilhåndtering
#TODO: Restrukturering, oppsettet med funksjoner blei rot..

import arcpy
import sys
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


def profil(inputfc, terreng, outputfc='punkter'):
    #Langer punkt kvar 1 meter langs input profil
    arcpy.GeneratePointsAlongLines_management(inputfc_profil, outputfc1, 'DISTANCE',
                                          Distance='1 Meters')
    #Legger z koordinater på punkter ut frå valt raster surface
    arcpy.ddd.AddSurfaceInformation(outputfc1, terreng, 'Z', 'BILINEAR')
    #Tar ut koordinatlister frå puntkter feature class 
    with arcpy.da.SearchCursor(outputfc, ["SHAPE", 'Z']) as cursor:
        x_list = []
        y_list = []
        z_list = []
        for row in cursor:
            x, y = row[0]
            x_list.append(x)
            y_list.append(y)
            z_list.append(row[1])
    #Etablerer Pandas dataframe for forenkling av vidare databehandling
    df = pd.DataFrame(list(zip(x_list, y_list, z_list)), columns =['X', 'Y', 'Z'])
    
    #Regner ut distanse mellom punkter (muligens unødvening, kan kanskje bruke index
    #til punkter istadenfor? Sidan kvar punkt er etalbert per meter?.
    df['DIST'] = np.sqrt(((df.X - df.X.shift(1))**2)+((df.Y - df.Y.shift(1))**2))
    df['M'] = df.DIST + df.DIST.shift(1)
    df.loc[0, 'M'] = 0
    df.loc[0, 'DIST'] = 0
    df.loc[1, 'M'] = df.loc[1, 'DIST']
    df.loc[0, 'H'] = 0
    
    #Regner ut lengden basert på avstand mellom punkter
    for i in range(2, len(df)):
        df.loc[i, 'M'] = df.loc[i-1, 'DIST'] + df.loc[i-1, 'M']

    #Runder av meterverdien 
    df['M'] = df['M'].round(0)
    
    
    return df


def poly_skredbane(df, polynom=2):
    '''
    input
    df: dataframe med profil, fra profil funksjonen
    polynom: størrelseorden på polynomet, standard 2. grads
    
    output
    returnerer ein dataframe med kolonnene POLY og H_DEG, som er tilpassa profil og helling i grader
    '''
    #Tar meterverdi og høgdeverdi for meterverdi (altså høgdeprofilet) og regner
    #om til ein tilpasse andregradspolynom (forenkler geometrien)
    p = np.poly1d(np.polyfit(df['M'], df['Z'], polynom))
    
    #Etablerer pandas kollonen poly for representerer det forenkla polynomet
    df['POLY'] = p(df['M'])

    #Regner ut hellingen langs det forenkla profilet
    for i in range(1, len(df)):
        df.loc[i, 'H'] = ((df.loc[i, 'POLY'] - df.loc[i -1, 'POLY'])/(df.loc[i, 'M'] - df.loc[i - 1, 'M']))

    #Regner om frå hellingstall til vinkel
    df['H_DEG'] = np.degrees(np.arctan(df['H']))
    
    return df, p

def betapunkt(df, skredtype='sno', avik=0.2):
    '''
    snøskred: 10 graders helling, steinsprang: 23 graders helling, jordskred: 20 grader helling
    input
    df: dataframe frå poly_skredbane
    vinkel: vinkel for betapunkt
    avik: aviket frå betapunktvinkel som blir godtatt
    

    output
    
    '''
    if skredtype == 'sno':
        vinkel = -10
       
    elif skredtype == 'stein':
        vinkel = -23

    else:
        vinkel = -20

    #Finner 10 graders vinkelpunkt

    df10 = df.loc[(df['H_DEG'] <= (vinkel+avik)) & (df['H_DEG'] >= (vinkel-avik))]
    index_label = df[(df['H_DEG'] <= (vinkel+avik)) & (df['H_DEG'] >= (vinkel-avik))].index.tolist()
    skjeringspunkt = index_label[0]

    #Finner "koordinater" for 10 graderspunktet og topp punkt av grafen for å regne beta vinkel

    m_topp = df.loc[0, 'M']
    poly_topp = df.loc[0, 'POLY']
    m_10 = df10.loc[skjeringspunkt, 'M']
    poly_10 = df10.loc[skjeringspunkt, 'POLY']

    #Koordinatpunkt for betapunkt
    beta_utlop_x = df.loc[m_10, 'X']
    beta_utlop_y = df.loc[m_10, 'Y']

    #Regner ut beta vinkelen
    beta_helning = abs((poly_10 - poly_topp)/(m_10 - m_topp))
    beta_vinkel_radianer = np.arctan(beta_helning)
    beta_vinkel_grader = np.degrees(beta_vinkel_radianer)    
    
    return (beta_helning, beta_vinkel_radianer, beta_vinkel_grader, beta_utlop_x, beta_utlop_y)

def alfa_vinkel(beta_vinkel_grader, skredtype='sno'):
    if skredtype == 'sno':
        standardavik = 2.3
        beta = 0.96
        justering = -1.4
        
    elif skredtype == 'stein':
        standardavik = 2.16
        beta = 0.77
        justering = 3.9

    else:
        standardavik = 1.5
        beta = 0.96
        justering = -4.0
    
    alfa_vinkel = (beta * beta_vinkel_grader) + justering  
    alfa_vinkel_SD1 = alfa_vinkel - standardavik
    alfa_vinkel_SD2 = alfa_vinkel - (standardavik * 2)
    alfa_vinkel_SD3 = alfa_vinkel - (standardavik * 3) 
    alfa_helning = np.tan(np.radians(alfa_vinkel))
    alfa_helning_SD1 = np.tan(np.radians(alfa_vinkel_SD1))
    alfa_helning_SD2 = np.tan(np.radians(alfa_vinkel_SD2))
        
    return alfa_vinkel, alfa_vinkel_SD1, alfa_vinkel_SD2, alfa_vinkel_SD3, alfa_helning, alfa_helning_SD1, alfa_helning_SD2
    

    
def skredutlop(df, alfa_helning, p):
    '''
    input
    df: dataframe frå poly_skredbane
    
    output:
    utlop_x
    utlop_y
    '''
    try:
        #Lager til data for til hjelp for å beregne skjæringspunkt mellom skredbane og linje med alfa vinkel
        liste_meterverdi = [0, df.loc[len(df)-1, 'M']]
        liste_alfa = [df.loc[0, 'POLY'], df.loc[0, 'POLY'] - df.loc[len(df)-1, 'M']*alfa_helning[0]]
        liste_alfa_sigma1 = [df.loc[0, 'POLY'], df.loc[0, 'POLY'] - df.loc[len(df)-1, 'M']*alfa_helning[1]]
        liste_alfa_sigma2 = [df.loc[0, 'POLY'], df.loc[0, 'POLY'] - df.loc[len(df)-1, 'M']*alfa_helning[2]]

        q_alfa = np.polyfit(liste_meterverdi, liste_alfa, 1)
        q_alfa_sigma1 = np.polyfit(liste_meterverdi, liste_alfa_sigma1, 1)
        q_alfa_sigma2 = np.polyfit(liste_meterverdi, liste_alfa_sigma2, 1)

        #Finer rotpunktet mellom skredbanen (p) (andregradspolynom) og alfa vinkelplanet (q)
        x_0_alfa = np.roots(p - q_alfa)
        alfa_verdi = int(x_0_alfa.max().round())

        x_0_alfa_sigma1 = np.roots(p - q_alfa_sigma1)
        alfa_sigma1_verdi = int(x_0_alfa_sigma1.max().round())

        x_0_alfa_sigma2 = np.roots(p - q_alfa_sigma2)
        alfa_sigma2_verdi = int(x_0_alfa_sigma2.max().round())
       
        #Regner ut koordinatpunkt for uløpspunktet
        alfa_utlop_x = df.loc[alfa_verdi, 'X']
        alfa_utlop_y = df.loc[alfa_verdi, 'Y']

        alfa_sigma1_utlop_x = df.loc[alfa_sigma1_verdi, 'X']
        alfa_sigma1_utlop_y = df.loc[alfa_sigma1_verdi, 'Y']

        alfa_sigma2_utlop_x = df.loc[alfa_sigma2_verdi, 'X']
        alfa_sigma2_utlop_y = df.loc[alfa_sigma2_verdi, 'Y']
        
    except ValueError:
        arcpy.AddError("Utløpet går utenfor lengden på profilet")
        sys.exit(1)
        
    return [(alfa_utlop_x, alfa_utlop_y), (alfa_sigma1_utlop_x, alfa_sigma1_utlop_y), (alfa_sigma2_utlop_x, alfa_sigma2_utlop_y)]

def utlop_feature(alfa, beta, skredtype):
    #alfa_punkt = arcpy.Point(alfa[0], alfa[1])
    #beta_punkt = arcpy.Point(beta[0], beta[1])
    #sigma1_punkt = arcpy.Point(alfa[0], alfa[1])
    #sigma2_punkt =  arcpy.Point(alfa[0], alfa[1])

    fc = 'alfabetapunkt_'+skredtype
    arcpy.CreateFeatureclass_management(fgdb, fc, "Point", "", "", "", 25833)
    arcpy.AddField_management(fc, 'NAVN', "TEXT")

    with arcpy.da.InsertCursor(fc, ["NAVN", "SHAPE@"]) as cursor:
        cursor.insertRow(('Alfa', (alfa[0][0], alfa[0][1])))
        cursor.insertRow(('Beta', (beta[0], beta[1])))
        cursor.insertRow(('Sigma1', (alfa[1][0], alfa[1][1])))
        cursor.insertRow(('Sigma2', (alfa[2][0], alfa[2][1])))

    #features = [fc_alfa, fc_beta, fc_sigma1, fc_sigma2]
    
    
    data = fgdb + "\\" + fc
    aprx = arcpy.mp.ArcGISProject("CURRENT")
    aprxMap = aprx.listMaps(aprx.activeMap.name)[0] 
    aprxMap.addDataFromPath(data)

    arcpy.ApplySymbologyFromLayer_management(data, "alfabetapunkt.lyrx")

#Setter workspace
fgdb = arcpy.GetParameterAsText(3)
arcpy.env.workspace = fgdb

#Inputdata til beregning
inputfc_profil = arcpy.GetParameterAsText(0)
input_skredtype = arcpy.GetParameterAsText(2)
outputfc1 = 'punkter'
insurface = arcpy.GetParameterAsText(1)

skredprofil = profil(inputfc_profil, insurface)
skredbane, p = poly_skredbane(skredprofil)
beta = betapunkt(skredbane, input_skredtype)
alfa = alfa_vinkel(beta[2])
utlop = skredutlop(skredbane, (alfa[4], alfa[5], alfa[6]), p)
betapunkt = (beta[3], beta[4])

utlop_feature(utlop, betapunkt, input_skredtype)

