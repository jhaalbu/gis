#ArcGIS pro script for å rekne ut alfa-beta utløp av skredtyper
#Utviklet av Jan Helge Aalbu, jan.aalbu@asplanviak.no
#Under utvikling, 18.10.2020
#Brukes gjennom script tool i ArcGIS pro

#TODO: Betapunkt/linje definert av bruker. Bruke arcgis til å finne skjæringspunkt?? Problem med å treffe meterverdi? Nærmeste punkt??
#TODO: Rekne på profil, og ikkje polylinje? Eller berre bruke høggrads polynom?



import arcpy
import sys
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import random
from datetime import datetime

class Profil:
    def __init__(self, inputfc, terreng, tempfc_punkter):
        self.inputfc = inputfc
        self.terreng = terreng
        self.tempfc_punkter = tempfc_punkter
        #Langer punkt kvar 1 meter langs input profil
        arcpy.GeneratePointsAlongLines_management(inputfc, tempfc_punkter, 'DISTANCE',
                                            Distance='1 Meters')
        #Legger z koordinater på punkter ut frå valt raster surface
        arcpy.ddd.AddSurfaceInformation(tempfc_punkter, terreng, 'Z', 'BILINEAR')
        #Tar ut koordinatlister frå puntkter feature class 
        with arcpy.da.SearchCursor(tempfc_punkter, ["SHAPE", 'Z']) as cursor:
            x_list = []
            y_list = []
            z_list = []
            for row in cursor:
                x, y = row[0]
                x_list.append(x)
                y_list.append(y)
                z_list.append(row[1])
        #Etablerer Pandas dataframe for forenkling av vidare databehandling
        self.df = pd.DataFrame(list(zip(x_list, y_list, z_list)), columns =['X', 'Y', 'Z'])

        #Regner ut distanse mellom punkter (muligens unødvening, kan kanskje bruke index
        #til punkter istadenfor? Sidan kvar punkt er etalbert per meter?.
        self.df['DIST'] = np.sqrt(((self.df.X - self.df.X.shift(1))**2)+((self.df.Y - self.df.Y.shift(1))**2))
        self.df['M'] = self.df.DIST + self.df.DIST.shift(1)
        self.df.loc[0, 'M'] = 0
        self.df.loc[0, 'DIST'] = 0
        self.df.loc[1, 'M'] = self.df.loc[1, 'DIST']
        self.df.loc[0, 'H'] = 0
        
        #Regner ut lengden basert på avstand mellom punkter
        for i in range(2, len(self.df)):
            self.df.loc[i, 'M'] = self.df.loc[i-1, 'DIST'] + self.df.loc[i-1, 'M']

        #Runder av meterverdien 
        self.df['M'] = self.df['M'].round(0)

    def poly(self, polynom=2):
        '''
        input
        df: dataframe med profil, fra profil funksjonen
        polynom: størrelseorden på polynomet, standard 2. grads
        
        output
        returnerer ein dataframe med kolonnene POLY og H_DEG, som er tilpassa profil og helling i grader
        '''
        self.polynom = polynom
        #Tar meterverdi og høgdeverdi for meterverdi (altså høgdeprofilet) og regner
        #om til ein tilpasse polynom (forenkler geometrien)
        self.p = np.poly1d(np.polyfit(self.df.index, self.df['Z'], self.polynom))
        #Etablerer pandas kollonen poly for representerer det forenkla polynomet
        self.df['POLY'] = self.p(self.df.index)
        #Regner ut hellingen langs det forenkla profilet
        for i in range(1, len(self.df)):
            self.df.loc[i, 'H'] = ((self.df.loc[i, 'POLY'] - self.df.loc[i -1, 'POLY'])/(self.df.loc[i, 'M'] - self.df.loc[i - 1, 'M']))

        #Regner om frå hellingstall til vinkel
        self.df['H_DEG'] = np.degrees(np.arctan(self.df['H']))
        
    def get_profil(self):
        return self.df        

    def plot_profil(self):
        plt.plot(self.df['M'], self.df['Z'])
        plt.title(f'Høydeprofil {self.inputfc}')
        plt.xlabel("Lengde (m)")
        plt.ylabel("Høyde (m)")
        plt.grid(True)
        plt.show()

class Skred:
    '''
    Eit skredobjekt kan berre ha ein skredtype
    Profil input må vere ei dataframe frå Profil objekt
    Fra get_profil metoden
    '''

    def __init__(self, profil):
        self.profil = profil
        self.df = profil.df
        self.x_start = self.df['X'][0]
        self.y_start = self.df['Y'][0]
        self.z_topp = self.df['Z'][0]
        self.m_topp = self.df.loc[0, 'M']
        self.poly_topp = self.df.loc[0, 'POLY']
        self.m_ende = self.df['M'].iloc[-1]
        self.z_ende = self.df['Z'].iloc[-1]

    def skyggevinkel(self, analysetype):
        if analysetype == 'F':
            self.vinkler = [28,30,32,35]
        elif analysetype == 'M':
            self.vinkler = [25,28,30,32]
        
        profil = self.profil
        df = self.df

        hellingsliste = []
        for i in self.vinkler:
            hellingsliste.append(math.tan(math.radians(i)))
  
        #Liste med verdiene for lengdemeter fra 
        liste_meterverdi = [0, self.m_ende]
        self.koordinater = []
        self.plotverdier = []
                                 
        for helling in hellingsliste:
            liste_utlop = [self.poly_topp, self.poly_topp - self.m_ende*helling]
            print(liste_utlop)

            q = np.polyfit(liste_meterverdi, liste_utlop, 1)
            p = profil.p

            nullpunkt = np.roots(p - q)
            nullpunktsliste = list(nullpunkt)
            nullpunkt_sortert = sorted(nullpunktsliste)
            #print(f'før pop: {nullpunkt_sortert}')

            #while max(nullpunkt_sortert) > len(df):
            #    nullpunkt_sortert.pop(-1)  

            #print(f'etter pop: {nullpunkt_sortert}')
            #print(max(nullpunkt_sortert))    

            filtered = [i for i in nullpunkt_sortert if i.imag == 0]
            print(filtered)
            if max(filtered) > len(df):
                filtered.pop(-1)
            print(filtered)
            nullpunkt_m = int(max(filtered).real)
            print(f'nullpunkt_m:{nullpunkt_m}')
            #print(df)
            utlop_x = self.df.loc[nullpunkt_m, 'X']
            utlop_y = self.df.loc[nullpunkt_m, 'Y']
            m = self.df.loc[nullpunkt_m, 'M']
            poly = self.df.loc[nullpunkt_m, 'POLY']
            self.koordinater.append((utlop_x, utlop_y))
            self.plotverdier.append((m, poly))

        return self.koordinater, self.plotverdier

    def runout(self, sigma=5):
        skredtype = self.skredtype
        alfaparameter = {'sno':[2.3, 0.96, 1.4], 'stein':[2.16, 0.77, -3.9], 'jord':[1.5, 0.96, 4.0]}
        sd = alfaparameter[self.skredtype][0]
        tilpassing = alfaparameter[self.skredtype][1]
        justering = alfaparameter[self.skredtype][2]
        print(f'tilpassing er {tilpassing}, justering er {justering}, standardavik er {sd}, betavinkel er {self.beta_vinkel_grader}')
        self.alfa_vinkelliste = []
        self.alfa_hellingsliste = []
        for i in range(5):
            alfa_vinkel = (tilpassing * self.beta_vinkel_grader) - justering  - (sd * i)
            alfa_helning = np.tan(np.radians(alfa_vinkel))
            self.alfa_vinkelliste.append(alfa_vinkel)
            self.alfa_hellingsliste.append(alfa_helning)
        liste_meterverdi = [0, self.df.loc[len(self.df)-1, 'M']]
        sigmateller = 0
        self.alfa_koordinater = []
        self.alfa_plotverdier = []
        
        while True and sigmateller <= sigma:
            try:
                liste_alfa = [self.df.loc[0, 'POLY'], self.df.loc[0, 'POLY'] - self.df.loc[len(self.df)-1, 'M']*self.alfa_hellingsliste[sigmateller]]

                #Finer rotpunktet mellom skredbanen (p) (andregradspolynom) og alfa vinkelplanet (q)
                q_alfa = np.polyfit(liste_meterverdi, liste_alfa, 1)
                x_0_alfa = np.roots(self.profil.p - q_alfa)
                x_0_alfa_list = list(x_0_alfa)
                x0_sorted = sorted(x_0_alfa_list)
                print('x0_sorted', x0_sorted)
                print('max x0:, ', max(x0_sorted))
                print('len df', len(self.df))
                #Logikk for å håndtere fleire rotpunkter ved høgaregradspolynom
                # if max(x0_sorted) > len(self.df) and max(x0_sorted) <= 0:                                
                #     while max(x0_sorted) > len(self.df):
                #         x0_sorted.pop(-1)
                #     alfa_verdi = int(max(x0_sorted).round())
                # else:
                #     alfa_verdi = int(max(x0_sorted).round())
                                               
                #while max(x0_sorted) > len(self.df):

                 #   x0_sorted.pop(-1)
                
                #if max(x0_sorted) > len(self.df) and max(x0_sorted) > self.alfa_plotverdier[-1][0]: 
                 #   alfa_verdi = int(max(x0_sorted).round())
                #else:
                 #   alfa_verdi = int(max(x0_sorted).round())

                filtered = [i for i in x0_sorted if i.imag == 0]
                alfa_verdi = int(max(filtered).round())

                alfa_utlop_x = self.df.loc[alfa_verdi, 'X']
                alfa_utlop_y = self.df.loc[alfa_verdi, 'Y']
                alfa_m = self.df.loc[alfa_verdi, 'M']
                alfa_poly = self.df.loc[alfa_verdi, 'POLY']
                self.alfa_koordinater.append((alfa_utlop_x, alfa_utlop_y))
                self.alfa_plotverdier.append((alfa_m, alfa_poly))
                
            except:
                break
            sigmateller += 1
        #print(f'alfa koorinater = {self.alfa_koordinater} .. alfa_plotverdier = {self.alfa_plotverdier}')
        print(f'alfa_plotverdier = {self.alfa_plotverdier}')
        return self.alfa_koordinater, self.alfa_plotverdier



def lag_featurepunkt(skred, fgdb, profilnavn):
    now = datetime.now()
    tid = now.strftime("%H%M%S")
    fc = profilnavn + '_steinsprang' + tid
    arcpy.CreateFeatureclass_management(fgdb, fc, "Point", "", "", "", 25833)
    arcpy.AddField_management(fc, 'NAVN', "TEXT")
    arcpy.AddField_management(fc, 'Punkttype', "TEXT")


    with arcpy.da.InsertCursor(fc, ["Punkttype","NAVN", "SHAPE@"]) as cursor:
        cursor.insertRow(('0','Start', (skred.x_start, skred.y_start)))
        for i in range(len(skred.koordinater)):
            cursor.insertRow((str(i+1),skred.vinkler[i], (skred.koordinater[i][0], skred.koordinater[i][1])))
    
    return fc
                             
def feature_paa_kart(fgdb, fc): 
    data = fgdb + "\\" + fc
    aprx = arcpy.mp.ArcGISProject("CURRENT")
    aprxMap = aprx.listMaps(aprx.activeMap.name)[0] 
    aprxMap.addDataFromPath(data)
    #fc.symbology = "C:/Users/jan.aalbu/Documents/Koding/alfabeta/Alfa-beta_punkter.lyrx"
    #arcpy.ApplySymbologyFromLayer_management(data, symbologyLayer)


def plot_alfa(profil, skred):
    df = profil.df
    
    # arcpy.AddMessage(str(alfa))
    fig, ax = plt.subplots(figsize=(12,8))
    ax.plot(df['M'], df['Z'], label='Høgdeprofil') #Høgdeprofilet
    #ax.plot(df['M'], df['POLY'], label=f'Tilpasset profil {profil.polynom}. grads') #Forenkla høgdeprofil
    #ax.scatter(beta[5], beta[6], color='r', linewidth='1', label='Punkt med 10 grader helling') # 10 graders punkter
    #ax.plot([df['M'][0], skred.m_beta], [df['POLY'][0], skred.z_beta], label=f'Beta {round(skred.beta_vinkel_grader, 1)}\xb0', linestyle='--') #Beta 
    #ax.plot([df['M'][0], skred.alfa_plotverdier[0][0]], [df['POLY'][0], skred.alfa_plotverdier[0][1]], label=f'Alfa {round(skred.alfa_vinkelliste[0], 1)}\xb0') #Må plotte til skjæeringspunkt
    #loc = ticker.MultipleLocator(base=100.0)
    #ax.xaxis.set_major_locator(loc)
    start, end = ax.get_xlim()
    start = 0
    end = round(end, -1)
    ax.xaxis.set_ticks(np.arange(start, end, 20))
    start, end = ax.get_ylim()
    start = round(start, -1)
    end = round(end, -1)
    ax.yaxis.set_ticks(np.arange(start, end, 20))
    ax.set_aspect('equal')
    ax.set_xlabel('Avstand (m)')
    ax.set_ylabel('Høyde (m.o.h.)')
    #ax.grid(True)
    ax.legend()
    return ax

def plot_alfa_png(profil, skred, folder, fc):
    df = profil.df
    
    # arcpy.AddMessage(str(alfa))
    fig, ax = plt.subplots(figsize=(12,8))
    ax.plot(df['M'], df['Z'], label='Høgdeprofil') #Høgdeprofilet
    ax.plot(df['M'], df['POLY'], label=f'Tilpasset profil {profil.polynom}. grads') #Forenkla høgdeprofil
    #ax.scatter(beta[5], beta[6], color='r', linewidth='1', label='Punkt med 10 grader helling') # 10 graders punkter
    ax.plot([df['M'][0], skred.m_beta], [df['POLY'][0], skred.z_beta], label=f'Beta {round(skred.beta_vinkel_grader, 1)}\xb0', linestyle='--') #Beta 
    ax.plot([df['M'][0], skred.alfa_plotverdier[0][0]], [df['POLY'][0], skred.alfa_plotverdier[0][1]], label=f'Alfa {round(skred.alfa_vinkelliste[0], 1)}\xb0') #Må plotte til skjæeringspunkt
    #loc = ticker.MultipleLocator(base=100.0)
    #ax.xaxis.set_major_locator(loc)
    start, end = ax.get_xlim()
    start = 0
    end = round(end, -1)
    ax.xaxis.set_ticks(np.arange(start, end, 20))
    start, end = ax.get_ylim()
    start = round(start, -1)
    end = round(end, -1)
    ax.yaxis.set_ticks(np.arange(start, end, 20))
    ax.set_aspect('equal')
    ax.set_xlabel('Avstand (m)')
    ax.set_ylabel('Høyde (m.o.h.)')
    #ax.grid(True)
    ax.legend()
    print(folder)
    plt.savefig(str(folder)+'/'+fc, dpi=600) 
    plt.close()
    del fig
    
#Setter parameter fra ArcGIS
fgdb = arcpy.GetParameterAsText(0)
inputfc_profil = arcpy.GetParameterAsText(1)
analysetype = arcpy.GetParameterAsText(2)
insurface = arcpy.GetParameterAsText(3)

plott = ' '
polynom = 8

#C:\Users\jan.aalbu\OneDrive - Asplan Viak\Documents\ArcGIS\Projects\Vestnes\Vestnes.gdb
#Setter parameter for testing
    
##fgdb = "C:/Users/jan.aalbu/OneDrive - Asplan Viak/Documents/ArcGIS/Projects/Vestnes/Vestnes.gdb"
##input_skredtype = 'sno'
##inputfc_profil = 'skredlinje_3'
##insurface = 'C:/Users/jan.aalbu/OneDrive - Asplan Viak/Documents/ArcGIS/Projects/Vestnes/SjeggstadHøydemodell.tif'
##plott = "C:/Users/jan.aalbu/OneDrive - Asplan Viak/Documents/ArcGIS/Projects/Vestnes/Plot"
##standardavik = 2
##polynom = 5
##analysetype = 'F'


outputfc1 = 'alfa_temp_punkter_xvb1'
desc = arcpy.Describe(fgdb +'/'+ inputfc_profil)
profilnavn = desc.name
#Setter workspace
arcpy.env.workspace = fgdb

arcpy.management.Delete(fgdb +'/'+ outputfc1)

#linjer = []
# cursor = arcpy.da.SearchCursor(fgdb+'/'+inputfc_profil, ["SHAPE@"])
# for row in cursor:
#     print(row)
#     linjer.append(row)

#Lager fyste profil
#profil1 = Profil(inputfc_profil, insurface, outputfc1)
#profil1 = Profil(linjer[0], insurface, outputfc1)
#Etablerer tilpassa skredbane
#profil1.poly(polynom)
#Lager skred
#skred = Skred(profil1, input_skredtype)

# print(skred.x_beta, snoskred.y_beta)
#skred.runout(standardavik)
# print(skred.alfa_plotverdier)
#fc = lag_featurepunkt(skred, fgdb, profilnavn)

#feature_paa_kart(fgdb, fc)

# if plott:
#     plot_alfa(profil1, skred)
#     plt.show()


def alfa_beta(fgdb, inputfc_profil, insurface, outputfc1, polynom, analysetype, plott):
    arcpy.management.Delete(fgdb +'/'+ outputfc1)
    linjer = []
    cursor = arcpy.da.SearchCursor(fgdb+'/'+inputfc_profil, ["SHAPE@"])
    for row in cursor:
        print(row)
        linjer.append(row)
    
    for linje in linjer:
        arcpy.management.Delete(fgdb +'/'+ outputfc1)
        profil = Profil(linje, insurface, outputfc1)
        profil.poly(polynom)
        skred = Skred(profil)
        skred.skyggevinkel(analysetype)
        fc = lag_featurepunkt(skred, fgdb, profilnavn)
        #plot_alfa_png(profil, skred, plott, fc)
        feature_paa_kart(fgdb, fc)




alfa_beta(fgdb, inputfc_profil, insurface, outputfc1, polynom, analysetype, plott)


