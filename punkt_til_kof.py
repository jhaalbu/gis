# Python script: punkt_til_kof.py
# Author: Jan Helge Aalbu, Asplan Viak
# Scriptet tar inn ein feature class med punkter, og skrive ut ei .kof fil med
# navnet <featureclass navn>_export.kof

import arcpy
print("Arcpy importert")
inputfc = arcpy.GetParameterAsText(0)
output_folder = arcpy.GetParameterAsText(1)
#Hardkoder inn typenr og sosikoden
typenr = '05'
sosikode = '2200'
z = '0.0'

with arcpy.da.SearchCursor(inputfc, ["SHAPE@XY", "ID"]) as cursor:
    with open(output_folder + '\\' + inputfc + '_export.kof', 'w') as f: #Setter navn på fil til navn paa featureclass + _export.kof
        for row in cursor:      
            x,y = row[0]
            #z = row[2]
            punktnr = row[1]
            print(f"{typenr:^4}{punktnr:<11}{sosikode:<12}{y:<13.6f} {x:<13.6f} {z}", file=f) #printer linje for linje til fil
