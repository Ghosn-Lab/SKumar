# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 11:50:15 2020

@author: scsac
"""

import csv, json, os

def circosFileGenerator(inFileName, outDir):
    filename = inFileName
    
    isotypeGroups = set()
    jsonContentUnSorted = {}
    
    with open(filename, 'r') as csvH:
        csvReader = csv.reader(csvH)
        for lineageData in csvReader:
            if lineageData[1] not in jsonContentUnSorted:
                jsonContentUnSorted[lineageData[1]] = []
                jsonContentUnSorted[lineageData[1]].append(lineageData[2])
            else:
                jsonContentUnSorted[lineageData[1]].append(lineageData[2])
            isotypeGroups.add(lineageData[1])
        
    isotypeGroups = sorted(list(isotypeGroups)[1:])
    #del jsonContentUnSorted['isotype']
    
    jsonContentSorted = []
    for isotype in isotypeGroups:
        jsonContentSorted.append({'name': isotype, 'matrix': [0]*len(isotypeGroups)})
    
    sharedLineageIDs = {}
    i = 0
    while i < len(isotypeGroups):
        for isotype in isotypeGroups[i+1:]:
            sharedLI = set(jsonContentUnSorted[isotypeGroups[i]]).intersection(set(jsonContentUnSorted[isotype]))
            if len(sharedLI) != 0:
                if isotypeGroups[i] not in sharedLineageIDs:
                    sharedLineageIDs[isotypeGroups[i]] = list(sharedLI)
                else:
                    sharedLineageIDs[isotypeGroups[i]] += list(sharedLI)
                iso1count = 0
                iso2count = 0
                for lineage in sharedLI:
                    iso1count += jsonContentUnSorted[isotypeGroups[i]].count(lineage)
                    iso2count += jsonContentUnSorted[isotype].count(lineage)
                isoItem1 = list(isoT['matrix'] for isoT in jsonContentSorted if isoT['name'] == isotypeGroups[i])
                isoItem2 = list(isoT['matrix'] for isoT in jsonContentSorted if isoT['name'] == isotype)
                isoItem1[0][isotypeGroups.index(isotype)] = iso1count
                isoItem2[0][isotypeGroups.index(isotypeGroups[i])] = iso2count
        i += 1
    
    for isotype in isotypeGroups:
        lineageCount = 0
        if isotype in sharedLineageIDs:
            totalNonSharedLIds = set(jsonContentUnSorted[isotype]) - set(sharedLineageIDs[isotype])
            lineageCount = len(totalNonSharedLIds)
        else:
            lineageCount = len(set(jsonContentUnSorted[isotype]))
        isoItem = list(isoT['matrix'] for isoT in jsonContentSorted if isoT['name'] == isotype)
        isoItem[0][isotypeGroups.index(isotype)] = lineageCount
    
    with open(os.path.join(outDir, "circos.json"), 'w') as wFile:
        wFile.write(json.dumps(jsonContentSorted))