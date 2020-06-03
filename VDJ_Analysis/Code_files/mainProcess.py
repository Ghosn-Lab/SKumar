# -*- coding: utf-8 -*-
"""
@author: scsac
"""

import sys, os
from shutil import copy
import BM3_vdj_integrate as intVDJ
import cdr3_aa_viz as cdr3AA
import circos_Mat_Generator as circosG
import csvToDict_v5 as mutProc

def main():
    if len(sys.argv) == 1:
        print("Generate files for VDJ visualizations")
        print("Usage: python3 mainProcess.py [1234] [filename(s)] [out_directory]")
        print("1:CDR3 AA, 2:ADT integration, 3:Circos, 4:CDR3 Mutations")
        print("Inputs for the modes in the order specified:")
        print("Mode 1: VDJ filtered matrix and output directory")
        print("Mode 2: ADT normalized matrix, IMGT processed mutation file, IMGT processed VDJ matrix and output directory")
        print("Mode 3: IMGT processed VDJ matrix and output directory")
        print("Mode 4: IMGT processed mutation file, IMGT processed VDJ matrix and output directory")
    else:
        mode = int(sys.argv[1])
        if mode == 1:
            inFName = sys.argv[2]
            outDir = sys.argv[3]
            cdr3AA.cdr3AAData(inFName, outDir)
        elif mode == 2:
            inFName = sys.argv[2]
            inFNames = []
            inFNames.append(sys.argv[3])
            inFNames.append(sys.argv[4])
            outDir = sys.argv[5]
            mutFNames = mutProc.mutationDataProcessing(inFNames, outDir)
            outFName = os.path.basename(inFName)
            copy(inFName, os.path.join(outDir, outFName))
            intVDJ.integrateMutations(inFName, sorted(mutFNames), os.path.join(outDir, outFName))
        elif mode == 3:
            inFName = sys.argv[2]
            outDir = sys.argv[3]
            circosG.circosFileGenerator(inFName, outDir)
        elif mode == 4:
            inFNames = []
            inFNames.append(sys.argv[2])
            inFNames.append(sys.argv[3])
            outDir = sys.argv[4]
            mutProc.mutationDataProcessing(inFNames, outDir)
        else:
            print("Wrong mode of operation! Provide 1/2/3/4")
    
main()