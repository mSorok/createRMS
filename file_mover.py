#!/usr/bin/python -tt

import os
import os.path
import re
import time
import sys
import commands
import shutil



def main():
    
    
    #option 1 : repertoire ou se trouvent les fichiers mol et les fichiers de resultats
    #option 2: le type de scan (sscan, scan ou fscan)
    #option 3: o ou n (si oui on pas fichiers aromatiques)
    
    molfilesdir = sys.argv[1]
    scantype = sys.argv[2]
    aroma = sys.argv[3]
    

    if aroma=='o':
        scandir='a'+scantype
    else:
        scandir=scantype
    

    rayons = [0,2,4,6,8,10]
    path = ''
    # creation de repertoires pour chaque diametre
    try:
        for i in rayons:
            #shutil.rmtree('mol-sig-results/sscan'+str(i))
            os.mkdir(path+'mol-sig-results/'+scandir+str(i))
    except OSError :
        sys.stderr.write('Dir already exists')

    
    source = os.listdir(path+molfilesdir+'/')
    for files in source:
        for i in rayons:
            if files.endswith("."+scantype+str(i)):
                shutil.move(path+molfilesdir+'/'+files,path+'mol-sig-results/'+scandir+str(i))
    
    
    


if __name__ == "__main__":
    main()


