#!/usr/bin/python -tt

import os
import os.path
import re
import time
import sys
import commands
import shutil



def main():
    #commands.getstatusoutput("alias molconvert='/env/cns/proj/agc/tools/LINUX64/MarvinBeans/bin/molconvert' ") 
	
	
	molconvert = sys.argv[1]
	
    
    path = ''
    indir = 'MetaCyc-MOLfiles'
    outdirP = 'MolFiles_FULL'
    outdirPA = 'MolFiles_FULL_aroma'
    
    if not os.path.exists(outdirP):
        os.makedirs(outdirP)
    if not os.path.exists(outdirPA):
        os.makedirs(outdirPA)


    #protonnation
    inFileList = os.listdir(path+''+indir)
    inFileList = [os.path.join(path+''+indir, f) 
        for f in inFileList 
            if os.path.splitext(f)[1] == '.mol']
            
    i=0       
    for molfile in inFileList:
        t = molfile.split('/')
        n=t[1]
        
        outfile = os.path.join(path+''+outdirP,n)
        commande = ''.join([molconvert,' mol:H ',molfile,' -o ',outfile])
        i=i+1
        print "protonation ",n,i,"/",len(inFileList)
        commands.getstatusoutput(commande)
    
    
    #aromatisation
    inFileList2 = os.listdir(path+''+outdirP)
    inFileList2 = [os.path.join(path+''+outdirP, f) 
        for f in inFileList2 
            if os.path.splitext(f)[1] == '.mol']
    
    i=0
    for molfile in inFileList2:

        t=molfile.split('/')
        n=t[1]
        
        outfile = os.path.join(path+''+outdirPA, n)
        commande = ''.join([molconvert,' mol:a ',molfile,' -o ',outfile])
        i=i+1
        print "aromatization ",n,i,"/",len(inFileList)
        commands.getstatusoutput(commande)
        
    
    
    
if __name__ == "__main__":
    main()



