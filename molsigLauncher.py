#!/usr/bin/python -tt

import os
import os.path
import re
import time
import sys
import commands
import shutil



def main():
    
    #MAIN PROGRAMM
	#for each file in MOLfile directory launches molsig on five heights

    
    #arguments
    #1: MolFiles_FULL / MolFiles_FULL_aroma (MOLfiles directory)
    #2: scan/sscan/fsscan (molsig type)
    #3: o/n (with or without aromatization)
    
    
    path = ''
    dir = sys.argv[1] 
    scantype = sys.argv[2]
    aroma = sys.argv[3]
    
    rayons = [0,2,4,6,8,10]
    
    fileList = os.listdir(path+''+dir)
    fileList = [os.path.join(path+''+dir, f) 
        for f in fileList 
            if os.path.splitext(f)[1] == '.mol']
    
    
    
    # launch molsig on 5 heights
    count=0
    for molfile in fileList :
        for i in rayons :
            commande = path+'molsig-code/bin/sscan '+molfile+' '+scantype+str(i)
            print commande
            commands.getstatusoutput(commande)
        #print "done "+molfile
        count+=1
    
    #check
    print count
    
    
    print "moving files to correct directory"
    # launch of the programm which will put in order the created files
    (status, output) = commands.getstatusoutput('python file_mover.py '+dir+' '+scantype+' '+aroma)
    


if __name__ == "__main__":
    main()





