
#!/usr/bin/python -tt

import os
import os.path
import re
import time
import sys
import commands




def main():

    # "scan" type (scan, sscan or fscan)
    if sys.argv[1] != '':
        molSigType = sys.argv[1]
    else:
        molSigType = 'sscan'
     
    
    if sys.argv[2] != '':
        aroma = sys.argv[2]
    else:
        aroma = 'n'
    
    
    molsigdir = molSigType
    if aroma == 'o':
        molsigdir = 'a'+molSigType
        
    
    sys.stderr.write(molsigdir+"\n")
    
    rayons = [0,2,4,6,8,10]
    
    sys.stderr.write("beguin\n")
    
    # transformes the .dat file to a tabulated file with only one reaction by line
    commands.getstatusoutput('python reactions_normalizer.py')    
           

    
    reaFile = open('reactions_tab','r')
    
    orphanCompounds = set([])
    reactionsPassed = set([])
    
    sys.stderr.write("reactions file opened\n")
    for line in reaFile:
        data = line.split("\t")
        frame = data[0]
        right = data[3].split('$')
        left = data[1].split('$')
        direction = data[2]
        
        #sys.stderr.write("line "+frame+"\n")
        
        
        #if all compounds are known
        if  data[1].find('UNKNOWN')== -1 and data[3].find('UNKNOWN') == -1 :
            
            substratesTab = []
            productsTab = []
            
            if direction == 'LEFT-TO-RIGHT':
                substratesTab = left
                productsTab = right
            elif direction == 'PHYSIOL-LEFT-TO-RIGHT':
                direction = 'LEFT-TO-RIGHT'
                substratesTab = left
                productsTab = right
            elif direction == 'IRREVERSIBLE-LEFT-TO-RIGHT':
                direction = 'LEFT-TO-RIGHT'
                substratesTab = left
                productsTab = right
            elif direction == 'RIGHT-TO-LEFT':
                substratesTab = right
                productsTab = left
            elif direction == 'IRREVERSIBLE-RIGHT-TO-LEFT':
                direction = 'RIGHT-TO-LEFT'
                substratesTab = right
                productsTab = left
            elif direction == 'PHYSIOL-RIGHT-TO-LEFT':
                direction = 'RIGHT-TO-LEFT'
                substratesTab = right
                productsTab = left
            elif direction == 'REVERSIBLE':
                substratesTab = left
                productsTab = right
            elif direction == '':
                direction = 'REVERSIBLE'
                substratesTab = left
                productsTab = right
            else:
                substratesTab = left
                productsTab = right
            
            
            
            substratesDic = {}
            productsDic = {}
            
            
            #substrates curation
            
            halfmolS = False
            
            for st in substratesTab:
                
                st = st.replace('\n','')
                st = st.replace('|','')
                elem = st.split('*')
                c = elem[0] #the compound
                c = c.replace('\n','')

                
                coef = 1
                if len(elem)> 1 :
                    coef = elem[1]
                if str(coef).find('N')>-1:
                    coef = 3    
                if coef=='m+q':
                    coef = 3 #generic value   
                if coef=='n':
                    coef = 3 #generic value
                if coef == '2n':
                    coef = 6
                if coef == '3n':
                    coef = 9
                if coef == '4n':
                    coef = 12
                if coef == '(n+1)':
                    coef = 4
                
                    
                    
                if coef == 0.5:  #half-oxygen molecule case
                    halfmolS =True
                    
                #coef = int(coef)

                #adjustments
                if c == 'NAD-P-OR-NOP':
                    c='NADP'
                c = c.replace('|','')
                

                
                if c != 'E-':
                    substratesDic[c]=coef
                
            
            
            
            #products curation
            
            halfmolP = False
            
            for pt in productsTab:
                pt = pt.replace('\n','')
                pt = pt.replace('|','')
                elem = pt.split('*')
                c = elem[0] 
                c = c.replace('\n','')
                
                
                coef = 1
                if len(elem)> 1 :
                    coef = elem[1]
                if str(coef).find('N')>-1:
                    coef = 3    
                if coef=='m+q':
                    coef = 3 #generic value    
                if coef=='n':
                    coef = 3 #generic value
                if coef == '2n':
                    coef = 6
                if coef == '3n':
                    coef = 9
                if coef == '4n':
                    coef = 12
                if coef == '(n+1)':
                    coef = 4
                
                    
                if coef == 0.5: #hal-oxygen molecule case
                    halfmolP =True
                    
                #coef = int(coef)
                    
                #adjustemnts
                if c == 'NAD-P-OR-NOP':
                    c='NADP'
                    
                c = c.replace('|','')
                
                if c != 'E-':
                    productsDic[c]=coef
                
            # coefficient regularization in case of half molecules -> everything is multiplied by 2
            if halfmolS == True or halfmolP == True:
                for c in substratesDic.keys():
                    substratesDic[c] = substratesDic[c]*2
                    substratesDic[c] = int(substratesDic[c])
                for c in productsDic.keys():
                    productsDic[c] = productsDic[c]*2
                    productsDic[c] = int(productsDic[c])
                
            
            # determination if all compound files exist 
            allFilesExist = True
            for k in substratesDic.keys():
                if not os.path.exists('mol-sig-results/'+molsigdir+'10/'+k+'.'+molSigType+'10'):
                    allFilesExist = False
                   

            for k in productsDic.keys():
                if not os.path.exists('mol-sig-results/'+molsigdir+'10/'+k+'.'+molSigType+'10'):
                    allFilesExist= False
            
            
            
            if allFilesExist and len(substratesDic.keys()) != 0 and len(productsDic.keys()) != 0:
                for i in rayons:
                    
                    #path where are all the molsig files at the given rayon
                    source = 'mol-sig-results/'+molsigdir+str(i)
                    substratesSignature = {}
                    productsSignature = {}
                
                
                    # molsigs recuperation for substrates 
                    for c in substratesDic.keys():
                        filename = source+'/'+c+'.'+molSigType+str(i)
                        try:
                            f = open(filename,"r")
                            
                            for line in f:
                                match = re.search(r"^0.0.*",line)
                                if not match :
                                    
                                    l = line.split(' ') 
                                    motif = l[1]
                                    motif = motif[:-1]
                                    
                                    motif = motif.replace('-1','')
                                    motif = motif.replace('+1','')
                                    
                                    occurency = float(l[0])*float(substratesDic[c])
                                    
                                    if motif in substratesSignature.keys():
                                        substratesSignature[motif] = substratesSignature[motif]+occurency
                                    else:
                                        substratesSignature[motif] = occurency
                                        
                        except IOError:
                            sys.stderr.write(c+"\n")
            
                        
                    
                    
                    
                    # recovery of molsigs for products
                    for c in productsDic.keys():
                        filename = source+'/'+c+'.'+molSigType+str(i)
                        try:
                            f = open(filename,"r")
                            
                            for line in f:
                                match = re.search(r"^0.0.*",line)
                                if not match :
                                    
                                    l = line.split(' ')
                                    motif = l[1]
                                    motif = motif[:-1]
                                    
                      
                                    motif = motif.replace('-1','')
                                    motif = motif.replace('+1','')
                                    
                                    occurency = float(l[0])*float(productsDic[c])
                                    
                                    if motif in productsSignature.keys():
                                        productsSignature[motif] = productsSignature[motif]+occurency
                                    else:
                                        productsSignature[motif] = occurency
                                        
                        except IOError:
                            sys.stderr.write(c+"\n")
                    
                    
                    
                    #RMS computation
                    # substraction: products - substrates
                    reactionSignature = {}
                    #3 cases for subparts:
                    # -presence in both dictionnaries
                    # - present in roducts and absent in substrates
                    # - present in substrates and absent in products
                    
                    for p in productsSignature.keys():
                        if p in substratesSignature.keys():
                            #first case : 
                            reactionSignature[p] =productsSignature[p]-substratesSignature[p]
                        else:
                            #second case :
                            reactionSignature[p] =productsSignature[p]
                
                    #third cases :
                    for p in substratesSignature.keys():
                        if p not in productsSignature.keys():
                            reactionSignature[p]= 0 - substratesSignature[p]
                
                    # ordering of keys for the reaction signatures so they are in the same order, even for different reactions

                    parts = sorted(reactionSignature.keys())
                
                    lineOut = ''

                    lineOut= str(i/2)+"\t"+str(frame)+"\t"
                    for j in parts:
                        if reactionSignature[j] != 0.0 :
                            lineOut = lineOut+str(reactionSignature[j])+"*"+str(j)+'$'
                    lineOut= lineOut+"0.0"
                
                    print lineOut
                    
            #else:
                 #sys.stderr.write("all files missing\n")
                    
                    


    reaFile.close()

 
if __name__ == "__main__":
    main()





