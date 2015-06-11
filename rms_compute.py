
#!/usr/bin/python -tt

import os
import os.path
import re
import time
import sys
import commands




def main():

    #argument precisant le type de molsig
    #au choix : scan, sscan, fscan
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
    #transforme le fichier .dat en un fichier tabule avec une seule reaction par ligne
    commands.getstatusoutput('python reactions_normalizer.py')    
           

    
    reaFile = open('/env/cns/proj/agc/msorokina/RMS/createRMS/reactions_tab','r')
    
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
        
        
        #si tous les composants sont connus
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
            
            
            #curation des substrats
            
            halfmolS = False
            
            for st in substratesTab:
                
                st = st.replace('\n','')
                st = st.replace('|','')
                elem = st.split('*')
                c = elem[0] #le composant en soi
                c = c.replace('\n','')

                
                coef = 1
                if len(elem)> 1 :
                    coef = elem[1]
                if str(coef).find('N')>-1:
                    coef = 3    
                if coef=='m+q':
                    coef = 3 #valeur generique        
                if coef=='n':
                    coef = 3 #valeur generique
                if coef == '2n':
                    coef = 6
                if coef == '3n':
                    coef = 9
                if coef == '4n':
                    coef = 12
                if coef == '(n+1)':
                    coef = 4
                
                    
                    
                if coef == 0.5:  #cas de demi-molecules d'oxygene
                    halfmolS =True
                    
                #coef = int(coef)

                #ajustements connus sur les noms d'elements chimiques
                if c == 'NAD-P-OR-NOP':
                    c='NADP'
                c = c.replace('|','')
                

                
                if c != 'E-':
                    substratesDic[c]=coef
                
            
            
            
            #curation des produits
            
            halfmolP = False
            
            for pt in productsTab:
                pt = pt.replace('\n','')
                pt = pt.replace('|','')
                elem = pt.split('*')
                c = elem[0] #le composant en soi
                c = c.replace('\n','')
                
                
                coef = 1
                if len(elem)> 1 :
                    coef = elem[1]
                if str(coef).find('N')>-1:
                    coef = 3    
                if coef=='m+q':
                    coef = 3 #valeur generique        
                if coef=='n':
                    coef = 3 #valeur generique
                if coef == '2n':
                    coef = 6
                if coef == '3n':
                    coef = 9
                if coef == '4n':
                    coef = 12
                if coef == '(n+1)':
                    coef = 4
                
                    
                if coef == 0.5: #cas de demi-molecules d'oxygene
                    halfmolP =True
                    
                #coef = int(coef)
                    
                #ajustements connus sur les noms d'elements chimiques
                if c == 'NAD-P-OR-NOP':
                    c='NADP'
                    
                c = c.replace('|','')
                
                if c != 'E-':
                    productsDic[c]=coef
                
            #regularisation des coefficients si presence de demi-molecules -> tout est multiplie par 2
            if halfmolS == True or halfmolP == True:
                for c in substratesDic.keys():
                    substratesDic[c] = substratesDic[c]*2
                    substratesDic[c] = int(substratesDic[c])
                for c in productsDic.keys():
                    productsDic[c] = productsDic[c]*2
                    productsDic[c] = int(productsDic[c])
                
            
            #determination si tous les fichiers des composants existent sans le repertoire sscan0 (par exemple)
            allFilesExist = True
            for k in substratesDic.keys():
                if not os.path.exists('/env/cns/proj/agc/msorokina/RMS/createRMS/mol-sig-results/'+molsigdir+'10/'+k+'.'+molSigType+'10'):
                    allFilesExist = False
                   

            for k in productsDic.keys():
                if not os.path.exists('/env/cns/proj/agc/msorokina/RMS/createRMS/mol-sig-results/'+molsigdir+'10/'+k+'.'+molSigType+'10'):
                    allFilesExist= False
            
            
            
            if allFilesExist and len(substratesDic.keys()) != 0 and len(productsDic.keys()) != 0:
                for i in rayons:
                    
                    #path ou se trouvent tous les fichiers molsig a ce rayon
                    source = '/env/cns/proj/agc/msorokina/RMS/createRMS/mol-sig-results/'+molsigdir+str(i)
                    substratesSignature = {}
                    productsSignature = {}
                
                
                    #recuperation des molsig pour les substrats
                    for c in substratesDic.keys():
                        filename = source+'/'+c+'.'+molSigType+str(i)
                        try:
                            f = open(filename,"r")
                            
                            for line in f:
                                match = re.search(r"^0.0.*",line)
                                if not match :
                                    
                                    l = line.split(' ') #en 0 on a donc le chiffre et en 1 les atomes
                                    motif = l[1]
                                    motif = motif[:-1]
                                    
                                    #enleve les charges
                                    motif = motif.replace('-1','')
                                    motif = motif.replace('+1','')
                                    
                                    occurency = float(l[0])*float(substratesDic[c])
                                    
                                    if motif in substratesSignature.keys():
                                        substratesSignature[motif] = substratesSignature[motif]+occurency
                                    else:
                                        substratesSignature[motif] = occurency
                                        
                        except IOError:
                            sys.stderr.write(c+"\n")
            
                        
                    
                    
                    
                    #recuperation des molsig pour les produits
                    for c in productsDic.keys():
                        filename = source+'/'+c+'.'+molSigType+str(i)
                        try:
                            f = open(filename,"r")
                            
                            for line in f:
                                match = re.search(r"^0.0.*",line)
                                if not match :
                                    
                                    l = line.split(' ') #en 0 on a donc le chiffre et en 1 les atomes
                                    motif = l[1]
                                    motif = motif[:-1]
                                    
                                    #enleve les charges
                                    motif = motif.replace('-1','')
                                    motif = motif.replace('+1','')
                                    
                                    occurency = float(l[0])*float(productsDic[c])
                                    
                                    if motif in productsSignature.keys():
                                        productsSignature[motif] = productsSignature[motif]+occurency
                                    else:
                                        productsSignature[motif] = occurency
                                        
                        except IOError:
                            sys.stderr.write(c+"\n")
                    
                    
                    
                    #calcul des RMS
                    # maintenant il faut soustraire : products-substrates
                    reactionSignature = {}
                    #3 cas pour les sous parties:
                    # -presence dans les 2 dict
                    # -present dans produts et absent dans substracts
                    # -present dans substrates et absent dans products
                    
                    for p in productsSignature.keys():
                        if p in substratesSignature.keys():
                            #premier cas : 
                            reactionSignature[p] =productsSignature[p]-substratesSignature[p]
                        else:
                            #deuxieme cas :
                            reactionSignature[p] =productsSignature[p]
                
                    #3e cas :
                    for p in substratesSignature.keys():
                        if p not in productsSignature.keys():
                            reactionSignature[p]= 0 - substratesSignature[p]
                
                    #maintenant in faut ordonner les cles de reactionSignature pour qu'elles soient dans le meme ordre meme pour des
                    #reactions differentes (ce qu'on cherche en fait)

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





