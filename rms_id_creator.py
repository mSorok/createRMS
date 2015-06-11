
#!/usr/bin/python -tt

import os
import os.path
import re
import time
import sys
import commands
from Node import Node




def main():
    
    #ouverture du fichier MR_RMSf_chain.txt contenant les identifiants de reactions et la chaine correspondante de RMSf, ordonnee
    if sys.argv[1] != '':
        fileName = sys.argv[1]
    else:
        fileName = "MR_RMSf_chain.txt"
    
    f = open(fileName,"r")
    
    rms_rmsid_list = {} #cle=rmsf, valeur= rms_id
    
    rms_chains = {} #cle:MR_id, valeur:tableau de rmsf
    
    root = None
    
    
    #pour chaque ligne dans ce fichier
    for line in f:
        #les lidgnes sont organisees de facon suivante:
        #MR_id    chaine de RMSf separes de '$' et tries dans l'ordre du 
        lineTab = line.split("\t") #0 : MR_id ; 1 : RMSf chain
        
        MR_id = lineTab[0]
        rmsfstring = lineTab[1]
        rmsfstring = rmsfstring.replace("\n","")
        
        
        rmsf = rmsfstring.split('$') #l'indice du rmsf dans ce tableau correspond au diametre du rmsf
        rms_chains[MR_id] = rmsf
        
        #creaction de la racine de l'arbre si elle n'existe pas deja
        if root==None:
            root = Node(rmsf[0],0)
        
        lastParent = root #on definit le dernier parent
        
        for i in range(1,len(rmsf)):
            #i est le diametre
            
            currentNode = None
            
            found = False
            for n in lastParent.children:
                if n.id_rmsf == rmsf[i]:
                    found = True
                    currentNode = n
                    
            
            if found==False:
                #creation du noeud
                currentNode = Node(rmsf[i],i)
                currentNode.defineParent(lastParent)
                lastParent.addChild(currentNode)
            
            lastParent = currentNode
            
            
    #l'arbre est cree, il est represente par root
    
    #definition des ids definitifs
    
    
    root.setRmsId("RMS-0")
    parcours(root)
    
   
    
    
    
    
            
        
  
        
        
def parcours(tree):
    if tree != None:
        #on considere que le RMS_id de tree est deja defini
        count = 1
        for n in tree.children:
            nid = tree.RMS_id+'.'+str(count)
            n.setRmsId(nid)
            parcours(n)
            count+=1
            

    print tree.id_rmsf+"\t"+str(tree.diametre)+"\t"+tree.RMS_id   


 
if __name__ == "__main__":
    main()



