
#!/usr/bin/python -tt

import os
import os.path
import re
import time
import sys
import commands
from Node import Node




def main():
    
    # openinf of the MR_RMSf_chain.txt file containing the identifiers of tall the reactions and the correspondinf RMSf chain
    if sys.argv[1] != '':
        fileName = sys.argv[1]
    else:
        fileName = "MR_RMSf_chain.txt"
    
    f = open(fileName,"r")
    
    rms_rmsid_list = {} #key=rmsf, value= rms_id
    
    rms_chains = {} #key:MR_id, value: list of RMSf
    
    root = None
    
    

    for line in f:

        lineTab = line.split("\t") #0 : MR_id ; 1 : RMSf chain
        
        MR_id = lineTab[0]
        rmsfstring = lineTab[1]
        rmsfstring = rmsfstring.replace("\n","")
        
        
        rmsf = rmsfstring.split('$')
        rms_chains[MR_id] = rmsf
        
        # creation of the tree root if desn't exist yet
        if root==None:
            root = Node(rmsf[0],0)
        
        lastParent = root #last parent definition
        
        for i in range(1,len(rmsf)):
            #i is the diameter
            
            currentNode = None
            
            found = False
            for n in lastParent.children:
                if n.id_rmsf == rmsf[i]:
                    found = True
                    currentNode = n
                    
            
            if found==False:
                #node creation
                currentNode = Node(rmsf[i],i)
                currentNode.defineParent(lastParent)
                lastParent.addChild(currentNode)
            
            lastParent = currentNode
            
            
    #the tree is created, it is represented by root
    
    #definition of the final ids
    
    
    root.setRmsId("RMS-0")
    parcours(root)
    
   
    
    
    
    
            
        
  
        
        
def parcours(tree):
    if tree != None:
        count = 1
        for n in tree.children:
            nid = tree.RMS_id+'.'+str(count)
            n.setRmsId(nid)
            parcours(n)
            count+=1
            

    print tree.id_rmsf+"\t"+str(tree.diametre)+"\t"+tree.RMS_id   


 
if __name__ == "__main__":
    main()



