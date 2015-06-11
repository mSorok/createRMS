#!/usr/bin/python -tt

import os
import os.path
import re
import time
import sys
import commands




def main():
    
    
    #par defaut, on a une nouvelle reaction
    newReactionFlag=1
    
    leftCompounds = []
    rightCompounds = []
    
    lastLcompound = ''
    lastRcompound = ''
    
    frame = ''
    
    reactionDirection = ''
    
    
    leftCflag = 0
    rightCflag = 0
    
    
    reaFile = open('/env/cns/proj/agc/msorokina/RMS/createRMS/reactions.dat','r')
    
    reaFileTab = open('/env/cns/proj/agc/msorokina/RMS/createRMS/reactions_tab','w')
    
    for line in reaFile:
        if not line.startswith('#'):
            #on ne lit pas les commentaires
            
            if line.startswith('//'):
                
                #on a termine une reaction - on print tout et on passe tous les flag a 0
                newReactionFlag = 0
                leftC = '$'.join(leftCompounds)
                rightC = '$'.join(rightCompounds)
                
                
                leftCflag=0
                rightCflag=0
                
                reaFileTab.write(frame[:-1]+"\t"+leftC+"\t"+reactionDirection+"\t"+rightC+"\n")
                
                
            # frame de la reaction
            if line.startswith('UNIQUE-ID'):
                tabf = line.split(' - ')
                frame=tabf[1]
                newReactionFlag = 1
                
                reactionDirection = ''
                leftCflag=0
                rightCflag=0
                                
                
                
            #composants de gauche  
            if line.startswith('LEFT'):
                if leftCflag==0:
                    tabc = line.split(' - ')
                    leftCflag = 1
                    leftCompounds = []
                    c = tabc[1]
                    leftCompounds.append(c[:-1])
                    lastLcompound = c[:-1]
            
                elif leftCflag==1:
                    tabc = line.split(' - ')
                    c = tabc[1]
                    leftCompounds.append(c[:-1])
                    lastLcompound = c[:-1]
                
            #coefficients de stoechiometrie de gauche
            if line.startswith('^COEFFICIENT') and leftCflag == 1 and rightCflag == 0:
                tabc = line.split(' - ')
                leftCompounds.remove(lastLcompound)
                c = tabc[1]
                lastLcompound = lastLcompound+'*'+c[:-1]
                leftCompounds.append(lastLcompound)
                
                
                
                
            #composants de droite   
            if line.startswith('RIGHT'): 
                if rightCflag==0:
                    tabc = line.split(' - ')
                    rightCflag = 1
                    rightCompounds = []
                    c = tabc[1]
                    rightCompounds.append(c[:-1])
                    lastRcompound = c[:-1]
                    leftCflag = 0 #on le met bien a 0 pour pas que ca puisse gener au niveau des coefficients
                elif rightCflag==1:
                    tabc = line.split(' - ')
                    c = tabc[1]
                    rightCompounds.append(c[:-1])
                    lastRcompound = c[:-1]
                
            #coefficients de stoechiometrie de droite
            if line.startswith('^COEFFICIENT') and leftCflag == 0 and rightCflag == 1:
                tabc = line.split(' - ')
                c = tabc[1]
                rightCompounds.remove(lastRcompound)
                lastRcompound = lastRcompound+'*'+c[:-1]
                rightCompounds.append(lastRcompound)
                
            
            #direction de la reaction   
            if line.startswith('REACTION-DIRECTION'):
                tabD = line.split(' - ')
                d = tabD[1]
                reactionDirection = d[:-1]
            
    
    

    
    reaFileTab.close()
    reaFile.close()

 
if __name__ == "__main__":
    main()





