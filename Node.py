
class Node:
    """ classe representant le noeud-RMS_id dans l'arborescence generale """
    """ identifiant du noeud = id a assigner une fois l'arbre a ete contruit. il s'agit du RMS_id definitif"""
    """ diametre pour lequel il a ete calcule"""
    """ Noeud parent"""
    """ Liste de noeuds enfants """
    
    
    
    
    
    def __init__(self,idtmp,diametre):
        #constructeur
        self.id_rmsf = idtmp
        self.parent = None
        self.diametre = diametre
        self.RMS_id = ""
        self.children = []
        
        
    def defineParent(self,parent):
        self.parent = parent
        
    def addChild(self,child):
        self.children.append(child)
   
    
    def setRmsId(self,def_id):
        """setteur pour l'id du RMS definitif """
        self.RMS_id = def_id
        
        
        

