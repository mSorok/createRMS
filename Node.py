
class Node:
    """ class representing the RMS node in the general tree """
    """ id : node indentifier, the id to assign once the tree is constructed. It's the final RMSid """
    """ diametre : the diameter for which the tree was computed """"
    """ parent node """
    """ List of child nodes """
    
    
    
    
    
    def __init__(self,idtmp,diametre):
        #constructor
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
        """ setter for the final RMSid """
        self.RMS_id = def_id
        
        
        

