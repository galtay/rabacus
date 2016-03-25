""" A module supplying containers for groups of ions.  
""" 

class AbsorbingSpecies:
    r""" Photoabsorbing species in primordial chemistry models,
    :math:`\{ {\scriptstyle{\rm HI, \, HeI, \, HeII }} \}`. """ 
    def __init__( self ):
        self.H1 = None
        self.He1 = None
        self.He2 = None

class RecombiningSpecies:
    r""" Recombining species in primoridial chemistry model,
    :math:`\{ {\scriptstyle{\rm HII, \, HeII, \, HeIII }} \}`. """ 
    def __init__( self ):
        self.H2 = None
        self.He2 = None
        self.He3 = None
