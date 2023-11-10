from parcels import (JITParticle, Variable)
from operator import attrgetter
import numpy as np


### Forward
class Iceberg_FW(JITParticle):
    # Ocean - at surface
    tempS = Variable('tempS', initial=0, to_write=False)   # temperature [°C]
    uveloS = Variable('uveloS', initial=0, to_write=False) # u-velocity [m/s]
    vveloS = Variable('vveloS', initial=0, to_write=False) # v-velocity [m/s]
    
    # Wind - at surface
    uvela = Variable('uvela', initial=0, to_write=False)   # u-velocity [m/s]
    vvela = Variable('vvela', initial=0, to_write=False)   # v-velocity [m/s]

    # Ocean - at iceberg base
    tempB = Variable('tempB', initial=0, to_write=False)   # temperature [°C]
    uveloB = Variable('uveloB', initial=0, to_write=False) # u-velocity [m/s]
    vveloB = Variable('vveloB', initial=0, to_write=False) # v-velocity [m/s]
    
    # Ocean - depth-integrated along iceberg draft
    tempD = Variable('tempD', initial=0, to_write=False)   # temperature [°C]
    uveloD = Variable('uveloD', initial=0, to_write=False) # u-velocity [m/s]
    vveloD = Variable('vveloD', initial=0, to_write=False) # v-velocity [m/s]
    
    # Iceberg properties
    L = Variable('L', initial=0, to_write=False)          # Iceberg length [m]
    W = Variable('W', initial=0, to_write=False)          # Iceberg width [m]
    T = Variable('T', initial=0, to_write=False)          # Iceberg thickness [m]
    M = Variable('M', initial=0)                          # Iceberg mass [kg]
    
    prev_L = Variable('prev_L', to_write=False, initial=attrgetter('L'))
    prev_W = Variable('prev_W', to_write=False, initial=attrgetter('W'))
    prev_T = Variable('prev_T', to_write=False, initial=attrgetter('T'))
    
    # Melt terms
    Mvr = Variable('Mvr', initial=0)                      # buoyant convection [m/s]
    Mbr = Variable('Mbr', initial=0)                      # basal melt [m/s]
    Mer = Variable('Mer', initial=0)                      # wave erosion [m/s]
      
    prev_Mvr = Variable('prev_Mvr', to_write=False, initial=attrgetter('Mvr'))
    prev_Mbr = Variable('prev_Mbr', to_write=False, initial=attrgetter('Mbr'))
    prev_Mer = Variable('prev_Mer', to_write=False, initial=attrgetter('Mer'))
    
    # Others
    bath = Variable('bath', initial=0, to_write=False)    # bathymetry [m]
    
    distance = Variable('distance', initial=1000000)      # distance to ODP Site 696
    prev_dist = Variable('prev_dist', to_write=False, initial=attrgetter('distance'))
    
    check = Variable('check', to_write=False, initial=-1) # check iceberg existance
    prev_check = Variable('prev_check', to_write=False, initial=attrgetter('check'))
    
    
    
### Backward
class Iceberg_BW(JITParticle):
    # Ocean - at surface
    tempS = Variable('tempS', initial=0, to_write=False)   # temperature [°C]
    uveloS = Variable('uveloS', initial=0, to_write=False) # u-velocity [m/s]
    vveloS = Variable('vveloS', initial=0, to_write=False) # v-velocity [m/s]
    
    # Wind - at surface
    uvela = Variable('uvela', initial=0, to_write=False)   # u-velocity [m/s]
    vvela = Variable('vvela', initial=0, to_write=False)   # v-velocity [m/s]

    # Ocean - at iceberg base
    tempB = Variable('tempB', initial=0, to_write=False)   # temperature [°C]
    uveloB = Variable('uveloB', initial=0, to_write=False) # u-velocity [m/s]
    vveloB = Variable('vveloB', initial=0, to_write=False) # v-velocity [m/s]
    
    # Ocean - depth-integrated along iceberg draft
    tempD = Variable('tempD', initial=0, to_write=False)   # temperature [°C]
    uveloD = Variable('uveloD', initial=0, to_write=False) # u-velocity [m/s]
    vveloD = Variable('vveloD', initial=0, to_write=False) # v-velocity [m/s]
    
    # Iceberg properties
    L = Variable('L', initial=0, to_write=False)          # Iceberg length [m]
    W = Variable('W', initial=0, to_write=False)          # Iceberg width [m]
    T = Variable('T', initial=0, to_write=False)          # Iceberg thickness [m]
    M = Variable('M', initial=0)                          # Iceberg mass [kg]
    
    prev_L = Variable('prev_L', to_write=False, initial=attrgetter('L'))
    prev_W = Variable('prev_W', to_write=False, initial=attrgetter('W'))
    prev_T = Variable('prev_T', to_write=False, initial=attrgetter('T'))
    
    # Melt terms
    Mvr = Variable('Mvr', initial=0)                      # buoyant convection [m/s]
    Mbr = Variable('Mbr', initial=0)                      # basal melt [m/s]
    Mer = Variable('Mer', initial=0)                      # wave erosion [m/s]
      
    prev_Mvr = Variable('prev_Mvr', to_write=False, initial=attrgetter('Mvr'))
    prev_Mbr = Variable('prev_Mbr', to_write=False, initial=attrgetter('Mbr'))
    prev_Mer = Variable('prev_Mer', to_write=False, initial=attrgetter('Mer'))
    
    # Others
    bath = Variable('bath', initial=0, to_write=False)    # bathymetry [m]
    reg = Variable('reg', initial=0, dtype=np.int32)      # region Carter
    
    check = Variable('check', to_write=False, initial=-1) # check iceberg existance
    prev_check = Variable('prev_check', to_write=False, initial=attrgetter('check')) 



### Forward daily
class Iceberg_FW_day(JITParticle):
    # Ocean - at surface
    tempS = Variable('tempS', initial=0)   # temperature [°C]
    uveloS = Variable('uveloS', initial=0) # u-velocity [m/s]
    vveloS = Variable('vveloS', initial=0) # v-velocity [m/s]
    
    # Wind - at surface
    uvela = Variable('uvela', initial=0)   # u-velocity [m/s]
    vvela = Variable('vvela', initial=0)   # v-velocity [m/s]

    # Ocean - at iceberg base
    tempB = Variable('tempB', initial=0)   # temperature [°C]
    uveloB = Variable('uveloB', initial=0) # u-velocity [m/s]
    vveloB = Variable('vveloB', initial=0) # v-velocity [m/s]
    
    # Ocean - depth-integrated along iceberg draft
    tempD = Variable('tempD', initial=0)   # temperature [°C]
    uveloD = Variable('uveloD', initial=0) # u-velocity [m/s]
    vveloD = Variable('vveloD', initial=0) # v-velocity [m/s]
    
    # Iceberg properties
    L = Variable('L', initial=0, to_write=False)          # Iceberg length [m]
    W = Variable('W', initial=0, to_write=False)          # Iceberg width [m]
    T = Variable('T', initial=0, to_write=False)          # Iceberg thickness [m]
    M = Variable('M', initial=0)                          # Iceberg mass [kg]
    
    prev_L = Variable('prev_L', to_write=False, initial=attrgetter('L'))
    prev_W = Variable('prev_W', to_write=False, initial=attrgetter('W'))
    prev_T = Variable('prev_T', to_write=False, initial=attrgetter('T'))
    
    # Melt terms
    Mvr = Variable('Mvr', initial=0)                      # buoyant convection [m/s]
    Mbr = Variable('Mbr', initial=0)                      # basal melt [m/s]
    Mer = Variable('Mer', initial=0)                      # wave erosion [m/s]
      
    prev_Mvr = Variable('prev_Mvr', to_write=False, initial=attrgetter('Mvr'))
    prev_Mbr = Variable('prev_Mbr', to_write=False, initial=attrgetter('Mbr'))
    prev_Mer = Variable('prev_Mer', to_write=False, initial=attrgetter('Mer'))
    
    # Others
    bath = Variable('bath', initial=0, to_write=False)    # bathymetry [m]
    
    distance = Variable('distance', initial=1000000)      # distance to ODP Site 696
    prev_dist = Variable('prev_dist', to_write=False, initial=attrgetter('distance'))
    
    check = Variable('check', to_write=False, initial=-1) # check iceberg existance
    prev_check = Variable('prev_check', to_write=False, initial=attrgetter('check'))



### Surface
class Iceberg_FW_Surf(JITParticle):
    # Ocean - at surface
    temp = Variable('temp', initial=0, to_write=False)   # temperature [°C]
    uvelo = Variable('uvelo', initial=0, to_write=False) # u-velocity [m/s]
    vvelo = Variable('vvelo', initial=0, to_write=False) # v-velocity [m/s]
    
    # Wind - at surface
    uvela = Variable('uvela', initial=0, to_write=False)   # u-velocity [m/s]
    vvela = Variable('vvela', initial=0, to_write=False)   # v-velocity [m/s]
    
    # Iceberg properties
    L = Variable('L', initial=0, to_write=False)          # Iceberg length [m]
    W = Variable('W', initial=0, to_write=False)          # Iceberg width [m]
    T = Variable('T', initial=0, to_write=False)          # Iceberg thickness [m]
    M = Variable('M', initial=0)                          # Iceberg mass [kg]
    
    prev_L = Variable('prev_L', to_write=False, initial=attrgetter('L'))
    prev_W = Variable('prev_W', to_write=False, initial=attrgetter('W'))
    prev_T = Variable('prev_T', to_write=False, initial=attrgetter('T'))
    
    # Melt terms
    Mvr = Variable('Mvr', initial=0)                      # buoyant convection [m/s]
    Mbr = Variable('Mbr', initial=0)                      # basal melt [m/s]
    Mer = Variable('Mer', initial=0)                      # wave erosion [m/s]
      
    prev_Mvr = Variable('prev_Mvr', to_write=False, initial=attrgetter('Mvr'))
    prev_Mbr = Variable('prev_Mbr', to_write=False, initial=attrgetter('Mbr'))
    prev_Mer = Variable('prev_Mer', to_write=False, initial=attrgetter('Mer'))
    
    # Others
    bath = Variable('bath', initial=0, to_write=False)    # bathymetry [m]
    
    distance = Variable('distance', initial=1000000)      # distance to ODP Site 696
    prev_dist = Variable('prev_dist', to_write=False, initial=attrgetter('distance'))
    
    check = Variable('check', to_write=False, initial=-1) # check iceberg existance
    prev_check = Variable('prev_check', to_write=False, initial=attrgetter('check'))



### Backward
class Iceberg_BW_Surf(JITParticle):
    # Ocean - at surface
    temp = Variable('temp', initial=0, to_write=False)     # temperature [°C]
    uvelo = Variable('uvelo', initial=0, to_write=False)   # u-velocity [m/s]
    vvelo = Variable('vvelo', initial=0, to_write=False)   # v-velocity [m/s]
    
    # Wind - at surface
    uvela = Variable('uvela', initial=0, to_write=False)   # u-velocity [m/s]
    vvela = Variable('vvela', initial=0, to_write=False)   # v-velocity [m/s]
    
    # Iceberg properties
    L = Variable('L', initial=0, to_write=False)           # Iceberg length [m]
    W = Variable('W', initial=0, to_write=False)           # Iceberg width [m]
    T = Variable('T', initial=0, to_write=False)           # Iceberg thickness [m]
    M = Variable('M', initial=0)                           # Iceberg mass [kg]
    
    prev_L = Variable('prev_L', to_write=False, initial=attrgetter('L'))
    prev_W = Variable('prev_W', to_write=False, initial=attrgetter('W'))
    prev_T = Variable('prev_T', to_write=False, initial=attrgetter('T'))
    
    # Melt terms
    Mvr = Variable('Mvr', initial=0)                       # buoyant convection [m/s]
    Mbr = Variable('Mbr', initial=0)                       # basal melt [m/s]
    Mer = Variable('Mer', initial=0)                       # wave erosion [m/s]
      
    prev_Mvr = Variable('prev_Mvr', to_write=False, initial=attrgetter('Mvr'))
    prev_Mbr = Variable('prev_Mbr', to_write=False, initial=attrgetter('Mbr'))
    prev_Mer = Variable('prev_Mer', to_write=False, initial=attrgetter('Mer'))
    
    # Others
    bath = Variable('bath', initial=0, to_write=False)     # bathymetry [m]
    reg = Variable('reg', initial=0, dtype=np.int32)       # region Carter
    
    check = Variable('check', to_write=False, initial=-1)  # check iceberg existance
    prev_check = Variable('prev_check', to_write=False, initial=attrgetter('check')) 