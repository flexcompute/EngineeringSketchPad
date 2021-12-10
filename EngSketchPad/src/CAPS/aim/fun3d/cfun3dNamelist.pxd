# This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

cdef extern from "cfdTypes.h":
    
    cdef enum cfdSurfaceTypeEnum:
        UnknownBoundary, Inviscid, Viscous, Farfield, Extrapolate, Freestream, BackPressure, 
        Symmetry, SubsonicInflow, SubsonicOutflow, MassflowIn, MassflowOut,
        FixedInflow, FixedOutflow, MachOutflow

    ctypedef struct cfdSurfaceStruct:
        char   *name
    
        cfdSurfaceTypeEnum surfaceType # "Global" boundary condition types
    
        int    bcID                 # ID of boundary
    
        # Wall specific properties
        int    wallTemperatureFlag  # Temperature flag
        double wallTemperature      # Temperature value -1 = adiabatic  >0 = isothermal
        double wallHeatFlux         # Wall heat flux. to use Temperature flag should be true and wallTemperature < 0
    
        # Symmetry plane
        int symmetryPlane        # Symmetry flag / plane
    
        # Stagnation quantities
        double totalPressure    # Total pressure
        double totalTemperature # Total temperature
        double totalDensity     # Total density
    
        # Static quantities
        double staticPressure   # Static pressure
        double staticTemperature# Static temperature
        double staticDensity    # Static temperature
    
        # Velocity components
        double uVelocity  # x-component of velocity
        double vVelocity  # y-component of velocity
        double wVelocity  # z-component of velocity
        double machNumber # Mach number
    
        # Massflow
        double massflow # Mass flow through a boundary


    # Structure to hold CFD boundary condition information
    ctypedef struct cfdBoundaryConditionStruct:
        char             *name   # Name of BCsStruct
        int              numSurfaceProp # Number of unique BC ids
        cfdSurfaceStruct *surfaceProp   # Surface properties for each bc - length of numSurfaceProp