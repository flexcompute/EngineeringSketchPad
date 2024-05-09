##============================================================##
##                                                            ##
##                    sierraInputGenerator.py                 ##
##                                                            ##
## This function generates an example input file for Sierra   ##
## Mechanics                                                  ##
##                                                            ##
## Inputs:                                                    ##
## *********                                                  ##
##                                                            ##
## Outputs:                                                   ##
## *********                                                  ##
## The output of this function is an input.i file for use in  ##
## Sierra Mechanics.                                          ##
##                                                            ##
## ---------------------------------------------------------- ##
##                                                            ##
## Written by  : Jack Rossetti                                ##
## Date        : 20231004                                     ##
##                                                            ##
## Edits/Notes :                                              ##
## =============================                              ##
## 20231004 JSR |                                             ##
## **************                                             ##
## >> Code has been written                                   ##
##                                                            ##
##============================================================##


## [importPrint]
from __future__ import print_function
## [importPrint]

# Import os module
## [import]
try:
    import os
except:
    print ("Unable to import os module")
    raise SystemError
import argparse
import numpy as np
import time
## [import]

def write_materialDef(fpo, material):

    fpo.write("  begin material %s\n" % material["name"])
    fpo.write("    density = %8.2f\n" % material["rho"])
    fpo.write("    begin parameters for model %s\n" % material["model"])
    fpo.write("      youngs modulus = %8.2e\n" % material["yngmod"])
    fpo.write("      poissons ratio = %8.2f\n" % material["posrat"])
    fpo.write("    end parameters for model %s\n" % material["model"])
    fpo.write("  end material %s\n" % material["name"])
    fpo.write("\n")

####============================================================

def write_elementBlock(fpo, block):

    if(len(block["blockID"]) == 1):
    ## begin if
        fpo.write("    begin parameters for block block_%d\n" % block["blockID"][0])
        fpo.write("      material = %s\n" % block["material"]["name"])
        fpo.write("      model    = %s\n" % block["material"]["model"])
        fpo.write("      section  = %s\n" % block["section"])
        fpo.write("    end parameters for block block_%d\n" % block["blockID"][0])
        fpo.write("\n")
    elif(len(block["blockID"]) != 1):
        for ii in range(0,len(block["blockID"])):
        ## begin for ii
            fpo.write("    begin parameters for block block_%d\n" % block["blockID"][ii])
            fpo.write("      material = %s\n" % block["material"][ii]["name"])
            fpo.write("      model    = %s\n" % block["material"][ii]["model"])
            fpo.write("      section  = %s\n" % block["section"][ii])
            fpo.write("    end parameters for block block_%d\n" % block["blockID"][ii])
            fpo.write("\n")
        ## end for ii
    ## end if

####============================================================

def write_contactDefinition(fpo, friction_model):
    fpo.write("      begin contact definition test\n")
    fpo.write("        CONTACT SURFACE surf17 CONTAINS sideset_17\n")
    fpo.write("        CONTACT SURFACE surf18 CONTAINS sideset_18\n")
    fpo.write("        CONTACT SURFACE surf19 CONTAINS sideset_19\n")
    fpo.write("        CONTACT SURFACE surf20 CONTAINS sideset_20\n")
    fpo.write("        CONTACT SURFACE surf21 CONTAINS sideset_21\n")
    fpo.write("        CONTACT SURFACE surf22 CONTAINS sideset_22\n")
    fpo.write("        CONTACT SURFACE surf23 CONTAINS sideset_23\n")
    fpo.write("        CONTACT SURFACE surf24 CONTAINS sideset_24\n")
    fpo.write("        CONTACT SURFACE surf25 CONTAINS sideset_25\n")
    fpo.write("        CONTACT SURFACE surf26 CONTAINS sideset_26\n")
    fpo.write("        CONTACT SURFACE surf27 CONTAINS sideset_27\n")
    fpo.write("        CONTACT SURFACE surf28 CONTAINS sideset_28\n")
    fpo.write("        CONTACT SURFACE surf29 CONTAINS sideset_29\n")
    fpo.write("        CONTACT SURFACE surf30 CONTAINS sideset_30\n")
    fpo.write("        CONTACT SURFACE surf31 CONTAINS sideset_31\n")
    fpo.write("        CONTACT SURFACE surf32 CONTAINS sideset_32\n")
    fpo.write("        CONTACT SURFACE surf33 CONTAINS sideset_33\n")
    fpo.write("        CONTACT SURFACE surf34 CONTAINS sideset_34\n")
    fpo.write("        CONTACT SURFACE surf35 CONTAINS sideset_35\n")
    fpo.write("        CONTACT SURFACE surf36 CONTAINS sideset_36\n")
    fpo.write("        CONTACT SURFACE surf37 CONTAINS sideset_37\n")
    fpo.write("        CONTACT SURFACE surf38 CONTAINS sideset_38\n")
    fpo.write("        CONTACT SURFACE surf39 CONTAINS sideset_39\n")
    fpo.write("        CONTACT SURFACE surf40 CONTAINS sideset_40\n")
    fpo.write("        CONTACT SURFACE surf41 CONTAINS sideset_41\n")
    fpo.write("        CONTACT SURFACE surf42 CONTAINS sideset_42\n")
    fpo.write("        CONTACT SURFACE surf43 CONTAINS sideset_43\n")
    fpo.write("        CONTACT SURFACE surf44 CONTAINS sideset_44\n")
    fpo.write("        CONTACT SURFACE surf45 CONTAINS sideset_45\n")
    fpo.write("        CONTACT SURFACE surf46 CONTAINS sideset_46\n")
    fpo.write("        CONTACT SURFACE surf47 CONTAINS sideset_47\n")
    fpo.write("        CONTACT SURFACE surf48 CONTAINS sideset_48\n")
    fpo.write("        CONTACT SURFACE surf49 CONTAINS sideset_49\n")
    fpo.write("        CONTACT SURFACE surf50 CONTAINS sideset_50\n")
    fpo.write("        CONTACT SURFACE surf51 CONTAINS sideset_51\n")
    fpo.write("        CONTACT SURFACE surf52 CONTAINS sideset_52\n")
    fpo.write("        CONTACT SURFACE surf53 CONTAINS sideset_53\n")
    fpo.write("        CONTACT SURFACE surf54 CONTAINS sideset_54\n")
    fpo.write("        CONTACT SURFACE surf55 CONTAINS sideset_55\n")
    fpo.write("        CONTACT SURFACE surf56 CONTAINS sideset_56\n")
    fpo.write("        CONTACT SURFACE surf57 CONTAINS sideset_57\n")
    fpo.write("        CONTACT SURFACE surf58 CONTAINS sideset_58\n")
    fpo.write("        CONTACT SURFACE surf59 CONTAINS sideset_59\n")
    fpo.write("        CONTACT SURFACE surf60 CONTAINS sideset_60\n")
    fpo.write("        CONTACT SURFACE surf61 CONTAINS sideset_61\n")
    fpo.write("        CONTACT SURFACE surf62 CONTAINS sideset_62\n")
    fpo.write("        CONTACT SURFACE surf63 CONTAINS sideset_63\n")
    fpo.write("        CONTACT SURFACE surf64 CONTAINS sideset_64\n")
    fpo.write("        CONTACT SURFACE surf65 CONTAINS sideset_65\n")
    fpo.write("        CONTACT SURFACE surf66 CONTAINS sideset_66\n")
    fpo.write("        CONTACT SURFACE surf67 CONTAINS sideset_67\n")
    fpo.write("        CONTACT SURFACE surf68 CONTAINS sideset_68\n")
    fpo.write("        CONTACT SURFACE surf69 CONTAINS sideset_69\n")
    fpo.write("        CONTACT SURFACE surf70 CONTAINS sideset_70\n")
    fpo.write("        CONTACT SURFACE surf71 CONTAINS sideset_71\n")
    fpo.write("        CONTACT SURFACE surf72 CONTAINS sideset_72\n")
    fpo.write("        CONTACT SURFACE surf73 CONTAINS sideset_73\n")
    fpo.write("        CONTACT SURFACE surf74 CONTAINS sideset_74\n")
    fpo.write("        CONTACT SURFACE surf75 CONTAINS sideset_75\n")
    fpo.write("        CONTACT SURFACE surf76 CONTAINS sideset_76\n")
    fpo.write("        CONTACT SURFACE surf77 CONTAINS sideset_77\n")
    fpo.write("        CONTACT SURFACE surf78 CONTAINS sideset_78\n")
    fpo.write("        CONTACT SURFACE surf79 CONTAINS sideset_79\n")
    fpo.write("        CONTACT SURFACE surf80 CONTAINS sideset_80\n")
    fpo.write("        CONTACT SURFACE surf81 CONTAINS sideset_81\n")
    fpo.write("        CONTACT SURFACE surf82 CONTAINS sideset_82\n")
    fpo.write("        CONTACT SURFACE surf83 CONTAINS sideset_83\n")
    fpo.write("        CONTACT SURFACE surf84 CONTAINS sideset_84\n")
    fpo.write("        CONTACT SURFACE surf85 CONTAINS sideset_85\n")
    fpo.write("        CONTACT SURFACE surf86 CONTAINS sideset_86\n")
    fpo.write("        CONTACT SURFACE surf87 CONTAINS sideset_87\n")
    fpo.write("        CONTACT SURFACE surf88 CONTAINS sideset_88\n")
    fpo.write("        CONTACT SURFACE surf01 CONTAINS sideset_1\n")
    fpo.write("        CONTACT SURFACE surf03 CONTAINS sideset_3\n")
    fpo.write("        CONTACT SURFACE surf05 CONTAINS sideset_5\n")
    fpo.write("        CONTACT SURFACE surf06 CONTAINS sideset_6\n")
    fpo.write("        CONTACT SURFACE surf13 CONTAINS sideset_13\n")
    fpo.write("        CONTACT SURFACE surf14 CONTAINS sideset_14\n")
    fpo.write("        SKIN ALL BLOCKS = ON\n")    
    fpo.write("\n");
    fpo.write("        begin constant friction model smooth\n")
    fpo.write("          friction coefficient = 0.100000\n")
    fpo.write("        end constant friction model smooth\n")    
    fpo.write("\n");
    fpo.write("        begin constant friction model rough\n")
    fpo.write("          friction coefficient = 0.600000\n")
    fpo.write("        end constant friction model rough\n")    
    fpo.write("\n");
    fpo.write("        begin interaction defaults\n")
    fpo.write("          friction model  = smooth\n")
    fpo.write("          general contact = on\n")
    fpo.write("        end interaction defaults\n")    
    fpo.write("\n");
    fpo.write("        begin interaction bearingsSides\n")
    fpo.write("          SURFACES = surf18 surf21 surf24 \$ \n")
    fpo.write("                     surf27 surf30 surf33 \$ \n")
    fpo.write("                     surf36 surf39 surf42 \$ \n")
    fpo.write("                     surf45 surf48 surf51 \$ \n")
    fpo.write("                     surf54 surf57 surf60 \$ \n")
    fpo.write("                     surf63 surf66 surf69 \$ \n")
    fpo.write("                     surf72 surf75 surf78 \$ \n")
    fpo.write("                     surf81 surf84 surf87\n")
    fpo.write("          friction model  = frictionless\n")
    fpo.write("        end interaction bearingsSides\n")    
    fpo.write("\n");
    fpo.write("        begin interaction bearingsRaces\n")
    fpo.write("          SURFACES = surf19 surf22 surf25 \$\n")
    fpo.write("                     surf28 surf31 surf34 \$\n")
    fpo.write("                     surf37 surf40 surf43 \$\n")
    fpo.write("                     surf46 surf49 surf52 \$\n")
    fpo.write("                     surf53 surf56 surf59 \$\n")
    fpo.write("                     surf62 surf65 surf68 \$\n")
    fpo.write("                     surf71 surf74 surf77 \$\n")
    fpo.write("                     surf80 surf83 surf86 \$\n")
    fpo.write("                     surf03 surf13\n")
    fpo.write("          friction model  = tied\n")
    fpo.write("        end interaction bearingsRaces\n")    
    fpo.write("\n");
    fpo.write("        begin interaction bearingsInterfaces\n")
    fpo.write("          SURFACES = surf17 surf20 surf23 \$\n")
    fpo.write("                     surf26 surf29 surf32 \$\n")
    fpo.write("                     surf35 surf38 surf41 \$\n")
    fpo.write("                     surf44 surf47 surf50 \$\n")
    fpo.write("                     surf55 surf58 surf61 \$\n")
    fpo.write("                     surf64 surf67 surf70 \$\n")
    fpo.write("                     surf73 surf76 surf79 \$\n")
    fpo.write("                     surf82 surf85 surf88 \$\n")
    fpo.write("                     surf01 surf05 surf06 \$\n")
    fpo.write("                     surf14\n")
    fpo.write("          friction model  = rough\n")
    fpo.write("        end interaction bearingsInterfaces\n")
    fpo.write("\n");
    #fpo.write("        skin all blocks = on\n")
    #fpo.write("\n")
    #fpo.write("        begin constant friction model %s\n" % friction_model["name"])
    #fpo.write("          friction coefficient = %f\n" % friction_model["coeff"])
    #fpo.write("        end constant friction model %s\n" % friction_model["name"])
    #fpo.write("\n")
    #fpo.write("        begin interaction defaults\n")
    #fpo.write("          friction model  = %s\n" % friction_model["name"])
    #fpo.write("          #friction model  = %s\n" % 'tied')
    #fpo.write("          general contact = on\n")
    #fpo.write("        end interaction defaults\n")
    #fpo.write(" \n")
    fpo.write("      end contact definition test\n")
    fpo.write("\n")
    
####============================================================
    
def set_BoundaryConditions(fpo, bcs):
        
    for ibc in bcs:
    ## begin for ibc
        if(ibc["bcType"] == "fixed displacement"):
        ## begin if
            fpo.write("      begin fixed displacement %s\n" % ibc["name"])
            fpo.write("        node set = %s\n"             % ibc["nodeset"])
            fpo.write("        components = %s\n"           % ibc["components"])
            fpo.write("      end fixed displacement %s\n"   % ibc["name"])
            fpo.write("\n")
        elif(ibc["bcType"] == "pressure"):
            fpo.write("      begin pressure %s\n"    % ibc["name"])
            fpo.write("        surface       = %s\n" % ibc["sideset"])
            fpo.write("        read variable = pressure\n")
            fpo.write("        object type   = node\n")
            fpo.write("        scale factor  = %e\n" % ibc["scale"])
            fpo.write("      end pressure %s\n"      % ibc["name"])
            fpo.write("\n")
        ## end if
    ## end for ibc
####============================================================

def sierra_writeInputs(case_dir           , 
                       materials          , 
                       components         , 
                       friction_model     ,
                       boundary_conditions ):

    fpo = open("%s/input.i" % case_dir, "w");
    fpo.write("begin sierra input\n")
    fpo.write("\n")
    fpo.write("  title mutant aerostructural example\n")
    fpo.write("\n")
    fpo.write("  define direction x with vector 1.0 0.0 0.0\n")
    fpo.write("  define direction y with vector 0.0 1.0 0.0\n")
    fpo.write("  define direction z with vector 0.0 0.0 1.0\n")
    fpo.write("\n")
    fpo.write("  begin function gravity_accel\n")
    fpo.write("    type is constant\n")
    fpo.write("    begin values\n")
    fpo.write("      1.0\n")
    fpo.write("    end\n")
    fpo.write("  end function gravity_accel\n")
    fpo.write("\n")
    
    #
    # Write material properties:
    #
    for imat in materials:
    ## begin for imat
        write_materialDef(fpo, imat);
    ## end for imat
    
    fpo.write("  begin solid section ug\n")
    fpo.write("  end solid section ug\n")
    fpo.write("\n")
    fpo.write("  begin finite element model fem\n")
    fpo.write("    database name = tetgenMesh.g\n")
    fpo.write("    database type = exodusII\n")
    fpo.write("\n")
    
    #
    # Write element block information:
    #
    for icomp in components:
    ## begin for icomp
        write_elementBlock(fpo, icomp);
    ## end for icomp

    fpo.write("  end finite element model fem\n")
    fpo.write("\n")
    fpo.write("  begin adagio procedure agio_procedure\n")
    fpo.write("\n")
    fpo.write("    begin time control\n")
    fpo.write("      begin time stepping block p0\n")
    fpo.write("        start time = 0.0\n")
    fpo.write("        begin parameters for adagio region agio_region\n")
    fpo.write("          time increment = 1.0\n")
    fpo.write("        end parameters for adagio region agio_region\n")
    fpo.write("      end time stepping block p0\n")
    fpo.write("      termination time = 1.0\n")
    fpo.write("    end time control\n")
    fpo.write("\n")
    fpo.write("    begin adagio region agio_region\n")
    fpo.write("\n")
    fpo.write("      begin adaptive time stepping\n")
    fpo.write("      end adaptive time stepping\n")
    fpo.write("\n")

    write_contactDefinition(fpo, friction_model);
    
    fpo.write("      use finite element model fem\n")
    fpo.write("\n")

    fpo.write("      begin gravity down\n")
    fpo.write("        function = gravity_accel\n")
    fpo.write("        direction = z\n")
    fpo.write("        gravitational constant = 9.81\n")
    fpo.write("      end gravity down\n")
    fpo.write("\n")

    #fpo.write("      begin gravity gForce\n")
    #fpo.write("        function = gravity_accel\n")
    #fpo.write("        direction = x\n")
    #fpo.write("        gravitational constant = -9.81\n")
    #fpo.write("      end gravity gForce\n")
    #fpo.write("\n")
    
    set_BoundaryConditions(fpo, boundary_conditions);

    fpo.write("      begin results output agio_region_output\n")
    fpo.write("        database name = results.e\n")
    fpo.write("        at time 1.0 increment = 1.0\n")
    fpo.write("        nodal variables   = displacement\n")
    fpo.write("        nodal variables   = force_external\n")
    fpo.write("        nodal variables   = force_internal\n")
    fpo.write("        nodal variables   = force_inertial\n")
    fpo.write("        nodal variables   = reaction\n")
    fpo.write("        nodal variables   = mass\n")
    fpo.write("        nodal variables   = residual\n")
    fpo.write("        element variables = log_strain\n")
    fpo.write("        element variables = principal_stresses\n")
    fpo.write("        element variables = min_principal_stress\n")
    fpo.write("        element variables = max_principal_stress\n")
    fpo.write("        element variables = stress\n")
    fpo.write("        element variables = von_mises\n")
    fpo.write("        element variables = strain_energy\n")
    fpo.write("        element variables = strain_energy_density\n")
    fpo.write("        element variables = element_mass\n")
    fpo.write("        element variables = element_shape\n")
    fpo.write("        element variables = centroid\n")
    fpo.write("        element variables = volume\n")
    fpo.write("      end results output agio_region_output\n")
    fpo.write("\n")
    fpo.write("      begin solver\n")
    fpo.write("        begin control contact\n")
    fpo.write("          target relative residual = 5e-4\n")
    fpo.write("          maximum iterations       = 150\n")
    fpo.write("          level = 1\n")
    fpo.write("        end control contact\n")
    fpo.write("\n")
    fpo.write("        begin loadstep predictor\n")
    fpo.write("          type = scale_factor\n")
    fpo.write("          scale factor = 0.0 0.0\n")
    fpo.write("        end loadstep predictor\n")
    fpo.write("        level 1 predictor = none\n")
    fpo.write("\n")
    fpo.write("        begin cg\n")
    fpo.write("          reference = Belytschko\n")
    fpo.write("          acceptable relative residual = 1.0e10\n")
    fpo.write("          target relative residual     = 5e-5\n")
    fpo.write("          maximum iterations           = 250\n")
    fpo.write("          begin full tangent preconditioner\n")
    fpo.write("            minimum smoothing iterations = 100\n")
    fpo.write("          end full tangent preconditioner\n")
    fpo.write("        end cg\n")
    fpo.write("      end solver\n")
    fpo.write("\n")
    fpo.write("    end adagio region agio_region\n")
    fpo.write("  end adagio procedure agio_procedure\n")
    fpo.write("end sierra input\n")

    fpo.close();
    
####============================================================

####
## begin main function/script:
####

start = time.time();

## [argparse]
# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Sierra Input File Example',
                                 prog = 'sierraInputGenerator.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

caseDir       = '.';

#
#  Create element block definitions:
##
## Define materials
##
materials  = [];

aluminum    ={ "name"  : "aluminum",
               "model" : "elastic",
               "rho"   : 2810,
               "yngmod": 71.7e9,
               "posrat": 0.3  };

steel      ={ "name"  : "steel",
              "model" : "elastic",
              "rho"   : 7810,
              "yngmod": 200e9,
              "posrat": 0.3  };

unobtain   ={ "name"  : "unobtainium",
              "model" : "elastic",
              "rho"   : 99999,
              "yngmod": 500e9,
              "posrat": 1.0  };
    
materials.append(aluminum);
materials.append(steel);
materials.append(unobtain);

components = [];
component  = {"name"     : "component_1",
              "blockID"  : [1],
              "sideSetID": [1,2],
              "material" : materials[0],
              "section"  : "ug"  }
components.append(component);

component  = {"name"     : "component_2",
              "blockID"  : [2],
              "sideSetID": [3,4,5],
              "material" : materials[0],
              "section"  : "ug"  }
components.append(component);

component  = {"name"     : "component_3",
              "blockID"  : [3],
              "sideSetID": [6,7],
              "material" : materials[0],
              "section"  : "ug"  }
components.append(component);

component  = {"name"     : "component_4",
              "blockID"  : [4,5,6],
              "sideSetID": [8,9,10,11,12,13,14,15,16],
              "material" : [materials[0], materials[2], materials[0]],
              "section"  : ["ug", "ug", "ug"]  }
components.append(component);
#print("material: %s" % component["material"][0]["name"])
nBlocks = 0;
nSide   = 0;
for icomp in components:
## begin for icomp
    #print(len(icomp["blockID"]));
    for ii in range(0,len(icomp["blockID"])):
    ## begin for ii
        print("Component -> %s" % icomp["name"]);
        print("\tblockID : %d" % icomp["blockID"][ii]);
        if(len(icomp["blockID"]) == 1):
            print(icomp["material"]["name"]);
        elif(len(icomp["blockID"]) != 1):
            print(icomp["material"][ii]["name"]);
        ## end if
    ## end for ii

    nBlocks = nBlocks + len(icomp["blockID"]);
    nSide   = nSide   + len(icomp["sideSetID"]);
    #print(nBlocks);
## end for icomp

nBase = len(components);
nComp = nBlocks;

friction_model = {"name" : "smooth",
                  "coeff": 0.1}

fixedBC = { "bcType"     : "fixed displacement",
            "name"       : "fixedEnd_BC",
            "nodeset"    : "nodelist_11",
            "components" : "x y z"}

pressBC = { "bcType" : "pressure",
            "name"   : "aero_pressure",
            "sideset": "sideset_8 sideset_9 sideset_10",
            "scale"  :  1.0 }

boundary_conditions = [fixedBC, pressBC] 

sierra_writeInputs(caseDir            , 
                   materials          , 
                   components         , 
                   friction_model     ,
                   boundary_conditions )

end = time.time();
tot = end-start;
hh  = int((tot-tot%3600)/3600);
mm  = int((tot-tot%60)/60 - 60*hh);
ss  = int(tot - tot%1 - 60*mm - 3600*hh);

print("Sierra Input Writer Timing:");
print("    Time Elapsed:\t %f" % tot);
print("    Time Elapsed:\t %02d hh:\t%02d mm:\t%02d ss" % (hh,mm,ss))

exit()
