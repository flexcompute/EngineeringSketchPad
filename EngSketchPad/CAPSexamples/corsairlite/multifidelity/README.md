# Multi-fidelity Wing Design Using CAPS Phases

This is an example of Corsairlite and CAPS Phasing used to optimize a wing in phases by evolving aerodynamic and geometric fidelity. 

## Problem Statement
For this demonstration, we consider a design optimization problem originally proposed in Hoburg and Abbeel 2014. The problem seeks to size a UAV flying an outbound cruise segment, followed by a return cruise segment, with a third sprint segment being used to size the engine. The goal is to minimize the total weight of consumed fuel over the three segments.

## Setup
If running the examples in the IDE, one environment variable that will be helpful is the `ESP_PREFIX` environment. This populates the text field when trying to load a pyCAPS script to save time and repetition. For this project, all pyCAPS scripts live in `pyCAPS/`.
```
$> export ESP_PREFIX="pyCAPS/"
```

To run the phasing scripts in the terminal as normal pyCAPS scripts, please run them from the root `multifidelity/` directory so that all the appropriate dependent files can be found.

## File Structure
This directory is organized as follows.
```
.
|--- csm
|--- pyCAPS
|    |---pickle
|--- hoburg
|--- README.md
```
The `` directory contains all the csm and pyCAPS files in their respective directories. The `hoburg/` directory is automatically created when serveESP is run and the first phase has been started. 

Within the `pyCAPS/` folder, the common variables and constraints to all phases are defined in the `defineFormulation.py`. This is imported in each phase as the base formulation, then variables and constraints corresponding to the given phase are appended to the formulation. The `pyCAPS/pickle` directory contains the results of the optimization for each phase. The results of the previous phase are seeded in as the initial guess for the next phase.


## Design Phases
### Phase GPSize
The first phase, `GPSize`, solves the original geometric program formulation proposed in the paper. The values for wing thickness, area, aspect ratio, and taper are solved and updated in the CAPS geometry. The airfoil is constrained to a NACA 24XX.

#### Files for this phase
`csm/naca.csm` and `pyCAPS/GPSize.py`

### Phase MSES
The second phase, `MSES`, increases the aerodynamic model fidelity by replacing the fitted XFOIL data for profile drag with a runtime call to MSES. 

#### Files for this phase
`csm/naca.csm` and `pyCAPS/MSES.py`

### Phase Camber
The third phase, `Camber`, introduces new geometric variables. In addition to those from phase 2, the airfoil camber and camber location are also optimized, albeit constrained to a NACA 4-digit airfoil.

#### Files for this phase
`csm/naca.csm` and `pyCAPS/Camber.py`

### Phase CMConstraint
The optimizer in phase `Camber` exploits the lack of a thickness constraint on the trailing edge of the airfoil, resulting in an unrealistic aft camber and aerodynamic load. To correct for this, we impose a constraint on the moment coefficient that limits the amount of lift concentrated in the thin, aft section of the airfoil. `CMConstraint` starts again from `MSES`, the last successful design phase.

#### Files for this phase
`csm/naca.csm` and `pyCAPS/CMConstraint.py`

### Phase Kulfan2
Phase `Kulfan2` again introduces new geometric variables by expanding the airfoil design space from a NACA 4-digit family of airfoils to a general two-parameter Kulfan airfoil.

#### Files for this phase
`csm/kulfan2.csm` and `pyCAPS/Kulfan2.py`

### Phase FlowTrip
The optimizer in phase `Kulfan2` performs optimally at the design condition, but has poor performance outside the cruise condition. The airfoil has a low $C_{L,max}$ and optimistically assumes the flow stays laminar for most of the chord. For phase `FlowTrip`, we redesign the aircraft with a constraint in MSES that trips the flow at 35% to prevent this. `FlowTrip` starts again from `CMConstraint`, the last successful design phase.

#### Files for this phase
`csm/kulfan2.csm` and `pyCAPS/FlowTrip.py`

### Phase Kulfan4
Phase `Kulfan4` expands the design space again by using a four-parameter Kulfan airfoil.

#### Files for this phase
`csm/kulfan4.csm` and `pyCAPS/Kulfan4.py`

## To Note
- This is a work in progress and is used as a part of a submission for the 2023 AIAA SciTech. 
- Other than phase `GPSize`, the optimizations take on the order of 5 to 20 minutes to solve depending on the computing machine capabilities.
- The phases with Kulfan airfoils (phases 5-7) require the Python `torch` package to work correctly.
- Feel free to email dongjoon@mit.edu with any questions