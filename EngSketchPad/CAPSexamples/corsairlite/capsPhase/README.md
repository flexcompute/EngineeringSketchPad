# CAPS Phasing with Optimization Demo

This is a demonstration of the CAPS phasing capability with optimization using Corsair. The demo aims to size a rubber band powered airplane carrying golf balls. 

## Setup
One enviornment variable that will be helpful is the `ESP_PREFIX` environment. This populates the text field when trying to load a pyCAPS script to save time and repetition. For this project, all pyCAPS scripts live in `pyCAPS/`.
```
$> export ESP_PREFIX="pyCAPS/"
```

## File Structure
This directory is organized as follows.
```
.
|--- csm
|--- pyCAPS
|--- README.md
```
The `csm/` and `pyCAPS/` directories contains all the csm and pyCAPS files respectively.

## Design Phases

### Phase 1.1: Size Rectangular Wing
The first phase sizes a rectangular wing, optimizing for lift-to-drag ratio. The wing geometry for this phase is defined and constrained by `csm/wing1.csm` and the optimization formulation is contained in `pyCAPS/sizeWing.py`. 

To start up the CAPS environment, run the following command in the terminal
```
$> serveESP csm/wing1.csm
```
This brings up the CAPS environment in the browser with the CSM file loaded. To start the phasing mode, click on Tool -> Caps. This will prompt a project name and an intent for the phase.

This demo project will be called "ostrich" but can realistically be called anything else. The project name determines the name of the folder created containing all the phasing files. The intent of the first phase is to complete initial sizing of the rectangular wing. 

Once in CAPS mode, we can run the associated pyCAPS script by clicking on Tool -> Pyscript. The text will be populated with the `ESP_PREFIX` environment variable. Enter `pyCAPS/sizeWing.py` to bring up the Python text editor. Click Save and run at the top to run the script. This will solve the corsairlite optimization formulation contained in the Python file and update the wing area and aspect ratio with the optimized geometry.

To end up the phase, click Caps -> Commit Phase. This will exit the CAPS mode and save the current state.

### Phase 1.2: Optimize Taper
The second phase improves upon the first phase by introducing a model for taper. First, the CSM script needs to evolve in order to be able to represent a taper. `csm/wing2.csm` has a design parameter `wing:taper` which represents the taper ratio.

There are two ways to load the new CSM file. The first way is staying in the ESP window and clicking on File -> Open, and entering the path to `wing2.csm`. The second way is closing the window and restarting the session in the terminal as before.
```
$> serveESP csm/wing2.csm
```
To begin the phase once again, click Tool -> Caps. Entering either "ostrich" or "ostrich:1.1" as the project name will begin phase 1.2. Click Tool -> Pyscript and open up `pyCAPS/optTaperAvl.py`. This script runs AVL as a blackbox tool to compute the Oswald efficiency factor $e$ and seeks to maximize $e$. The geometry will be updated with the optimized taper.

Click Caps -> Commit Phase to exit the CAPS mode and save the current state.
