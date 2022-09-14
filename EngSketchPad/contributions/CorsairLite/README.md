# Corsairlite Optimization Examples
This repository contains a number of examples that demonstrate the uses of Corsairlite in a variety of optimization problems ranging from simple test problems to full aircraft optimization problems. 

## Corsairlite Documentation
The Corsairlite documentation is provided in `corsairlite.md`. 

## Test Problems
The following files demonstrate the use of Corsairlite on simple test optimization problems.
```
qp.py - Quadratic Program
gp.py - Geometric Program
sp.py - Signomial Program
```

`boydBox.py` defines the geometric program test problem given by Boyd: optimization of a box for maximum volume given constraints on the areas of the box floor and walls as well as the aspect ratio of various sides.

## Aircraft Optimization Problems
The following files demonstrate the use of Corsairlite on aircraft optimization problems.
```
hoburg_gp.py
hoburg_blackbox.py
```
These problems optimize an aircraft for minimum fuel weight on three segments as defined by Hoburg and Abbeel in "Geometric Programming for Aircraft Design Optimization". 

The first file, `hoburg_gp.py` defines all of the constraints as a geometric program as described in the original paper. The results are saved in a pickle file called `hoburg.pl`.

The second file, `hoburg_blackbox.py` has all of the same definitions, except the drag model fitted to XFOIL data is instead moved into a Blackbox constraint. The constraint is set using a `RuntimeConstraint` class and is evaluated at runtime for each iteration. The optimization is solved using Logspace Sequential Quadratic Programming and the initial guess is taken from the results of the geometric programming optimzation from `hoburg.pl`. Note that `hoburg_blackbox.py` requires `scipy` which is not distributed with ESP but can be installed with `pip install scipy`.

## ESP Phasing with Corsairlite
There is one example of phasing with Corsairlite. `capsPhase/` directory contains the files and instructions to design and optimize a rubber band powered aircraft carrying golf balls using the ESP phasing architecture and Corsairlite. Read the `capsPhase/README.md` for more details.
