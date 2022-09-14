# CorsairLite Quickstart Guide
This file contains a number of tips for successfully using the CorsairLite software package.  

## Overview
CorsairLite is a light weight version of Corsair included in the ESP distribution.  CorsairLite allows users to define an optimization problem in a simple readable syntax, and then solve it using one of a number of provided optimization solvers.  This repository includes a number of example files that users can explore, though this file should contain all you need to get started.  

## The Header
The standard header in CorsairLite will be as follows:
```
from corsairlite import units
from corsairlite.optimization import Variable, Constant, Formulation
from corsairlite.optimization.solve import solve
```
These three statements import the most critical classes for developing optimization problems in CorsairLite, and should cover the needs of most new users.

Some additional import options can be seen in the `hoburg.gp` file that deal with storing data via python's `dill` package.

## Variables
A variable is a parameter that will be modified by the optimizer to find the optimal solution.  A variable `x` is defined by the statement:
```
x = Variable(name = 'x', guess = 1.0, units = 'm', description = 'The Variable X')
```
The `Variable()` class constructor requires inputs name, guess, and units, though a description is highly recommended.  Best practice is to included the keyword flag, but notation can be condensed into:
```
x = Variable('x', 1.0, 'm', 'The Variable X')
```
if so desired.  Name must be a string.  The guess either a float or an instance of a CorsairLite quantity (that is a float/units pair), but a quantity input must have units that are consistent with the defined units in the variable constructor (that is a quantity with units of 'feet' can be used in a variable constructor with 'm' for meters, but using a quantity with units 'seconds' would throw an error).  The units value must be a string that matches with the CorsairLite units definition.  These are the same as the `pint` python package.  The description value must be a string.  If a variable has no units, the strings '' and '-' are both parsed as a dimensionless quantity.

Variables in an optimization problem must have unique names.  Best practice for sequential variables would be to use `x_1`, `x_2` etc.  Variables must be scalars in CorsairLite, though support for vector and matrix variables will come in a future version.

## Constants
A constant is a parameter that a user wishes to use symbolically, but is not modified by the optimizer in determining the optimal solution.  In post-processing, some types of optimization solvers support the ability to determine the sensitivity of the optimal solution with respect to these constants (LP, QP, GP).  

A constant is defined the same way as a variable, but using the `Constant()` constructor:
```
g = Constant(name = 'g', guess = 9.81, units = 'm/s^2', description = 'Gravitational Constant')
```

Though variables and constants can be defined in any order, best practice is to split the Variables and Constants into two separate blocks.  Within each block, best practice is to order variables and constants either in order of appearance in the optimization problem (recommended for small problems) or in alphabetical order (recommended for large problems).

## Expressions
Variables and Constants combine to form `Expressions` using standard arithmetic operators: `+`, `-`, `*`, `/`, and `**`.  Note that the current version does not support raising to the power of variables (ex: `2**x`), though this will be supported in future versions.

An expression can be constructed as follows:
```
x   = Variable(name = 'x', guess = 1, units = '-', description = 'The X variable')
y   = Variable(name = 'y', guess = 1, units = '-', description = 'The Y variable')
z   = Variable(name = 'z', guess = 1, units = '-', description = 'The Z variable')

c   = Constant(name = 'c', guess = 3, units = '-', description = 'A Constant')

expr = x**2 * y**3 * z + c * x / ( y * z )
```

## The Objective
The objective function can be any valid CorsairLite expression.  The current version does not support arbitrary python functions as objective functions, though some version of this feature will be included in future versions.  The objective is defined when constructing a formulation (see below).

## The Constraints
CorsairLite generally expects constraints to be defined in the form of a python list.  The simplest form of constraint is simply defined by using an expression in combination with a relational operator, ie:
```
x**2 + y**2 <= 1
```

Constraints can also be defined using two expressions:
```
x + y == c*u + w
```

Best practice is generally to define an empty list, and then append constraints using the `+=` operator:
```
constraints = []

constraints += [x <= 1]
constraints += [x**2 == y, z >= y]
```

More complex constraints can be constructed when an explicit relationship is not available.  These are called `RuntimeConstraints` and act as black box functions in the optimization problem.  An example of this is shown in `hoburg_blackbox.py`, but users interested in using this capability are encouraged to contact the developers.

## Defining a Formulation
A formulation (ie, the optimization problem to be solved) is constructed using the `Formulation()` constructor.  The `Formulation()` constructor takes two inputs:
```
formulation = Formulation(objective, constraints)
```

## Selecting an Optimization Algorithm
CorsairLite contains logic that will select the best optimization algorithm for the formulation that has been defined.  However, this logic does take some time to process, and so best practice is for the user to select the desired optimization algorithm.  The algorithm is set via the attribute in the formulation option:
```
formulation.solverOptions.solveType
```

Options available in CorsairLite are:
```
formulation.solverOptions.solveType = 'lp'
formulation.solverOptions.solveType = 'qp'
formulation.solverOptions.solveType = 'gp'
formulation.solverOptions.solveType = 'sp'
formulation.solverOptions.solveType = 'pccp'
formulation.solverOptions.solveType = 'sqp'
formulation.solverOptions.solveType = 'lsqp'
formulation.solverOptions.solveType = 'slcp'
```

For most users, the 'sqp' flag is probably most appropriate.  However, users interested in other algorithms are encouraged to contact the developers.  Using the correct algorithm can significantly alter how long the optimizer takes to solve.  

The type of solver can also be changed using:
```
formulation.solverOptions.solver = 'cvxopt'
```
In CorsairLite, CVXOPT is the only available solver.  However, the full version of Corsair supports more options.

## Setting Solver Options
Many other options are available in the `formulation.solverOptions` object.  New users are unlikely to need to adjust any of these parameters, though all of the possible options can be inspected at the bottom of the `corsairlite.optimization.formulation` file.

## Solving the Optimization Problem
The optimization problem is solved simply by calling:
```
res = solve(formulation)
```

## Printing and Extracting the Results
The `res` (short for 'result') object is an instance of a custom python class that contains a large amount of data associated with the optimization solve.  Calling `print(res.result(5))` will print a summary of the major details, where the input to the function is the number of decimals printed.  

The optimal variables are extracted using `res.variables['x']` where the input to the dictionary object is the `name` associated with the variable.

## Optional Dependencies
CorsairLite has been stripped down to limit the size of the ESP distribution.  Functionality can be enhanced by installing the `torch` and `pandas` python packages

## Upgrading to the Full Corsair Package
Users who wish to utilize more advanced features are encouraged to contact the developers to upgrade to the full Corsair package.  This package is currently distributed via invitation to a private GitHub repository, but will be made public in the near future.

## Contact Information
For any additional questions, feel free to contact the developers:
Cody Karcher, Primary Corsair Developer, ckarcher@mit.edu
Dongjoon Lee, CorsairLite Support, dongjoon@mit.edu

