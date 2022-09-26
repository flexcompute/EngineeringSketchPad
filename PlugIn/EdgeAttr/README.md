# Installation
Skip this step if using [pre-pulit ESP with edgeAttr included](https://github.com/flexcompute/EngineeringSketchPad/releases/tag/ESP121-Pre-Built-edgeAttr).

EdgeAttr folder can (and should) be located outside the ESP root folder.
It will build a plug-in that is not part of the normal release and so it should be in a place that is independent.

With a terminal/shell that has the ESP environment, change directories into EdgeAttr.
Then simply type `make` (or `nmake -f NMakefile` in Windows).
This builds the plug-in and puts the resultant shared object/DLL in `$ESP_ROOT/lib`.

# Usage and Example
Set an attribute for the edge going through set of points (x1,y1,z1), (x2,y2,z2), (x3,y3,z3) ...
```
udparg  edgeAttr  attrname  $edgeName attrstr $leadingEdge
udprim  edgeAttr  xyzs "0;0;0; 0;1;0; 0;2;0"
```
