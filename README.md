# Engineering Sketch Pad (ESP)
User Manual:  
https://flexcompute.github.io/EngineeringSketchPad/EngSketchPad/ESP/ESP-help.html  
Valid CSM statements:  
https://flexcompute.github.io/EngineeringSketchPad/EngSketchPad/ESP/ESP-help.html#sec5.4

# Staging prebuild artifacts for commit

```
[edit download.sh and stage.sh to confirm versions]
./download.sh
./stage.sh
```
which rebuids libegads.so to resolve indirect OpenCASCADE libraries.
Commit new, removed, and modfied files.
