Unpack the g'zipped tar image attached. This can (and should
be done outside the ESP distribution). It will build a plug-in that
is not part of the normal release and so it should be in a place
that is independent.

With a terminal/shell that has the ESP environment, change
directories into EdgeAttr (the directory made by unpacking the tar
image). Then simply type "make" (or "nmake -f NMakefile" in Windows).
This builds the plug-in and puts the resultant shared object/DLL in
$ESP_ROOT/lib.