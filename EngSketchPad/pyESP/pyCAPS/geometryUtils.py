
## View or take a screen shot of the geometry configuration. The use of this function to save geometry requires the \b matplotlib module.
# \b <em>Important</em>: If both showImage = True and filename is not None, any manual view changes made by the user
# in the displayed image will be reflected in the saved image.
#
# \param **kwargs See below.
#
# Valid keywords:
# \param viewerType What viewer should be used (default - "capsViewer").
#                   Options: "capsViewer" or "matplotlib" (options are case insensitive).
#                   Important: if $filename is not None, the viewer is changed to matplotlib.
# \param portNumber Port number to start the server listening on (default - 7681).
# \param title Title to add to each figure (default - None).
# \param filename Save image(s) to file specified (default - None). Note filename should
# not contain '.' other than to indicate file type extension (default type = *.png).
# 'file' - OK, 'file2.0Test' - BAD, 'file2_0Test.png' - OK, 'file2.0Test.jpg' - BAD.
# \param directory Directory path were to save file. If the directory doesn't exist
# it will be made. (default - current directory).
# \param viewType  Type of view for the image(s). Options: "isometric" (default),
# "fourview", "top" (or "-zaxis"), "bottom" (or "+zaxis"), "right" (or "+yaxis"),
# "left" (or "-yaxis"), "front" (or "+xaxis"), "back" (or "-xaxis").
# \param combineBodies Combine all bodies into a single image (default - False).
# \param ignoreBndBox Ignore the largest body (default - False).
# \param showImage Show image(s) (default - False).
# \param showAxes Show the xyz axes in the image(s) (default - False).
# \param showTess Show the edges of the tessellation (default - False).
# \param dpi Resolution in dots-per-inch for the figure (default - None).
# \param tessParam Custom tessellation paremeters, see EGADS documentation for makeTessBody function.
# values will be scaled by the norm of the bounding box for the body (default - [0.0250, 0.0010, 15.0]).
def _viewGeometryMatplotlib(self, bodies, **kwargs):

     try:
         import matplotlib.pyplot as plt
         from mpl_toolkits.mplot3d import Axes3D
     except:
         print("Error: Unable to import matplotlib - viewing the geometry with matplotlib is not possible")
         raise ImportError

     def createFigure(figure, view, xArray, yArray, zArray, triArray, linewidth, edgecolor, transparent, facecolor):
         figure.append(plt.figure())

         if view == "fourview":
             for axIndex in range(4):
                 figure[-1].add_subplot(2, 2, axIndex+1, projection='3d')

         else:
             figure[-1].gca(projection='3d')

         for ax in figure[-1].axes:
             ax.plot_trisurf(xArray,
                             yArray,
                             zArray,
                             triangles=triArray,
                             linewidth=linewidth,
                             edgecolor=edgecolor,
                             alpha=transparent,
                             color=facecolor)

     # Plot parameters from keyword Args
     facecolor = "#D3D3D3"
     edgecolor = facecolor
     transparent = 1.0
     linewidth = 0

     viewType = ["isometric", "fourview",
                 "top"   , "-zaxis",
                 "bottom", "+zaxis",
                 "right",  "+yaxis",
                 "left",   "-yaxis",
                 "front",  "+xaxis",
                 "back",   "-xaxis"]

     title = kwargs.pop("title", None)
     combineBodies = kwargs.pop("combineBodies", False)
     ignoreBndBox = kwargs.pop("ignoreBndBox", False)
     showImage =  kwargs.pop("showImage", False)
     showAxes =  kwargs.pop("showAxes", False)

     showTess = kwargs.pop("showTess", False)
     dpi = kwargs.pop("dpi", None)

     tessParam = kwargs.pop("tessParam", [0.0250, 0.0010, 15.0])

     if showTess:
         edgecolor = "black"
         linewidth = 0.2

     filename = kwargs.pop("filename", None)
     directory = kwargs.pop("directory", None)

     if (filename is not None and "." not in filename):
         filename += ".png"

     view = str(kwargs.pop("viewType", "isometric")).lower()
     if view not in viewType:
         view = viewType[0]
         print("Unrecongized viewType, defaulting to " + str(view))

     # How many bodies do we have
     numBody = len(bodies)

     if len(bodies) < 1:
         raise caps.CAPSError(caps.CAPS_BADVALUE, msg = "while saving screen shot of geometry"
                                                      + "\nThe number of bodies in the problem is less than 1!!!" )

     # If ignoring the bounding box body, determine which one it is
     if numBody > 1 and ignoreBndBox == True:
         for bodyIndex, body in enumerate(bodies):
             box = body.getBoundingBox()

             size = sqrt((box[0][0]-box[1][0])**2 +
                         (box[0][1]-box[1][1])**2 +
                         (box[0][2]-box[1][2])**2)
             if size > sizeBndBox:
                 sizeBndBox = size
                 bndBoxIndex = bodyIndex
     else:
         ignoreBndBox = False

     # Loop through bodies 1 by 1
     for bodyIndex in range(numBody):

         # If ignoring the bounding box - skip this body
         if ignoreBndBox and bodyIndex == bndBoxIndex:
             continue

         body = bodies[bodyIndex]
         box = body.getBoundingBox()

         size = sqrt((box[0][0]-box[1][0])**2 +
                     (box[0][1]-box[1][1])**2 +
                     (box[0][2]-box[1][2])**2)


         params[0] = tessParam[0]*size
         params[1] = tessParam[1]*size
         params[2] = tessParam[2]

         tess = body.makeTessBody(params)
         faces = body.getBodyTopos(egads.FACE)

         # If not combining bodies erase arrays
         if not combineBodies:
             del xArray[:]
             del yArray[:]
             del zArray[:]
             del triArray[:]

         del gIDArray[:]

         # Loop through face by face and get tessellation information
         for faceIndex in range(numFace):
             points, uv, ptype, pindex, tris, triNeighbor = tess.getTessFace(faceIndex+1)

             plen = len(points)
             tlen = len(tris)

             for i in range(plen):
                  globalID = tess.localToGlobal(faceIndex+1, i+1)

                  tri = []
                  if globalID not in gIDArray:
                      xArray.append(points[3*i  ])
                      yArray.append(points[3*i+1])
                      zArray.append(points[3*i+2])
                      gIDArray.append(globalID)

             for i in range(tlen):
                 tri = []
                 globalID = tess.localToGlobal(faceIndex+1, tris[3*i + 0])
                 tri.append(globalID-1 + numPointTot)

                 globalID = tess.localToGlobal(faceIndex+1, tris[3*i + 1])
                 tri.append(globalID-1 + numPointTot)

                 globalID = tess.localToGlobal(faceIndex+1, tris[3*i + 2])
                 tri.append(globalID-1 + numPointTot)

                 triArray.append(tri)

         if combineBodies:
             numPointTot += len(gIDArray)
         else:
             createFigure(figure, view,
                          xArray, yArray, zArray, triArray,
                          linewidth, edgecolor, transparent, facecolor)

     if combineBodies:
         createFigure(figure, view,
                      xArray, yArray, zArray, triArray,
                      linewidth, edgecolor, transparent, facecolor)

     for fig in figure:

         if len(fig.axes) > 1:
             fig.subplots_adjust(left=2.5, bottom=0.0, right=97.5, top=0.95,
                                 wspace=0.0, hspace=0.0)
         if title is not None:
             fig.suptitle(title)

         if dpi is not None:
             fig.set_dpi(dpi)

         for axIndex, ax in enumerate(fig.axes):
             xLimits = ax.get_xlim()
             yLimits = ax.get_ylim()
             zLimits = ax.get_zlim()

             center = []
             center.append((xLimits[1] + xLimits[0])/2.0)
             center.append((yLimits[1] + yLimits[0])/2.0)
             center.append((zLimits[1] + zLimits[0])/2.0)

             scale = max(abs(xLimits[1] - xLimits[0]),
                         abs(yLimits[1] - yLimits[0]),
                         abs(zLimits[1] - zLimits[0]))

             scale = 0.6*(scale/2.0)

             ax.set_xlim(-scale + center[0], scale + center[0])
             ax.set_ylim(-scale + center[1], scale + center[1])
             ax.set_zlim(-scale + center[2], scale + center[2])

             if showAxes:
                 ax.set_xlabel("x", size="medium")
                 ax.set_ylabel("y", size="medium")
                 ax.set_zlabel("z", size="medium")
             else:
                 ax.set_axis_off() # Turn the axis off

             #ax.invert_xaxis() # Invert x- axis
             if axIndex == 0:

                 if view in ["top", "-zaxis"]:
                     ax.view_init(azim=-90, elev=90) # Top

                 elif view in ["bottom", "+zaxis"]:
                     ax.view_init(azim=-90, elev=-90) # Bottom

                 elif view in ["right", "+yaxis"]:
                     ax.view_init(azim=-90, elev=0) # Right

                 elif view in ["left", "-yaxis"]:
                     ax.view_init(azim=90, elev=0) # Left

                 elif view in ["front", "+xaxis"]:
                     ax.view_init(azim=180, elev=0) # Front

                 elif view in ["back", "-xaxis"]:
                     ax.view_init(azim=0, elev=0) # Back

                 else:
                     ax.view_init(azim=120, elev=30) # Isometric

             if axIndex == 1:
                 ax.view_init(azim=180, elev=0) # Front

             elif axIndex == 2:
                 ax.view_init(azim=-90, elev=90) # Top

             elif axIndex == 3:
                 ax.view_init(azim=-90, elev=0) # Right

     if showImage:
         plt.show()

     if filename is not None:

         for figIndex, fig in enumerate(figure):
             if len(figure) > 1:
                 fname = filename.replace(".", "_" + str(figIndex) + ".")
             else:
                 fname = filename

             if directory is not None:
                 if not os.path.isdir(directory):
                     print("Directory ( " + directory + " ) does not currently exist - it will be made automatically")
                     os.makedirs(directory)

                 fname = os.path.join(directory, fname)

             print("Saving figure - ", fname)
             fig.savefig(fname, dpi=dpi)