from pyCAPS import capsViewer

# Initialize viewer
viewer = capsViewer()

# Define nodes
nodes =  [ [0, 0, 0], 
           [1, 1, 0], 
           [0, 1, 0], 
           [2, 0, 0],
           [2, 2, 0], 
           [3, 2, 0], 
           [5, 0, 0], 
           [7, 3, 0]]

# Define element connectivity
connectivity = [ [1, 2, 3], 
                 [1, 4, 2],
                 [4, 5, 2], 
                 [6, 5, 4], 
                 [4, 7, 8, 6]]

# Define scalar data at nodes
color = [1, 1, 2, 2, 4, 2, 3, 4]

# Add nodes
viewer.addVertex( nodes)

# Add element connectivity 
viewer.addIndex( connectivity)

# Add color/scalar data values at nodes 
viewer.addColor( color)

# Add second set of values 
color = [0, 30, 40, 50, 60, 90, 10, 100]
viewer.addColor( color, name = "Color 2")

# Add "mesh" connectivity - enables ploting the elemental mesh
viewer.addLineIndex(connectivity)

# Set mesh color - a single value set a monotone color
viewer.addLineColor(1, colorMap = "grey")

# Take all current "items" that have been added and create a triangle graphic primitive
viewer.createTriangle()

# Start the webserver 
viewer.startServer()

# Optionally during each "add" function call we could have kept track of the "items" 
# being created and explicitly used them to create the triangle graphic primitive.
#
# For example: 
# 
# Setup empty list to hold wvItems
#itemList = []

#itemList.append(viewer.addVertex( nodes )

#itemList.append(viewer.addIndex( connectivity ))

#itemList.append(viewer.addColor( color ))

#itemList.append(viewer.addLineIndex( connectivity ))

#itemList.append(viewer.addLineColor( 1, "grey" ))

#viewer.createTriangle(itemList)
