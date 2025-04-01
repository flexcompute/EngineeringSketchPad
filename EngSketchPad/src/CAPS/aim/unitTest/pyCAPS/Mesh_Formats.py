
import os
import sys

# Get library extension
if sys.platform.startswith('darwin') or sys.platform.startswith('linux'):
    ext = '.so'
elif sys.platform.startswith('win32'):
    ext = '.dll'

# Get the path to the meshWriter directory
meshWriter = os.path.join(os.path.dirname(__file__), '..', '..', 'meshWriter')

# List all directories in meshWriter
Writers = [ os.path.basename(f.path) for f in os.scandir(meshWriter) if f.is_dir() ]

# bdf is special
Writers.append('bdfsmallWriter')
Writers.append('bdflargeWriter')
Writers.append('bdffreeWriter')

# Check if the shared library exists
Mesh_Formats = []
for writer in Writers:
    if os.path.isfile( os.path.join(os.environ["ESP_ROOT"],"lib",writer+ext) ):
        Mesh_Formats += [writer.replace('Writer','')]

# Remove file formats not supporting quads
Mesh_Formats_Quad = Mesh_Formats.copy()
if "fast"   in Mesh_Formats_Quad: Mesh_Formats_Quad.remove("fast")
if "exodus" in Mesh_Formats_Quad: Mesh_Formats_Quad.remove("exodus")
