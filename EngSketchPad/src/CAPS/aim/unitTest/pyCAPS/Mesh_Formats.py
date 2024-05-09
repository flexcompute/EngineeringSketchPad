
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

# Check if the shared library exists
Mesh_Formats = []
for writer in Writers:
    if os.path.isfile( os.path.join(os.environ["ESP_ROOT"],"lib",writer+ext) ):
        Mesh_Formats += [writer.replace('Writer','')]