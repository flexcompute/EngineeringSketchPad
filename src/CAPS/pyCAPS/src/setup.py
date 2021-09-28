# To execute type:
# python setup.py sdist bdist_wheel
# bdist_wheel - build a wheel
# sdist - build a tar ball of the source

try:
    print "Using setuptools"
    from setuptools import setup
#     from setuptools import Command
    from setuptools import Extension
except ImportError: # Try getting them from the distutils
    print "Using distutils"
    
    from distutils.core import setup
#     from distutils.core import Command
    from distutils.extension import Extension

import os
import platform
import sys
from re import sub as reSub

# Set version number
version = "2.2.1"

# macOS and Linux -> -lcaps -locsm -legads -ludunits2 -lwsserver $(PYTHONLIB) -lm -lz 
# Windows -> caps.lib wsserver.lib ws2_32.lib ocsm.lib egads.lib udunits2.lib z.lib $(PYTHONLIB) 

# libs = ["caps", "ocsm", "egads", "udunits2", "wsserver"]
libs = ["capsstatic", "ocsm", "egadstatic", "udunits2", "wsserver"]

if platform.system() == "Windows":
    libs += ["ws2_32"] #"legacy_stdio_definitions", "Advapi32" - not sure about these
    libs += ["udunits2"]
	   
# Get enviroment variable
ESP_ROOT = os.environ["ESP_ROOT"]

compile_args = []
link_args = []

# Grab extra "flags" from the command line
temp = sys.argv[:]
for i in temp:

    if "--CC=" in i:
        os.environ["CC"] = reSub(' +',' ',i).replace("--CC=", "")
        sys.argv.remove(i)
        
    if "--COPTS=" in i:
        compile_args = i.replace("--COPTS=", "").split()
        sys.argv.remove(i)
        
    if "--SOFLGS=" in i:
        link_args = i.replace("--SOFLGS=", "").split()
        sys.argv.remove(i)

capsModule = Extension("pyCAPS", 
                       sources = ["pyCAPS.c"],
                       libraries = libs,
                       library_dirs = [ESP_ROOT+ "/lib"],
                       include_dirs = [ESP_ROOT+'/include'],
                       extra_compile_args = compile_args, 
                       extra_link_args = link_args)

# class CleanCommand(Command):
#     """Custom clean command to tidy up the project root."""
#     user_options = []
#     def initialize_options(self):
#         pass
#     def finalize_options(self):
#         pass
#     def run(self):
#         os.system('rm -vrf ./build ./*.egg-info')

setup(name = 'pyCAPS',
      version = version,
      description = 'Python Extension Module to CAPS',
      long_description = '''
pyCAPS is a Python extension module to interact with Computational Aircraft Prototype Syntheses (CAPS) routines in the 
Python environment. Written in Cython, pyCAPS natively handles all type conversions/casting, while logically grouping CAPS function 
calls together to simplify a user's experience. Additional functionality not directly available through the CAPS API 
is also provided.

Place clearance statement here
''',
      author='Ryan Durscher',
      url = 'https://acdl.mit.edu/ESP',
      platforms = ["macOS", "Linux", "Windows"],
      license = "",
      classifiers=[
                   # How mature is this project? Common values are
                   #   3 - Alpha
                   #   4 - Beta
                   #   5 - Production/Stable
                   'Development Status :: 4 - Beta',
                   'Programming Language :: Python :: 2.6',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3.3',
                   'Programming Language :: Python :: 3.4',
                   'Programming Language :: Python :: 3.5', 
                   'Programming Language :: Python :: 3.6', 
                   'Programming Language :: Python :: 3.7'],
      
      #  These are all part of the standard Python dist
      #install_requires=["cypython.version", "libc.stdlib", "libs.stdio", "json", "os", "time", "math", "re", "webbrowser", "struct"],
      
      ext_modules = [capsModule],
#       include_package_data=True,
#       package_data={'': package}
#       cmdclass = { 'clean': CleanCommand},
)
