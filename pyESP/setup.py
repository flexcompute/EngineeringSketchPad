
from setuptools import setup, find_packages

setup(
    name='pyESP',
    version='1.28.0',
    description='Engineering Sketch Pad',
    author='Mike Park',
    author_email='mike@flexcompte.com',
    packages=find_packages(),  # Automatically finds 'pyEGADS' and 'pyCAPS'
    install_requires=[
            'ctypes',
            'weakref',
            'atexit',
        ],
)