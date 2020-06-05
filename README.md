# Mesh generator for arbitrarily shaped fractures
Python script to automatically generate the finite-volume mesh for any given complex fracture network geometry.

An open-source platform SALOME (https://www.salome-platform.org/) is used for creating the geometry and mesh which is conforming to the complex fracture network. 
SALOME supports CAD (computer-aided design) modeling which is very useful for complex shapes. 

This script writes the generated mesh files in OpenFOAM format. Salome to OpenFOAM mesh conversion function is taken from https://github.com/nicolasedh/salomeToOpenFOAM.
