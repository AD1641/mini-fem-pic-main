# Generating input files

Input files for this prototype application can be done using the `SALOME`
program.

1. Generate geometry using salome's `geometry` tools
2. Select `Create groups` after right-clicking the geometry objects in the left
   sidebar
3. Create groups of faces: 
    -  "inlet" faces,
    -  wall or "fixed" faces
4.  Select the `meshes` tools
5. Generate and compute the meshes for the geometry
6. Right click on the geometry in the sidebar and select "export as -> DAT"
7. Repeat 6. for the "inlet" and "fixed" faces.



