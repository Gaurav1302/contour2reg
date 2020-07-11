# Contour to Registartion (2D-3D Registration)

## Prerequisites

1. Setup your Python environment using the `requirement.txt` file. The code was tested on Python `3.7.5`.
2. Install the following softwares:
  Blender (`v2.83`)

## Generate Random Meshes

To generate random meshes, set up the configurations in `Generate/config.py` file and run the following script.

```
cd Generate/
python blender_generate.py
```

A Random organ-like mesh (`randomMesh.stl`) is created in `config.Data_path/` folder with a different folder within this folder for each random mesh. 
