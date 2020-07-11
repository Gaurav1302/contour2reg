import config as cfg
import os
import subprocess
import random


if cfg.MeshExists:
    # Check if the mesh is an'.stl' file
    if not cfg.MeshPath.endswith( ".stl" ):
        raise IOError("Given mesh must be a .stl file.")
    else:
        print("==========================================")
        print("Mesh Found! Skipping mesh generation.")
        print("==========================================")
else:
    # Create Random Meshes
    print("==========================================")
    print("Generating ", cfg.num, " Random Meshes")
    print("==========================================")
    for i in range(cfg.num):
        # Create folder for this sample
        outputDir = os.path.join(cfg.Data_path,str(i+1))
        if (not os.path.isdir(outputDir)):
            os.mkdir(outputDir)
        else:
            print("Directory already exists..")

        # Call blender script (generateRandomMesh.py) to generate random meshes
        a = ["blender"]
        a += ["--background"]
        a += ["-noaudio"]
        a += ["--python"]
        a += ["generateRandomMesh.py"]
        a += ["--"]
        a += ["--outdir"]
        a += [outputDir]
        # blender --background -noaudio --python generateRandomMesh.py -- --outdir "../Data/1/"
        p = subprocess.call( a, shell=False, timeout=None )

        print("=================")
        print("Completed ", i+1, "/", cfg.num)
        print("=================")

    print("==========================================")
    print("Generation Complete")
    print("==========================================")
