import config as cfg
import os
import subprocess
import random
import shutil

pwd = os.getcwd()

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
        p1 = subprocess.call( a, shell=False, timeout=None )


        print("==========================================")
        print("Input File:  ", outputDir + "/randomMesh.stl")
        print("Saving File: ", outputDir + "/input.vtk")
        # print("Running ", cfg.num_simulations, " simulation(s)...")
        print("==========================================")
        os.chdir(outputDir)
        # print(os.getcwd())
        shutil.copy(cfg.randomMesh_geo_path, "./randomMesh.geo")

        b = ["gmsh"]
        b += ["randomMesh.geo"]
        b += ["-o"]
        b += ["input.vtk"]
        b += ["-3"]
        # gmsh randomMesh.geo -o input.vtk -3
        p2 = subprocess.call(b, shell=False, timeout=None)

        os.chdir(pwd)

        print("=================")
        print("Completed ", i+1, "/", cfg.num)
        print("=================")

    print("==========================================")
    print("Generation Complete")
    print("==========================================")
