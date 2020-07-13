import config as cfg
import os
import subprocess
import random
import shutil

# random.seed(30)

pwd = os.getcwd()

if cfg.MeshExists_ps:
    # Check if the mesh is an '.stl' file
    if not cfg.MeshPath_ps.endswith( ".stl" ):
        raise IOError("Given mesh must be a .stl file. Aborting...")
        exit()
    else:
        print("==========================================")
        print("Mesh Found! Skipping mesh generation.")
        print("==========================================")
        MeshPath = cfg.MeshPath_ps
else:
    # Create Random Mesh
    print("==========================================")
    print("Generating Random Mesh.")
    print("==========================================")
    # Call blender script (generateRandomMesh.py) to generate random meshes
    a = "blender --background -noaudio --python generateRandomMesh.py -- --outdir " + cfg.Data_path_ps
    os.system(a)
    # Save the path of this random mesh
    MeshPath = os.path.join(cfg.Data_path_ps, "randomMesh.stl")


for i in range(cfg.startNum_simulations_ps, cfg.num_simulations_ps):
    # Create folder for this simulation
    outputDir = os.path.join(cfg.Data_path_ps,str(i+1))
    if (not os.path.isdir(outputDir)):
        os.mkdir(outputDir)

    # Copy the random mesh to this folder for generating its simulation
    shutil.copy(MeshPath, outputDir + "/randomMesh.stl")

    print("==========================================")
    print("Input File:  ", outputDir + "/randomMesh.stl")
    print("Saving File: ", outputDir + "/input.vtk")
    # print("Running ", cfg.num_simulations, " simulation(s)...")
    print("==========================================")
    os.chdir(outputDir)
    shutil.copy(cfg.randomMesh_geo_path, "./randomMesh.geo")
    os.system("gmsh randomMesh.geo -o input.vtk -3")


    # Set up random boundary conditions for this sample
    print("==========================================")
    print("Set up boundary conditions (simulationGen)")
    print("==========================================")
    os.system("simulationGen . i+1")

    # Generate simulation files:
    print("==========================================")
    print("Generate simulation files (vtk2elmer)")
    print("==========================================")
    os.system("vtk2elmer .")

    # Run the simulation:
    print("==========================================")
    print("Running Simulations (ElmerSolver)")
    print("==========================================")
    os.chdir("./simulation/")
    if not os.path.isdir("./result/"):
        os.mkdir("./result/")
    os.system("ElmerSolver case.sif")

    # Voxelize the input
    print("==========================================")
    print("Running Voxelize on output (voxelize)")
    print("==========================================")
    os.chdir(pwd) # Go back to "Generate/" folder
    b = "voxelize " + os.path.join(outputDir, "simulation/result/case_t0001.vtu") + " " + str(cfg.cubeSideLength) + " 0.025 " + outputDir + " -s 0.3 -b " + os.path.join(outputDir, "simulation.vtu") + " -m 0.1"
    os.system(b)


    print("=================")
    print("Completed ", i+1, "/", cfg.num_simulations_ps)
    print("=================")

print("==========================================")
print("Generation Complete")
print("==========================================")
