import argparse
import os
import sys
import subprocess
import shutil
import surface2Volume
import setupRandomSimulation
import vtk2Elmer
import voxelize
import voxelizeDisplacementField
import makePartialSurface
import traceback
import gc
from utils import *
#import tracemalloc
#tracemalloc.start()

##################################################
## Argument passing
parser = argparse.ArgumentParser(description="Generate a number of random simulations")
parser.add_argument("--num", type=int, help="Number of simulations to create.")
parser.add_argument("--outdir", type=path, default="../Data/", help="Folder in which to save the simulation files.")
parser.add_argument("--mesh", type=filepath, help="Surface mesh (stl) on which to run the simulation. If empty, will generate random meshes.")
parser.add_argument("--startNum", type=int, default=0, help="Optionally start from this simulation number")
parser.add_argument("--size", type=float, default=0.3, help="Size of simulation/voxelization domain")
parser.add_argument("--grid_size", type=float, default=64, help="Number of voxels per dimension (results in grid_size^3 voxels)")
parser.add_argument("--rigid_displacement", action="store_true", help="Optionally add random rigid displacement to intraoperative data. Note: this option was off when creating the results shown in the publication." )
parser.add_argument("--retry_with_less_force", action="store_true", help="If simulation fails, decrease force magnitude and retry.")

skip_group = parser.add_argument_group("pipeline skipping", description="The following options allow to skip steps of the pipeline if they have already been performed earlier (mainly intended for debugging a single step of the pipeline)." )
skip_group.add_argument("--skip_mesh_generation", action="store_true", help="Optionally do not generate meshes (if they were already generated beforehand)")
skip_group.add_argument("--skip_meshing", action="store_true", help="Optionally skip creating volumes from surfaces (if they were already generated beforehand)")
skip_group.add_argument("--skip_simulation_setup", action="store_true", help="Optionally skip setting up simulation (if they were already generated beforehand)")
skip_group.add_argument("--skip_simulations", action="store_true", help="Optionally skip simulation (if they were already generated beforehand)")
skip_group.add_argument("--skip_voxelization", action="store_true", help="Optionally skip voxelization" )

args = parser.parse_args()

##################################################
## Generate random mesh if needed:
print("==========================================")
print("Mesh setup")
if not args.skip_mesh_generation:
    if not args.mesh:

        print("Generating {:d} random meshes".format(args.num))
        #blender --background -noaudio --python generateRandomMesh.py -- [other arguments]
        a = ["blender"]
        a += ["--background"]
        a += ["-noaudio"]
        a += ["--python"]
        a += ["generateRandomMesh.py"]
        a += ["--"]
        a += ["--num"]
        a += [str(args.num)]
        a += ["--startNum"]
        a += [str(args.startNum)]
        a += ["--outdir"]
        a += [args.outdir]
        p = subprocess.call( a, shell=False, timeout=None )
    else:

        if not args.mesh.endswith( ".stl" ):
            raise IOError("Given mesh must be a .stl file.")

        # If a mesh is given, simply copy it to every folder:
        print("Duplicating given mesh {:d} times".format(args.num))
        #for i in range(0,args.num):
        for i in range(args.startNum,args.startNum+args.num):
            curDir = directoryForIndex( i, args.outdir )
            shutil.copyfile( args.mesh, os.path.join( curDir, "surface.stl" ) )


for i in range(args.startNum,args.startNum+args.num):
    print("Mesh number:", i)
    curDir = directoryForIndex( i, args.outdir )
    curSeed = i+1

    #snapshot1 = tracemalloc.take_snapshot()
    ##################################################
    ## Generate volume mesh if needed:
    print("==========================================")
    print("Convert surfaces to volume meshes", i)
    if not args.skip_meshing:
        if args.mesh and i != args.startNum:   # If a volume for this mesh already exists, simply duplicate

            firstDir = directoryForIndex( args.startNum, args.outdir )
            firstVolumeFile = os.path.join( firstDir, "volume.vtk" )
            newVolumeFile = os.path.join( curDir, "volume.vtk" )
            print("Duplicating volume mesh to:", newVolumeFile)
            shutil.copyfile( firstVolumeFile, newVolumeFile )

        else:   # Otherwise try to generate a volume from the surface:
            try:
                surface2Volume.surface2Volume( os.path.join( curDir, "surface.stl" ), "volume.vtk" )
            except KeyboardInterrupt:
                exit()
            except:
                print("\tFailed to create volume for mesh {}. Skipping.".format(i))
                print(sys.exc_info())
                continue


    gc.collect()

    ##################################################
    ## Generate random boundary conditions (and material properties)
    print("==========================================")
    print("Generate simulations", i)
    if not args.skip_simulation_setup:
        volumeFile = os.path.join( curDir, "volume.vtk" )
        outputDir = os.path.join( curDir, "sim" )
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        try:
            setupRandomSimulation.mesh2simulation( volumeFile, curSeed, outputDir )
        except KeyboardInterrupt:
            exit()
        except:
            print("\tFailed to create simulation files for mesh {}. Skipping.".format(i))
            print(sys.exc_info())
            continue

    ##################################################
    ## Convert the files to Elmer format and run the Elmer simulation:
    print("==========================================")
    print("Run simulations", i)
    if not args.skip_simulations:
        try:
            simDir = os.path.join( curDir, "sim" )
            simulationFile = os.path.join( simDir, "simulation.vtk" )

            simulationSuccess = False
            maxForces = 3
            forcePrefactor = 1
            gaveUp = False
            while not simulationSuccess and not gaveUp:

                vtk2Elmer.convert( simulationFile, curSeed, maxForces=maxForces, forcePrefactor=forcePrefactor )

                outputDir = os.path.join( simDir, "result" )
                if not os.path.exists( outputDir ):
                    os.makedirs( outputDir )
                # Try to delete any previous results:
                resultFileName = os.path.join( outputDir, "case_t0001.vtu" )
                if os.path.exists( resultFileName ):
                    os.remove( resultFileName )

                # Run the ElmerSolver to solve the given deformation problem:
                p = subprocess.call(["ElmerSolver", "case.sif"], cwd=simDir)

                # If the result file has been created, consider this a success:
                if os.path.exists( resultFileName ):
                    simulationSuccess = True
                else:
                    if forcePrefactor > 0.01 and args.retry_with_less_force:
                        forcePrefactor = forcePrefactor*0.3
                        print("Could not solve simulation. Trying a with a smaller force:",forcePrefactor)
                    else:
                        gaveUp = True

        except KeyboardInterrupt:
            exit()
        except:
            print("\tFailed to run simulation {}. Skipping.".format(i))
            print(sys.exc_info())
            continue


    ##################################################
    ## Voxelize the simulation input and output:
    print("==========================================")
    print("Voxelization", i)
    if not args.skip_voxelization:
        volumeFile = os.path.join( curDir, "volume.vtk" )
        simDir = os.path.join( curDir, "sim" )
        simulationResultFile = os.path.join( simDir, "result", "case_t0001.vtu" )

        try:

            if not os.path.isfile( volumeFile ):
                raise IOError("Could not find {}".format(volumeFile))

            if not os.path.isfile( simulationResultFile ):
                raise IOError("Could not find {}".format(simulationResultFile))

            # Create the pre-oprative signed distance field:
            reader = vtkUnstructuredGridReader()
            reader.SetFileName( volumeFile )
            reader.Update()
            mesh = vtkUnstructuredGrid()
            mesh.DeepCopy(reader.GetOutput())
            del reader
            if mesh.GetNumberOfPoints() <= 0:    # Loaded successfully?
                raise IOError("Could not read {}".format(volumeFile))
            mesh = unstructuredGridToPolyData( mesh )

            grid = voxelize.createGrid( args.size, args.grid_size )

            surface = extractSurface( mesh )
            voxelize.distanceField( surface, grid, "preoperativeSurface", signed=True )

            # Create the intra-operative distance field:
            reader = vtkXMLUnstructuredGridReader()
            reader.SetFileName( simulationResultFile )
            reader.Update()
            mesh = vtkUnstructuredGrid()
            mesh.DeepCopy(reader.GetOutput())
            if mesh.GetNumberOfPoints() <= 0:    # Loaded successfully?
                raise IOError("Could not read {}".format(simulationResultFile))
            mesh = unstructuredGridToPolyData( mesh )

            surface = extractSurface( mesh )
            fullSurfaceArea = surfaceArea( surface )

            maxTranslation = 0
            maxRotation = 0
            if args.rigid_displacement:
                maxTranslation=0.02
                maxRotation=5

            surfaceModified, randomTransform = makePartialSurface.makePartialSurface( surface, curSeed,
                    maxTranslation, maxRotation, path=curDir )
            partialSurfaceArea = surfaceArea( surfaceModified )

            maxNoise = random.random()*1e-2
            print("Adding noise of maximum magnitude:",  maxNoise)
            surfaceModified = makePartialSurface.applyNoise( surfaceModified, curSeed, maxNoise )

            writerSurface = vtkXMLPolyDataWriter()
            writerSurface.SetFileName( os.path.join( curDir, "partialSurface.vtp" ) )
            writerSurface.SetInputData( surfaceModified )
            writerSurface.Update()

            voxelize.distanceFieldFromCloud( surfaceModified, grid, "intraoperativeSurface" )

            # Create the ground truth displacement from the deformed mesh:
            mesh = voxelizeDisplacementField.undeform( mesh )

            grid = voxelizeDisplacementField.interpolateToGrid( mesh, grid, args.size/(args.grid_size-1) )
            makePartialSurface.applyRigidTransform( grid, randomTransform )

            visibleAmount = partialSurfaceArea/fullSurfaceArea

            # Make sure we have at least a bit of surface which is visible:
            if visibleAmount < 0.01:
                raise Exception("Less than one percent of surface visible. Ignoring sample")

            visibleSurfaceAmountArr = makeSingleFloatArray( "visibleSurfaceAmount", visibleAmount )
            grid.GetFieldData().AddArray( visibleSurfaceAmountArr )
            print("Visible surface amount:", visibleAmount)

            arraysToKeep = ["preoperativeSurface"]
            arraysToKeep += ["intraoperativeSurface"]
            arraysToKeep += ["displacement"]
            arraysToKeep += ["rigidDisplacement"]
            removeUnwantedArrays( grid.GetCellData(), [] )
            removeUnwantedArrays( grid.GetPointData(), arraysToKeep )

            # Calculate statistics and add them as field data:
            getMaximumOfArray( grid.GetPointData().GetArray("displacement"), grid.GetFieldData())
            getMaximumOfArray( grid.GetPointData().GetArray("rigidDisplacement"), grid.GetFieldData())
            totalDisplacement = sumOfArrays( grid.GetPointData().GetArray("displacement"),
                    grid.GetPointData().GetArray("rigidDisplacement"), "totalDisplacement" )
            print("Max total displacement:", getMaximumOfArray( totalDisplacement, grid.GetFieldData()))

            # Append the max noise value to the field data
            fieldArr = makeSingleFloatArray( "maxNoise", maxNoise )
            grid.GetFieldData().AddArray( fieldArr )

            resultFile = os.path.join( curDir, "voxelized.vts" )
            print("Writing to {}".format( resultFile ))
            writer = vtkXMLStructuredGridWriter()
            writer.SetFileName( resultFile )
            writer.SetInputData( grid )
            writer.Update()

            del totalDisplacement, surface, grid, visibleAmount, visibleSurfaceAmountArr
            del fullSurfaceArea, surfaceModified, partialSurfaceArea


        except KeyboardInterrupt:
            exit()
        except:
            print("\tFailed to voxelize mesh {}. Skipping.".format(i))
            print(sys.exc_info())
            traceback.print_exc()

    gc.collect()

    #snapshot2 = tracemalloc.take_snapshot()
    #top = snapshot2.compare_to(snapshot1,"lineno")
    #for stat in top:
    ##    print(stat)
