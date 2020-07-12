import vtk
import numpy
import torch
import os
import random
import math
import time
from torch.autograd import Variable
from torch.utils.data import Dataset
from const import *
from util import *
from shutil import copyfile
import traceback

import matplotlib as mpl
mpl.use('Agg')
mpl.use('PS')   # generate postscript output by default

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs

from vtk.util import numpy_support


preloadedSamples = {}

def loadSample( path, i, loadSurfaceArray=False, visPercentage=100, loadMeshes=False ):

    if i in preloadedSamples:
        #print( "Sample {:d} already loaded.".format(i) )
        return preloadedSamples[i]
    #else:
        #print( "Sample {:d} not found. Loading...".format(i) )

    # The input file
    filenameInput = path + "/{:d}/voxels_preprocessed.vtu".format(i)
    #filenameTarget = path + "/{:d}/voxels_result.vtu".format(i)

    samplePath = path + "/{:d}/".format(i)

    if not os.path.isfile(filenameInput):
        return None

    try:

        startTime = time.time()
        # -------------------------------------------
        # Read input:
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filenameInput)
        reader.Update()
        voxelsInput = reader.GetOutput()

        cellData = voxelsInput.GetCellData()
        sdfNumpy = numpy_support.vtk_to_numpy( cellData.GetArray("signedDistance") )
        zeroDisplacementNumpy = numpy_support.vtk_to_numpy( cellData.GetArray("zeroDisplacement") )
        displacementNumpy = numpy_support.vtk_to_numpy( cellData.GetArray("displacement") )

        sdf = numpy.reshape( sdfNumpy, (1, gridSize, gridSize, gridSize) )
        zeroDisplacement = numpy.reshape( zeroDisplacementNumpy, (1, gridSize, gridSize, gridSize) )
        displacement = numpy.reshape( displacementNumpy, (gridSize, gridSize, gridSize, 3) )
        displacement = numpy.transpose( displacement, (3,0,1,2) )

        if loadSurfaceArray:
            surfaceNumpy = numpy_support.vtk_to_numpy( cellData.GetArray("surface") )
            surface = numpy.reshape( surfaceNumpy, (1, gridSize, gridSize, gridSize) )

        #visibleAmounts = [0.25, 0.5, 0.75, 1.0]

        #visibleDisplacement = []
        #for j in range(0,len(visibleAmounts)):
        #    arrName = "visibleDisplacement{:d}".format( int(visibleAmounts[j]*100 ) )
        #    dispNumpy = numpy_support.vtk_to_numpy( cellData.GetArray(arrName) )
        #    disp = numpy.reshape( dispNumpy, (gridSize, gridSize, gridSize, 3) )
        #    disp = numpy.transpose( disp, (3,0,1,2) )
        #    visibleDisplacement.append( disp )
        arrName = "visibleDisplacement{:d}".format( visPercentage )
        try:
            dispNumpy = numpy_support.vtk_to_numpy( cellData.GetArray(arrName) )
        except AttributeError as err:
            print("Could not find array \"{:s}\", will try to load \"visibleDisplacement\" instead.".format( arrName ))
            dispNumpy = numpy_support.vtk_to_numpy( cellData.GetArray("visibleDisplacement") )
        visibleDisplacement = numpy.reshape( dispNumpy, (gridSize, gridSize, gridSize, 3) )
        visibleDisplacement = numpy.transpose( visibleDisplacement, (3,0,1,2) )
        #visibleDisplacement.append( disp )
        arrName = "isVisible{:d}".format( visPercentage )
        try:
            isVisibleNumpy = numpy_support.vtk_to_numpy( cellData.GetArray(arrName) )
        except AttributeError as err:
            print("Could not find array \"{:s}\", will try to load \"isVisible\" instead.".format( arrName ))
            isVisibleNumpy = numpy_support.vtk_to_numpy( cellData.GetArray("isVisible") )

        isVisible = numpy.reshape( isVisibleNumpy, (gridSize, gridSize, gridSize) )


        mask = sdf < 0
        #mask = mask.expand( 3,-1,-1,-1 )
        displacement = displacement*mask

        meshes = []
        #if loadMeshes:
        #    meshesPath = path + "/{:d}/meshes".format(i)
        #    for filename in os.listdir(meshesPath):
        #        if filename.endswith(".vtu"):
        #            print("Loading mesh {} from {}".format( os.path.basename(filename), filename ) )
        #            meshReader = vtk.vtkXMLUnstructuredGridReader()
        #            meshReader.SetFileName(filename)
        #            meshReader.Update()
        #            mesh = meshReader.GetOutput()
        #            meshes.append( { "filename":os.path.basename(filename), "mesh":mesh } )
        #

        # -------------------------------------------
        # Put together training sample:
        sample = {
                "mesh":numpy.concatenate( (sdf, zeroDisplacement) ),
                "visDispl":visibleDisplacement,
                "isVisible":isVisible.astype( numpy.uint8 ),
                #"visDispl25":visibleDisplacement[0],
                #"visDispl50":visibleDisplacement[1],
                #"visDispl75":visibleDisplacement[2],
                #"visDispl100":visibleDisplacement[3],
                "target":displacement,
                "path":samplePath,
                "index":i,
                "meshes":meshes,
                #"visibleAreaSize":visPercentage,
                }
        #print("vis percentage:", visPercentage)
        if loadSurfaceArray:
            sample["surface"] = surface
            print("Loading surface data", surface.shape )

        #print( "Disk load time:", time.time() - startTime )
        #preloadedSamples[i] = sample
        return sample
    except KeyboardInterrupt as err:
        print("KeyboardInterrupt")
        raise( err )
    except Exception as err:
        print(traceback.format_exc())
        print(err)
        return None


def loadSimData( path, num, startIndex=0, maxDefThreshold=0.1 ):

    data = []
    for i in range( startIndex, num + startIndex ):
        sample = loadSample( path, i )
        if sample is not None:
            maxDeformation = numpy.amax( numpy.linalg.norm( sample["target"], ord=2, axis=0 ) )
            #maxForce = numpy.amax( numpy.linalg.norm( sample["force"], ord=2, axis=0 ) )
            if maxDeformation < maxDefThreshold:# and not ( maxDeformation > 0 and maxForce < 1e-10 ):
                data.append( sample )
            #else:
                #print("Warning: Skipping sample {:d} due to large deformation:\n\tDeformation: {:e}, Threshold: {:e}".format(i,maxDeformation,maxDefThreshold))
        else:
            print("Warning: Could not load sample {:d}!".format(i))

        printProgressBar(i+1-startIndex, num, prefix='Loading', suffix='', decimals=2)

    print("Loaded {:d} samples. Skipped {:d} samples.".format(len(data),num-len(data)) )
    return data

# Remove any samples where the deformation is above a certain threshold
def removeLargeDeformations( data, threshold ):
    data[:] = [x for x in data if numpy.amax( numpy.linalg.norm( x["target"], ord=2, axis=0 ) ) < threshold ]
    return data

# Check if the sample has a deformation but no force. If so, return True. This sample will be discarded.
def deformationMissingForce( sample ):
    maxDeformation = numpy.amax( numpy.linalg.norm( sample["target"], ord=2, axis=0 ) )
    maxForce = numpy.amax( numpy.linalg.norm( sample["force"], ord=2, axis=0 ) )
    if maxDeformation > 0 and maxForce < 1e-10:
        return True
    return False

# If there is no force, but there is deformation, remove the sample. This can happen
# if there's an error in discretization (voxelization), for example:
def removeZeroDeformations( data ):
    data[:] = [x for x in data if not deformationMissingForce(x)]
    return data

def flipSample( d, dim ):
    f = {}
    if type(d["mesh"]) == numpy.ndarray:    # Numpy tensors?
        # Flip data along an axis:
        f["mesh"] = numpy.ascontiguousarray( numpy.flip( d["mesh"], dim+1 ) )
        #f["force"] = numpy.ascontiguousarray( numpy.flip( d["force"], dim+1 ) )
        f["visDispl"] = numpy.ascontiguousarray( numpy.flip( d["visDispl"], dim+1 ) )
        #f["visDispl25"] = numpy.ascontiguousarray( numpy.flip( d["visDispl25"], dim+1 ) )
        #f["visDispl50"] = numpy.ascontiguousarray( numpy.flip( d["visDispl50"], dim+1 ) )
        #f["visDispl75"] = numpy.ascontiguousarray( numpy.flip( d["visDispl75"], dim+1 ) )
        #f["visDispl100"] = numpy.ascontiguousarray( numpy.flip( d["visDispl100"], dim+1 ) )
        f["target"] = numpy.ascontiguousarray( numpy.flip( d["target"], dim+1 ) )
        f["isVisible"] = numpy.ascontiguousarray( numpy.flip( d["isVisible"], dim ) )
    else:   # Torch tensors
        # Flip data along an axis:
        f["mesh"] = flipTorch( d["mesh"], dim+1 )
        f["visDispl"] = flipTorch( d["visDispl"], dim+1 )
        #f["visDispl25"] = flipTorch( d["visDispl25"], dim+1 )
        #f["visDispl50"] = flipTorch( d["visDispl50"], dim+1 )
        #f["visDispl75"] = flipTorch( d["visDispl75"], dim+1 )
        #f["visDispl100"] = flipTorch( d["visDispl100"], dim+1 )
        f["target"] = flipTorch( d["target"], dim+1 )
        f["isVisible"] = flipTorch( d["isVisible"], dim )


    # Flip sign of the vector component of force and mesh along the flipped axis:
    f["visDispl"][dim,:,:,:] = -f["visDispl"][dim,:,:,:]
    #f["visDispl25"][dim,:,:,:] = -f["visDispl25"][dim,:,:,:]
    #f["visDispl50"][dim,:,:,:] = -f["visDispl50"][dim,:,:,:]
    #f["visDispl75"][dim,:,:,:] = -f["visDispl75"][dim,:,:,:]
    #f["visDispl100"][dim,:,:,:] = -f["visDispl100"][dim,:,:,:]
    f["target"][dim,:,:,:] = -f["target"][dim,:,:,:]

    # Copy path:
    # Note: Saving will not work on augmented samples at this point!
    f["path"] = d["path"]
    f["index"] = d["index"]
    f["meshes"] = d["meshes"]
    return f


def augment( data ):
    augmented = []
    for d in data:
        augmented += augmentSample( d )

    return augmented

# TODO check gravity. If it is applied, we should not flip in the direction in which gravity pulls!
def augmentSample( d ):
    # Flip tensors:
    x = flipSample( d, 0 )
    y = flipSample( d, 1 )
    z = flipSample( d, 2 )
    xy = flipSample( x, 1 )
    xz = flipSample( x, 2 )
    yz = flipSample( y, 2 )
    xyz = flipSample( xy, 2 )

    augmented = [d,x,y,z,xy,xz,yz,xyz]
    return augmented


def sampleToTorch( d, cuda=True ):
    sample = {}
    if cuda:
        for v in ["mesh", "visDispl25", "visDispl50", "visDispl75", "visDispl100", "visDispl", "target", "surface" ]:
            if v in d:
                startTime = time.time()
                sample[v] = torch.Tensor( d[v] ).cuda()
                #if v == "visDispl":
                    #print("Upload time:", time.time() - startTime )
        sample["isVisible"] = torch.Tensor( d["isVisible"] ).cuda()
    else:
        for v in ["mesh", "visDispl25", "visDispl50", "visDispl75", "visDispl100", "visDispl", "target", "surface" ]:
            if v in d:
                sample[v] = torch.Tensor( d[v] )
        sample["isVisible"] = torch.Tensor( d["isVisible"] )

    sample["path"] = d["path"]
    sample["index"] = d["index"]
    sample["meshes"] = d["meshes"]
    if "visibleAreaSize" in d:
        sample["visibleAreaSize"] = d["visibleAreaSize"]
    return sample


def toTorch( data, cuda=True ):
    dataTorch = []
    for d in data:
        sample = sampleToTorch( d )
        dataTorch.append( sample )

    return dataTorch

def saveResults( data, outPath, startIndex=0 ):
    errors = []
    errorsRel = []

    binsDist = 30
    binsDeform = 10
    maxDist = 0.3
    maxDeform = 0.1
    for i in range( 0, binsDist ):
        errorsDef = []
        errorsDefRel = []
        for j in range( 0, binsDeform ):
            errorsDef.append( [] )
            errorsDefRel.append( [] )
        errors.append( errorsDef )
        errorsRel.append( errorsDefRel )

    metricsPerVoxel = []

    max_0_0_err = 0
    max_0_0_err_index = -1

    reader = vtk.vtkXMLUnstructuredGridReader()
    writer = vtk.vtkXMLUnstructuredGridWriter()
    for i in range( 0, len(data) ):

        d = data[i]

        origIndex = d["index"]
        print("Processing result for sample:",origIndex )
        #random.seed(origIndex)
        #for j in range( 0, 10 ):
            #random.randint()

        reader.SetFileName( d["path"] + "/voxels_preprocessed.vtu" )
        reader.Update()
        voxelsTarget = reader.GetOutput()

        cellData = voxelsTarget.GetCellData()
        if cellData.GetArray("displacement") is not None:
            if cellData.GetArray("displacementOrig") is None:
                cellData.GetArray("displacement").SetName("displacementOrig")

        target = d["target"].cpu().numpy().squeeze()
        targetMag = numpy.linalg.norm( target, 2, 0 ).squeeze()

        output = d["output"].squeeze()
        outputMag = numpy.linalg.norm( output, 2, 0 ).squeeze()

        meshNumpy = d["mesh"].cpu().numpy().squeeze()
        mesh = numpy.reshape( meshNumpy[0,:,:,:], (gridSize, gridSize, gridSize) )

        isVisible = d["isVisible"].cpu()
        isVisibleVoxelList = numpy.transpose( numpy.nonzero( isVisible ) )
        #print( isVisibleVoxelList )

        if isVisibleVoxelList.shape[1] > 0:

          # For every voxel, calculate the distance to visible surface:
          dist2Visible = distanceFunction( mesh, isVisibleVoxelList )

          for x in range( 0, gridSize ):
              for y in range( 0, gridSize ):
                  for z in range( 0, gridSize ):
                      if mesh[x,y,z] < 0:
                          dist = dist2Visible[x,y,z]
                          err = d["displacementError"][x,y,z]
                          errRel = d["displacementErrorRel"][x,y,z]
                          displTargetMag = targetMag[x,y,z]
                          displOutputMag = outputMag[x,y,z]
                          distBin = int(dist/maxDist*binsDist)
                          displTargetBin = int(displTargetMag/maxDeform*binsDeform)

                          #voxelMetrics = dict()
                          #voxelMetrics["err"] = err
                          #voxelMetrics["errRel"] = errRel
                          #voxelMetrics["displacement"] = displ
                          #voxelMetrics["dist2visible"] = dist
                          metricsPerVoxel.append( (displOutputMag,displTargetMag,dist,err) )
                          #print(dist,displ,err)
                          #print("\t",distBin,displBin)
                          if distBin >= 0 and distBin < binsDist and displTargetBin >= 0 and displTargetBin < binsDeform:
                              errors[distBin][displTargetBin].append( err )
                              errorsRel[distBin][displTargetBin].append( errRel )
                              if distBin == 0 and displTargetBin == 0:
                                  if err > max_0_0_err:
                                      max_0_0_err = err
                                      max_0_0_err_index = i

          d["dist2Visible"] = dist2Visible

        write = True
        if write:

            if "output" in d:
                disp = numpy.transpose( d["output"], (1,2,3,0) )
                disp = numpy.ravel( disp )
                vtkArray = numpy_support.numpy_to_vtk( disp, deep=True, array_type=vtk.VTK_DOUBLE )
                vtkArray.SetName("displacement")
                vtkArray.SetNumberOfComponents( 3 )
                vtkArray.SetNumberOfTuples( gridSize*gridSize*gridSize )
                cellData.AddArray( vtkArray )
            if "displacementError" in d:
                #disp = numpy.transpose( d["displacementError"], (1,2,3,0) )
                disp = numpy.ravel( d["displacementError"] )
                vtkArray = numpy_support.numpy_to_vtk( disp, deep=True, array_type=vtk.VTK_DOUBLE )
                vtkArray.SetName("displacementError")
                vtkArray.SetNumberOfComponents( 1 )
                vtkArray.SetNumberOfTuples( gridSize*gridSize*gridSize )
                cellData.AddArray( vtkArray )
            if "displacementErrorRel" in d:
                #disp = numpy.transpose( d["displacementError"], (1,2,3,0) )
                disp = numpy.ravel( d["displacementErrorRel"] )
                vtkArray = numpy_support.numpy_to_vtk( disp, deep=True, array_type=vtk.VTK_DOUBLE )
                vtkArray.SetName("displacementErrorRel")
                vtkArray.SetNumberOfComponents( 1 )
                vtkArray.SetNumberOfTuples( gridSize*gridSize*gridSize )
                cellData.AddArray( vtkArray )
            if "visDispl" in d:
                disp = numpy.transpose( d["visDispl"], (1,2,3,0) )
                disp = numpy.ravel( disp )
                vtkArray = numpy_support.numpy_to_vtk( disp, deep=True, array_type=vtk.VTK_DOUBLE )
                vtkArray.SetName("visDispl")
                vtkArray.SetNumberOfComponents( 3 )
                vtkArray.SetNumberOfTuples( gridSize*gridSize*gridSize )
                cellData.AddArray( vtkArray )
            if "dist2Visible" in d:
                #disp = numpy.transpose( d["displacementError"], (1,2,3,0) )
                disp = numpy.ravel( d["dist2Visible"] )
                vtkArray = numpy_support.numpy_to_vtk( disp, deep=True, array_type=vtk.VTK_DOUBLE )
                vtkArray.SetName("dist2Visible")
                vtkArray.SetNumberOfComponents( 1 )
                vtkArray.SetNumberOfTuples( gridSize*gridSize*gridSize )
                cellData.AddArray( vtkArray )
            if "zeroDisplacement" in d:
                #disp = numpy.transpose( d["displacementError"], (1,2,3,0) )
                disp = numpy.ravel( d["zeroDisplacement"] )
                vtkArray = numpy_support.numpy_to_vtk( disp, deep=True, array_type=vtk.VTK_DOUBLE )
                vtkArray.SetName("zeroDisplacement")
                vtkArray.SetNumberOfComponents( 1 )
                vtkArray.SetNumberOfTuples( gridSize*gridSize*gridSize )
                cellData.AddArray( vtkArray )
            if "tmpSurface" in d:
                #disp = numpy.transpose( d["displacementError"], (1,2,3,0) )
                disp = numpy.ravel( d["tmpSurface"] )
                vtkArray = numpy_support.numpy_to_vtk( disp, deep=True, array_type=vtk.VTK_DOUBLE )
                vtkArray.SetName("tmpSurface")
                vtkArray.SetNumberOfComponents( 1 )
                vtkArray.SetNumberOfTuples( gridSize*gridSize*gridSize )
                cellData.AddArray( vtkArray )

            if "visibleAreaSize" in d:
                print("Writing visible area:", d["visibleAreaSize"])
                vtkArray = vtk.vtkIntArray()
                vtkArray.SetNumberOfTuples(1)
                vtkArray.SetNumberOfComponents(1)
                vtkArray.SetName( "visibleAreaSize" )
                vtkArray.SetTuple1( 0, d["visibleAreaSize"] )
                voxelsTarget.GetFieldData().AddArray( vtkArray )



            folder = outPath + "/{:d}/".format(d["index"])
            if not os.path.exists(folder):
                os.makedirs(folder)
            filename = folder + "voxels_output.vtu"
            if "suffix" in d:
                filename = folder + "voxels_output" + str(d["suffix"]) + ".vtu"
            print("\t", filename)
            writer.SetFileName( filename )
            writer.SetInputData( voxelsTarget )
            writer.Write()
            print("Saved:", folder)


    print("Found the maximum error for a voxel with 0 distance and 0 deformation in sample: {:d} (Error: {:e})".format( max_0_0_err_index, max_0_0_err ) )
        #copyfile( d["path"] + "/voxels_input.vtu", folder + "/voxels_input.vtu")
        #copyfile( d["path"] + "/voxels_result.vtu", folder + "/voxels_result.vtu")
        #copyfile( d["path"] + "/voxels_result.vtu", folder + "/voxels_result.vtu")

    #dists = numpy.linspace( 0, 0.3, 30 )
    minErr = numpy.zeros( (binsDist, binsDeform) )
    maxErr = numpy.zeros( (binsDist, binsDeform) )
    meanErr = numpy.zeros( (binsDist, binsDeform) )
    medianErr = numpy.zeros( (binsDist, binsDeform) )
    minErrRel = numpy.zeros( (binsDist, binsDeform) )
    maxErrRel = numpy.zeros( (binsDist, binsDeform) )
    meanErrRel = numpy.zeros( (binsDist, binsDeform) )
    medianErrRel = numpy.zeros( (binsDist, binsDeform) )

    samplesPerBin = numpy.zeros( (binsDist, binsDeform) )
    for i in range( 0, binsDist ):
        for j in range( 0, binsDeform ):
            arr = numpy.array( errors[i][j] )
            #print(arr.shape, arr.min(), arr.max() )
            if arr.shape[0] > 0:
                minErr[i,j] = arr.min()
                maxErr[i,j] = arr.max()
                meanErr[i,j] = arr.mean()
                medianErr[i,j] = numpy.median( arr )

            arrRel = numpy.array( errorsRel[i][j] )
            if arrRel.shape[0] > 0:
                minErrRel[i,j] = arrRel.min()
                maxErrRel[i,j] = arrRel.max()
                meanErrRel[i,j] = arrRel.mean()
                medianErrRel[i,j] = numpy.median( arrRel )

            samplesPerBin[i,j] = len( errors[i][j] )

    numpy.save( outPath + "/errorMetricsMin", minErr )
    numpy.save( outPath + "/errorMetricsMax", maxErr )
    numpy.save( outPath + "/errorMetricsMean", meanErr )
    numpy.save( outPath + "/errorMetricsMedian", medianErr )

    numpy.save( outPath + "/errorMetricsMinRel", minErrRel )
    numpy.save( outPath + "/errorMetricsMaxRel", maxErrRel )
    numpy.save( outPath + "/errorMetricsMeanRel", meanErrRel )
    numpy.save( outPath + "/errorMetricsMedianRel", medianErrRel )
    numpy.save( outPath + "/metrics", metricsPerVoxel )


    fig = plt.figure()

    plt.clf()
    ax = fig.gca()
    #im = ax.imshow( samplesPerBin, norm=LogNorm(vmin=0.01,vmax=1))
    im = ax.imshow( samplesPerBin, norm=LogNorm())
    ax.set_ylabel("dist from visible (cm)")
    ax.set_xlabel("deformation (cm)")


    #axcolor = f.add_axes([0.90, 0.02, 0.03, 0.79])
    #im = ax.matshow(C, cmap=cm.gray_r, norm=LogNorm(vmin=0.01, vmax=1))
    ticks = []
    t = 1
    while t < samplesPerBin.max():
        ticks.append( t )
        t *= 10
    #tickLabels
    cbar = fig.colorbar( im, format='$%.0f$', ticks=ticks)
    #cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented colorbar
    plt.savefig( outPath + "/samplesPerBin" )



def overwriteArray( arrayName, sample, cellData ):
    print("Adding " + arrayName, cellData.GetNumberOfArrays() )
    if cellData.GetArray( arrayName ) is not None:
        cellData.RemoveArray( arrayName )
        print("\tremoved old")
    arr = numpy.transpose( sample[arrayName].data.cpu().numpy(), (1,2,3,0) )
    print(arr.shape, arr.shape[3])
    if arr.shape[3] == 1 or arr.shape[3] == 3:
        arr = numpy.ravel( arr )
        vtkArray = numpy_support.numpy_to_vtk( arr, deep=True, array_type=vtk.VTK_DOUBLE )
        vtkArray.SetName(arrayName)
        vtkArray.SetNumberOfComponents( sample[arrayName].size(0) )
        vtkArray.SetNumberOfTuples( gridSize*gridSize*gridSize )
        cellData.AddArray( vtkArray )
        print("\tAdded " + vtkArray.GetName(), cellData.GetNumberOfArrays() )
    elif arr.shape[3] == 2:
        arr0 = numpy.ravel( arr[:,:,:,0] )
        vtkArray0 = numpy_support.numpy_to_vtk( arr0, deep=True, array_type=vtk.VTK_DOUBLE )
        vtkArray0.SetName(arrayName + "0")
        vtkArray0.SetNumberOfComponents( 1 )
        vtkArray0.SetNumberOfTuples( gridSize*gridSize*gridSize )
        cellData.AddArray( vtkArray0 )
        print("\tAdded " + vtkArray0.GetName(), cellData.GetNumberOfArrays() )
        arr1 = numpy.ravel( arr[:,:,:,1] )
        vtkArray1 = numpy_support.numpy_to_vtk( arr1, deep=True, array_type=vtk.VTK_DOUBLE )
        vtkArray1.SetName(arrayName + "1")
        vtkArray1.SetNumberOfComponents( 1 )
        vtkArray1.SetNumberOfTuples( gridSize*gridSize*gridSize )
        cellData.AddArray( vtkArray1 )
        print("\tAdded " + vtkArray1.GetName(), cellData.GetNumberOfArrays() )



def save( data, outPath, startIndex=0 ):
    reader = vtk.vtkXMLUnstructuredGridReader()
    writer = vtk.vtkXMLUnstructuredGridWriter()
    for i in range( 0, len(data) ):
        d = data[i]

        reader.SetFileName( d["path"] + "/voxels_preprocessed.vtu" )
        reader.Update()
        voxelsTarget = reader.GetOutput()

        cellData = voxelsTarget.GetCellData()

        if "mesh" in d:
            overwriteArray( "mesh", d, cellData )
        if "force" in d:
            overwriteArray( "force", d, cellData )
        if "target" in d:
            overwriteArray( "target", d, cellData )

        folder = outPath + "/{:d}/".format(i+startIndex)
        if not os.path.exists(folder):
            os.makedirs(folder)
        writer.SetFileName( folder + "voxels_output.vtu" )
        writer.SetInputData( voxelsTarget )
        writer.Write()


class SimulatedDataSet( Dataset ):

    def __init__( self, path, num, startIndex=0, maxDefThreshold=0.1, augment=True,
            visPercentage=-1, loadSurfaceArray=False, loadMeshes=False ):

        self.validSamples = []

        self.drawingWeights = []

        #self.drawn = []

        visibleDisplacementCounter = {
            100: 0,
            75: 0,
            50: 0,
            25: 0,
        }

        # Preload data and check for valid samples:
        for i in range( startIndex, num + startIndex ):
            printProgressBar( i+1-startIndex, num, "Checking Data", decimals=2 )

            if visPercentage == -1:
                if i % 4 == 0:
                    visP = 100
                elif i % 4 == 1:
                    visP = 75
                elif i % 4 == 2:
                    visP = 50
                else:
                    visP = 25
            sample = loadSample( path, i, loadSurfaceArray, visP, loadMeshes=loadMeshes )
            if sample is not None:
                maxDeformation = numpy.amax( numpy.linalg.norm( sample["target"], ord=2, axis=0 ) )
                #maxForce = numpy.amax( numpy.linalg.norm( sample["force"], ord=2, axis=0 ) )
                # If this sample does not have a deformation above maxDefThreshold and
                # it does not have a deformation without a force (could happen due to bad voxelization)
                if maxDeformation < maxDefThreshold: #and not ( maxDeformation > 0 and maxForce < 1e-10 ):
                    # ... then remember this is a valid sample.
                    self.validSamples.append(i)
                    # Calculate the probability of this sample to be drawn:
                    drawingWeight = max( maxDeformation,1e-3 )
                    # Store the drawing weight. If augmented, we need the weight 8 times
                    if not augment:
                        self.drawingWeights.append(drawingWeight)
                    else:
                        for i in range(0,8):
                            self.drawingWeights.append(drawingWeight)

                visibleDisplacementCounter[visP] += 1

        print( "Number of Samples by visibleDisplacement:\n", visibleDisplacementCounter )

        self.path = path
        self.augment = augment
        self.visPercentage = visPercentage
        self.loadMeshes = loadMeshes

    def __len__(self):
        if self.augment:
            return len( self.validSamples*8 )
        else:
            return len( self.validSamples )

    def __getitem__( self, i ):

        try:
            sample = None
            baseIndex = None
            weight = None
            if self.augment:
                baseIndex = self.validSamples[int(i/8)]
                weight = self.drawingWeights[int(i/8)]
            else:
                baseIndex = self.validSamples[i]
                weight = self.drawingWeights[i]

            #self.drawn.append(baseIndex)
            visPercentage = self.visPercentage
            if visPercentage == -1:
                if baseIndex % 4 == 0:
                    visPercentage = 100
                elif baseIndex % 4 == 1:
                    visPercentage = 75
                elif baseIndex % 4 == 2:
                    visPercentage = 50
                else:
                    visPercentage = 25

            sample = loadSample( self.path, baseIndex, visPercentage=visPercentage, loadMeshes=self.loadMeshes )
            #sample = sampleToTorch( sample, cuda=True )
            sample = sampleToTorch( sample, cuda=False )

            if self.augment:

                augmentedIndex = i % 8

                if augmentedIndex == 1:   # Flip x
                    sample = flipSample( sample, 0 )
                elif augmentedIndex == 2:   # Flip y
                    sample = flipSample( sample, 1 )
                elif augmentedIndex == 3:   # Flip z
                    sample = flipSample( sample, 2 )
                elif augmentedIndex == 4:   # Flip xy
                    sample = flipSample( flipSample(sample,0), 1 )
                elif augmentedIndex == 5:   # Flip xz
                    sample = flipSample( flipSample(sample,0), 2 )
                elif augmentedIndex == 6:   # Flip yz
                    sample = flipSample( flipSample(sample,1), 2 )
                elif augmentedIndex == 7:   # Flip xyz
                    sample = flipSample( flipSample( flipSample(sample,0),1 ),2 )
            return sample
            #return sampleToTorch( sample, cuda=True )
        except KeyboardInterrupt as err:
            raise( err )
        except Exception as err:
            print(traceback.format_exc())
            print(err)
            newIndex = random.randint( 0, len( self )-1 )
            print( "Warning. There was an error while loading sample {:d}. Will try loading a random sample ({:d}) instead.".format( i, newIndex ) )
            return self[newIndex]


##########################################################
# Generate a displacement boundary condition:

def isBoundary( sdf, x, y, z ):
    if sdf[x,y,z] > 0:
        return False

    if x == 0 or x == gridSize-1:
        return True
    if y == 0 or y == gridSize-1:
        return True
    if z == 0 or z == gridSize-1:
        return True

    if sdf[x+1,y,z] > 0 or sdf[x-1,y,z] > 0:
        return True
    if sdf[x,y+1,z] > 0 or sdf[x,y-1,z] > 0:
        return True
    if sdf[x,y,z+1] > 0 or sdf[x,y,z-1] > 0:
        return True
    return False

def checkIfDistIsSmaller( distance, v, vNext, offset ):
    if vNext[0] < 0 or vNext[0] > gridSize-1:
        return False
    if vNext[1] < 0 or vNext[1] > gridSize-1:
        return False
    if vNext[2] < 0 or vNext[2] > gridSize-1:
        return False

    distAtSource = distance[v]
    distAtTarget = distance[vNext]
    if distAtSource + offset < distAtTarget:
        distance[vNext[0],vNext[1],vNext[2]] = distAtSource + offset
        return True
    return False

def signedDistanceFunction( mesh, surface ):

    distance = numpy.empty( mesh.shape, numpy.float )
    distance.fill( numpy.inf )
    front = []
    for voxel in surface:
        v = (voxel[0],voxel[1],voxel[2])
        front.append( v )
        distance[ v ] = 0
    offsetStraight = 1
    offsetDiagonal = math.sqrt( 2 )

    while len( front ) > 0:
        v = front.pop(0)
        for dx in range(-1,2):
            for dy in range(-1,2):
                for dz in range(-1,2):
                    if not (dx==0 and dy==0 and dz==0):
                        vNext = ( v[0]+dx, v[1]+dy, v[2]+dz )
                        offset = offsetStraight
                        if abs(dx)+abs(dy)+abs(dz) > 1:
                            offset = offsetDiagonal
                        if checkIfDistIsSmaller( distance, v, vNext, offset ):
                            front.append( vNext )

    distance[mesh==0] *= -1
    return distance/(2*gridSize)

def validVoxel( x, y, z, size ):
    if x < 0 or y < 0 or z < 0 or \
        x >= size[0] or y >= size[1] or z >= size[2]:
        return False
    return True

def getNeighbours( x, y, z ):
    return [(x-1,y,z),(x+1,y,z),(x,y-1,z),(x,y+1,z),(x,y,z-1),(x,y,z+1)]
def getNeighboursDiag( x, y, z ):
    return [(x-1,y-1,z),(x+1,y-1,z),(x-1,y+1,z),(x+1,y+1,z),
            (x,y-1,z-1),(x,y+1,z-1),(x,y-1,z+1),(x,y+1,z+1),
            (x-1,y,z-1),(x+1,y,z-1),(x-1,y,z+1),(x+1,y,z+1)]
def getNeighboursDiag3( x, y, z ):
    return [(x-1,y-1,z-1),(x+1,y-1,z-1),(x-1,y+1,z-1),(x-1,y-1,z+1),
            (x+1,y+1,z+1),(x+1,y+1,z-1),(x-1,y+1,z+1),(x+1,y-1,z+1)]

def dist2( dx,dy,dz ):
    return dx*dx + dy*dy + dz*dz

def distanceFunction( mesh, surface ):
    distance = numpy.empty( mesh.shape, numpy.float )
    distance.fill( -1 )

    front = []
    for voxel in surface:
        v = (voxel[0].item(),voxel[1].item(),voxel[2].item())
        front.append( v )

    for x in range(0,distance.shape[0]):
        for y in range(0,distance.shape[1]):
            for z in range(0,distance.shape[2]):
                minDist = 9999999
                for v in front:
                    d2 = dist2( v[0]-x,v[1]-y,v[2]-z )
                    if d2 < minDist:
                        minDist = d2
                distance[x,y,z] = math.sqrt( minDist )


    return distance*voxel2Meters

def distanceFunctionBkup2( mesh, surface ):
    distance = numpy.empty( mesh.shape, numpy.float )
    distance.fill( -1 )
    front = []
    for voxel in surface:
        v = (voxel[0],voxel[1],voxel[2])
        front.append( v )

    offsetStraight = 1
    offsetDiagonal = math.sqrt(2)
    offsetDiagonal3 = math.sqrt(3)

    gridSize = distance.shape

    while len( front ) > 0:
        newFront = []
        for v in front:
            vDist = distance[ v ]
            n1 = getNeighbours( v[0],v[1],v[2] )
            n2 = getNeighboursDiag( v[0],v[1],v[2] )
            n3 = getNeighboursDiag3( v[0],v[1],v[2] )
            for n in n1:
                if validVoxel( n[0],n[1],n[2], gridSize ):
                    if distance[n] == -1:
                        distance[n] = vDist + offsetStraight
                        newFront.append( n )
            for n in n2:
                if validVoxel( n[0],n[1],n[2], gridSize ):
                    if distance[n] == -1:
                        distance[n] = vDist + offsetDiagonal
                        newFront.append( n )
            for n in n3:
                if validVoxel( n[0],n[1],n[2], gridSize ):
                    if distance[n] == -1:
                        distance[n] = vDist + offsetDiagonal3
                        newFront.append( n )

        front = newFront


def distanceFunctionBkup( mesh, surface ):

    distance = numpy.empty( mesh.shape, numpy.float )
    distance.fill( numpy.inf )
    front = []
    for voxel in surface:
        v = (voxel[0],voxel[1],voxel[2])
        front.append( v )
        distance[ v ] = 0
    offsetStraight = 1
    offsetDiagonal = math.sqrt( 2 )

    while len( front ) > 0:
        v = front.pop(0)
        for dx in range(-1,2):
            for dy in range(-1,2):
                for dz in range(-1,2):
                    if not (dx==0 and dy==0 and dz==0):
                        vNext = ( v[0]+dx, v[1]+dy, v[2]+dz )
                        offset = offsetStraight
                        if abs(dx)+abs(dy)+abs(dz) > 1:
                            offset = offsetDiagonal
                        if checkIfDistIsSmaller( distance, v, vNext, offset ):
                            front.append( vNext )

    return distance*voxel2Meters


def preprocess( path, i ):

    try:
        #filenameInput = path + "/{:d}/voxels_input.vtu".format(i)
        filenameTarget = path + "/{:d}/voxels_result.vtu".format(i)
        samplePath = path + "/{:d}/".format(i)

        #if not os.path.isfile(filenameInput) or not os.path.isfile(filenameTarget):
        if not os.path.isfile(filenameTarget):
            return None

        random.seed( i )

        # -------------------------------------------
        # Load the original input data:
        #reader = vtk.vtkXMLUnstructuredGridReader()
        #reader.SetFileName(filenameInput)
        #reader.Update()
        #voxelsInput = reader.GetOutput()

        #cellData = voxelsInput.GetCellData()
        #meshNumpy = numpy_support.vtk_to_numpy( cellData.GetArray("mesh") )
        #mesh = numpy.reshape( meshNumpy, (gridSize, gridSize, gridSize) )

        # -------------------------------------------
        # Read target:
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filenameTarget)
        reader.Update()
        voxels = reader.GetOutput()

        cellData = voxels.GetCellData()

        displacementNumpy = numpy_support.vtk_to_numpy( cellData.GetArray("displacement") )
        displacement = numpy.reshape( displacementNumpy, (gridSize, gridSize, gridSize, 3) )

        sdfNumpy = numpy_support.vtk_to_numpy( cellData.GetArray("signedDistance") )
        sdf = numpy.reshape( sdfNumpy, (gridSize, gridSize, gridSize) )

        # Find all the voxels which are visible from the outside (i.e. boundary voxels):
        visibleVoxels = numpy.zeros( sdf.shape )
        for x in range(0,gridSize):
            for y in range(0,gridSize):
                for z in range(0,gridSize):
                    #if abs(sdf[x,y,z]) <= voxel2Meters:
                    if isBoundary( sdf, x, y, z ):
                        visibleVoxels[x,y,z] = 1
        visibleVoxelsList = numpy.transpose(numpy.nonzero( visibleVoxels ))

        # Add the list of visible voxels:
        vtkArray = numpy_support.numpy_to_vtk( numpy.ravel(visibleVoxels), deep=True, array_type=vtk.VTK_DOUBLE )
        vtkArray.SetName("visibleVoxels")
        vtkArray.SetNumberOfComponents( 1 )
        vtkArray.SetNumberOfTuples( gridSize*gridSize*gridSize )
        cellData.AddArray( vtkArray )

        # Add the sdf:
        #sdf = signedDistanceFunction( mesh, visibleVoxelsList )
        #vtkArray = numpy_support.numpy_to_vtk( numpy.ravel(sdf), deep=True, array_type=vtk.VTK_DOUBLE )
        #vtkArray.SetName("SDF")
        #vtkArray.SetNumberOfComponents( 1 )
        #vtkArray.SetNumberOfTuples( gridSize*gridSize*gridSize )
        #cellData.AddArray( vtkArray )

        # Random sphere:
        center = random.randint( 0, visibleVoxelsList.shape[0]-1 )
        centerPoint = visibleVoxelsList[center]
        #radius = random.random()*0.05 + 0.03
        #radiusInVoxels = radius*meters2Voxels

        # Calculate how far away the furthest surface voxel is from the sphere center:
        maxDistFromSphereCenter = 0
        for surfacePoint in visibleVoxelsList:
            dist = numpy.linalg.norm( surfacePoint - centerPoint )
            maxDistFromSphereCenter = max( dist, maxDistFromSphereCenter )

        # List of sphere sizes:
        visibleAmounts = [0.25, 0.5, 0.75, 1.0]

        visibleDisplacement = []
        for i in range(0,len(visibleAmounts)):
            visibleDisplacement.append( numpy.zeros( displacement.shape ) )

        isVisible = []
        for i in range(0,len(visibleAmounts)):
            isVisible.append( numpy.zeros( (gridSize, gridSize, gridSize), dtype=bool ) )


        for point in visibleVoxelsList:
            dist = numpy.linalg.norm( point - centerPoint, 2 )

            # For each amount, check if the distance is small enough to be added to the visible sphere:
            for i in range(0,len(visibleAmounts)):
                if dist <= visibleAmounts[i]*maxDistFromSphereCenter:
                    visibleDisplacement[i][point[0],point[1],point[2],:] = \
                            displacement[point[0],point[1],point[2],:]
                    isVisible[i][point[0],point[1],point[2]] = True


        for i in range(0,len(visibleAmounts)):
            disp = numpy.ravel( visibleDisplacement[i] )
            vtkArray = numpy_support.numpy_to_vtk( disp, deep=True, array_type=vtk.VTK_DOUBLE )
            vtkArray.SetName("visibleDisplacement{:d}".format( int(visibleAmounts[i]*100) ) )
            vtkArray.SetNumberOfComponents( 3 )
            vtkArray.SetNumberOfTuples( gridSize*gridSize*gridSize )
            cellData.AddArray( vtkArray )

            isVis = numpy.ravel( isVisible[i] )
            vtkArray = numpy_support.numpy_to_vtk( isVis, deep=True, array_type=vtk.VTK_CHAR )
            vtkArray.SetName("isVisible{:d}".format( int(visibleAmounts[i]*100) ) )
            vtkArray.SetNumberOfComponents( 1 )
            vtkArray.SetNumberOfTuples( gridSize*gridSize*gridSize )
            cellData.AddArray( vtkArray )



        arraysToKeep = ["displacement", "signedDistance", "visibleVoxels", "force", "zeroDisplacement"]
        for i in range(0,len(visibleAmounts)):
            arraysToKeep.append("visibleDisplacement{:d}".format( int(visibleAmounts[i]*100) ) )
            arraysToKeep.append("isVisible{:d}".format( int(visibleAmounts[i]*100) ) )

        # Remove any unused arrays:
        for i in range( cellData.GetNumberOfArrays(), 0, -1 ):
            name = cellData.GetArrayName(i-1)
            if not name in arraysToKeep:
                cellData.RemoveArray(i-1)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName( samplePath + "/voxels_preprocessed.vtu" )
        writer.SetInputData( voxels )
        writer.Write()

        return True
    except KeyboardInterrupt as err:
        print("KeyboardInterrupt")
        raise( err )
    #except Exception as err:
    #    print(err)
    #    return False
