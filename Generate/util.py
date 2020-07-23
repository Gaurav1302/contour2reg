# from const import *
import config as cfg
import sys
import torch
import gc
from vtk import *

gridSize = cfg.gridSize
gridSizeInMeters = cfg.gridSizeInMeters
meters2Voxels = gridSize/gridSizeInMeters
voxel2Meters = gridSizeInMeters/gridSize


def unstructuredGridToPolyData( ug ):
    geometryFilter = vtkGeometryFilter()
    geometryFilter.SetInputData( ug )
    geometryFilter.Update()
    return geometryFilter.GetOutput()


# Print iterations progress
def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s %s %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()


def flipTorch(x, dim):
    xsize = x.size()
    dim = x.dim() + dim if dim < 0 else dim
    x = x.view(-1, *xsize[dim:])
    x = x.view(x.size(0), x.size(1), -1)[:, getattr(torch.arange(x.size(1)-1,
                      -1, -1), ('cpu','cuda')[x.is_cuda])().long(), :]
    return x.view(xsize)

def printNumberOfCudaTensors():
    num = 0
    for obj in gc.get_objects():
        try:
            if torch.is_tensor(obj):
                if obj.is_cuda:
                    #print(type(obj), obj.size())
                    num += 1
            elif (hasattr(obj, 'data') and torch.is_tensor(obj.data)):
                if obj.is_cuda or ob.data.is_cuda:
                    #print(type(obj), obj.size())
                    num += 1
        except Exception as ex:
            pass
    print("Cuda Tensors:", num)
