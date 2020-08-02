import numpy as np
import vtk
from vtk import *
from vtk.util.numpy_support import vtk_to_numpy
from vtk import vtkStructuredPointsReader
import matplotlib
from util import unstructuredGridToPolyData
matplotlib.use('agg')
import matplotlib.pyplot as plt
import cv2
import os
import config as cfg

def generate_contour(data_path, id):

    colors = vtk.vtkNamedColors()

    # ren = vtk.vtkRenderer()
    # renWin = vtk.vtkRenderWindow()
    # renWin.AddRenderer(ren)

    # # create a renderwindowinteractor
    # iren = vtk.vtkRenderWindowInteractor()
    # iren.SetRenderWindow(renWin)

    # read .vtu file
    simulationResultFile = data_path + "/" + str(id) + "/simulation/result/case_t0001.vtu"
    mesh_2d_path = data_path + "/" + str(id) + "/mesh_2d.png"
    contour_path = data_path + "/" + str(id) + "/contour_input.png"

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(simulationResultFile)
    reader.Update()
    output = reader.GetOutput()

    mesh = unstructuredGridToPolyData(output)
    output = mesh

    # show edges
    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.SetInputData(output)
    featureEdges.BoundaryEdgesOn()
    # featureEdges.FeatureEdgesOff()
    # featureEdges.ManifoldEdgesOff()
    # featureEdges.NonManifoldEdgesOff()
    featureEdges.Update()

    # mapper
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(output)
    mapper.ScalarVisibilityOff()

    # Visualise edges
    edgeMapper = vtk.vtkPolyDataMapper()
    edgeMapper.SetInputConnection(featureEdges.GetOutputPort())
    edgeActor = vtk.vtkActor()
    edgeActor.SetMapper(edgeMapper)

    # Create the Actor
    actor = vtk.vtkActor()
    # actor.GetProperty().EdgeVisibilityOn()
    # actor.GetProperty().SetSpecular(0.6)
    # actor.GetProperty().SetSpecularPower(30)
    # actor.GetProperty().SetLineWidth(2.0)
    actor.SetMapper(mapper)
    # backface = vtk.vtkProperty()
    # backface.SetColor(colors.GetColor3d("green"))
    # actor.SetBackfaceProperty(backface)

    # Camera setting
    # camera = vtk.vtkCamera()
    # camera.SetViewUp(-0.00168938124049409, 0.995662775211076, -0.0930203421301063)
    # camera.SetViewAngle(30)
    # camera.SetClippingRange(1.12612880757745, 1.4342182273413)
    # camera.SetPosition(0.219592083767933, -0.137774627753506, -1.21051641507013)
    # camera.SetFocalPoint(0.237973213540144, -0.0247719586331365, -0.00130243207377139)
    # camera.SetParallelScale(0.259807631459159)
    # camera.SetParallelProjection(0)

    # # camera.Set2DManipulators(1,3,2,2,2,2,3,1,4)
    # # camera.Set3DManipulators(4,1,2,3,4,1,2,4,2)

    # print(camera.GetFocalPoint(),
    # camera.GetPosition(),
    # camera.GetViewUp(),
    # camera.GetViewAngle(),
    # camera.GetClippingRange())

    # Create the Renderer
    renderer = vtk.vtkRenderer()
    # renderer.AddActor(actor)
    # renderer.AddActor(edgeActor)
    renderer.AddActor(actor)
    # renderer.SetBackground(1, 1, 1)  # Set background to white
    renderer.SetBackground(colors.GetColor3d("White"))
    # renderer.SetActiveCamera(camera)

    '''
    renderer = ren
    renderer_window = renWin
    interactor = iren
    '''

    # Create the RendererWindow
    renderer_window = vtk.vtkRenderWindow()
    renderer_window.AddRenderer(renderer)
    renderer_window.Render()

    # Create the RendererWindowInteractor and display the vtk_file
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderer_window)

    # screenshot code:
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(renderer_window)
    w2if.SetInputBufferTypeToRGB()
    w2if.ReadFrontBufferOff()
    w2if.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(mesh_2d_path)
    writer.SetInputConnection(w2if.GetOutputPort())
    writer.Write()

    ## enable user interface interactor
    # interactor.Initialize()
    # interactor.Start()

    image = cv2.imread(mesh_2d_path)
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    edged = cv2.Canny(gray, 150, 200)
    # edged = 255-edged
    contours, hierarchy = cv2.findContours(edged, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    contour_mask = np.zeros(image.shape) + 255
    cv2.drawContours(contour_mask, contours, -1, (0, 0, 0), 3)
    # cv2.imwrite(contour_path, edged)
    cv2.imwrite(contour_path, contour_mask)



if __name__ == '__main__':

    data_path = cfg.Data_path_ps
    num_samples = cfg.num_simulations_ps
    start_num = cfg.startNum_simulations_ps
    for i in range(start_num+1, num_samples+1):
        generate_contour(data_path, i)
        print(i, "/", num_samples, " done.")
