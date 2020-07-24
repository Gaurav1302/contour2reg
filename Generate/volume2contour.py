import numpy as np
import vtk
from vtk import *
from vtk.util.numpy_support import vtk_to_numpy
from vtk import vtkStructuredPointsReader
import matplotlib
from util import unstructuredGridToPolyData
# matplotlib.use('TkAgg')
matplotlib.use('agg')
import matplotlib.pyplot as plt


def main():

    colors = vtk.vtkNamedColors()

    # ren = vtk.vtkRenderer()
    # renWin = vtk.vtkRenderWindow()
    # renWin.AddRenderer(ren)

    # # create a renderwindowinteractor
    # iren = vtk.vtkRenderWindowInteractor()
    # iren.SetRenderWindow(renWin)




    # create source/read.vtu file
    simulationResultFile = "../Data/1/simulation/result/case_t0001.vtu"

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

    save_filename = "../../TestScreenshot.png"
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(save_filename)
    writer.SetInputConnection(w2if.GetOutputPort())
    writer.Write()


    # enable user interface interactor
    # interactor.Initialize()
    # interactor.Start()

    import cv2
    import numpy as np

    image = cv2.imread(save_filename)

    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    edged = cv2.Canny(gray, 150, 200)
    # cv2.waitKey(0)

    contours, hierarchy = cv2.findContours(edged, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)


    # plt.imshow(edged)
    # plt.show()
    # cv2.imshow('Canny Edges After Contouring', edged)
    # cv2.waitKey(0)

    cv2.imwrite('../../edge.png', edged)

    cv2.drawContours(image, contours, -1, (0, 255, 0), 3)

    # cv2.imshow('Contours', image)
    cv2.imwrite('../../cont.png', edged)
    # cv2.waitKey(0)
    # cv2.destroyAllWindows()


if __name__ == '__main__':
    main()
