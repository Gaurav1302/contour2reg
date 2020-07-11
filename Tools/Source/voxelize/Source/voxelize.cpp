#include <vtkVersion.h>
#include <vtkProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkLine.h>
#include <vtkImageData.h>
#include <vtkProbeFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkFloatArray.h>
#include <vtkWarpScalar.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkWarpVector.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkGeometryFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkIntersectionPolyDataFilter.h>
//include <vtkFillHolesFilter.h>
#include <vtkAppendFilter.h>
#include <vtkCellCenters.h>
#include <vtkPointInterpolator.h>
#include <vtkPointLocator.h>
#include <vtkGaussianKernel.h>
#include <vtkCellDataToPointData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPolyDataNormals.h>
#include <vtkCleanPolyData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkMath.h>
#include <sstream>


void writeCSV( vtkSmartPointer<vtkUnstructuredGrid> ugrid, std::string filename );
void GetCellCenter(vtkUnstructuredGrid* poly, const unsigned int cellId, double center[3]);
int checkSideOfPlane( double testPoint[3], double closestPoint[3], double normal[3] );

bool isBorderCell( int x, int y, int z );

int main(int argc, char* argv[])
{
    //parse command line arguments
    if(argc < 2)
    {
        std::cerr << "Usage: " << argv[0]
            << " Filename(.vtu) [Resolution [Padding [OutputDir]]] [-b boundaryFile] [-s cubeSize] [-m maxDisplacement ]" << std::endl;
        return EXIT_FAILURE;
    }
    std::string filename = argv[1];

    // The number of voxels in each direction of the output mesh
    float outputResolution = 20;
    // Padding added in each direction around the output mesh:
    float padding = 0;

    if( argc > 2 )
    {
        std::stringstream sstr;
        sstr.str( argv[2] );
        sstr >> outputResolution;
    }
    if( argc > 3 )
    {
        std::stringstream sstr;
        sstr.str( argv[3] );
        sstr >> padding;
    }
    std::cout << "Resolution: " << outputResolution << ", Padding: " << padding << std::endl;

    std::string outputFolder = ".";
    if(argc > 4)
    {
        outputFolder = argv[4];
    }

    std::string boundaryFile = "";
    for( int i = 0; i < argc-1; i++ )
    {
        if( strcmp(argv[i], "-b") == 0)
        {
            if( i + 1 >= argc )
            {
                std::cerr << "Wrong usage of '-b'" << std::endl;
                return EXIT_FAILURE;
            }
            boundaryFile = argv[i+1];
            std::cout << "Will add boundary conditions from file: " << boundaryFile << std::endl;
        }
    }

    bool useFixedSize = false;
    float fixedSize = 0.3;
    for( int i = 0; i < argc-1; i++ )
    {
        if( strcmp(argv[i], "-s") == 0)
        {
            if( i + 1 >= argc )
            {
                std::cerr << "Wrong usage of '-s'" << std::endl;
                return EXIT_FAILURE;
            }
            std::stringstream sstr;
            sstr.str( argv[i+1] );
            sstr >> fixedSize;
            std::cout << "Will voxelize the cube of size " << fixedSize << " around the origin." << std::endl;
            useFixedSize = true;
        }
    }

    bool maxDisplacementCheck = false;
    double maxDisplacementThreshold = -1;
    for( int i = 0; i < argc-1; i++ )
    {
        if( strcmp(argv[i], "-m") == 0)
        {
            if( i + 1 >= argc )
            {
                std::cerr << "Wrong usage of '-m'" << std::endl;
                return EXIT_FAILURE;
            }
            std::stringstream sstr;
            sstr.str( argv[i+1] );
            sstr >> maxDisplacementThreshold;
            std::cout << "Will check for maximum displacement of " << maxDisplacementThreshold << "m" << std::endl;
            maxDisplacementCheck = true;
        }
    }

    // -----------------------------------------------------
    // Read the input file. This file should be the result of a finite element simulation and
    // should have a 3D "displacement" vector for every point. It should be an unstructured
    // grid (i.e. arbitary points).
    // This could, for example, be generated using the ElmerGUI simulation software.
    vtkUnstructuredGrid* ugrid;
    vtkSmartPointer<vtkXMLUnstructuredGridReader> readerXML =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    vtkSmartPointer<vtkUnstructuredGridReader> reader =
        vtkSmartPointer<vtkUnstructuredGridReader>::New();
    if( filename.substr(filename.find_last_of(".") + 1) == "vtu" )
    {
        std::cout << "File found. Opening:" << std::endl;
        readerXML->SetFileName(filename.c_str());
        std::cout << "Opening: " << readerXML->GetFileName() << std::endl;
        readerXML->Update();
        ugrid = readerXML->GetOutput();
        std::cout << "\tRead." << std::endl;
    } else if( filename.substr(filename.find_last_of(".") + 1) == "vtk" ) {
        std::cout << "File found. Opening:" << std::endl;
        reader->SetFileName(filename.c_str());
        std::cout << "Opening: " << reader->GetFileName() << std::endl;
        reader->Update();
        ugrid = reader->GetOutput();
        std::cout << "\tRead." << std::endl;
    } else {
        throw("Input File must be .vtu or .vtk!");
    }

    // Remove unused data:
    if( ugrid->GetPointData()->HasArray("vonmises") )
        ugrid->GetPointData()->RemoveArray("vonmises");
    if( ugrid->GetPointData()->HasArray("stress 1") )
        ugrid->GetPointData()->RemoveArray("stress 1");
    if( ugrid->GetPointData()->HasArray("stress 2") )
        ugrid->GetPointData()->RemoveArray("stress 2");
    if( ugrid->GetPointData()->HasArray("stress 3") )
        ugrid->GetPointData()->RemoveArray("stress 3");
    if( ugrid->GetPointData()->HasArray("stress 4") )
        ugrid->GetPointData()->RemoveArray("stress 4");
    if( ugrid->GetPointData()->HasArray("stress 5") )
        ugrid->GetPointData()->RemoveArray("stress 5");
    if( ugrid->GetPointData()->HasArray("stress 6") )
        ugrid->GetPointData()->RemoveArray("stress 6");


    // -----------------------------------------------------
    // If there is an array called "displacement" then we assume this mesh is the result of a
    // ElmerSolver FEM simulation.
    // Elmer applies the deformation (i.e. displacement) to the points, meaning it moves the points.
    // We need to un-apply the deformation, to get back to the reference frame, i.e. to determine
    // where the points were initially. Note that we do keep the "displacement" information in
    // the seperate array, we don't throw it away. We just move the points back to where they came
    // from. This can be done by using the vtkWarpVector algorithm.
    bool isSimulationResult = false;
    if( ugrid->GetPointData()->HasArray("displacement") )
    {
        // If there's a "displacement" array, assume this is the result of an Elmer simulation:
        isSimulationResult = true;

        vtkDataArray* warpData = ugrid->GetPointData()->GetArray("displacement");

        // Check the maximum displacement. If it is larger than threshold, refuse to voxelize:
        if( maxDisplacementCheck )
        {
            double maxDisplacement = warpData->GetMaxNorm();
            std::cout << "Found maximum displacement: " << maxDisplacement << std::endl;
            if( maxDisplacement > maxDisplacementThreshold )
            {
                std::cerr << "Maximum displacement in this sample is greater than the allowed threshold "
                          << maxDisplacementThreshold << std::endl;
                std::cerr << "Aborting voxelization." << std::endl;
                return EXIT_FAILURE;
            }
        }

        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(ugrid->GetPoints());
        polydata->GetPointData()->AddArray(warpData);
        polydata->GetPointData()->SetActiveVectors(warpData->GetName());

        // WarpVector will use the array marked as active vector in polliver.vtuydata.
        // It has to be a 3 component array with the same number of tuples as points in polydata.
        vtkSmartPointer<vtkWarpVector> warpVector =
          vtkSmartPointer<vtkWarpVector>::New();
        warpVector->SetScaleFactor(-1);
        warpVector->SetInputData(polydata);
        warpVector->Update();

        // Update the points in the unstructured grid to move them back to the undeformed positions:
        ugrid->SetPoints( warpVector->GetPolyDataOutput()->GetPoints() );

        // Optionally write the unwarped mesh, for debugging purposes:
        /*vtkSmartPointer<vtkXMLUnstructuredGridWriter> ugridWriter =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        ugridWriter->SetFileName("unwarped.vtu");
        ugridWriter->SetInputData( ugrid );
        ugridWriter->Write();*/
    } else if( ugrid->GetPointData()->HasArray("rest_position") )   // SOFA simulation? Also remove distortion
    {
        // If there's a "displacement" array, assume this is the result of an Elmer simulation:
        isSimulationResult = true;

        vtkDataArray* restPositions = ugrid->GetPointData()->GetArray("rest_position");
        vtkSmartPointer<vtkDoubleArray> warpData = vtkSmartPointer<vtkDoubleArray>::New();
        warpData->SetNumberOfComponents(3);
        warpData->SetNumberOfTuples(ugrid->GetNumberOfPoints());
        warpData->SetName("displacement");
        for( int i = 0; i < ugrid->GetNumberOfPoints(); i++ )
        {
          double* displaced;
          double* original;
          displaced = ugrid->GetPoint(i);
          original = restPositions->GetTuple3(i);
          double dx = displaced[0] - original[0];
          double dy = displaced[1] - original[1];
          double dz = displaced[2] - original[2];
          warpData->SetTuple3(i, dx, dy, dz);
          /*std::cout << original[0] << " ";
          std::cout << original[1] << " ";
          std::cout << original[2] << " -> ";
          std::cout << displaced[0] << " ";
          std::cout << displaced[1] << " ";
          std::cout << displaced[2] << " Diff:";
          std::cout << dx << " ";
          std::cout << dy << " ";
          std::cout << dz << std::endl;*/

          ugrid->GetPoints()->SetPoint(i, original);
        }

        // Check the maximum displacement. If it is larger than threshold, refuse to voxelize:
        if( maxDisplacementCheck )
        {
            double maxDisplacement = warpData->GetMaxNorm();
            std::cout << "Found maximum displacement: " << maxDisplacement << std::endl;
            if( maxDisplacement > maxDisplacementThreshold )
            {
                std::cerr << "Maximum displacement in this sample is greater than the allowed threshold "
                          << maxDisplacementThreshold << std::endl;
                std::cerr << "Aborting voxelization." << std::endl;
                return EXIT_FAILURE;
            }
        }

        /*vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(ugrid->GetPoints());
        polydata->GetPointData()->AddArray(warpData);
        polydata->GetPointData()->SetActiveVectors(warpData->GetName());

        // WarpVector will use the array marked as active vector in polliver.vtuydata.
        // It has to be a 3 component array with the same number of tuples as points in polydata.
        vtkSmartPointer<vtkWarpVector> warpVector =
          vtkSmartPointer<vtkWarpVector>::New();
        warpVector->SetScaleFactor(-1);
        warpVector->SetInputData(polydata);
        warpVector->Update();

        // Update the points in the unstructured grid to move them back to the undeformed positions:
        ugrid->SetPoints( warpVector->GetPolyDataOutput()->GetPoints() );*/

        ugrid->GetPointData()->AddArray(warpData);

        // Optionally write the unwarped mesh, for debugging purposes:
        /*std::string unwarpedFileName = outputFolder + "/unwarped.vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> ugridWriter =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        ugridWriter->SetFileName(unwarpedFileName.c_str());
        ugridWriter->SetInputData( ugrid );
        ugridWriter->Write();*/
    }

    for( int i = 0; i < ugrid->GetCellData()->GetNumberOfArrays(); i++ )
    {
        std::cout << "Cell array: " << ugrid->GetCellData()->GetArrayName(i) << std::endl;
    }
    for( int i = 0; i < ugrid->GetPointData()->GetNumberOfArrays(); i++ )
    {
        std::cout << "Point array: " << ugrid->GetPointData()->GetArrayName(i) << std::endl;
    }


    // -----------------------------------------------------
    // Add an array representing the material:
    /*vtkSmartPointer<vtkFloatArray> material = vtkSmartPointer<vtkFloatArray>::New();
    material->SetNumberOfComponents(1);
    material->SetName("material");
    for( unsigned int i = 0; i < ugrid->GetNumberOfPoints(); i++ )
        material->InsertNextTuple1(1);
    ugrid->GetPointData()->AddArray( material );*/


    // -----------------------------------------------------
    // Next, we want to interpolate the data onto a structured, rectiliniar grid.
    // First, we create the grid, then we use the vtkProbeFilter to resample the
    // unstructured grid onto the structured grid.

    // Create a grid of points to interpolate over
    vtkSmartPointer<vtkPoints> gridPoints =
            vtkSmartPointer<vtkPoints>::New();
    double bounds[6];
    double center[3];
    ugrid->GetPoints()->GetBounds( bounds );
    center[0] = (bounds[0] + bounds[1])*0.5;
    center[1] = (bounds[2] + bounds[3])*0.5;
    center[2] = (bounds[4] + bounds[5])*0.5;
    std::cout << "Input Center:" << std::endl;
    std::cout << "\t " << center[0] << " " << center[1] << " " << center[2] << std::endl;
    std::cout << "Input Bounds:" << std::endl;
    std::cout << "\tx: " << bounds[0] << " - " << bounds[1] << std::endl;
    std::cout << "\ty: " << bounds[2] << " - " << bounds[3] << std::endl;
    std::cout << "\tz: " << bounds[4] << " - " << bounds[5] << std::endl;

    // Make cube large enough to contain full model, determine the longest side of the mesh:
    double outputSize = bounds[1] - bounds[0];
    if( useFixedSize )
    {
        outputSize = fixedSize;
        center[0] = 0;
        center[1] = 0;
        center[2] = 0;
    } else {
        if( bounds[3] - bounds[2] > outputSize )
            outputSize = bounds[3] - bounds[2];
        if( bounds[5] - bounds[4] > outputSize )
            outputSize = bounds[5] - bounds[4];
        outputSize = outputSize + 2*padding;
    }
    double stepSize = (outputSize)/(outputResolution);
    double minX = center[0] - 0.5*outputSize;
    double maxX = center[0] + 0.5*outputSize;
    double minY = center[1] - 0.5*outputSize;
    double maxY = center[1] + 0.5*outputSize;
    double minZ = center[2] - 0.5*outputSize;
    double maxZ = center[2] + 0.5*outputSize;
    double epsilon = 1e-10;
    for ( double x = minX; x < maxX + epsilon; x += stepSize )
    {
        for ( double y = minY; y < maxY + epsilon; y += stepSize )
        {
            for ( double z = minZ; z < maxZ + epsilon; z += stepSize )
            {
                gridPoints->InsertNextPoint ( x, y, z );
            }
        }
    }

    // Create a dataset from the grid points
    vtkSmartPointer<vtkPolyData> gridPolyData = vtkSmartPointer<vtkPolyData>::New();
    gridPolyData->SetPoints(gridPoints);


    // -----------------------------------------------------
    // Generate a structured Grid.
    vtkSmartPointer<vtkStructuredGrid> sgrid = vtkSmartPointer<vtkStructuredGrid>::New();
    sgrid->SetDimensions( outputResolution+1, outputResolution+1, outputResolution+1 );
    // Apply the points to the structured grid. Now we have the cells/surfaces etc. setup.
    sgrid->SetPoints(gridPoints);

    vtkSmartPointer<vtkCellCenters> centers = vtkSmartPointer<vtkCellCenters>::New();
    centers->SetInputData( sgrid );
    centers->Update();
    vtkPolyData* cellCenters = centers->GetOutput();

    // Perform the interpolation
    vtkSmartPointer<vtkPointInterpolator> interpolator = vtkSmartPointer<vtkPointInterpolator>::New();
    vtkSmartPointer<vtkGaussianKernel> gaussian = vtkSmartPointer<vtkGaussianKernel>::New();
    gaussian->SetRadius(stepSize*15);
    gaussian->SetSharpness(10);
    interpolator->SetSourceData( ugrid );
    interpolator->SetInputData( cellCenters );
    interpolator->SetKernel( gaussian );
    interpolator->SetNullPointsStrategy(vtkPointInterpolator::CLOSEST_POINT);
    interpolator->Update();
    //vtkDataSet* interpolated = interpolator->GetOutput();

    /*vtkSmartPointer<vtkProbeFilter> probeFilter = vtkSmartPointer<vtkProbeFilter>::New();
    probeFilter->SetValidPointMaskArrayName("mesh");
    probeFilter->SetSourceData( ugrid );
    probeFilter->SetInputData( cellCenters );      // Interpolate 'Source' at these points
    probeFilter->Update();*/

    // GMSH and Elmer work with unstructured grids, so even though our data is inside a structured grid,
    // we need to write it to an unstructured grid:
    vtkSmartPointer<vtkUnstructuredGrid> voxelGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    voxelGrid->SetPoints( sgrid->GetPoints() );
    vtkSmartPointer<vtkCellArray> arr = vtkSmartPointer<vtkCellArray>::New();
    for( unsigned int i = 0; i < sgrid->GetNumberOfCells(); i++ )
    {
        arr->InsertNextCell(sgrid->GetCell(i));
        //std::cout << sgrid->GetCell(i) << std::endl;
    }
    voxelGrid->SetCells(12, arr);
    voxelGrid->Modified();
    voxelGrid->GetCellData()->AddArray(cellCenters->GetCellData()->GetArray("centers"));

    //std::cout << sgrid->GetNumberOfCells() << " " << voxelGrid->GetNumberOfCells() << " cells" << std::endl;
    //std::cout << sgrid->GetNumberOfPoints() << " " << voxelGrid->GetNumberOfPoints() << " points" << std::endl;

    // Copy all the resampled arrays to the voxel grid:
    unsigned int numberOfArrays = interpolator->GetOutput()->GetPointData()->GetNumberOfArrays();
    std::cout << "Copying arrays:" << std::endl;
    for(unsigned int i = 0; i < numberOfArrays; i++)
    {
        std::cout << "\t" << interpolator->GetOutput()->GetPointData()->GetArray(i)->GetName() << std::endl;
        //std::cout << "\t\t" << interpolator->GetOutput()->GetPointData()->GetArray(i)->GetNumberOfTuples() << std::endl;
        voxelGrid->GetCellData()->AddArray(interpolator->GetOutput()->GetPointData()->GetArray(i));
        //sgrid->GetCellData()->AddArray(interpolator->GetOutput()->GetPointData()->GetArray(i));
    }
    /*unsigned int numberOfArrays = probeFilter->GetOutput()->GetPointData()->GetNumberOfArrays();
    std::cout << "Copying arrays:" << std::endl;
    for(unsigned int i = 0; i < numberOfArrays; i++)
    {
        std::cout << "\t" << probeFilter->GetOutput()->GetPointData()->GetArray(i)->GetName() << std::endl;
        voxelGrid->GetCellData()->AddArray(probeFilter->GetOutput()->GetPointData()->GetArray(i));
        sgrid->GetCellData()->AddArray(probeFilter->GetOutput()->GetPointData()->GetArray(i));
    }*/


    // -----------------------------------------------------
    // If there's a boundary file given, also take the data from that file:
    if( boundaryFile.size() != 0 )
    {
        std::cout << "Calculating which voxels are on the surface of the mesh:" << std::endl;

        vtkSmartPointer<vtkXMLUnstructuredGridReader> boundaryConditionsReader =
                    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
        boundaryConditionsReader->SetFileName(boundaryFile.c_str());
        boundaryConditionsReader->Update();
        vtkUnstructuredGrid* boundaryCond = boundaryConditionsReader->GetOutput();

        // Extract the surface:
        vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter =
            vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        surfaceFilter->SetInputData(boundaryCond);
        surfaceFilter->Update();
        vtkPolyData* surfacePolyData = surfaceFilter->GetOutput();
        std::cout << "Number of extracted surface elements: " << surfacePolyData->GetNumberOfCells() << std::endl;


        // Write Cleaned (DEBUG ONLY):
        std::string oFilename = outputFolder + "cleaned.vtp";

        vtkSmartPointer<vtkCellDataToPointData> c2p = vtkSmartPointer<vtkCellDataToPointData>::New();
        c2p->SetInputData( boundaryCond );
        c2p->Update();

        vtkSmartPointer<vtkPointInterpolator> surfInterpolator = vtkSmartPointer<vtkPointInterpolator>::New();
        vtkSmartPointer<vtkGaussianKernel> surfGaussian = vtkSmartPointer<vtkGaussianKernel>::New();
        surfGaussian->SetRadius(stepSize*2);
        surfGaussian->SetSharpness(10);
        surfInterpolator->SetSourceData( c2p->GetOutput() );
        surfInterpolator->SetInputData( cellCenters );
        surfInterpolator->SetKernel( surfGaussian );
        surfInterpolator->SetNullPointsStrategy(vtkPointInterpolator::MASK_POINTS);
        surfInterpolator->Update();

        // Copy all the resampled arrays to the voxel grid:
        unsigned int numberOfArrays = surfInterpolator->GetOutput()->GetPointData()->GetNumberOfArrays();
        std::cout << "Copying arrays: " << numberOfArrays << std::endl;
        for(unsigned int i = 0; i < numberOfArrays; i++)
        {
            std::cout << "\t" << surfInterpolator->GetOutput()->GetPointData()->GetArray(i)->GetName() << std::endl;
            voxelGrid->GetCellData()->AddArray(surfInterpolator->GetOutput()->GetPointData()->GetArray(i));
            //sgrid->GetCellData()->AddArray(surfInterpolator->GetOutput()->GetPointData()->GetArray(i));
        }
        numberOfArrays = surfInterpolator->GetOutput()->GetCellData()->GetNumberOfArrays();
        std::cout << "Copying cell arrays: " << numberOfArrays << std::endl;
        for(unsigned int i = 0; i < numberOfArrays; i++)
        {
            std::cout << "\t" << surfInterpolator->GetOutput()->GetCellData()->GetArray(i)->GetName() << std::endl;
            voxelGrid->GetCellData()->AddArray(surfInterpolator->GetOutput()->GetCellData()->GetArray(i));
            //sgrid->GetCellData()->AddArray(interpolator->GetOutput()->GetCellData()->GetArray(i));
        }

        vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
        normalGenerator->SetInputData( surfacePolyData );
        normalGenerator->ComputePointNormalsOn();
        normalGenerator->ComputeCellNormalsOff();
        normalGenerator->AutoOrientNormalsOn();
        normalGenerator->SplittingOff();
        normalGenerator->Update();
        //std::cout << "\tPoints:" << normalGenerator->GetOutput()->GetNumberOfPoints() << std::endl;
        //std::cout << "\tCells:" << normalGenerator->GetOutput()->GetNumberOfCells() << std::endl;

        // Get the generated Normals:
        vtkDataArray* pointNormals = normalGenerator->GetOutput()->GetPointData()->GetArray("Normals");
        //for( unsigned int i = 0; i < normalGenerator->GetOutput()->GetPointData()->GetNumberOfArrays(); i ++ )
        //	std::cout << "array "  << normalGenerator->GetOutput()->GetPointData()->GetArray(i)->GetName() << std::endl;

        vtkSmartPointer<vtkDoubleArray> sdf = vtkSmartPointer<vtkDoubleArray>::New();
        sdf->SetNumberOfComponents(1);
        sdf->SetNumberOfTuples( voxelGrid->GetNumberOfCells() );
        sdf->SetName("signedDistance");

        // Create the tree
        vtkSmartPointer<vtkPointLocator> pointLocator =
            vtkSmartPointer<vtkPointLocator>::New();
        pointLocator->SetDataSet( surfacePolyData );
        pointLocator->BuildLocator();

        std::cout << "Calculating SDF" << std::endl;

        for( unsigned int i = 0; i < cellCenters->GetNumberOfPoints(); i++ )
        {
            double testPoint[3];
            cellCenters->GetPoint( i, testPoint );

            //Find the closest points to TestPoint
            vtkIdType pointID = pointLocator->FindClosestPoint(testPoint);
            double closestPoint[3];
            surfacePolyData->GetPoint( pointID, closestPoint );
            double dist = sqrt(vtkMath::Distance2BetweenPoints( testPoint, closestPoint ));
            double normal[3];
            pointNormals->GetTuple( pointID, normal );
            if( checkSideOfPlane( testPoint, closestPoint, normal ) == -1 )
            {
                sdf->SetTuple1( i, -dist );
            } else {
                sdf->SetTuple1( i, dist );
            }
        }
        //std::cout << "Calculated SDF" << std::endl;

        // Calculate signed distance function:
        vtkSmartPointer<vtkImplicitPolyDataDistance> distFunction = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
        distFunction->SetInput( surfacePolyData );
        vtkSmartPointer<vtkDoubleArray> signedDistance = vtkSmartPointer<vtkDoubleArray>::New();
        //signedDistance->SetNumberOfComponents(1);
        //signedDistance->SetNumberOfTuples(cellCenters->GetNumberOfPoints());
        signedDistance->SetName("signedDistance");
        distFunction->FunctionValue( cellCenters->GetPoints()->GetData(), signedDistance );

        //voxelGrid->GetCellData()->AddArray(signedDistance);
        voxelGrid->GetCellData()->AddArray(sdf);

    }


    // Remove unused data:

    if( voxelGrid->GetCellData()->HasArray("Normals") )
        voxelGrid->GetCellData()->RemoveArray("Normals");
    if( voxelGrid->GetCellData()->HasArray("vtkValidPointMask") )
        voxelGrid->GetCellData()->RemoveArray("vtkValidPointMask");


    // -----------------------------------------------------
    // Output:
    std::string outBaseFileName = outputFolder + "/voxels_" + (isSimulationResult ? "result" : "input");

    // Write file
    std::string oFilename = outBaseFileName + ".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(oFilename.c_str());
    writer->SetInputData(voxelGrid);
    writer->Write();
    std::cout << "Wrote " << oFilename << std::endl;

    // Write again as csv file:
    //oFilename = outBaseFileName + ".csv";
    //writeCSV( voxelGrid, oFilename );

    return EXIT_SUCCESS;
}

void writeCSV( vtkSmartPointer<vtkUnstructuredGrid> ugrid, std::string filename )
{
    vtkCellData* cData = ugrid->GetCellData();

    vtkSmartPointer<vtkCellCenters> centers = vtkSmartPointer<vtkCellCenters>::New();
    centers->SetInputData( ugrid );
    centers->Update();
    vtkPolyData* cellCenters = centers->GetOutput();

    std::ofstream oFile(filename.c_str());
    int numberOfCellArrays = cData->GetNumberOfArrays();
    // Header:
    for(int i = 0; i < numberOfCellArrays; i++)
    {
        if( cData->GetArray(i)->GetNumberOfComponents() > 1 ) {
            for( int k = 0; k < cData->GetArray(i)->GetNumberOfComponents(); ++k )
            {
                oFile << cData->GetArrayName(i) << "(" << k << ") ";
            }
        } else {
            oFile << cData->GetArrayName(i) << " ";
        }
    }
    oFile << "pos(0) pos(1) pos(2)" << std::endl;

    // Content:
    double pos[3];
    for( int j = 0; j < ugrid->GetNumberOfCells(); ++j )
    {
        //ugrid->GetPoint(j, pos);
        cellCenters->GetPoint(j, pos );
        for(int i = 0; i < numberOfCellArrays; i++)
        {
            //double* vals = pData->GetArray(i)->GetTuple(j);
            //vtkTuple tuple = pData->GetArray(i)->GetT
            for( int k = 0; k < cData->GetArray(i)->GetNumberOfComponents(); ++k )
                oFile << cData->GetArray(i)->GetComponent(j, k) << " ";
        }
        oFile << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    }

    oFile.close();
    std::cout << "Wrote " << filename << std::endl;
}

void GetCellCenter(vtkUnstructuredGrid* poly, const unsigned int cellId, double center[3])
{
  double pcoords[3] = {0,0,0};
  double weights[poly->GetMaxCellSize()];

  vtkCell* cell = poly->GetCell(cellId);
  int subId = cell->GetParametricCenter(pcoords);
  cell->EvaluateLocation(subId, pcoords, center, weights);
}

int checkSideOfPlane( double testPoint[3], double closestPoint[3], double normal[3] )
{
	double pThisSide[3];
	double pOtherSide[3];
	pThisSide[0] = closestPoint[0] + normal[0];
	pThisSide[1] = closestPoint[1] + normal[1];
	pThisSide[2] = closestPoint[2] + normal[2];
	pOtherSide[0] = closestPoint[0] - normal[0];
	pOtherSide[1] = closestPoint[1] - normal[1];
	pOtherSide[2] = closestPoint[2] - normal[2];
	double dist2ThisSide = vtkMath::Distance2BetweenPoints( testPoint, pThisSide );
	double dist2OtherSide = vtkMath::Distance2BetweenPoints( testPoint, pOtherSide );
	if( dist2ThisSide <= dist2OtherSide )
		return 1;
	else
		return -1;

}



