#include "writeElmerFiles.h"

#include <vtkSurfaceReconstructionFilter.h>
#include <vtkPolyData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkStructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkPointLocator.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPointSet.h>
#include <vtkPointDataToCellData.h>
#include <vtkArrayCalculator.h>
#include <vtkDoubleArray.h>
#include <fstream>
#include <map>

#define BOUNDARY_ZERO_FORCE 1
#define BOUNDARY_ZERO_DISPLACEMENT 2
#define BOUNDARY_FIXED_DISPLACEMENT 3

vtkIdType parentForCell( vtkCell* cell, vtkUnstructuredGrid* ugrid );

double manipulateX = 0;
double manipulateY = 0;
double manipulateZ = 0;

void writeElmerFiles( vtkSmartPointer<vtkUnstructuredGrid> ugrid, std::string path )
{
    std::map<int,int> numberOfElements;

    // Compute the surface of the mesh:
    /*vtkSmartPointer<vtkDataSetSurfaceFilter> surf =
            vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surf->SetInputData(ugrid);
    surf->Update();*/

    // Compute the surface of the mesh:
    /*vtkSmartPointer<vtkStructuredGridGeometryFilter> surf =
            vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
    surf->SetExtent(10, 20, 0, 10, 0, 0);
    surf->SetInputData(sgrid);
    surf->Update();*/
    //vtkPolyData* surface = surf->GetOutput();

    vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
    locator->SetDataSet(ugrid);

    // ---------------------------------------------------------------------
    // Write the nodes file:
    std::ofstream nodeFile( (path + "/mesh.nodes").c_str() );
    for( unsigned int i = 0; i < ugrid->GetNumberOfPoints(); i++ )
    {
        double pos[3];// = ugrid->GetPoint(i);
        ugrid->GetPoint(i, pos);
        nodeFile << i+1 << " -1 " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    }
    std::cout << "Wrote mesh.nodes" << std::endl;

    // ---------------------------------------------------------------------
    // Write the elements file:
    std::ofstream elemFile( (path + "/mesh.elements").c_str() );
    // Helper for interpolation:
    /*vtkSmartPointer<vtkPointDataToCellData> p2c =
            vtkSmartPointer<vtkPointDataToCellData>::New();
    vtkSmartPointer<vtkStructuredGrid> pointSet = vtkSmartPointer<vtkStructuredGrid>::New();
    pointSet->SetPoints( ugrid->GetPoints() );
    pointSet->SetFieldData( ugrid->GetPointData() );
    p2c->SetInputData( ugrid );
    p2c->Update();
    vtkDataSet* res = p2c->GetOutput();
    vtkDataArray* material = res->GetCellData()->GetArray("material");*/

    vtkSmartPointer<vtkIdTypeArray> volumeElementIDs = vtkSmartPointer<vtkIdTypeArray>::New();
    volumeElementIDs->SetNumberOfComponents(1);
    volumeElementIDs->SetNumberOfTuples(ugrid->GetNumberOfCells());

    int body = 1;
    int numberOfVolumeElements = 0;
    for( unsigned int i = 0; i < ugrid->GetNumberOfCells(); i++ )
    {
        vtkCell* cell = ugrid->GetCell(i);
        /*if( material->GetTuple1(i) > 0.1 )
            body = 2;
        else
            body = 1;*/

        int cellType = typeForCell( cell );
        if( numberOfElements.find(cellType) == numberOfElements.end() )
            numberOfElements[cellType] = 0;
        numberOfElements[cellType] = numberOfElements[cellType] + 1;

        if( cell->GetCellType() == VTK_TETRA )
        {
            numberOfVolumeElements ++;
            volumeElementIDs->SetTuple1( i, numberOfVolumeElements );
            elemFile << numberOfVolumeElements << " " << body << " " << cellType;
            for( unsigned int j = 0; j < cell->GetNumberOfPoints(); j++)
            {
                elemFile << " " << cell->GetPointId(j) + 1;
            }
            elemFile << "\n";
        } else {
            volumeElementIDs->SetTuple1( i, 0 );
        }
    }
    std::cout << "Wrote mesh.elements" << std::endl;


    // ---------------------------------------------------------------------
    // Write the boundary file:
    std::ofstream boundaryFile( (path + "/mesh.boundary").c_str() );
    vtkDataArray* zeroDisplacement = ugrid->GetCellData()->GetArray("zeroDisplacement");
    vtkDataArray* manipulation = ugrid->GetCellData()->GetArray("manipulation");
    /*std::cout << force->GetNumberOfTuples() << std::endl;
    vtkSmartPointer<vtkArrayCalculator> normCalculator =
            vtkSmartPointer<vtkArrayCalculator>::New();
    //normCalculator->SetAttributeModeToUseCellData();
    normCalculator->SetAttributeType( VTK_ATTRIBUTE_MODE_USE_CELL_DATA );
    normCalculator->SetInputData(surface);
    normCalculator->AddVectorArrayName("force");
    normCalculator->SetFunction("mag(force)");
    normCalculator->SetResultArrayName("forceNorm");
    normCalculator->Update();

    std::cout << "Cells: " << normCalculator->GetDataSetOutput()->GetNumberOfCells() << std::endl;
    std::cout << "Points: " << normCalculator->GetDataSetOutput()->GetNumberOfPoints() << std::endl;
    //std::cout << normCalculator->GetOutput()

    for( unsigned int i = 0; i < normCalculator->GetDataSetOutput()->GetCellData()->GetNumberOfArrays(); i++ )
    {
        std::cout <<"array: "<<normCalculator->GetDataSetOutput()->GetCellData()->GetArrayName(i)<< std::endl;
    }
    for( unsigned int i = 0; i < normCalculator->GetDataSetOutput()->GetPointData()->GetNumberOfArrays(); i++ )
    {
        std::cout <<"array2: "<<normCalculator->GetDataSetOutput()->GetPointData()->GetArrayName(i)<< std::endl;
    }*/

    //vtkDataArray* forceNorm = normCalculator->GetDataSetOutput()->GetCellData()->GetArray("forceNorm");

    vtkSmartPointer<vtkDoubleArray> manipulationNorm = vtkSmartPointer<vtkDoubleArray>::New();
    manipulationNorm->SetNumberOfComponents(1);
    manipulationNorm->SetNumberOfTuples(manipulation->GetNumberOfTuples());
    for( unsigned int i = 0; i < manipulation->GetNumberOfTuples(); i++ )
    {
        double f[3];
        manipulation->GetTuple(i,f);
        double norm = sqrt( f[0]*f[0] + f[1]*f[1] + f[2]*f[2] );
        manipulationNorm->SetTuple1(i,norm);
    }
    int numberOfSurfaceElements = 0;
    for( unsigned int i = 0; i < ugrid->GetNumberOfCells(); i++ )
    {
        vtkCell* cell = ugrid->GetCell(i);
        int cellType = typeForCell( cell );
        if( cell->GetCellType() == VTK_TRIANGLE )
        {
            vtkIdType parent = parentForCell( cell, ugrid );
            vtkIdType parentAsVolumeIndex = volumeElementIDs->GetTuple1( parent );       // As volume nodes have new indexing, use this.
            int boundaryIndex = BOUNDARY_ZERO_FORCE;
            if( zeroDisplacement->GetTuple1(i) > 0 )
            {
                boundaryIndex = BOUNDARY_ZERO_DISPLACEMENT;

            } else if( manipulationNorm->GetTuple1(i) > 0 ) {
                double* f;
                f = manipulation->GetTuple3(i);
                manipulateX = f[0];
                manipulateY = f[1];
                manipulateZ = f[2];
                boundaryIndex = BOUNDARY_FIXED_DISPLACEMENT;
            }
            boundaryFile << i+1 << " " << boundaryIndex << " " << parentAsVolumeIndex << " 0 " << cellType;
            cell = ugrid->GetCell(i);
            for( unsigned int j = 0; j < cell->GetNumberOfPoints(); j++)
            {
                //double pos[3];
                //cell->GetPoints()->GetPoint(2-j, pos);
                //boundaryFile << " " << locator->FindClosestPoint( pos ) + 1;
                boundaryFile << " " << cell->GetPointId(j) + 1;
            }
            boundaryFile << "\n";
            numberOfSurfaceElements ++;
        }
    }

    std::cout << "Wrote mesh.boundary" << std::endl;

    // ---------------------------------------------------------------------
    // Write the header file: vtkSmartPointer<vtkUnstructuredGrid> ugrid, vtkSmartPointer<vtkUnstructuredGrid> surface,
    std::ofstream headerFile( (path + "/mesh.header").c_str() );
    headerFile << ugrid->GetNumberOfPoints() << " ";
    headerFile << numberOfVolumeElements << " ";
    headerFile << numberOfSurfaceElements << "\n";
    headerFile << numberOfElements.size() << "\n";
    for( std::map<int,int>::iterator it = numberOfElements.begin(); it != numberOfElements.end(); ++it )
        headerFile << it->first << " " << it->second << "\n";
    //headerFile << "808 " << ugrid->GetNumberOfCells() << "\n";
    //headerFile << "404 " << surface->GetNumberOfCells() << "\n";
    std::cout << "Wrote mesh.header" << std::endl;
}

vtkIdType parentForCell( vtkCell* cell, vtkUnstructuredGrid* ugrid )
{
    vtkIdType pID0 = cell->GetPointId(0);
    vtkIdType pID1 = cell->GetPointId(1);
    vtkIdType pID2 = cell->GetPointId(2);
    for( int i = 0; i < ugrid->GetNumberOfCells(); i++ )
    {
        vtkCell* c = ugrid->GetCell(i);
        if( c->GetCellType() != VTK_TRIANGLE )
        {
            int matchingPoints = 0;
            for( int j = 0; j < c->GetNumberOfPoints(); j++)
            {
                vtkIdType pID = c->GetPointId(j);
                if( pID0 == pID || pID1 == pID || pID2 == pID )
                    matchingPoints ++;
            }
            //if( matchingPoints > 0 )
                //std::cout << matchingPoints << " points match " << std::endl;
            if( matchingPoints == 3 )
            {
                return i;
            }
        }
    }
    return 0;
}

int vtkCellType2ElmerCellType( int cellType )
{
    if( cellType == VTK_TRIANGLE  )
        return 303;
    if( cellType == VTK_QUAD )
        return 404;
    if( cellType == VTK_HEXAHEDRON )
        return 808;
    if( cellType == VTK_TETRA )
        return 504;
    return 0;
}

int typeForCell( vtkCell* cell )
{
    if( cell->GetCellType() == VTK_TRIANGLE  )
        return 303;
    if( cell->GetCellType() == VTK_QUAD )
        return 404;
    if( cell->GetCellType() == VTK_HEXAHEDRON )
        return 808;
    if( cell->GetCellType() == VTK_TETRA )
        return 504;
    return 0;
}

void writeSimulation( std::string path, float youngsModulus )
{

    std::ofstream simfile( (path + "/case.sif").c_str() );
    // Write header:
    simfile << "Header\n" <<
            "\tCHECK KEYWORDS Warn\n" <<
            "\tMesh DB \".\" \".\"\n" <<
            "\tInclude Path \"\"\n" <<
            "\tResults Directory \"result\"\n" <<
            "End\n\n";

    // Simulation:
    simfile <<
            "Simulation\n" <<
            "\tMax Output Level = 5\n" <<
            "\tCoordinate System = Cartesian\n" <<
            "\tCoordinate Mapping(3) = 1 2 3\n" <<
            "\tSimulation Type = Steady state\n" <<
            "\tSteady State Max Iterations = 1\n" <<
            "\tSolver Input File = case.sif\n" <<
            "\tPost File = case.vtu\n" <<
            "End\n\n";

    // Constants:
    /*simfile <<
            "Constants\n" <<

            "\tStefan Boltzmann = 5.67e-08\n" <<
            "\tPermittivity of Vacuum = 8.8542e-12\n" <<
            "\tBoltzmann Constant = 1.3807e-23\n" <<
            "\tUnit Charge = 1.602e-19\n" <<
            "End\n\n";*/

            //"\tGravity(4) = 0 -1 0 9.82\n" <<

    // Solver:
    simfile <<
            "Solver 1\n" <<
            "\tEquation = Non Linear elasticity\n" <<
            "\tProcedure = \"ElasticSolve\" \"ElasticSolver\"\n" <<
            "\tVariable = -dofs 3 Displacement\n" <<
            "\tExec Solver = Always\n" <<
            "\tStabilize = True\n" <<
            "\tBubbles = False\n" <<
            "\tLumped Mass Matrix = False\n" <<
            "\tOptimize Bandwidth = True\n" <<
            "\tSteady State Convergence Tolerance = 1.0e-5\n" <<
            "\tNonlinear System Convergence Tolerance = 1.0e-7\n" <<
            "\tNonlinear System Max Iterations = 20\n" <<
            "\tNonlinear System Newton After Iterations = 3\n" <<
            "\tNonlinear System Newton After Tolerance = 1.0e-3\n" <<
            "\tNonlinear System Relaxation Factor = 1\n" <<
            "\tLinear System Solver = Iterative\n" <<
            "\tLinear System Iterative Method = GCR\n" <<
            "\tLinear System Max Iterations = 500\n" <<
            "\tLinear System Convergence Tolerance = 1.0e-10\n" <<
            "\tBiCGstabl polynomial degree = 2\n" <<
            "\tLinear System Preconditioning = ILU1\n" <<
            "\tLinear System ILUT Tolerance = 1.0e-3\n" <<
            "\tLinear System Abort Not Converged = False\n" <<
            "\tLinear System Residual Output = 1\n" <<
            "\tLinear System Precondition Recompute = 1\n" <<
            "End\n\n";

    /*imfile <<
            "Solver 2\n"
            "\tExec Solver = after timestep\n"
            "\tEquation = \"result output\"\n" <<
            "\tProcedure = \"ResultOutputSolve\" \"ResultOutputSolver\"\n" <<
            "\tOutput File Name = \"result/case\"\n" <<
            "\tOutput Format = String \"vtu\"\n" <<:w
            "\tBinary Output = True\n" <<
            "End\n\n";*/

    // Equation:
    simfile <<
            "Equation 1\n" <<
            "\tName = \"Elasticity\"\n" <<
            //"\tCalculate Stresses = True\n" <<
            "\tActive Solvers(1) = 1\n" <<
            "End\n\n" ;

    simfile <<
            "Body 1\n" <<
            "\tTarget Bodies(1) = 1\n" <<
            "\tEquation = 1\n" <<
            "\tMaterial = 1\n" <<
            "End\n\n" ;
    simfile <<
            "Material 1\n" <<
            "\tName = \"Liver (Human, approx)\"\n" <<
            "\tDensity = 1060.0\n" <<
            "\tPoisson ratio = 0.35\n" <<
            "\tYoungs modulus = " << youngsModulus << "\n" <<
            "End\n\n" ;

    /*simfile <<
            "Boundary Condition 1\n" <<
            "\tTarget Boundaries(1) = 1\n" <<
            "\tName = \"ZeroF\"\n" <<
            "\tForce 3 = 0\n" <<
            "\tForce 2 = 0\n" <<
            "\tForce 1 = 0\n" <<
            "End\n\n" ;*/

    simfile <<
            "Boundary Condition 2\n" <<
            "\tTarget Boundaries(1) = 2\n" <<
            "\tName = \"ZeroU\"\n" <<
            "\tDisplacement 3 = 0\n" <<
            "\tDisplacement 2 = 0\n" <<
            "\tDisplacement 1 = 0\n" <<
            "End\n\n" ;

    simfile <<
            "Boundary Condition 3\n" <<
            "\tTarget Boundaries(1) = 3\n" <<
            "\tName = \"Pull\"\n" <<
            "\t!A 'Force' in a boundary condition is really a pressure\n" <<
            "\t!This means that the values given here are in Pascal\n" <<
            "\t!Don't ask me why they call it a 'Force'... :/\n" <<
            "\tForce 1 = " << manipulateX << "\n" <<
            "\tForce 2 = " << manipulateY << "\n" <<
            "\tForce 3 = " << manipulateZ << "\n" <<
            "End\n\n" ;

}
