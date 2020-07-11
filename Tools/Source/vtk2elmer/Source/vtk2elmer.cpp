/******************************************************
 * This program takes a .vtu unstructured volume mesh as input.
 * It then outputs everything necessary to perform a finit element
 * simulation, optionally adding
 * The result files can be passed to the program vtk2elmer to convert
 * to mesh.* files which can be simulated in Elmer.
 *******************************************************/

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
#include <vtkCellLocator.h>
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
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkWarpVector.h>
#include <sstream>
#include <vtkDirectory.h>

#include "writeElmerFiles.h"

void runSimulation( std::string path );

inline bool fileExists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

int main(int argc, char* argv[])
{
	//parse command line arguments
    if(argc < 2)
    {
      std::cerr << "Usage: " << argv[0]
              << " Foldername [YoungsModulus]" << std::endl;
      return EXIT_FAILURE;
    }

    std::string foldername = argv[1];

    float youngsModulus = 3000.0;
    if(argc > 2)
    {
      std::stringstream sstr;
      sstr << argv[2];
      sstr >> youngsModulus;
    }

    // -----------------------------------------------------
    // Read the input file. This file should contain an unstructured mesh:
    vtkUnstructuredGrid* ugrid;
    vtkSmartPointer<vtkXMLUnstructuredGridReader> readerXML =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    vtkSmartPointer<vtkUnstructuredGridReader> reader =
        vtkSmartPointer<vtkUnstructuredGridReader>::New();
    if( fileExists( (foldername + "/simulation.vtu").c_str() ) )
    {
        std::cout << "File input.vtu found. Opening:" << std::endl;
        readerXML->SetFileName((foldername + "/simulation.vtu").c_str());
        readerXML->Update();
        ugrid = readerXML->GetOutput();
        std::cout << "\tRead." << std::endl;
    } else if( fileExists( (foldername + "/simulation.vtk").c_str() ) ) {
        std::cout << "File input.vtk found. Opening:" << std::endl;
        reader->SetFileName((foldername + "/simulation.vtk").c_str());
        std::cout << reader->GetFileName() << std::endl;
        reader->Update();
        ugrid = reader->GetOutput();
        std::cout << "\tRead." << std::endl;
    } else {
        throw("Cannot find simulation.vtu or simulation.vtk in given input folder!");
    }


    // -----------------------------------------------------
    // Read the surface input file. This file should contain the surface
    // of the same mesh, along with arrays which represent the boundary conditions.
    /*vtkSmartPointer<vtkXMLUnstructuredGridReader> readerSurf =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    readerSurf->SetFileName((foldername + "/surface.vtu").c_str());
    readerSurf->Update();
    vtkSmartPointer<vtkUnstructuredGrid> surface = readerSurf->GetOutput();*/

    // -----------------------------------------------------
    // Write mesh files
    std::string outputFolder = foldername + "/simulation";
    vtkDirectory::MakeDirectory( outputFolder.c_str() );
    writeElmerFiles( ugrid, outputFolder );

    // -----------------------------------------------------
    // Write simulation files
    writeSimulation( outputFolder, youngsModulus );

    // -----------------------------------------------------
    // Run the simulation
    //runSimulation( outputFolder );

    return EXIT_SUCCESS;
}

/*void runSimulation( std::string path )
{
    std::cout << "Running simulation: " << path << std::endl;
    std::stringstream cmd;
    cmd << "mkdir -p " << path << "/result";
    std::cout << "\t" << cmd.str() << std::endl;

    int result = system( cmd.str().c_str() );
    std::cout << "Result: " << result << std::endl;

    cmd.str("");
    cmd << "cd " << path << "; ElmerSolver case.sif";// > simLog.txt";
    std::cout << "\t" << cmd.str() << std::endl;

    result = system( cmd.str().c_str() );
    std::cout << "Result: " << result << std::endl;

    //os.execute('cd ' .. self.mesh.path .. '; ElmerSolver case.sif > /dev/null' ) then
}*/
