/******************************************************
 * This program takes a folder as input and expects there to be an
 * input.vtu file inside it, which contains an unstructured volume mesh.
 * It then outputs everything necessary to perform a finite element
 * simulation, optionally adding random values.
 * The result files can be passed to the program vtk2elmer to convert
 * to mesh.* files which can be simulated in Elmer.
 *******************************************************/

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkAppendFilter.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>
#include <vtkPolyDataNormals.h>
#include <vtkDoubleArray.h>
#include <vtkMinimalStandardRandomSequence.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkSphere.h>
#include <vtkIdFilter.h>
#include <vtkDirectory.h>
#include <vtkCleanPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkVectorOperators.h>
#include <sstream>
#include <map>

#define PI 3.14159265

void GetCellCenter(vtkPolyData *poly, const unsigned int cellId, double center[]);
std::vector<vtkIdType> GetRandomSurfaceRegion( vtkUnstructuredGrid* ugrid, double radius,
                                               vtkSmartPointer<vtkMinimalStandardRandomSequence> rnd,
                                               int minimumRequiredCells,
                                               std::vector<vtkIdType> mask = std::vector<vtkIdType>() );

vtkVector3d getCellCenter( vtkCell* cell );

inline bool fileExists (const std::string& name) {
  ifstream f(name.c_str());
  return f.good();
}

vtkVector3d randomPerpendicularVector( vtkVector3d in, vtkSmartPointer<vtkMinimalStandardRandomSequence> rnd );

int main(int argc, char* argv[])
{

  //parse command line arguments
  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0]
              << " Folder [Seed]" << std::endl;
    return EXIT_FAILURE;
  }

  std::string foldername = argv[1];

  int seed = 1;
  if(argc > 2)
  {
    std::stringstream sstr;
    sstr << argv[2];
    sstr >> seed;
  }
  std::cout << "Using random seed: " << seed << std::endl;


  // -----------------------------------------------------
  // Read the input file. This file should contain an unstructured mesh:
  vtkUnstructuredGrid* ugrid;
  vtkSmartPointer<vtkXMLUnstructuredGridReader> readerXML =
      vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  vtkSmartPointer<vtkUnstructuredGridReader> reader =
      vtkSmartPointer<vtkUnstructuredGridReader>::New();
  if( fileExists( (foldername + "/input.vtu").c_str() ) )
  {
    std::cout << "File input.vtu found. Opening:" << std::endl;
    readerXML->SetFileName((foldername + "/input.vtu").c_str());
    readerXML->Update();
    ugrid = readerXML->GetOutput();
    std::cout << "\tRead." << std::endl;
  } else if( fileExists( (foldername + "/input.vtk").c_str() ) ) {
    std::cout << "File input.vtk found. Opening:" << std::endl;
    reader->SetFileName((foldername + "/input.vtk").c_str());
    reader->Update();
    ugrid = reader->GetOutput();
    std::cout << "\tRead." << std::endl;
  } else {
    throw("Cannot find input.vtu or input.vtk in given input folder!");
  }


  vtkSmartPointer<vtkMinimalStandardRandomSequence> rnd = vtkSmartPointer<vtkMinimalStandardRandomSequence>::New();
  rnd->SetSeed(seed);

  // -----------------------------------------------------
  // Calculate size of mesh:
  double bounds[6];
  ugrid->GetPoints()->GetBounds( bounds );
  double inputMeshSize = std::max( std::max( bounds[1]-bounds[0], bounds[3]-bounds[2] ),
      bounds[5]-bounds[4]);

  // -----------------------------------------------------
  // Add zero-displacement border:
  vtkSmartPointer<vtkFloatArray> zeroDisplacement = vtkSmartPointer<vtkFloatArray>::New();
  zeroDisplacement->SetNumberOfComponents(1);
  zeroDisplacement->SetNumberOfTuples( ugrid->GetNumberOfCells() );
  zeroDisplacement->SetName("zeroDisplacement");
  for( unsigned int i = 0; i < zeroDisplacement->GetNumberOfTuples(); i++ )
    zeroDisplacement->SetTuple1( i, 0 );
  // Limit the size of the zero-displacement border to be (at maximum) slightly larger than
  // half the mesh size:
  double radius = 0.025 + 0.03*rnd->GetValue()*inputMeshSize;
  int minimumRequiredCells = 3;
  std::vector<vtkIdType> IDsDisplacement = GetRandomSurfaceRegion( ugrid, radius, rnd, minimumRequiredCells);
  std::cout << "Assigning zero displacement boundary condition to " << IDsDisplacement.size() << " surface cells." << std::endl;
  for( int i = 0; i < IDsDisplacement.size(); i++)
  {
    vtkIdType id = IDsDisplacement[i];
    zeroDisplacement->SetTuple1( id, 1 );
  }
  ugrid->GetCellData()->AddArray( zeroDisplacement );

  // -----------------------------------------------------
  // Add random force (or random displacement) to some surface cells:
  vtkSmartPointer<vtkFloatArray> manipulation = vtkSmartPointer<vtkFloatArray>::New();
  manipulation->SetNumberOfComponents(3);
  manipulation->SetNumberOfTuples( ugrid->GetNumberOfCells() );
  manipulation->SetName("manipulation");
  for( unsigned int i = 0; i < manipulation->GetNumberOfTuples(); i++ )
    manipulation->SetTuple3(i, 0,0,0 );


  radius = 0.015 + 0.01*rnd->GetValue(); rnd->Next();
  minimumRequiredCells = 1;
  std::vector<vtkIdType> IDsManipulation = GetRandomSurfaceRegion( ugrid, radius, rnd, minimumRequiredCells, IDsDisplacement);

  double x = rnd->GetValue(); rnd->Next();
  double y = rnd->GetValue(); rnd->Next();
  double z = rnd->GetValue(); rnd->Next();
  vtkVector3d manipulationVec(x,y,z);

  // Calculate area of the random surface region:
  double area = 0;
  for( int i = 0; i < IDsManipulation.size(); i++)
  {
    vtkTriangle* t = (vtkTriangle*)ugrid->GetCell( IDsManipulation[i] );
    area += t->ComputeArea();
  }
  std::cout << "Area of manipulation surface region: ~" << area << " m^2" << std::endl;

  // rescale the manipulation
  double r = rnd->GetValue(); rnd->Next();
  double magnitudeForce = 1*r;
  double magnitude = magnitudeForce/area;    // random magnitude

  std::cout << "Random Force: " << magnitudeForce << " N" << std::endl;
  std::cout << "Resulting Pressure: " << magnitude << " Pa" << std::endl;

  manipulationVec = manipulationVec.Normalized()*magnitude;


  std::cout << "Assigning manipulation boundary condition to " << IDsManipulation.size() << " surface cells." << std::endl;
  for( int i = 0; i < IDsManipulation.size(); i++)
  {
    vtkIdType id = IDsManipulation[i];
    manipulation->SetTuple3( id, manipulationVec[0], manipulationVec[1], manipulationVec[2] );
  }
  ugrid->GetCellData()->AddArray( manipulation );

  /*for( int i = 0; i < region->GetNumberOfCells(); i++ )
       region->GetCell(i)->Print(std::cout);*/

  // -----------------------------------------------------
  // Write result:
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
      vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName((foldername + "/simulation.vtu").c_str());
  //writer->SetDataMode( vtkXMLWriter::Ascii );
  writer->SetInputData( ugrid );
  writer->Update();
  std::cout << "Wrote to " << (foldername + "/simulation.vtu") << std::endl;

  // Convert surface to unstructured grid using the appendFilter:
  /*vtkSmartPointer<vtkAppendFilter> appendFilter =
        vtkSmartPointer<vtkAppendFilter>::New();
    appendFilter->AddInputData(surface);
    appendFilter->Update();

    // -----------------------------------------------------
    // Write result:
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writerSurf =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writerSurf->SetFileName((foldername + "/surface.vtu").c_str());
    writerSurf->SetDataMode( vtkXMLWriter::Ascii );
    writerSurf->SetInputData( appendFilter->GetOutput() );
    writerSurf->Update();
    std::cout << "Wrote to " << (foldername + "/surface.vtu") << std::endl;*/

  return EXIT_SUCCESS;
}

void GetCellCenter(vtkPolyData* poly, const unsigned int cellId, double center[3])
{
  double pcoords[3] = {0,0,0};
  double weights[poly->GetMaxCellSize()];
  vtkCell* cell = poly->GetCell(cellId);
  int subId = cell->GetParametricCenter(pcoords);
  cell->EvaluateLocation(subId, pcoords, center, weights);
}

bool isInMask( const vtkIdType ID, std::vector<vtkIdType>* mask )
{
  for( int i = 0; i < mask->size(); i++ )
  {
    if(mask->at(i) == ID)
    {
      return true;
    }
  }
  return false;
}

std::vector<vtkIdType> GetRandomSurfaceRegion(vtkUnstructuredGrid* ugrid, double radius, vtkSmartPointer<vtkMinimalStandardRandomSequence> rnd, int minimumRequiredCells, std::vector<vtkIdType> mask )
{
  // Find the surface cells:
  //std::map<vtkIdType,Point> triangleCenters;
  std::vector<vtkIdType> triangleIDs;
  for( int i = 0; i < ugrid->GetNumberOfCells(); i++ )
  {
    vtkCell* cell = ugrid->GetCell(i);
    if( cell->GetCellType() == VTK_TRIANGLE )
    {
      //Point p = getCellCenter( cell );
      //triangleCenters[i] = p;
      triangleIDs.push_back( i );
    }
  }

  std::vector<vtkIdType> regionIDs;
  while( regionIDs.size() < minimumRequiredCells )
  {
    regionIDs.clear();

    double squaredRadius = radius*radius;

    unsigned int randomTriangleID;
    bool foundInMask = false;
    int retries = 0;
    // Make sure the randomTriangleID is NOT part of the mask
    do
    {
      randomTriangleID = rnd->GetValue()*triangleIDs.size(); rnd->Next();
      foundInMask = isInMask( randomTriangleID, &mask );
      retries ++;
    } while( foundInMask == true && retries < 100 );


    //Point center = triangleCenters[randomTriangleID];
    vtkCell* centerCell = ugrid->GetCell( triangleIDs[randomTriangleID] );
    vtkVector3d center = getCellCenter( centerCell );
    // Make sure the first cell in the region is the center cell:
    regionIDs.push_back(triangleIDs[randomTriangleID]);
    // Search for all closeby triangles:
    for( int i = 0; i < triangleIDs.size(); i++ )
    {
      if( i == randomTriangleID )
        continue;

      vtkCell* cell = ugrid->GetCell( triangleIDs[i] );
      vtkPoints* points = cell->GetPoints();
      double pos[3];
      bool outside = false;
      for( int j = 0; j < points->GetNumberOfPoints(); j++ )
      {
        points->GetPoint(j,pos);
        double dx = (pos[0]-center[0]);
        double dy = (pos[1]-center[1]);
        double dz = (pos[2]-center[2]);
        double squaredDist = dx*dx + dy*dy + dz*dz;
        if( squaredDist > squaredRadius )
        {
          outside = true;
          break;
        }
      }
      if( !outside )
      {
        regionIDs.push_back(triangleIDs[i]);
      }
    }
    /*std::map<vtkIdType, Point>::iterator it;
        for ( it = triangleCenters.begin(); it != triangleCenters.end(); it++ )
        {
            // Calculate the center of the triangle:
            Point curCenter = it->second;
            double dx = (curCenter.x-center.x);
            double dy = (curCenter.y-center.y);
            double dz = (curCenter.z-center.z);
            double squaredDist = dx*dx + dy*dy + dz*dz;
            if( squaredDist < squaredRadius )
            {
                regionIDs.push_back( it->first );
            }
        }*/
    if( regionIDs.size() < minimumRequiredCells )
    {
      std::cout << "Found " << regionIDs.size() << " cells, which is less than the minimum required (" << minimumRequiredCells << "). Retrying with larger radius." << std::endl;
      radius = radius + 0.01;
    }
  }

  return regionIDs;
}

vtkVector3d getCellCenter( vtkCell* cell )
{
  double paramCoords[2];
  cell->GetParametricCenter(paramCoords);
  double center[3];
  double weights[3];
  int subID = 0;
  cell->EvaluateLocation(subID,paramCoords,center,weights);
  vtkVector3d p(center[0],center[1],center[2]);
  return p;
}


vtkVector3d randomPerpendicularVector( vtkVector3d in, vtkSmartPointer<vtkMinimalStandardRandomSequence> rnd )
{
  in.Normalize();
  vtkVector3d unitX( 1, 0, 0 );
  vtkVector3d c = in.Cross(unitX);
  double ang = 2*PI*rnd->GetValue(); rnd->Next();
  // Rodrigues' rotation formula (partial, because last term is always zero in our case):
  vtkVector3d cRot = c*std::cos(ang) + in.Cross(c)*std::sin(ang);
  return cRot;
}
