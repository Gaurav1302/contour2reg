#ifndef WRITEELMERFILE_H
#define WRITEELMERFILE_H

#include <string>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

int typeForCell( vtkCell* cell );
int vtkCellType2ElmerCellType( int cellType );
void writeElmerFiles( vtkSmartPointer<vtkUnstructuredGrid> ugrid, std::string path );
void writeSimulation( std::string path, float youngsModulus );

#endif // WRITEELMERFILE_H
