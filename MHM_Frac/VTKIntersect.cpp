//
//  VTKIntersect.cpp
//  Benchmark0a
//
//  Created by Philippe Devloo on 21/09/18.
//

#include "VTKIntersect.h"

/*=========================================================================
 
 Program:   Visualization Toolkit
 Module:    TestBooleanOperationPolyDataFilter.cxx
 
 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.
 =========================================================================*/

#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkReverseSense.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkDataSetMapper.h>
#include "vtkCamera.h"

#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"

void GenerateVTKInput()
{
    {
        TPZGmshReader gmsh;
        gmsh.fPZMaterialId[3]["domain"] = 1;

        TPZGeoMesh *gmesh = 0;
#ifdef MACOSX
        gmesh = gmsh.GeometricGmshMesh("../cube.msh");
#else
        gmesh = gmsh.GeometricGmshMesh("cube.msh");
#endif
        std::ofstream out("cube.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        delete gmesh;
    }
    {
        TPZGmshReader gmsh;
        gmsh.fPZMaterialId[2]["domain"] = 1;

        TPZGeoMesh *gmesh = 0;
#ifdef MACOSX
        gmesh = gmsh.GeometricGmshMesh("../plane.msh");
#else
        gmesh = gmsh.GeometricGmshMesh("plane.msh");
#endif
        std::ofstream out("plane.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        delete gmesh;
    }

}

vtkActor* GetBooleanOperationActor( double x, int operation )
{

    operation = vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION;
    double centerSeparation = .25;
    
    vtkSmartPointer<vtkUnstructuredGridReader> sphere1 =
    vtkSmartPointer<vtkUnstructuredGridReader>::New();
    sphere1->SetFileName("cube.vtk");
    
    
    vtkSmartPointer<vtkUnstructuredGrid> ugrd1 = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ugrd1 = sphere1->GetOutput();
    
    vtkGeometryFilter *gf1 = vtkGeometryFilter::New();
    gf1->SetInputData(sphere1->GetOutput());
    
    vtkSmartPointer<vtkUnstructuredGridReader> sphere2 =
    vtkSmartPointer<vtkUnstructuredGridReader>::New();
    sphere2->SetFileName("plane.vtk");
    
    vtkGeometryFilter *gf2 = vtkGeometryFilter::New();
    gf2->SetInputData(sphere2->GetOutput());
    
    vtkSmartPointer<vtkIntersectionPolyDataFilter> intersection =
    vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();
    
    intersection->SetInputData(0, gf1->GetOutput());
    intersection->SetInputData( 1, gf2->GetOutput());
    
//    intersection->Update();
    
    //intersection->Get
    
    vtkSmartPointer<vtkDistancePolyDataFilter> distance =
    vtkSmartPointer<vtkDistancePolyDataFilter>::New();
    distance->SetInputConnection( 0, intersection->GetOutputPort( 1 ) );
    distance->SetInputConnection( 1, intersection->GetOutputPort( 2 ) );
    
    vtkSmartPointer<vtkThreshold> thresh1 =
    vtkSmartPointer<vtkThreshold>::New();
    thresh1->AllScalarsOn();
    thresh1->SetInputArrayToProcess
    ( 0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Distance" );
    thresh1->SetInputConnection( distance->GetOutputPort( 0 ) );
    
    vtkSmartPointer<vtkThreshold> thresh2 =
    vtkSmartPointer<vtkThreshold>::New();
    thresh2->AllScalarsOn();
    thresh2->SetInputArrayToProcess
    ( 0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "Distance" );
    thresh2->SetInputConnection( distance->GetOutputPort( 1 ) );
    
    if ( operation == vtkBooleanOperationPolyDataFilter::VTK_UNION )
    {
        thresh1->ThresholdByUpper( 0.0 );
        thresh2->ThresholdByUpper( 0.0 );
    }
    else if ( operation == vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION )
    {
        thresh1->ThresholdByLower( 0.0 );
        thresh2->ThresholdByLower( 0.0 );
    }
    else // Difference
    {
        thresh1->ThresholdByUpper( 0.0 );
        thresh2->ThresholdByLower( 0.0 );
    }
    
    vtkSmartPointer<vtkDataSetSurfaceFilter> surface1 =
    vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surface1->SetInputConnection( thresh1->GetOutputPort() );
    
    vtkSmartPointer<vtkDataSetSurfaceFilter> surface2 =
    vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surface2->SetInputConnection( thresh2->GetOutputPort() );
    
    vtkSmartPointer<vtkReverseSense> reverseSense =
    vtkSmartPointer<vtkReverseSense>::New();
    reverseSense->SetInputConnection( surface2->GetOutputPort() );
    if ( operation == 2 ) // difference
    {
        reverseSense->ReverseCellsOn();
        reverseSense->ReverseNormalsOn();
    }
    
    vtkSmartPointer<vtkAppendPolyData> appender =
    vtkSmartPointer<vtkAppendPolyData>::New();
    appender->SetInputData(0, gf1->GetOutput());
    
//
//    appender->SetInputConnection( surface1->GetOutputPort() );
//    if ( operation == 2)
//    {
//        appender->AddInputConnection( reverseSense->GetOutputPort() );
//    }
//    else
//    {
//        appender->AddInputConnection( surface2->GetOutputPort() );
//    }
    
    vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkDataSetMapper * mppr = vtkDataSetMapper::New();
    mppr->SetInputData(ugrd1);
    
    mapper->ScalarVisibilityOff();
    
    vtkActor *actor = vtkActor::New();
    actor->SetMapper( mppr);//mapper );
    
    return actor;
}

int TestBooleanOperationPolyDataFilter2(int, char *[])
{
    vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
    
    vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer( renderer );
    
    vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renWinInteractor->SetRenderWindow( renWin );
    
    vtkActor *unionActor =
    GetBooleanOperationActor( -2.0, vtkBooleanOperationPolyDataFilter::VTK_UNION );
    renderer->AddActor( unionActor );
    unionActor->Delete();
    
    vtkActor *intersectionActor =
    GetBooleanOperationActor(  0.0, vtkBooleanOperationPolyDataFilter::VTK_INTERSECTION );
    renderer->AddActor( intersectionActor );
    intersectionActor->Delete();
    
    vtkActor *differenceActor =
    GetBooleanOperationActor(  2.0, vtkBooleanOperationPolyDataFilter::VTK_DIFFERENCE );
    renderer->AddActor( differenceActor );
    differenceActor->Delete();
    
    renWin->Render();
    renWinInteractor->Start();
    
    return EXIT_SUCCESS;
}

void VTKWindow()
{
    // The renderer generates the image
    // which is then displayed on the render window.
    // It can be thought of as a scene to which the actor is added
    vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
    vtkActor * act = GetBooleanOperationActor(0, 1);
    renderer->AddActor(act);
    renderer->SetBackground(0.1, 0.2, 0.4);
    // Zoom in a little by accessing the camera and invoking its "Zoom" method.
    renderer->ResetCamera();
    renderer->GetActiveCamera()->Zoom(1.5);
    
    // The render window is the actual GUI window
    // that appears on the computer screen
    vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetSize(2000, 2000);
    renderWindow->AddRenderer(renderer);
    
    // The render window interactor captures mouse events
    // and will perform appropriate camera or actor manipulation
    // depending on the nature of the events.
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    
    // This starts the event loop and as a side effect causes an initial render.
    renderWindowInteractor->Start();
    
    
}
