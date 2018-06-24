//
//  TPZFractureIntersection.cpp
//  Benchmark0a
//
//  Created by Philippe Devloo on 22/06/18.
//

#include <stdio.h>
#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
#include "TPZGmshReader.h"
#include "TPZVTKGeoMesh.h"
#include "FractureIntersectionConfig.h"

#include "pzgmesh.h"
#include "TPZPointCloud.h"


int main(int argc, char *argv[]) {
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    FractureIntersectionConfig config;
    {
        TPZGmshReader gmsh;
        gmsh.fPZMaterialId[2]["inflow"] = -1;
        gmsh.fPZMaterialId[2]["outflow"] = -2;
        gmsh.fPZMaterialId[2]["top"] = -3;
        gmsh.fPZMaterialId[2]["lateral"] = -4;
        gmsh.fPZMaterialId[3]["domain"] = config.volumematid;
        config.materialids.insert(1);
        config.bcmaterialids.insert(-1);
        config.bcmaterialids.insert(-2);
        config.bcmaterialids.insert(-3);
        config.bcmaterialids.insert(-4);
        TPZGeoMesh *gmesh = 0;
#ifdef MACOSX
        gmesh = gmsh.GeometricGmshMesh("../BoxReservoir.msh");
#else
        gmesh = gmsh.GeometricGmshMesh("BoxReservoir.msh");
#endif
        config.gmesh = gmesh;
        config.pointcloud.InsertGmesh(gmesh);
    }
    config.CreateSurfaceMesh(2);
    config.CreateSurfaceMesh(1);
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(config.gmesh, out);
    }

    return 0;
}
