//
//  FractureIntersectionConfig.hpp
//  Benchmark0a
//
//  Created by Philippe Devloo on 22/06/18.
//

#ifndef FractureIntersectionConfig_hpp
#define FractureIntersectionConfig_hpp

#include <stdio.h>
#include <set>
#include "TPZPointCloud.h"
#include "pzgeoelside.h"

class TPZGeoMesh;

struct FractureIntersectionConfig
{
    std::set<int> materialids;
    std::set<int> bcmaterialids;
    
    TPZGeoMesh *gmesh;
    
    TPZPointCloud pointcloud;
    
    int volumematid = 1;
    
    int MHMatId = 101;
    
    int RibMatId = 102;
    
    /// create a "dimension" element neighbour of higher dimensional elements
    void CreateSurfaceMesh(int dimension);
    
    /// return all the sides of gel of a given dimension
    static void GetSides(TPZGeoEl *gel, int dimension, TPZStack<int> &sides);
    
    /// verify of an element of the side dimension is neighbour of the element/side
    static bool HasElementNeighbour(TPZGeoElSide &gelside);
    
    /// compute the material id of the element
    int ComputeMaterialId(TPZGeoElSide &gelside);
};
#endif /* FractureIntersectionConfig_hpp */
