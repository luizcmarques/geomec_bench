//
//  FractureIntersectionConfig.cpp
//  Benchmark0a
//
//  Created by Philippe Devloo on 22/06/18.
//

#include "FractureIntersectionConfig.h"
#include "pzgeoel.h"
#include "pzgeoelbc.h"

void FractureIntersectionConfig::CreateSurfaceMesh(int dimension)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->Dimension() <= dimension) {
            continue;
        }
        TPZStack<int> sides;
        GetSides(gel, dimension, sides);
        for(auto side:sides)
        {
            TPZGeoElSide gelside(gel,side);
            if(HasElementNeighbour(gelside)) continue;
            int locmatid = ComputeMaterialId(gelside);
            TPZGeoElBC gbc(gelside,locmatid);
        }
    }
}

/// return all the sides of gel of a given dimension
void FractureIntersectionConfig::GetSides(TPZGeoEl *gel, int dimension, TPZStack<int> &sides)
{
    int nsides = gel->NSides();
    for (int side = 0; side<nsides; side++) {
        if (gel->SideDimension(side) == dimension) {
            sides.Push(side);
        }
    }
}


/// verify of an element of the side dimension is neighbour of the element/side
bool FractureIntersectionConfig::HasElementNeighbour(TPZGeoElSide &gelside)
{
    if (gelside.Element()->Dimension() == gelside.Dimension()) {
        return true;
    }
    TPZGeoElSide neighbour = gelside.Neighbour();
    while (neighbour != gelside) {
        if (neighbour.Element()->Dimension() == neighbour.Dimension()) {
            return true;
        }
        neighbour = neighbour.Neighbour();
    }
    return false;
}

/// compute the material id of the element
int FractureIntersectionConfig::ComputeMaterialId(TPZGeoElSide &gelside)
{
    TPZGeoMesh *gmesh = gelside.Element()->Mesh();
    int sidedim = gelside.Dimension();
    int dim = gmesh->Dimension();
    if (sidedim == dim-1) {
        return MHMatId;
    }
    if (gelside.Element()->Dimension() == sidedim+1) {
        int matid = gelside.Element()->MaterialId();
        if (bcmaterialids.find(matid) != bcmaterialids.end()) {
            return matid;
        }
    }
    TPZGeoElSide neighbour = gelside.Neighbour();
    while (neighbour != gelside) {
        int matid = neighbour.Element()->MaterialId();
        if (bcmaterialids.find(matid) != bcmaterialids.end()) {
            return matid;
        }
        neighbour = neighbour.Neighbour();
        
    }
    return MHMatId;
}
