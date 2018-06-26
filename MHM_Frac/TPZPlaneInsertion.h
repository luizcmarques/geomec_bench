//
//  TPZPlaneInsertion.hpp
//  FractureIntersection
//
//  Created by Philippe Devloo on 26/06/18.
//

#ifndef TPZPlaneInsertion_hpp
#define TPZPlaneInsertion_hpp

#include <stdio.h>
#include "pzgmesh.h"
#include "pzfmatrix.h"
#include "FractureIntersectionConfig.h"

class TPZPlaneIntersection
{
    TPZGeoMesh LocalMesh;
    
    FractureIntersectionConfig *Config;
    
    /// indices of the nodes in the global mesh
    TPZStack<int64_t> GlobalIndex;
    
    /// map of the indices of the global mesh into the local mesh
    std::map<int64_t, int64_t> GlobalToLocal;
    
    int PlotCount = 0;
    
    /// compute the number of points in ksi and eta to cover the node cloud
    void ComputeNumPoints(TPZVec<int> &npoints);
    
    /// insert the corner points of the fracture plane, projecting them to ribs if necessary
    void InsertCornerPoints();
public:
    
    ///
    TPZPlaneIntersection(TPZFMatrix<REAL> &cornerco, FractureIntersectionConfig *config) : Config(config)
    {
        LocalMesh.SetDimension(config->gmesh->Dimension()-1);
        LocalMesh.NodeVec().Resize(4);
        for(int i=0; i<4; i++)
        {
            TPZManVector<REAL,3> co(3);
            for(int j=0; j<3; j++) co[j] = cornerco(j,i);
            LocalMesh.NodeVec()[i].Initialize(co, LocalMesh);
        }
        TPZManVector<int64_t,4> indices(4);
        for (int i=0; i<4; i++) indices[i] = i;
        int64_t index;
        LocalMesh.CreateGeoElement(EQuadrilateral, indices, 1, index);
    }
    
    /// Find points in the original mesh that are close to the fracture plane
    void FindClosePoints();
    
    /// Find the points on one dimensional elements that cut the plane
    void FindRibIntersections();
    
    /// Find fracture border intersections
    void ComputeFractureBorderIntersection();
};
#endif /* TPZPlaneInsertion_hpp */
