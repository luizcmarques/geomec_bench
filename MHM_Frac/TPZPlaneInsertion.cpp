//
//  TPZPlaneInsertion.cpp
//  FractureIntersection
//
//  Created by Philippe Devloo on 26/06/18.
//

#include "TPZPlaneInsertion.h"
#include "pzvec_extras.h"
#include "pzgeoel.h"


/// compute the number of points in ksi and eta to cover the node cloud
void TPZPlaneIntersection::ComputeNumPoints(TPZVec<int> &npoints)
{
    TPZManVector<TPZManVector<REAL,3>,4> coords(4);
    for (int i=0; i<4; i++) {
        coords[i].Resize(3, 0.);
        LocalMesh.NodeVec()[i].GetCoordinates(coords[i]);
    }
    REAL distksi = std::max(dist(coords[0],coords[1]),dist(coords[3],coords[2]));
    REAL disteta = std::max(dist(coords[0],coords[3]),dist(coords[1],coords[2]));
    npoints[0] = distksi/Config->pointcloud.delxmin;
    npoints[1] = disteta/Config->pointcloud.delxmin;
}

/// Find points in the original mesh that are close to the fracture plane
void TPZPlaneIntersection::FindClosePoints()
{
    TPZManVector<int,2> npoints;
    ComputeNumPoints(npoints);
    for(int iksi=0; iksi<=npoints[0]; iksi++)
        for (int ieta=0; ieta<=npoints[1]; ieta++) {
            TPZManVector<REAL,2> point(2);
            point[0] = -1.+iksi*2./npoints[0];
            point[1] = -1.+ieta*2./npoints[1];
            TPZManVector<REAL,3> xco(3);
            LocalMesh.Element(0)->X(point,xco);
            int64_t close = Config->pointcloud.FindClosePoint(xco);
            if(close != -1 && GlobalToLocal.find(close) == GlobalToLocal.end())
            {
                // we need to insert the point
                TPZManVector<REAL,3> globalco(3);
                Config->gmesh->NodeVec()[close].GetCoordinates(globalco);
                int64_t index = LocalMesh.NodeVec().AllocateNewElement();
                GlobalToLocal[close] = index;
                LocalMesh.NodeVec()[index].Initialize(globalco, LocalMesh);
                GlobalIndex.Resize(index+1,-1);
                GlobalIndex[index] = close;
            }
        }
}

/// insert the corner points of the fracture plane, projecting them to ribs if necessary
void TPZPlaneIntersection::InsertCornerPoints()
{
    TPZManVector<REAL,3> localco(3);
    GlobalIndex.Resize(4, -1);
    for (int corner=0; corner<4; corner++) {
        LocalMesh.NodeVec()[corner].GetCoordinates(localco);
        int64_t globalindex;
        globalindex = Config->pointcloud.FindClosePoint(localco);
        if (globalindex != -1) {
            GlobalIndex[corner] = globalindex;
            GlobalToLocal[globalindex] = corner;
            continue;
        }
        TPZManVector<REAL,3> qsi(3);
        int64_t initialindex = -1;
        int targetdim = 3;
        TPZGeoEl *gelglob = Config->gmesh->FindElement(localco, qsi, initialindex, targetdim);
        /// check if qsi is interior in the element or close to a boundary
        /// insert the point in the localmesh and global mesh
        /// divide the element (volume, face or rib)
    }
}
