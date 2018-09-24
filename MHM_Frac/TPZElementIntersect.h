//
//  TPZElementIntersect.hpp
//  FractureIntersection
//
//  Created by Philippe Devloo on 11/09/18.
//

#ifndef TPZElementIntersect_hpp
#define TPZElementIntersect_hpp

#include <stdio.h>
#include "pzmanvector.h"
#include "pzgeoelside.h"

/// this class is dedicated to compute the intersection of a line with an element
class TPZElementIntersect
{
public:
    /// points that define the line
    TPZManVector<REAL,3> first,last;
    
    TPZElementIntersect()
    {
        
    }
    
    TPZElementIntersect(const TPZElementIntersect &copy) : first(copy.first), last(copy.last)
    {
        
    }
    
    TPZElementIntersect &operator=(const TPZElementIntersect &copy)
    {
        first = copy.first;
        last = copy.last;
        return *this;
    }
    
    /// find the side which intersects the line starting from the given point
    bool FindSideIntersection(TPZVec<REAL> &init, TPZGeoElSide &intersect, TPZVec<REAL> &sidepar);
};

#endif /* TPZElementIntersect_hpp */
