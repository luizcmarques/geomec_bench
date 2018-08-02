//
//  TPZFractureNeighborData.h
//  Benchmark0a
//  Class that stores fracture neighbor data in terms of geometric element indexes
//  Created by Pablo Carvalho on 02/08/18.
//

#ifndef TPZFractureNeighborData_h
#define TPZFractureNeighborData_h

#include <stdio.h>
#include <iostream>
#include <set>
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZVTKGeoMesh.h"

class TPZFractureNeighborData {
  
private:
    
    std::set<int64_t> m_pivot_indexes;
    
    std::set<int64_t> m_non_pivot_indexes;

    std::set<int64_t> m_fracture_indexes;
    
    std::set<int64_t> m_gel_left_indexes;
    
    std::set<int64_t> m_gel_right_indexes;
    
    /// Insert left or Right neighbor element index
    void InsertNextNeighbor(TPZGeoEl * gel);
    
public:

    /// Default constructor
    TPZFractureNeighborData();
    
    /// Default desconstructor
    ~TPZFractureNeighborData();
    
    /// Copy constructor
    TPZFractureNeighborData(TPZFractureNeighborData & other);
    
    /// Constructor based on a computational mesh and fracture material id
    TPZFractureNeighborData(TPZCompMesh & cmesh, int fracture_id);

    /// Set node pivots and non pivots
    void SetPivotIndexes(std::set<int64_t> & pivots, std::set<int64_t> & non_pivots);
    
    /// Set geometric fracture indexes
    void SetFractureIndexes(std::set<int64_t> & fracture_indexes);
    
    /// Set geometric indexes for left and right sides
    void SetLeftRightIndexes(std::set<int64_t> & left_indexes, std::set<int64_t> & right_indexes);
    
    /// Get node pivots
    std::set<int64_t> & GetPivotIndexes();

    /// Get node non pivots
    std::set<int64_t> & GetNonPivotIndexes();

    /// Get geometric fracture indexes
    std::set<int64_t> & GetFractureIndexes();
    
    /// Get geometric indexes for left
    std::set<int64_t> & GetLeftIndexes();
    
    /// Get geometric indexes for right
    std::set<int64_t> & GetRightIndexes();
    
};


#endif /* TPZFractureNeighborData_h */
