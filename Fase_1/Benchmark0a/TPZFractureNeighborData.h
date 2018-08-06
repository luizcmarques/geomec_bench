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
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"

class TPZFractureNeighborData {
  
private:
    
    int m_fracture_id;
    
    TPZGeoMesh * m_geometry;
    
    std::vector<TPZGeoElSide> m_pivot_indexes;
    
    std::vector<TPZGeoElSide> m_non_pivot_indexes;

    std::vector<int64_t> m_fracture_indexes;
    
    std::vector<int64_t> m_gel_left_indexes;
    
    std::vector<int64_t> m_gel_right_indexes;
    
    /// Insert left or Right neighbor element index
    void InsertNextNeighbor(TPZGeoEl * gel);
    
    /// For the fracture material identifier computes the pivots and non pivots nodes
    void IdentifyPivotsAndNonPivots();
    
    /// Compute and insert the first left and right elements
    void ComputeLeftAndRightSeeds();
    
    /// Insert left neighbor element index
    void InsertLeftSideIndexes(TPZStack<TPZGeoElSide> & neighbors);
    
    /// Insert right neighbor element index
    void InsertRightSideIndexes(TPZStack<TPZGeoElSide> & neighbors);
    
    /// Verify if the geometric element candidate is a neighbor by a face of the geometric target element
    bool AreFaceNeighbors(TPZGeoEl * gel_target, TPZGeoEl * gel_candidate);
    
public:

    /// Default constructor
    TPZFractureNeighborData();
    
    /// Default desconstructor
    ~TPZFractureNeighborData();
    
    /// Copy constructor
    TPZFractureNeighborData(TPZFractureNeighborData & other);
    
    /// Constructor based on a computational mesh and fracture material id
    TPZFractureNeighborData(TPZGeoMesh * geometry, int fracture_id);

    /// Set fracture Identifier
    void SetFractureIdentifier(int fracture_id);
    
    /// Set node pivots and non pivots
    void SetPivotIndexes(std::vector<TPZGeoElSide> & pivots, std::vector<TPZGeoElSide> & non_pivots);
    
    /// Set geometric fracture indexes
    void SetFractureIndexes(std::vector<int64_t> & fracture_indexes);
    
    /// Set geometric indexes for left and right sides
    void SetLeftRightIndexes(std::vector<int64_t> & left_indexes, std::vector<int64_t> & right_indexes);
    
    /// Get fracture Identifier
    int & GetFractureIdentifier();
    
    /// Get node pivots
    std::vector<TPZGeoElSide> & GetPivotIndexes();

    /// Get node non pivots
    std::vector<TPZGeoElSide> & GetNonPivotIndexes();

    /// Get geometric fracture indexes
    std::vector<int64_t> & GetFractureIndexes();
    
    /// Get geometric indexes for left
    std::vector<int64_t> & GetLeftIndexes();
    
    /// Get geometric indexes for right
    std::vector<int64_t> & GetRightIndexes();
    
};


#endif /* TPZFractureNeighborData_h */
