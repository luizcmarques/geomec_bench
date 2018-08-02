//
//  TPZFractureNeighborData.cpp
//  Benchmark0a
//
//  Created by Pablo Carvalho on 02/08/18.
//

#include "TPZFractureNeighborData.h"

/// Default constructor
TPZFractureNeighborData::TPZFractureNeighborData(){
    
}

/// Default desconstructor
TPZFractureNeighborData::~TPZFractureNeighborData(){
    
}

/// Copy constructor
TPZFractureNeighborData::TPZFractureNeighborData(TPZFractureNeighborData & other){
    m_pivot_indexes         = other.m_pivot_indexes;
    m_non_pivot_indexes     = other.m_non_pivot_indexes;
    m_gel_left_indexes      = other.m_gel_left_indexes;
    m_gel_right_indexes     = other.m_gel_right_indexes;
}


/// Set node pivots and non pivots
void TPZFractureNeighborData::SetPivotIndexes(std::set<int64_t> & pivots, std::set<int64_t> & non_pivots){
    m_pivot_indexes     = pivots;
    m_non_pivot_indexes =   non_pivots;
}

/// Set geometric fracture indexes
void TPZFractureNeighborData::SetFractureIndexes(std::set<int64_t> & fracture_indexes){
    
}

/// Set geometric indexes for left and right sides
void TPZFractureNeighborData::SetLeftRightIndexes(std::set<int64_t> & left_indexes, std::set<int64_t> & right_indexes){
    m_gel_left_indexes  = left_indexes;
    m_gel_right_indexes =  right_indexes;
}

/// Get node pivots
std::set<int64_t> & TPZFractureNeighborData::GetPivotIndexes(){
    return m_pivot_indexes;
}

/// Get node non pivots
std::set<int64_t> & TPZFractureNeighborData::GetNonPivotIndexes(){
    return m_non_pivot_indexes;
}

/// Get geometric fracture indexes
std::set<int64_t> & TPZFractureNeighborData::GetFractureIndexes(){
    return m_fracture_indexes;
}

/// Get geometric indexes for left
std::set<int64_t> & TPZFractureNeighborData::GetLeftIndexes(){
    return m_gel_left_indexes;
}

/// Get geometric indexes for right
std::set<int64_t> & TPZFractureNeighborData::GetRightIndexes(){
    return m_gel_left_indexes;
}

/// Constructor based on a computational mesh and fracture material id
TPZFractureNeighborData::TPZFractureNeighborData(TPZCompMesh & cmesh, int fracture_id){
    
    TPZGeoMesh *gmesh = cmesh.Reference();
    gmesh->ResetReference();
    cmesh.LoadReferences();
    int64_t ncel = cmesh.NElements();
    
    for (int64_t el=0; el<ncel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        
        // Filtering elements by fracture material identifier
        TPZGeoEl *gel = cel->Reference();
        if(gel->MaterialId() != fracture_id || gel->Dimension()!=gmesh->Dimension()-1){
            continue;
        }
        
        m_fracture_indexes.insert(gel->Index());
        
        unsigned int n_corner_sides = gel->NCornerNodes();
        for (unsigned int i_side = 0; i_side < n_corner_sides; i_side++) {
            
            TPZStack<TPZCompElSide> cel_side_vec;
            TPZGeoElSide gel_corner_side(gel,i_side);
            gel_corner_side.EqualLevelCompElementList(cel_side_vec, 0, 0);
            
            unsigned int fracture_neighbors_counter = 0;
            unsigned int d_minus_one_neighbors_counter = 0;
            
            unsigned int n_neighbors = cel_side_vec.size();
            for (unsigned int i_neighbor = 0; i_neighbor < n_neighbors; i_neighbor++) {
                TPZGeoEl * gel_neighbor = cel_side_vec[i_neighbor].Element()->Reference();
#ifdef PZDEBUG
                if (!gel_neighbor) {
                    DebugStop();
                }
#endif
                if (gel_neighbor->MaterialId() == fracture_id) {
                    fracture_neighbors_counter++;
                    m_fracture_indexes.insert(gel_neighbor->Index());
                }
                
                if (gel_neighbor->Dimension() == gmesh->Dimension()-1) {
                    d_minus_one_neighbors_counter++;
                }
                
                // inserting volumetric neighbors elements
                if (gel_neighbor->Dimension() == gmesh->Dimension()) {
                    if(m_gel_left_indexes.size() == 0 || m_gel_right_indexes.size() == 0){
                        if (m_gel_left_indexes.size() != 0 && m_gel_left_indexes.count(gel_neighbor->Index()) == 0 ) {
                            m_gel_right_indexes.insert(gel_neighbor->Index());
                        }
                        else{
                            m_gel_left_indexes.insert(gel_neighbor->Index());
                        }
                    }else{
                        bool insert_element_Q = m_gel_left_indexes.count(gel_neighbor->Index()) == 0 &&  m_gel_right_indexes.count(gel_neighbor->Index()) == 0;
                        if (insert_element_Q) {
                            InsertNextNeighbor(gel_neighbor);
                        }
 
                    }
                }
                
                
            }
            
            bool is_non_pivot_Q = (fracture_neighbors_counter == 0 && d_minus_one_neighbors_counter == 1);
            if (is_non_pivot_Q) {
                m_non_pivot_indexes.insert(gel->NodeIndex(i_side));
            }else{
                m_pivot_indexes.insert(gel->NodeIndex(i_side));
            }
        }
    }
    
    int left_id = 100;
    for (auto it=m_gel_left_indexes.begin(); it != m_gel_left_indexes.end(); ++it){
        int64_t iel = *it;
        gmesh->Element(iel)->SetMaterialId(left_id);
    }
    
    int right_id = 200;
    for (auto it=m_gel_right_indexes.begin(); it != m_gel_right_indexes.end(); ++it){
        int64_t iel = *it;
        gmesh->Element(iel)->SetMaterialId(right_id);
    }
    
    std::ofstream filegvtk("Geometry_labels.vtk"); //Impressão da malha geométrica (formato vtk)
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
    
    
}

void TPZFractureNeighborData::InsertNextNeighbor(TPZGeoEl * gel){
    
    int n_corners = gel->NCornerNodes();
    
    for (unsigned int i_edge_side = n_corners; i_edge_side < gel->NSides()-1; i_edge_side++) {
        TPZGeoElSide edge_neighbor = gel->Neighbour(i_edge_side);
        TPZGeoEl * gel_edge_neighbor = edge_neighbor.Element();
        if(m_gel_left_indexes.count(gel_edge_neighbor->Index()) == 1){
            m_gel_left_indexes.insert(gel->Index());
        }
        if(m_gel_right_indexes.count(gel_edge_neighbor->Index()) == 1){
            m_gel_right_indexes.insert(gel->Index());
        }
    
    }
    
}


