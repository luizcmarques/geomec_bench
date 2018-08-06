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
    m_fracture_id           = other.m_fracture_id;
    m_geometry              = other.m_geometry;
    m_pivot_indexes         = other.m_pivot_indexes;
    m_non_pivot_indexes     = other.m_non_pivot_indexes;
    m_gel_left_indexes      = other.m_gel_left_indexes;
    m_gel_right_indexes     = other.m_gel_right_indexes;
}

/// Set fracture Identifier
void TPZFractureNeighborData::SetFractureIdentifier(int fracture_id){
    m_fracture_id = fracture_id;
}

/// Set node pivots and non pivots
void TPZFractureNeighborData::SetPivotIndexes(std::vector<TPZGeoElSide> & pivots, std::vector<TPZGeoElSide> & non_pivots){
    m_pivot_indexes     = pivots;
    m_non_pivot_indexes =   non_pivots;
    DebugStop();
}

/// Set geometric fracture indexes
void TPZFractureNeighborData::SetFractureIndexes(std::vector<int64_t> & fracture_indexes){
    DebugStop();
}

/// Set geometric indexes for left and right sides
void TPZFractureNeighborData::SetLeftRightIndexes(std::vector<int64_t> & left_indexes, std::vector<int64_t> & right_indexes){
    m_gel_left_indexes  = left_indexes;
    m_gel_right_indexes =  right_indexes;
}

/// Get fracture Identifier
int & TPZFractureNeighborData::GetFractureIdentifier(){
    return m_fracture_id;
}

/// Get node pivots
std::vector<TPZGeoElSide> & TPZFractureNeighborData::GetPivotIndexes(){
    return m_pivot_indexes;
}

/// Get node non pivots
std::vector<TPZGeoElSide> & TPZFractureNeighborData::GetNonPivotIndexes(){
    return m_non_pivot_indexes;
}

/// Get geometric fracture indexes
std::vector<int64_t> & TPZFractureNeighborData::GetFractureIndexes(){
    return m_fracture_indexes;
}

/// Get geometric indexes for left
std::vector<int64_t> & TPZFractureNeighborData::GetLeftIndexes(){
    return m_gel_left_indexes;
}

/// Get geometric indexes for right
std::vector<int64_t> & TPZFractureNeighborData::GetRightIndexes(){
    return m_gel_left_indexes;
}

/// Constructor based on a computational mesh and fracture material id
TPZFractureNeighborData::TPZFractureNeighborData(TPZGeoMesh * geometry, int fracture_id){
    m_geometry = geometry;
    m_fracture_id = fracture_id;
    
    std::ofstream filegvtk_b("Geometry_labels.vtk"); //Impressão da malha geométrica (formato vtk)
    TPZVTKGeoMesh::PrintGMeshVTK(m_geometry, filegvtk_b,true);
    
    IdentifyPivotsAndNonPivots();
    
    ComputeLeftAndRightSeeds();
    
    int dimension = m_geometry->Dimension();
    /// Computations based on left data
    for (auto it=m_pivot_indexes.begin(); it != m_pivot_indexes.end(); ++it){
        
        TPZGeoElSide node_side = *it;
        
//        TPZStack<TPZGeoElSide> egde_neighbors;
//        node_side.Element()->AllHigherDimensionSides(node_side.Side(), dimension, egde_neighbors);
        
        TPZStack<TPZGeoElSide> all_neighbors;
        node_side.AllNeighbours(all_neighbors);
        
//        for (int i = 0 ; i < all_neighbors.size(); i++) {
//            std::cout << "i element = "<< i << std::endl;
//            all_neighbors[i].Element()->Print();
//            std::cout << std::endl;
//        }
        
//        InsertLeftSideIndexes(all_neighbors);
        InsertRightSideIndexes(all_neighbors);
        int aka = 0;
//        break;
    }
    

//    int64_t ncel = geometry->NElements();
//    for (int64_t el=0; el<ncel; el++) {
//
//        // Filtering elements by fracture material identifier
//        TPZGeoEl *gel = geometry->Element(el);
//        if(gel->MaterialId() != m_fracture_id || gel->Dimension()!=geometry->Dimension()-1){
//            continue;
//        }
//
//        m_fracture_indexes.insert(gel->Index());
//
//        /// identify left and right seeds
//        TPZStack<TPZCompElSide> vols_side_vec;
//        TPZGeoElSide fracture_side(gel,gel->NSides()-1);
//        fracture_side.EqualLevelCompElementList(vols_side_vec, 0, 0);
//
//        if (vols_side_vec.size() != 2) {
//            DebugStop();
//        }
//
//        TPZGeoEl * gel_left = vols_side_vec[0].Element()->Reference();
//        TPZGeoEl * gel_right = vols_side_vec[1].Element()->Reference();
//        m_gel_left_indexes.insert(gel_left->Index());
//        m_gel_right_indexes.insert(gel_right->Index());
//
//
//        unsigned int n_corner_sides = gel->NCornerNodes();
//        for (unsigned int i_side = 0; i_side < n_corner_sides; i_side++) {
//
//
//            TPZStack<TPZCompElSide> cel_side_vec;
//            TPZGeoElSide gel_corner_side(gel,i_side);
//            gel_corner_side.EqualLevelCompElementList(cel_side_vec, 0, 0);
//
//            unsigned int fracture_neighbors_counter = 0;
//            unsigned int d_minus_one_neighbors_counter = 0;
//
//            unsigned int n_neighbors = cel_side_vec.size();
//            for (unsigned int i_neighbor = 0; i_neighbor < n_neighbors; i_neighbor++) {
//                TPZGeoEl * gel_neighbor = cel_side_vec[i_neighbor].Element()->Reference();
//#ifdef PZDEBUG
//                if (!gel_neighbor) {
//                    DebugStop();
//                }
//#endif
//                if (gel_neighbor->MaterialId() == fracture_id) {
//                    fracture_neighbors_counter++;
//                    m_fracture_indexes.insert(gel_neighbor->Index());
//                }
//
//                if (gel_neighbor->Dimension() == geometry->Dimension()-1) {
//                    d_minus_one_neighbors_counter++;
//                }
//
//                if (gel_neighbor->Dimension() == geometry->Dimension()-2) {
//                    break;
//                }
//
//
//
//                // inserting volumetric neighbors elements
//                if (gel_neighbor->Dimension() == geometry->Dimension()) {
//                    if(m_gel_left_indexes.size() == 0 || m_gel_right_indexes.size() == 0){
//                        if (m_gel_left_indexes.size() != 0 && m_gel_left_indexes.count(gel_neighbor->Index()) == 0 ) {
//                            m_gel_right_indexes.insert(gel_neighbor->Index());
//                        }
//                        else{
//                            m_gel_left_indexes.insert(gel_neighbor->Index());
//                        }
//
//                    }else{
//                        bool insert_element_Q = m_gel_left_indexes.count(gel_neighbor->Index()) == 0 &&  m_gel_right_indexes.count(gel_neighbor->Index()) == 0;
//                        if (insert_element_Q) {
//                            InsertNextNeighbor(gel_neighbor);
//                        }
//                    }
//                }
//
//            }
//
//        }
//    }
    
    int left_id = 100;
    for (auto it=m_gel_left_indexes.begin(); it != m_gel_left_indexes.end(); ++it){
        int64_t iel = *it;
        geometry->Element(iel)->SetMaterialId(left_id);
    }
    
    int right_id = 200;
    for (auto it=m_gel_right_indexes.begin(); it != m_gel_right_indexes.end(); ++it){
        int64_t iel = *it;
        geometry->Element(iel)->SetMaterialId(right_id);
    }
    
    std::ofstream filegvtk("Geometry_labels.vtk"); //Impressão da malha geométrica (formato vtk)
    TPZVTKGeoMesh::PrintGMeshVTK(geometry, filegvtk,true);
    
}

void TPZFractureNeighborData::IdentifyPivotsAndNonPivots(){
    
    int64_t n_gel = m_geometry->NElements();
    for (int64_t iel=0; iel< n_gel; iel++) {
        
        // Filtering elements by fracture material identifier
        TPZGeoEl *gel = m_geometry->Element(iel);
        if(gel->MaterialId() != m_fracture_id || gel->Dimension()!=m_geometry->Dimension()-1){
            continue;
        }
        bool has_point_neighbor_Q = false;
        unsigned int n_corner_sides = gel->NCornerNodes();
        for (unsigned int i_side = 0; i_side < n_corner_sides; i_side++) {
            
            TPZStack<TPZGeoElSide> all_neighbors;
            TPZGeoElSide gel_corner_side(gel,i_side);
            gel_corner_side.AllNeighbours(all_neighbors);
            
            unsigned int d_minus_one_neighbors_counter = 0;
            
            unsigned int n_neighbors = all_neighbors.size();
            for (unsigned int i_neighbor = 0; i_neighbor < n_neighbors; i_neighbor++) {
                TPZGeoEl * gel_neighbor = all_neighbors[i_neighbor].Element();
#ifdef PZDEBUG
                if (!gel_neighbor) {
                    DebugStop();
                }
#endif
                
                if (gel_neighbor->Dimension() == m_geometry->Dimension()-1) {
                    d_minus_one_neighbors_counter++;
                }
                
                if (gel_neighbor->Dimension() == m_geometry->Dimension()-2) {
                    has_point_neighbor_Q = true;
                }
                
            }
            bool is_non_pivot_Q = (d_minus_one_neighbors_counter == 0 && has_point_neighbor_Q);
            if (is_non_pivot_Q) {
                m_non_pivot_indexes.push_back(gel_corner_side);
            }else{
                m_pivot_indexes.push_back(gel_corner_side);
            }
        }
    }
    
}

void TPZFractureNeighborData::ComputeLeftAndRightSeeds(){
    

    int64_t n_gel = m_geometry->NElements();
    for (int64_t iel=0; iel < n_gel; iel++) {
        
        // Filtering elements by fracture material identifier
        TPZGeoEl *gel = m_geometry->Element(iel);
        if(gel->MaterialId() != m_fracture_id || gel->Dimension()!=m_geometry->Dimension()-1){
            continue;
        }
        
        m_fracture_indexes.push_back(gel->Index());
        
        /// identify left and right seeds
        TPZStack<TPZGeoElSide> all_neigh;
        TPZGeoElSide fracture_side(gel,gel->NSides()-1);
        fracture_side.AllNeighbours(all_neigh);
        if (all_neigh.size() != 2) {
            DebugStop();
        }
        
        TPZGeoEl * gel_left = all_neigh[0].Element();
        TPZGeoEl * gel_right = all_neigh[1].Element();
        m_gel_left_indexes.push_back(gel_left->Index());
        m_gel_right_indexes.push_back(gel_right->Index());
        break;
    }
    
    
}

/// Insert left neighbor element index
void TPZFractureNeighborData::InsertLeftSideIndexes(TPZStack<TPZGeoElSide> & neighbors){

#ifdef PZDEBUG
    if (neighbors.size() == 0 && m_gel_left_indexes.size() == 0) {
        DebugStop();
    }
#endif

    for (int i = 0; i < m_gel_left_indexes.size(); i++) {
        int64_t gel_left_index = m_gel_left_indexes[i];
        TPZGeoEl * gel_target = m_geometry->Element(gel_left_index);
        
        unsigned int n_neighbors = neighbors.size();
        for (unsigned int i_neighbor = 0; i_neighbor < n_neighbors; i_neighbor++) {
            TPZGeoEl * gel_candidate = neighbors[i_neighbor].Element();
            bool are_neighbors_Q = AreFaceNeighbors(gel_target, gel_candidate);
            if (are_neighbors_Q) {
                bool is_not_set_memeber_Q = true;
                for (int i = 0; i < m_gel_left_indexes.size(); i++) {
                    if(m_gel_left_indexes[i] == gel_candidate->Index()){
                        is_not_set_memeber_Q =  false;
                    }
                }
                if (is_not_set_memeber_Q) {
                    m_gel_left_indexes.push_back(gel_candidate->Index());
                }
            }
        }
    }
    
}

/// Insert left neighbor element index
void TPZFractureNeighborData::InsertRightSideIndexes(TPZStack<TPZGeoElSide> & neighbors){
    
#ifdef PZDEBUG
    if (neighbors.size() == 0 && m_gel_right_indexes.size() == 0) {
        DebugStop();
    }
#endif
    
    for (int i = 0; i < m_gel_right_indexes.size(); i++) {
        int64_t gel_right_index = m_gel_right_indexes[i];
        TPZGeoEl * gel_target = m_geometry->Element(gel_right_index);
        
        unsigned int n_neighbors = neighbors.size();
        for (unsigned int i_neighbor = 0; i_neighbor < n_neighbors; i_neighbor++) {
            TPZGeoEl * gel_candidate = neighbors[i_neighbor].Element();
            bool are_neighbors_Q = AreFaceNeighbors(gel_target, gel_candidate);
            if (are_neighbors_Q) {
                bool is_not_set_memeber_Q = true;
                for (int i = 0; i < m_gel_right_indexes.size(); i++) {
                    if(m_gel_right_indexes[i] == gel_candidate->Index()){
                        is_not_set_memeber_Q =  false;
                    }
                }
                if (is_not_set_memeber_Q) {
                    m_gel_right_indexes.push_back(gel_candidate->Index());
                }
            }
        }
    }

    
    
}

bool TPZFractureNeighborData::AreFaceNeighbors(TPZGeoEl * gel_target, TPZGeoEl * gel_candidate){

    bool are_neighbors_Q = false;
    unsigned int n_corner_sides = gel_candidate->NCornerNodes();
    unsigned int gel_candidate_side = gel_candidate->NSides() - 1;
    
    if (gel_target->Dimension() != gel_candidate->Dimension()) {
        return  are_neighbors_Q;
    }
    
    
    for (unsigned int i_edge_side = n_corner_sides; i_edge_side < gel_candidate_side; i_edge_side++) {
        TPZGeoElSide edge_neighbor = gel_candidate->Neighbour(i_edge_side);
        TPZGeoEl * gel_edge_neighbor = edge_neighbor.Element();
        if(gel_edge_neighbor->Index() == gel_target->Index() && gel_edge_neighbor->Dimension() == gel_target->Dimension()){
            are_neighbors_Q = true;
        }
    }
    
    return are_neighbors_Q;
}

void TPZFractureNeighborData::InsertNextNeighbor(TPZGeoEl * gel){
    
//    int n_corners = gel->NCornerNodes();
//
//    for (unsigned int i_edge_side = n_corners; i_edge_side < gel->NSides()-1; i_edge_side++) {
//        TPZGeoElSide edge_neighbor = gel->Neighbour(i_edge_side);
//        TPZGeoEl * gel_edge_neighbor = edge_neighbor.Element();
//        if(m_gel_left_indexes.count(gel_edge_neighbor->Index()) == 1){
//            m_gel_left_indexes.insert(gel->Index());
//        }
//        if(m_gel_right_indexes.count(gel_edge_neighbor->Index()) == 1){
//            m_gel_right_indexes.insert(gel->Index());
//        }
//
//    }
    DebugStop();
}


