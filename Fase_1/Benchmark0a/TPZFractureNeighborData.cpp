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
    m_fracture_id               = other.m_fracture_id;
    m_geometry                  = other.m_geometry;
    m_boundaries_material_ids   = other.m_boundaries_material_ids;
    m_pivot_indexes             = other.m_pivot_indexes;
    m_non_pivot_indexes         = other.m_non_pivot_indexes;
    m_gel_left_indexes          = other.m_gel_left_indexes;
    m_gel_right_indexes         = other.m_gel_right_indexes;
}

/// Set fracture Identifier
void TPZFractureNeighborData::SetFractureIdentifier(int fracture_id){
    m_fracture_id = fracture_id;
}

/// Get fracture Identifier
int & TPZFractureNeighborData::GetFractureMaterialId(){
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
std::set<int64_t> & TPZFractureNeighborData::GetLeftIndexes(){
    return m_gel_left_indexes;
}

/// Get geometric indexes for right
std::set<int64_t> & TPZFractureNeighborData::GetRightIndexes(){
    return m_gel_left_indexes;
}

/// Constructor based on a computational mesh and fracture material id
TPZFractureNeighborData::TPZFractureNeighborData(TPZGeoMesh * geometry, int fracture_id, std::set<int> & boundaries_ids){
    m_geometry = geometry;
    m_fracture_id = fracture_id;
    m_boundaries_material_ids = boundaries_ids;
    ClassifyNeighboursofPivots();
}

void TPZFractureNeighborData::BuildFractureElements(){
    
    int64_t n_gel = m_geometry->NElements();
    for (int64_t iel=0; iel < n_gel; iel++) {

        // Filtering elements by fracture material identifier
        TPZGeoEl *gel = m_geometry->Element(iel);
        if(gel->MaterialId() != m_fracture_id || gel->Dimension()!=m_geometry->Dimension()-1 || gel->HasSubElement() == 1){
            continue;
        }
        
//        std::cout << "gel->HasSubElement() = " << gel->HasSubElement() << std::endl;
        m_fracture_indexes.push_back(gel->Index());
    }
    
    if(m_fracture_indexes.size() == 0){
        std::cout << "No fracture elements with id = " << m_fracture_id << "  were detected." << std::endl;
        DebugStop();
    }
    
}

void TPZFractureNeighborData::BuildPivotDataStructure(){
    
    for (auto ifrac : m_fracture_indexes) {
        
        // Filtering elements by fracture material identifier
        TPZGeoEl *gel = m_geometry->Element(ifrac);
        unsigned int n_corner_sides = gel->NCornerNodes();
        for (unsigned int i_side = 0; i_side < n_corner_sides; i_side++) {
            
            TPZStack<TPZGeoElSide> all_neighbors;
            TPZGeoElSide gel_corner_side(gel,i_side);
            gel_corner_side.AllNeighbours(all_neighbors);
            
            unsigned int d_minus_bc_neighbors_counter = 0;
            unsigned int d_minus_other_fracture_neighbors_counter = 0;
            unsigned int d_minus_same_fracture_neighbors_counter = 0;
            
            bool has_point_neighbor_Q = false;
            unsigned int n_neighbors = all_neighbors.size();
            for (unsigned int i_neighbor = 0; i_neighbor < n_neighbors; i_neighbor++) {
                TPZGeoEl * gel_neighbor = all_neighbors[i_neighbor].Element();
#ifdef PZDEBUG
                if (!gel_neighbor) {
                    DebugStop();
                }
#endif
                
                if (gel_neighbor->HasSubElement() == 1) {
                    continue;
                }
                
                if (gel_neighbor->Dimension() == m_geometry->Dimension()-1) {
                    int material_id = gel_neighbor->MaterialId();
                    bool is_boundary_memeber_Q = m_boundaries_material_ids.find(material_id) != m_boundaries_material_ids.end();
                    if (is_boundary_memeber_Q) {
                        d_minus_bc_neighbors_counter++;
                    }else{
                        if (material_id == m_fracture_id) {
                            d_minus_same_fracture_neighbors_counter++;
                        }else{
                            d_minus_other_fracture_neighbors_counter++;
                        }
                        
                    }
                }
                
                if (gel_neighbor->Dimension() == m_geometry->Dimension()-2) {
                    
//                    TPZStack<TPZGeoElSide> gel_sides;
//                    TPZGeoElSide gel_corner_side(gel_neighbor,gel_neighbor->NSides()-1);
//                    gel_corner_side.AllNeighbours(gel_sides);
//                    unsigned int fracture_counter = 0;
//                    for (auto gel_side: gel_sides) {
//
//                        if (gel_side.Element()->HasSubElement() == 1) {
//                            continue;
//                        }
//                        if (gel_side.Element()->Dimension() != m_geometry->Dimension() - 1) {
//                            continue;
//                        }
//
//                        if(gel_side.Element()->MaterialId() == m_fracture_id){
//                            fracture_counter++;
//                        }
//                    }
//
//                    if (fracture_counter == 1) {
//                        has_point_neighbor_Q = true;
//                    }

                    has_point_neighbor_Q = true;
                }
                
            }
            
            bool is_boundary_neighbour_Q = (d_minus_bc_neighbors_counter != 0 && has_point_neighbor_Q);
            bool is_the_end_of_fractrue_Q = d_minus_same_fracture_neighbors_counter == 0 && has_point_neighbor_Q;
            bool is_internal_side_Q = !is_the_end_of_fractrue_Q;
            
            bool is_pivot_Q = (is_boundary_neighbour_Q && is_the_end_of_fractrue_Q) || is_internal_side_Q;
            
            if (is_pivot_Q) {
                bool is_already_inserted_Q = HasInsertedPivot(gel_corner_side);
                if (!is_already_inserted_Q) {
                    m_pivot_indexes.push_back(gel_corner_side);
                }
            }else{
                bool is_already_inserted_Q = HasInsertedNonPivot(gel_corner_side);
                if (!is_already_inserted_Q) {
                    m_non_pivot_indexes.push_back(gel_corner_side);
                }
            }
        }
    }
    
    if(m_pivot_indexes.size() == 0){
        std::cout << "No pivots were detected." << std::endl;
        DebugStop();
    }
    
}

/// Verify if the non pivot has been inserted on non pivot structure
bool TPZFractureNeighborData::HasInsertedNonPivot(TPZGeoElSide pivot_side){
    
    TPZGeoNode node_candidate = pivot_side.Element()->Node(pivot_side.Side());
    
    bool has_inserted_non_pivot_Q = false;
    for (auto non_pivot: m_non_pivot_indexes) {
        TPZGeoNode node_target = non_pivot.Element()->Node(non_pivot.Side());
        if (node_candidate.Id() == node_target.Id()) {
            has_inserted_non_pivot_Q = true;
        }
    }

    return has_inserted_non_pivot_Q;
    
}

/// Verify if the pivot has been inserted on pivot structure
bool TPZFractureNeighborData::HasInsertedPivot(TPZGeoElSide pivot_side){
    
    TPZGeoNode node_candidate = pivot_side.Element()->Node(pivot_side.Side());
    
    bool has_inserted_pivot_Q = false;
    for (auto pivot: m_pivot_indexes) {
        TPZGeoNode node_target = pivot.Element()->Node(pivot.Side());
        if (node_candidate.Id() == node_target.Id()) {
            has_inserted_pivot_Q = true;
        }
    }
    
    return has_inserted_pivot_Q;
    
}

std::set<int64_t> TPZFractureNeighborData::PivotNeighbours(TPZGeoElSide pivotside){
    
    TPZStack<TPZGeoElSide> all_neighbors;
    pivotside.AllNeighbours(all_neighbors);
    std::set<int64_t> neigh_indexes;
    for (int i = 0; i < all_neighbors.size(); i++) {
        if (all_neighbors[i].HasSubElement() == 1) {
            continue;
        }
        if(all_neighbors[i].Element()->Dimension() == m_geometry->Dimension()){
            neigh_indexes.insert(all_neighbors[i].Element()->Index());
        }
        
    }
    return neigh_indexes;
}

/// Find the 2 elements that will initiate the element sets
// m_gel_left_indexes and m_gel_right_indexes have to be empty
// take a neighbouring fracture element - insert its two neighbours along the face
void TPZFractureNeighborData::InsertFractureNeighbours(std::set<int64_t> pivot_neighbours){
    
    
    int64_t fractture_neighbour_index = m_fracture_indexes[0];
    TPZGeoEl * fracture_neighbour = m_geometry->Element(fractture_neighbour_index);
    
    /// identify left and right seeds
    TPZStack<TPZGeoElSide> all_neigh;
    TPZGeoElSide fracture_side(fracture_neighbour,fracture_neighbour->NSides()-1);
    fracture_side.AllNeighbours(all_neigh);
    if (all_neigh.size() != 2) {
        DebugStop();
    }
    
    TPZGeoEl * gel_left = all_neigh[0].Element();
    TPZGeoEl * gel_right = all_neigh[1].Element();
    m_gel_left_indexes.insert(gel_left->Index());
    m_gel_right_indexes.insert(gel_right->Index());
    
}


void TPZFractureNeighborData::ClassifyElements(std::set<int64_t> pivotneighbours)
{
    while(pivotneighbours.size())
    {
        std::set<int64_t> found;
        for (auto elindex : pivotneighbours) {
            bool has_left_index_neighbor_Q = HasLeftIndexNeighbour(elindex);
            bool has_right_index_neighbor_Q = HasRightIndexNeighbour(elindex);
            if(has_left_index_neighbor_Q) m_gel_left_indexes.insert(elindex);
            if(has_right_index_neighbor_Q) m_gel_right_indexes.insert(elindex);
            if (has_left_index_neighbor_Q || has_right_index_neighbor_Q) {
                found.insert(elindex);
            }
        }
        for (auto elindex: found) {
            pivotneighbours.erase(elindex);
        }
    }
}


/// look if there is a fracture element neighbour (ok)
int64_t TPZFractureNeighborData::FractureNeighbourIndex(TPZGeoElSide gelside){

#ifdef PZDEBUG
    {
        int meshdim = m_geometry->Dimension();
        if (gelside.Dimension() != meshdim -1) {
            DebugStop();
        }
    }
#endif
    
    if (gelside.Element()->HasSubElement() == 1) {
        return -1;
    }
    
    TPZGeoElSide neighbour = gelside.Neighbour();
    while(neighbour != gelside)
    {
        if(neighbour.Element()->MaterialId() == m_fracture_id) return neighbour.Element()->Index();
        neighbour = neighbour.Neighbour();
    }
    return -1;
    
}

bool TPZFractureNeighborData::HasLeftIndexNeighbour(int64_t gel_index){
    
    TPZGeoEl * gel = m_geometry->Element(gel_index);
    unsigned int n_corner_sides = gel->NCornerNodes();
    unsigned int gel_candidate_side = gel->NSides() - 1;
    
    if (gel->HasSubElement() == 1) {
        return false;
    }
    
    if (m_gel_left_indexes.find(gel_index) != m_gel_left_indexes.end()) {
        return true;
    }
    
    for(int side=n_corner_sides; side<gel_candidate_side; side++)
    {
        TPZGeoElSide gelside(gel,side);
        if (FractureNeighbourIndex(gelside) != -1) {
            continue;
        }
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            int64_t neighindex = neighbour.Element()->Index();
            if (m_gel_left_indexes.find(neighindex) != m_gel_left_indexes.end()) {
                return true;
            }
            neighbour = neighbour.Neighbour();
        }
    }
    return false;
    
}

bool TPZFractureNeighborData::HasRightIndexNeighbour(int64_t gel_index){
    
    TPZGeoEl * gel = m_geometry->Element(gel_index);
    unsigned int n_corner_sides = gel->NCornerNodes();
    unsigned int gel_candidate_side = gel->NSides() - 1;
    
    if (gel->HasSubElement() == 1) {
        return false;
    }
    
    if (m_gel_right_indexes.find(gel_index) != m_gel_right_indexes.end()) {
        return true;
    }
    
    for(int side=n_corner_sides; side<gel_candidate_side; side++)
    {
        TPZGeoElSide gelside(gel,side);
        if (FractureNeighbourIndex(gelside) != -1) {
            continue;
        }
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            int64_t neighindex = neighbour.Element()->Index();
            if (m_gel_right_indexes.find(neighindex) != m_gel_right_indexes.end()) {
                return true;
            }
            neighbour = neighbour.Neighbour();
        }
    }
    return false;
    
}

bool TPZFractureNeighborData::HasClassifiedNeighbour(TPZGeoElSide & pivot_side){
    
    std::set<int64_t> neighbors = PivotNeighbours(pivot_side);
    
    bool has_classified_Q = false;
    bool has_left_index_Q;
    bool has_right_index_Q;
    for (auto index: neighbors) {
        if(HasLeftIndexNeighbour(index)){
            has_left_index_Q  = true;
        };
        if(HasRightIndexNeighbour(index)){
            has_right_index_Q  = true;
        };
    }
    if (has_left_index_Q && has_right_index_Q) {
        has_classified_Q = true;
    }
    
    return has_classified_Q;
    
}


/// Classify the neighbouring elements of the pivots
void TPZFractureNeighborData::ClassifyNeighboursofPivots(){

    std::ofstream filegvtk_b("Geometry_labels.vtk"); //Impressão da malha geométrica (formato vtk)
    TPZVTKGeoMesh::PrintGMeshVTK(m_geometry, filegvtk_b,true);
    
    BuildFractureElements();
    BuildPivotDataStructure();
    
    // Computes the seeds
    TPZGeoElSide pivot = m_pivot_indexes[0];
    std::set<int64_t> neighbours =  PivotNeighbours(pivot);
    InsertFractureNeighbours(neighbours);
    
    // Auxiliary list for account pivots that are processed
    std::set<TPZGeoElSide> pivots;
    for (auto gel_side: m_pivot_indexes) {
        pivots.insert(gel_side);
    }
    
    while (pivots.size()) {
        std::set<TPZGeoElSide> to_be_deleted;
        for (auto pivot: pivots) {
            bool has_classified_neighbor_Q = HasClassifiedNeighbour(pivot);
            if (has_classified_neighbor_Q) {
                std::set<int64_t> neighbours =  PivotNeighbours(pivot);
                ClassifyElements(neighbours);
                to_be_deleted.insert(pivot);
            }
        }
        for (auto pivot: to_be_deleted) {
            pivots.erase(pivot);
        }
    }
    
    
    
    int left_id = 100;
    for (auto it=m_gel_left_indexes.begin(); it != m_gel_left_indexes.end(); ++it){
        int64_t iel = *it;
        m_geometry->Element(iel)->SetMaterialId(left_id);
    }
    
    int right_id = 200;
    for (auto it=m_gel_right_indexes.begin(); it != m_gel_right_indexes.end(); ++it){
        int64_t iel = *it;
        m_geometry->Element(iel)->SetMaterialId(right_id);
    }
    
    std::ofstream filegvtk("Geometry_labels.vtk"); //Impressão da malha geométrica (formato vtk)
    TPZVTKGeoMesh::PrintGMeshVTK(m_geometry, filegvtk,true);
    
}

