//
//  TPMRSSegregatedAnalysis.cpp
//  PMRS
//
//  Created by Omar Durán on 9/13/18.
//

#include "TPZSegregatedAnalysisDFN.h"


TPZSegregatedAnalysisDFN::TPZSegregatedAnalysisDFN(){
    m_simulation_data       = NULL;
    m_elastoplast_analysis  = NULL;
    m_darcy_analysis    = NULL;
}

TPZSegregatedAnalysisDFN::~TPZSegregatedAnalysisDFN(){
    
}

TPZSegregatedAnalysisDFN::TPZSegregatedAnalysisDFN(const TPZSegregatedAnalysisDFN & other){
    m_simulation_data       = other.m_simulation_data;
    m_elastoplast_analysis  = other.m_elastoplast_analysis;
    m_darcy_analysis    = other.m_darcy_analysis;
}

void TPZSegregatedAnalysisDFN::ApplyMemoryLink(){
    
    if (!m_simulation_data) {
        DebugStop();
    }
 
//    TPZManVector<std::pair<int, TPZManVector<int,12>>,12>  material_ids = m_simulation_data->MaterialIds();
//    TPZManVector<int,10> volumetric_mat_id(1);
//    int matid = material_ids[0].first;
//    volumetric_mat_id[0] = matid;
    
        int matid_elast = m_simulation_data->Get_elasticity_matid();
        int matid_darcy = m_simulation_data->Get_darcy_matid();

        TPZMaterial * material_geo = m_elastoplast_analysis->Mesh()->FindMaterial(matid_elast);
        TPZMaterial * material_res = m_darcy_analysis->Mesh()->FindMaterial(matid_darcy);
        if (!material_geo || !material_res) {
            DebugStop();
        }
    
        TPZMatWithMem<TPZMemoryDFN> * mat_with_memory_geo = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> * >(material_geo);
        TPZMatWithMem<TPZMemoryDFN> * mat_with_memory_res = dynamic_cast<TPZMatWithMem<TPZMemoryDFN> * >(material_res);
        if (!mat_with_memory_geo || !mat_with_memory_res) {
            DebugStop();
        }
        mat_with_memory_res->SetMemory(mat_with_memory_geo->GetMemory());
    
}

void TPZSegregatedAnalysisDFN::ConfigurateAnalysis(DecomposeType decompose_geo, DecomposeType decompose_res, TPZSimulationData * simulation_data,  TPZCompMesh * cmesh_geomechanics, TPZCompMesh * cmesh_reservoir, TPZManVector<TPZCompMesh * , 2> & mesh_vec, TPZStack<std::string> & post_pro_var_res, TPZStack<std::string> & post_pro_var_geo){
    
    if (!simulation_data || !cmesh_geomechanics || !cmesh_reservoir) {
        DebugStop();
    }
    
    this->SetSimulationData(simulation_data);
    bool mustOptimizeBandwidth = true;
    
    // The Geomechanics Simulator
    m_elastoplast_analysis = new TPZMatElastoPlasticAnalysis;
    m_elastoplast_analysis->SetCompMesh(cmesh_geomechanics,mustOptimizeBandwidth);
    m_elastoplast_analysis->ConfigurateAnalysis(decompose_geo, m_simulation_data);
    
    // The Flux Simulator
    m_darcy_analysis = new TPZDarcyAnalysis;
    m_darcy_analysis->SetCompMesh(cmesh_reservoir,mustOptimizeBandwidth);
    m_darcy_analysis->ConfigurateAnalysis(decompose_res, mesh_vec, m_simulation_data,post_pro_var_res);
    
    this->ApplyMemoryLink();
    
}

void TPZSegregatedAnalysisDFN::ExecuteOneTimeStep(bool must_accept_solution_Q){
    m_darcy_analysis->ExecuteOneTimeStep(must_accept_solution_Q);
    m_elastoplast_analysis->ExecuteOneTimeStep(must_accept_solution_Q);
}

void TPZSegregatedAnalysisDFN::PostProcessTimeStep(std::string & geo_file, std::string & res_file){
    m_darcy_analysis->PostProcessTimeStep(res_file);
    m_elastoplast_analysis->PostProcessTimeStep(geo_file);
}

void TPZSegregatedAnalysisDFN::ExecuteTimeEvolution(){
    
//    std::string file_res("ReservoirFlow.vtk");
//    std::string file_geo("Geomechanic.vtk");
//    
//    int n_max_fss_iterations = 10; // @TODO:: MS, please to xml file structure
//    int n_enforced_fss_iterations = 5; // @TODO:: MS, please to xml file structure
//    int n_time_steps = m_simulation_data->ReportingTimes().size();
//    REAL r_norm = m_simulation_data->epsilon_res();
//    REAL dx_norm = m_simulation_data->epsilon_cor();
//    
//    bool error_stop_criterion_Q = false;
//    bool dx_stop_criterion_Q = false;
//    for (int it = 0; it < n_time_steps; it++) {
//        for (int k = 1; k <= n_max_fss_iterations; k++) {
//            this->ExecuteOneTimeStep(false);
//            error_stop_criterion_Q = (m_darcy_analysis->Get_error() < r_norm) && (m_elastoplast_analysis->Get_error() < r_norm);
//            dx_stop_criterion_Q = (m_darcy_analysis->Get_dx_norm() < dx_norm) && (m_elastoplast_analysis->Get_dx_norm() < dx_norm);
//            
//            if ((error_stop_criterion_Q && (k > n_enforced_fss_iterations)) || dx_stop_criterion_Q) {
//                this->ExecuteOneTimeStep(true);
//                this->PostProcessTimeStep(file_geo, file_res);
//                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with residue norm for res = " << m_darcy_analysis->Get_error() << std::endl;
//                std::cout << "TPMRSSegregatedAnalysis:: Iterative process converged with residue norm for geo = " << m_elastoplast_analysis->Get_error() << std::endl;
//                break;
//            }
//        }
//    }
    
}
