//
//  TPZDarcyAnalysis.cpp
//  Benchmark0a
//
//  Created by Pablo Carvalho on 14/09/18.
//

#include "TPZDarcyAnalysis.h"

TPZDarcyAnalysis::TPZDarcyAnalysis() : TPZAnalysis(){
    
    m_simulation_data = NULL;
    m_X_n.Resize(0, 0);
    m_X.Resize(0, 0);
    m_mesh_vec.Resize(0);
    m_error = 0;
    m_dx_norm = 0;
    m_k_iterations = 0;
    m_post_processor = NULL;
    m_var_names.resize(0);
    m_vec_var_names.resize(0);
    
}

TPZDarcyAnalysis::~TPZDarcyAnalysis(){
    
}

TPZDarcyAnalysis::TPZDarcyAnalysis(const TPZDarcyAnalysis & other){
    
    m_simulation_data   = other.m_simulation_data;
    m_X_n               = other.m_X_n;
    m_X                 = other.m_X;
    m_mesh_vec          = other.m_mesh_vec;
    m_error             = other.m_error;
    m_dx_norm           = other.m_dx_norm;
    m_k_iterations      = other.m_k_iterations;
    m_post_processor    = other.m_post_processor;
    m_var_names         = other.m_var_names;
    m_vec_var_names     = other.m_vec_var_names;
    
}

void TPZDarcyAnalysis::ConfigurateAnalysis(DecomposeType decomposition, TPZManVector<TPZCompMesh * , 2> & mesh_vec, TPZSimulationData * simulation_data, TPZVec<std::string> & var_names){
    
    SetSimulationData(simulation_data);
    TPZStepSolver<STATE> step;
    unsigned int n_threads = m_simulation_data->Get_n_threads();

    if(!Mesh()){
        std::cout << "Call SetCompMesh method." << std::endl;
        DebugStop();
    }
    
    m_mesh_vec = mesh_vec;
    switch (decomposition) {
        case ELU:
        {
            TPZSkylineNSymStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(n_threads);
            this->SetStructuralMatrix(struct_mat);
        }
            break;
        case ELDLt:
        {
            TPZSymetricSpStructMatrix struct_mat(Mesh());
            struct_mat.SetNumThreads(n_threads);
            this->SetStructuralMatrix(struct_mat);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    step.SetDirect(decomposition);
    this->SetSolver(step);
    this->Solution().Resize(Mesh()->Solution().Rows(), 1);
    m_X.Resize(Mesh()->Solution().Rows(), 1);
    m_X_n.Resize(Mesh()->Solution().Rows(), 1);

    m_post_processor = new TPZPostProcAnalysis;
    m_post_processor->SetCompMesh(Mesh());

    int n_vols = m_simulation_data->Get_volumetric_material_id().size();
    TPZManVector<int,10> post_mat_id(n_vols);
    for (int ivol = 0; ivol < n_vols; ivol++)
    {
        int matid = m_simulation_data->Get_volumetric_material_id()[ivol];
        post_mat_id[ivol] = matid;
    }

    for (auto i : var_names) {
        m_var_names.Push(i);
    }

    m_post_processor->SetPostProcessVariables(post_mat_id, m_var_names);
    TPZFStructMatrix structmatrix(m_post_processor->Mesh());
    structmatrix.SetNumThreads(n_threads);
    m_post_processor->SetStructuralMatrix(structmatrix);

}

void TPZDarcyAnalysis::ExecuteNewtonInteration(){
    this->Assemble();
    this->Rhs() *= -1.0;
    this->Solve();
}

void TPZDarcyAnalysis::ExecuteOneTimeStep(bool must_accept_solution_Q){
    
    LoadLastState();
    AcceptTimeStepSolution();

    LoadCurrentState();
    this->AcceptTimeStepSolution();

    TPZFMatrix<STATE> dx;
    bool residual_stop_criterion_Q = false;
    bool correction_stop_criterion_Q = false;
    REAL norm_res, norm_dx;
    REAL r_norm = m_simulation_data->Get_epsilon_res();
    REAL dx_norm = m_simulation_data->Get_epsilon_cor();
    int n_it = m_simulation_data->Get_n_iterations();

    for (int i = 1; i <= n_it; i++) {
        this->ExecuteNewtonInteration();
        dx = Solution();
        norm_dx  = Norm(dx);
        m_X_n += dx;

//        this->AcceptTimeStepSolution();
        LoadCurrentState();
        this->AssembleResidual();
        norm_res = Norm(Rhs());
        residual_stop_criterion_Q   = norm_res < r_norm;
        correction_stop_criterion_Q = norm_dx  < dx_norm;

        m_k_iterations = i;
        m_error = norm_res;
        m_dx_norm = norm_dx;

        if (residual_stop_criterion_Q ||  correction_stop_criterion_Q) {
#ifdef PZDEBUG
            std::cout << "TPZDarcyAnalysis:: Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "TPZDarcyAnalysis:: Number of iterations = " << i << std::endl;
            std::cout << "TPZDarcyAnalysis:: Correction norm = " << norm_dx << std::endl;
#endif
            if (must_accept_solution_Q) {
                m_X = m_X_n;
            }
            break;
        }
    }

    if (residual_stop_criterion_Q == false) {
        std::cout << "TPZDarcyAnalysis:: Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
}

void TPZDarcyAnalysis::PostProcessTimeStep(std::string & file, bool is_stantdard_post_pro_Q){
    
    if (is_stantdard_post_pro_Q) {
        this->StandardPostProcessTimeStep(file);
        return;
    }
    
    int dim = Mesh()->Dimension();
    int div = 0;
    TPZStack< std::string> vecnames;
    m_post_processor->TransferSolution();
    m_post_processor->DefineGraphMesh(dim,m_var_names,vecnames,file);
    m_post_processor->PostProcess(div,dim);
}

void TPZDarcyAnalysis::StandardPostProcessTimeStep(std::string & file){
    
    int postProcessResolution = m_simulation_data->Get_vtk_resolution();
    int dim = Mesh()->Reference()->Dimension();
    this->DefineGraphMesh(dim,m_var_names,m_vec_var_names,file);
    this->PostProcess(postProcessResolution,dim);
    
}

void TPZDarcyAnalysis::AcceptTimeStepSolution(){
    
//    bool state = m_simulation_data->IsCurrentStateQ();
//    if (state) {
//        m_simulation_data->Set_must_accept_solution_Q(true);
//        LoadCurrentState();
//        AssembleResidual();
//        m_simulation_data->Set_must_accept_solution_Q(false);
//    }else{
//        m_simulation_data->Set_must_accept_solution_Q(true);
//        LoadLastState();
//        AssembleResidual();
//        m_simulation_data->Set_must_accept_solution_Q(false);
//    }
}


void TPZDarcyAnalysis::LoadCurrentState(){
    LoadSolution(m_X_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_mesh_vec, Mesh());
}

void TPZDarcyAnalysis::LoadLastState(){
    LoadSolution(m_X);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_mesh_vec, Mesh());
}

