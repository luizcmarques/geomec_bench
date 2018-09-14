//
//  TPZDarcyAnalysis.h
//  Benchmark0a
//
//  Created by Pablo Carvalho on 14/09/18.
//

#ifndef TPZDarcyAnalysis_h
#define TPZDarcyAnalysis_h

#include <stdio.h>
#include "pzanalysis.h"
#include "TPZSimulationData.h"
#include "pzpostprocanalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"

class TPZDarcyAnalysis : public TPZAnalysis {
    
private:
    
    /// Pointer of Simulation data object
    TPZSimulationData * m_simulation_data;
    
    /// Solution at n+1 state
    TPZFMatrix<STATE> m_X_n;
    
    /// Solution at n (past) state
    TPZFMatrix<STATE> m_X;
    
    /// Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZManVector<TPZCompMesh * , 2> m_mesh_vec;
    
    /// Residue error
    STATE m_error;
    
    /// Correction variation
    STATE m_dx_norm;
    
    /// number of Newton iterations
    int m_k_iterations;
    
    /// Post-processor object
    TPZPostProcAnalysis * m_post_processor;
    
    /// Variables being postprocessed as scalar
    TPZStack<std::string> m_var_names;
    
    /// Variables being postprocessed as vector on standard post-process
    TPZStack<std::string> m_vec_var_names;
    
public:
    
    /// Default constructor
    TPZDarcyAnalysis();
    
    /// Destructor
    ~TPZDarcyAnalysis();
    
    /// Copy constructor
    TPZDarcyAnalysis(const TPZDarcyAnalysis & other);
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPZSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }
    
    /// Configurate the solver being used to compute the approximation
    void ConfigurateAnalysis(DecomposeType decomposition, TPZManVector<TPZCompMesh * , 2> & mesh_vec, TPZSimulationData * simulation_data, TPZVec<std::string> & var_names);
    
    /// Execute a single newton iteration
    void ExecuteNewtonInteration();
    
    /// Execute the evolution for a single time step
    void ExecuteOneTimeStep(bool must_accept_solution_Q = true);
    
    /// Post-processing the variables for a single time step from memory (is_stantdard_post_pro_Q = false)
    void PostProcessTimeStep(std::string & file, bool is_stantdard_post_pro_Q = true);
    
    /// Post-processing the variables for a single time step from DoF
    void StandardPostProcessTimeStep(std::string & file);
    
    /// Update the memory with the converged time step solution
    void AcceptTimeStepSolution();
    
    /// Load the current state for the hdiv and 2 meshes
    void LoadCurrentState();
    
    /// Load the last state for the hdiv and 2 meshes
    void LoadLastState();
    
    /** @brief Set Residue error */
    void Set_error(STATE error)
    {
        m_error = error;
    }
    
    /** @brief Get Residue error */
    STATE Get_error()
    {
        return m_error;
    }
    
    /** @brief Set Correction variation */
    void Set_dx_norm(STATE dxnorm)
    {
        m_dx_norm = dxnorm;
    }
    
    /** @brief Get Correction variation */
    STATE Get_dx_norm()
    {
        return m_dx_norm;
    }
    
    /** @brief Set number of Newton iterations */
    void Set_k_iterations(int kiterations)
    {
        m_k_iterations = kiterations;
    }
    
    /** @brief Get number of Newton iterations */
    int Get_k_iterations()
    {
        return m_k_iterations;
    }
    
    /** @brief Set variables being postprocessed as vector on standard post-process */
    void Set_vec_var_names(TPZStack<std::string> & vec_var_names)
    {
        m_vec_var_names = vec_var_names;
    }
    
    /** @brief Get variables being postprocessed as vector on standard post-process */
    TPZStack<std::string> & Get_vec_var_names()
    {
        return m_vec_var_names;
    }
    
};

#endif /* TPMRSMonoPhasicAnalysis_h */
