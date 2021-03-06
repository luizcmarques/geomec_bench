//
//  TPMRSGeomechanicAnalysis.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#ifndef TPZMatElastoPlasticAnalysis_h
#define TPZMatElastoPlasticAnalysis_h

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

class TPZMatElastoPlasticAnalysis : public TPZAnalysis {
    
private:
    
    /// Pointer of Simulation data object
    TPZSimulationData * m_simulation_data;
    
    /// Solution at n+1 state
    TPZFMatrix<STATE> m_X_n;
    
    /// Solution at n (past) state
    TPZFMatrix<STATE> m_X;
    
    /// Residue error
    STATE m_error;
    
    /// Correction variation
    STATE m_dx_norm;
    
    /// number of Newton iterations
    int m_k_iterations;
    
    /// Post-processor object
    TPZPostProcAnalysis * m_post_processor;
    
    /// Variables being postprocessed
    TPZStack<std::string> m_var_names;
    
public:
    
    /// Default constructor
    TPZMatElastoPlasticAnalysis();
    
    /// Destructor
    ~TPZMatElastoPlasticAnalysis();
    
    /// Copy constructor
    TPZMatElastoPlasticAnalysis(const TPZMatElastoPlasticAnalysis & other);
    
    /// Set the pointer of Simulation data object
    void SetSimulationData(TPZSimulationData * simulation_data){
        m_simulation_data = simulation_data;
    }
    
    /// Configurate the solver being used to compute the approximation
    void ConfigurateAnalysis(DecomposeType decomposition, TPZSimulationData * simulation_data);
    
    /// Execute a single newton iteration
    void ExecuteNewtonInteration();
    
    /// Execute the evolution for a single pseudo time step
    void ExecuteOneTimeStep(bool must_accept_solution_Q = true);
    
    /// Post-processing the variables for a single pseudo time step
    void PostProcessTimeStep(std::string & file);
    
    /// Update the memory with the converged pseudo time step solution
    void AcceptPseudoTimeStepSolution();
    
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
    
    
};

#endif /* TPZMatElastoPlasticAnalysis_h */
