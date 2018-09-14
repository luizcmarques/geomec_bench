//
//  TPZSimulationData.h
//  PZ
//
//  Created by Omar on 8/28/18.
//
//

#ifndef TPZSimulationData_h
#define TPZSimulationData_h

#include <stdio.h>
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzstack.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include "TPZGmshReader.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckgeom.h"


/** @brief Object conatining several kind of informations being used anytime and anywhere */
class TPZSimulationData
{
    
protected:
    
    /** @brief Spatial refinemenet level */
    int m_h_level;
    
    /** @brief Polynomial order for elasticity component */
    int m_elasticity_order;
    
    /** @brief Polynomial order for diffusion component */
    int m_darcy_order;
    
    /** @brief Physical dimension of the domain */
    int m_dimesion;
    
    /** @brief Number of iteration */
    int m_n_iterations;
    
    /** @brief Residue overal tolerance */
    REAL m_epsilon_res;
    
    /** @brief Correction overal tolerance */
    REAL m_epsilon_cor;
    
    /** @brief Number of thread */
    int m_n_threads;
    
    /** @brief Name for the Gmsh geometry file being used */
    std::string m_geometry_file;
    
    /** @brief Name for the vtk files being postprocessed */
    std::string m_vtk_file;
    
    /** @brief Number of vtk resolution during postprocessing */
    int m_vtk_resolution;
    
public:
    
    
    /** @brief default constructor */
    TPZSimulationData();
    
    /** @brief default constructor */
    TPZSimulationData(const TPZSimulationData & other);
    
    /** @brief default constructor */
    TPZSimulationData &operator=(const TPZSimulationData &other);
    
    /** @brief destructor */
    ~TPZSimulationData();
    
    /** @brief Print object attributes */
    void Print();
    
    /// Access methods
    
    /** @brief Set the spatial refinemenet level */
    void Set_h_level(int h_level){
        m_h_level = h_level;
    }
    
    /** @brief Get the spatial refinemenet level */
    int Get_h_level(){
        return m_h_level;
    }
    
    /** @brief Set the polynomial order for the elasticity approximation */
    void Set_elasticity_order(int elasticity_order){
        m_elasticity_order = elasticity_order;
    }
    
    /** @brief Get the polynomial order for the elasticity approximation  */
    int Get_elasticity_order(){
        return m_elasticity_order;
    }
    
    /** @brief Set the polynomial order for the Darcy's approximation */
    void Set_darcy_order(int darcy_order){
        m_darcy_order = darcy_order;
    }
    
    /** @brief Get the polynomial order for the Darcy's approximation  */
    int Get_darcy_order(){
        return m_darcy_order;
    }
    
    /** @brief Set the problem dimension */
    void Set_dimesion(int dimesion){
        m_dimesion = dimesion;
    }
    
    /** @brief Get the problem dimension  */
    int Get_dimesion(){
        return m_dimesion;
    }
    
    /** @brief Set Newton iterations */
    void Set_n_iteraions(int n_iterations){
        m_n_iterations = n_iterations;
    }
    
    /** @brief Get Newton iterations  */
    REAL Get_n_iteraions(){
        return m_n_iterations;
    }
    
    /** @brief Set residue tolerance */
    void Set_epsilon_res(REAL epsilon_res){
        m_epsilon_res = epsilon_res;
    }
    
    /** @brief Get residue tolerance  */
    REAL Get_epsilon_res(){
        return m_epsilon_res;
    }
    
    /** @brief Set correction tolerance */
    void Set_epsilon_cor(REAL epsilon_cor){
        m_epsilon_cor = epsilon_cor;
    }
    
    /** @brief Get correction tolerance  */
    REAL Get_epsilon_cor(){
        return m_epsilon_cor;
    }
    
    /** @brief Get number of threads being used  */
    int Get_n_threads(){
        return m_n_threads;
    }
    
    /** @brief Set number of threads being used */
    void Set_n_threads(int n_threads){
        m_n_threads = n_threads;
    }
    
    /** @brief Get name for the Gmsh geometry file being used  */
    std::string Get_geometry_file(){
        return m_geometry_file;
    }
    
    /** @brief Set name for the Gmsh geometry file being used */
    void Set_n_threads(std::string geometry_file){
        m_geometry_file = geometry_file;
    }
    
    /** @brief Get name for the vtk files being postprocessed  */
    std::string Get_vtk_file(){
        return m_vtk_file;
    }
    
    /** @brief Set name for the vtk files being postprocessed */
    void Set_vtk_file(std::string vtk_file){
        m_vtk_file = vtk_file;
    }
    
    /** @brief Get number of vtk resolution during postprocessing  */
    int Get_vtk_resolution(){
        return m_vtk_resolution;
    }
    
    /** @brief Set number of vtk resolution during postprocessing */
    void Set_vtk_resolution(int vtk_resolution){
        m_vtk_resolution = vtk_resolution;
    }

};

#endif /* TPZSimulationData_h */
