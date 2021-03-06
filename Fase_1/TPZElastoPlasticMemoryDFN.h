//
//  TPZElastoPlasticMemoryDFN.h
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#ifndef TPZElastoPlasticMemoryDFN_h
#define TPZElastoPlasticMemoryDFN_h

#include <stdio.h>
#include <iostream>
#include <string>
#include "TPZSavable.h"
#include "Hash/TPZHash.h"
#include "TPZStream.h"
#include "TPZPersistenceManager.h"
#include "TPZTensor.h"
#include "TPZPlasticState.h"

class TPZElastoPlasticMemoryDFN{
    
private:
    
    /// Current stress state
    TPZTensor<REAL> m_sigma_n;
    
    /// Current plastic_strain state
    TPZPlasticState<REAL> m_plastic_strain_n;
    
    /// Current displacement field
    TPZManVector<REAL,3> m_u_n;
    
    /// Last stress state
    TPZTensor<REAL> m_sigma;
    
    /// Last plastic_strain state
    TPZPlasticState<REAL> m_plastic_strain;
    
    /// Last displacement field
    TPZManVector<REAL,3> m_u;
    
    /// Initial stress state
    TPZTensor<REAL> m_sigma_0;
    
    /// Initial plastic_strain state
    TPZPlasticState<REAL> m_plastic_strain_0;
    
    /// Initial displacement field
    TPZManVector<REAL,3> m_u_0;
    
    
public:
    
    /// Default constructor
    TPZElastoPlasticMemoryDFN();
    
    /// Copy constructor
    TPZElastoPlasticMemoryDFN(const TPZElastoPlasticMemoryDFN & other);
    
    /// Assignement constructor
    const TPZElastoPlasticMemoryDFN & operator=(const TPZElastoPlasticMemoryDFN & other);
    
    /// Desconstructor
    virtual ~TPZElastoPlasticMemoryDFN();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPZElastoPlasticMemoryDFN & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;
    
    
    /// Set the current stress state
    void SetSigma_n(TPZTensor<REAL> & sigma_n){
        m_sigma_n = sigma_n;
    }
    
    /// Get the current stress state
    TPZTensor<REAL> & GetSigma_n(){
        return m_sigma_n;
    }
    
    /// Set the current plastic_strain state
    void SetPlasticState_n(TPZPlasticState<REAL> & plastic_strain_n){
        m_plastic_strain_n = plastic_strain_n;
    }
    
    /// Get the current plastic_strain state
    TPZPlasticState<REAL> & GetPlasticState_n(){
        return m_plastic_strain_n;
    }
    
    /// Set the current displacement field
    void Setu_n(TPZManVector<REAL,3> & u_n){
        m_u_n = u_n;
    }
    
    /// Get the current displacement field
    TPZManVector<REAL,3> & Getu_n(){
        return m_u_n;
    }
    
    /// Set the last stress state
    void SetSigma(TPZTensor<REAL> & sigma){
        m_sigma = sigma;
    }
    
    /// Get the last stress state
    TPZTensor<REAL> & GetSigma(){
        return m_sigma;
    }
    
    /// Set the last plastic_strain state
    void SetPlasticState(TPZPlasticState<REAL> & plastic_strain){
        m_plastic_strain = plastic_strain;
    }
    
    /// Get the last plastic_strain state
    TPZPlasticState<REAL> & GetPlasticState(){
        return m_plastic_strain;
    }
    
    /// Set the last displacement field
    void Setu(TPZManVector<REAL,3> & u){
        m_u = u;
    }
    
    /// Get the last displacement field
    TPZManVector<REAL,3> & Getu(){
        return m_u;
    }
    
    /// Set the initial stress state
    void SetSigma_0(TPZTensor<REAL> & sigma_0){
        m_sigma_0 = sigma_0;
    }
    
    /// Get the initial stress state
    TPZTensor<REAL> & GetSigma_0(){
        return m_sigma_0;
    }
    
    /// Set the initial plastic_strain state
    void SetPlasticState_0(TPZPlasticState<REAL> & plastic_strain_0){
        m_plastic_strain_0 = plastic_strain_0;
    }
    
    /// Get the initial plastic_strain state
    TPZPlasticState<REAL> & GetPlasticState_0(){
        return m_plastic_strain_0;
    }
    
    /// Set the initial displacement field
    void Setu_0(TPZManVector<REAL,3> & u_0){
        m_u_0 = u_0;
    }
    
    /// Get the initial displacement field
    TPZManVector<REAL,3> & Getu_0(){
        return m_u_0;
    }
    
};

#endif /* TPZElastoPlasticMemoryDFN_h */
