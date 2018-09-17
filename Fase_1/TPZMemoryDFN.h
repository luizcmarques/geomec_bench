//
//  TPZMemoryDFN.cpp
//  Benchmark0a
//
//  Created by Pablo Carvalho on 17/09/18.
//


#ifndef TPZMemoryDFN_h
#define TPZMemoryDFN_h

#include <stdio.h>
#include "TPZMemoryDFN.h"
#include "TPZElastoPlasticMemoryDFN.h"
#include "TPZMonoPhasicMemoryDFN.h"

class TPZMemoryDFN : public TPZMonoPhasicMemoryDFN, public TPZElastoPlasticMemoryDFN {
    
private:
    
    /// Biot-Willis coefficient
    REAL m_alpha;
    
    /// Constrained specific storage at constant strain
    REAL m_Se;
    
public:
    
    /// Default constructor
    TPZMemoryDFN();
    
    /// Copy constructor
    TPZMemoryDFN(const TPZMemoryDFN & other);
    
    /// Assignement constructor
    const TPZMemoryDFN & operator=(const TPZMemoryDFN & other);
    
    /// Desconstructor
    virtual ~TPZMemoryDFN();
    
    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPZMemoryDFN & memory ){
        memory.Print(out);
        return out;
    }
    
    virtual int ClassId() const;
    
    void SetAlpha(REAL alpha){
        m_alpha  = alpha;
    }
    
    REAL GetAlpha(){
        return m_alpha;
    }
    
    void SetSe(REAL Se){
        m_Se  = Se;
    }
    
    REAL GetSe(){
        return m_Se;
    }
    
};

#endif /* TPZMemoryDFN_h */

