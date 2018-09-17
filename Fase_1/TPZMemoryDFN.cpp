//
//  TPZMemoryDFN.cpp
//  PMRS
//
//  Created by Omar and Manouchehr on 9/11/18.
//

#include "TPZMemoryDFN.h"


TPZMemoryDFN::TPZMemoryDFN() : TPZMonoPhasicMemoryDFN() , TPZElastoPlasticMemoryDFN() {
    m_alpha = 0.5;
    m_Se = 1.0*0.0000145038;
}

TPZMemoryDFN::TPZMemoryDFN(const TPZMemoryDFN & other): TPZMonoPhasicMemoryDFN(other), TPZElastoPlasticMemoryDFN(other) {
    m_alpha = other.m_alpha;
    m_Se    = other.m_Se;
}

const TPZMemoryDFN & TPZMemoryDFN::operator=(const TPZMemoryDFN & other) {
    
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_alpha = other.m_alpha;
    m_Se    = other.m_Se;
    
    return *this;
}

TPZMemoryDFN::~TPZMemoryDFN(){
    
}

const std::string TPZMemoryDFN::Name() const {
    return "TPZMemoryDFN";
}

void TPZMemoryDFN::Write(TPZStream &buf, int withclassid) const {
    TPZMonoPhasicMemoryDFN::Write(buf, withclassid);
    TPZElastoPlasticMemoryDFN::Write(buf, withclassid);
}

void TPZMemoryDFN::Read(TPZStream &buf, void *context){
    TPZMonoPhasicMemoryDFN::Read(buf, context);
    TPZElastoPlasticMemoryDFN::Read(buf, context);
}

void TPZMemoryDFN::Print(std::ostream &out) const {
    TPZMonoPhasicMemoryDFN::Print(out);
    TPZElastoPlasticMemoryDFN::Print(out);
}

int TPZMemoryDFN::ClassId() const {
    return Hash("TPZMemoryDFN");
}
