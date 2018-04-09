
#include "pzgmesh.h"
#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"

#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"

#include "TPZDarcy1DMaterial.h"
#include "TPZDarcy2DMaterial.h"
#include "TPZLagrangeMultiplier.h"
#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"

//------------------Benchmarks for Geomec------------------------


/**
 * @brief Funcao para criar a malha geometrica do problema a ser simulado
 * @note A malha sera unidimensional formada por nel elementos de tamanho elsize
 * @param uNDiv number of divisions ortogonal to the plates performed on the domain
 * @param vNDiv number of divisions parallel to the plates performed on the domain
 * @param nel numero de elementos
 * @param elsize tamanho dos elementos
 */
TPZGeoMesh *CreateGMesh();

void UniformRefine(TPZGeoMesh* gmesh, int nDiv);

/**
 * @brief Funcao para criar a malha computacional da velocidade a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional da pressão a ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder);

/**
 * @brief Funcao para criar a malha computacional multi-fisica ser simulado
 * @note Responsavel pela criacao dos espacos de aproximacao do problema
 * @param gmesh malha geometrica
 * @param pOrder ordem polinomial de aproximacao
 */
TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder);


//Função para criar interface entre elmentos:

TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);


//Variáveis globais do problema:

const int dim = 2; //Dimensão do problema
const int dimFrac = 1; //Dimensão do problema
const int matID = 1; //Materia do elemento volumétrico
//const int matFrac = 11; //Materia do elemento 1D, fratura
const int matBCbott = 2, matBCtop = 3, matBCright = 4, matBCleft = 5, matFrac = 6, matPointLeft = 7, matPointRight = 8; //Materiais das condições de contorno
const int matInterface = 17; //Material do elemento de interface

const int matFluxWrap = 21;
//const int matPoint =-5;//Materia de um ponto
int dirichlet = 0, neumann = 1, penetration = 2, pointtype=5, dirichletPress=6; //Condições de contorno do problema ->default Dirichlet na esquerda e na direita
const REAL visco=1., perm=1., theta=-1.; //Coeficientes: viscosidade, fator simetria

void AddMultiphysicsInterfaces(TPZCompMesh &cmesh);
void BreakConnectivity(TPZCompMesh &cmesh, int matId);

using namespace std;

// definition of sol analytic
void sol_exact1(const TPZVec<REAL> & x, TPZVec<STATE>& f, TPZFMatrix<STATE> &dsol);

//Função principal do programa:

int main(int argc, char *argv[])
{
    
    TPZMaterial::gBigNumber = 1.e16;
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    //Dados do problema:
    
    HDivPiola = 1;
    
    int pOrder = 2; //Ordem polinomial de aproximação
    
    //Gerando malha geométrica:
    
    TPZGeoMesh *gmesh = CreateGMesh(); //Função para criar a malha geometrica
    
    
#ifdef PZDEBUG
    std::ofstream fileg("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
    std::ofstream filegvtk("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
    gmesh->Print(fileg);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
#endif
    
    //Gerando malha computacional:
    
    TPZCompMesh *cmesh_v = CMesh_v(gmesh, pOrder); //Função para criar a malha computacional da velocidade
    TPZCompMesh *cmesh_p = CMesh_p(gmesh, pOrder); //Função para criar a malha computacional da pressão
    
    
    {
        std::ofstream filecv("MalhaC_v.txt"); //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecp("MalhaC_p.txt"); //Impressão da malha computacional da pressão (formato txt)
        cmesh_v->Print(filecv);
        cmesh_p->Print(filecp);
    }
    
    BreakConnectivity(*cmesh_v, matFrac);

    {
        std::ofstream filecv("MalhaC_v.txt"); //Impressão da malha computacional da velocidade (formato txt)
        cmesh_v->Print(filecv);
    }
    
    TPZCompMesh *cmesh_m = CMesh_m(gmesh, pOrder); //Função para criar a malha computacional multifísica

    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_v;
    meshvector[1] = cmesh_p;
    
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    cmesh_m->LoadReferences();
    
    {
        std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
        cmesh_m->Print(filecm);
    }

    
    AddMultiphysicsInterfaces(*cmesh_m);
    
    {
        std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
        cmesh_m->Print(filecm);
    }
    


    std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
    gmesh->Print(fileg1);
    

    
    //Resolvendo o Sistema:
    int numthreads = 0;
    
    bool optimizeBandwidth = false; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    TPZSkylineNSymStructMatrix matskl(cmesh_m); //caso nao simetrico ***
    
    matskl.SetNumThreads(numthreads);

//    std::set<int> matids;
//    matids.insert(matID);
//    matids.insert(matFrac);
//    matids.insert(matBCbott);
//    matids.insert(matBCright);
//    matids.insert(matBCtop);
//    matids.insert(matBCleft);
//    matids.insert(matIntFrac);
//    matskl.SetMaterialIds(matids);
    
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an.SetSolver(step);
    an.Assemble(); //Assembla a matriz de rigidez (e o vetor de carga) global
    
    
//#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    {
        std::ofstream filestiff("stiffness.txt");
        an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
        
        std::ofstream filerhs("rhs.txt");
        an.Rhs().Print("R = ",filerhs,EMathematicaInput);
        
        std::ofstream fileAlpha("alpha.txt");
        an.Solution().Print("Alpha = ",fileAlpha,EMathematicaInput);
    }
//#endif
    
    an.Solve();
    

    //Imprimindo vetor solução:
    {
        TPZFMatrix<STATE> solucao=cmesh_m->Solution();//Pegando o vetor de solução, alphaj
        std::ofstream solout("sol.txt");
        solucao.Print("Sol",solout,EMathematicaInput);//Imprime na formatação do Mathematica
    }

    
    //Calculo do erro
    
    //    TPZManVector<REAL,3> Errors;
    //    ofstream ErroOut("Erro.txt");
    //    an.SetExact(sol_exact1);
    //    an.PostProcessError(Errors,ErroOut);
    
    //Pós-processamento (paraview):
    
    std::string plotfile("Benchmark_0a_DarcyTest.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("P");
    vecnames.Push("V");
    scalnames.Push("f");
    vecnames.Push("V_exact");
    scalnames.Push("P_exact");
    
    int postProcessResolution = 2; //  keep low as possible
    int dim = gmesh->Dimension();
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(postProcessResolution,dim);
    
    std::cout << "FINISHED!" << std::endl;
    
    return 0;
}


// definition of v analytic
void sol_exact1(const TPZVec<REAL> & x, TPZVec<STATE>& f, TPZFMatrix<STATE> &dsol){
    
    f.resize(3);
    
    STATE xv = x[0];
    STATE yv = x[1];
    STATE theta = atan(yv/xv);
    STATE r=sqrt(xv*xv+yv*yv);
    
    STATE v_x =  -r*sin(theta);
    STATE v_y =  r*cos(theta);
    STATE p =   0.;
    
    f[0] = v_x; // x direction
    f[1] = v_y; // y direction
    f[2] = p; //
}



TPZGeoMesh *CreateGMesh()
{
    
    int64_t id, index;
    
    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:
    
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    //Aqui é implementado um método para malhas criadas no GMSH
    
    //std::string dirname = PZSOURCEDIR;
    std::string grid;
    
    grid = "/Users/pablocarvalho/Documents/GitHub/geomec_bench/Fase_1/Benchmark0a/gmsh/msh/GeometryBench.msh";
    
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    gmesh = Geometry.GeometricGmshMesh(grid);
    

    //if (neighbour.Element()->Dimension() == 2 && neighbour.Element()->MaterialId() == matID && gelside.Element()->MaterialId() == matFrac)
    
    // Criando e inserindo elemento de interface:
    
//    long nel = gmesh->NElements();
//    for (long el = 0; el<nel; el++) {
//        TPZGeoEl *gel = gmesh->Element(el);
//
//        if (gel->Dimension() == 1) {
//
//
//            int nsides = gel->NSides();
//
//            TPZGeoElSide gelside(gel,nsides-1);
//            TPZGeoElSide neighbour = gelside.Neighbour();
//            while (neighbour != gelside) {
//
//                if (neighbour.Element()->Dimension() == gmesh->Dimension() - 1) {
//
//                    break;
//
//                }
//                neighbour = neighbour.Neighbour();
//
//            }
//
//            int mat_id = gel->MaterialId() + 10; // tagging interfaces materials on boundaries
//            if (neighbour == gelside) {
//                TPZGeoElBC(gelside, mat_id);
//            }
//
//
//        }
//
//        if (gel->Dimension() != gmesh->Dimension()) {
//
//            continue;
//
//        }
//
//        int nsides = gel->NSides();
//        for (int is = 0; is<nsides; is++) {
//            if (gel->SideDimension(is) != gmesh->Dimension() - 1) {
//                continue;
//            }
//
//            TPZGeoElSide gelside(gel,is);
//
//            TPZGeoElSide neighbour = gelside.Neighbour();
//            while (neighbour != gelside) {
//
//                if (neighbour.Element()->Dimension() == gmesh->Dimension() - 1) {
//
//                    break;
//
//                }
//                neighbour = neighbour.Neighbour();
//
//            }
//
//            if (neighbour == gelside) {
//                TPZGeoElBC(gelside, matInterface);
//            }
//        }
//    }
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();
    
//    int n_div = 0;
//    UniformRefine(gmesh,n_div);
    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;
    
}

void UniformRefine(TPZGeoMesh* gmesh, int nDiv)
{
    for(int D = 0; D < nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZVec< TPZGeoEl * > filhos;
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            gel->Divide(filhos);
        }
    }
    gmesh->BuildConnectivity();
}

TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);

    return NULL;
}


TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);//Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim);//Insere dimensão do modelo
    
    
    //Definição do espaço de aprximação:
    
    TPZMat2dLin *material = new TPZMat2dLin(matID); //Criando material que implementa a formulação fraca do problema modelo
    TPZMat1dLin *materialFrac = new TPZMat1dLin(matFrac);
    
    cmesh->InsertMaterialObject(material); //Insere material na malha
    cmesh->InsertMaterialObject(materialFrac); //Insere material na malha

    
    //cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1:
    //cmesh->ApproxSpace().CreateDisconnectedElements(true); //Criando elementos desconectados (descontínuo)
    
    cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
 
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xbin(1,1,0.), xfin(1,1,0.);
    material->SetMaterial(xkin, xcin, xfin);
    materialFrac->SetMaterial(xkin, xcin, xbin, xfin);
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(3,1,0.);
    
    TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    TPZMaterial * BCond2 = material->CreateBC(material, matBCright, dirichlet, val1, val2);//Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    TPZMaterial * BCond3 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    //Mat Frac:
    
    TPZMaterial * BCond4 = materialFrac->CreateBC(materialFrac, matPointRight , dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond4); //Insere material na malha
    
    TPZMaterial * BCond5 = materialFrac->CreateBC(materialFrac, matPointLeft , dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond5); //Insere material na malha
    
    //Criando material para FluxWrap
    
    TPZBndCond *FluxWrapBC;
    FluxWrapBC = material->CreateBC(material,matFluxWrap,dirichlet,val1,val2);
    cmesh->InsertMaterialObject(FluxWrapBC);

    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:

    
    
    int ncel = cmesh->NElements();
    
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    std::set<int> matids;
    matids.insert(matID);
    matids.insert(matBCbott);
    matids.insert(matBCright);
    matids.insert(matBCtop);
    matids.insert(matBCleft);
    matids.insert(matFluxWrap);
    
    cmesh->AutoBuild(matids);

    cmesh->SetDimModel(dimFrac);
    cmesh->SetAllCreateFunctionsHDiv();
    gmesh->ResetReference();
    matids.clear();
    matids.insert(matFrac);
    matids.insert(matPointLeft);
    matids.insert(matPointRight);

    cmesh->AutoBuild(matids);
//    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->ExpandSolution();
    
    return cmesh;
    
}

TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder)
{
    
    // @omar::
    
    //pOrder--; // Space restriction apapapa
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação

    //Definição do espaço de aprximação:
    
    TPZMat2dLin *material = new TPZMat2dLin(matID);//criando material que implementa a formulacao fraca do problema modelo
    TPZMat2dLin *materialFrac = new TPZMat2dLin(matFrac);//criando material que implementa a formulacao fraca do problema modelo
    cmesh->InsertMaterialObject(material); //Insere material na malha
    cmesh->InsertMaterialObject(materialFrac); //Insere material na malha
    
    cmesh->SetAllCreateFunctionsDiscontinuous(); //Criando funções H1
 //   cmesh->ApproxSpace().CreateDisconnectedElements(true);

    //Dimensões do material (para H1 e descontínuo):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    material->SetMaterial(xkin, xcin, xfin);
    materialFrac->SetMaterial(xkin, xcin, xfin);
 
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(3,1,0.);
    
//    TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
//    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
//
//    TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
//    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
//    TPZMaterial * BCond2 = material->CreateBC(material, matBCright, dirichlet, val1, val2);//Cria material que implementa a condicao de contorno direita
//    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
//
//    TPZMaterial * BCond3 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
//    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    //Mat Frac:
    
//    TPZMaterial * BCond4 = materialFrac->CreateBC(materialFrac, matPointRight, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
//    cmesh->InsertMaterialObject(BCond4); //Insere material na malha
//
//    TPZMaterial * BCond5 = materialFrac->CreateBC(materialFrac, matPointLeft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
//    cmesh->InsertMaterialObject(BCond5); //Insere material na malha
    
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    
//    int ncel = cmesh->NElements();
//    for(int i =0; i<ncel; i++){
//        TPZCompEl * compEl = cmesh->ElementVec()[i];
//        if(!compEl) continue;
//        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
//        if(facel)DebugStop();
//
//    }
    
    
    std::set<int> matids;
    matids.insert(matID);
    cmesh->SetDimModel(dim); //Insere dimensão do modelo

    cmesh->AutoBuild(matids);

    gmesh->ResetReference();
    matids.clear();
    matids.insert(matFrac);
    cmesh->SetDimModel(dimFrac); //Insere dimensão do modelo

    cmesh->AutoBuild(matids);

    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    
    return cmesh;
    
}

TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder)
{
    
    //Criando malha computacional:
    int bc_inte_order = 10;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    
    // Criando material:
    

    TPZDarcy2DMaterial *material = new TPZDarcy2DMaterial(matID,dim,1,theta);//criando material que implementa a formulacao fraca do problema modelo
    TPZFMatrix<REAL> K(dim,dim),invK(dim,dim);
    
    
    invK(0,0)=1./(3.38801);
    invK(1,1)=1./(2.566*1.e-4);
    K(0,0)=3.38801;
    K(1,1)=2.566*1.e-4;
//    invK(0,0)=1./(343.29);
//    invK(1,1)=1./(0.026);

    
    material->SetPermeabilityTensor(K, invK);
    
    TPZDarcy2DMaterial *materialFrac = new TPZDarcy2DMaterial(matFrac,dimFrac,1,theta);//criando material que implementa a formulacao fraca do problema modelo
    REAL kf = 4.68789*1.e3;
    REAL Dyf = 6.5*1.e-5;
    materialFrac->SetPermeability(kf*Dyf);
    
    // Inserindo material na malha
    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (sol_exact1);
    
    material->SetForcingFunctionExact(solp);
    
    cmesh->InsertMaterialObject(material);
    cmesh->InsertMaterialObject(materialFrac);
    
    
    //Condições de contorno:
    
    STATE Pjusante = 54.9;
    STATE Pmontante = 55.0;
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(3,1,0.);
    
    val1(0,0) = 0.; //botton
    
    TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, neumann , val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    //BCond0->SetForcingFunction(sol_exact1,bc_inte_order);
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    val1(0,0) = 0.; //top
    
    TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, neumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond1->SetForcingFunction(solucao_exact2,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    val1(0,0) = Pjusante; // right
    
    TPZMaterial * BCond2 = material->CreateBC(material, matBCright, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond3->SetForcingFunction(solucao_exact1,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    val1(0,0) = Pmontante; // left
    
    TPZMaterial * BCond3 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond2->SetForcingFunction(solucao_exact2,bc_inte_order);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    
    val1(0,0) =  Pjusante; // right
    
    TPZMaterial * BCond4 = materialFrac->CreateBC(materialFrac, matPointRight , dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond4); //Insere material na malha
    
    val1(0,0) = Pmontante; // left
    
    TPZMaterial * BCond5 = materialFrac->CreateBC(materialFrac, matPointLeft , dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond5); //Insere material na malha
    
    
    val1(0,0) = 0.0;
    
    TPZMaterial *MatLagrange = new TPZLagrangeMultiplier(matInterface,dimFrac,1); //
    cmesh->InsertMaterialObject(MatLagrange);
    
    TPZBndCond *FluxWrapBC = material->CreateBC(material,matFluxWrap,dirichlet,val1,val2);
    cmesh->InsertMaterialObject(FluxWrapBC);
    
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
    }
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

void BreakConnectivity (TPZCompMesh &cmesh, int matId)
{

    TPZGeoMesh *gmesh = cmesh.Reference();
    gmesh->ResetReference();
    cmesh.LoadReferences();
    cmesh.SetDimModel(dim);
    cmesh.SetAllCreateFunctionsHDiv();
    int64_t ncel = cmesh.NElements();
    
    
    for (int64_t el=0; el<ncel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if(gel->MaterialId()!=matId){
            continue;
        }
        if(gel->Dimension()!=dim-1){
            DebugStop();
        }
        
        TPZStack<TPZCompElSide> neigh;
        int nsides = gel->NSides();
        
        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gelside.Neighbour();
        
        gelside.EqualLevelCompElementList(neigh, 0, 0);
        
        if(neigh.size()!=2){
            DebugStop();
        }
        gel->ResetReference();
        neigh[0].Element()->Reference()->ResetReference();
        neigh[1].Element()->Reference()->ResetReference();
        
        //working on element 0
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(neigh[0].Element());
            
            if(!intel){
                DebugStop();
            }
            
            intel->LoadElementReference();
            
            int locindex = intel->MidSideConnectLocId(neigh[0].Side());
            TPZConnect &midsideconnect = intel->MidSideConnect(neigh[0].Side());
            if(midsideconnect.NElConnected() != 2)
            {
                DebugStop();
            }
            
            //Duplica um connect
            int64_t index = cmesh.AllocateNewConnect(midsideconnect.NShape(), midsideconnect.NState(), midsideconnect.Order());
            
            intel->SetConnectIndex(locindex, index);
            midsideconnect.DecrementElConnected();
            cmesh.ConnectVec()[index].IncrementElConnected();
            intel->SetSideOrient(neigh[0].Side(), 1);
            
            
            
            TPZGeoElBC bc(intel->Reference(),neigh[0].Side(),matFluxWrap);
            cmesh.CreateCompEl(bc.CreatedElement(), index);
            
            TPZCompEl *var = cmesh.Element(index);
            var->Reference()->ResetReference();
            intel->Reference()->ResetReference();
            
            
        }
        
        // working on element 1
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(neigh[1].Element());
            
            if(!intel){
                DebugStop();
            }
            
            intel->LoadElementReference();
            
            intel->SetSideOrient(neigh[1].Side(), 1);
            
            int64_t index;
            
            TPZGeoElBC bc(intel->Reference(),neigh[1].Side(),matFluxWrap);
            cmesh.CreateCompEl(bc.CreatedElement(), index);
            TPZCompEl *var = cmesh.Element(index);
            var->Reference()->ResetReference();

            intel->Reference()->ResetReference();
            
        }
        
    }

    cmesh.ExpandSolution();
}

void AddMultiphysicsInterfaces(TPZCompMesh &cmesh)
{
    
    TPZGeoMesh *gmesh = cmesh.Reference();
    std::set<int> velmatid;
    velmatid.insert(matFrac);

    
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int matid = gel->MaterialId();
        if(velmatid.find(matid) == velmatid.end())
        {
            continue;
        }
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            
            TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
            LeftElIndices[0]=0;
            RightElIndices[0]=1;
            
            if (neighbour.Element()->Dimension() == 1 && neighbour.Element()->MaterialId() == matFluxWrap) { //oioioi IDFlux -> ID
                // create an interface element
                TPZCompElSide celside = gelside.Reference();
                TPZCompElSide Wrapneigh = neighbour.Reference();
                if (!celside || !Wrapneigh) {
                    DebugStop();
                }
                std::cout << "Created an element between volumetric element " << neighbour.Element()->Index() <<
                " side " << neighbour.Side() <<
                " and interface element " << gelside.Element()->Index() << std::endl;
                TPZGeoElBC gelbc(gelside,matInterface);
                int64_t index;
                TPZMultiphysicsInterfaceElement *intf = new
                TPZMultiphysicsInterfaceElement(cmesh,gelbc.CreatedElement(),index,Wrapneigh,celside);
                intf->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
        
            }
            neighbour = neighbour.Neighbour();
        }
    }
    
}


