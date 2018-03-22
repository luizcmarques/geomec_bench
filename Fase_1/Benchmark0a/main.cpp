
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

#include "TPZDarcyPMaterial.h"
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
const int matID = 1; //Materia do elemento volumétrico
//const int matFrac = 11; //Materia do elemento 1D, fratura
const int matBCbott = 2, matBCtop = 3, matBCright = 4, matBCleft = 5, matFrac = 6; //Materiais das condições de contorno
const int matInterface = 17; //Material do elemento de interface
const int matIntBCbott = matBCbott+10, matIntBCtop=matBCtop+10,  matIntBCright=matBCright+10, matIntBCleft=matBCleft+10, matIntFrac=matFrac+10; //Materiais das condições de contorno (elementos de interface)
//const int matPoint =-5;//Materia de um ponto
int dirichlet = 0, neumann = 1, penetration = 2, pointtype=5, dirichletPress=6; //Condições de contorno do problema ->default Dirichlet na esquerda e na direita
const REAL visco=1., perm=1., theta=-1.; //Coeficientes: viscosidade, fator simetria

void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget);

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
    
    HDivPiola = 0;
    
    int pOrder = 1; //Ordem polinomial de aproximação
    
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
    TPZCompMesh *cmesh_m = CMesh_m(gmesh, pOrder); //Função para criar a malha computacional multifísica
    
#ifdef PZDEBUG
    {
        std::ofstream filecv("MalhaC_v.txt"); //Impressão da malha computacional da velocidade (formato txt)
        std::ofstream filecp("MalhaC_p.txt"); //Impressão da malha computacional da pressão (formato txt)
        cmesh_v->Print(filecv);
        cmesh_p->Print(filecp);
    }
#endif
    
    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_v;
    meshvector[1] = cmesh_p;
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh_m);
    cmesh_m->LoadReferences();
    
#ifndef DarcyCase
    
    AddMultiphysicsInterfaces(*cmesh_m,matInterface,matID);
    AddMultiphysicsInterfaces(*cmesh_m,matIntBCbott,matBCbott);
    AddMultiphysicsInterfaces(*cmesh_m,matIntBCtop,matBCtop);
    AddMultiphysicsInterfaces(*cmesh_m,matIntBCright,matBCright);
    AddMultiphysicsInterfaces(*cmesh_m,matIntBCleft,matBCleft);
    AddMultiphysicsInterfaces(*cmesh_m,matIntFrac,matFrac);
    
#endif
    
#ifdef PZDEBUG
    std::ofstream fileg1("MalhaGeo2.txt"); //Impressão da malha geométrica (formato txt)
    gmesh->Print(fileg1);
    
    std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
    cmesh_m->Print(filecm);
#endif
    
    //Resolvendo o Sistema:
    int numthreads = 0;
    
    bool optimizeBandwidth = true; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)
    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    TPZSkylineNSymStructMatrix matskl(cmesh_m); //caso nao simetrico ***
    matskl.SetNumThreads(numthreads);
    
    std::set<int> matids;
    matids.insert(matID);
    matids.insert(matBCbott);
    matids.insert(matBCright);
    matids.insert(matBCtop);
    matids.insert(matBCleft);
    matids.insert(matIntFrac);
    matskl.SetMaterialIds(matids);
    
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an.SetSolver(step);
    an.Assemble(); //Assembla a matriz de rigidez (e o vetor de carga) global
    
    
#ifdef PZDEBUG
    //Imprimir Matriz de rigidez Global:
    {
        std::ofstream filestiff("stiffness.txt");
        an.Solver().Matrix()->Print("K1 = ",filestiff,EMathematicaInput);
        
        std::ofstream filerhs("rhs.txt");
        an.Rhs().Print("R = ",filerhs,EMathematicaInput);
        
        std::ofstream fileAlpha("alpha.txt");
        an.Solution().Print("Alpha = ",fileAlpha,EMathematicaInput);
    }
#endif
    
    an.Solve();
    
#ifdef PZDEBUG
    //Imprimindo vetor solução:
    {
        TPZFMatrix<STATE> solucao=cmesh_m->Solution();//Pegando o vetor de solução, alphaj
        std::ofstream solout("sol.txt");
        solucao.Print("Sol",solout,EMathematicaInput);//Imprime na formatação do Mathematica
    }
#endif
    
    //Calculo do erro
    
    //    TPZManVector<REAL,3> Errors;
    //    ofstream ErroOut("Erro.txt");
    //    an.SetExact(sol_exact1);
    //    an.PostProcessError(Errors,ErroOut);
    
    //Pós-processamento (paraview):
    
    std::string plotfile("DarcyCurve.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("P");
    vecnames.Push("V");
    scalnames.Push("f");
    vecnames.Push("V_exact");
    scalnames.Push("P_exact");
    
    int postProcessResolution = 3; //  keep low as possible
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
    
    
    // Criando e inserindo elemento de interface:
    
    long nel = gmesh->NElements();
    for (long el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        
        if (gel->Dimension() == 1) {
            
            
            int nsides = gel->NSides();
            
            TPZGeoElSide gelside(gel,nsides-1);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                
                if (neighbour.Element()->Dimension() == gmesh->Dimension() - 1) {
                    
                    break;
                    
                }
                neighbour = neighbour.Neighbour();
                
            }
            
            int mat_id = gel->MaterialId() + 10; // tagging interfaces materials on boundaries
            if (neighbour == gelside) {
                TPZGeoElBC(gelside, mat_id);
            }
            
            
        }
        
        if (gel->Dimension() != gmesh->Dimension()) {
            
            continue;
            
        }
        
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++) {
            if (gel->SideDimension(is) != gmesh->Dimension() - 1) {
                continue;
            }
            
            TPZGeoElSide gelside(gel,is);
            
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                
                if (neighbour.Element()->Dimension() == gmesh->Dimension() - 1) {
                    
                    break;
                    
                }
                neighbour = neighbour.Neighbour();
                
            }
            
            if (neighbour == gelside) {
                TPZGeoElBC(gelside, matInterface);
            }
        }
    }
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    int n_div = 0;
    UniformRefine(gmesh,n_div);
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
    TPZMat2dLin *material2 = new TPZMat2dLin(matFrac);
    
    cmesh->InsertMaterialObject(material); //Insere material na malha
    cmesh->InsertMaterialObject(material2); //Insere material na malha

    
    //cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1:
    //cmesh->ApproxSpace().CreateDisconnectedElements(true); //Criando elementos desconectados (descontínuo)
    
    cmesh->SetAllCreateFunctionsHDiv(); //Criando funções HDIV:
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    material->SetMaterial(xkin, xcin, xfin);
    material2->SetMaterial(xkin, xcin, xfin);
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(3,1,0.);
    
    TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, dirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno superior
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    TPZMaterial * BCond3 = material->CreateBC(material, matBCright, dirichlet, val1, val2);//Cria material que implementa a condicao de contorno direita
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    TPZMaterial * BCond2 = material->CreateBC(material, matBCleft, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
  //  TPZMaterial * BCond4 = material2->CreateBC(material2, matFrac, dirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
  //  cmesh->InsertMaterialObject(BCond4); //Insere material na malha
    
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
    cmesh->AutoBuild(matids);
    
    gmesh->ResetReference();
    matids.clear();
    matids.insert(matFrac);

    cmesh->AutoBuild(matids);
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
    
}

TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder)
{
    
    // @omar::
    
    //pOrder--; // Space restriction apapapa
    
    //Criando malha computacional:
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    
    // @omar::
    //cmesh->SetAllCreateFunctionsDiscontinuous();
    
    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    
    //Criando material:
    //Criando material cujo nSTATE = 2 ou seja linear
    
    TPZMat2dLin *material = new TPZMat2dLin(matID);//criando material que implementa a formulacao fraca do problema modelo
    
    cmesh->InsertMaterialObject(material); //Insere material na malha
    
    //Dimensões do material (para H1 e descontínuo):
    TPZFMatrix<STATE> xkin(1,1,0.), xcin(1,1,0.), xfin(1,1,0.);
    material->SetMaterial(xkin, xcin, xfin);
    
    //Condições de contorno:
    
    
    TPZFMatrix<REAL> val1(1,1,0.), val2(2,1,0.);
    
    
    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha
    
    int ncel = cmesh->NElements();
    for(int i =0; i<ncel; i++){
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if(!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *>(compEl);
        if(facel)DebugStop();
        
    }
    
    
    std::set<int> materialids;
    materialids.insert(matID);
    cmesh->AutoBuild(materialids);
    cmesh->LoadReferences();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild();
    
    
    // @omar::
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
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
    
    
    TPZDarcyPMaterial *material = new TPZDarcyPMaterial(matID,dim,1,visco,perm,theta);//criando material que implementa a formulacao fraca do problema modelo
    
    // Inserindo material na malha
    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (sol_exact1);
    
    material->SetForcingFunctionExact(solp);
    
    cmesh->InsertMaterialObject(material);
    
    
    //Condições de contorno:
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(3,1,0.);
    
    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0
    
    TPZMaterial * BCond0 = material->CreateBC(material, matBCbott, neumann , val1, val2); //Cria material que implementa a condição de contorno inferior
    //BCond0->SetForcingFunction(p_exact1, bc_inte_order);
    BCond0->SetForcingFunction(sol_exact1,bc_inte_order);
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha
    
    TPZMaterial * BCond1 = material->CreateBC(material, matBCtop, neumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    //BCond1->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond1->SetForcingFunction(solucao_exact2,bc_inte_order);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha
    
    
    TPZFMatrix<STATE> val2rg(3,1,0.);
    val2rg(0,0)=1.;
    
    TPZMaterial * BCond3 = material->CreateBC(material, matBCright, neumann, val1, val2); //Cria material que implementa a condicao de contorno direita
    //BCond3->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond3->SetForcingFunction(solucao_exact1,bc_inte_order);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha
    
    TPZFMatrix<STATE> val2lf(3,1,0.);
    val2lf(0,0)=1.;
    
    TPZMaterial * BCond2 = material->CreateBC(material, matBCleft, neumann, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    //BCond2->SetForcingFunction(p_exact1,bc_inte_order);
    //BCond2->SetForcingFunction(solucao_exact2,bc_inte_order);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha
    
    
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


void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget)
{
    
    TPZGeoMesh *gmesh = cmesh.Reference();
    std::set<int> velmatid;
    velmatid.insert(matFrac);

    
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int matid = gel->MaterialId();
        if(velmatid.find(matid) != velmatid.end())
        {
            int nsides = gel->NSides();
            TPZGeoElSide gelside(gel,nsides-1);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                
                TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
                LeftElIndices[0]=0;
                RightElIndices[0]=0;
                
                if (neighbour.Element()->Dimension() == 2 && neighbour.Element()->MaterialId() == matID && gelside.Element()->MaterialId() == matFrac) { //oioioi IDFlux -> ID
                    // create an interface element
                    TPZCompElSide celside = gelside.Reference();
                    TPZCompElSide celneigh = neighbour.Reference();
                    if (!celside || !celneigh) {
                        DebugStop();
                    }
                    std::cout << "Created an element between volumetric element " << neighbour.Element()->Index() <<
                    " side " << neighbour.Side() <<
                    " and interface element " << gelside.Element()->Index() << std::endl;
                    TPZGeoElBC gelbc(gelside,matID);
                    int64_t index;
                    TPZMultiphysicsInterfaceElement *intf = new
                    TPZMultiphysicsInterfaceElement(cmesh,gelbc.CreatedElement(),index,celneigh,celside);
                    intf->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
            
                }
                neighbour = neighbour.Neighbour();
            }
        }
        
    }
    
}


