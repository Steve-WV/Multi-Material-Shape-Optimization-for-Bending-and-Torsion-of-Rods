#ifndef __PDE_HPP__
#define __PDE_HPP__

#include <string>
#include "mesh.hpp"

#include <thread>
#include <mutex>

#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>

class ImplicitPDESystem {
public:
	
  typedef Eigen::VectorXd Vec;
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMat;
  typedef Eigen::Triplet<double> T;
  typedef Eigen::SimplicialLDLT<SpMat> SpSolver; 
  typedef Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> CGSolver;
  typedef Eigen::SimplicialLLT<SpMat> SPDsolver;
    
  inline ImplicitPDESystem( Mesh* mesh ) : m(mesh) { }
  
  void Init(std::string, bool);
  void InitU();
  void Step();
  void Solve();
  void Export(std::string) const;
  void Export_Data(std::string) const; 
  void InitPhi(); 
  

private:
	void VLoop(double*, Mesh::TriIt, Mesh:: TriIt);
	void CalcV();
	void PrepKM();
	void KMLoop( SpMat*, SpMat*, SpMat*, Vec*, Mesh::TriIt, Mesh::TriIt);
	void PoissonLoop(SpMat*,Vec*,Mesh::TriIt,Mesh::TriIt);
	void PrepPoisson(); 
	void CalcF();
	void ForceLoop(Mesh::TriIt, Mesh::TriIt);
	void F_TorsionLoop(Mesh::TriIt, Mesh::TriIt);
	void CalcF_Torsion();
	void ELoop(double*, Mesh::TriIt,Mesh::TriIt);
	void CalcE();
	void phatLoop(Vec2*, double*, Mesh::TriIt ti, Mesh::TriIt te );
	void uxyLoop(Vec2*,Mesh::TriIt ti, Mesh::TriIt te);
	void DxyLoop(double*, double*, double*, Mesh::TriIt ti, Mesh:: TriIt te);
	void IntBasisLoop(Mesh:: TriIt, Mesh:: TriIt);
	void CalcIntBasis();
	void Calcphat();
	void Calcuxy();
	void CalcDxy();
	void LoadOptions(std::string);
	
  // Access to the basis functions, etc.
  Mesh* m;
    
  // basis functions (that are not constrained) and elements
  int N_dof, N_tri;
  
 
  Vec U, U_Db, Phi, Phi_Db,ones;
  Vec F, F0, F_bend1,F_bend2,F_bend3,delta_RM, B,F_Torsion,F_Phi,B1;
  Vec G0, G1,G2;
  SpMat K, L, M, S,S_Phi,S1;
  double energy, energy1,energy_torsion2,energy_torsion1,total_energy, total_energy1,energy_bending1,energy_bending2;
  double energy_phase_abnahme,energy_torsion_abnahme,energy_bending_abnahme;
  double RM,H1_0;
 double energy_bendingtestung,energy_torsiontesttung,energyphase_testung;    
 
 double eps,tau,delta,omega1,theta1,omega2,omega3;
 double E, G; 
 int opt_case;
 
  Vec2 phat,uxyhat;
  double uav,Dxhat, Dyhat,Dxyhat;
  double v;
    
  // linear solver
  SpSolver solver;
  CGSolver cg; 
  SPDsolver poissonsolver; 
  
  // mutexes and threading
  mutable std::vector<std::mutex> idx_mutex;
  mutable std::mutex mut;
  mutable std::mutex mut1;
  std::vector<std::pair<Mesh::TriIt, Mesh::TriIt> > tri_th;
  std :: vector<std::pair<Mesh::TriIt,Mesh::TriIt>> tri_tl;
 
  // general params
  const int numThreads = 4; // number of threads
};


#endif
