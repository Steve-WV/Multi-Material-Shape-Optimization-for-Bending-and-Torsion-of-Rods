#include <iostream>
#include <fstream>
using namespace std; 

#include<random>
#include "pde.hpp"
#include "mesh.hpp"
#include "vertex.hpp"
#include "tmv.hpp"

#include "boost/tokenizer.hpp"

#include "pde.hpp"



const double c0 = sqrt(2.0)/12.0;
//const double delta=0.25;
//const double eps = 0.003;
//const double tau = (1e-7);
//const double omega1 = 3.0;  // 1= minimising torsion; -1= maximising torsion, 0= no weigthing
//const double theta1 = 0.1; // regularization parameters for torsion

const double theta= 0.01; // regularization parameters for bending 
//const double omega2 = -3.0; // 1= minimising bending; -1= maximizing bending, 0= no weigthing (minimal bending rigidity D_min)
//const double omega3 = 0.0 ;  // 1= minimising bending; -1= maximizing bending, 0= no weigthing (maximal bending rigidity D_min)
//const double E = 71.0;   // Young's Modulus 
//const double G = 26.0;    // Shear Modulus G= E/2(1+poissonnumber)


// Double Well Potential W + Derivative DW

inline double W( double s ) {
	return 1.0/4.0 * sqr(s)*sqr(s-1);
}

inline double DW( double s ) {
	return 1.0/2.0*s*(1.0 - 3.0*s + 2.0*sqr(s));
}


// factor in the calculation of the prandtl's stress function:   integral (zeta(u)*nabla phi *nabla v) = 2 integral (v) for all v

inline double zeta(double s,double c)
{
 
return (1.0/((s*(1.0-c)+c))); 
} 

// derivative of zeta 

inline double Dzeta(double s,double c)
{ 
	
return (-(1.0-c)/(sqr(s*(1.0-c)+c))); 
} 

// set u to 0, if u< 0 
inline double eta(double s)
{
	if (s<0.0) {
			return 0.0;
		   }
	
	else return s;  
}

                      
                      
                      
                      
// initialize u
void ImplicitPDESystem::InitU() {
	auto vi = m->VertexBegin(), ve = m->VertexEnd();
	for (; vi!=ve; ++vi) {
		
		if (!(*vi)->Boundary()) {
		   double x = (*vi)->x();
		   double y = (*vi)->y();
		   double r = std::sqrt(sqr(x)+sqr(y));
		   
		   if (opt_case == 1) {(*vi)-> u() = sqr(0.5*(1-tanh((r-0.44618)/(sqrt(2)*sqr(eps)))));}
		   else {(*vi)-> u() = ((rand() % 100)*0.01);}
		   
	   				 } 
	   
	    else {
		   (*vi)->u() = 0.0;
		 }		
	}
}




// Phi is set to zero

void ImplicitPDESystem::InitPhi() {
	auto vi = m->VertexBegin(), ve = m->VertexEnd();
	for (; vi!=ve; ++vi) {
		
		double x = (*vi)->x();
		double y = (*vi)->y();
		(*vi)-> Phi() = 0.0;
	}
}

// Initialize the system: triangulation, boundary points etc.

void ImplicitPDESystem::Init(std::string config_file, bool vis_only=false) {
	
	std::cerr << "Setting boundary... ";
	auto bvi = m->BdryVertexBegin(), bve = m->BdryVertexEnd();
	for (; bvi!=bve; ++bvi) {
		(*bvi)->MarkDirBdry();
	}
	std::cerr << "done; ";
	
	// now we can index...
	std::cerr << "Indexing " << m->VertexSize() << " vertices, ";
	auto vi = m->VertexBegin(), ve = m->VertexEnd();
	auto vi_p = m-> VertexBegin(), ve_p = m-> VertexEnd();
	N_dof = 0;
	// count basis functions
	for (; vi!=ve; ++vi) if ( !(*vi)->DirBdry() ) ++N_dof;
	idx_mutex = std::vector<std::mutex>(N_dof);
		
	std::cerr << N_dof << " dof, " << m->VertexSize() - N_dof << " on Dirichlet boundary... ";
    
	int idx = 0, idx_bdry = 0;
	vi = m->VertexBegin();
	for (; vi!=ve; ++vi) {
		if ( !(*vi)->DirBdry() ) {
			(*vi)->SetIndex(idx);
			++idx;
		} else 	{
			(*vi)->SetIndex(idx_bdry+N_dof);
			++idx_bdry;
		}
	}

	std::cerr << "done; ";
		
	std::cerr << "Preparing " << m->TriSize() << " triangles... ";
	N_tri = 0;
	Mesh::TriIt ti = m->TriBegin(), te = m->TriEnd();
	for (; ti!=te; ++ti) {
		(*ti)->Init_area();
		(*ti)->Init_grad_s();
		++N_tri;
	}
		
		
	std::cerr << "Preparing " << numThreads << " threads... ";
	size_t elpth = ceil( (double)N_tri/numThreads);
	ti = m->TriBegin();
	for ( int i=0; i<numThreads; ++i ) {
		std::pair< Mesh::TriIt, Mesh::TriIt > bounds;
		bounds.first = ti;
		int j = 0;
		while (j<elpth && ti!=te) {++j; ++ti;}
		bounds.second = ti;
		tri_th.push_back(bounds);
	}
#ifdef EIGEN_HAS_OPENMP
	std::cerr << "Eigen is running in parallel... ";
	Eigen::setNbThreads(numThreads);
	Eigen::initParallel();
#endif 
	std::cerr << "ok; ";

	std::cerr << "zustandsvektor... ";
	

	U = Vec(N_dof+1);	
	Phi =Vec(N_dof+1);
        
	Phi_Db=Vec(m->VertexSize()-N_dof);
	U_Db=Vec(m->VertexSize()-N_dof);

	vi = m->VertexBegin();
	
	for (; vi!=ve; ++vi) {
		if (!(*vi)->DirBdry()) {
			(*vi)->Attach( &U[(*vi)->Index()]);
			(*vi)->Attach1( &Phi[(*vi)->Index()]);
		}
		else { 
			(*vi)->Attach( &U_Db[(*vi)->Index()-N_dof] );
			(*vi)->Attach1( &Phi_Db[(*vi)->Index()-N_dof]);
		}
	}
	
	std::cerr << "done." << std::endl;
	
	LoadOptions(config_file);
	//std::cerr << E << G << eps << theta1 << opt_case << omega1 << omega2 << omega3 << tau << delta << std::endl;
	
	InitU();

	PrepKM();
	ones = Vec::Ones(N_dof+1);

	CalcV();
	std::cerr << "Int u of initial condition = " << v << std::endl;
	
	InitPhi();
	PrepKM();
	CalcIntBasis();
        
	// saddle point matrix for the semi-implizit sheme is computed and factorized. 
	//just one factorization necessary
	S = M + tau*(delta/c0)*K + L;
	solver.compute(S);

	CalcE();
}


// calculation of the volume by VLoop and CalcV

void ImplicitPDESystem::VLoop(double* v, Mesh::TriIt ti, Mesh::TriIt te ) {
	double u[NumIntPts];
	for (; ti!=te; ++ti) {
		(*ti)->Calc_u( u );
		for (int k = 0; k<NumIntPts; ++k) {
			//*v += (u[k]*(1-theta1)+theta1) * GaussWeights[k]*(*ti)->Area();
			  *v += u[k]*GaussWeights[k]*(*ti)->Area();
		}
	}
}

void ImplicitPDESystem::CalcV() {
	std::vector<double> th_v;
	for ( int j = 0; j < numThreads; ++j ) {
		th_v.push_back( 0.0 );
	}
	// distribute elements on threads
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j )
		threads.push_back( std::thread( &ImplicitPDESystem::VLoop, this, &th_v[j], tri_th[j].first, tri_th[j].second)  );

	// join threads
	for (auto &thread : threads) thread.join();
	energy = 0.0;
	// add up all volume from all threads
	
	for (int j = 0; j < numThreads; ++j) v += th_v[j];
}






// stiffness matrix for penalized Poisson's euqtaion is computed, use S_P= int (zeta(u)*grad_sj*grad_si)
void ImplicitPDESystem::PoissonLoop( SpMat* th_s1, Vec* th_b1, Mesh::TriIt ti, Mesh::TriIt te ) {
	
	std::vector<T> vals_s1;
	*th_b1 = Vec::Zero(N_dof+1);
	double u[NumIntPts];
	
	Vec2 grad_sj, grad_si;

	for (; ti!=te; ++ti) {
		(*ti)-> Calc_u(u);
                  
		for (int bj=0; bj<3; ++bj) {
			
			Vertex* vert_j = (*ti)->v(static_cast<VertexName>(bj));
			int idx_j = vert_j->Index();
			
			(*ti)->Calc_grad_s(bj, grad_sj);
			
			for (int bi=0; bi<3; ++bi) {
				Vertex* vert_i = (*ti)->v(static_cast<VertexName>(bi));
				int idx_i = vert_i->Index();
					
				(*ti)->Calc_grad_s(bi, grad_si);
				
				double k = dot(grad_sj, grad_si) * (*ti)->Area();
				double s= 0.0;
				
				for (int i = 0; i< 3; ++i) {
					s += zeta(eta(u[i]),theta1)*k*GaussWeights[i];
				}
				
				if (!vert_j->DirBdry() && !vert_i->DirBdry()) {
					vals_s1.push_back( T(idx_i, idx_j, s) );
				} 
				
				else if (!vert_i->DirBdry()) {
					double c_j = vert_j->Phi();
					(*th_b1)[ idx_i ] -= c_j * s;
				}
			}
		}
	}
	(*th_s1).setFromTriplets( vals_s1.begin(), vals_s1.end() );		
}

// compute stiffness matrix for penalized poisson' problem
void ImplicitPDESystem::PrepPoisson(){

	std::vector<SpMat> th_s1;
	std::vector<Vec> th_b1;

	for (int j= 0; j< numThreads; ++j){      
		th_s1.push_back(SpMat(N_dof+1,N_dof+1));
		th_b1.push_back( Vec() );
	}
 
	std:: vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j )
		threads.push_back( std::thread( &ImplicitPDESystem::PoissonLoop, this, &th_s1[j],&th_b1[j],tri_th[j].first, tri_th[j].second));

	// join threads 
	for (auto &thread:threads)thread.join();
	//std::cerr<< "La gente esta muy loca"; 

	// using N_dof+1 x N_dof+1 matrices, so that the matrix in the Poisson problem has the same size as the matrix in the semi-implicit procedure.
	S1=SpMat(N_dof+1,N_dof+1);
	B1=Vec(N_dof+1);
	for (int j=0; j<numThreads; ++j){
		S1+=th_s1[j];
		B1+=th_b1[j];
	}
	S1.coeffRef(N_dof,N_dof) = 1.0;
}

// matrices for the semi-implicit procedure

void ImplicitPDESystem::KMLoop ( SpMat* th_k, SpMat* th_m, SpMat* th_l, Vec* th_b, Mesh::TriIt ti, Mesh::TriIt te ) {
	
	std::vector<T> vals_k, vals_m, vals_l;
	*th_b = Vec::Zero(N_dof+1);
    
	Vec2 grad_sj, grad_si;
	
	for (; ti!=te; ++ti) {
		for (int bj=0; bj<3; ++bj) {
			
			Vertex* vert_j = (*ti)->v(static_cast<VertexName>(bj));
			int idx_j = vert_j->Index();
			
			(*ti)->Calc_grad_s(bj, grad_sj);
			
			for (int bi=0; bi<3; ++bi) {
				
				Vertex* vert_i = (*ti)->v(static_cast<VertexName>(bi));
				int idx_i = vert_i->Index();
					
				(*ti)->Calc_grad_s(bi, grad_si);
				
				double k = dot(grad_sj, grad_si) * (*ti)->Area();
				double m = (bi == bj ? 1.0/6.0 : 1.0/12.0) * (*ti)->Area();
                                
					
				if (!vert_j->DirBdry() && !vert_i->DirBdry()) {
					vals_k.push_back( T(idx_i, idx_j, k) );
					vals_m.push_back( T(idx_i, idx_j, m) );
                                        
				} else if (!vert_i->DirBdry()) {
					double c_j = vert_j->u();
					(*th_b)[ idx_i ] -= c_j * k;
				}
			
			}
			// lumped matrix 
			if (!vert_j->DirBdry()) {
				double l= (1.0/3.0)* (*ti)->Area(); 
				vals_l.push_back( T(N_dof, idx_j, l) );
				vals_l.push_back( T(idx_j, N_dof, l) );
			}
			
		}
	}
	(*th_k).setFromTriplets( vals_k.begin(), vals_k.end() );	
	(*th_m).setFromTriplets( vals_m.begin(), vals_m.end() );
	(*th_l).setFromTriplets( vals_l.begin(), vals_l.end() );
}

// compute matrices for the semi-implicit procedure
void ImplicitPDESystem::PrepKM() {
	std::vector<SpMat> th_k, th_m, th_l;
	std::vector<Vec> th_b;

	for ( int j = 0; j < numThreads; ++j ) {
		th_k.push_back( SpMat(N_dof+1, N_dof+1) );
		th_l.push_back( SpMat(N_dof+1, N_dof+1) );
		th_m.push_back( SpMat(N_dof+1, N_dof+1) );
		th_b.push_back( Vec() );
	}

	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j ) 
		threads.push_back( std::thread( &ImplicitPDESystem::KMLoop, this, 
	&th_k[j],
	&th_m[j],
	&th_l[j],
	&th_b[j],
	tri_th[j].first, tri_th[j].second)  );

           
	// join threads
	for (auto &thread : threads) thread.join();
	
	K = SpMat(N_dof+1,N_dof+1);
	L = SpMat(N_dof+1,N_dof+1);
	M = SpMat(N_dof+1,N_dof+1);
	B = Vec::Zero(N_dof+1);
   
	for (int j = 0; j < numThreads; ++j) {
		K += th_k[j];
		L += th_l[j];
		M += th_m[j];
		B += th_b[j];
	}
}

// loop for integral(s_j*x) and Integral(s_j*y) for basis functions s_j
void ImplicitPDESystem::IntBasisLoop(Mesh::TriIt ti, Mesh:: TriIt te){ 

	double sj[NumIntPts];
	Vec2 p[NumIntPts];

	for (; ti!=te; ++ti){
   
		(*ti)-> Calc_Coord( p );
		
		for (int bj=0; bj<3; ++bj ){
			Vertex* vert_j= (*ti)-> v(static_cast<VertexName>(bj));
      
      
			(*ti)-> Calc_s(bj,sj); 

			double g0=0.0, g1=0.0; 
			for (int k=0; k< NumIntPts; ++k){
        
				double tempg0 = (p[k].x())*sj[k]; 
				double tempg1 = (p[k].y())*sj[k];
           
				g0+= tempg0*GaussWeights[k]*(*ti)-> Area();
				g1+= tempg1 * GaussWeights[k]*(*ti) -> Area();
			}
        
        
			if (!vert_j -> DirBdry()){ 
				int idx = vert_j -> Index(); 
				std:: lock_guard<std::mutex> lock( idx_mutex[idx] ) ;
				G0[idx]+=g0 ;
				G1[idx]+=g1;
	
			}
	
        
		}

	}	
        
}

// compute of integral(s_j*x) and integral(s_j*y)
void ImplicitPDESystem:: CalcIntBasis() { 

	G0 = Vec::Zero(N_dof+1);
	G1 = Vec::Zero(N_dof+1);
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j ) 
		threads.push_back( std::thread( &ImplicitPDESystem::IntBasisLoop, this, tri_th[j].first, tri_th[j].second)  );
 
	// join threads
	for (auto &thread : threads) thread.join();
}


// Force for the right hand side 
void ImplicitPDESystem::ForceLoop(Mesh::TriIt ti, Mesh::TriIt te) {
	
	double sj[NumIntPts];
	double u[NumIntPts];
	Vec2 p[NumIntPts];

	for (; ti!=te; ++ti ) {
		(*ti)->Calc_Coord( p );
		(*ti)->Calc_u( u );
		
		
		for (int bj=0; bj<3; ++bj) {
			Vertex* vert_j = (*ti)->v(static_cast<VertexName>(bj));
			
			(*ti)->Calc_s(bj, sj);
			double f=0.0, f0 = 0.0, f1 = 0.0,f2=0.0,f3=0.0;
			
			for (int k = 0; k<NumIntPts; ++k) {
						
				double temp0 = -(delta/(c0*eps)) * DW(u[k]) * sj[k];
					
				
				double temp1= (sqr((p[k]-phat).x()))*sj[k];
				double temp2= (sqr((p[k]-phat).y()))*sj[k];
				double temp3= (((p[k]-phat).x())*((p[k]-phat).y()))*sj[k];
				
				
				
				f0 += temp0 * GaussWeights[k]*(*ti)->Area();
				
				f1 += temp1 * GaussWeights[k]*(*ti)->Area();
				f2+= temp2*GaussWeights[k]*(*ti)->Area();
				f3 += temp3 *GaussWeights[k]*(*ti)->Area(); 
				 
				
			}
			
			
			
			
			
			if (!vert_j->DirBdry()) {
				int idx = vert_j->Index();
				std::lock_guard<std::mutex> lock( idx_mutex[idx] );
				
				F0[idx] += f0;
				F_bend1[idx] += f1;
				
				F_bend2[idx]+= f2;
				
				F_bend3[idx]+= f3;
				
				
			}

		}
		

	}


		
} 

// computing the force 
void ImplicitPDESystem::CalcF() {
	
	F0 = Vec::Zero(N_dof+1);
	F_bend1 = Vec::Zero(N_dof+1);
	F_bend2 = Vec::Zero(N_dof+1);
	F_bend3 = Vec::Zero(N_dof+1);
	
	
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j ) 
		threads.push_back( std::thread( &ImplicitPDESystem::ForceLoop, this, tri_th[j].first, tri_th[j].second)  );
 
	// join threads
	for (auto &thread : threads) thread.join();
}


void ImplicitPDESystem::ELoop(double* e, Mesh::TriIt ti, Mesh::TriIt te ) {
	double u[NumIntPts];
	Vec2 grad_u;
	for (; ti!=te; ++ti) {
		(*ti)->Calc_u( u );
		(*ti)->Calc_grad_u( grad_u );

		for (int k = 0; k<NumIntPts; ++k) {
			*e += delta/(eps) * W(u[k]) * GaussWeights[k]*(*ti)->Area();
		}
		*e += ( (delta*eps)/(2.0) * grad_u.normsqr()  ) * (*ti)->Area();
	}
}

void ImplicitPDESystem::CalcE() {
	std::vector<double> th_e;
	for ( int j = 0; j < numThreads; ++j ) {
		th_e.push_back( 0.0 );
	}
	// distribute elements on threads
	std::vector<std::thread> threads1;
	for ( int j = 0; j < numThreads; ++j )
		threads1.push_back( std::thread( &ImplicitPDESystem::ELoop,
	this, &th_e[j], tri_th[j].first, tri_th[j].second)  );

	// join threads
	for (auto &thread : threads1) thread.join();
	energy = 0.0;
	// add up all energies from all threads
	
	for (int j = 0; j < numThreads; ++j) energy += (1.0/c0)*th_e[j];
	
}



 
void ImplicitPDESystem::phatLoop(Vec2* p, double* u_av, Mesh::TriIt ti, Mesh::TriIt te ) {
	double u[NumIntPts];
	Vec2 coord[NumIntPts];

	for (; ti!=te; ++ti) {
		(*ti)->Calc_u( u );
		(*ti)->Calc_Coord( coord );
		
		for (int k = 0; k<NumIntPts; ++k) {
			
			*p += (u[k]*(1.0-theta1)+theta1)*coord[k] * GaussWeights[k]*(*ti)->Area();
			*u_av += (u[k]*(1.0-theta1)+theta1) * GaussWeights[k]*(*ti)->Area();
		}
	}
}



void ImplicitPDESystem::Calcphat() {
	std::vector<Vec2> th_p;
	std::vector<double> th_uav;
	
	for ( int j = 0; j < numThreads; ++j ) {
		th_p.push_back( Vec2(0.0,0.0) );
		th_uav.push_back( 0.0 );
	}
	// distribute elements on threads
	std::vector<std::thread> threads;
	
	for ( int j = 0; j < numThreads; ++j )
		threads.push_back( std::thread( &ImplicitPDESystem::phatLoop,
	this, &th_p[j], &th_uav[j], tri_th[j].first, tri_th[j].second)  );

	// join threads
	for (auto &thread : threads) thread.join();
	phat = Vec2(0.0, 0.0);
	uav = 0.0;
	
	// add up all energies from all threads
	for (int j = 0; j < numThreads; ++j) {
		phat += th_p[j];
		uav += th_uav[j];
	}
	phat /= uav;
}




void ImplicitPDESystem:: uxyLoop(Vec2* uxy,Mesh:: TriIt ti, Mesh:: TriIt te){


	double u[NumIntPts];
	Vec2 coord1[NumIntPts];
	for (; ti!=te; ++ti ) {
		(*ti) ->Calc_Coord(coord1);
		(*ti)-> Calc_u(u);
		for (int k = 0; k<NumIntPts; ++k) {	
			 
	*uxy+=2.0*(u[k]*(1.0-theta1)+theta1)*(coord1[k]-phat)*GaussWeights[k]*(*ti)->Area();
	                 
		}
                                
	}

}





void ImplicitPDESystem:: Calcuxy(){
	std::vector<Vec2> th_uxy;
	for ( int j = 0; j < numThreads; ++j ) {
		
		th_uxy.push_back(Vec2(0.0,0.0));
	}
	// distribute elements on threads
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j )
		threads.push_back( std::thread( &ImplicitPDESystem::uxyLoop,
	this,&th_uxy[j], tri_th[j].first, tri_th[j].second)  );

	// join threads
	for (auto &thread : threads) thread.join();
	uxyhat = Vec2(0.0,0.0);
	// add up all energies from all threads
	for (int j = 0; j < numThreads; ++j) {
		uxyhat += th_uxy[j];
	}
	uxyhat /= uav;
}



void ImplicitPDESystem::DxyLoop (double* Dx, double* Dy, double* Dxy, Mesh:: TriIt ti, Mesh:: TriIt te){

	double u[NumIntPts];
	Vec2 coord[NumIntPts];

	for (; ti!=te; ++ti) {
		(*ti)->Calc_u( u );
		(*ti)->Calc_Coord( coord );
		for (int k = 0; k<NumIntPts; ++k) {
			
			*Dx += (u[k]*(1.0-theta1)+theta1)*(sqr((coord[k]-phat).x())) * GaussWeights[k]*(*ti)->Area();
			*Dy += (u[k]*(1.0-theta1)+theta1)*(sqr((coord[k]-phat).y())) * GaussWeights[k]*(*ti)->Area();
			*Dxy += (u[k]*(1.0-theta1)+theta1)*((coord[k]-phat).x())*((coord[k]-phat).y()) * GaussWeights[k]*(*ti)->Area();
		   
		   
		         
		                                  }
		
		
                	}
		
		
		
    }



void ImplicitPDESystem::CalcDxy() {
	//std::vector<Vec2> th_p;
	std::vector<double> th_Dx,th_Dy,th_Dxy;
	
	for ( int j = 0; j < numThreads; ++j ) {
		//th_p.push_back( Vec2(0.0,0.0) );
		th_Dx.push_back( 0.0 );
		th_Dy.push_back( 0.0 );
		th_Dxy.push_back( 0.0 );
	}
	// distribute elements on threads
	std::vector<std::thread> threads;
	
	for ( int j = 0; j < numThreads; ++j )
		threads.push_back( std::thread( &ImplicitPDESystem::DxyLoop,
	this, &th_Dx[j], &th_Dy[j],&th_Dxy[j], tri_th[j].first, tri_th[j].second)  );

	// join threads
	for (auto &thread : threads) thread.join();
	Dxhat = 0.0;
	Dyhat = 0.0;
	Dxyhat= 0.0;
	// add up all energies from all threads
	for (int j = 0; j < numThreads; ++j) {
		Dxhat += th_Dx[j];
		Dyhat += th_Dy[j];
		Dxyhat+= th_Dxy[j];
	}
	
}



void ImplicitPDESystem::F_TorsionLoop(Mesh::TriIt ti, Mesh::TriIt te){

	double sj[NumIntPts];
	double u[NumIntPts];
	Vec2 grad_Phi;

	for (; ti!=te; ++ti ) {
	
		(*ti)->Calc_gradPhi(grad_Phi);
		(*ti)->Calc_u(u);
	
		
		for (int bj=0; bj<3; ++bj) {
			Vertex* vert_j = (*ti)->v(static_cast<VertexName>(bj));
			
			(*ti)->Calc_s(bj, sj);
			
			double f = 0.0;
			for (int k = 0; k<NumIntPts; ++k) {
				double temp = -grad_Phi.normsqr() * Dzeta(eta(u[k]),theta1) * sj[k];
				f += temp * GaussWeights[k]*(*ti)->Area();
			}
			if (!vert_j->DirBdry()) {
				int idx = vert_j->Index();
				std::lock_guard<std::mutex> lock( idx_mutex[idx] );
				F_Torsion[idx] += f;
			}

		}
	}
		
} 

void ImplicitPDESystem::CalcF_Torsion() {
	
	F_Torsion= Vec::Zero(N_dof+1); 
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j ) 
		threads.push_back( std::thread( &ImplicitPDESystem::F_TorsionLoop, this, tri_th[j].first, tri_th[j].second));
 
	// join threads
	for (auto &thread : threads) thread.join();
  
}



void ImplicitPDESystem::Step() {
	
	std::cerr << "Step." << std::endl;

	Calcphat();
	Calcuxy();
	CalcDxy();
	
	// calculate Prandtl's stress function by solving penalized poisson problem
	PrepPoisson();  
	S_Phi=(1.0/G)*S1;
	poissonsolver.compute(S_Phi);
	F_Phi= 2.0*M*ones;
	Phi=poissonsolver.solve(F_Phi);
	
	//poisson problem solved 
	CalcE(); 
	RM=sqrt( (sqr(Dxhat-Dyhat)/4)+sqr(Dxyhat));
	
        Vec test3=M*U; 

	// calculate right hand side
	CalcF();
	
	F_bend1-= G0*uxyhat.x();
	F_bend2-= G1* uxyhat.y();
	F_bend3-= 0.5*(G0*uxyhat.x()+G1*uxyhat.y());
	
        //variation of RM 
	delta_RM=((Dxhat-Dyhat)*(F_bend1-F_bend2)+4.0*Dxyhat*F_bend3)/(2.0*(sqrt(sqr(Dxhat-Dyhat)+4.0*sqr(Dxyhat)+sqr(theta)))); 
	
	CalcF_Torsion();
		
	// final right hand side
	F_bend1*=(1.0-theta1);
	F_bend2*=(1.0-theta1);
	F = (tau/(eps))*(F0-E*omega2*(0.5*F_bend1+0.5*F_bend2)-E*omega3*(0.5*F_bend1+0.5*F_bend2+delta_RM) - (1/G)*omega1*F_Torsion) + M*U + tau*B;
	
	F(N_dof) = v;

	Vec test2=M*ones; 
	
	energy_bending2=E*(RM-0.5*(Dxhat+Dyhat));
	energy_torsion2=2.0*Phi.dot(test2); 

	//calculate L^infty error
	Vec U_test = U;
	U = solver.solve(F);
	//CalcV();
	//std::cerr<< "Konvergenz? " << ((U_test-U).norm())/N_dof<<std::endl;
}
	
void ImplicitPDESystem::Solve() {
	
	int k = 0;
	Export("out" + std::to_string(k));
	while (k< 4*1e5) {
		Step();
		++k;
		
		if  ((k%5000 ==0)) { 
			 
	Export("out" + std::to_string(k));	
	std::fstream dataablage1; 
	dataablage1.open("daten1.dat",ios::app);
	dataablage1 << "Perimeter ="<< energy<< std::endl; 
	dataablage1 << "Energie Biegung ="<< (1/E)*energy_bending2<< std::endl; 
	dataablage1 << "Energie Torsion ="<< (1/G)*energy_torsion2 << std::endl; 
		
        
	                        }
		
		
                  }
	
}



// load options regarding shear and bulk modulus etc.
void ImplicitPDESystem::LoadOptions(std::string fname) {
	std::ifstream is;
	is.open(fname.c_str());

	
	std::string line;
	while( std::getline(is, line) ){
		
		boost::char_separator<char> sep(" ");
		boost::tokenizer< boost::char_separator<char> > tok(line, sep);
		boost::tokenizer< boost::char_separator<char> >::iterator ti = tok.begin(), te = tok.end();
		
		if (ti==te) continue; //empty line.
		
		if ( *ti == "eps" ) eps = std::stod(*(++ti));
		if ( *ti == "tau" ) tau = std::stod(*(++ti));
		if ( *ti == "delta" ) delta = std::stod(*(++ti));
		if ( *ti == "theta1" ) theta1 = std::stod(*(++ti));
		if ( *ti == "omega1" ) omega1 = std::stod(*(++ti));
		if ( *ti == "omega2" ) omega2 = std::stod(*(++ti));
		if ( *ti == "omega3" ) omega3 = std::stod(*(++ti));
		if ( *ti == "E" ) E = std::stod(*(++ti));
		if ( *ti == "G" ) G = std::stod(*(++ti));
		if ( *ti == "opt_case" ) opt_case = std::stoi(*(++ti));
	}
}



