#include <iostream>
#include <array>
#include <cmath>

#include "cuddObj.hh"
#include "RungeKutta4.hh"
#include "SymbolicSet.hh"
#include "SymbolicModelUnbalanced.hh"
#include "TicToc.hh"
#include "dimensionize.hh"
#include "extras.hh"
#include "system.hh"

using namespace std;

// NCS: assuming delays of N=Nmax=Nmin=4
#define N 4

// NCS: state/input types
typedef array<double,ssDIM+isDIM> ncsICState_t;
typedef array<double,N*ssDIM+isDIM> ncsState_t;


// GLOBALS
Cudd cuddManager;
TicToc tt;

// NCS: This function recursively constructs the normal ncs-state given any initial x0 and u0
void constructNCSStates(size_t depth, sysState_t& x0, sysInput_t& u, vector<sysState_t>& currentPath, vector<ncsState_t>& pathsPool, sysInput_t& u_add){

	sysState_t x_current;
	(depth == 0)?x_current = x0:x_current=currentPath[depth-1];
	vector<sysState_t> post_states = sys_posts(x_current,u);


	for(size_t i=0; i<post_states.size(); i++){
		currentPath[depth]=post_states[i];
		if(depth == N-1){
			ncsState_t tmp;
			for(size_t j=0; j<N; j++)
				for(size_t k=0; k<ssDIM; k++)
					tmp[j*ssDIM+k] = currentPath[j][k];

			for(size_t j=0; j<isDIM; j++)
				tmp[N*ssDIM+j]=u_add[j];

			pathsPool.push_back(tmp);
		}
		else
			constructNCSStates(depth+1,x0,u,currentPath, pathsPool,u_add);
	}
}


// NCS: Post function for initial states
// (x0,u0) -- u --> (x1,....,xn,u1)
bool ncsIC_post(ncsICState_t& x, sysInput_t& u, vector<ncsState_t>& ManyxNcs){
  
  sysInput_t u0;
  sysState_t xs;

  // get x0
  for(size_t i=0; i<ssDIM; i++)
    xs[i]  = x[i];
  
  // get u0
  for(size_t i=0; i<isDIM; i++)
    u0[i]  = x[ssDIM+i];

  vector<sysState_t> tmp(N);
  constructNCSStates(0,xs,u0,tmp,ManyxNcs,u);

  return true;
}

// NCS: Normal post function
// (x11,....,x1n,u1) -- u --> (x21,....,x2n,u2)
bool ncs_post(ncsState_t& x, sysInput_t& u, vector<ncsState_t>& ManyxNcs){
  
  sysInput_t u1;
  sysState_t xs;

  // get x1n
  for(size_t i=0; i<ssDIM; i++)
    xs[i]  = x[ssDIM*(N-1)+i];
  
  // get u1
  for(size_t i=0; i<isDIM; i++)
    u1[i]  = x[ssDIM*N+i];

  vector<sysState_t> tmp(N);
  constructNCSStates(0,xs,u1,tmp,ManyxNcs,u);

  return true;
}

// NCS: This function constructs the space of the Ncs-Initial-State formed by (x0,u0) tuple
// where x0 \in X and u0 \in U
scots::SymbolicSet createNcsICStateSpace(Cudd &mgr) {

  double  lb[ssDIM+isDIM];  
  double  ub[ssDIM+isDIM]; 
  double eta[ssDIM+isDIM];
  
  size_t k=0;
  for(size_t i=0; i<ssDIM; i++, k++){
    lb[k]  = sysSSLb[i];
    ub[k]  = sysSSUb[i];
    eta[k] = sysSSeta[i];
  }

  for(size_t i=0; i<isDIM; i++, k++){
    lb[k]  = sysISLb[i];
    ub[k]  = sysISUb[i];
    eta[k] = sysISeta[i];
  }

  scots::SymbolicSet ss(mgr, ssDIM+isDIM, lb, ub, eta);
  return ss;
}

// NCS: This function constructs the space of the Ncs-State formed by (x1,....xn,u1) tuple
// where xi \in X and u1 \in U
scots::SymbolicSet createNcsStateSpace(Cudd &mgr) {

  double  lb[ssDIM*N+isDIM];  
  double  ub[ssDIM*N+isDIM]; 
  double eta[ssDIM*N+isDIM];
  
  size_t k=0;
  for(size_t i=0; i<ssDIM*N; i++, k++){
    lb[k]  = sysSSLb[i%ssDIM];
    ub[k]  = sysSSUb[i%ssDIM];
    eta[k] = sysSSeta[i%ssDIM];
  }

  for(size_t i=0; i<isDIM; i++, k++){
    lb[k]  = sysISLb[i];
    ub[k]  = sysISUb[i];
    eta[k] = sysISeta[i];
  }

  scots::SymbolicSet ss(mgr, ssDIM*N+isDIM, lb, ub, eta);
  return ss;
}

// NCS: This function constructs the normal input space U
scots::SymbolicSet createSysInputSpace(Cudd &mgr) {
  scots::SymbolicSet ss(mgr, isDIM, sysISLb, sysISUb, sysISeta);
  return ss;
}

int main(){
  // creating the ncs-intial-state space and filling all points in it
  scots::SymbolicSet ncsICSpace = createNcsICStateSpace(cuddManager);
  ncsICSpace.addGridPoints();
  PrintSymSetInfo(ncsICSpace, "ncsICSpace");

  // creating input space
  scots::SymbolicSet isSpace = createSysInputSpace(cuddManager);
  isSpace.addGridPoints();
  PrintSymSetInfo(isSpace, "isSpace");

  // creating the ncs-state/ ncs-post-state spaces
  scots::SymbolicSet ncsStateSpace = createNcsStateSpace(cuddManager);
  PrintSymSetInfo(ncsStateSpace, "ncsStateSpace");
  scots::SymbolicSet ncsPostStateSpace = createNcsStateSpace(cuddManager);
  PrintSymSetInfo(ncsPostStateSpace, "ncsPostStateSpace");

  // creating a transition relation of the initial transitions
  scots::SymbolicModelUnbalanced<ncsICState_t, sysInput_t, ncsState_t> ic_transitions(&ncsICSpace, &isSpace, &ncsStateSpace);
  cout << "Computing the initial_transitions ... ";
  tt.tic();
  ic_transitions.computeDirectedTransitionRelation(ncsIC_post);
  cout << "% ";
  tt.toc();

  // projecting the posts of the initial-ransitions as pres

  BDD ic_trans_bdd = ProjectBDD(cuddManager, ic_transitions.getTransitionRelation(), getSymSetBddVars(ncsStateSpace));
  ncsStateSpace.setSymbolicSet(ic_trans_bdd);

  // creating a transition relation of the full transitions
  scots::SymbolicModelUnbalanced<ncsState_t, sysInput_t, ncsState_t> full_transitions(&ncsStateSpace, &isSpace, &ncsPostStateSpace);
  cout << " Computing the normal_transitions ... ";
  tt.tic();
  full_transitions.computeDirectedTransitionRelation(ncs_post);
  cout << "% ";
  tt.toc();

  // saving bdds
  SaveBDD(ic_transitions.getTransitionRelation(), "ic_transitions.bdd");
  SaveBDD(full_transitions.getTransitionRelation(), "full_transitions.bdd");

  return 0;
}



