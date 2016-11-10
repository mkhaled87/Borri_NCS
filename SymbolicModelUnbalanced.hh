/*
 * SymbolicModel.hh
 *
 *  created on: 03.11.2016
 *  author: M. Khaled
 *  remarks: based on the original class SymbolicModel created on: 09.10.2015 by the author: rungger
 */

#ifndef SYMBOLICMODELUNBALANCED_HH_
#define SYMBOLICMODELUNBALANCED_HH_

#include <iostream>
#include <stdexcept>

#include "cuddObj.hh"
#include "dddmp.h"

#include "SymbolicSet.hh"


namespace scots {
/*
 * class: SymbolicModelUnbalanced
 * 
 * it stores the bdd of the transition relation and its bdd variable information
 * it provides an iterator to loop over all elements in the stateSpace and inputSpace
 *
 * it is different from the normal Symbolic Model as it supports post-states with different structures
 *
 */

template<class preStateType_, class inputType_, class postStateType_>
class SymbolicModelUnbalanced {
protected:
  /* var: ddmgr_ */
  Cudd *ddmgr_;

  /* var: stateSpace_ */
  SymbolicSet *stateSpace_;
  /* var: inputSpace_ */
  SymbolicSet *inputSpace_;
  /* var: stateSpacePost_ */
  SymbolicSet *stateSpacePost_;

  /* var: nssVars_ 
   * number of state space bdd variables (= number of pre variables)  */
  size_t nssVars_;
  /* var: preVars_ 
   * array of indices of the state space pre bdd variables  */
  size_t* preVars_;

  /* var: npostssVars_ 
   * number of post state space bdd variables (= number of post variables)  */
  size_t npostssVars_;
  /* var: postVars_ 
   * array of indices of the state space post bdd variables  */
  size_t* postVars_;

  /* var: nisVars_ 
   * number of input space bdd variables */
  size_t nisVars_;
  /* var: inpVars_ 
   * array of indices of the input space bdd variables  */
  size_t* inpVars_;


  /* var: iterator_ 
   * CuddMintermIterator to iterate over all elements in stateSpace_ and
   * inputSpace_ */
  CuddMintermIterator *iterator_;
  BDD it;

  /* var: transitionRelation_ 
   * the bdd representation of the transition relation  X x U x X
   * the bdd variables of the transition relation are 
   * given by preVars_ x inpVars_ x postVars_ */
  BDD transitionRelation_;
  friend class FixedPoint;

public:
  /* constructor: SymbolicModelUnbalanced 
   *
   * Representation of the transition relation as BDD in the space 
   *
   *   preX  x U x postX
   *  
   * provide SymbolicSet for preX 
   * provide SymbolicSet for U
   * provide SymbolicSet for postX 
   *
   * the SymbolicSet for preX and postX need to be identical, except the
   * BDD variable IDs need to differ
   *
   * 
   */
  SymbolicModelUnbalanced(SymbolicSet *stateSpace,  SymbolicSet *inputSpace,  SymbolicSet *stateSpacePost):
	  	  	    stateSpace_(stateSpace), inputSpace_(inputSpace), stateSpacePost_(stateSpacePost) {
    if(stateSpace_->ddmgr_!=inputSpace_->ddmgr_ || stateSpacePost_->ddmgr_!=inputSpace_->ddmgr_ ) {
      std::ostringstream os;
      os << "Error: scots::SymbolicModel: stateSpace and inputSpace need to have the same dd manager.";
      throw std::invalid_argument(os.str().c_str());
    }

    int differ = 0;
    /* check if stateSpace and stateSpacePost have different BDD IDs */
    for(size_t i=0; i<stateSpace_->dim_; i++) 
      for(size_t j=0; j<stateSpace_->nofBddVars_[i]; j++) 
        if(stateSpace_->indBddVars_[i][j]==stateSpacePost_->indBddVars_[i][j])
          differ=1;
    if(differ) {
      std::ostringstream os;
      os << "Error: scots::SymbolicModel: stateSpace and stateSpacePost are not allowed to have the same BDD IDs.";
      throw std::invalid_argument(os.str().c_str());
    }


    ddmgr_=stateSpace_->ddmgr_;

    nssVars_=0;
    for(size_t i=0; i<stateSpace_->dim_; i++) 
      for(size_t j=0; j<stateSpace_->nofBddVars_[i]; j++) 
       nssVars_++;

    nisVars_=0;
    for(size_t i=0; i<inputSpace_->dim_; i++) 
      for(size_t j=0; j<inputSpace_->nofBddVars_[i]; j++) 
       nisVars_++;

    npostssVars_ = 0;
    for(size_t i=0; i<stateSpacePost_->dim_; i++)
      for(size_t j=0; j<stateSpacePost_->nofBddVars_[i]; j++)
    	  npostssVars_++;

    /* set the preVars_ to the variable indices of stateSpace_ 
     * set the postVars_ to the variable indices of stateSpacePost_ */
    preVars_ = new size_t[nssVars_];
    postVars_ = new size_t[npostssVars_];

    for(size_t k=0, i=0; i<stateSpace_->dim_; i++) {
      for(size_t j=0; j<stateSpace_->nofBddVars_[i]; k++, j++) {
       preVars_[k]=stateSpace_->indBddVars_[i][j];
      }
    }

    for(size_t k=0, i=0; i<stateSpacePost_->dim_; i++) {
      for(size_t j=0; j<stateSpacePost_->nofBddVars_[i]; k++, j++) {
       postVars_[k]=stateSpacePost_->indBddVars_[i][j];
      }
    }


    inpVars_ = new size_t[nisVars_];
    for(size_t k=0, i=0; i<inputSpace_->dim_; i++) {
      for(size_t j=0; j<inputSpace_->nofBddVars_[i]; k++, j++) {
       inpVars_[k]=inputSpace_->indBddVars_[i][j];
      }
    }

    /* initialize the transition relation */
    transitionRelation_=stateSpace_->symbolicSet_*inputSpace_->symbolicSet_;
    it=transitionRelation_;
  };

  ~SymbolicModelUnbalanced(void) {
    delete[] preVars_;
    delete[] postVars_;
    delete[] inpVars_;
  };


  /* function:  getSize 
   * get the number of elements in the transition relation */
  inline double getSize(void) {
    return transitionRelation_.CountMinterm(nssVars_+nisVars_+npostssVars_);
  };


  /* function:  getTransitionRelation 
   * get the SymbolicSet which represents transition relation in X x U x X */
  inline BDD getTransitionRelation(void) const {
    return transitionRelation_;
  }; 


  /* function:  setTransitionRelation 
   * set the transitionRelation_ BDD to transitionRelation */
  void setTransitionRelation(BDD transitionRelation) {
    transitionRelation_=transitionRelation;
  };  /* iterator methods */


  /* function: begin
   * initilize the iterator and compute the first element */
  inline void begin(void) {
    size_t nvars_=nssVars_+nisVars_;
    std::vector<size_t> ivars_;
    ivars_.reserve(nvars_);
    for(size_t i=0; i<nssVars_; i++) 
      ivars_.push_back(preVars_[i]);
    for(size_t i=0; i<nisVars_; i++) 
      ivars_.push_back(inpVars_[i]);
    /* set up iterator */
    iterator_ = new CuddMintermIterator(it,ivars_,nvars_);
  }

  /* function: next
   * compute the next minterm */
  inline void next(void) {
    ++(*iterator_);
  }

  /* function: progress
   * print progess of iteratin in percent */
  inline void progress(void) const {
    iterator_->printProgress();
  }  /* function: done


   * changes to one if iterator reached the last element */
  inline int done(void) {
    if(iterator_->done()) {
      delete iterator_;
      iterator_=NULL;
      return 1;
    } else 
    return 0;
  }


  /* function: currentMinterm
   * returns the pointer to the current minterm in the iteration */
  inline const int* currentMinterm(void) const {
    if (iterator_)
      return iterator_->currentMinterm();
    else
      return NULL;
  }

  /* function:  computeTransitionRelation
     *
     * provide the solution of the system at sampling time and
     * provide the solution of the linear system associated with the growth bound
     * at sampling time
     *
     * see the example directory for the specific format
     *
     */
    template<class F1>
    void computeDirectedTransitionRelation(F1 &system_post) {

      /* create the BDD's with numbers 0,1,2,.., #gridPoints */
      //size_t     ssDim=stateSpace_->getDimension();
      size_t postssDim=stateSpacePost_->getDimension();

      const size_t* nvars= stateSpacePost_->getNofBddVars();
      BDD **bddVars = new BDD*[postssDim];
      for(size_t n=0, i=0; i<postssDim; i++) {
        bddVars[i]= new BDD[nvars[i]];
        for(size_t j=0; j<nvars[i]; j++)  {
          bddVars[i][j]=ddmgr_->bddVar(postVars_[n+j]);
        }
        n+=nvars[i];
      }

      const size_t* ngp= stateSpacePost_->getNofGridPoints();
      BDD **num = new BDD*[postssDim];
      for(size_t i=0; i<postssDim; i++) {
        num[i] = new BDD[ngp[i]];
        int *phase = new int[nvars[i]];
        for(size_t j=0;j<nvars[i];j++)
          phase[j]=0;
        for(size_t j=0;j<ngp[i];j++) {
          int *p=phase;
          int x=j;
          for (; x; x/=2) *(p++)=0+x%2;
          num[i][j]= ddmgr_->bddComputeCube(bddVars[i],(int*)phase,nvars[i]);
        }
        delete[] phase;
        delete[] bddVars[i];
      }
      delete[] bddVars;


      /* bdd nodes in pre and input variables */
      DdManager *mgr = ddmgr_->getManager();
      size_t ndom=nssVars_+nisVars_;
      int*  phase = new int[ndom];
      DdNode**  dvars = new DdNode*[ndom];
      for(size_t i=0;i<nssVars_; i++)
        dvars[i]=Cudd_bddIthVar(mgr,preVars_[i]);
      for(size_t i=0;i<nisVars_; i++)
        dvars[nssVars_+i]=Cudd_bddIthVar(mgr,inpVars_[i]);

      /* initialize cell radius
       * used to compute the growth bound */
      double eta[postssDim];
      stateSpacePost_->copyEta(&eta[0]);

      double first[postssDim];
      stateSpacePost_->copyFirstGridPoint(&first[0]);

      transitionRelation_=ddmgr_->bddZero();
      const int* minterm;

      /** big loop over all state elements and input elements **/
      for(begin(); !done(); next()) {
        progress();
        minterm=currentMinterm();

        /* current state */
        preStateType_ x;
        std::vector<postStateType_> ManyXdash;
        inputType_ u;

        stateSpace_->mintermToElement(minterm,&x[0]);

        /* current input */
        inputSpace_->mintermToElement(minterm,&u[0]);

        bool to_include = system_post(x, u, ManyXdash);
        if(!to_include)
      	  continue;

        if(ManyXdash.size() == 0)
        	continue;

        for(size_t xi=0; xi < ManyXdash.size(); xi++){

        	postStateType_ xdash = ManyXdash[xi];

			/* check if the post gets out of the bounds */
			double Lb[postssDim];
			double Ub[postssDim];
			stateSpacePost_->copyFirstGridPoint(&Lb[0]);
			stateSpacePost_->copyLastGridPoint(&Ub[0]);

			bool out_of_bounds = false;
			for(size_t i=0; i<postssDim; i++){
				if (xdash[i]>Ub[i] || xdash[i]<Lb[i])
					out_of_bounds = true;
			}
			if(out_of_bounds)
				continue;


			BDD post=ddmgr_->bddOne();
			for(size_t i=0; i<postssDim; i++) {
			  BDD zz=ddmgr_->bddZero();
			  int j = std::lround(((xdash[i]-first[i])/eta[i]));
			  zz|=num[i][j];
			  post &= zz;
			}

			if(!(post==ddmgr_->bddZero())) {
				/* compute bdd for the current x and u element and add x' */
				for(size_t i=0;i<nssVars_; i++)
				  phase[i]=minterm[preVars_[i]];

				for(size_t i=0;i<nisVars_; i++)
				  phase[nssVars_+i]=minterm[inpVars_[i]];

				BDD current(*ddmgr_,Cudd_bddComputeCube(mgr,dvars,phase,ndom));
				current&=post;

				transitionRelation_ +=current;
			}
        }
      }

      for(size_t i=0; i<postssDim; i++)
        delete[] num[i];

      delete[] num;
      delete[] dvars;
      delete[] phase;
    }


}; /* close class def */
} /* close namespace */

#endif /* SYMBOLICMODELUNBALANCED_HH_ */
