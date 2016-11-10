#ifndef EXTRAS_HH_
#define EXTRAS_HH_

#include <iostream>
#include <fstream>
#include <vector>
#include "cuddObj.hh"

using namespace std;

vector<size_t> getSymSetBddVars(scots::SymbolicSet& ss){
	vector<size_t> ret;
	for(size_t i=0; i<ss.getDimension(); i++)
		for(size_t j=0; j<ss.getNofBddVars()[i]; j++)
			ret.push_back(ss.getIndBddVars()[i][j]);

	return ret;
}

void PrintSymSetInfo(scots::SymbolicSet& ss, const char* name){
	cout << "Info for " << name << " : " << endl;
	cout << "    size: " << ss.getSize() << endl;
	cout << "     dim: " << ss.getDimension() << endl;
	cout << " bddVars: ";

	vector<size_t> bddVars = getSymSetBddVars(ss);
	for(size_t i=0; i<bddVars.size(); i++)
		cout << bddVars[i] << " ";

	cout << endl;
}

BDD ProjectBDD(Cudd& cuddManager, const BDD& srcBDD, vector<size_t> projVars){
	size_t nAllBddVars = cuddManager.ReadSize();
	vector<BDD> otherBDDs;

	for(size_t i=0; i<nAllBddVars; i++)
		if (!(find(projVars.begin(),projVars.end(),i) != projVars.end()))
			otherBDDs.push_back(cuddManager.bddVar(i));

	BDD otherBdDDsCube = cuddManager.bddComputeCube(otherBDDs.data(), NULL, otherBDDs.size());
	return srcBDD.ExistAbstract(otherBdDDsCube);
}

void SaveBDD(const BDD& srcBDD, const char* filename){
	FILE *file = fopen (filename,"w");
	if (file == NULL){
		std::ostringstream os;
		os << "Error: Unable to open file for writing." << filename << "'.";
		throw std::runtime_error(os.str().c_str());
	}

	/* before we save the BDD to file, we save it to another manager,
	 * because the current manager is classified as ADD manager */
	Cudd mdest;
	BDD tosave = srcBDD.Transfer(mdest);

	int storeReturnValue = Dddmp_cuddBddStore(
	  mdest.getManager(),
	  NULL,
	  tosave.getNode(),
	  //(char**)varnameschar, // char ** varnames, IN: array of variable names (or NULL)
	  NULL, // char ** varnames, IN: array of variable names (or NULL)
	  NULL,
	  DDDMP_MODE_BINARY,
	  // DDDMP_VARNAMES,
	  DDDMP_VARIDS,
	  NULL,
	  file
	);

	fclose(file);
	if (storeReturnValue!=DDDMP_SUCCESS)
	  throw "Error: Unable to write BDD to file.";
	else
	  cout << "BDD saved to file: "<< filename << std::endl;
}

#endif
