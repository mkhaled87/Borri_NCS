#ifndef SYSTEM_HH_
#define SYSTEM_HH_

#include <array>
#include <vector>

using namespace std;

#define VEHICLE_SYSTEM

#ifdef DINT_SYSTEM

// SYS: system dimensions, types and bounds for input and state sets
#define ssDIM 2
#define isDIM 1
typedef array<double,ssDIM> sysState_t;
typedef array<double,isDIM> sysInput_t;
const double  sysSSLb[ssDIM] = {0.60,0.30};
const double  sysSSUb[ssDIM] = {3.40,3.30};
const double sysSSeta[ssDIM] = {0.20,0.30};
const double  sysISLb[isDIM] = {0.00};
const double  sysISUb[isDIM] = {1.00};
const double sysISeta[isDIM] = {0.30};

#define tau 0.3

/* SYS: system's ode RHS */
auto system_rhs =[](sysState_t& dxdt,  const sysState_t &x, sysInput_t &u) {
    const double b[2]    =  {0, 1};
    const double a[2][2] = {{0, 1},
							{0, 0}};

    dxdt[0] = a[0][0]*x[0]+a[0][1]*x[1] + b[0]*u[0];
    dxdt[1] = a[1][0]*x[0]+a[1][1]*x[1] + b[1]*u[0];
};

/* SYS: computation of the growth bound (the result is stored in r)  */
auto radius_post = [](sysState_t &r, sysInput_t &u) {
	  /* the ode to determine the radius of the cell which over-approximates the
	   * attainable set see: http://arxiv.org/abs/1503.03715v1 */
	  auto growth_bound_ode = [](sysState_t &drdt,  const sysState_t &r, const sysInput_t &u) {
		const double b[2]    =  {0, 1};
		const double a[2][2] = {{0, 1},
								{0, 0}};

		drdt[0] = a[0][0]*r[0]+a[0][1]*r[1] + b[0]*u[0];
		drdt[1] = a[1][0]*r[0]+a[1][1]*r[1] + b[1]*u[0];
	  };
	  const int nint=5;
	  OdeSolver ode_solver(ssDIM,nint,tau);
	  ode_solver(growth_bound_ode,r,u);
};

bool sys_constaints(sysState_t &x, sysInput_t &u){
	return true;
}

#endif

#ifdef VEHICLE_SYSTEM

// SYS: system dimensions, types and bounds for input and state sets
#define ssDIM 3
#define isDIM 2
typedef array<double,ssDIM> sysState_t;
typedef array<double,isDIM> sysInput_t;
const double  sysSSLb[ssDIM] = { 0, 0,-M_PI-0.4};
const double  sysSSUb[ssDIM] = { 6, 5, M_PI+0.4};
const double sysSSeta[ssDIM] = {.2,.2,.2};
const double  sysISLb[isDIM] = {-1, -1};
const double  sysISUb[isDIM] = { 1,  1};
const double sysISeta[isDIM] = {0.3,0.3};

#define tau 0.30

/* SYS: system's ode RHS */
auto system_rhs =[](sysState_t& dxdt,  const sysState_t &x, sysInput_t &u) {
      double alpha=atan(tan(u[1])/2.0);
      dxdt[0] = u[0]*cos(alpha+x[2])/cos(alpha);
      dxdt[1] = u[0]*sin(alpha+x[2])/cos(alpha);
      dxdt[2] = u[0]*tan(u[1]);
};

/* SYS: computation of the growth bound (the result is stored in r)  */
auto radius_post = [](sysState_t &r, sysInput_t &u) {
    double c = std::abs(u[0]*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1));
    r[0] = r[0]+c*r[2]*0.3;
    r[1] = r[1]+c*r[2]*0.3;
};

bool sys_constaints(sysState_t &x, sysInput_t &u){
	double xy_obsts[3][4] = {
			{1.2, 1.4, 2.4, 5.0},
			{2.8, 3.0, 0.0, 2.6},
			{4.4, 4.6, 2.4, 5.0}
	};
	for(size_t i=0; i<3; i++)
		if(x[0]>=xy_obsts[i][0] &&
		   x[0]<=xy_obsts[i][1] &&
		   x[1]>=xy_obsts[i][2] &&
		   x[2]<=xy_obsts[i][3])
			return false;

	return true;
}

#endif

// SYS: Post function that takes care of the growth of the states
vector<sysState_t> sys_posts(sysState_t& x, sysInput_t& u){

	const int nint=5;
	OdeSolver ode_solver(ssDIM,nint,tau);

	vector<sysState_t> Posts;
	if(!sys_constaints(x,u))
		return Posts;

	// The evolution of the state
	sysState_t xsource = x;
	ode_solver(system_rhs, xsource, u);

	// The growth radius
	sysState_t r;
    for(size_t i=0; i<ssDIM; i++)
      r[i]=sysSSeta[i]/2.0;
	radius_post(r, u);

	double xBoxLb[ssDIM];
	double xBoxUb[ssDIM];
	size_t xBoxPointsCount[ssDIM];
	vector<vector<double>> xBoxDataPerDim(ssDIM);
	for(size_t i=0; i<ssDIM; i++){
		xBoxLb[i] = xsource[i] - r[i];
		xBoxUb[i] = xsource[i] + r[i];
		xBoxPointsCount[i] = (xBoxUb[i]-xBoxLb[i])/sysSSeta[i];

		double data = xBoxLb[i] + r[i];
		if(xBoxPointsCount[i] <= 0)
			xBoxDataPerDim[i].push_back(data);
		else{
			for(size_t j=0; j<xBoxPointsCount[i]; j++){
				xBoxDataPerDim[i].push_back(data);
				data+=sysSSeta[i];
			}
		}
	}

	Dimensionize<double> DMZ(ssDIM, xBoxDataPerDim);
	vector<vector<double>> PostsVectors = DMZ.DoDimensionize();

	for(size_t i=0; i<PostsVectors.size(); i++){
		sysState_t tmp;
		for(size_t j=0; j<PostsVectors[i].size(); j++){
			tmp[j] = PostsVectors[i][j];
		}
		Posts.push_back(tmp);
	}

	return Posts;
}

#endif
