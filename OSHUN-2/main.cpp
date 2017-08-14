///////////////////////////////////////////////////////////
//   Contributing authors :  Michail Tzoufras
//   Last Modified:          May 18th 2010
///////////////////////////////////////////////////////////
//   
//   Main program, including the main loop. 
//
///////////////////////////////////////////////////////////

// Standard libraries 
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <stdarg.h>
 

#include <map>
#include <string>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cstring>


// My libraries 
#include "lib-array.h"
#include "lib-algorithms.h"

// Misc Declerations
#include "state.h"
#include "formulary.h"
#include "setup.h"
#include "vlasov.h"
#include "export.h"



//**************************************************************
//  This Clock controls an iteration loop, given an initial
//  time tout_start*dt_out it evaluates the number of time-
//  steps it takes to get to (tout_start+1)*dt_out and can 
//  be incremented; 

    class Clock {
    public:
//      Constructor
        Clock(int tout_start, double dt_out, double CFL) {
            _hstep = 0;
            _t_start = double(tout_start)*dt_out;
            _numh  = size_t(static_cast<int>(dt_out/CFL))+1;
            _h     = dt_out/static_cast<double>(_numh);
        }

//      Clock readings
        double h()     const {return _h;}                    // Step size
        size_t numh()  const {return _numh;}                 // # steps 
        size_t tick()  const {return _hstep;}                // Current step
        double time()  const {return tick()*h() + _t_start;} // Current time
             
//      Increment time
        Clock& operator++() { ++_hstep; return *this;}  

    private:
        double _t_start;   // Initial output timestep
        size_t _numh;      // # of steps derived from this CFL
        double _h;         // resulting time-step from CFL
        size_t _hstep;
    };
//--------------------------------------------------------------
//**************************************************************

void periodic_boundaries(State1D& Y, size_t Nbc);
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
int main(int argc, char** argv) {

//  INPUT PARAMETERS
    vector< string > oTags;
    oTags.push_back("Time");
    oTags.push_back("Space");
    oTags.push_back("px-x_0");
    oTags.push_back("px-x_1");
    oTags.push_back("Ex");
    oTags.push_back("n");
    //oTags.push_back("T_eV");

    vector<double> qs; qs.push_back(1.0); qs.push_back(-1.0);
    vector<double> ms; ms.push_back(1.0); ms.push_back(10.0);

//  Number of hamornics
    size_t my_ls[] = { 40, 28};   vector<size_t> ls(my_ls,my_ls+sizeof(my_ls)/sizeof(size_t));
                                  vector<size_t> Nm;

//  Grid in p-space
    size_t my_ps[] = { 122, 160}; vector<size_t> ps(my_ps,my_ps+sizeof(my_ps)/sizeof(size_t));
    vector<double> pmax; pmax.push_back(1.2); pmax.push_back(2.6);
    vector<double> pth;   pth.push_back(0.3); pth.push_back(0.15);

//  Output p-axis
    vector<size_t> Npx; Npx.push_back(122); Npx.push_back(118);

//  Output x-axis
    vector< size_t > Nx;   Nx.push_back(144);
    vector< double > xmin; xmin.push_back(-10.0);
    vector< double > xmax; xmax.push_back(10.0);

//  Move all of the data into a single container 
    Grid_Info grid(ls, Nm, ms, qs, 
                   xmin, xmax, Nx,
                   xmin, xmax, Nx,
                         pmax, ps,
                               Npx);

//  INITIALIZATION
    Export_Files::Xport out( grid.axis, oTags);     
    Export_Files::Restart_Facility Re;
    State1D Y( grid.axis.Nx(0), grid.Nl, grid.Np, grid.charge, grid.mass);
    Y = 0.0;    

//  Set up gaussians in the middle of the grid
    for (size_t s(0); s < qs.size(); ++s) {
        valarray<double> setGauss( Gaussian( grid.axis.Np(s),  grid.axis.pmin(s), grid.axis.pmax(s), pth[s]) );
        valarray<double> xax( grid.axis.x(0) );
        for (size_t i(0); i < grid.axis.Np(s); ++i) {
            for (size_t j(0); j < grid.axis.Nx(0); ++j) {
                Y.SH(s,0)(i,j) = max(0.00,exp( (-0.5) *double(xax[j])*double(xax[j]) ) )* setGauss[i];
            }
        }
    }


//  ITERATION LOOP TESTS-->
    int tout_start(0);
    double dt_out(0.2);
    double CLF(0.1);
    double numh  = size_t(static_cast<int>(dt_out/CLF))+1; 
    double h     = dt_out/static_cast<double>(numh);
    Algorithms::RK3<State1D> RK(Y);
    RKFunctor1D rkF(ls, pmax, ps, grid.axis.xmin(0), grid.axis.xmax(0), grid.axis.Nx(0));   
//    Re.Read(111,1,Y);
 
    Output_Data::Output_Preprocessor_1D  output( grid, oTags); 
    output( Y, tout_start);

//  ITERATION LOOP 
    for (size_t  t_out(tout_start+1); t_out< 101; ++t_out) {
        for (Clock W(t_out-1,dt_out,CLF); W.tick() < W.numh(); ++W) {
            Y = RK(Y, W.h(), &rkF);
            periodic_boundaries(Y,3);
            cout << "Time = " <<  W.time()<< "\n";

            // THE FOKKER-PLANCK OPERATOR GOES HERE
        }

        cout << "Output...\n";
        output( Y, t_out);
    }
//    Re.Write(111,2,Y);

}
//--------------------------------------------------------------------------- 




// Periodic boundary for a single node
void periodic_boundaries(State1D& Y,size_t  Nbc){

    // Harmonics:x0 "Right-Bound ---> Left-Guard" 
    for(size_t s(0); s < Y.Species(); ++s) {
        for(int i(0); i < Y.DF(s).dim(); ++i) {
            copy(Y.SH(s,i).xRB(Nbc), Y.SH(s,i).xRB(Nbc).end(), Y.SH(s,i).xLG(Nbc)); 
         }
    }
    // Fields:   x0 "Right-Bound ---> Left-Guard"
    for (size_t i(0); i < Nbc; ++i)  Y.Ex()(i) = Y.Ex()(Y.Ex().numx()-2*Nbc+i) ;

    // Harmonics:x0 "Left-Bound ---> Right-Guard" 
    for(size_t s(0); s < Y.Species(); ++s) {
        for(size_t i(0); i < Y.DF(s).dim(); ++i) {
            copy(Y.SH(s,i).xLB(Nbc), Y.SH(s,i).xLB(Nbc).end(),  Y.SH(s,i).xRG(Nbc));
        }
    } 
    // Fields:   x0 "Left-Bound ---> Right-Guard"
    for (size_t i(0); i < Nbc; ++i)  Y.Ex()(Y.Ex().numx()-Nbc+i) = Y.Ex()(i+Nbc);
}
//--------------------------------------------------------------------------- 

