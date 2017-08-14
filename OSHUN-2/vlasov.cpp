///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	May 17 2013
///////////////////////////////////////////////////////////

//   
//   This cpp file contains the definitions for the functions
//   required to export the data
///////////////////////////////////////////////////////////
//
// 
//   class Export_Formatted_Data::
//
//   This class receives the output matrices and saves the 
//   data in txt files with recognizable name after it attaches
//   a small header with information necessary for plotting. 
//
// 
//   class Restart_Facility::
//
//   This class writes restart files from each node, and 
//   reads restart files for each node.  
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <algorithm>
    #include <cstdlib>

    #include <math.h>
    #include <map>

//  My libraries
    #include "lib-array.h"
    #include "lib-algorithms.h"

//  Declerations
    #include "state.h"
    #include "formulary.h"
    #include "setup.h"
    #include "vlasov.h"


//**************************************************************
//--------------------------------------------------------------
//  Current
Current_1D::Current_1D( double pmin, double pmax, size_t Np, 
                        size_t Nx ) 
  : Jx(Nx),
    pr(Algorithms::MakeAxis(pmin,pmax,Np)), 
       invg(pr) {

       for (size_t i(0); i < pr.size(); ++i) { 
           invg[i] = 1.0 / (sqrt(1.0+pr[i]*pr[i]));
       }
};
//--------------------------------------------------------------

void Current_1D::operator()(const DistFunc1D& Din, Field1D& Exh) {
    SHarmonic1D f1( *Din(1) );
    f1 = f1.mpaxis( invg );
    for (size_t i(0); i < Jx.numx(); ++i) {
        Jx(i) = Algorithms::moment( f1.xVec(i), pr, 3);
    }
    Jx *= Din.q()/Din.mass() * 4.0 * M_PI / 3.0;
    Exh += Jx;
}
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Electric_Field_1D::Electric_Field_1D(size_t Nl, 
                         double pmin, double pmax, size_t Np, 
                         double xmin, double xmax, size_t Nx) 
   : A1(Nl), A2(Nl), Hp0(Nl),
     H(Np,Nx), G(Np,Nx),
     pr(Algorithms::MakeAxis(pmin,pmax,Np)),
     invpr(pr) {
//      - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Inverted momentum axis
         for (size_t i(0); i < pr.size(); ++i) { 
             invpr[i] = 1.0/pr[i]; 
         }
         double idp = (-1.0)/ (2.0*(pmax-pmin)/double(Np-1)); // -1/(2dp) 

//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       Calculate the A1 * l/(l+1) * 1/(2Dp), A2 * (-1)/(2Dp) parameters
         for (size_t l(1); l < Nl; ++l){
             double ld(l);
             A1[l] = idp *(-1.0) *  (ld+1.0)/(2.0*ld+1.0)  * ld/(ld+1.0);
             A2[l] = idp *                ld/(2.0*ld+1.0);
         }
         A1[0] = idp;
         // A2[0] is not used

//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       -2*Dp*H at the 0 momentum cell
         Hp0[0] = 1.0 / pr[0];
         for (size_t l(1); l < Nl; ++l) {
             double ld(l);
             Hp0[l] = Hp0[l-1] * (pr[0]/pr[1]) * (2.0*ld+1.0)/(2.0*ld-1.0);
         }
         Hp0 *= (-2.0)*(pr[1]-pr[0]);

    }
//--------------------------------------------------------------

void Electric_Field_1D::operator()(const DistFunc1D& Din, const Field1D& FEx, DistFunc1D& Dh) {

    valarray<double> Ex(FEx.array());
    Ex *= Din.q();
    size_t l0(A1.size()-1);

    MakeG00(*Din(0));
    Ex *= A1[0];  *Dh(1) += G.mxaxis(Ex);

//  1 < l < l0
    for (size_t l(1); l < l0; ++l){
        MakeGH(*Din(l),l);
        Ex *= A2[l] / A1[l-1];   *Dh(l-1) += H.mxaxis(Ex);
        Ex *= A1[l] / A2[l];     *Dh(l+1) += G.mxaxis(Ex);
    }

//  m = 0,  l = l0
    MakeGH(*Din(l0),l0);
    Ex *= A2[l0] / A1[l0-1];  *Dh(l0-1) += H.mxaxis(Ex);
}
//--------------------------------------------------------------

//  Make derivatives 2*Dp*(l+1/l)*G and -2*Dp*H for a given f
    void Electric_Field_1D::MakeGH(SHarmonic1D& f, size_t el){
        valarray<double> invpax(invpr); 
        double ld(el); 

        invpax *= (-2.0)*(ld+1.0) * (pr[1]-pr[0]);
       
        G = f;                   H = f;
        G = G.Dp();
                                 H  = H.mpaxis(invpax);
                                 H += G;
        G *= -(2.0*ld+1.0)/ld;
        G += H;

        for (size_t i(0); i < G.numx(); ++i) G(0,i) = 0.0;
        for (size_t i(0); i < H.numx(); ++i) H(0,i) = f(1,i) * Hp0[el];
    }
//--------------------------------------------------------------

//  Calculation of G00 = -2*Dp* df/dp(p0)
    void Electric_Field_1D::MakeG00(SHarmonic1D& f) {
        G = f; G = G.Dp();                            

        double p0p1_sq( pr[0]*pr[0]/(pr[1]*pr[1]) ),
               inv_mp0p1_sq( 1.0/(1.0-p0p1_sq) ),
               g_r = -4.0*(pr[1]-pr[0]) * pr[0]/(pr[1]*pr[1]),
               f00;

        for (size_t i(0); i < f.numx(); ++i) {
             f00    = ( f(0,i) - f(1,i) * p0p1_sq) * inv_mp0p1_sq;	
             G(0,i) = ( f(1,0) - f00) * g_r;
        }
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Spatial_Advection_1D::Spatial_Advection_1D(size_t Nl, 
                         double pmin, double pmax, size_t Np, 
                         double xmin, double xmax, size_t Nx) 
   : A1(Nl), A2(Nl),
     vr(Algorithms::MakeAxis(pmin,pmax,Np)) {
//      - - - - - - - - - - - - - - - - - - - - - - - - - - -

         for (size_t i(0); i < vr.size(); ++i) { 
             vr[i] = vr[i]/(sqrt(1.0+vr[i]*vr[i]));
         }

         double idx = (-1.0) / (2.0*(xmax-xmin)/double(Nx-1)); // -1/(2dx) 

//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
//       Calculate the "A1, A2" parameters
         for (size_t l(0); l < Nl; ++l){
             double ld(l);
             A1[l] = idx *(-1.0) * (ld+1.0) / (2.0*ld+1.0);
             A2[l] = idx *(-1.0) *  ld      / (2.0*ld+1.0);    
         }
         A2[0] = 1.0;

    }
//--------------------------------------------------------------

//   Advection in x
void Spatial_Advection_1D::operator()(const DistFunc1D& Din, DistFunc1D& Dh) {

        valarray<double> v_A1(vr); v_A1 *= 1.0/Din.mass();        
        valarray<double> v_A2(v_A1);        
        size_t l0(A1.size()-1);

        SHarmonic1D sh1( *Din(0) ), sh2(sh1);     
        sh1 = sh1.Dx(); // One copy of SH(0)
        v_A1 *= A1[0];      *Dh(1) += sh1.mpaxis(v_A1);

        for (size_t l(1); l < l0; ++l) {
            sh1 = *Din(l);
            sh2 = sh1.Dx(); // Two copies of SH(l)
            v_A1 *= A1[l]/A1[l-1];  *Dh(l+1) += sh1.mpaxis(v_A1);  
            v_A2 *= A2[l]/A2[l-1];  *Dh(l-1) += sh2.mpaxis(v_A2);  
        }

        sh2 = *Din(l0);            
        sh2 = sh2.Dx(); // One copy of SH(l0)
        v_A2 *= A2[l0]/A2[l0-1];  *Dh(l0-1) += sh2.mpaxis(v_A2);  

}
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Functor to be used in the Runge-Kutta methods 

//  Constructor
RKFunctor1D::RKFunctor1D(vector<size_t> Nl,vector<double> pmax, vector<size_t> Np, 
                         double xmin, double xmax, size_t Nx) {
    for (size_t s(0); s < Nl.size(); ++s){ 
        double pmin( pmax[s] / ( double(Np[s] * 2 - 1)) );
        SA.push_back( Spatial_Advection_1D(Nl[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );
        EF.push_back( Electric_Field_1D(Nl[s], pmin, pmax[s], Np[s], xmin, xmax, Nx) );
        JX.push_back( Current_1D(pmin, pmax[s], Np[s], Nx) );
    }
}

//  Collect all of the terms
void RKFunctor1D::operator()(const State1D& Yin, State1D& Yslope){
   Yslope = 0.0;
   for (size_t s(0); s < Yin.Species(); ++s) {
       SA[s](Yin.DF(s),Yslope.DF(s));
       EF[s](Yin.DF(s),Yin.Ex(),Yslope.DF(s));
       JX[s](Yin.DF(s),Yslope.Ex());
       Yslope.DF(s) = Yslope.DF(s).Filterp();
   }
}
//--------------------------------------------------------------
//**************************************************************
