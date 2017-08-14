///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//   Modified:  	Sep  12th 2009
//   Last Modified:	July 26th 2011
///////////////////////////////////////////////////////////

//   
//   This file contains the definitions for the  
//   classes ElectricField, MagneticField, Current, 
//   Spatial_Advection and MaxwellEq.
///////////////////////////////////////////////////////////
//
//   classes :
//
//   1. class ElectricField :
//      The class ElectricField  stores a reference to the slope Yh 
//      and steps it calls Exyz(Yin), where Yin is the current
//      state, in order to calculate the effect of the 
//      electric field on the distribution function.
//
//   2. class MagneticField :
//      The class Magneticfield stores a reference to the slope Yh 
//      and steps it calls Bxyz(Yin), where Yin is the current
//      state, in order to calculate the effect of the 
//      magnetic field on the distribution function.
//
//   3. class Current :
//      The class Current stores a reference to the slope Yh and
//      in subsequent steps it calls Jx(Yin), where Yin is the 
//      state, in order  to calculate effect of the current
//      on the Electric field.
//
//   4. Spatial_Advection :
//      The class Spatial_Advection  stores a reference to the slope Yh and
//      in subsequent steps it calls Axy(Yin), where Yin is the 
//      state, in order to calculate effect of the spatial 
//      advection on the distribution function. 
//
//   5. class MaxwellEq :
//      The class MaxwellEq stores a reference to the slope Yh and
//      in subsequent steps it calls Ampere(Yin) and Faraday(Yin),
//      where Yin is the state, in order to update the electromagnetic
//      fields 
// 
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

//  Standard Libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <math.h>
    #include <float.h>

//  My Libraries
    #include "matrices.h"

//  Declerations
    #include "decl-input.h"
    #include "decl-state.h"
    #include "decl-vlasovmax.h"




//**************************************************************
//**************************************************************
//   Definition for the Electric Field in the x direction
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    ElectricField::ElectricField(Stat& Yslope) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yslope),
         l0(Inputdata::IN().inp().l0), m0(Inputdata::IN().inp().m0), 
         A1(l0+1,m0+1), A2(l0+1,m0+1),
         B1(l0+1),      B2(l0+1),
         C1(l0+1),      C3(l0+1),
         C2(l0+1,m0+1), C4(l0+1,m0+1),
         H(Yslope.SH(0,0)), G(Yslope.SH(0,0)),TMP(Yslope.SH(0,0)),
         Hp0(l0+1),
         pr(Inputdata::IN().inp().pr.dim(), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(0)), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(Inputdata::IN().inp().pr.dim()-1))),
         invpr(pr)
         {
         complex<double> lc, mc;

//       Calculate the "A1, A2" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(1); l<l0+1; ++l)
             for (size_t m=0; m<((m0<l)?m0:l)+1; ++m){
                 lc = static_cast< complex<double> >(l);
                 mc = static_cast< complex<double> >(m);
                 A1(l,m) = (lc + 1.0 - mc) * lc / ( (2.0*lc + 1.0)*(lc+1.0) );
                 A2(l,m) = -(lc + mc) / (2.0*lc + 1.0);    
             }
         A1(0,0) = -1.0;
         A1 *= 0.5 / pr.dx(); A2 *= 0.5 / pr.dx();
//         A2(0,0) = 1.0; This has been commented out because you never get to use A2(m,m)
//         for (size_t m(1); m < m0+1; ++m) A2(m,m) = A2(l0,m-1);
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "B1, B2" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(1); l<l0+1; ++l){
             lc = static_cast< complex<double> >(l);
             B1[l] =(-0.5/pr.dx()) * lc* lc     /(2.0*lc+1.0);
             B2[l] =(-0.5/pr.dx()) * lc*(lc+1.0)/(2.0*lc+1.0);
         }
         B2[0] = 1.0;
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
//       Calculate the "C1, C3" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(1); l<l0+1; ++l){
             lc = static_cast< complex<double> >(l);
             C1[l] =(0.25/pr.dx()) * lc /((2.0*lc+1.0)*(lc+1.0));
             C3[l] =(0.25/pr.dx()) /(2.0*lc+1.0);
         }        
         C1[0] = -0.25/pr.dx();
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "C2, C4" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(2); l<l0+1; ++l){
             for (size_t m(2); m<((m0<l)?m0:l)+1; ++m){
                 lc = static_cast< complex<double> >(l);
                 mc = static_cast< complex<double> >(m);
                 C2(l,m) =(-0.25/pr.dx()) * lc * (lc-mc+2.0)*(lc-mc+1.0)/((2.0*lc+1.0)*(lc+1.0));
                 C4(l,m) =(-0.25/pr.dx()) * (lc+mc-1.0)*(lc+mc)/(2.0*lc+1.0);
             }
         }   
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -


//       Calculate the "H" parameters at the p0 cell
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Hp0 = 1.0 / pr(0); 
         for (size_t l(1); l < l0+1; ++l) { Hp0[l] = Hp0[l-1]*(pr(0)/pr(1)) * (2.0*l+1.0)/(2.0*l-1);}
         Hp0 *= -2.0 * pr.dx(); 
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Create inverted momentum axis   
         for (size_t i(0); i<invpr.dim(); ++i) { invpr(i) = 1.0/pr(i); } 
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& ElectricField::Implicit_Ex(Stat& Yin, Matrix2D< complex<double> >& Ex){ 
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in x,
//  which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

//      m = 0, l = 0
        G  = G00(Yin.SH(0,0)); 
        Ex *= A1(0,0);             Yh.SH(1,0) += G.mxy_matrix(Ex);

//      m = 0, l = 1
        MakeGH(1,Yin.SH(1,0));
        Ex *= A2(1,0) / A1(0,0);   Yh.SH(0,0) += H.mxy_matrix(Ex);
        Ex *= A1(1,0) / A2(1,0);   Yh.SH(2,0) += G.mxy_matrix(Ex);

//      m = 0, l = 2
        MakeGH(2,Yin.SH(2,0));
        Ex *= A2(2,0) / A1(1,0);   Yh.SH(1,0) += H.mxy_matrix(Ex);

//      m = 0, l = 3
        MakeGH(3,Yin.SH(3,0));
        Ex *= A2(3,0) / A2(2,0);   Yh.SH(2,0) += H.mxy_matrix(Ex);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(1,Yin.SH(1,1));
        Ex *= A1(1,1) / A2(3,0);   Yh.SH(2,1) += G.mxy_matrix(Ex); 

//      m = 1, l = 2
        MakeGH(2,Yin.SH(2,1));
        Ex *= A2(2,1) / A1(1,1);   Yh.SH(1,1) += H.mxy_matrix(Ex);

//      m = 1, l = 3
        MakeGH(3,Yin.SH(3,1));
        Ex *= A2(3,1) / A2(2,1);   Yh.SH(2,1) += H.mxy_matrix(Ex);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        return Yh;
    }
//--------------------------------------------------------------



//--------------------------------------------------------------
    Stat& ElectricField::Implicit_Ex1D(Stat& Yin, Matrix2D< complex<double> >& Ex){ 
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in x,
//  which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

//      m = 0, l = 0
        G  = G00(Yin.SH(0,0)); 
        Ex *= A1(0,0);             Yh.SH(1,0) += G.mxy_matrix(Ex);

//      m = 0, l = 1
        MakeGH(1,Yin.SH(1,0));
        Ex *= A2(1,0) / A1(0,0);   Yh.SH(0,0) += H.mxy_matrix(Ex);
        Ex *= A1(1,0) / A2(1,0);   Yh.SH(2,0) += G.mxy_matrix(Ex);

//      m = 0, l = 2
        MakeGH(2,Yin.SH(2,0));
        Ex *= A2(2,0) / A1(1,0);   Yh.SH(1,0) += H.mxy_matrix(Ex);

//      m = 0, l = 3
        MakeGH(3,Yin.SH(3,0));
        Ex *= A2(3,0) / A2(2,0);   Yh.SH(2,0) += H.mxy_matrix(Ex);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        return Yh;
    }
//--------------------------------------------------------------



//--------------------------------------------------------------
    Stat& ElectricField::Implicit_Eyz(Stat& Yin, Matrix2D< complex<double> >& Em, 
                                                 Matrix2D< complex<double> >& Ep){ 
//--------------------------------------------------------------
//  This is the calculation for the implicit electric field in y
//  and z, which is to say Ex as it acts on the first few harmonics.
//--------------------------------------------------------------

//      m = 0, l = 0
        G  = G00(Yin.SH(0,0)); 
        Em *= C1[0];            Yh.SH(1,1) += G.mxy_matrix(Em);

//      m = 0, l = 1
        MakeGH(1,Yin.SH(1,0));
        Em *= C1[1]   / C1[0];  Yh.SH(2,1) += G.mxy_matrix(Em);

//      m = 0, 1 < l < l0
        MakeGH(2,Yin.SH(2,0));
        Em *= C3[2]   / C1[1];  Yh.SH(1,1) += H.mxy_matrix(Em);

//      m = 0,  l = 3 
        MakeGH(3,Yin.SH(3,0));
        Em *= C3[3]   / C3[2];  Yh.SH(2,1) += H.mxy_matrix(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(1,Yin.SH(1,1));
        Ep *= B2[1];             H = H.mxy_matrix(Ep); Yh.SH(0,0) += H.Re();
        Em *= C1[1] / C3[3];     TMP = G;              Yh.SH(2,2) += TMP.mxy_matrix(Em); 
        Ep *= B1[1] / B2[1];     G = G.mxy_matrix(Ep); Yh.SH(2,0) += G.Re();

//      m = 1, l = 2
        MakeGH(2,Yin.SH(2,1));
        Ep *= B2[2]   / B1[1];   H = H.mxy_matrix(Ep); Yh.SH(1,0) += H.Re();

//      m = 1, l = 3
        MakeGH(3,Yin.SH(3,1));
        Em *= C3[3]   / C1[1];   TMP = H;              Yh.SH(2,2) += TMP.mxy_matrix(Em);
        Ep *= B2[3]   / B2[2];   H = H.mxy_matrix(Ep); Yh.SH(2,0) += H.Re();
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 2, l = 2
        MakeGH(2,Yin.SH(2,2));
        Ep *= C4(2,2) / B2[3];        Yh.SH(1,1) += H.mxy_matrix(Ep); 

//      m = 2, l = 3
        MakeGH(3,Yin.SH(3,2));
        Ep *= C4(3,2) / C4(2,2);      Yh.SH(2,1)  += H.mxy_matrix(Ep);

        if ( m0 > 2) { 
//          m = 3, l = 3
            MakeGH(3,Yin.SH(3,3));
            Ep *= C4(3,3) / C4(3,2);  Yh.SH(2,2)  += H.mxy_matrix(Ep);
        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        return Yh;
    }
//--------------------------------------------------------------



//--------------------------------------------------------------
    Stat& ElectricField::Exyz1D(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for the electric field
//--------------------------------------------------------------
        complex<double> ii(0.0,1.0);
        Matrix2D< complex<double> > Ex(Yin.EMF().Ex().matrix());
        Matrix2D< complex<double> > Em(Yin.EMF().Ez().matrix()); Em *= (-1.0)*ii; Em += Yin.EMF().Ey().matrix();
        Matrix2D< complex<double> > Ep(Yin.EMF().Ez().matrix()); Ep *=  ii;       Ep += Yin.EMF().Ey().matrix();


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 0
        G  = G00(Yin.SH(0,0)); 
        Ex *= A1(0,0); TMP = G; Yh.SH(1,0) += G.mxy_matrix(Ex);
//        Em *= C1[0];            Yh.SH(1,1) += TMP.mxy_matrix(Em);

//      m = 0, l = 1
        MakeGH(1,Yin.SH(1,0));
        Ex *= A2(1,0) / A1(0,0);          Yh.SH(0,0) += H.mxy_matrix(Ex);
        Ex *= A1(1,0) / A2(1,0); TMP = G; Yh.SH(2,0) += G.mxy_matrix(Ex);
//        Em *= C1[1]   / C1[0]  ;          Yh.SH(2,1) += TMP.mxy_matrix(Em);

//      m = 0, 1 < l < l0
        for (size_t l(2); l < l0; ++l){
            MakeGH(l,Yin.SH(l,0));
            Ex *= A2(l,0) / A1(l-1,0); TMP = H;  Yh.SH(l-1,0) += H.mxy_matrix(Ex);
//            Em *= C3[l]   / C1[l-1];             Yh.SH(l-1,1) += TMP.mxy_matrix(Em);
            Ex *= A1(l,0) / A2(l,0);   TMP = G;  Yh.SH(l+1,0) += G.mxy_matrix(Ex);
//            Em *= C1[l]   / C3[l];               Yh.SH(l+1,1) += TMP.mxy_matrix(Em);
        }
//      m = 0,  l = l0
        MakeGH(l0,Yin.SH(l0,0));
        Ex *= A2(l0,0) / A1(l0-1,0); TMP = H;  Yh.SH(l0-1,0) += H.mxy_matrix(Ex);
//        Em *= C3[l0]   / C1[l0-1];             Yh.SH(l0-1,1) += TMP.mxy_matrix(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// This is the end for the 1D version
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


        return Yh;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& ElectricField::Exyz(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for the electric field
//--------------------------------------------------------------
        complex<double> ii(0.0,1.0);
        Matrix2D< complex<double> > Ex(Yin.EMF().Ex().matrix());
        Matrix2D< complex<double> > Em(Yin.EMF().Ez().matrix()); Em *= (-1.0)*ii; Em += Yin.EMF().Ey().matrix();
        Matrix2D< complex<double> > Ep(Yin.EMF().Ez().matrix()); Ep *=  ii;       Ep += Yin.EMF().Ey().matrix();


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 0, l = 0
        G  = G00(Yin.SH(0,0)); 
        Ex *= A1(0,0); TMP = G; Yh.SH(1,0) += G.mxy_matrix(Ex);
        Em *= C1[0];            Yh.SH(1,1) += TMP.mxy_matrix(Em);

//      m = 0, l = 1
        MakeGH(1,Yin.SH(1,0));
        Ex *= A2(1,0) / A1(0,0);          Yh.SH(0,0) += H.mxy_matrix(Ex);
        Ex *= A1(1,0) / A2(1,0); TMP = G; Yh.SH(2,0) += G.mxy_matrix(Ex);
        Em *= C1[1]   / C1[0]  ;          Yh.SH(2,1) += TMP.mxy_matrix(Em);

//      m = 0, 1 < l < l0
        for (size_t l(2); l < l0; ++l){
            MakeGH(l,Yin.SH(l,0));
            Ex *= A2(l,0) / A1(l-1,0); TMP = H;  Yh.SH(l-1,0) += H.mxy_matrix(Ex);
            Em *= C3[l]   / C1[l-1];             Yh.SH(l-1,1) += TMP.mxy_matrix(Em);
            Ex *= A1(l,0) / A2(l,0);   TMP = G;  Yh.SH(l+1,0) += G.mxy_matrix(Ex);
            Em *= C1[l]   / C3[l];               Yh.SH(l+1,1) += TMP.mxy_matrix(Em);
        }
//      m = 0,  l = l0
        MakeGH(l0,Yin.SH(l0,0));
        Ex *= A2(l0,0) / A1(l0-1,0); TMP = H;  Yh.SH(l0-1,0) += H.mxy_matrix(Ex);
        Em *= C3[l0]   / C1[l0-1];             Yh.SH(l0-1,1) += TMP.mxy_matrix(Em);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      m = 1, l = 1
        MakeGH(1,Yin.SH(1,1));
        Ep *= B2[1];              H = H.mxy_matrix(Ep); Yh.SH(0,0) += H.Re();
        Ex *= A1(1,1) / A2(l0,0); TMP = G;              Yh.SH(2,1) += G.mxy_matrix(Ex); 
        Em *= C1[1] / C3[l0];     G = TMP;              Yh.SH(2,2) += TMP.mxy_matrix(Em); 
        Ep *= B1[1] / B2[1];      G = G.mxy_matrix(Ep); Yh.SH(2,0) += G.Re();

//      m = 1, l = 2
        MakeGH(2,Yin.SH(2,1));
        Ex *= A2(2,1) / A1(1,1); TMP = H;              Yh.SH(1,1) += TMP.mxy_matrix(Ex);
        Ep *= B2[2]   / B1[1];   H = H.mxy_matrix(Ep); Yh.SH(1,0) += H.Re();
        Ex *= A1(2,1) / A2(2,1); TMP = G;              Yh.SH(3,1) += G.mxy_matrix(Ex);
        Em *= C1[2]   / C1[1];   G = TMP;              Yh.SH(3,2) += TMP.mxy_matrix(Em);
        Ep *= B1[2]   / B2[2];   G = G.mxy_matrix(Ep); Yh.SH(3,0) += G.Re();

//      m = 1, 1 < l < l0
        for (size_t l(3); l < l0; ++l){
            MakeGH(l,Yin.SH(l,1));
            Ex *= A2(l,1) / A1(l-1,1); TMP = H;              Yh.SH(l-1,1) += H.mxy_matrix(Ex);
            Em *= C3[l]   / C1[l-1];   H = TMP;              Yh.SH(l-1,2) += TMP.mxy_matrix(Em);
            Ep *= B2[l]   / B1[l-1];   H = H.mxy_matrix(Ep); Yh.SH(l-1,0) += H.Re();
            Ex *= A1(l,1) / A2(l,1);   TMP = G;              Yh.SH(l+1,1) += G.mxy_matrix(Ex);
            Em *= C1[l]   / C3[l];     G = TMP;              Yh.SH(l+1,2) += TMP.mxy_matrix(Em);
            Ep *= B1[l]   / B2[l];     G = G.mxy_matrix(Ep); Yh.SH(l+1,0) += G.Re();
         }
//       m = 1,  l = l0
         MakeGH(l0,Yin.SH(l0,1));
         Ex *= A2(l0,1) / A1(l0-1,1); TMP = H;              Yh.SH(l0-1,1) += H.mxy_matrix(Ex);
         Em *= C3[l0]   / C1[l0-1];   H = TMP;              Yh.SH(l0-1,2) += TMP.mxy_matrix(Em);
         Ep *= B2[l0]   / B1[l0-1];   H = H.mxy_matrix(Ep); Yh.SH(l0-1,0) += H.Re();
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        C4(l0,1) = B2[l0];
        for (size_t m(2); m < m0; ++m){
//          m > 1 , l = m
            MakeGH(m,Yin.SH(m,m));
            Ep *= C4(m,m) / C4(l0,m-1);              Yh.SH(m-1,m-1) += H.mxy_matrix(Ep); 
            Ex *= A1(m,m) / A2(l0,m-1); TMP = G;     Yh.SH(m+1,m  ) += G.mxy_matrix(Ex); 
            Em *= C1[m]   / C3[l0];     G = TMP;     Yh.SH(m+1,m+1) += TMP.mxy_matrix(Em); 
            Ep *= C2(m,m) / C4(m,m);                 Yh.SH(m+1,m-1) += G.mxy_matrix(Ep);

//          m > 1 , l = m+1
            MakeGH(m+1,Yin.SH(m+1,m));
            Ex *= A2(m+1,m) / A1(m,m);   TMP = H;     Yh.SH(m  ,m  ) += TMP.mxy_matrix(Ex);
            Ep *= C4(m+1,m) / C2(m,m);                Yh.SH(m  ,m-1) += H.mxy_matrix(Ep);
            if ( m+1 < l0) { //always true except when m = m0-1 = l0-1 
                Ex *= A1(m+1,m) / A2(m+1,m); TMP = G;     Yh.SH(m+2,m  ) += G.mxy_matrix(Ex);
                Em *= C1[m+1]   / C1[m];     G = TMP;     Yh.SH(m+2,m+1) += TMP.mxy_matrix(Em);
                Ep *= C2(m+1,m) / C4(m+1,m);              Yh.SH(m+2,m-1) += G.mxy_matrix(Ep);

//              m > 1, 1 < l < l0
                for (size_t l(m+2); l < l0; ++l){
                    MakeGH(l,Yin.SH(l,m));
                    Ex *= A2(l,m) / A1(l-1,m); TMP = H;   Yh.SH(l-1,m  ) += H.mxy_matrix(Ex);
                    Em *= C3[l]   / C1[l-1];   H = TMP;   Yh.SH(l-1,m+1) += TMP.mxy_matrix(Em);
                    Ep *= C4(l,m) / C2(l-1,m);            Yh.SH(l-1,m-1) += H.mxy_matrix(Ep);
                    Ex *= A1(l,m) / A2(l,m);   TMP = G;   Yh.SH(l+1,m  ) += G.mxy_matrix(Ex);
                    Em *= C1[l]   / C3[l];     G = TMP;   Yh.SH(l+1,m+1) += TMP.mxy_matrix(Em);
                    Ep *= C2(l,m) / C4(l,m);              Yh.SH(l+1,m-1) += G.mxy_matrix(Ep); 
                 }
//               m > 1,  l = l0
                 MakeGH(l0,Yin.SH(l0,m));
                 Ex *= A2(l0,m) / A1(l0-1,m); TMP = H;    Yh.SH(l0-1,m  ) += H.mxy_matrix(Ex);
                 Em *= C3[l0]   / C1[l0-1];   H = TMP;    Yh.SH(l0-1,m+1) += TMP.mxy_matrix(Em);
                 Ep *= C4(l0,m) / C2(l0-1,m);             Yh.SH(l0-1,m-1) += H.mxy_matrix(Ep);
             }
         }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         MakeGH(m0,Yin.SH(m0,m0));
         Ep *= C4(m0,m0) / C4(l0,m0-1);              Yh.SH(m0-1,m0-1) += H.mxy_matrix(Ep); 
//       m = m0, l0 > l > m0
         if ( m0 < l0) { 
            Ex *= A1(m0,m0) / A2(l0,m0-1); TMP = G;  Yh.SH(m0+1,m0  ) += TMP.mxy_matrix(Ex); 
            Ep *= C2(m0,m0) / C4(m0,m0);             Yh.SH(m0+1,m0-1) += G.mxy_matrix(Ep);

//          m = m0 , l = m0+1
            MakeGH(m0+1,Yin.SH(m0+1,m0));
            Ex *= A2(m0+1,m0) / A1(m0,m0); TMP = H;        Yh.SH(m0,m0  )  += TMP.mxy_matrix(Ex);
            Ep *= C4(m0+1,m0) / C2(m0,m0);                 Yh.SH(m0,m0-1)  += H.mxy_matrix(Ep);
            if ( m0+1 < l0) { 
                Ex *= A1(m0+1,m0) / A2(m0+1,m0); TMP = G;  Yh.SH(m0+2,m0  )+= TMP.mxy_matrix(Ex);
                Ep *= C2(m0+1,m0) / C4(m0+1,m0);           Yh.SH(m0+2,m0-1)+= G.mxy_matrix(Ep);

//              m = m0, m0+2 < l < l0
                for (size_t l(m0+2); l < l0; ++l){
                    MakeGH(l,Yin.SH(l,m0));
                    Ex *= A2(l,m0) / A1(l-1,m0); TMP = H;  Yh.SH(l-1,m0)   += TMP.mxy_matrix(Ex);
                    Ep *= C4(l,m0) / C2(l-1,m0);           Yh.SH(l-1,m0-1) += H.mxy_matrix(Ep);
                    Ex *= A1(l,m0) / A2(l,  m0); TMP = G;  Yh.SH(l+1,m0  ) += TMP.mxy_matrix(Ex);
                    Ep *= C2(l,m0) / C4(l,  m0);           Yh.SH(l+1,m0-1) += G.mxy_matrix(Ep); 
                 }

//               m > 1,  l = l0
                 MakeGH(l0,Yin.SH(l0,m0));
                 Ex *= A2(l0,m0) / A1(l0-1,m0); TMP = H;    Yh.SH(l0-1,m0)   += TMP.mxy_matrix(Ex);
                 Ep *= C4(l0,m0) / C2(l0-1,m0);             Yh.SH(l0-1,m0-1) += H.mxy_matrix(Ep);
             } 
        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


        return Yh;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void ElectricField::MakeGH(size_t el, SHarmonic& f){
//--------------------------------------------------------------
        GSlice_iter< complex<double> >  it1(f.p0(1));
        Axis< complex<double> > invpax(invpr); 
        complex<double> lc(static_cast< complex<double> >(el)); 

        invpax *= -2.0*(lc+1.0)*pr.dx();
       
        G = f;                   H = f;
        G = G.Dp();
                                 H  = H.mpaxis(invpax.array());
                                 H += G;
        G *= -(2.0*lc+1.0)/lc;
        G += H;

        for ( GSlice_iter< complex<double> > itG(G.p0(0)); itG!=itG.end(); ++itG)
            *itG = 0.0;
        for ( GSlice_iter< complex<double> > itH(H.p0(0)); itH!=itH.end(); ++itH)
            *itH = (*it1++)*Hp0[el]; 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    SHarmonic& ElectricField::G00(SHarmonic& f){ 
//--------------------------------------------------------------
//  Calculation of G00
//--------------------------------------------------------------
        complex<double> p0p1_sq(pr(0)*pr(0)/(pr(1)*pr(1))),
                        inv_mp0p1_sq(1.0/(1.0-p0p1_sq)),
                        g_r = -4.0*pr.dx()*pr(0)/(pr(1)*pr(1)),
                        f00;
        GSlice_iter< complex<double> > it0(f.p0(0)), it1(f.p0(1));

        G = f; G = G.Dp();                            
        for ( GSlice_iter< complex<double> > itG(G.p0(0)); itG!=itG.end(); ++itG){
             f00    =  ( (*it0++) - (*it1)*p0p1_sq)*inv_mp0p1_sq;	
             *itG =   ((*it1++)-f00)*g_r;
        }
        return G;
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Definition for the Magnetic Field
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    MagneticField::MagneticField(Stat& Yslope) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yslope),
         l0(Inputdata::IN().inp().l0), m0(Inputdata::IN().inp().m0), 
         A1(m0+1),      B1(l0+1),
         A2(l0+1,m0+1),
         A3(0.5),
         FLM(Yslope.SH(0,0)) 
         {
         complex<double> lc, mc;
         complex<double> c01(0.0,1.0);

//       Calculate the "A1" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t m=0; m<m0+1; ++m){
             mc = static_cast< complex<double> >(m);
             A1[m] = (-1.0)*c01*mc;
         }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "A2" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(1); l<l0+1; ++l)
             for (size_t m=0; m<((m0<l)?m0:l)+1; ++m){
                 lc = static_cast< complex<double> >(l);
                 mc = static_cast< complex<double> >(m);
                 A2(l,m) = (-0.5)*(lc+1.0-mc)*(lc+mc);
             }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "A3" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         A3 = 0.5;
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "B1" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(0); l<l0+1; ++l){
            lc = static_cast< complex<double> >(l);
            B1[l] = (-1.0)*lc*(lc+1.0); 
         }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& MagneticField::Bxyz(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for the magnetic field
//--------------------------------------------------------------
        complex<double> ii(0.0,1.0);
        Matrix2D< complex<double> > Bx(Yin.EMF().Bx().matrix());
        Matrix2D< complex<double> > Bm(Yin.EMF().By().matrix()); Bm *= (-1.0)*ii; Bm += Yin.EMF().Bz().matrix();
        Matrix2D< complex<double> > Bp(Yin.EMF().By().matrix()); Bp *=  ii; Bp += Yin.EMF().Bz().matrix();

//      m = 0, 1 < l < l0+1
        Bp *= A3;
        for (size_t l(1); l < l0+1; ++l){
            FLM = Yin.SH(l,0);      Yh.SH(l,1) += FLM.mxy_matrix(Bp);
        }

//      m = 1, l = 1
        FLM = Yin.SH(1,1); Bx *= A1[1];                           Yh.SH(1,1) += FLM.mxy_matrix(Bx); 
        FLM = Yin.SH(1,1); Bm *= B1[1]; FLM = FLM.mxy_matrix(Bm); Yh.SH(1,0) += FLM.Re(); 
//      m = 1, l > 1
        for (size_t l(2); l < l0+1; ++l){
            FLM = Yin.SH(l,1);                                                Yh.SH(l,2) += FLM.mxy_matrix(Bp); 
            FLM = Yin.SH(l,1);                                                Yh.SH(l,1) += FLM.mxy_matrix(Bx); 
            FLM = Yin.SH(l,1); Bm *= B1[l]/B1[l-1]; FLM = FLM.mxy_matrix(Bm); Yh.SH(l,0) += FLM.Re(); 
        }        
        Bm *= 1.0/B1[l0];        

//      m > 1, l = m
        for (size_t m(2); m < m0; ++m){
            FLM = Yin.SH(m,m); Bx *= A1[m]/A1[m-1];                   Yh.SH(m,m  ) += FLM.mxy_matrix(Bx); 
            FLM = Yin.SH(m,m); Bm *= A2(m,m);                         Yh.SH(m,m-1) += FLM.mxy_matrix(Bm); 
            for (size_t l(m+1); l < l0+1; ++l){
                FLM = Yin.SH(l,m);                                    Yh.SH(l,m+1) += FLM.mxy_matrix(Bp); 
                FLM = Yin.SH(l,m);                                    Yh.SH(l,m  ) += FLM.mxy_matrix(Bx); 
                FLM = Yin.SH(l,m); Bm *= A2(l,m)/A2(l-1,m);           Yh.SH(l,m-1) += FLM.mxy_matrix(Bm); 
            }
            Bm *= 1.0/A2(l0,m);
        }

//      m = m0, l >= m0
        FLM = Yin.SH(m0,m0); Bx *= A1[m0]/A1[m0-1];                   Yh.SH(m0,m0)   += FLM.mxy_matrix(Bx); 
        FLM = Yin.SH(m0,m0); Bm *= A2(m0,m0)/*/A2(l0,m0-1)*/;             Yh.SH(m0,m0-1) += FLM.mxy_matrix(Bm);
        for (size_t l(m0+1); l < l0+1; ++l){
            FLM = Yin.SH(l,m0);                                     Yh.SH(l,m0  )  += FLM.mxy_matrix(Bx); 
            FLM = Yin.SH(l,m0); Bm *= A2(l,m0)/A2(l-1,m0);          Yh.SH(l,m0-1)  += FLM.mxy_matrix(Bm); 
        }

        return Yh;
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Definition for the Current in the x direction
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    Current::Current(Stat& Yslope) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yslope),
         jayx(Yslope.FLD(0)), 
         jayy(Yslope.FLD(0)), 
         jayz(Yslope.FLD(0)), 
         p3og(Inputdata::IN().inp().pr.dim(), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(0)), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(Inputdata::IN().inp().pr.dim()-1))),
         Delta_p(static_cast< complex<double> >(Inputdata::IN().inp().pr.dx())),
         small(static_cast< complex<double> >(Inputdata::IN().inp().pr(0)))
         {
             for (size_t i(0); i<p3og.dim(); ++i) // p^3/gamma 
                p3og(i) = (p3og(i)*p3og(i))*(p3og(i)/sqrt(1.0+p3og(i)*p3og(i))); 

             small *= small; small *= small; small *= Inputdata::IN().inp().pr(0); 
             small *= 0.2; small *= 1.0/(Inputdata::IN().inp().pr(1)); 
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& Current::Jxyz(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for the electric field
//--------------------------------------------------------------
        complex<double> c01(0.0,1.0);
        SHarmonic f10(Yin.SH(1,0)); 
        SHarmonic f11(Yin.SH(1,1)); 
        f10  = f10.mpaxis(p3og.array());
        f11  = f11.mpaxis(p3og.array());

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//      Basic integration
        for (size_t ix(0); ix < f10.numx(); ++ix){
            for (size_t iy(0); iy < f10.numy(); ++iy){
                jayx(ix,iy)  = f10(0,ix,iy);
                jayx(ix,iy) += f10(f10.nump()-1,ix,iy);
                jayx(ix,iy) *= 0.5;
                for (size_t ip (1); ip < f10.nump()-1; ++ip) jayx(ix,iy)+= f10(ip,ix,iy);
            }
        }
        jayx *= Delta_p;

//      Add tiny contribution from 0-p0
        for (size_t ix(0); ix < f10.numx(); ++ix){
            for (size_t iy(0); iy < f10.numy(); ++iy){
                jayx(ix,iy) += small*f10(1,ix,iy);
             }
        }

//      Calculate the final minus current
        jayx *= 4.0 * M_PI / 3.0; 

        Yh.EMF().Ex() += jayx;// (minus minus ---> +current)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//      Basic integration
        for (size_t ix(0); ix < f11.numx(); ++ix){
            for (size_t iy(0); iy < f11.numy(); ++iy){
                jayz(ix,iy)  = f11(0,ix,iy);
                jayz(ix,iy) += f11(f11.nump()-1,ix,iy);
                jayz(ix,iy) *= 0.5;
                for (size_t ip (1); ip < f11.nump()-1; ++ip) jayz(ix,iy)+= f11(ip,ix,iy);
            }
        }
        jayz *= Delta_p;

//      Add tiny contribution from 0-p0
        for (size_t ix(0); ix < f11.numx(); ++ix){
            for (size_t iy(0); iy < f11.numy(); ++iy){
                jayz(ix,iy) += small*f11(1,ix,iy);
             }
        }

//      Calculate the final minus complex current 
        jayz *= 8.0 * M_PI / 3.0;

        jayy = jayz;      jayy  = jayy.Re();
        jayz -= jayy;     jayz *= c01;
        
//      minus minus ---> +current        
        Yh.EMF().Ey() += jayy;
        Yh.EMF().Ez() += jayz;

        return Yh;
        
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Definition for the Spatial Advection 
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    Spatial_Advection::Spatial_Advection(Stat& Yslope) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yslope),
         fd1(Yslope.SH(0,0)), fd2(Yslope.SH(0,0)),
         l0(Inputdata::IN().inp().l0), m0(Inputdata::IN().inp().m0), 
         A1(l0+1,m0+1), A2(l0+1,m0+1),
         B1(l0+1),      B2(l0+1),
         C1(l0+1),      C3(l0+1),
         C2(l0+1,m0+1), C4(l0+1,m0+1),
         vr(Inputdata::IN().inp().pr.dim(), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(0)), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(Inputdata::IN().inp().pr.dim()-1))),
         idx(static_cast< complex<double> >(Inputdata::IN().inp().x.dx())),
         idy(static_cast< complex<double> >(Inputdata::IN().inp().y.dx()))
         {

         for (size_t i(0); i<vr.dim(); ++i) vr(i) = vr(i)/(sqrt(1.0+vr(i)*vr(i)));
         idx = -1.0/(2.0*idx); 
         idy = -1.0/(2.0*idy); 

         complex<double> lc, mc;

//       Calculate the "A1, A2" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(0); l<l0+1; ++l){
             for (size_t m=0; m<((m0<l)?m0:l)+1; ++m){
                 lc = static_cast< complex<double> >(l);
                 mc = static_cast< complex<double> >(m);
                 A1(l,m) = idx *(-1.0) * (lc-mc+1.0) / (2.0*lc+1.0);
                 A2(l,m) = idx *(-1.0) * (lc+mc)     / (2.0*lc+1.0);    
             }
         }

//         A1(0,0) = -1.0;
//         A1 *= 0.5 / pr.dx(); A2 *= 0.5 / pr.dx();
//         A2(0,0) = 1.0;
//         for (size_t m(1); m < m0+1; ++m) A2(m,m) = A2(l0,m-1);
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "B1, B2" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(0); l<l0+1; ++l){
             lc = static_cast< complex<double> >(l);
             B1[l] = idy * (lc + 1.0) * lc / (2.0*lc + 1.0);
             B2[l] = (-1.0)*B1[l];
         }

//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "C1, C3" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(0); l<l0+1; ++l){
             lc = static_cast< complex<double> >(l);
             C1[l] = (-0.5) * idy / (2.0*lc + 1.0);
             C3[l] = (-1.0) * C1[l];
         }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

//       Calculate the "C2, C4" parameters
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t l(0); l<l0+1; ++l){
             for (size_t m=0; m<((m0<l)?m0:l)+1; ++m){
                 lc = static_cast< complex<double> >(l);
                 mc = static_cast< complex<double> >(m);
                 C2(l,m) = idy * 0.5 * (lc + 2.0 - mc)*(lc - mc + 1.0) / (2.0*lc + 1.0);
                 C4(l,m) = idy * (-0.5) * (lc + mc - 1.0)*(lc + mc) / (2.0*lc + 1.0);
             }
         }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& Spatial_Advection::Axy(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for the spatial advection
//--------------------------------------------------------------
        Axis< complex<double> > vt(vr);        

//      Advection in x
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for (size_t m(0); m < ((m0<l0)?(m0+1):m0); ++m){
            fd1 = Yin.SH(m,m);     fd1 = fd1.Dx();
            vt *= A1(m,m);                          Yh.SH(m+1,m) += fd1.mpaxis(vt.array());
            for (size_t l(m+1); l < l0; ++l) {
               fd1 = Yin.SH(l,m);          fd1 = fd1.Dx();  
               vt *= A2(l,m)/A1(l-1,m);    fd2 = fd1;  Yh.SH(l-1,m) += fd1.mpaxis(vt.array());  
               vt *= A1(l,m)/A2(l  ,m);                Yh.SH(l+1,m) += fd2.mpaxis(vt.array());  
            }
            fd1 = Yin.SH(l0,m);            fd1 = fd1.Dx();
            vt *= A2(l0,m)/A1(l0-1,m);                 Yh.SH(l0-1,m) += fd1.mpaxis(vt.array());  
            vt *= 1.0/A2(l0,m); 
        }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
            
//       m = 0, advection in y
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         fd1 = Yin.SH(0,0);     fd1 = fd1.Dy();
         vt *= C1[0];                            Yh.SH(1,1) += fd1.mpaxis(vt.array());
         fd1 = Yin.SH(1,0);     fd1 = fd1.Dy();
         vt *= C1[1]/C1[0];                      Yh.SH(2,1) += fd1.mpaxis(vt.array());
         for (size_t l(2); l < l0; ++l) {
             fd1 = Yin.SH(l,0);     fd1 = fd1.Dy();
             vt *= C3[l]/C1[l-1];   fd2 = fd1;   Yh.SH(l-1,1) += fd1.mpaxis(vt.array());
             vt *= C1[l]/C3[l];                  Yh.SH(l+1,1) += fd2.mpaxis(vt.array());
         }
         fd1 = Yin.SH(l0,0);    fd1 = fd1.Dy();
         vt *= C3[l0]/C1[l0-1];              Yh.SH(l0-1,1) += fd1.mpaxis(vt.array());
         vt *= 1.0   /C3[l0];
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -

         
//       m = 1, advection in y
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         fd1 = Yin.SH(1,1);     fd1 = fd1.Dy();
         vt *= C1[1];           fd2 = fd1;                                 Yh.SH(2,2) += fd1.mpaxis(vt.array());
         vt *= B2[1]/C1[1];     fd1 = fd2;  fd2 = fd2.mpaxis(vt.array());  Yh.SH(0,0) += fd2.Re();
         vt *= B1[1]/B2[1];                 fd1 = fd1.mpaxis(vt.array());  Yh.SH(2,0) += fd1.Re();

         fd1 = Yin.SH(2,1);     fd1 = fd1.Dy();
         vt *= C1[2]/B1[1];     fd2 = fd1;                                 Yh.SH(3,2) += fd1.mpaxis(vt.array());
         vt *= B2[2]/C1[2];     fd1 = fd2;  fd2 = fd2.mpaxis(vt.array());  Yh.SH(1,0) += fd2.Re();
         vt *= B1[2]/B2[2];                 fd1 = fd1.mpaxis(vt.array());  Yh.SH(3,0) += fd1.Re();
       
         for (size_t l(3); l < l0; ++l){
             fd1 = Yin.SH(l,1);     fd1 = fd1.Dy();                           
             vt *= C3[l]/B1[l-1];   fd2 = fd1;                                Yh.SH(l-1,2) += fd1.mpaxis(vt.array());
             vt *= C1[l]/C3[l];     fd1 = fd2;                                Yh.SH(l+1,2) += fd2.mpaxis(vt.array());
             vt *= B2[l]/C1[l];     fd2 = fd1;  fd1 = fd1.mpaxis(vt.array()); Yh.SH(l-1,0) += fd1.Re();
             vt *= B1[l]/B2[l];                 fd2 = fd2.mpaxis(vt.array()); Yh.SH(l+1,0) += fd2.Re();
         }
         fd1 = Yin.SH(l0,1);     fd1 = fd1.Dy();                           
         vt *= C3[l0]/B1[l0-1];  fd2 = fd1;                                Yh.SH(l0-1,2) += fd1.mpaxis(vt.array());
         vt *= B2[l0]/C3[l0];                fd2 = fd2.mpaxis(vt.array()); Yh.SH(l0-1,0) += fd2.Re();
         vt *= 1.0   /B2[l0];
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         
//       m > 1, advection in y
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         for (size_t m(2); m < m0; ++m){
//           m > 1, l = m
             fd1 = Yin.SH(m,m);       fd1 = fd1.Dy();
             vt *= C4(m,m);           fd2 = fd1;     Yh.SH(m-1,m-1) += fd1.mpaxis(vt.array());
             vt *= C2(m,m)/C4(m,m);   fd1 = fd2;     Yh.SH(m+1,m-1) += fd2.mpaxis(vt.array());
             vt *= C1[m]  /C2(m,m);                  Yh.SH(m+1,m+1) += fd1.mpaxis(vt.array());

//           m > 1, l = m+1
             fd1 = Yin.SH(m+1,m);       fd1 = fd1.Dy();
             vt *= C4(m+1,m)/C1[m];     fd2 = fd1;     Yh.SH(m  ,m-1) += fd1.mpaxis(vt.array());
             if (m+1 < l0) {
                 vt *= C2(m+1,m)/C4(m+1,m);   fd1 = fd2;   Yh.SH(m+2,m-1) += fd2.mpaxis(vt.array());
                 vt *= C1[m+1]  /C2(m+1,m);                Yh.SH(m+2,m+1) += fd1.mpaxis(vt.array());

//               m > 1, 3 < l < l0
                 for (size_t l(m+2); l < l0; ++l){
                     fd1 = Yin.SH(l,m);       fd1 = fd1.Dy();
                     vt *= C4(l,m)/C1[m+1];   fd2 = fd1;     Yh.SH(l-1,m-1) += fd1.mpaxis(vt.array());
                     vt *= C2(l,m)/C4(l,m);   fd1 = fd2;     Yh.SH(l+1,m-1) += fd2.mpaxis(vt.array());
                     vt *= C3[l]  /C2(l,m);   fd2 = fd1;     Yh.SH(l-1,m+1) += fd1.mpaxis(vt.array());
                     vt *= C1[l]  /C3[l];                    Yh.SH(l+1,m+1) += fd2.mpaxis(vt.array());
                 }
                 
                 fd1 = Yin.SH(l0,m);       fd1 = fd1.Dy();
                 vt *= C4(l0,m)/C1[l0-1];  fd2 = fd1;       Yh.SH(l0-1,m-1) += fd1.mpaxis(vt.array());
                 vt *= C3[l0]/C4(l0,m);                     Yh.SH(l0-1,m+1) += fd2.mpaxis(vt.array());
                 vt *= 1.0/C3[l0];
             }
             else {
                 vt *= 1.0/C4(m+1,m);
             }
         }
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 
//       - - - - - - - - - - - - - - - - - - - - - - - - - - -
         fd1 = Yin.SH(m0,m0);       fd1 = fd1.Dy();
         vt *= C4(m0,m0);           fd2 = fd1;     Yh.SH(m0-1,m0-1) += fd1.mpaxis(vt.array());
         if (m0 < l0) {
             vt *= C2(m0,m0)/C4(m0,m0);            Yh.SH(m0+1,m0-1) += fd2.mpaxis(vt.array());
             for (size_t l(m0+1); l < l0; ++l){
                 fd1 = Yin.SH(l,m0);          fd1 = fd1.Dy();
                 vt *= C4(l,m0)/C2(l-1,m0);   fd2 = fd1;     Yh.SH(l-1,m0-1) += fd1.mpaxis(vt.array());
                 vt *= C2(l,m0)/C4(l,  m0);                  Yh.SH(l+1,m0-1) += fd2.mpaxis(vt.array());
             }
             fd1 = Yin.SH(l0,m0);           fd1 = fd1.Dy();
             vt *= C4(l0,m0)/C2(l0-1,m0);   Yh.SH(l0-1,m0-1) += fd1.mpaxis(vt.array());
         }

        return Yh;    

    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Definition for Maxwell's Equations
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    MaxwellEq::MaxwellEq(Stat& Yslope) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : Yh(Yslope),
         tmpEB(Yslope.FLD(0)), 
         idx(static_cast< complex<double> >(Inputdata::IN().inp().x.dx())),
         idy(static_cast< complex<double> >(Inputdata::IN().inp().y.dx()))
         {
         idx = -1.0/(2.0*idx); 
         idy = -1.0/(2.0*idy);  
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Stat& MaxwellEq::Faraday(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for Faraday's Law 
//--------------------------------------------------------------

//      dBx/dt += - dEz/dy  
        tmpEB          = Yin.EMF().Ez(); 
        tmpEB         *= (-1.0) * idy;
        Yh.EMF().Bx() += tmpEB.Dy();

//      dBy/dt +=   dEz/dx       
        tmpEB          = Yin.EMF().Ez(); 
        tmpEB         *= idx;
        Yh.EMF().By() += tmpEB.Dx();        

//      dBz/dt +=   dEx/dy       
        tmpEB          = Yin.EMF().Ex(); 
        tmpEB         *= idy;
        Yh.EMF().Bz() += tmpEB.Dy();    

//      dBz/dt += - dEy/dx       
        tmpEB          = Yin.EMF().Ey(); 
        tmpEB         *= (-1.0) * idx;
        Yh.EMF().Bz() += tmpEB.Dx();   

        return Yh;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Stat& MaxwellEq::Ampere(Stat& Yin){ 
//--------------------------------------------------------------
//  This is the core calculation for Ampere's Law 
//--------------------------------------------------------------

//      dEx/dt +=   dBz/dy       
        tmpEB          = Yin.EMF().Bz(); 
        tmpEB         *= idy;
        Yh.EMF().Ex() += tmpEB.Dy();

//      dEy/dt +=  - dBz/dx       
        tmpEB          = Yin.EMF().Bz(); 
        tmpEB         *= (-1.0) * idx;
        Yh.EMF().Ey() += tmpEB.Dx();        

//      dEz/dt +=  - dBx/dy       
        tmpEB          = Yin.EMF().Bx(); 
        tmpEB         *= (-1.0) * idy;
        Yh.EMF().Ez() += tmpEB.Dy();    

//      dEz/dt += dBy/dx       
        tmpEB          = Yin.EMF().By(); 
        tmpEB         *=  idx;
        Yh.EMF().Ez() += tmpEB.Dx();   

        return Yh;
    }
//--------------------------------------------------------------
//**************************************************************
