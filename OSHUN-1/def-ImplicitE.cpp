///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	First Created:	June  2nd 2011
//	Last Modified:	Aug   1st 2011
///////////////////////////////////////////////////////////

//   
//   This file contains the definitions for the 
//   classes required to calculate the effect of the implicit
//   electric field and collisions for the first order harmonic
///////////////////////////////////////////////////////////
//
//   class Implicit_Efield::
//
// 
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

//  Standard Libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <limits>

    #include <math.h>
    #include <stdio.h>
    #include <float.h>

//  My Libraries
    #include "matrices.h"

//  Declerations
    #include "decl-nmethods.h"
    #include "decl-input.h"
    #include "decl-state.h"
    #include "decl-vlasovmax.h"
    #include "decl-actions.h"
    #include "decl-RK.h"
    #include "decl-fokkerplanck.h"
    #include "decl-implicitE.h"
//**************************************************************


//**************************************************************
//**************************************************************
//*******THE CLASS FOR THE CURRENT DEFINITION ******************
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    Electric_Field_Methods::Current_xyz::Current_xyz(Stat& Yin) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : jayx(Yin.FLD(0)), 
         jayy(Yin.FLD(0)), 
         jayz(Yin.FLD(0)), 
         p3og(Inputdata::IN().inp().pr.dim(), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(0)), 
            static_cast< complex<double> >(Inputdata::IN().inp().pr(Inputdata::IN().inp().pr.dim()-1))),
         Delta_p(static_cast< complex<double> >(Inputdata::IN().inp().pr.dx())),
         small(static_cast< complex<double> >(Inputdata::IN().inp().pr(0)))
         {
             for (size_t i(0); i<p3og.dim(); ++i) { // calculate p^3/g
                p3og(i) = (p3og(i)*p3og(i))*(p3og(i)/sqrt(1.0+p3og(i)*p3og(i))); 
             }

             small *= small; small *= small; small *= Inputdata::IN().inp().pr(0); 
             small *= 0.2; small *= 1.0/(Inputdata::IN().inp().pr(1)); 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Current_xyz::J(int component) {
//--------------------------------------------------------------
//  Return the current component 1 --> x, 2 --> y, 3 --> z
//--------------------------------------------------------------

    switch (component) {
        case 1: { 
            return Jx(); 
            break;
        }
        case 2: { 
            return Jy(); 
            break;
        }
        case 3: { 
            return Jz(); 
            break;
        }
        default: {
                cout << "There is no such component for the current! \n";
                exit(1);
                break;
        }
    }
    return Jx();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Current_xyz::Jx() {
       return jayx; 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Current_xyz::Jy() {
        return jayy;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Current_xyz::Jz() {
        return jayz;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Electric_Field_Methods::
    Current_xyz::calculate_J(Stat& Yin) {
//--------------------------------------------------------------
//  Update the total current 
//--------------------------------------------------------------
        complex<double> c01(0.0,1.0);
        SHarmonic f10(Yin.SH(1,0)); 
        SHarmonic f11(Yin.SH(1,1)); 
        f10  = f10.mpaxis(p3og.array());
        f11  = f11.mpaxis(p3og.array());

//   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
//      Basic integration
        for (size_t ix(0); ix < f10.numx(); ++ix){
            for (size_t iy(0); iy < f10.numy(); ++iy){
                jayx(ix,iy)  = f10(0,ix,iy);
                jayx(ix,iy) += f10(f10.nump()-1,ix,iy);
                jayx(ix) *= 0.5;
                for (size_t ip(1); ip < f10.nump()-1; ++ip) jayx(ix,iy)+= f10(ip,ix,iy);
            }
        }
        jayx *= Delta_p;

//      Add tiny contribution from 0-p0
        for (size_t ix(0); ix < f10.numx(); ++ix){
            for (size_t iy(0); iy < f10.numy(); ++iy){
                jayx(ix,iy) += small*f10(1,ix,iy);
             }
        }

//      Calculate the final current (minus is for electrons)
        jayx *= (-4.0) * M_PI / 3.0; 
//   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -

//      Basic integration of the complex "current"
        for (size_t ix(0); ix < f11.numx(); ++ix){
            for (size_t iy(0); iy < f11.numy(); ++iy){
                jayz(ix,iy)  = f11(0,ix,iy);
                jayz(ix,iy) += f11(f11.nump()-1,ix,iy);
                jayz(ix) *= 0.5;
                for (size_t ip(1); ip < f11.nump()-1; ++ip) jayz(ix,iy)+= f11(ip,ix,iy);
            }
        }
        jayz *= Delta_p;

//      Add tiny contribution from 0-p0
        for (size_t ix(0); ix < f11.numx(); ++ix){
            for (size_t iy(0); iy < f11.numy(); ++iy){
                jayz(ix,iy) += small*f11(1,ix,iy);
             }
        }

//      Calculate the final complex current 
        jayz *= (-8.0) * M_PI / 3.0;

//      Calculate Jy
        jayy = jayz;      jayy = jayy.Re();

//      Calculate Jz
        jayz *= c01;      jayz = jayz.Re();

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Electric_Field_Methods::
    Current_xyz::calculate_J(int component, Stat& Yin) {
//--------------------------------------------------------------
//  Update the current component 1 --> x, 2 --> y, 3 --> z
//--------------------------------------------------------------

    switch (component) {
        case 1: { 
            calculate_Jx(Yin); 
            break;
        }
        case 2: { 
            calculate_Jy(Yin); 
            break;
        }
        case 3: { 
            calculate_Jz(Yin); 
            break;
        }
        default: {
                cout << "There is no such component for the current! \n";
                exit(1);
                break;
        }
    }
          
}
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Electric_Field_Methods::
    Current_xyz::calculate_Jx(Stat& Yin) {
//--------------------------------------------------------------
//  Update the current in x-direction 
//--------------------------------------------------------------
        SHarmonic f10(Yin.SH(1,0)); 
        f10  = f10.mpaxis(p3og.array());

//      Basic integration
        for (size_t ix(0); ix < f10.numx(); ++ix){
            for (size_t iy(0); iy < f10.numy(); ++iy){
                jayx(ix,iy)  = f10(0,ix,iy);
                jayx(ix,iy) += f10(f10.nump()-1,ix,iy);
                jayx(ix) *= 0.5;
                for (size_t ip(1); ip < f10.nump()-1; ++ip) jayx(ix,iy)+= f10(ip,ix,iy);
            }
        }
        jayx *= Delta_p;

//      Add tiny contribution from 0-p0
        for (size_t ix(0); ix < f10.numx(); ++ix){
            for (size_t iy(0); iy < f10.numy(); ++iy){
                jayx(ix,iy) += small*f10(1,ix,iy);
             }
        }

//      Calculate final current (minus is for electrons)
        jayx *= (-4.0) * M_PI / 3.0; 
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Electric_Field_Methods::
    Current_xyz::calculate_Jy(Stat& Yin) {
//--------------------------------------------------------------
//  Update the current in y-direction 
//--------------------------------------------------------------
        SHarmonic f11(Yin.SH(1,1)); 
        f11  = f11.mpaxis(p3og.array());

//   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -

//      Basic integration of the complex "current"
        for (size_t ix(0); ix < f11.numx(); ++ix){
            for (size_t iy(0); iy < f11.numy(); ++iy){
                jayy(ix,iy)  = f11(0,ix,iy);
                jayy(ix,iy) += f11(f11.nump()-1,ix,iy);
                jayy(ix) *= 0.5;
                for (size_t ip(1); ip < f11.nump()-1; ++ip) jayy(ix,iy)+= f11(ip,ix,iy);
            }
        }
        jayy *= Delta_p;

//      Add tiny contribution from 0-p0
        for (size_t ix(0); ix < f11.numx(); ++ix){
            for (size_t iy(0); iy < f11.numy(); ++iy){
                jayy(ix,iy) += small*f11(1,ix,iy);
             }
        }

//      Calculate the final complex current 
        jayy *= (-8.0) * M_PI / 3.0;

//      Calculate Jy
        jayy = jayy.Re();

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Electric_Field_Methods::
    Current_xyz::calculate_Jz(Stat& Yin) {
//--------------------------------------------------------------
//  Update the total current 
//--------------------------------------------------------------
        complex<double> c01(0.0,1.0);
        SHarmonic f11(Yin.SH(1,1)); 
        f11  = f11.mpaxis(p3og.array());

//   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -

//      Basic integration of the complex "current"
        for (size_t ix(0); ix < f11.numx(); ++ix){
            for (size_t iy(0); iy < f11.numy(); ++iy){
                jayz(ix,iy)  = f11(0,ix,iy);
                jayz(ix,iy) += f11(f11.nump()-1,ix,iy);
                jayz(ix) *= 0.5;
                for (size_t ip(1); ip < f11.nump()-1; ++ip) jayz(ix,iy)+= f11(ip,ix,iy);
            }
        }
        jayz *= Delta_p;

//      Add tiny contribution from 0-p0
        for (size_t ix(0); ix < f11.numx(); ++ix){
            for (size_t iy(0); iy < f11.numy(); ++iy){
                jayz(ix,iy) += small*f11(1,ix,iy);
             }
        }

//      Calculate the final complex current 
        jayz *= (-8.0) * M_PI / 3.0;

//      Calculate Jz
        jayz *= c01;      jayz = jayz.Re();

    }
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
    Electric_Field_Methods::Efield_xyz::Efield_xyz(Stat& Yin) 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : efieldx(Yin.FLD(0)), 
         efieldy(Yin.FLD(1)), 
         efieldz(Yin.FLD(2)) {
     }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Efield_xyz::E(int component) {
//--------------------------------------------------------------
//  Return the current component 1 --> x, 2 --> y, 3 --> z
//--------------------------------------------------------------

    switch (component) {
        case 1: { 
            return Ex(); 
            break;
        }
        case 2: { 
            return Ey(); 
            break;
        }
        case 3: { 
            return Ez(); 
            break;
        }
        default: {
                cout << "There is no such component for the current! \n";
                exit(1);
                break;
        }
    }
    return Ex();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Efield_xyz::Ex() {
       return efieldx; 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Efield_xyz::Ey() {
        return efieldy;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Field2D& Electric_Field_Methods::Efield_xyz::Ez() {
        return efieldz;
    }
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
//  Definition of the pure virtual destructor for the abstract
//  class Explicit_Method
//--------------------------------------------------------------
    Electric_Field_Methods::Efield_Method::
    ~Efield_Method(){}
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Electric_Field_Methods::Explicit_E_Field::
    Explicit_E_Field(Stat& Yin, int tout_start)
        : t(static_cast<double>(tout_start) *
          Inputdata::IN().cont().dt_out),         // Initialize time 
          Y(Yin) {
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Electric_Field_Methods::Explicit_E_Field::
    advance(Explicit_Method* rk, Euler_Backward* eb){ }
//--------------------------------------------------------------

//--------------------------------------------------------------
    bool Electric_Field_Methods::Explicit_E_Field::
    implicitE() const { 
        return false;
    } 
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//**************************************************************
//   Definition for the implicit electric field 
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Electric_Field_Methods::Implicit_E_Field::
    Implicit_E_Field(Stat& Yin, int tout_start)
//--------------------------------------------------------------
//  Constructor for the implicit electric field
//--------------------------------------------------------------
       : t(static_cast<double>(tout_start) *
         Inputdata::IN().cont().dt_out),         // Initialize time 
         Y(Yin),
         JN(Yin),   J0(Yin),                     // New current, Current from E = 0
         J_Ex(Yin), J_Ey(Yin), J_Ez(Yin),        // Current due to the effect of Ex, Ey, Ez
         EN(Yin),   E0(Yin),   DE(Yin),          // New E, old E, perturbation E
         f00(Y.SH(0,0)),
         f10(f00), f11(f00),
         f20(f00), f21(f00), f22(f00),
         if_implicitES(Inputdata::IN().list().if_implicitES),
         if_implicit1D(Inputdata::IN().list().if_implicit1D),
         l0(Inputdata::IN().inp().l0), 
         m0(Inputdata::IN().inp().m0),
         idx(static_cast< complex<double> >(Inputdata::IN().inp().x.dx())),
         idy(static_cast< complex<double> >(Inputdata::IN().inp().y.dx())) {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         idx = -1.0/(2.0*idx); 
         idy = -1.0/(2.0*idy);  

        using Inputdata::IN; 
        Nbc = IN().list().RKLevel;
        szx = IN().inp().x.dim(); 
        szy = IN().inp().y.dim();


    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Electric_Field_Methods::Implicit_E_Field::
    advance(Explicit_Method* rk, Euler_Backward* eb){ 
//--------------------------------------------------------------
//  Calculate the implicit electric field
//--------------------------------------------------------------

        int zeros_in_det(1);      // This counts the number of zeros in the determinant 
        int execution_attempt(0); // This counts the number of attempts to find invert the E-field

        f00 = Y.SH(0,0); 
        f10 = Y.SH(1,0); f11 = Y.SH(1,1); 
        f20 = Y.SH(2,0); f21 = Y.SH(2,1); f22 = Y.SH(2,2); 

        FindDE();                           //  Reset DE

// - - - - - - - - - - - - - - - - - - - - - -
        while ( (zeros_in_det > 0) && ( execution_attempt < 4) ) {  // Execute this loop at most twice
            zeros_in_det = 0;                                       // Count the zeros of the determinant
            ++execution_attempt;                                // Count the execusion attempts
// - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - -
                                                // Effect of E = 0 on f00, f10, f11
            (*eb).advance_1((*rk).time());      // Collisions for f10, f11
            J0.calculate_J(Y);
            Y.SH(0,0) = f00;
            Y.SH(1,0) = f10; Y.SH(1,1) = f11;
// - - - - - - - - - - - - - - - - - - - - - -

            Y.EMF().Ex() = DE.Ex();             // Ex = DEx
            (*rk).advance_IEx();                // Effect of DEx on Y00, Y10, Y11, Y20, Y21, Y22
            (*eb).advance_1((*rk).time());      // Collisions for Y10, Y11
            J_Ex.calculate_J(Y);                // Evaluate J(DEx)
            Y.SH(0,0) = f00;
            Y.SH(1,0) = f10; Y.SH(1,1) = f11;
            Y.SH(2,0) = f20; Y.SH(2,1) = f21; Y.SH(2,2) = f22;
// - - - - - - - - - - - - - - - - - - - - - -

            if ( !(if_implicit1D) ) {
                Y.EMF().Ey() = DE.Ey();             // Ey = DEy
                (*rk).advance_IEy();                // Effect of DEy on Y00, Y10, Y11, Y20, Y21, Y22
                (*eb).advance_1((*rk).time());      // Collisions for Y10, Y11
                J_Ey.calculate_J(Y);                // Evaluate J(DEy)
                Y.SH(0,0) = f00;
                Y.SH(1,0) = f10; Y.SH(1,1) = f11;
                Y.SH(2,0) = f20; Y.SH(2,1) = f21; Y.SH(2,2) = f22;
// - - - - - - - - - - - - - - - - - - - - - -

                Y.EMF().Ez() = DE.Ez();             // Ez = DEz
                (*rk).advance_IEz();                // Effect of DEz on Y00, Y10, Y11, Y20, Y21, Y22
                (*eb).advance_1((*rk).time());      // Collisions for Y10, Y11
                J_Ez.calculate_J(Y);                // Evaluate J(DEz)
                Y.SH(0,0) = f00;
                Y.SH(1,0) = f10; Y.SH(1,1) = f11;
                Y.SH(2,0) = f20; Y.SH(2,1) = f21; Y.SH(2,2) = f22;
// - - - - - - - - - - - - - - - - - - - - - -

                if ( !(if_implicitES) ) {
                    Ampere();                           // Calculate JN
                }
            }
// - - - - - - - - - - - - - - - - - - - - - -

            //  calculate new EN 

            Matrix2D< complex<double> > sgm(3,3);     
            valarray< complex<double> > clm(3);

            if ( !(if_implicit1D) ) {
                for (size_t iy(0); iy < szy; ++iy){
                    for (size_t ix(0); ix < szx; ++ix){
                        sgm(0,0) =  ( J_Ex.Jx()(ix,iy) - J0.Jx()(ix,iy) ) / DE.Ex()(ix,iy);  // sxx = dJ(E_x)_x/DE_x
                        sgm(0,1) =  ( J_Ey.Jx()(ix,iy) - J0.Jx()(ix,iy) ) / DE.Ey()(ix,iy);  // sxy = dJ(E_y)_x/DE_y
                        sgm(0,2) =  ( J_Ez.Jx()(ix,iy) - J0.Jx()(ix,iy) ) / DE.Ez()(ix,iy);  // sxz = dJ(E_z)_x/DE_z

                        sgm(1,0) =  ( J_Ex.Jy()(ix,iy) - J0.Jy()(ix,iy) ) / DE.Ex()(ix,iy);  // syx = dJ(E_x)_y/DE_x
                        sgm(1,1) =  ( J_Ey.Jy()(ix,iy) - J0.Jy()(ix,iy) ) / DE.Ey()(ix,iy);  // syy = dJ(E_y)_y/DE_y
                        sgm(1,2) =  ( J_Ez.Jy()(ix,iy) - J0.Jy()(ix,iy) ) / DE.Ez()(ix,iy);  // syz = dJ(E_z)_y/DE_z

                        sgm(2,0) =  ( J_Ex.Jz()(ix,iy) - J0.Jz()(ix,iy) ) / DE.Ex()(ix,iy);  // szx = dJ(E_x)_z/DE_x
                        sgm(2,1) =  ( J_Ey.Jz()(ix,iy) - J0.Jz()(ix,iy) ) / DE.Ey()(ix,iy);  // szy = dJ(E_y)_z/DE_y
                        sgm(2,2) =  ( J_Ez.Jz()(ix,iy) - J0.Jz()(ix,iy) ) / DE.Ez()(ix,iy);  // szz = dJ(E_z)_z/DE_z
  
                        clm[0] =  JN.Jx()(ix,iy) - J0.Jx()(ix,iy);
                        clm[1] =  JN.Jy()(ix,iy) - J0.Jy()(ix,iy);
                        clm[2] =  JN.Jz()(ix,iy) - J0.Jz()(ix,iy);
                
                    // Solve the 3 by 3 system of equations
                        complex<double> D_sgm( Det33(sgm) );            // The Determinant of the conductivity tensor
                        if ( abs( D_sgm.real() ) > 6.0*DBL_MIN ) {
                            EN.Ex()(ix,iy) = Detx33(clm, sgm) / D_sgm;
                            EN.Ey()(ix,iy) = Dety33(clm, sgm) / D_sgm;
                            EN.Ez()(ix,iy) = Detz33(clm, sgm) / D_sgm;
                        }
                        else {
                            ++zeros_in_det;                                 // Try again with a new perturbation as below
                            DE.Ex()(ix,iy) *= 2;                        // DEx is 2 times larger than before
                            DE.Ey()(ix,iy) *= 1.3;                      // DEy is 1.3 times larger  than before

                            EN.Ex()(ix,iy) = 0.0; 
                            EN.Ey()(ix,iy) = 0.0;
                            EN.Ez()(ix,iy) = 0.0;
                        }
                    }
                }
            }
            else {
                for (size_t iy(0); iy < szy; ++iy){
                    for (size_t ix(0); ix < szx; ++ix){
                        // Solution for alculation JN = 0
                        EN.Ex()(ix,iy) =  DE.Ex()(ix,iy) * J0.Jx()(ix,iy) 
                                       /( J0.Jx()(ix,iy) - J_Ex.Jx()(ix,iy) );
                        EN.Ey()(ix,iy) = 0.0;// DE.Ey()(ix,iy) * J0.Jy()(ix,iy) 
                        EN.Ez()(ix,iy) = 0.0;// DE.Ey()(ix,iy) * J0.Jy()(ix,iy) 
                        // EN.Ey()(ix,iy) =  DE.Ey()(ix,iy) * J0.Jy()(ix,iy) 
                        //               /( J0.Jy()(ix,iy) - J_Ey.Jy()(ix,iy) );
                    }
                }
            }
// - - - - - - - - - - - - - - - - - - - - - -

        } // End the while statement
// - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - -
        Y.EMF().Ex() = EN.Ex();
        Y.EMF().Ey() = EN.Ey();
        Y.EMF().Ez() = EN.Ez();

        if ( zeros_in_det > 0) {
                        cout << "WARNING, Det = 0 in "<<zeros_in_det <<"locations \n";
                        if ( zeros_in_det > 8) { exit(1);}
        }

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    bool Electric_Field_Methods::Implicit_E_Field::
    implicitE() const { 
        return true;
    } 
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Electric_Field_Methods::Implicit_E_Field::FindDE(){
//--------------------------------------------------------------
//  Reset DE
//--------------------------------------------------------------
        double Eps(16.0*numeric_limits<double>::epsilon());  
        double LargeEps(sqrt(Eps));  

        DE.Ex() = Y.EMF().Ex();
        DE.Ey() = Y.EMF().Ey();
        DE.Ez() = Y.EMF().Ez();

        // Calculate DEx = LargeEps*(|E|)+Eps
        for (size_t iy(0); iy < szy; ++iy){
            for (size_t ix(0); ix < szx; ++ix){
                DE.Ex()(ix,iy) *=  DE.Ex()(ix,iy);
                DE.Ey()(ix,iy) *=  DE.Ey()(ix,iy);
                DE.Ez()(ix,iy) *=  DE.Ez()(ix,iy);

                DE.Ex()(ix,iy) +=  DE.Ey()(ix,iy);
                DE.Ex()(ix,iy) +=  DE.Ez()(ix,iy);
                DE.Ex()(ix,iy)  =  sqrt(DE.Ex()(ix,iy)); 
                DE.Ex()(ix,iy) *=  LargeEps;
                DE.Ex()(ix,iy) +=  Eps;
            }
        }
        DE.Ey() = DE.Ex(); 
        DE.Ez() = DE.Ex(); 

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Electric_Field_Methods::Implicit_E_Field::Ampere(){
//--------------------------------------------------------------
//  Use Ampere's law to calculate JN
//--------------------------------------------------------------
        Field2D tmpJi(Y.EMF().Bz());

//      Jx =  dBz/dy       
//      tmpJi    = Yin.EMF().Bz(); (same as assignment above)
        tmpJi   *= idy;
        JN.Jx()  = tmpJi.Dy();

//      Jy = -dBz/dx       
        tmpJi    = Y.EMF().Bz(); 
        tmpJi   *= (-1.0) * idx;
        JN.Jy()  = tmpJi.Dx();        

//      Jz = -dBx/dy       
        tmpJi    = Y.EMF().Bx(); 
        tmpJi   *= (-1.0) * idy;
        JN.Jz()  = tmpJi.Dy();    

//      Jz += dBy/dx       
        tmpJi    = Y.EMF().By(); 
        tmpJi   *=  idx;
        JN.Jz() += tmpJi.Dx();   
    }
//--------------------------------------------------------------



//**************************************************************

