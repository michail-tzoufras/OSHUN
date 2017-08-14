//   Contributing authors :	Michail Tzoufras
//
//	Modified:	Apr  14 2011
//	Last Modified:	July 25 2012
///////////////////////////////////////////////////////////

//   
//   This file contains the definitions for the  
//   class PLaserSource
///////////////////////////////////////////////////////////
//
//   class PLaserSource:
//      The class PLaserSource stores a reference to the slope Yh 
//      and steps it calls Source(Yin), where Yin is the current
//      state, in order to calculate the effect of phenomenological
//      laser source on the distribution function.
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
    #include <algorithm>

//  My Libraries
    #include "matrices.h"

//  Declerations
    #include "decl-input.h"
    #include "decl-state.h"
    #include "decl-plsource.h"



//*******************************************************************
//-------------------------------------------------------------------
     double ZLOGei_v2(double ne, double Te) {
//-------------------------------------------------------------------
//   Calculate the Coulomb logarithm * Zeta for electron-ion collisions
//-------------------------------------------------------------------
//      Note: that the results here assume the distribution functions
//      are nonrelativistic, as is the case for the rest of the F-P 
//      part of the code.

        static double lnei;

//      if the density is positive
        if (ne > 0.0000001) {
            Te /= (3.0*ne);
            Te *= 511000; // Temperature in eV
            ne *= Inputdata::IN().list().density_np;

            /*debug*/ //cout << "Temperature = " << Te << ",   Density = "<< ne <<"\n";
            Te = log(Te); 
            ne = log(ne);

            lnei = 24.0 - 0.5*ne + Te;
            if (lnei > 2.0) return lnei * Inputdata::IN().list().Zeta;
        }
 
        return  2.0 * Inputdata::IN().list().Zeta;
    }
//-------------------------------------------------------------------
//*******************************************************************

//**************************************************************
//--------------------------------------------------------------
    PLaserSource::PLaserSource(Stat& Yin, int tout_start)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : t1(tout_start*Inputdata::IN().cont().dt_out), 
         Y00(Yin.SH(0,0)),
         nmin(0.00001), 
         acoeff(2.08050325*1.0e+12/(Inputdata::IN().inp().xglob.dx()*
                                    Inputdata::IN().list().density_np*
                                    Inputdata::IN().list().lambda_0)),
         linear_Ivst(Inputdata::IN().list().linear_Ivst),
         T0(1000.0 * sqrt(1.0e-14*Inputdata::IN().list().lambda_0*Inputdata::IN().list().lambda_0)),
         polynomial_Ivst(Inputdata::IN().list().polynomial_Ivst),
         pr(Inputdata::IN().inp().pr), 
         gamma(Inputdata::IN().inp().pr),
         fp4(pr),
         IL(Inputdata::IN().list().I_0 * Inputdata::IN().list().ab_eta), 
          //  (Inputdata::IN().list().ab_depth * sqrt(M_PI))),
         IL_sqrt(sqrt(IL)),
         IL_min(0.00001*IL),
         ILx(Inputdata::IN().inp().x.dim()),
         ILy(Inputdata::IN().inp().y.dim()),
         Ty(Inputdata::IN().inp().y.dim()),
         IL_map(Inputdata::IN().inp().x.dim(), Inputdata::IN().inp().y.dim()),
         rise_time(Inputdata::IN().list().rise_time),
         fall_time(Inputdata::IN().list().fall_time),
         flat_time(Inputdata::IN().list().flat_time){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//       The intensity profile is: I0(t)*e^(-2*(r/w0)^2) and absorption takes place over
//       depth L so that I0_r(t) ~ I0(t)*e^(-((x-x0)/L)^2). The idea is that if you 
//       integrate over x you get I0(t)
        
        using Inputdata::IN;

        // Calculate the intensity normalization coefficient
        double sumI(0.0);
        for (size_t i(0); i < IN().inp().xglob.dim(); ++i){
	    sumI += exp((-1.0)*pow((IN().inp().xglob(i)-IN().list().cntr_x0)/IN().list().ab_depth,2));
        }
        
        // Construct the intensity profile map
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (size_t i(0); i < IN().inp().x.dim(); ++i){
            ILx[i] = exp((-1.0)*pow((IN().inp().x(i)-IN().list().cntr_x0)/IN().list().ab_depth,2));
        }
        ILx *= 1.0/sumI;
        for (size_t j(0); j < IN().inp().y.dim(); ++j){
            double tmp2 = ( IN().inp().y(j) - IN().list().cntr_r0 ) / IN().list().spot_w0;
            ILy[j] =  exp(-2.0*tmp2*tmp2);
        }
        
        // Make the gamma axis
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (size_t ip(0); ip < pr.dim(); ++ip) gamma(ip) = sqrt(1.0+pr(ip)*pr(ip)); 


        // Calculate the base temperature profile
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Ty = T0 * sqrt(ILy*IL);
 
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void PLaserSource::Laser_Source(const double& tnew){ 
//--------------------------------------------------------------
//  Calculate the effect of the laser source
//--------------------------------------------------------------
        
        double I0_t = t_Profile(tnew);  // This is I0(t)*dt
        double dt = tnew -t1;
        double ILt = I0_t/dt;
        t1 = tnew;
        if (I0_t < 0) return;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//      Calculate n for this distribution function
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Matrix2D< double > n0(Y00.numx(),Y00.numy());
        double p0p1_sq(pr(0)*pr(0)/(pr(1)*pr(1))),
               inv_mp0p1_sq(1.0/(1.0-p0p1_sq));
        complex<double> f00;

        for (size_t iy(0); iy < Y00.numy(); ++iy){
            for (size_t ix(0); ix < Y00.numx(); ++ix){

                n0(ix,iy)  = Y00(0,ix,iy).real()*pr(0)*pr(0);
                n0(ix,iy) += Y00(Y00.nump()-1,ix,iy).real()
                              *pr(pr.dim()-1)*pr(pr.dim()-1);
                n0(ix,iy) *= 0.5;
                for (size_t ip(1); ip < Y00.nump()-1; ++ip) 
                    n0(ix,iy)+= pr(ip)*pr(ip)*Y00(ip,ix,iy).real();

                f00  =  Y00(1,ix,iy)/5.0*p0p1_sq;
                f00 += (Y00(0,ix,iy) - Y00(1,ix,iy)*p0p1_sq)*inv_mp0p1_sq
                     *(1.0/3.0-1.0/5.0*p0p1_sq); 
                f00 *= pr(0)*(pr(0)/pr.dx())*pr(0); 

                n0(ix,iy) += f00.real();
            }
        }
        n0 *= 4.0 * M_PI * pr.dx();
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//      Calculate n<gamma> for this distribution function
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Matrix2D< double > Vsq(Y00.numx(),Y00.numy());

        for (size_t iy(0); iy < Y00.numy(); ++iy){
            for (size_t ix(0); ix < Y00.numx(); ++ix){

                Vsq(ix,iy)  = Y00(0,ix,iy).real()*pow(pr(0),4)/pow(gamma(0),2);
                Vsq(ix,iy) += Y00(Y00.nump()-1,ix,iy).real()
                               * pow(pr(pr.dim()-1),4) / pow(gamma(pr.dim()-1),2);
                Vsq(ix,iy) *= 0.5;
                for (size_t ip(1); ip < Y00.nump()-1; ++ip){ 
                    Vsq(ix,iy)+= pow(pr(ip),4)/pow(gamma(ip),2)*Y00(ip,ix,iy).real();
                }
                f00  =  Y00(1,ix,iy)/7.0*p0p1_sq;
                f00 += (Y00(0,ix,iy) - Y00(1,ix,iy)*p0p1_sq)*inv_mp0p1_sq
                     *(1.0/5.0-1.0/7.0*p0p1_sq); 
                f00 *= pow(pr(0),4)*(pr(0)/pr.dx());                

                Vsq(ix,iy) += f00.real();
            }
        }
        Vsq *= 4.0 * M_PI /3.0 * pr.dx();
        
        for (size_t iy(0); iy < Y00.numy(); ++iy){
            for (size_t ix(0); ix < Y00.numx(); ++ix){
                Vsq(ix,iy) /= n0(ix,iy);
             }
        }
        // Convert to temperature in eV
        Vsq *= (3000.0/4.19)*(3000.0/4.19);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        valarray<double> TL_hot(Ty);
        TL_hot *= sqrt(ILt); 

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//      Modify the distribution function 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        for (size_t yy(0); yy < Y00.numy(); ++yy){
            for (size_t xx(0); xx < Y00.numx(); ++xx){ 
                
                double T_hot = Ty[yy] + Vsq(xx,yy);
                Make_Gaussian_f(fp4, pr, gamma, T_hot);
//                cout << "Thot_new = " << T_hot << ",     Thot_old = " << Ty[yy] + Vsq(xx,yy);

                if ((n0(xx,yy) > nmin) && (ILy[yy]*ILx[xx]*ILt*IL > IL_min) && (T_hot > Vsq(xx,yy)) ) {

                    double alpha = acoeff * IL_sqrt * dt / n0(xx,yy);
                           alpha *= ILx[xx] * sqrt(ILy[yy]*ILt);       // a *= I(x) * ( I(t)*I(y) )^(1/2)
                    // cout << alpha << "\n";
                    if ( (alpha < 0) || (alpha > 0.99) ) {
                        cout << " alpha = " << alpha << "\n";
                        exit(1);
                    }
                    for (size_t ip(0); ip < Y00.nump(); ++ip) {
                        Y00(ip,xx,yy) *= (1.0-alpha);
                        Y00(ip,xx,yy) += alpha*fp4(ip)*n0(xx,yy);
                    }
              //      double sumIx(0.0);
              //      double suma(0.0);
              //      cout <<  "x,y = " << xx << ","<< yy <<",        Thot =  " << T_hot << 
              //               ",        Ix(x) =  " << ILx[xx] << ",        alpha = " << alpha << "\n";
              //      sumIx += ILx[xx];
              //      suma += alpha;
                }
            }
            //cout <<  "Ty =  " << Ty[yy] << ",          ";
            //cout << "Total Ix = " << sumIx<<",          ";
            //cout << "Total alpha = " << suma<<"\n";
         }
        // exit(1);
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    }            
//--------------------------------------------------------------
                    
//--------------------------------------------------------------
    double PLaserSource::t_Profile(const double& t2){ 
//--------------------------------------------------------------
//  Calculate the effect of the laser profile
//--------------------------------------------------------------
        double t, CI;

        if ( t2 < rise_time ) {
            t = t2 / rise_time;
            if ( polynomial_Ivst ) {
                return (10*pow(t,3)-15*pow(t,4)+6*pow(t,5))*(t2-t1);
            }
            if ( linear_Ivst ) {
                return t * (t2-t1);
            }
            return 0;
        } else {
            if ( t2 < rise_time + flat_time) { 
                return t2-t1;
            } else {
                if ( t2 < rise_time + flat_time + fall_time) { 
                    t = (rise_time + flat_time + fall_time - t2) / fall_time;
                    if ( polynomial_Ivst ) {
                        return (10*pow(t,3)-15*pow(t,4)+6*pow(t,5))*(t2-t1);
                    }
                    if ( linear_Ivst ) {
                        return t *(t2-t1);
                    }
                    return 0;
                 } 
                 else return -1.0;
             } 
         }
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//**************************************************************
// An external electric field is introduced, the electric field 
// is sinusoidal Ex = v_os* sin(w0(t-x/c));
//**************************************************************
//**************************************************************
//--------------------------------------------------------------
    ExternalEfield::ExternalEfield(Stat& Yin, int tout_start)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : t1(tout_start*Inputdata::IN().cont().dt_out),
         omega_0(3.0e+10*2.0*M_PI/(1.0e-4*Inputdata::IN().list().lambda_0)),
         omega_p(5.64 * 1.0e+4*sqrt(Inputdata::IN().list().density_np)),
         w0overwp(omega_0/omega_p), 
         polarization(Inputdata::IN().list().polarization_direction), 
         LaserFields(Yin.EMF()),
//       -----
         linear_Ivst(Inputdata::IN().list().linear_Ivst),
         polynomial_Ivst(Inputdata::IN().list().polynomial_Ivst),
//       -----
         vos(Inputdata::IN().list().lambda_0 * sqrt(7.3e-19*Inputdata::IN().list().I_0)),
         IL_xprofile(Inputdata::IN().inp().x.dim()),
         IL_yprofile(Inputdata::IN().inp().y.dim()),
         xaxis(Inputdata::IN().inp().x),
         rise_time(Inputdata::IN().list().rise_time),
         fall_time(Inputdata::IN().list().fall_time),
         flat_time(Inputdata::IN().list().flat_time){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//       The intensity profile is: I0(t)*e^(-2*(r/w0)^2) and absorption takes place over
//       depth L so that I0_r(t) ~ I0(t)*e^(-((x-x0)/L)^2). The idea is that if you 
//       integrate over x you get I0(t)
        
        using Inputdata::IN;

        // Calculate the intensity normalization coefficient
        // double sumI(0.0);
        // for (size_t i(0); i < IN().inp().xglob.dim(); ++i){
	//    sumI += exp((-1.0)*pow((IN().inp().xglob(i)-IN().list().cntr_x0)/IN().list().ab_depth,2));
        // }
        
        // Construct the intensity profile map
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (size_t i(0); i < IN().inp().x.dim(); ++i){
            IL_xprofile[i] = exp((-1.0)*pow((IN().inp().x(i)-IN().list().cntr_x0)/IN().list().ab_depth,2));
        }
        //  IL_xprofile *= 1.0/sumI;

        // debug:  double suum(0);
        // debug:  for (size_t i(0); i < IN().inp().x.dim(); ++i){
        // debug:      suum += IL_xprofile[i];
        // debug:      cout << IL_xprofile[i] <<",     sum = " << suum << "\n";
        // debug:  }
        
        for (size_t j(0); j < IN().inp().y.dim(); ++j){
            double tmp2 = ( IN().inp().y(j) - IN().list().cntr_r0 ) / IN().list().spot_w0;
            IL_yprofile[j] =  exp(-2.0*tmp2*tmp2);
        }
        
        if ( abs(polarization) > 3 ) {
            polarization = 0;
        }

        // Calculate the base temperature profile
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //Ty = T0 * sqrt(ILy*IL);
 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    int ExternalEfield:: POL()  const {return polarization;} 
//--------------------------------------------------------------

//--------------------------------------------------------------
    void ExternalEfield::Laser_Efield(const double& tnew){ 
//--------------------------------------------------------------
//  Calculate the effect of the laser source
//--------------------------------------------------------------
        
        //cout << "omega_0 = " << omega_0 <<"\n";
        //cout << "omega_p = " << omega_p <<"\n";
        //cout << "np/nc = " << pow(omega_p / omega_0,2) <<"\n";
        // cout << "vos = " << vos <<"\n";
        // cout << "polarization = " << POL() <<"\n";

        double dt = tnew -t1;
        double Envelope_Value = t_Profile(tnew)/dt; // This is Envelope(t)*dt/dt = Envelope(t)

        // debug: cout << "Envelope = " << Envelope_Value << /*",  Carrier Signal = " << Carrier_Signal(tnew) <<*/"\n"; 

        switch (POL()) {
            case 1: {
               for (size_t i(0); i < Inputdata::IN().inp().x.dim(); ++i){
                    for (size_t j(0); j < Inputdata::IN().inp().y.dim(); ++j){
                        LaserFields.Ex()(i,j) = vos * IL_xprofile[i] * IL_yprofile[j] 
                                * Envelope_Value * Carrier_Signal(tnew-xaxis(i));
                    }
                }
                break;
            }
            case 2: {
                for (size_t i(0); i < Inputdata::IN().inp().x.dim(); ++i){
                    for (size_t j(0); j < Inputdata::IN().inp().y.dim(); ++j){
                        LaserFields.Ey()(i,j) = vos  * IL_xprofile[i] * IL_yprofile[j]
                                 * Envelope_Value * Carrier_Signal(tnew-xaxis(i));
                    }
                }
                break;
            }
            case 3: {
                for (size_t i(0); i < Inputdata::IN().inp().x.dim(); ++i){
                    for (size_t j(0); j < Inputdata::IN().inp().y.dim(); ++j){
                        LaserFields.Ez()(i,j) = vos * IL_xprofile[i] * IL_yprofile[j] 
                                * Envelope_Value * Carrier_Signal(tnew-xaxis(i));
                    }
                }
                break;
            } 
            default:
                break;
        }
        t1 = tnew;
        
    }            
//--------------------------------------------------------------
                    
//--------------------------------------------------------------
    double ExternalEfield::Carrier_Signal(const double& phase){ 
//--------------------------------------------------------------
//  The carrier signal of the laser
//--------------------------------------------------------------
        return sin(w0overwp*phase);
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    double ExternalEfield::t_Profile(const double& t2){ 
//--------------------------------------------------------------
//  Calculate the effect of the laser profile
//--------------------------------------------------------------
        double t;

        if ( t2 < rise_time ) {
            t = t2 / rise_time;
            if ( polynomial_Ivst ) {
                return (10*pow(t,3)-15*pow(t,4)+6*pow(t,5))*(t2-t1);
            }
            if ( linear_Ivst ) {
                return t * (t2-t1);
            }
            return 0;
        } else {
            if ( t2 < rise_time + flat_time) { 
                return t2-t1;
            } else {
                if ( t2 < rise_time + flat_time + fall_time) { 
                    t = (rise_time + flat_time + fall_time - t2) / fall_time;
                    if ( polynomial_Ivst ) {
                        return (10*pow(t,3)-15*pow(t,4)+6*pow(t,5))*(t2-t1);
                    }
                    if ( linear_Ivst ) {
                        return t *(t2-t1);
                    }
                    return 0;
                 } 
                 else return -1.0;
             } 
         }
    }
//--------------------------------------------------------------
//**************************************************************


//*******************************************************************
//*******************************************************************
// Inverse Bremsstrahlung operator
//*******************************************************************
//*******************************************************************

//*******************************************************************
//--------------------------------------------------------------
    IB_Operator::IB_Operator(valarray<double>& fslope)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       : fh(fslope),
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//       Velocity
         vr(0.0, Inputdata::IN().list().nump), 
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//       Constants for Integrals
         U4(0.0, Inputdata::IN().list().nump), 
         U4m1(0.0, Inputdata::IN().list().nump), 
         U2(0.0, Inputdata::IN().list().nump), 
         U2m1(0.0, Inputdata::IN().list().nump), 
//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         Inv_Uav6(0.0, Inputdata::IN().list().nump), 
         gn(0.0, Inputdata::IN().list().nump), 
         Qn(0.0, Inputdata::IN().list().nump), 
         Pn(0.0, Inputdata::IN().list().nump)
         {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

         //classical electron radius
         double re(2.8179402894e-13);           
         double kp(sqrt(4.0*M_PI*(Inputdata::IN().list().density_np)*re));
         double omega_0(3.0e+10*2.0*M_PI/(1.0e-4*Inputdata::IN().list().lambda_0));
         double omega_p(5.64 * 1.0e+4*sqrt(Inputdata::IN().list().density_np));

         vw_coeff_cube = omega_p/omega_0 * kp*re;
         Qn_coeff = kp*re/6.0;

//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         // Determine vr
         for (size_t i(0); i < Inputdata::IN().inp().pr.dim(); ++i) {
             vr[i] = Inputdata::IN().inp().pr(i);
             vr[i] = vr[i] / (sqrt(1.0+vr[i]*vr[i]));         //you can comment this line to get to the right T and n
         }

//       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         // Determine U4, U4m1, U2, U2m1, U1, U1m1
         for (size_t i(1); i < U4.size(); ++i) {
             U4[i]   = 0.5 * pow(vr[i],4)     * (vr[i]-vr[i-1]);         
             U4m1[i] = 0.5 * pow(vr[i-1],4)   * (vr[i]-vr[i-1]);         
             U2[i]   = 0.5 * vr[i]*vr[i]      * (vr[i]-vr[i-1]);         
             U2m1[i] = 0.5 * vr[i-1]*vr[i-1]  * (vr[i]-vr[i-1]);         
         }

         // Determine <vr>
         Inv_Uav6_nm1 = pow( 2.0/vr[0], 6); 
         for (size_t i(0); i < Inputdata::IN().inp().pr.dim(); ++i) {
             Inv_Uav6[i] = pow( 2.0/(vr[i+1]+vr[i]), 6); 
         }

         // Determine Qn 
         Qn[0] = 1.0 / ((vr[0]*vr[0]*vr[1])/2.0);
         for (size_t i(1); i < Qn.size()-1; ++i) {
             Qn[i] = 1.0 / (vr[i]*vr[i]*(vr[i+1]-vr[i-1])/2.0);
         }        
         // Determine Pn 
         Pnm1 = 2.0/(vr[0]*vr[0]);
         for (size_t i(0); i < Pn.size()-1; ++i) {
             Pn[i] = 1.0 / ((vr[i+1]-vr[i])/2.0*(vr[i+1]+vr[i]));
         }    

    }
//--------------------------------------------------------------

//-------------------------------------------------------------------
    valarray<double>& IB_Operator::Getslope(const valarray<double>& fin, const double vos) {
//-------------------------------------------------------------------

//-------------------------------------------------------------------
       double ZLn_ei;
       double I4, I2;
       double vw_cube, qn_c;
       double p0overp1_sq(Inputdata::IN().inp().pr(0)/ Inputdata::IN().inp().pr(1));
       p0overp1_sq *= p0overp1_sq;
/*debug*/ // cout << p0overp1_sq<<"\n";

//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       double f00((fin[0]-fin[1]*p0overp1_sq)/(1.0-p0overp1_sq));
/*debug*/ // cout << "(p0/p1)^2 =" << p0overp1_sq << ",   f[1] = " << fin[1] << ",   f[0] = " << fin[0]<< ",   f[-1] = " << f00 <<"\n";
/*debug*/ // cout << "P[-1] = " << Pnm1 << ",    Inv_U_av-1 = " << Inv_Uav6_nm1 << "\n";


//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     Evaluate the integrals in a standard manner
       I4 = 0;
       for (int n(1); n < U4.size(); ++n) {
           I4 += U4[n]*fin[n]+U4m1[n]*fin[n-1]; 
       }

       I2 = 0;
       for (int n(1); n < U2.size(); ++n) {
           I2 += U2[n]*fin[n]+U2m1[n]*fin[n-1]; 
       }

       ZLn_ei = ZLOGei_v2(4.0*M_PI*I2,4.0*M_PI*I4);
       vw_cube = vw_coeff_cube * ZLn_ei * 4.0*M_PI*I2;
//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       for (size_t i(0); i < gn.size(); ++i) { //last location initialized with "1"
            gn[i] = Inv_Uav6[i]*vw_cube*vw_cube;
            gn[i] = 1.0/(1.0+gn[i]);
/*debug*/ //  cout << "g["<< vr[i]<<"] = " << gn[i]<< "\n"; 
       } 

/*debug*/// cout << "vos = "<< vos << ",    n/np = " << 4.0*M_PI*I2<<  ",    ZLOG = " << ZLn_ei <<",     vw = "<< pow(vw_cube,1.0/3.0) << "\n"; 

       qn_c = Qn_coeff*vos*vos*ZLn_ei*4.0*M_PI*I2; // where 4.0*M_PI*I2 is ne/np

       fh[0] = qn_c * Qn[0] * (  Pn[0]   * gn[0]   * (fin[1]-fin[0]) 
                                    - Pnm1 /(1.0 + Inv_Uav6_nm1*vw_cube*vw_cube) * (fin[0]-f00));
/*debug*/ // cout << "fh["<< vr[0]<<"] = " << fh[0]<< "\n"; 
       for (size_t i(1); i < fh.size()-1; ++i) { 
            fh[i] = qn_c * Qn[i] * (  Pn[i]   * gn[i]   * (fin[i+1]-fin[i]) 
                                    - Pn[i-1] * gn[i-1] * (fin[i]-fin[i-1])); 
/*debug*/ // cout << "fh["<< vr[i]<<"] = " << fh[i]<< "\n"; 
       }
/*debug*/ // exit(1);
/*debug*/ // cout  << "-------------------------\n";
/*debug*/ // cout  << "-------------------------\n";
/*debug*/ // cout  << "-------------------------\n";

       return fh;
    }
//-------------------------------------------------------------------
//*******************************************************************


//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

//*******************************************************************
//*******************************************************************
//  Runge Kutta loop for the electron-electron collisions
//*******************************************************************
//*******************************************************************


//*******************************************************************
//-------------------------------------------------------------------
    RungeKutta4_ib::RungeKutta4_ib(valarray<double>& fin, int tout_start)
//-------------------------------------------------------------------
//  Constructor
//-------------------------------------------------------------------
       :t(static_cast<double>(tout_start) *
          Inputdata::IN().cont().dt_out),  // Initialize time = 0
        Tout(t + Inputdata::IN().cont().dt_out),                                   // The next output time
        fh(0.0,Inputdata::IN().list().nump), 
        f0(0.0,Inputdata::IN().list().nump), 
        f1(0.0,Inputdata::IN().list().nump),    // 3 local "State" Variables
        f(fin),                                 // Initialize a reference to "State" Y
        Inversebremsstrahlung(fh){              // Initialize "Actions"
    
    }
//--------------------------------------------------------------

//  Output time
    double& RungeKutta4_ib::tout() {return Tout;}

//  Output time
    size_t& RungeKutta4_ib::numh() {return num_h;}

//  Output time
    double& RungeKutta4_ib::th() {return h;}

//  Real time of the simulation (can only be modified in the RK)
    double& RungeKutta4_ib::time()  {return t;}

//  Call Advection Actions 
    valarray<double>& RungeKutta4_ib::F(const valarray<double>& fin, const double vosc) { 
        return Inversebremsstrahlung.Getslope(fin,vosc);
    }

//--------------------------------------------------------------
    RungeKutta4_ib& RungeKutta4_ib::advance(const double vosc){    
//--------------------------------------------------------------
//  Take a step using RK4
//--------------------------------------------------------------

//      Initialization
        f0 = f; f1 = f; 

//      Step 1
        F(f1,vosc);                          // fh = F(f1)
        fh *= (0.5*h);   f1 += fh;      // f1 = f1 + (h/2)*fh
        fh *= (1.0/3.0); f  += fh;      // f  = f  + (h/6)*fh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        F(f1,vosc);     f1  = f0;            // fh = F(f1)
        fh *= (0.5*h);   f1 += fh;      // f1 = f0 + (h/2)*fh
        fh *= (2.0/3.0); f  += fh;      // f  = f  + (h/3)*fh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 3
        F(f1,vosc);                          // fh = F(f1)
        fh *= h;          f0 += fh;     // f1 = f0 + h*Yh
        fh *= (1.0/3.0);  f  += fh;     // f  = f  + (h/3)*fh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
//      Step 4
        F(f0,vosc);                          // fh = F(f0)
        fh *= (h/6.0);    f += fh;      // f  = f  + (h/6)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

        t += h;
       
        return *this;
    
    }
//-------------------------------------------------------------------
//*******************************************************************





//**************************************************************
//**************************************************************
// Inverse Brehmsstrahlung heating
//**************************************************************
//**************************************************************
//--------------------------------------------------------------
    InverseBremsstrahlung::InverseBremsstrahlung(Stat& Yin, int tout_start)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       :t(static_cast<double>(tout_start) *
          Inputdata::IN().cont().dt_out),         // Initialize time = tout_start * dt_out
         tout(t + Inputdata::IN().cont().dt_out),  // Initialize tout = time + dt_out 
         fc(0.0, Inputdata::IN().list().nump),
         rk4_ib(fc, tout_start),
//       ------------------------
         InvBremsstrahlung(Inputdata::IN().list().inverse_bremsstrahlung), 
         Y(Yin),
//       ------------------------
         omega_0(3.0e+10*2.0*M_PI/(1.0e-4*Inputdata::IN().list().lambda_0)),
         omega_p(5.64 * 1.0e+4*sqrt(Inputdata::IN().list().density_np)),
         w0overwp(omega_0/omega_p), 
//       -----
         linear_Ivst(Inputdata::IN().list().linear_Ivst),
         polynomial_Ivst(Inputdata::IN().list().polynomial_Ivst),
//       -----
         U2(0.0, Inputdata::IN().list().nump), 
         U2m1(0.0, Inputdata::IN().list().nump), 
//       -----
         vos(Inputdata::IN().list().lambda_0 * sqrt(7.3e-19*Inputdata::IN().list().I_0)),
         IL_xprofile(Inputdata::IN().inp().x.dim()),
         IL_yprofile(Inputdata::IN().inp().y.dim()),
         xaxis(Inputdata::IN().inp().x),
         rise_time(Inputdata::IN().list().rise_time),
         fall_time(Inputdata::IN().list().fall_time),
         flat_time(Inputdata::IN().list().flat_time){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//       The intensity profile is: I0(t)*e^(-2*(r/w0)^2) and absorption takes place over
//       depth L so that I0_r(t) ~ I0(t)*e^(-((x-x0)/L)^2). The idea is that if you 
//       integrate over x you get I0(t)
        
        using Inputdata::IN;
        Nbc = IN().list().RKLevel;
        szx = IN().inp().x.dim();      // size of useful x axis 
        szy = IN().inp().y.dim();      // size of useful y axis

        // Calculate the intensity normalization coefficient
        // double sumI(0.0);
        // for (size_t i(0); i < IN().inp().xglob.dim(); ++i){
	//    sumI += exp((-1.0)*pow((IN().inp().xglob(i)-IN().list().cntr_x0)/IN().list().ab_depth,2));
        // }

       for (size_t i(1); i < U2.size(); ++i) {
             U2[i]   = 0.5 * IN().inp().pr(i)*IN().inp().pr(i)      * (IN().inp().pr(i)-IN().inp().pr(i-1));         
             U2m1[i] = 0.5 * IN().inp().pr(i-1)*IN().inp().pr(i-1)      * (IN().inp().pr(i)-IN().inp().pr(i-1));         
         }
 
        // Construct the intensity profile map
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        for (size_t i(0); i < IN().inp().x.dim(); ++i){
            IL_xprofile[i] = exp((-1.0)*pow((IN().inp().x(i)-IN().list().cntr_x0)/IN().list().ab_depth,2));
        }
        //  IL_xprofile *= 1.0/sumI;

        
        for (size_t j(0); j < IN().inp().y.dim(); ++j){
            double tmp2 = ( IN().inp().y(j) - IN().list().cntr_r0 ) / IN().list().spot_w0;
            IL_yprofile[j] =  exp(-2.0*tmp2*tmp2);
        }
        
        if ( abs(InvBremsstrahlung) > 1 ) {
            InvBremsstrahlung = 0;
        }

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    int InverseBremsstrahlung:: IBSOURCE()  const {return InvBremsstrahlung;} 
//--------------------------------------------------------------

//--------------------------------------------------------------
    void InverseBremsstrahlung::loop(const double& tnew){ 
//--------------------------------------------------------------
//  Calculate the effect of the laser source
//--------------------------------------------------------------
        
        //cout << "omega_0 = " << omega_0 <<"\n";
        //cout << "omega_p = " << omega_p <<"\n";
        //cout << "np/nc = " << pow(omega_p / omega_0,2) <<"\n";
        // cout << "vos = " << vos <<"\n";
        // cout << "polarization = " << POL() <<"\n";

        double dt = tnew -t;
        double Envelope_Value = t_Profile(tnew)/dt; // This is Envelope(t)*dt/dt = Envelope(t)
        /*debug*/ // cout << "vos("<<tnew<<") = " << Envelope_Value << "\n";

        // debug: cout << "Envelope = " << Envelope_Value << /*",  Carrier Signal = " << Carrier_Signal(tnew) <<*/"\n"; 

        /*switch (POL()) {
            case 1: {
               for (size_t i(0); i < Inputdata::IN().inp().x.dim(); ++i){
                    for (size_t j(0); j < Inputdata::IN().inp().y.dim(); ++j){
                        LaserFields.Ex()(i,j) = vos * IL_xprofile[i] * IL_yprofile[j] 
                                * Envelope_Value * Carrier_Signal(tnew-xaxis(i));
                    }
                }
                break;
            }
            case 2: {
                for (size_t i(0); i < Inputdata::IN().inp().x.dim(); ++i){
                    for (size_t j(0); j < Inputdata::IN().inp().y.dim(); ++j){
                        LaserFields.Ey()(i,j) = vos  * IL_xprofile[i] * IL_yprofile[j]
                                 * Envelope_Value * Carrier_Signal(tnew-xaxis(i));
                    }
                }
                break;
            }
            case 3: {
                for (size_t i(0); i < Inputdata::IN().inp().x.dim(); ++i){
                    for (size_t j(0); j < Inputdata::IN().inp().y.dim(); ++j){
                        LaserFields.Ez()(i,j) = vos * IL_xprofile[i] * IL_yprofile[j] 
                                * Envelope_Value * Carrier_Signal(tnew-xaxis(i));
                    }
                }
                break;
            } 
            default:
                break;
        }*/

        Matrix2D< double > Vos(Inputdata::IN().inp().x.dim(),Inputdata::IN().inp().y.dim());
        for (size_t iy(0); iy < Inputdata::IN().inp().y.dim(); ++iy){
            for (size_t ix(0); ix < Inputdata::IN().inp().x.dim(); ++ix){
                   Vos(ix,iy) = vos * sqrt(IL_xprofile[ix] * IL_yprofile[iy]); 
            }
        }

//      Initialization
        rk4_ib.tout()  = tnew;
        rk4_ib.numh()  = static_cast<size_t>(static_cast<int>((tnew-t)/(Inputdata::IN().list().smaller_dt)))+1; 
        rk4_ib.th()    = (tnew-t)/static_cast<double>(rk4_ib.numh());


        for (size_t iy(0); iy < szy-2*Nbc; ++iy){
            for (size_t ix(0); ix < szx-2*Nbc; ++ix){
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                // Copy data for a specific location in space to valarray
                for (size_t ip(0); ip < fc.size(); ++ip){
                    fc[ip] = (Y.SH(0,0)(ip,ix+Nbc,iy+Nbc)).real();
                }
                double I2_before(0.0);
                for (int n(1); n < U2.size(); ++n) {
                   I2_before+= U2[n]*fc[n]+U2m1[n]*fc[n-1]; 
                 }

                // Time loop: Update the valarray
                rk4_ib.time() = t;
                for (size_t h_step(0); h_step < rk4_ib.numh(); ++h_step){
                    rk4_ib.advance(Vos(ix+Nbc,iy+Nbc)); 
                    /*debug*/ // cout << "Vos("<<ix<<","<<iy<<") = "<< Vos(ix+Nbc,iy+Nbc) << "\n";
                }
                double I2_after(0.0);
                for (int n(1); n < U2.size(); ++n) {
                   I2_after += U2[n]*fc[n]+U2m1[n]*fc[n-1]; 
                 }
                fc *= I2_before/I2_after;
                // Return updated data to the harmonic
                for (int ip(0); ip < fc.size(); ++ip){
                    Y.SH(0,0)(ip,ix+Nbc,iy+Nbc) = fc[ip];
                }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            }
        }
        t = tnew;
        
    }            
//--------------------------------------------------------------

//--------------------------------------------------------------
    double InverseBremsstrahlung::t_Profile(const double& t2){ 
//--------------------------------------------------------------
//  Calculate the effect of the laser profile
//--------------------------------------------------------------
        double t;

        if ( t2 < rise_time ) {
            t = t2 / rise_time;
            if ( polynomial_Ivst ) {
                return (10*pow(t,3)-15*pow(t,4)+6*pow(t,5))*(t2-t);
            }
            if ( linear_Ivst ) {
                return t * (t2-t);
            }
            return 0;
        } else {
            if ( t2 < rise_time + flat_time) { 
                return t2-t;
            } else {
                if ( t2 < rise_time + flat_time + fall_time) { 
                    t = (rise_time + flat_time + fall_time - t2) / fall_time;
                    if ( polynomial_Ivst ) {
                        return (10*pow(t,3)-15*pow(t,4)+6*pow(t,5))*(t2-t);
                    }
                    if ( linear_Ivst ) {
                        return t *(t2-t);
                    }
                    return 0;
                 } 
                 else return -1.0;
             } 
         }
    }
//--------------------------------------------------------------
//**************************************************************







//**************************************************************
//**************************************************************
//****** NUMERICS INFORMATION **********************************
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    NumericsInfo::NumericsInfo(int tout_start)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       // : t1(tout_start*Inputdata::IN().cont().dt_out), 
       : pr(Inputdata::IN().inp().pr), 
         gamma(Inputdata::IN().inp().pr),
         mem_block(23),
         num_checks(0) {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
        for (size_t ip(0); ip < pr.dim(); ++ip) {
            gamma(ip) = sqrt(1.0+pr(ip)*pr(ip));
        } 

         
        // Convert energy density to Joules/cm^2
        double wp, cwp;
        wp  = 5.64 * 1.0e+4*sqrt(Inputdata::IN().list().density_np);
        cwp = 3.0  * 1.0e+10/wp;
        energy_denorm = (1.602 * 1.0e-19 * Inputdata::IN().list().density_np) *
                         (cwp * Inputdata::IN().inp().x.dx()) *
                         (cwp * Inputdata::IN().inp().y.dx()) ;

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void NumericsInfo::Collect_Numerics(Stat& Yin, const double& tnew, 
                                        const size_t& h_step){    
//--------------------------------------------------------------
//  Collect the numerics information
//--------------------------------------------------------------

//      This ensures that if the energy tracking diagnostic is off the
//      Energy_Info buffer will always stay empty. As a result, output will 
//      never happen.
        if  ( !(Inputdata::IN().list().o_EHist) ) return;

//      This is to specify how many distinct items there are in the list, i.e.
//      how many times the function is called in each iteration loop. At restart 
//      this is recalculated, but it is irrelevant because it will not be exported
//      as the file already exists.
        if  ( !(num_checks > 0) ){
            if  (h_step == 0) {
                --num_checks;  
            }
            else {
                num_checks *= -1;
            }
        }        

        double DeltaE =  Total_Energy(Yin.SH(0,0)) - Energy_in_box;
        Energy_Info.push_back( static_cast<float>(DeltaE) );
        Energy_in_box += DeltaE;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    bool NumericsInfo::Export_Numerics() const
        { return (Energy_Info.size() > mem_block) ; }   
//--------------------------------------------------------------

//--------------------------------------------------------------
    int  NumericsInfo::List_Size() const 
        { return abs(num_checks); }
//--------------------------------------------------------------

//--------------------------------------------------------------
    double NumericsInfo::Total_Energy(SHarmonic& Y00){ 
//--------------------------------------------------------------
//  Calculate the Energy 
//--------------------------------------------------------------
        
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//      Calculate n<gamma> for this distribution function
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Matrix2D< double > Vsq(Y00.numx(),Y00.numy());
        double p0p1_sq(pr(0)*pr(0)/(pr(1)*pr(1))),
               inv_mp0p1_sq(1.0/(1.0-p0p1_sq));
        complex<double> f00;

        for (size_t iy(0); iy < Y00.numy(); ++iy){
            for (size_t ix(0); ix < Y00.numx(); ++ix){

                Vsq(ix,iy)  = Y00(0,ix,iy).real()*pow(pr(0),4)/pow(gamma(0),2);
                Vsq(ix,iy) += Y00(Y00.nump()-1,ix,iy).real()
                               * pow(pr(pr.dim()-1),4) / pow(gamma(pr.dim()-1),2);
                Vsq(ix,iy) *= 0.5;
                for (size_t ip(1); ip < Y00.nump()-1; ++ip){ 
                    Vsq(ix,iy)+= pow(pr(ip),4)/pow(gamma(ip),2)*Y00(ip,ix,iy).real();
                }
                f00  =  Y00(1,ix,iy)/7.0*p0p1_sq;
                f00 += (Y00(0,ix,iy) - Y00(1,ix,iy)*p0p1_sq)*inv_mp0p1_sq
                     *(1.0/5.0-1.0/7.0*p0p1_sq); 
                f00 *= pow(pr(0),4)*(pr(0)/pr.dx());                

                Vsq(ix,iy) += f00.real();
            }
        }
        Vsq *= 4.0 * M_PI /3.0 * pr.dx();
        
        // Convert to temperature in eV
        Vsq *= (3000.0/4.19)*(3000.0/4.19);
 
        // OUTPUT
        double TotalEnergy(0.0);
        int    rkl(Inputdata::IN().list().RKLevel);

        for (int ix(rkl); ix < Y00.numx()-rkl; ++ix){
            for (int iy(rkl); iy < Y00.numy()-rkl; ++iy){
                TotalEnergy += Vsq(ix,iy);
            }
        }
         
        // Convert energy density to Joules/cm
        TotalEnergy *= energy_denorm;

        return TotalEnergy;
    }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//**************************************************************





//**************************************************************
//--------------------------------------------------------------
    void Make_Gaussian_f(Axis<double>& f, Axis<double>& p,
                         Axis<double>& g,      double&  Thot){
//--------------------------------------------------------------
//  Generate a Gaussian distribution with temperature T_hot 
//  normalized to 1
//  f   : This is where the distribution will be stored
//  p   : The momentum axis
//  g   : The gamma axis
//  Thot: The hot temperature in eV. At output we calculate the 
//        new Thot which may differ slightly.
//--------------------------------------------------------------

      double PTH(4.19/3000.0*sqrt(Thot));
      if (PTH > p(1)){ 
        

          // Generate the Gaussian with PTH 
          //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          f  = p; 
          f *= p; 
          f *=  -0.5/(PTH*PTH); 
          for (size_t ip(0); ip < f.dim(); ++ip) f(ip) = exp(f(ip));
          f *= 1.0 /pow(sqrt(2*M_PI)*PTH,3);
  
          // renormalize to make sure it yields n0 = 1;
          //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          double n0  = f(0)*p(0)*p(0);
                 n0 += f(f.dim()-1) * p(p.dim()-1)*p(p.dim()-1);
                 n0 *= 0.5;
                 for (size_t ip(1); ip < f.dim()-1; ++ip){ 
                      n0 += f(ip) * p(ip)*p(ip);
                 }

          double p0p1_sq(p(0)*p(0)/(p(1)*p(1))),
                 inv_mp0p1_sq(1.0/(1.0-p0p1_sq));
          double f00  =  f(1)/5.0*p0p1_sq;
                 f00 += (f(0) - f(1)*p0p1_sq)*inv_mp0p1_sq
                      *(1.0/3.0-1.0/5.0*p0p1_sq); 
                 f00 *= p(0)*(p(0)/p.dx())*p(0); 

                 n0 += f00;
                 n0 *= 4.0 * M_PI * p.dx();

                 f *= 1.0/n0;

          // Calculate vt^2 --> Th 
          //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          double vtsq = f(0) * pow(p(0),4) / pow(g(0),2);
                 vtsq  += f(f.dim()-1) * pow(p(p.dim()-1),4)/ pow(g(p.dim()-1),2);
                 vtsq  *= 0.5;                              
                 for (size_t ip(1); ip < f.dim()-1; ++ip){ 
                     vtsq += f(ip) * pow(p(ip),4) / pow(g(ip),2);
                 }
                 f00  =  f(1)/7.0*p0p1_sq;
                 f00 += (f(0) - f(1)*p0p1_sq)*inv_mp0p1_sq
                       *(1.0/5.0-1.0/7.0*p0p1_sq); 
                 f00 *= pow(p(0),4)*(p(0)/p.dx()); 
                 vtsq += f00;
  
                 vtsq *= 4.0 * M_PI / 3.0 * p.dx();
  
                 Thot = vtsq*(3000.0/4.19)*(3000.0/4.19);
             }
             else {
                 Thot = -1.0;
             }
     }       
//--------------------------------------------------------------


     

