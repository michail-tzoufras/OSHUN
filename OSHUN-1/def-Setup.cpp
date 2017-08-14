///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Modified:	May 17 2011
//	Last Modified:	Nov 28 2011
///////////////////////////////////////////////////////////

//   
//   This cpp file contains the definitions for the 
//   initialization routines.
///////////////////////////////////////////////////////////
//
// 
//   namespace Setup_Parameters::
//
//   At the moment this just creates the folders where the
//   data is to be saved;
//
//   namespace Setup_Y::
//
//   Initializes density and temperature profiles
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <fstream>
    #include <iomanip>
    #include <cstdlib>
    #include <sstream>
    #include <string>
    #include <cstring>
    #include <math.h>
    #include <ctime>
    #include <cstdlib>
    #include <float.h>

    #include <sys/stat.h>
    #include <sys/types.h>

//  My Libraries
    #include "matrices.h"

//  Declerations
    #include "decl-input.h"
    #include "decl-state.h"
    #include "decl-setup.h"


//**************************************************************
//**************************************************************
// Make the directory tree
//**************************************************************
//**************************************************************

//**************************************************************
//---------------------------------------------------------------
     int Setup_Parameters::makefolder(string _name){
//---------------------------------------------------------------
//   Create a folder
//---------------------------------------------------------------

         mode_t _permissions(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
         char*   foldername = new char [_name.size()];
         strcpy(foldername, _name.c_str());
         int   status(mkdir(foldername,_permissions));
         
         delete[] foldername;
         return status;
     }
//---------------------------------------------------------------


//---------------------------------------------------------------
     void Setup_Parameters::folders(){
//--------------------------------------------------------------
//   create the directory tree
//--------------------------------------------------------------
        using Inputdata::IN;

        if (makefolder("RESTART") != 0) cout << "Warning: Folder 'RESTART' exists\n";

        if (makefolder("OUTPUT") != 0) cout << "Warning: Folder 'OUTPUT' exists\n";

        if ( IN().list().o_EHist ) {
            if ( makefolder("OUTPUT/NUM") != 0) cout << "Warning: Folder 'OUTPUT/NUM' exists\n";
        }

        if ( IN().list().o_Ex ||  IN().list().o_Ey || IN().list().o_Ez ||
             IN().list().o_Bx ||  IN().list().o_By || IN().list().o_Bz )  {
            if (makefolder("OUTPUT/FLD") != 0) 
                cout<<"Warning: Folder 'OUTPUT/FLD' exists\n";

            if (IN().list().o_Ex) { 
                if (makefolder("OUTPUT/FLD/EX") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/FLD/EX' exists\n";
            }
            if (IN().list().o_Ey) { 
                if (makefolder("OUTPUT/FLD/EY") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/FLD/EY' exists\n";
            }
            if (IN().list().o_Ez) { 
                if (makefolder("OUTPUT/FLD/EZ") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/FLD/EZ' exists\n";
            }
            if (IN().list().o_Bx) { 
                if (makefolder("OUTPUT/FLD/BX") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/FLD/BX' exists\n";
            }
            if (IN().list().o_By) {
                if (makefolder("OUTPUT/FLD/BY") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/FLD/BY' exists\n";
            }
            if (IN().list().o_Bz) {
                if (makefolder("OUTPUT/FLD/BZ") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/FLD/BZ' exists\n";
            } 
        }

        if ( IN().list().o_x1x2 ||  IN().list().o_pth ||
             IN().list().o_G  ||
             IN().list().o_Px || IN().list().o_PxPx ||
             IN().list().o_Py || IN().list().o_PxPy || IN().list().o_PyPy ||
             IN().list().o_Pz || IN().list().o_PxPz || IN().list().o_PyPz || IN().list().o_PzPz || 
             IN().list().o_Vx || IN().list().o_VxVx ||
             IN().list().o_Vy || IN().list().o_VxVy || IN().list().o_VyVy ||
                                 IN().list().o_VxVz || IN().list().o_VyVz || IN().list().o_VzVz ||
             IN().list().o_Vsq || IN().list().o_Temperature || IN().list().o_Pressure || IN().list().o_Qx || IN().list().o_Qy || IN().list().o_p1x1 )  {

            if (makefolder("OUTPUT/MOM") != 0)    
                cout<<"Warning: Folder 'OUTPUT/MOM' exists\n";

//          Relativistic Energy - Momentum Tensor  
//          - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - 
            if (IN().list().o_x1x2) { 
                if (makefolder("OUTPUT/MOM/N") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/N' exists\n";
            }
            if (IN().list().o_pth) { 
                if (makefolder("OUTPUT/MOM/Pth") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/Pth' exists\n";
            }
            if (IN().list().o_G) { 
                if (makefolder("OUTPUT/MOM/Gam") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/Gam' exists\n";
            }
            if (IN().list().o_Px) { 
                if (makefolder("OUTPUT/MOM/Px") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/Px' exists\n";
            }
            if (IN().list().o_PxPx) { 
                if (makefolder("OUTPUT/MOM/PxPx") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/PxPx' exists\n";
            }
            if (IN().list().o_Py) { 
                if (makefolder("OUTPUT/MOM/Py") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/Py' exists\n";
            }
            if (IN().list().o_PxPy) { 
                if (makefolder("OUTPUT/MOM/PxPy") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/PxPy' exists\n";
            }
            if (IN().list().o_PyPy) { 
                if (makefolder("OUTPUT/MOM/PyPy") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/PyPy' exists\n";
            }
            if (IN().list().o_Pz) { 
                if (makefolder("OUTPUT/MOM/Pz") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/Pz' exists\n";
            }
            if (IN().list().o_PxPz) { 
                if (makefolder("OUTPUT/MOM/PxPz") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/PxPz' exists\n";
            }
            if (IN().list().o_PyPz) { 
                if (makefolder("OUTPUT/MOM/PyPz") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/PyPz' exists\n";
            }
            if (IN().list().o_PzPz) { 
                if (makefolder("OUTPUT/MOM/PzPz") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/PzPz' exists\n";
            }

//          Nonrelativistic output
//          - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - 
            if (IN().list().o_Vx) { 
                if (makefolder("OUTPUT/MOM/Vx") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/Vx' exists\n";
            }
            if (IN().list().o_VxVx) { 
                if (makefolder("OUTPUT/MOM/VxVx") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/VxVx' exists\n";
            }
            if (IN().list().o_Vy) { 
                if (makefolder("OUTPUT/MOM/Vy") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/Vy' exists\n";
            }
            if (IN().list().o_VxVy) { 
                if (makefolder("OUTPUT/MOM/VxVy") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/VxVy' exists\n";
            }
            if (IN().list().o_VyVy) { 
                if (makefolder("OUTPUT/MOM/VyVy") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/VyVy' exists\n";
            }
            if (IN().list().o_VxVz) { 
                if (makefolder("OUTPUT/MOM/VxVz") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/VxVz' exists\n";
            }
            if (IN().list().o_VyVz) { 
                if (makefolder("OUTPUT/MOM/VyVz") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/VyVz' exists\n";
            }
            if (IN().list().o_VzVz) { 
                if (makefolder("OUTPUT/MOM/VzVz") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/VzVz' exists\n";
            }
            if (IN().list().o_Vsq) { 
                if (makefolder("OUTPUT/MOM/Vsq") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/Vsq' exists\n";
            }
            if (IN().list().o_Temperature) { 
                if (makefolder("OUTPUT/MOM/T_eV") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/T_eV' exists\n";
            }
            if (IN().list().o_Pressure) { 
                if (makefolder("OUTPUT/MOM/P_Mbar") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/P_Mbar' exists\n";
            }
            if (IN().list().o_ND) { 
                if (makefolder("OUTPUT/MOM/ND") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/ND' exists\n";
            }
            if (IN().list().o_Nu) { 
                if (makefolder("OUTPUT/MOM/Nu") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/Nu' exists\n";
            }
            if (IN().list().o_Qx) { 
                if (makefolder("OUTPUT/MOM/Qx") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/Qx' exists\n";
            }
            if (IN().list().o_Qy) { 
                if (makefolder("OUTPUT/MOM/Qy") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/MOM/Qy' exists\n";
            }
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//      TWO STREAM STUFF
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//            if (IN().list().o_p1x1) { 
//                if (makefolder("OUTPUT/MOM/p1x1") != 0) 
//                    cout<<"Warning: Folder 'OUTPUT/MOM/p1x1' exists\n";
//            }
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

        }

        if ( IN().list().o_p1x1x2 ||  IN().list().o_p1p2p3 ||  IN().list().o_p1x1) {
            if (makefolder("OUTPUT/DISTR") != 0)    
                cout<<"Warning: Folder 'OUTPUT/DISTR' exists\n";

            if (IN().list().o_p1x1x2) { 
                if (makefolder("OUTPUT/DISTR/P1X1X2") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/DISTR/P1X1X2' exists\n";
            }
            if (IN().list().o_p1p2p3) { 
                if (makefolder("OUTPUT/DISTR/P1P2P3") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/DISTR/P1P2P3' exists\n";
            }
            if (IN().list().o_p1x1) { 
                if (makefolder("OUTPUT/DISTR/P1X1") != 0) 
                    cout<<"Warning: Folder 'OUTPUT/DISTR/P1X1' exists\n";
            }
        }
 

   }
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//**************************************************************
// INITIALIZATION
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
     void Setup_Y::initialize(Stat& Y){
//--------------------------------------------------------------
//   Initialization of the harmonics and the fields
//--------------------------------------------------------------
        
        using Inputdata::IN;

        complex<double> c01(0.0,-1.0);
        Y = 0.0;
        Y.SH(0,0) = 1.0;

//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

        Matrix2D<double> Temperature_map(IN().inp().x.dim(),IN().inp().y.dim());

//      Setup the density profile
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//      x - direction
//      - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^
        // Create a long axis that allows for additional cells for smoothing
        vector<double> long_x(IN().inp().x.dim()+2*IN().list().smooth_x);
        for (size_t i(0); i < IN().list().smooth_x; ++i) {           // cells left of left-guard
            long_x.at(i) = IN().inp().x(0) -                         // x(i) = x(0)-(x(smooth-i)-x(0)) = x(0)-(smooth-i)*dx
			               ( IN().inp().x(IN().list().smooth_x-i)-IN().inp().x(0) ); 
        }
        for (size_t i(0); i < IN().inp().x.dim(); ++i)      {        // regular + guard cells
            long_x.at(i+IN().list().smooth_x) = IN().inp().x(i);                  
        }
        for (size_t i(0); i < IN().list().smooth_x; ++i) {           // cells right of right-guard
            long_x.at(i+IN().list().smooth_x+IN().inp().x.dim()) =   // x(i+N) = x(N-1)+x(i+1)-x(0) = x(N-1)+(i+1)*dx
			            IN().inp().x(IN().inp().x.dim()-1) + ( IN().inp().x(i+1)-IN().inp().x(0) );
        }
//      - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^


        // Create density profile for the long axis
        vector<double> xdens_profile(IN().inp().x.dim()+2*IN().list().smooth_x);

	if ( IN().list().setup_poly_x ) {
            Polynomial_x( xdens_profile, long_x, IN().list().start_x0, 
                        IN().list().rise_x, IN().list().flat_x,
                        IN().list().fall_x, IN().list().nmin_x);
	}
        // if ( IN().list().setup_poly_x ) {
        //    Polynomial_x_har( Y.SH(0,0), IN().inp().x, IN().list().start_x0, 
        //                IN().list().rise_x, IN().list().flat_x,
        //                IN().list().fall_x, IN().list().nmin_x);
	// }
        if (IN().list().setup_pcwise_lin_x || IN().list().setup_exp_x) { 
            Piecewise_Linear_x(xdens_profile, long_x, 
                                IN().list().xloc, IN().list().xdens, 
                                IN().list().bndX);
        }
        if ( IN().list().setup_exp_x ) { 
            Exponential_x( xdens_profile, long_x, IN().list().exp_nmin, 
                        IN().list().exp_xmin, IN().list().density_np,
                        IN().list().lambda_0);
        }
        // if ( IN().list().setup_gauss_dens_x ) {
        //    Gaussian_x_har( Y.SH(0,0), IN().inp().x, IN().list().sigma_x, 
        //                IN().list().center_x, IN().list().ampl_x);
        // }
        Smooth121(xdens_profile, IN().list().smooth_x);
        for (size_t i(0); i < IN().inp().x.dim(); ++i) {
            for (GSlice_iter< complex<double> > it(Y.SH(0,0).x0(i)); it != it.end(); ++it) {
                *it *= xdens_profile.at(i+IN().list().smooth_x); 
            }
        }

        // for (size_t i(0); i < long_x.size(); ++i) {            // cells right of right-guard
        // cout << "i = "<< i << ",    x = "<< long_x.at(i) << ",     n = " << xdens_profile.at(i) << "\n";
        // long_x.at(i+IN().list().smooth_x+IN().inp().x.dim()) = IN().inp().x(IN().inp().x.dim()-1) + ( IN().inp().x(i+1)-IN().inp().x(0) );
        // }
        // exit(1);
        // Smooth121_x( Y.SH(0,0), IN().inp().x, IN().list().smooth_x, IN().list().bndX);

//      y - direction
//      - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^ - ^
        if (IN().list().setup_gauss_dens_y) 
            Gaussian_y( Y.SH(0,0), IN().inp().y, IN().list().sigma_y,
                        IN().list().center_y, IN().list().ampl_y);
        if (IN().list().setup_pcwise_lin_y) 
            Piecewise_Linear_y( Y.SH(0,0), IN().inp().y, 
                                IN().list().yloc, IN().list().ydens, IN().list().bndY);
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//      Setup the temperature profile
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Temperature_map = 1.0;

//      x - direction
        if (IN().list().setup_gauss_temp_x) 
            Gaussian_tx( Temperature_map , IN().inp().x, 
                   IN().list().sigma_tx, IN().list().center_tx, IN().list().ptx);
        if (IN().list().setup_pcwise_tmp_x) 
            Piecewise_Linear_tmpx( Temperature_map, 
                   IN().inp().x, IN().list().xtloc, IN().list().xtemp);
        if (IN().list().setup_sine_temp_x) 
            Sinusoidal_tmpx( Temperature_map, 
                   IN().inp().x, IN().list().pt_amb, IN().list().pt_sinmax, IN().list().xmax - IN().list().xmin);

//      y - direction
        if (IN().list().setup_gauss_temp_y) 
            Gaussian_ty( Temperature_map , IN().inp().y, 
                   IN().list().sigma_ty, IN().list().center_ty, IN().list().pty);
        if (IN().list().setup_pcwise_tmp_y) 
            Piecewise_Linear_tmpy( Temperature_map, 
                   IN().inp().y, IN().list().ytloc, IN().list().ytemp);
        if (IN().list().setup_sine_temp_y) 
            Sinusoidal_tmpy( Temperature_map, 
                   IN().inp().y, IN().list().pt_amb, IN().list().pt_sinmay, IN().list().ymax - IN().list().ymin); 
//      Set the ambient temperature
        for (size_t i(0); i < Temperature_map.dim(); ++i) { 
            if (Temperature_map(i) < IN().list().pt_amb) 
                Temperature_map(i) = IN().list().pt_amb;
        }

        Gaussian_p( Y.SH(0,0), IN().inp().pr, Temperature_map); 
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//      Setup the electromagnetic fields
//      - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Y.EMF() = 0.0;
        // Y.EMF().Ex() = 0.01; 

//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//      WEIBEL STUFF
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
/*        for (size_t i(0); i< IN().inp().x.dim();++i){
            for (size_t j(0); j< IN().inp().y.dim();++j){
                Y.EMF().Bx()(i,j) =  0.9 *cos(2.0*M_PI*IN().inp().y(j) / (10.0 +IN().inp().y(2)-IN().inp().y(1)))
                                           * sin(2.0*M_PI*IN().inp().x(i) / (10.0 +IN().inp().x(2)-IN().inp().x(1)));
             }
        } 
        for (size_t i(0); i< IN().inp().x.dim();++i){
            for (size_t j(0); j< IN().inp().y.dim();++j){
                Y.EMF().By()(i,j) = (-0.9)* sin(2.0*M_PI*IN().inp().y(j) / (10.0 +IN().inp().y(2)-IN().inp().y(1)))
                                           * cos(2.0*M_PI*IN().inp().x(i) / (10.0 +IN().inp().x(2)-IN().inp().x(1)));
             }
        } 


        SHarmonic f0(Y.SH(1,1));
        f0 = 0.00104*c01;
        //Gaussian_p_sinx( Y.SH(1,1), IN().inp().pr, 0.05*IN().list().sigma_p, 2*IN().list().center_p, IN().inp().x, 0.1);
        Gaussian_p_NOTsinxsiny( f0, IN().inp().pr, 0.05*IN().list().sigma_p, 2*IN().list().center_p, IN().inp().x, IN().inp().y, 0.1);

        // ---------  Pz  --------- 
        double moment;
        Axis<double>  pr(IN().inp().pr);
        Axis< double > gamma(pr), invgamma(pr);
        for (size_t ip(0); ip < pr.dim(); ++ip) gamma(ip)    = sqrt(1.0+pr(ip)*pr(ip)); 
        for (size_t ip(0); ip < pr.dim(); ++ip) invgamma(ip) = 1.0 / gamma(ip); 
        
        complex<double> p0p1_sq(pr(0)*pr(0)/(pr(1)*pr(1))),
                        inv_mp0p1_sq(1.0/(1.0-p0p1_sq)),
                        f00;


        moment = f0(0,0,0).imag() * pow(pr(0),3) * invgamma(0);
        moment += f0(f0.nump()-1,0,0).imag() * pow(pr(pr.dim()-1),3) * invgamma(pr.dim()-1);
        moment *= 0.5;
        for (size_t ip(1); ip < f0.nump()-1; ++ip) 
            moment += pow(pr(ip),3)*f0(ip,0,0).imag()*invgamma(ip);

        f00  =  f0(1,0,0);
        f00 *= pow(pr(0),5)/(pr.dx()*pr(1)*5.0); 
        moment += f00.imag();

        moment *= 8.0 * M_PI / 3.0 * pr.dx();
        // ---------------------------- 
        Y.EMF().Bx() *= moment;
        Y.EMF().By() *= moment;*/

//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><


     }       
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
    void Setup_Y::Smooth121(vector<double>& dens, size_t smooth_level){
//--------------------------------------------------------------
//  Smooth the "dens" axis using a 1-2-1 filter "smooth_level" times
//--------------------------------------------------------------

       for (size_t lvl(0); lvl < smooth_level; ++lvl){            
           double ol1(dens.at(lvl));
           double re(0.0);
           for (size_t i(lvl+1); i < dens.size()-lvl-1; ++i) {
               re = 0.25 * ( ol1 + 2.0 * dens.at(i) + dens.at(i+1) ) ; 
               ol1 = dens.at(i); 
               dens.at(i) = re;
           }
       }
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    void Setup_Y:: Gaussian_p(SHarmonic& h, const Axis<double>& p, Matrix2D<double>& pt){
//--------------------------------------------------------------
//  Place constant temperature in the harmonic
//--------------------------------------------------------------
        double tmp;
        Matrix2D<double> alpha(pt), coeff(pt);
        
        for (int j(0); j < h.numx(); ++j){
            for (int k(0); k < h.numy(); ++k){
                alpha(j,k) = 1.0/(2.0*pt(j,k)*pt(j,k));
                coeff(j,k) = 1.0/(sqrt(2.0*M_PI)*2.0*M_PI*pt(j,k)*pt(j,k)*pt(j,k)); 
            }
        }

        for (int i(0); i < p.dim(); ++i){
           tmp =  p(i)*p(i);
           for (int j(0); j < h.numx(); ++j){
               for (int k(0); k < h.numy(); ++k){
                   h(i,j,k) *= exp(-(alpha(j,k) * tmp));
                   h(i,j,k) *= coeff(j,k);
               }
           }
        }

        //for (int i(0); i < p.dim(); ++i){
        //   tmp =  p(i)*p(i);
        //   for (int j(0); j < h.numx(); ++j){
        //       for (int k(0); k < h.numy(); ++k){
        //           h(i,j,k) *= exp(-(alpha(j,k) * tmp)) + 0.6*exp(-(4000.0*(p(i) - 0.13)*(p(i)-0.13)));
        //           h(i,j,k) *= coeff(j,k);
        //       }
        //   }
        //}
     }       
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//                Setup of the Density Profile 
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    void Setup_Y:: Piecewise_Linear_x(vector<double>& smooth_xdens, 
                               vector<double>& long_x, vector<double>& xloc,
                               vector<double>& xdens, int bndX){
//--------------------------------------------------------------
//  Piecewise linear density profile in the x direction
//--------------------------------------------------------------

        for (size_t i(0); i < long_x.size(); ++i) {
            size_t t(1);
            double xx(long_x.at(i));

            if (long_x.at(i) < xloc.at(0) - FLT_EPSILON) { // left of the (entire) box
                if (bndX == 0) { // periodic
                    xx += xloc.at(xloc.size()-1) - xloc.at(0);
                } 
                else {           // mirror 
                    xx = xloc.at(0) + (xloc.at(0) - xx);
                }
            }
            if (long_x.at(i) > xloc.at(xloc.size()-1) + FLT_EPSILON) { // right of the (entire) box 
                if (bndX == 0) { // periodic 
                    xx += xloc.at(0)-xloc.at(xloc.size()-1);
                }
                else {           // mirror
                    xx =  xloc.at(xloc.size()-1) + (xloc.at(xloc.size()-1) - xx);
                }
            }
            while ((xx > xloc.at(t)) && (t < xloc.size()-1)) ++t; 
            smooth_xdens.at(i) = xdens.at(t-1) + (xdens.at(t)-xdens.at(t-1))/(xloc.at(t)-xloc.at(t-1)) 
                                  * (xx - xloc.at(t-1)); // n(x) = n(t-1) + [n(t)-n(t-1)]/[x(t)-x(t-1)] * x-x(t-1)
            //cout<< i<< "    " << xx <<  "    " << smooth_xdens.at(i) << "\n";
        }
        
     }       
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Setup_Y:: Piecewise_Linear_y(SHarmonic& h, const Axis<double>& y,
                               vector<double>& yloc, vector<double>& ydens, int bndY){
//--------------------------------------------------------------
//  Piecewise linear density profile in the y direction
//--------------------------------------------------------------
        double dens_y, yy;
        size_t t(0);
        for (size_t i(0); i < y.dim(); ++i) {
            t = 1;
            yy = y(i);

            if (y(i) < yloc.at(0)) { 
                if (bndY == 0) { 
                    yy += yloc.at(yloc.size()-1) - yloc.at(0);
                } 
                else {
                    yy = yloc.at(0) + (yloc.at(0) - yy);
                }
            }
                    
            if (y(i) > yloc.at(yloc.size()-1))  { 
                if (bndY == 0) { 
                    yy += yloc.at(0)-yloc.at(yloc.size()-1);
                }
                else {
                    yy =  yloc.at(yloc.size()-1) + (yloc.at(yloc.size()-1) - yy);
                }
            }

            if (y(i) < yloc.at(0))              yy += yloc.at(yloc.size()-1) - yloc.at(0);
            if (y(i) > yloc.at(yloc.size()-1))  yy += yloc.at(0)-yloc.at(yloc.size()-1) ;

            while ((yy > yloc.at(t)) && (t < yloc.size()-1)) ++t; 
            dens_y = ydens.at(t-1) + (ydens.at(t)-ydens.at(t-1))/(yloc.at(t)-yloc.at(t-1)) 
                                  * (yy - yloc.at(t-1));
            
            for (GSlice_iter< complex<double> > it(h.y0(i)); it != it.end(); ++it) *it *= dens_y; 
        }
        
     }       
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Setup_Y:: Polynomial_x(vector<double>& smooth_xdens, vector<double>& long_x,// const Axis<double>& x, 
                         double start_x0, double rise_x, double flat_x,
                         double fall_x, double nmin_x){
//--------------------------------------------------------------
//  Polynomial density profile in the x direction
//--------------------------------------------------------------
        size_t ic(0);
        double xt(0.0), nt(0.0);

        while ( (ic < long_x.size()) && (long_x.at(ic) <= start_x0) ){
            // for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = nmin_x;
            smooth_xdens.at(ic) = nmin_x;
            ++ic;
        }
        while ( (ic < long_x.size()) && (long_x.at(ic) <= rise_x) && (long_x.at(ic) > start_x0) ){
            xt = (long_x.at(ic) - start_x0) / (rise_x - start_x0);
            nt = 10*pow(xt,3)-15*pow(xt,4)+6*pow(xt,5);
            nt = nmin_x + (1.0-nmin_x)*nt;
            smooth_xdens.at(ic) = nt;
            // for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = nt;
            ++ic;
        }
        while ( (ic < long_x.size()) && (long_x.at(ic) <= flat_x) && (long_x.at(ic) > rise_x) ){
            // for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = 1.0;
            smooth_xdens.at(ic) = 1.0; // nmin_x;
            ++ic;
        }
        while ( (ic < long_x.size()) && (long_x.at(ic) <= fall_x) && (long_x.at(ic) > flat_x) ){
            xt = 1.0 - ((long_x.at(ic)-flat_x) / (fall_x-flat_x));
            nt = 10*pow(xt,3)-15*pow(xt,4)+6*pow(xt,5);
            nt = nmin_x + (1.0-nmin_x)*nt;
            // for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = nt;
            smooth_xdens.at(ic) = nt;
            ++ic;
        }
        while ( (ic < long_x.size()) && (long_x.at(ic) > fall_x)){
            //for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = nmin_x;
            smooth_xdens.at(ic) = nmin_x;
            ++ic;
        }
    
        

     }       
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Setup_Y:: Polynomial_x_har(SHarmonic& h, const Axis<double>& x, 
                         double start_x0, double rise_x, double flat_x,
                         double fall_x, double nmin_x){
//--------------------------------------------------------------
//  Polynomial density profile in the x direction
//--------------------------------------------------------------
        size_t ic(0);
        double xt(0.0), nt(0.0);

        while ((x(ic) <= start_x0) && (ic < x.dim())){
            for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = nmin_x;
            ++ic;
        }
        while ((x(ic) <= rise_x) && (x(ic) > start_x0) && (ic < x.dim())){
            xt = (x(ic) - start_x0) / (rise_x - start_x0);
            nt = 10*pow(xt,3)-15*pow(xt,4)+6*pow(xt,5);
            nt = nmin_x + (1.0-nmin_x)*nt;
            for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = nt;
            ++ic;
        }
        while ((x(ic) <= flat_x) && (x(ic) > rise_x) && (ic < x.dim())){
            for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = 1.0;
            ++ic;
        }
        while ((x(ic) <= fall_x) && (x(ic) > flat_x) && (ic < x.dim())){
            xt = 1.0 - ((x(ic)-flat_x) / (fall_x-flat_x));
            nt = 10*pow(xt,3)-15*pow(xt,4)+6*pow(xt,5);
            nt = nmin_x + (1.0-nmin_x)*nt;
            for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = nt;
            ++ic;
        }
        while ((x(ic) > fall_x) && (ic < x.dim())){
            for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = nmin_x;
            ++ic;
        }
    
        

     }       
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Setup_Y:: Exponential_x(vector<double>& smooth_xdens, vector<double>& long_x,
                        /* const Axis<double>& x,*/ double exp_nmin, double exp_xmin, 
                         double exp_nmax, double lambda_0){
//--------------------------------------------------------------
//  Exponential density profile in the x direction the critical
//  density is at x = 0; After the peak is reached in density we
//  use the piecewise linear profile. 
//--------------------------------------------------------------
        size_t ic(0);
        double xt(0.0), nt(exp_nmin), nc(1.0);
        nc = 2.0*M_PI/(lambda_0*0.0001)*3.0e+10; nc /= (5.64e+4); nc *= nc;
        double log_ncOnm(log(nc/exp_nmin));

        while ((long_x.at(ic) <= exp_xmin) && (ic < long_x.size()-1)){
            // for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = exp_nmin / exp_nmax; 
            smooth_xdens.at(ic) = exp_nmin / exp_nmax; 
            ++ic;
        }
        if ((ic == long_x.size()-1) && (long_x.at(ic) <= exp_xmin)) smooth_xdens.at(ic) = exp_nmin / exp_nmax;  
        xt = ( exp_xmin - long_x.at(ic)) / exp_xmin;
        nt = exp_nmin * exp( log_ncOnm * xt);
        while ((nt <= exp_nmax) && (long_x.at(ic) > exp_xmin) && (ic < long_x.size()-1)){
            // for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = nt / exp_nmax;
            smooth_xdens.at(ic) = nt / exp_nmax; 
            ++ic;
            xt = ( exp_xmin - long_x.at(ic)) / exp_xmin;
            nt = exp_nmin * exp( log_ncOnm * xt);
        }
        if ((ic == long_x.size()-1) && (long_x.at(ic) > exp_xmin) && (nt <= exp_nmax))  smooth_xdens.at(ic) = nt / exp_nmax;  

       // for (int i(0); x(i) < x.dim(); ++i) {
       //     cout << x(i) << ",       ";
             
        // Fill the rest of the box
        //while (ic < x.dim()){
        //    for (GSlice_iter< complex<double> > it(h.x0(ic)); it != it.end(); ++it) *it = 1.0;
        //    ++ic;
        //}
     }       
//--------------------------------------------------------------

//--------------------------------------------------------------
//  THE SMOOTHING FUNCTION BELOW IS DEFUNCT. SMOOTHING IS IMPLEMENTED
//  ELSEWHERE
//--------------------------------------------------------------
    void Setup_Y:: Smooth121_x(SHarmonic& h, const Axis<double>& x, 
                               size_t smooth_level, int bndX){ 
//--------------------------------------------------------------
//  Apply a 121 smoother to this harmonic in the x direction
//--------------------------------------------------------------
       
       for (size_t lvl(0); lvl < smooth_level; ++lvl){            // Make smooth_level passes with the smoother
           size_t ic(1);//Inputdata::IN().list().RKLevel+1);
           SHarmonic new_H(h);
           complex<double> prev, next;


           cout << "Return cell end = " << x.dim()-2*Inputdata::IN().list().RKLevel-2<< "\n";
           GSlice_iter< complex<double> > prev0(h.x0(x.dim()-2*Inputdata::IN().list().RKLevel-2)); 
           GSlice_iter< complex<double> > curr0(h.x0(0)); 
           GSlice_iter< complex<double> > next0(h.x0(1));
           
           GSlice_iter< complex<double> > currN(new_H.x0(0)); 

           while ( next0 != next0.end() ){ 
               *currN = (*prev0 + 2.0*(*curr0) + *next0 )* 0.25;
               ++currN; ++prev0;++curr0;++next0;
           }

           while ((ic < x.dim()-1/*-(Inputdata::IN().list().RKLevel+1)*/)){
               GSlice_iter< complex<double> > prev(h.x0(ic-1)); 
               GSlice_iter< complex<double> > curr(h.x0(ic)); 
               GSlice_iter< complex<double> > next(h.x0(ic+1));
               GSlice_iter< complex<double> > curN(new_H.x0(ic)); 
               
               while ( next != next.end() ){ 
                   *curN = (*prev + 2.0*(*curr) + *next )* 0.25;
                   ++curN; ++prev; ++curr; ++next;
               }
               ++ic;
           }

           cout << "Return cell beginning = " << 2*Inputdata::IN().list().RKLevel+1<< "\n";
           GSlice_iter< complex<double> > prevf(h.x0(x.dim()-2)); 
           GSlice_iter< complex<double> > currf(h.x0(x.dim()-1)); 
           GSlice_iter< complex<double> > nextf(h.x0(2*Inputdata::IN().list().RKLevel+1));
           
           GSlice_iter< complex<double> > curNf(new_H.x0(x.dim()-1)); 

           while ( nextf != nextf.end() ){ 
               *curNf = (*prevf + 2.0*(*currf) + *nextf )* 0.25;
               ++curNf; ++prevf;++currf;++nextf;
           }

           h = new_H;
       }         
   }       
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Setup_Y:: Gaussian_x_har(SHarmonic& h, const Axis<double>& x, 
                         double sigma, double center, double ampl){
//--------------------------------------------------------------
//  Gaussian density profile in the x direction
//--------------------------------------------------------------
        double alpha, tmp;
        alpha = 1.0/(sqrt(2.0)*sigma);
        for (int i=0; i<x.dim(); ++i){
            tmp = alpha * ( x(i) - center );
            for (int j=0; j < h.nump(); ++j)
                for (int k=0; k < h.numy();++k)
                    h(j,i,k) *= exp(-tmp*tmp);
            
        }
        h *= ampl;
     }       
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Setup_Y:: Gaussian_y(SHarmonic& h, const Axis<double>& x, 
                          double sigma, double center, double ampl){
//--------------------------------------------------------------
//  Gaussian density profile in the y direction
//--------------------------------------------------------------
        double alpha, tmp;
        alpha = 1.0/(sqrt(2.0)*sigma);
        for (int i=0; i < x.dim();++i){
            tmp = alpha * ( x(i) - center );
            for (int j=0; j<h.nump();++j)
                for (int k=0; k<h.numx();++k)
                    h(j,k,i) *= exp(-tmp*tmp);
            
        }
        h *= ampl;
     }       
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//                Setup of the Temperature Profile
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    void Setup_Y:: Piecewise_Linear_tmpx(Matrix2D<double>& Temp, 
          const Axis<double>& x, vector<double>& xtloc, vector<double>& xtemp){
//--------------------------------------------------------------
//  Piecewise linear temperature profile in the x direction
//--------------------------------------------------------------
        double tmp_x, xx;
        size_t t(0);
        for (size_t i(0); i < x.dim(); ++i) {
            t = 1;
            xx = x(i);
            if (x(i) < xtloc.at(0))               xx += xtloc.at(xtloc.size()-1) - xtloc.at(0);
            if (x(i) > xtloc.at(xtloc.size()-1))  xx += xtloc.at(0)-xtloc.at(xtloc.size()-1) ;

            while ((xx > xtloc.at(t)) && (t < xtloc.size()-1)) ++t; 
            tmp_x = xtemp.at(t-1) + (xtemp.at(t)-xtemp.at(t-1))/(xtloc.at(t)-xtloc.at(t-1)) 
                                  * (xx - xtloc.at(t-1));
            
            for (size_t j(0); j < Temp.dim2(); ++j)  Temp(i,j) *= tmp_x; 
        }
        
     }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Setup_Y:: Piecewise_Linear_tmpy(Matrix2D<double>& Temp, 
          const Axis<double>& y, vector<double>& ytloc, vector<double>& ytemp){
//--------------------------------------------------------------
//  Piecewise linear temperature profile in the y direction
//--------------------------------------------------------------
        double tmp_y, yy;
        size_t t(0);
        for (size_t j(0); j < y.dim(); ++j) {
            t = 1;
            yy = y(j);
            if (y(j) < ytloc.at(0))               yy += ytloc.at(ytloc.size()-1) - ytloc.at(0);
            if (y(j) > ytloc.at(ytloc.size()-1))  yy += ytloc.at(0)-ytloc.at(ytloc.size()-1) ;

            while ((yy > ytloc.at(t)) && (t < ytloc.size()-1)) ++t; 
            tmp_y = ytemp.at(t-1) + (ytemp.at(t)-ytemp.at(t-1))/(ytloc.at(t)-ytloc.at(t-1)) 
                                  * (yy - ytloc.at(t-1));
            
            for (size_t i(0); i < Temp.dim1(); ++i)  Temp(i,j) *= tmp_y; 
        }
        
     }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Setup_Y:: Sinusoidal_tmpx(Matrix2D<double>& Temp, 
          const Axis<double>& x, double pmin, double pmax, double box_size){
//--------------------------------------------------------------
//  Gaussian temperature profile in the x direction
//--------------------------------------------------------------
        double alpha, tmp;
        double umax(pmin);
        double Ampl, Aver(pmax);
        Ampl = 0.5*(umax*umax - Aver*Aver);
        Aver = 0.5*(umax*umax + Aver*Aver);
 
       //for (size_t i(0); i < x.dim(); ++i){ cout << x(i)<<"\n";}exit(1);
       //umax = 0.05;
       for (size_t i(0); i < x.dim(); ++i){
           for (size_t j(0); j < Temp.dim2(); ++j){
              Temp(i,j) = Aver+Ampl*sin(2.0*M_PI*(x(i))/(box_size+x(1)-x(0)));
              Temp(i,j) = sqrt(Temp(i,j));
           }
       }
       
       /* alpha = 1.0/(sqrt(2.0)*sigma_tx);
        for (size_t i(0); i < x.dim(); ++i){
            tmp = alpha * ( x(i) - center_tx );
            tmp = ptx * exp(-tmp*tmp);
            for (size_t j(0); j < Temp.dim2(); ++j) Temp(i,j) *= tmp;
        }*/
    }       
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Setup_Y:: Sinusoidal_tmpy(Matrix2D<double>& Temp, 
          const Axis<double>& y, double pmin, double pmax, double box_size){
//--------------------------------------------------------------
//  Gaussian temperature profile in the y direction
//--------------------------------------------------------------
        double alpha, tmp;
        double umax(pmin);
        double Ampl, Aver(pmax);
        Ampl = 0.5*(umax*umax - Aver*Aver);
        Aver = 0.5*(umax*umax + Aver*Aver);
 
       //for (size_t i(0); i < x.dim(); ++i){ cout << x(i)<<"\n";}exit(1);
       //umax = 0.05;
       for (size_t i(0); i < y.dim(); ++i){
           for (size_t j(0); j < Temp.dim1(); ++j){
              Temp(j,i) = Aver+Ampl*sin(2.0*M_PI*(y(i))/(box_size+y(1)-y(0)));
              Temp(j,i) = sqrt(Temp(j,i));
           }
       }
       
       /* alpha = 1.0/(sqrt(2.0)*sigma_tx);
        for (size_t i(0); i < x.dim(); ++i){
            tmp = alpha * ( x(i) - center_tx );
            tmp = ptx * exp(-tmp*tmp);
            for (size_t j(0); j < Temp.dim2(); ++j) Temp(i,j) *= tmp;
        }*/
    }       
//--------------------------------------------------------------
//--------------------------------------------------------------
    void Setup_Y:: Gaussian_tx(Matrix2D<double>& Temp, 
          const Axis<double>& x, double sigma_tx, double center_tx, double ptx){
//--------------------------------------------------------------
//  Gaussian temperature profile in the x direction
//--------------------------------------------------------------
        double alpha, tmp;
        alpha = 1.0/(sqrt(2.0)*sigma_tx);
        for (size_t i(0); i < x.dim(); ++i){
            tmp = alpha * ( x(i) - center_tx );
            tmp = ptx * exp(-tmp*tmp);
            for (size_t j(0); j < Temp.dim2(); ++j) Temp(i,j) *= tmp;
        }
    }       
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Setup_Y:: Gaussian_ty(Matrix2D<double>& Temp, 
          const Axis<double>& y, double sigma_ty, double center_ty, double pty){
//--------------------------------------------------------------
//  Gaussian temperature profile in the y direction
//--------------------------------------------------------------
        double alpha, tmp;
        alpha = 1.0/(sqrt(2.0)*sigma_ty);
        for (size_t j(0); j < y.dim(); ++j){
            tmp = alpha * ( y(j) - center_ty );
            tmp = pty * exp(-tmp*tmp);
            for (size_t i(0); i < Temp.dim1(); ++i) Temp(i,j) *= tmp;
        }
    }       
//--------------------------------------------------------------
//**************************************************************

//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//      TWO STREAM STUFF
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//--------------------------------------------------------------
    void Setup_Y:: Gaussian_p( SHarmonic& h, Axis<double>& x, double sigma, double center ){
        double tmp, pt(sqrt(sigma));
        double alpha, coeff;
        alpha = 0.5 / (pt *pt);
        coeff = 1.0/(sqrt(2.0*M_PI)*2.0*M_PI*pt*pt*pt);
       // cout << "center = " << center <<"\n";
        for (int i(0); i < x.dim(); ++i){
           tmp =  (x(i)-center)*(x(i)-center);
           for (int j(0); j < h.numx(); ++j){
               for (int k(0); k < h.numy(); ++k){
                   h(i,j,k) *= exp(-(alpha * tmp));
                   h(i,j,k) *= coeff;
               }
           }
        }
     }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Setup_Y:: Gaussian_p_sinx( SHarmonic& h, Axis<double>& p, double sigma, double center,
                                                  Axis<double>& x, double ampl ){
        double tmp, pt(sqrt(sigma));
        double alpha, coeff;
        alpha = 0.5 / (pt *pt);
        coeff = 1.0/(sqrt(2.0*M_PI)*2.0*M_PI*pt*pt*pt);
       // cout << "center = " << center <<"\n";
        for (int i(0); i < p.dim(); ++i){
           tmp =  (p(i)-center)*(p(i)-center);
           for (int j(0); j < h.numx(); ++j){
               for (int k(0); k < h.numy(); ++k){
                   h(i,j,k) *= exp(-(alpha * tmp));
                   h(i,j,k) *= ampl*sin(0.2*M_PI*x(j))*coeff;
               }
           }
        }
     }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Setup_Y:: Constant_p( SHarmonic& h, Axis<double>& x, double value, double center ){
        for (int i(0); x(i) < center; ++i){
           for (int j(0); j < h.numx(); ++j){
               for (int k(0); k < h.numy(); ++k){
                   h(i,j,k) += value; 
               }
           }
        }
     }
//--------------------------------------------------------------
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

//--------------------------------------------------------------
//**************************************************************

//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//     WEIBEL STUFF 
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

//--------------------------------------------------------------
    void Setup_Y:: Gaussian_p_sinxsiny( SHarmonic& h, Axis<double>& p, double sigma, double center,
                                                  Axis<double>& x, Axis<double>& y, double ampl ){
        double tmp, pt(sqrt(sigma));
        double alpha, coeff;
        alpha = 0.5 / (pt *pt);
        coeff = 1.0/(sqrt(2.0*M_PI)*2.0*M_PI*pt*pt*pt);
       // cout << "center = " << center <<"\n";
        for (int i(0); i < p.dim(); ++i){
           tmp =  (p(i)-center)*(p(i)-center);
           for (int j(0); j < h.numx(); ++j){
               for (int k(0); k < h.numy(); ++k){
                   h(i,j,k) *= exp(-(alpha * tmp));
                   h(i,j,k) *= ampl*sin(2.0*M_PI*x(j)/(10.0+x(2)-x(1)))*sin(2.0*M_PI*y(k)/(10.0+y(2)-y(1)))*coeff;
               }
           }
        }
     }
    void Setup_Y:: Gaussian_p_NOTsinxsiny( SHarmonic& h, Axis<double>& p, double sigma, double center,
                                                  Axis<double>& x, Axis<double>& y, double ampl ){
        double tmp, pt(sqrt(sigma));
        double alpha, coeff;
        alpha = 0.5 / (pt *pt);
        coeff = 1.0/(sqrt(2.0*M_PI)*2.0*M_PI*pt*pt*pt);
       // cout << "center = " << center <<"\n";
        for (int i(0); i < p.dim(); ++i){
           tmp =  (p(i)-center)*(p(i)-center);
           for (int j(0); j < h.numx(); ++j){
               for (int k(0); k < h.numy(); ++k){
                   h(i,j,k) *= exp(-(alpha * tmp));
                   h(i,j,k) *= ampl*coeff;
               }
           }
        }
     }
//--------------------------------------------------------------
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
