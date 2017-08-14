
//	Last Modified:	Dec  1 2011
///////////////////////////////////////////////////////////

//   
//   TEMPORARY INPUT
//   This cpp file contains the definitions for the input 
///////////////////////////////////////////////////////////

//
// 
//   namespace Inputdata::
//    The data in the input deck are collected in the class 
//    "Input". For an object obj of class Input the data can 
//    be accessed as follows: 
//
//--> obj.inp().l0  :: (size_t) the l0 number of harmonics
//       "     .m0  :: (size_t) the m0 number of harmonics
//       "     .x   :: Axis<double> for the spatial "x"
//       "     .pr  :: Axis<double> for the radial momentum
//       "     .CLF :: (double) the "Courant" condition
//
//--> obj.cont().n_out  :: (size_t) the number of output steps
//       "      .tstop  :: (double) the end time of the simulation
//       "      .dt_out :: (double) the time between output 
//
//--> obj.outp().x1 :: Axis<float> for the output x1 axis
//       "      .p1 :: Axis<float> for the output p1 axis
//       "      .p2 :: Axis<float> for the output p2 axis
//       "      .p3 :: Axis<float> for the output p3 axis
//
//   
//    Access to the entire list can be obtain from:
//--> obj.list().* 
//   
//    The input is made available to the rest of the program
//    with the function IN(), so that where "obj" above the
//    one will have to substitute "Inputdata::IN()"
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <fstream>
    #include <iomanip>
    #include <cstdlib>
    #include <sstream>
    #include <string>
    #include <math.h>
    #include <float.h>

//  My libraries
    #include "matrices.h"

//  Interface
    #include "decl-input.h"



//**************************************************************
//--------------------------------------------------------------
    Inputdata::input_list::input_list(){
//--------------------------------------------------------------
//  The constructor for the input_structure
//  works as a temporary input deck
//--------------------------------------------------------------

//      Parallel
        NnodesX = 16; NnodesY = 8; 

//      spherical harmonics (assumed l0 > 3, m0 > 1, l0 >= m0)
        l0 = 16;   m0 = 3; 
        nump = 300; numx_glob = 96; numy_glob = 48;

        xmin = -23780.328;  //  -400um
        xmax =  17835.246;  //   300um
        ymin = -44588.115;  //  -750um
        ymax =  44588.115;  //  -750um
        pmax = 0.6; pmin = pmax /(2.0*nump);

//      the ratio CLF/Dp 
        clf_dp = 30;
 
//      Algorithms
        if_implicitES = false;
        if_implicit1D = false;

//      Collisions for f00
        small_dt = 0.09;
        smaller_dt = 0.005;
        NB_algorithms = 6;
        if_tridiagonal = true;
        implicit_E     = true;

//      loop control 
        RKLevel = 2; 
        n_outsteps = 960; // every 250 fs
        t_stop =  2140299.52; //  240 ps

//      Restat files
        restart_time = 0; // start reading at the "restart_time"-th restart file
        n_restarts = 480;   // generate "n_restarts" restart files

//      boundary type (default periodic) 
//      0:periodic, 1:mirror 
        bndX = 1; bndY = 0;


//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//      TWO STREAM STUFF
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
        sigma_p = 0.22; center_p = 0.5;
//        sigma_x = 0.24; center_x = xmax/2.0;
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

//      Numerics output 
        o_EHist = false; 

//      Output options
        o_Ex =   true;  o_Ey  = true;  o_Ez = false;
        o_Bx =   false; o_By  = false; o_Bz = true;
        o_x1x2 = true; o_pth = false;

//      Output stress energy tensor
        o_G  = false;
        o_Px = false;  o_PxPx = false;
        o_Py = false;  o_PxPy = false; o_PyPy = false;
        o_Pz = false;  o_PxPz = false; o_PyPz = false; o_PzPz = false;

//      Nonrelativistic output 
        o_Vx = true;   o_VxVx = false;
        o_Vy = true;   o_VxVy = false; o_VyVy = false;
                       o_VxVz = false; o_VyVz = false; o_VzVz = false;
        o_Vsq = false; o_Qx   = true;  o_Qy  = true; 

        o_Temperature = true;
        o_Pressure = true;
        o_ND = false;
        o_Nu = false;

//      Output 
        o_p1x1 = false;
        o_p1x1x2 = false; 
        o_p1p2p3 = false;

//      number of output cells for cartesian data
        nump1 = 128; nump2 = 16; nump3 = 16; 
        numpx = 1000; 

//      Only generate data
        only_output = false;

//      Executable statements
        if (((NnodesX%2))!=0 && (NnodesX !=1)) 
            std::cout<<"The number of nodes "<<NnodesX<<" is not even\n";
        numx      = numx_glob / NnodesX;
        numx     += 2*RKLevel;        
        numx_glob = NnodesX * (numx-2*RKLevel);

        if (((NnodesY%2)!=0) && (NnodesY !=1)) 
            std::cout<<"The number of nodes "<<NnodesY<<" is not even\n";
        numy      = numy_glob / NnodesY;
        numy     += 2*RKLevel;        
        numy_glob = NnodesY * (numy-2*RKLevel);

//      Electron-ion collision parameters
        lnLambda = 10.3; Zeta = 4.0;

//      Density
        density_np = 1.0e+23;// in 1/cm^3

//      Density
        denormalize_fields = true;
        denormalize_length = true;

//      Phenomenological laser parameters
        polarization_direction = 0; // x = 1, y = 2, z = 3, 0 = nothing;
        inverse_bremsstrahlung = 0; // 
        linear_Ivst = true;
        polynomial_Ivst = false;
        spot_w0 = 14862.705; cntr_r0 = 0.0; cntr_x0 = -10701.0; // -100um 
        I_0 = 1.0e+16; lambda_0 = 0.351; //intensity in W/cm^2, laser wavelength in um 
        rise_time = 178352.46; //  10 ps
        flat_time = 1783524.6; // 100 ps 
        fall_time = 178352.46; //  10 ps
        ab_eta = 1.0; ab_depth = 7134.0; // 120um 

//-----------------------------------------------------------------------
//      Density profile

        // Polynomial
        setup_poly_x = false;
        nmin_x = 0.01;  nconp = 0.01;
        start_x0 = -10000.0, rise_x = 100000.0; flat_x = 120000.0; fall_x = 130000.0;

        // Exponential
        setup_exp_x = true;
        exp_nmin = 1.0e+21;// in 1/cm^3
        exp_xmin = -800.00;// -163um. This is the distance behind 0 which is
                                // where the critical density is, this is determined by the wavelength
                                // The maximum density is density_np

        // Gaussian
        setup_gauss_dens_x = false;
        //xmin = -200.0;   xmax = 512.5;
        ampl_x = 1.0; sigma_x = 96.0; center_x = 144.0;

        // X Piecewise Linear
        setup_pcwise_lin_x = false;
        //xloc.push_back(   -132936.07  );   xdens.push_back(   0.0005  ); 
        //xloc.push_back(   -112936.07  );   xdens.push_back(  /* 0.02); */0.0005  ); 
        //xloc.push_back(   -82936.07  );   xdens.push_back(  /* 0.02);*/ 0.0005  ); 
        //xloc.push_back(   1.00  );   xdens.push_back(  1.0); // 0.0005  );
        xloc.push_back(   -23780.328 );   xdens.push_back(   1.00  ); 
        xloc.push_back(   2.00  );   xdens.push_back(   1.0  );
        xloc.push_back(   4.0  );   xdens.push_back(   1.0  );
        xloc.push_back(    6.0  );   xdens.push_back(   1.0  ); 
        xloc.push_back(    8.0  );   xdens.push_back(   1.0  ); 
        xloc.push_back(    17835.246 );   xdens.push_back(   1.00  );
        // xloc.push_back(    210000.0  );   xdens.push_back(   0.0005  );
        // xloc.push_back(    260000.0  );   xdens.push_back(   0.0005  );
        // xloc.push_back(    305753.0  );   xdens.push_back(   0.0005  );
        //------> Need to check that xloc(0) < xloc(1) < ... etc

        xmin = xloc.at(0); 
        xmax = xloc.at(xloc.size()-1);
        if (xloc.size() != xdens.size()){ 
            std::cout<<"The density profile in 'x' was not declared properly\n";
            exit(1);
        }
        
        // Smoothing level
        smooth_x = 1;

        // Y Gaussian
        setup_gauss_dens_y = false;
        ampl_y = 1.0; sigma_y = 1.5; center_y = 0.0;
        //ymin = -570.0; ymax = 570.0;

        // Y Piecewise Linear
        setup_pcwise_lin_y = true;
        yloc.push_back( -44588.115); ydens.push_back(   1.0  ); 
        yloc.push_back(   2.0  );   ydens.push_back(   1.0  );
        yloc.push_back(   4.0  );   ydens.push_back(   1.0  );
        yloc.push_back(   6.00  );  ydens.push_back(   1.0  ); 
        yloc.push_back(   8.0  );   ydens.push_back(   1.0  ); 
        yloc.push_back(  44588.115); ydens.push_back(   1.0  );
        //------> Need to check that yloc(0) < yloc(1) < ... etc

        ymin = yloc.at(0); 
        ymax = yloc.at(yloc.size()-1);
        if (yloc.size() != ydens.size()){ 
            std::cout<<"The density profile in 'y' was not declared properly\n";
            exit(1);
        }
 //-----------------------------------------------------------------------

 //-----------------------------------------------------------------------
 //     Temperature profile

        // Minimum / ambient temperature
        pt_amb = 0.0312456; 

        // Gaussian
        setup_gauss_temp_x = false;
        sigma_tx = 150.0; center_tx = 0.0; ptx = 0.04625;

        // Sinusoidal Temperature
        setup_sine_temp_x = false;
        pt_sinmax = 0.025;

        // X
        setup_pcwise_tmp_x = true;
        xtloc.push_back(   -300000.1  );   xtemp.push_back(   0.03124562);  
        xtloc.push_back(   -23780.328 );  xtemp.push_back(    0.03124562); 
        xtloc.push_back(   -9690.48 );   xtemp.push_back(  0.03124562);    
        xtloc.push_back(         0.0 );    xtemp.push_back(   0.03124562);  
        xtloc.push_back(   9512.1312 );    xtemp.push_back(   0.03124562);   
//        xtloc.push_back(    425.0  );   xtemp.push_back(   0.0313  ); 
        xtloc.push_back(   17835.246  );    xtemp.push_back(   0.03124562);  // 0.15keV
        xtloc.push_back(   400000.0  );    xtemp.push_back(   0.03124562);   // 0.15keV


        if (xtloc.size() != xtemp.size()){ 
            std::cout<<"The temperature profile in 'x' was not declared properly\n";
            exit(1);
        }

        // Gaussian
        setup_gauss_temp_y = false;
        sigma_ty = 384.0; center_ty = 144.0; pty = 192.0;
        
        // Sinusoidal Temperature
        setup_sine_temp_y = false;
        pt_sinmay = 0.025;

        // Y
        setup_pcwise_tmp_y = true;
        ytloc.push_back(  -208709.64 );   ytemp.push_back(   1.00  ); 
        ytloc.push_back(  -44588.115  );   ytemp.push_back(   1.00  ); 
        ytloc.push_back(  -2.99 );   ytemp.push_back(   1.00  );
        ytloc.push_back(   2.99 );   ytemp.push_back(   1.00  ); 
        ytloc.push_back(   44566.115);   ytemp.push_back(   1.00  ); 
        ytloc.push_back(   208709.64 );   ytemp.push_back(   1.00  );

        if (ytloc.size() != ytemp.size()){ 
            std::cout<<"The temperature profile in 'y' was not declared properly\n";
            exit(1);
        }
 //-----------------------------------------------------------------------

        return;
    }
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------

//    Definition of the constructor for the input data struct
      Inputdata::decl_input:: decl_input(input_list in) :
          l0(in.l0), m0(in.m0), 
          pr(in.nump, in.pmin, in.pmax),
          x( in.numx, in.xmin, in.xmax), 
          xglob(in.numx_glob, in.xmin, in.xmax), 
          y( in.numy, in.ymin, in.ymax), 
          yglob(in.numy_glob, in.ymin, in.ymax), 
          CLF(/*pr.dx()**/in.clf_dp) { 

          x(0) = xglob(0)-in.RKLevel*xglob.dx();
          for (int i = 1; i< x.dim(); ++i) x(i) = x(i-1) + xglob.dx();
          y(0) = yglob(0)-in.RKLevel*yglob.dx();
          for (int i = 1; i< y.dim(); ++i) y(i) = y(i-1) + yglob.dx();
          return; 
       }
           

//--------------------------------------------------------------
//    Definition of the constructor for the loop control data struct
      Inputdata::decl_control::decl_control(input_list in) :
          n_out(in.n_outsteps), tstop(in.t_stop)
          { dt_out = tstop/double(n_out); return; }
//--------------------------------------------------------------

//--------------------------------------------------------------

//    Definition of the constructor for the output data struct
      Inputdata::decl_output::decl_output(input_list in) :
 //         x1(in.numx1,static_cast<float>(in.xmin), static_cast<float>(in.xmax)), 
          px(in.numpx,static_cast<float>(in.pmax)), 
          p1(in.nump1,static_cast<float>(in.pmax)), 
          p2(in.nump2,static_cast<float>(in.pmax)), 
          p3(in.nump3,static_cast<float>(in.pmax)) 
          { return; }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//  Definition of the Input Class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------

//  Constructor
    Inputdata::Input:: Input(){
        ilist = new input_list();
        in = new decl_input(*ilist); 
        control = new decl_control(*ilist); 
        out = new decl_output(*ilist);
    }
//  Destructor
    Inputdata::Input:: ~Input(){
        delete ilist;
        delete in;
        delete control;
        delete out;
    }

//--------------------------------------------------------------
//  Access 
//--------------------------------------------------------------
    Inputdata:: input_list& Inputdata::Input:: list() const  {return (*ilist);}
    Inputdata:: decl_input& Inputdata::Input:: inp() const {return (*in);}
    Inputdata:: decl_control& Inputdata::Input:: cont() const {return (*control);}
    Inputdata:: decl_output& Inputdata::Input:: outp() const {return (*out);}

//--------------------------------------------------------------
//**************************************************************


//--------------------------------------------------------------
//**************************************************************
//--------------------------------------------------------------
    Inputdata::Input& Inputdata::IN() {
        static Inputdata::Input in;
        return in;
    }

//--------------------------------------------------------------
//**************************************************************
//--------------------------------------------------------------
