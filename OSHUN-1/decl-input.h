///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Modified:	May 14th 2011
//	Last Modified:	Nov 28th 2011
///////////////////////////////////////////////////////////

//   
//   This header contains the declerations for the 
//   temporary input 
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

    #ifndef DECLERATION_INPUT_H
    #define DECLERATION_INPUT_H

//**************************************************************
//**************************************************************
//   Declerations for the Input namespace
//**************************************************************
//**************************************************************

//**************************************************************
    namespace Inputdata {  
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        struct input_list {
//--------------------------------------------------------------
//      decleration of input list (input deck)
//--------------------------------------------------------------
            size_t NnodesX, NnodesY;
            size_t l0, m0, nump, numx_glob, numy_glob, numx, numy;
            double xmin, xmax, ymin, ymax;
            double pmin, pmax;
            double clf_dp;
            double small_dt;
            double smaller_dt;
            int NB_algorithms;

//          Algorithms
            bool if_implicitES;
            bool if_implicit1D;
            bool if_tridiagonal;
            bool implicit_E;
            int RKLevel;
            int bndX, bndY;
            size_t n_outsteps; double t_stop;
            int restart_time;  int n_restarts;

//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//      TWO STREAM STUFF
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
            double sigma_p, center_p;
//            double sigma_x, center_x;
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

//          Output
            bool o_EHist;
            bool o_Ex, o_Ey, o_Ez, o_Bx, o_By, o_Bz, o_x1x2, o_pth, o_p1x1x2, o_p1p2p3;
            bool o_G, o_Px, o_PxPx, o_Py, o_PxPy, o_PyPy, o_Pz, o_PxPz, o_PyPz, o_PzPz;
            bool o_Vx, o_VxVx, o_Vy, o_VxVy, o_VyVy, o_VxVz, o_VyVz, o_VzVz;
            bool o_Vsq, o_Qx, o_Qy;
            bool o_Pressure, o_Temperature, o_ND, o_Nu, o_p1x1;
            bool only_output;
            size_t numpx, nump1, nump2, nump3;

//          Electron-ion collisions
            double lnLambda, Zeta, density_np;

//          Phenomenological laser parameters
            int polarization_direction;
            int inverse_bremsstrahlung;
            bool linear_Ivst, polynomial_Ivst;
            double spot_w0, cntr_r0, cntr_x0;
            double I_0, lambda_0;
            double rise_time, fall_time, flat_time;
            double ab_eta, ab_depth;

//          density initialization
            bool setup_poly_x, setup_exp_x, setup_gauss_dens_x, setup_gauss_dens_y;
            
            bool setup_pcwise_lin_x, setup_pcwise_lin_y;
            double nmin_x, nconp, start_x0, rise_x, flat_x, fall_x;
            double exp_nmin, exp_xmin;
            double ampl_x, sigma_x, center_x;
            double ampl_y, sigma_y, center_y;
            vector<double> xloc,  yloc,  xdens, ydens;
            size_t smooth_x;

//          temperature initialization
            bool setup_gauss_temp_x, setup_gauss_temp_y;
            bool setup_pcwise_tmp_x, setup_pcwise_tmp_y;
            bool setup_sine_temp_x, setup_sine_temp_y;
            double pt_amb;
            double pt_sinmax, pt_sinmay;
            double ptx, sigma_tx, center_tx;
            double pty, sigma_ty, center_ty;
            vector<double> xtloc, ytloc, xtemp, ytemp;

//          Denormalization option
            bool denormalize_fields, denormalize_length;

//          Struct constructor
            input_list( );

        };
//--------------------------------------------------------------

//--------------------------------------------------------------
        struct decl_input {
//--------------------------------------------------------------
//      decleration input information for "State"
//--------------------------------------------------------------
            size_t l0, m0;
            Axis<double>  pr;
            Axis<double>  x, xglob;
            Axis<double>  y, yglob;
            double CLF;

//          Struct constructor
            decl_input(input_list in);
        };
//--------------------------------------------------------------

//--------------------------------------------------------------
        struct decl_control {
//--------------------------------------------------------------
//      decleration of control information for the loop
//--------------------------------------------------------------
            size_t n_out;
            double tstop, dt_out;

//          Struct constructor
            decl_control(input_list in);
        };
//--------------------------------------------------------------

//--------------------------------------------------------------
        struct decl_output {
//--------------------------------------------------------------
//      decleration of output axis
//--------------------------------------------------------------
            Axis<float>  px, p1, p2, p3;//x1, 
            decl_output(input_list in);
        };
//--------------------------------------------------------------

//--------------------------------------------------------------
        class Input{
//--------------------------------------------------------------
//      Collection of fields and distribution function decleration
//--------------------------------------------------------------
        private:
            input_list *ilist;
            decl_input *in;
            decl_control *control;
            decl_output *out;

        public:
//          Constructors/Destructors
            Input(); ~Input();

//          Access
            input_list& list() const;
            decl_input& inp() const;
            decl_control& cont() const;
            decl_output& outp() const;
        };
//--------------------------------------------------------------

//--------------------------------------------------------------
       Input& IN();
//--------------------------------------------------------------

//**************************************************************

   }    

    #endif
