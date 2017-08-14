///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	May 18 2010
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the 
//   Collision module
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

    #ifndef DECLERATION_FOKKER_PLANCK_H
    #define DECLERATION_FOKKER_PLANCK_H


//**************************************************************
//--------------------------------------------------------------
        class EE_isotropic_conserva {
//--------------------------------------------------------------
//      Decleration of Energy/Number-Conserving Algorithm
//--------------------------------------------------------------
        private:
//          retain a reference to the slope
            valarray<double>&  fh;

//          Define the velocity axis
            valarray<double>  vr;

//          Parameters
            valarray<double>  U4, U4m1, U2, U2m1, U1, U1m1;
            valarray<double>  U3, Qn, Pn;

//          Define the integrals
            valarray<double>  J1, I2, I4;

//          Constant
            double c_kpre;
            int    NB;

//          Calculate G for low momentum cells 
            double G(const int& n, const valarray<double>& fin);

        public:
//          Constructors/Destructors
            EE_isotropic_conserva(valarray<double>& fslope); 
         
//          Explicit Advance
            valarray<double>& Getslope(const valarray<double>& fin);    
        };
//--------------------------------------------------------------
//**************************************************************

//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//--------------------------------------------------------------
        class EE_non_const {
//--------------------------------------------------------------
//      Decleration of Non-Conserving Algorithm (NOT IN USE)
//--------------------------------------------------------------
        private:
//          retain a reference to the slope
            valarray<double>&  fh;

//          Define the powers of the velocity
            valarray<double>  vr, df, ddf;

//          Parameters
            valarray<double>  U4, U4m1, U2, U2m1, U1, U1m1;
            valarray<double>  U3, Qn, Pn;
            valarray<double>  Inv_Dvh, Inv_Dvc, Inv_v, Inv_v2, Inv_v4;

//          Define the integrals
            valarray<double>  I1, I2, I4;

//          Constant
            double c_kpre, f0_c1, f0_c2, fpp0_c, f00, fpp00;
            int    NB;


        public:
//          Constructors/Destructors
            EE_non_const(valarray<double>& fslope); 
         
//          Explicit Advance
            valarray<double>& Getslope(const valarray<double>& fin);    
        };
//--------------------------------------------------------------
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

//**************************************************************
//   Decleration for the standard RK4 class for valarrays
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class RungeKutta4_ee {
//--------------------------------------------------------------
//      Decleration of the 4rth order Runge-Kutta Class
//--------------------------------------------------------------
        public:
//          Constructor
            RungeKutta4_ee(valarray<double>& fin, int tout_start);

//          Main function
            RungeKutta4_ee& advance();    
        
//          Access
            double&  tout();
            double&  time();
            size_t&  numh();
            double&  th();

        private:
//          Helper functions
            valarray<double>& F(const valarray<double>& fin);

//          Variables
	    valarray<double>  f0, f1, fh;
            valarray<double>&  f;

            EE_isotropic_conserva Collide;
            //EE_non_const Collide;

            size_t num_h;
            double h, t, Tout;
        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Definition for class that controls the explicit loop for
//   the electron-electron algorithm
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Explicit_EC_FP {
//--------------------------------------------------------------
//      Decleration of the Explicit Energy-Conserving Algorithm
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Explicit_EC_FP(Stat& Yin, int tout_start); 

//          Collisions 
            void loop(const double& tnew);  

        private:
//          Variables
            Stat& Y;

//          Algorithm
            RungeKutta4_ee      rk4_ee;

	    double              t, tout;
            int                 Nbc, szx, szy;

            valarray<double>    fc; 

        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//      Decleration of collisions with high order harmonics
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class Implicit_SC_step {
//--------------------------------------------------------------
//      Function class to facilitate taking a step 
//      for scattering with high order harmonics
//--------------------------------------------------------------
        private:

//          Define the powers of the velocity
            valarray<double>  vr;

//          Parameters
            valarray<double>  U4, U4m1, U2, U2m1, U1, U1m1;
            valarray<double>  s1;

//          Define the integrals
            valarray<double>  J1m, I2, I4, ss1;

//          Constant
            double c_kpre;

        public:
//          Constructors/Destructors
            Implicit_SC_step(); 
         
//          Calculate the coefficients
            void reset_coeff(const valarray<double>& f00, const double Dt);    

//          Explicit Advance
            void advance(valarray< complex<double> >& fin, const int el);    
        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Implicit_FLM_step {
//--------------------------------------------------------------
//      Function class to facilitate taking a step 
//      for collisions with high order harmonics
//--------------------------------------------------------------
        private:

//          Define the powers of the velocity
            valarray<double>  vr;

//          Parameters
            valarray<double>  U4, U4m1, U2, U2m1, U1, U1m1;

//          Define the integrals
            valarray<double>  J1m, I0, I2;

//          Define the derivatives of f0
            valarray<double>  df0, ddf0;

//          Variables for the matrix A, such that A*x = b
            valarray<double> Scattering_Term;
            Matrix2D<double> Alpha_Tri;
            bool if_tridiagonal;

//          Constant
            double I0_density, I2_temperature;
            double _ZLOGei, _LOGee, kpre, Dt;

        public:
//          Constructors/Destructors
            Implicit_FLM_step(); 
         
//          Calculate the coefficients
            void reset_coeff(const valarray<double>& f00, const double Delta_t);    

//          Explicit Advance
            void advance(valarray< complex<double> >& fin, const int el);    
        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Anisotropic_Collisions {
//--------------------------------------------------------------
//      Decleration of the Anisotropic Collisions
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Anisotropic_Collisions(); 

//          Advance
            Stat& f1_loop(Stat& Yin, const double Dt);  
            Stat& f1_loop1D(Stat& Yin, const double Dt);  
            Stat& flm_loop(Stat& Yin, const double Dt);  
            Stat& flm_loop1D(Stat& Yin, const double Dt);  

        private:
//          Boundary Cells
            int                         Nbc, szx, szy;

//          Number of Harmonics
            size_t                      l0, m0;

            valarray< complex<double> > fc;
            valarray< double >          f00;

            //Implicit_SC_step  implicit_step;
            Implicit_FLM_step  implicit_step;

        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Euler_Backward {
//--------------------------------------------------------------
//      Decleration of the Euler Backward Class
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Euler_Backward(Stat& Yin, int tout_start);

//          Advance
            void advance_1(const double& tnew);    
            void advance_lm(const double& tnew);    

        private:
//          Variables
            Stat& Y;
            double t;

            bool if_implicit1D;
            
//          Actions
            Anisotropic_Collisions   Flm_Coll;

        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//-------------------------------------------------------------------
     double LOGee(double ne, double Te); 
     double ZLOGei(double ne, double Te); 
//-------------------------------------------------------------------
//**************************************************************
    #endif
