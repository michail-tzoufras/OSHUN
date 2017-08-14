///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Nov 29rd 2011
//	Last Modified:	Aug 13   2012
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the  
//   class PLaserSource
///////////////////////////////////////////////////////////
//
//   class PLaserSource:
//      The class PLaserSource stores a reference to the slope Yh 
//      and steps it calls Source(Yin), where Yin is the current
//      state, in order to calculate the effect of phenomenological
//      laser source on the distribution function.
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

    #ifndef DECLERATION_PLASERSOURCE_H
    #define DECLERATION_PLASERSOURCE_H

//**************************************************************
//**************************************************************
//   Decleration of the Electric Field 
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class PLaserSource {
//--------------------------------------------------------------
//      Decleration of the Phenomenological Laser Source 
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            PLaserSource(Stat& Yin, int tout_start); 
         
//          Explicit Advance
            void Laser_Source(const double& tnew);    

//          Laser Temporal Evolution
            double t_Profile(const double& t2); 

//          Laser Temporal Evolution
            bool PLSOURCE() const { return (linear_Ivst || polynomial_Ivst);} 

//          Energy in box calculation
//            double Energy_h(const double& tnew); 

        private:
            SHarmonic& Y00;
            double t1;

            double Th;
            const double nmin, IL, IL_sqrt, IL_min, T0, acoeff;
            const double rise_time, fall_time, flat_time;

            Matrix2D< double > IL_map;
            valarray<double>     ILx, ILy, Ty;
            Axis< double >  pr, gamma, fp4;

            const bool linear_Ivst, polynomial_Ivst;
        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class ExternalEfield{
//--------------------------------------------------------------
//      Decleration of the Phenomenological Laser Source 
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            ExternalEfield(Stat& Yin, int tout_start); 
         
//          Explicit Advance
            void Laser_Efield(const double& tnew);    

//          Laser Temporal Evolution
            double t_Profile(const double& t2); 
            double Carrier_Signal(const double& phase); 

//          Laser polarization
            int POL() const; 

        private:
            EMF2D& LaserFields;
            double t1;

            int polarization;

            valarray<double>     IL_xprofile, IL_yprofile;

            Axis< double >  xaxis;

            const double rise_time, fall_time, flat_time;
            const double vos, omega_0, omega_p, w0overwp;
            const bool   linear_Ivst, polynomial_Ivst;
        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class IB_Operator{
//--------------------------------------------------------------
//      Decleration of Energy/Number-Conserving Algorithm
//--------------------------------------------------------------
        private:
//          retain a reference to the slope
            valarray<double>&  fh;

//          Define the velocity axis
            valarray<double>  vr;

//          Parameters
            valarray<double>  U4, U4m1, U2, U2m1;//, U1, U1m1;
            valarray<double>  Inv_Uav6, gn, Qn, Pn;

//          Define the integrals
            // valarray<double>  I2, I4; // J1, 

//          Constant
            double c_kpre, Qn_coeff, vw_coeff_cube;
            double Inv_Uav6_nm1, Pnm1;

        public:
//          Constructors/Destructors
            IB_Operator(valarray<double>& fslope); 
         
//          Explicit Advance
            valarray<double>& Getslope(const valarray<double>& fin, const double vos);    
        };
//--------------------------------------------------------------
//**************************************************************

//--------------------------------------------------------------
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

//**************************************************************
//   Decleration for the standard RK4 class for valarrays
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class RungeKutta4_ib {
//--------------------------------------------------------------
//      Decleration of the 4rth order Runge-Kutta Class
//--------------------------------------------------------------
        public:
//          Constructor
            RungeKutta4_ib(valarray<double>& fin, int tout_start);

//          Main function
            RungeKutta4_ib& advance(const double vosc);    
        
//          Access
            double&  tout();
            double&  time();
            size_t&  numh();
            double&  th();

        private:
//          Helper functions
            valarray<double>& F(const valarray<double>& fin, const double vosc);

//          Variables
	    valarray<double>  f0, f1, fh;
            valarray<double>&  f;

            IB_Operator Inversebremsstrahlung;

            size_t num_h;
            double h, t, Tout;
        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class InverseBremsstrahlung{
//--------------------------------------------------------------
//      Decleration of the Phenomenological Laser Source 
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            InverseBremsstrahlung(Stat& Yin, int tout_start); 
         
//          Explicit Advance
            void loop(const double& tnew);  

            valarray<double>  U2, U2m1;//, U1, U1m1;

//          Laser Temporal Evolution
            double t_Profile(const double& t2); 

//          Laser polarization
            int IBSOURCE() const; 

        private:

//          Variables
            Stat& Y;

//          Status 
            int                 InvBremsstrahlung;

//          Algorithm
            RungeKutta4_ib      rk4_ib;

	    double              t, tout;
            int                 Nbc, szx, szy;

//          harmonic data
            valarray<double>    fc;

//          Laser profile
            valarray<double>     IL_xprofile, IL_yprofile;

//          Laser profile
            Axis< double >  xaxis;

//          Define laser profile
            const double rise_time, fall_time, flat_time;
            const double vos, omega_0, omega_p, w0overwp;
            const bool   linear_Ivst, polynomial_Ivst;
        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class NumericsInfo {
//--------------------------------------------------------------
//      Decleration of the Phenomenological Laser Source 
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            NumericsInfo(int tout_start); 
         
//          Explicit Advance
            void Collect_Numerics(Stat& Yin, const double& tnew, const size_t& h_step);    

//          Check if the memory block is full
            bool Export_Numerics() const; 

            int List_Size() const; 

//          Laser Temporal Evolution
//            bool Check_Numerics(); 

            vector<float> Energy_Info;

        private:

            const int mem_block; 
            int       num_checks;
            double    energy_denorm, Energy_in_box;
            bool      disabled;

            Axis< double >  pr, gamma;

//          Energy in the box 
            double Total_Energy(SHarmonic& Y00); 
        };
//--------------------------------------------------------------
//**************************************************************


//--------------------------------------------------------------
     void Make_Gaussian_f(Axis<double>& f, Axis<double>& p, Axis<double>& g, double& Thot);
//--------------------------------------------------------------

//-------------------------------------------------------------------
     double ZLOGei_v2(double ne, double Te); 
//--------------------------------------------------------------


//--------------------------------------------------------------
    #endif
