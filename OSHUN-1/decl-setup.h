///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Apr 12 2011
///////////////////////////////////////////////////////////

//   
//   This header contains the definitions for the temporary
//   initialization routines. 
// 
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

    #ifndef DECLERATION_SETPUT_H
    #define DECLERATION_SETPUT_H

//**************************************************************
//**************************************************************
//   Declerations for the Input namespace
//**************************************************************
//*************************************************************

    namespace Setup_Parameters{  
        int makefolder(string _name);
        void folders();
    }


    namespace Setup_Y {  

//      Initialize the appropriate density and temperature profiles (from the list below)
        void initialize(Stat& Y);

//      Density profiles
        void Polynomial_x(vector<double>& smooth_xdens, vector<double>& long_x, // const Axis<double>& x, 
                         double start_x0, double rise_x, double flat_x,
                         double fall_x, double nmin_x);
        void Polynomial_x_har(SHarmonic& h, const Axis<double>& x, 
                         double start_x0, double rise_x, double flat_x,
                         double fall_x, double nmin_x);
        void Exponential_x(vector<double>& smooth_xdens, vector<double>& long_x, // const Axis<double>& x, 
                         double exp_nmin, double exp_xmin, 
                         double exp_nmax, double lambda_0);
        void Smooth121(vector<double>& dens, size_t smooth_level);
        void Smooth121_x(SHarmonic& h, const Axis<double>& x, size_t smooth_level, int bndX);
        void Gaussian_x_har(SHarmonic& h, const Axis<double>& x, double sigma, double center, double ampl);
        void Gaussian_y(SHarmonic& h, const Axis<double>& x, double sigma, double center, double ampl);
        void Piecewise_Linear_x(vector<double>& smooth_xdens, vector<double>& long_x, // const Axis<double>& x,
                                 vector<double>& xloc, vector<double>& xdens, int bndX);
        void Piecewise_Linear_y(SHarmonic& h, const Axis<double>& y, vector<double>& yloc, vector<double>& ydens, int bndY);

//      Temperature profiles
        void Gaussian_tx(Matrix2D<double>& Temp, const Axis<double>& x, double sigma_tx, double center_tx, double ptx);
        void Gaussian_ty(Matrix2D<double>& Temp, const Axis<double>& y, double sigma_ty, double center_ty, double pty);
        void Piecewise_Linear_tmpx(Matrix2D<double>& Temp, const Axis<double>& x, vector<double>& xtloc, vector<double>& xtemp);
        void Sinusoidal_tmpx(Matrix2D<double>& Temp, const Axis<double>& x, double pmin, double pmax, double box_size);
        void Sinusoidal_tmpy(Matrix2D<double>& Temp, const Axis<double>& y, double pmin, double pmax, double box_size);
        void Piecewise_Linear_tmpy(Matrix2D<double>& Temp, const Axis<double>& y, vector<double>& ytloc, vector<double>& ytemp);
        void Gaussian_p(SHarmonic& h, const Axis<double>& x, Matrix2D<double>& pt);

//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//      TWO STREAM STUFF
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
        void Gaussian_p(SHarmonic& h, Axis<double>& x, double sigma, double center);
        void Gaussian_p_sinx( SHarmonic& h, Axis<double>& p, double sigma, double center,
                                                  Axis<double>& x, double ampl );
        void Constant_p(SHarmonic& h, Axis<double>& x, double value, double center);
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
//      WEIBEL STUFF 
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><
        void Gaussian_p_sinxsiny( SHarmonic& h, Axis<double>& p, double sigma, double center,
                                                  Axis<double>& x, Axis<double>& y, double ampl );
        void Gaussian_p_NOTsinxsiny( SHarmonic& h, Axis<double>& p, double sigma, double center,
                                                  Axis<double>& x, Axis<double>& y, double ampl );
//><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><-><

    }



    #endif
