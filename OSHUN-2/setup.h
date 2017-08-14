///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Jun 4, 2013
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the 
//   operators:
//
//   1. Spatial Advection
//
//   2. Electric field
// 
//   3. Current 
// 
//   4. Runge-Kutta 1D Functor
// 
///////////////////////////////////////////////////////////

    #ifndef DECL_SETUP_H
    #define DECL_SETUP_H



//**************************************************************
        

class Grid_Info {
public:
//  Constructor
    Grid_Info(const vector<size_t> _Nl,   const vector<size_t> _Nm, 
              const vector<double> _mass, const vector<double> _charge,
//            local  spatial axes
              const vector<double> _xmin, const vector<double> _xmax,  const vector<size_t> _Nx,     
//            global spatial axes
              const vector<double> _xgmin,const vector<double> _xgmax, const vector<size_t> _Nxg,   
//            				  momentum axes
                                          const vector<double> _pmax,  const vector<size_t> _Np,     
//            							       output momentum axes
                                                                       const vector<size_t> _Npx) 
	      : Nl(_Nl), Nm(_Nm), Np(_Np), mass(_mass), charge(_charge),
                axis(_xmin, _xmax, _Nx, _xgmin, _xgmax, _Nxg, _pmax, _Np, _Npx){ }

//  Copy constructor
    Grid_Info(const Grid_Info& other): Nl(other.Nl), Nm(other.Nm), Np(other.Np),
                                        mass(other.mass), charge(other.charge),
                                        axis(other.axis){}
    ~Grid_Info(){};

    const vector<size_t> Nl;
    const vector<size_t> Nm;
    const vector<size_t> Np;
    const vector<double> mass;
    const vector<double> charge;

    const Algorithms::AxisBundle<double> axis;
};
//--------------------------------------------------------------
//**************************************************************

    #endif
