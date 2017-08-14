///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Jun 4, 2013
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the 
//   data containers:
//
//   1. Spherical Harmonics in 1D
//
//   2. Field in 1D
//
//   3. Distribution function in 1D
// 
//   4. State of the system in 1D
// 
///////////////////////////////////////////////////////////

    #ifndef DECL_STATE_H
    #define DECL_STATE_H

//**************************************************************
//--------------------------------------------------------------
//  Spherical Harmonic Declaration
    class SHarmonic1D {
    private:
        Array2D<double> *sh;

    public:
//      Constructors/Destructors
        SHarmonic1D(size_t nump, size_t numx);
        SHarmonic1D(const SHarmonic1D& other);
        ~SHarmonic1D();

//      Basic information
        Array2D<double>& array() const {return (*sh);}
        size_t dim()  const {return (*sh).dim();}
        size_t nump() const {return (*sh).dim1();}
        size_t numx() const {return (*sh).dim2();}

//      Access to the underlying array 
        double& operator()(size_t i, size_t j) {return (*sh)(i,j);} 
        double  operator()(size_t i, size_t j) const {return (*sh)(i,j);} 
        double& operator()(size_t i) {return (*sh)(i);}             //1D-style 
        double  operator()(size_t i) const {return (*sh)(i);}       //1D-style 

        vector<double> xVec(size_t j) const {return (*sh).d2c(j);}       //1D-style 

//      Access to boundary and guard cells
//      LB -> Left Boundary-cells, LG -> Left Guard-cells, cN -> Number of boundary cells
        GSlice_iter<double>  xLG(size_t cN) {return (*sh).d2c(0,cN-1);}
        GSlice_iter<double>  xLB(size_t cN) {return (*sh).d2c(cN,2*cN-1);}
        GSlice_iter<double>  xRG(size_t cN) {return (*sh).d2c(numx()-cN,numx()-1);}
        GSlice_iter<double>  xRB(size_t cN) {return (*sh).d2c(numx()-2*cN,numx()-cN-1);}

//      Operators
        SHarmonic1D& operator=(const double& d);
        SHarmonic1D& operator=(const SHarmonic1D& other);
        SHarmonic1D& operator*=(const double& d);
        SHarmonic1D& operator*=(const SHarmonic1D& shmulti);
        SHarmonic1D& operator+=(const double& d);
        SHarmonic1D& operator+=(const SHarmonic1D& shadd);
        SHarmonic1D& operator-=(const double& d);
        SHarmonic1D& operator-=(const SHarmonic1D& shmin);

//      Other Algebra
        SHarmonic1D& mpaxis(const valarray<double>& shmulti);
        SHarmonic1D& mxaxis(const valarray<double>& shmulti);

//      Derivatives
        SHarmonic1D& Dp();
        SHarmonic1D& Dx();

//      FilterP
        SHarmonic1D& Filterp(size_t N); 
    };
//--------------------------------------------------------------

// Field in 1D 
  class Field1D {

    private:
        valarray<double> *fi;

    public:
//      Constructors/Destructors
        Field1D(size_t numx);
        Field1D(const Field1D& other);
        ~Field1D();

//      Access to the underlying matrix
        valarray<double>& array() const {return (*fi);}
        size_t numx() const {return (*fi).size();}
        double& operator()(size_t i){ return (*fi)[i];}             //1D-style 
        double  operator()(size_t i) const {return (*fi)[i];} 

//      Operators
        Field1D& operator=(const double& d);
        Field1D& operator=(const valarray<double>& other);
        Field1D& operator=(const Field1D& other);
        Field1D& operator*=(const double& d);
        Field1D& operator*=(const valarray<double>& fimulti);
        Field1D& operator*=(const Field1D& fimulti);
        Field1D& operator+=(const double& d);
        Field1D& operator+=(const Field1D& fiadd);
        Field1D& operator-=(const double& d);
        Field1D& operator-=(const Field1D& fimin);

//      Derivatives
        Field1D& Dx();
    };
//--------------------------------------------------------------

//  1D Distribution function decleration
    class DistFunc1D {
    private:
        vector<SHarmonic1D> *df;
        size_t lmax; 
        double charge, ma;

        vector<SHarmonic1D*> _Neighbors;

    public:
        enum _Compus1D { N, S};  

//      Constructors/Destructors
        DistFunc1D(size_t l, size_t np, size_t nx, double q, double _ma);
        DistFunc1D(const DistFunc1D& other);
        ~DistFunc1D();

//      Basic info
        size_t dim()  const {return lmax+1;}
        size_t l0()   const {return lmax;  }
        double q()    const {return charge;}
        double mass() const {return ma;}

//      Access
        SHarmonic1D* operator()(int l);                   // Returns the address of a harmonic, NULL if out-of-bounds
        SHarmonic1D* operator()(int l) const;
        SHarmonic1D* Compus(size_t l, _Compus1D n) const; // Returns the location of some neighbor of l

//      Operators
        DistFunc1D& operator=(const double& d);
        DistFunc1D& operator=(const SHarmonic1D& h);
        DistFunc1D& operator=(const DistFunc1D& other);
        DistFunc1D& operator*=(const double& d);
        DistFunc1D& operator*=(const DistFunc1D& other);
        DistFunc1D& operator+=(const double& d);
        DistFunc1D& operator+=(const DistFunc1D& other);
        DistFunc1D& operator-=(const double& d);
        DistFunc1D& operator-=(const DistFunc1D& other);

//      Filter
        DistFunc1D& Filterp();
    };
//--------------------------------------------------------------

//  Collection of fields and distribution functions decleration
    class State1D {
    private:
        vector<DistFunc1D> *sp;
        Field1D *ex;

        size_t ns;

    public:
//      Constructors/Destructors
        State1D(size_t nx, vector<size_t> l0, vector<size_t> np, vector<double> q, vector<double> ma);
        State1D(const State1D& other);
        ~State1D();

//      Basic information
        size_t Species() const {return ns;}
        size_t Fields() const  {return 1;}

//      Access to underlying structures
//          Distributions
        DistFunc1D&  DF(size_t s)         {return (*sp)[s];}
        DistFunc1D&  DF(size_t s) const   {return (*sp)[s];}
//          Individual Harmonics
        SHarmonic1D& SH(size_t s, size_t lh)        {return *(((*sp)[s])(lh));} // Reference to spherical harmonic
        SHarmonic1D& SH(size_t s, size_t lh)  const {return *(((*sp)[s])(lh));}
        SHarmonic1D* SHp(size_t s, size_t lh)       {return   ((*sp)[s])(lh); } // Pointer to spherical harmonic
        SHarmonic1D* SHp(size_t s, size_t lh) const {return   ((*sp)[s])(lh); }
//          Fields
        Field1D& Ex() {return *ex;}
        Field1D& Ex() const {return *ex;}

//      Copy assignment Operator
        State1D& operator=(const State1D& other);
        State1D& operator=(const double& d);
        State1D& operator*=(const State1D& other);
        State1D& operator*=(const double& d);
        State1D& operator+=(const State1D& other);
        State1D& operator+=(const double& d);
        State1D& operator-=(const State1D& other);
        State1D& operator-=(const double& d);
    };
//--------------------------------------------------------------
//**************************************************************

    #endif
