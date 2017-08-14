///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Aug 31st 2009
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the data
//   structures that characterize the state of the 
//   system:
//
//   1. class Field2D :
//       This is a class decleration for one of the 
//       components of the fields: Ex, Ey, Ez, Bx, By or Bz.
//
//   2. class EMF2D :
//       This is the collection of the Electromagnetic fields.
// 
//   3. class SHarmonic :
//       This is a class decleration for a single spherical harmonic. It 
//       is a 3D Matrix (|p|,x,y).
//
//   4. class DFunc :
//       This is the class decleration for a collection of harmonics,
//       which is to say the distribution function. In essence 
//       it stores a "trapezoidal" array for the harmonics
//       Access to the individual harmonics for a DFunc F1 is 
//       attained by writing F1(i,j). If i,j are out of range the
//       attempt to access the data will fail. For lmax = 3, mmax = 2
//       the array would look like this:
//       0,0
//       1,0  1,1
//       2,0  2,1  2,2
//       3,0  3,1  3,2
//       where for each couple of numbers we have a "Harmonic class"
//
//   5. class Stat :
//       This is the class decleration for the state of the system,
//       including the distribution function and the fields.
//
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

    #ifndef DECL_STATE_H
    #define DECL_STATE_H

//**************************************************************
//**************************************************************
//  Decleration of the 2D Field class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    class Field2D {
//--------------------------------------------------------------
//  Spherical Harmonic Declaration
//--------------------------------------------------------------

    private:
        Matrix2D< complex<double> > *sh;

    public:
//      Constructors/Destructors
        Field2D(size_t numx, size_t numy);
        Field2D(const Field2D& other);
        ~Field2D();

//      Access to the underlying matrix
        Matrix2D< complex<double> >& matrix() const {return (*sh);}
        size_t numx() const {return (*sh).dim1();}
        size_t numy() const {return (*sh).dim2();}
        complex<double>& operator()(size_t i, size_t j){return (*sh)(i,j);} //Fortran-style
        complex<double>& operator()(size_t i){return (*sh)(i);}             //1D-style 

//      Access to slices
        Slice_iter< complex<double> >  x0(size_t i){return (*sh).d1c(i);}      // constant x0 
        GSlice_iter< complex<double> > x0(size_t i, size_t j)  {return (*sh).d1c(i,j);}
        Slice_iter< complex<double> >  y0(size_t i){return (*sh).d2c(i);}      // constant y0
        Slice_iter< complex<double> >  y0(size_t i, size_t j)  {return (*sh).d2c(i,j);}

//      Access to boundary and guard cells 
//      notation: x0 -> x coordinate is fixed
//                LB -> Left Boundary-cells, LG -> Left Guard-cells
        GSlice_iter< complex<double> > x0_LG(size_t cN) {return x0(0,cN-1);}
        GSlice_iter< complex<double> > x0_LB(size_t cN) {return x0(cN,2*cN-1);}
        GSlice_iter< complex<double> > x0_RG(size_t cN) {return x0(numx()-cN,numx()-1);}
        GSlice_iter< complex<double> > x0_RB(size_t cN) {return x0(numx()-2*cN,numx()-cN-1);}
        Slice_iter< complex<double> >  y0_LG(size_t cN) {return y0(0,cN-1);}
        Slice_iter< complex<double> >  y0_LB(size_t cN) {return y0(cN,2*cN-1);}
        Slice_iter< complex<double> >  y0_RG(size_t cN) {return y0(numy()-cN,numy()-1);}
        Slice_iter< complex<double> >  y0_RB(size_t cN) {return y0(numy()-2*cN,numy()-cN-1);}

        GSlice_iter< complex<double> > TrimGuards(size_t cN);
        GSlice_iter< complex<double> > Subarray(size_t st, size_t nx, size_t ny);

//      Operators
        Field2D& operator=(const complex<double>& d);
        Field2D& operator=(const Field2D& other);
        Field2D& operator*=(const complex<double>& d);
        Field2D& operator*=(const Field2D& shmulti);
        Field2D& operator+=(const complex<double>& d);
        Field2D& operator+=(const Field2D& shadd);
        Field2D& operator-=(const complex<double>& d);
        Field2D& operator-=(const Field2D& shmin);

//      Other Algebra
        Field2D& mxaxis(const valarray< complex<double> >& shmulti);
        Field2D& myaxis(const valarray< complex<double> >& shmulti);
        Field2D& Re();
        Field2D& Im();

//      Derivatives
        Field2D& Dx();
        Field2D& Dy();
    };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//  Decleration of the Electromagnetic Fields Class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    class EMF2D {
//--------------------------------------------------------------
//  Electromagnetic fields decleration
//--------------------------------------------------------------

    private:
        valarray<Field2D> *fie; 

    public:
//      Constructors/Destructors
        EMF2D(size_t nx, size_t ny);
        EMF2D(const EMF2D& other);
        ~EMF2D();

//      Access
        size_t dim()  const {return (*fie).size();}
        Field2D& operator()(size_t i) {return (*fie)[i];}      
        Field2D  operator()(size_t i) const {return (*fie)[i];} 

        Field2D& Ex() {return (*fie)[0];}      
        Field2D& Ey() {return (*fie)[1];}      
        Field2D& Ez() {return (*fie)[2];}      
        Field2D& Bx() {return (*fie)[3];}      
        Field2D& By() {return (*fie)[4];}      
        Field2D& Bz() {return (*fie)[5];}      

//      Operators
        EMF2D& operator=(const complex<double>& d);
        EMF2D& operator=(const Field2D& h);
        EMF2D& operator=(const EMF2D& other);
        EMF2D& operator*=(const complex<double>& d);
        EMF2D& operator*=(const EMF2D& other);
        EMF2D& operator+=(const complex<double>& d);
        EMF2D& operator+=(const EMF2D& other);
        EMF2D& operator-=(const complex<double>& d);
        EMF2D& operator-=(const EMF2D& other);

    };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//  Decleration of the Spherical Harmonic class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    class SHarmonic {
//--------------------------------------------------------------
//  Spherical Harmonic Declaration
//--------------------------------------------------------------

    private:
        Matrix3D< complex<double> > *sh;

    public:
//      Constructors/Destructors
        SHarmonic(size_t nump, size_t numx, size_t numy);
        SHarmonic(const SHarmonic& other);
        ~SHarmonic();

//      Access to the underlying matrix
        Matrix3D< complex<double> >& matrix() const {return (*sh);}
        size_t dim() const {return (*sh).dim();}
        size_t nump() const {return (*sh).dim1();}
        size_t numx() const {return (*sh).dim2();}
        size_t numy() const {return (*sh).dim3();}
        complex<double>& operator()(size_t i, size_t j, size_t k){return (*sh)(i,j, k);} 
        complex<double>& operator()(size_t i){return (*sh)(i);}             //1D-style 

//      Access to slices
        GSlice_iter< complex<double> > p0(size_t i){return (*sh).d1c(i,i);}
        GSlice_iter< complex<double> > p0(size_t b, size_t e) {return (*sh).d1c(b,e);}
        GSlice_iter< complex<double> > x0(size_t i){return (*sh).d2c(i,i);}
        GSlice_iter< complex<double> > x0(size_t b, size_t e) {return (*sh).d2c(b,e);}
        GSlice_iter< complex<double> > y0(size_t i){return (*sh).d3c(i,i);}
        GSlice_iter< complex<double> > y0(size_t b, size_t e) {return (*sh).d3c(b,e);}

//      Access to boundary and guard cells
//      notation: x0 -> x coordinate is fixed
//                LB -> Left Boundary-cells, LG -> Left Guard-cells
        GSlice_iter< complex<double> >  x0_LG(size_t cN) {return x0(0,cN-1);}
        GSlice_iter< complex<double> >  x0_LB(size_t cN) {return x0(cN,2*cN-1);}
        GSlice_iter< complex<double> >  x0_RG(size_t cN) {return x0(numx()-cN,numx()-1);}
        GSlice_iter< complex<double> >  x0_RB(size_t cN) {return x0(numx()-2*cN,numx()-cN-1);}
        GSlice_iter< complex<double> >  y0_LG(size_t cN) {return y0(0,cN-1);}
        GSlice_iter< complex<double> >  y0_LB(size_t cN) {return y0(cN,2*cN-1);}
        GSlice_iter< complex<double> >  y0_RG(size_t cN) {return y0(numy()-cN,numy()-1);}
        GSlice_iter< complex<double> >  y0_RB(size_t cN) {return y0(numy()-2*cN,numy()-cN-1);}

//      Operators
        SHarmonic& operator=(const complex<double>& d);
        SHarmonic& operator=(const SHarmonic& other);
        SHarmonic& operator*=(const complex<double>& d);
        SHarmonic& operator*=(const SHarmonic& shmulti);
        SHarmonic& operator+=(const complex<double>& d);
        SHarmonic& operator+=(const SHarmonic& shadd);
        SHarmonic& operator-=(const complex<double>& d);
        SHarmonic& operator-=(const SHarmonic& shmin);

//      Other Algebra
        SHarmonic& mpaxis(const valarray< complex<double> >& shmulti);
        SHarmonic& mxaxis(const valarray< complex<double> >& shmulti);
        SHarmonic& myaxis(const valarray< complex<double> >& shmulti);
        SHarmonic& mxy_matrix(Matrix2D< complex<double> >& mMatrix);
        SHarmonic& Re();

//      Derivatives
        SHarmonic& Dp();
        SHarmonic& Dx();
        SHarmonic& Dy();
    };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//  Decleration of the Distribution Function Class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    class DFunc {
//--------------------------------------------------------------
//  Distribution function decleration
//--------------------------------------------------------------

    private:
        valarray<SHarmonic> *mf;
        size_t lmax, mmax, sz; 

//      Indexing for "trapezoidal matrix"
        Matrix2D<int> ind;
        Matrix2D<int> indx() const {return ind;} 

    public:
//      Constructors/Destructors
        DFunc(size_t l, size_t m, size_t np, size_t nx, size_t ny);
        DFunc(const DFunc& other);
        ~DFunc();

//      Access
        size_t dim() const {return sz;}
        size_t l0()  const {return lmax;}
        size_t m0()  const {return mmax;}
        SHarmonic& operator()(size_t i) {return (*mf)[i];}      
        SHarmonic  operator()(size_t i) const {return (*mf)[i];} 
        SHarmonic& operator()(size_t il, size_t im) {return (*mf)[ind(il,im)];} 

//      Operators
        DFunc& operator=(const complex<double>& d);
        DFunc& operator=(const SHarmonic& h);
        DFunc& operator=(const DFunc& other);
        DFunc& operator*=(const complex<double>& d);
        DFunc& operator*=(const DFunc& other);
        DFunc& operator+=(const complex<double>& d);
        DFunc& operator+=(const DFunc& other);
        DFunc& operator-=(const complex<double>& d);
        DFunc& operator-=(const DFunc& other);
    };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//  Decleration of the State Class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    class Stat {
//--------------------------------------------------------------
//  Collection of fields and distribution function decleration
//--------------------------------------------------------------

    private:
        DFunc *ff;
        EMF2D *em;

    public:
//      Constructors/Destructors
        Stat();
        Stat(size_t l, size_t m, size_t np, size_t nx, size_t ny);
        Stat(const Stat& other);
        ~Stat();

//      Access to underlying structures
        DFunc& DF() const {return (*ff);}
        SHarmonic& SH(size_t lh,size_t mh) {return (*ff)(lh,mh);}
        SHarmonic& SH(size_t i) {return (*ff)(i);}
        EMF2D& EMF() const {return (*em);}
        Field2D& FLD(size_t ip) const {return (*em)(ip);}

//      Copy assignment Operator
        Stat& operator=(const Stat& other);
        Stat& operator=(const complex<double>& d);
        Stat& operator*=(const Stat& other);
        Stat& operator*=(const complex<double>& d);
        Stat& operator+=(const Stat& other);
        Stat& operator+=(const complex<double>& d);
        Stat& operator-=(const Stat& other);
        Stat& operator-=(const complex<double>& d);
    };
//--------------------------------------------------------------
//**************************************************************

    #endif
