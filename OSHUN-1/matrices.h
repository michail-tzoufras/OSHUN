///////////////////////////////////////////////////////////
//      Author:		Michail Tzoufras
//
//	Last Modified:	Jan 18th 2010
///////////////////////////////////////////////////////////

//
//   This header file contains the definitions for Matrix
//   and Axis types. These Matrix types allow us to perform
//   basic algebra such as "=, *=, +=, -=" as well as matrix
//   operations such as multiplication by a vector. The 
//   ability to do the central difference and iterate over 
//   generalized slices is also attached.
// 
//   1.  template<class T> class Slice_iter :
//       This is an iterator for the elements of a regular
//       matrix type. It utilizes the "slice" function.
// 
//   2.  template<class T> class GSlice_iter :
//       This is an iterator for the elements of a regular
//       matrix type. It utilizes the "gslice" function.
//
//   3. template<class T> class Matrix2D :
//       This is a 2D matrix with the associated algebra
//       and the ability to take the derivatives.
//
//   3. template<class T> class Matrix3D :
//       This is a 3D matrix with the associated algebra
//       and the ability to take the derivatives.
//
//   4. template<class T> class Axis :
//       This is the definition of an "Axis" class. Given
//       the number of points N it creates an Axis with N-1
//       intervals starting from xmin (or -xmax) and ending
//       at xmax. This "Axis" is basically a valarray.
//
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

    #ifndef MATRICES_N_AXIS_H
    #define MATRICES_N_AXIS_H

    using namespace std;

//**************************************************************
//**************************************************************
//  Non-Constant Slice Iterator
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
// forward Declarations to allow friend declarations
//--------------------------------------------------------------
    template<class T> class Slice_iter;
    template<class T> bool operator==(const Slice_iter<T>&, const Slice_iter<T>&);
    template<class T> bool operator!=(const Slice_iter<T>&, const Slice_iter<T>&);
    template<class T> bool operator< (const Slice_iter<T>&, const Slice_iter<T>&);
//--------------------------------------------------------------

//--------------------------------------------------------------
    template<class T> class Slice_iter{
//--------------------------------------------------------------
//  Iterator definition 
//--------------------------------------------------------------
        valarray<T>* v;
        slice s;
        size_t curr;     // index of current element
        T& ref(size_t i) const {return (*v)[s.start()+i*s.stride()];}

    public:
        // Type definitions 
        typedef ptrdiff_t             difference_type;
        typedef forward_iterator_tag  iterator_category;
        typedef T                     value_type;
        typedef T*                    pointer;
        typedef T&                    reference;

        Slice_iter(valarray<T>* vv,slice ss) :v(vv), s(ss), curr(0){ }
        
        Slice_iter end() const {
            Slice_iter t = *this;
            t.curr = s.size();    
            return t;
        }
        
        Slice_iter& operator++() {++curr; return *this;} //prefix
        Slice_iter  operator++(int) {Slice_iter t = *this; ++curr; return t;} //postfix

        T& operator[](size_t i) {return ref(i);}  // C style subscript
        T& operator()(size_t i) {return ref(i);}  // Fortran style subscript
        T& operator*() {return ref(curr);}        // Current element

	friend bool operator==<>(const Slice_iter& p, const Slice_iter& q);
	friend bool operator!=<>(const Slice_iter& p, const Slice_iter& q);
	friend bool operator< <>(const Slice_iter& p, const Slice_iter& q);
    };

//--------------------------------------------------------------
//  Comparison operators for this iterator:
//--------------------------------------------------------------
    template<class T> 
    bool operator==(const Slice_iter<T>& p, const Slice_iter<T>& q){
        return p.curr==q.curr && p.s.stride()==q.s.stride() && p.s.start()==q.s.start();
    }

    template<class T> 
    bool operator!=(const Slice_iter<T>& p, const Slice_iter<T>& q){
        return !(p==q);
    }
    
    template<class T> 
    bool operator< (const Slice_iter<T>& p, const Slice_iter<T>& q){
        return p.curr<q.curr && p.s.stride()==q.s.stride() && p.s.start()==q.s.start();
    }
//--------------------------------------------------------------

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// The Constant Slice Iterator is not currently in use
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//--------------------------------------------------------------
// forward declarations to allow friend declarations
//--------------------------------------------------------------
    template<class T> class Cslice_iter;
    template<class T> bool operator==(const Cslice_iter<T>&, const Cslice_iter<T>&);
    template<class T> bool operator!=(const Cslice_iter<T>&, const Cslice_iter<T>&);
    template<class T> bool operator< (const Cslice_iter<T>&, const Cslice_iter<T>&);
//--------------------------------------------------------------

//--------------------------------------------------------------
    template<class T> class Cslice_iter{
//--------------------------------------------------------------
//  Iterator definition 
//--------------------------------------------------------------
	valarray<T>* v;
	slice s;
	size_t curr; // index of current element
	const T& ref(size_t i) const { return (*v)[s.start()+i*s.stride()]; }

    public:
	Cslice_iter(valarray<T>* vv, slice ss): v(vv), s(ss), curr(0){}

	Cslice_iter end() const {
	    Cslice_iter t = *this;
	    t.curr = s.size(); // index of one plus last element
	    return t;
	}

	Cslice_iter& operator++() { curr++; return *this; }
	Cslice_iter operator++(int) { Cslice_iter t = *this; curr++; return t; }
	
	const T& operator[](size_t i) const { return ref(i); }
	const T& operator()(size_t i) const { return ref(i); }
	const T& operator*() const { return ref(curr); }

	friend bool operator==<>(const Cslice_iter& p, const Cslice_iter& q);
	friend bool operator!=<>(const Cslice_iter& p, const Cslice_iter& q);
	friend bool operator< <>(const Cslice_iter& p, const Cslice_iter& q);

    };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Comparison operators for this iterator:
//--------------------------------------------------------------
    template<class T> 
    bool operator==(const Cslice_iter<T>& p, const Cslice_iter<T>& q){
        return p.curr==q.curr && p.s.stride()==q.s.stride() && p.s.start()==q.s.start();
    }

    template<class T> 
    bool operator!=(const Cslice_iter<T>& p, const Cslice_iter<T>& q){
        return !(p==q);
    }
    
    template<class T> 
    bool operator< (const Cslice_iter<T>& p, const Cslice_iter<T>& q){
        return p.curr<q.curr && p.s.stride()==q.s.stride() && p.s.start()==q.s.start();
    }
//--------------------------------------------------------------
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//**************************************************************

//**************************************************************
//**************************************************************
//  Non-Constant Generalized Slice Iterator
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
// forward Declarations to allow friend declarations
//--------------------------------------------------------------
    template<class T> class GSlice_iter;
    template<class T> bool operator==(const GSlice_iter<T>&, const GSlice_iter<T>&);
    template<class T> bool operator!=(const GSlice_iter<T>&, const GSlice_iter<T>&);
    template<class T> bool operator< (const GSlice_iter<T>&, const GSlice_iter<T>&);
//--------------------------------------------------------------

//--------------------------------------------------------------
    template<class T> class GSlice_iter{
//--------------------------------------------------------------
//  Iterator definition 
//--------------------------------------------------------------
        valarray<T>* v;
        gslice gs;
        size_t curr;     // index of current element
        valarray<size_t> gsizes;
        T& ref(size_t i) const;

    public:
        // Type definitions 
        typedef ptrdiff_t             difference_type;
        typedef forward_iterator_tag  iterator_category;
        typedef T                     value_type;
        typedef T*                    pointer;
        typedef T&                    reference;

        // Constructor
        GSlice_iter(valarray<T>* vv,gslice gss);
        
        // Pointer to the end
        GSlice_iter end() const {
            GSlice_iter t = *this;
            t.curr = gsizes[0];    
            return t;
        }
        
        GSlice_iter& operator++() {++curr; return *this;} //prefix
        GSlice_iter  operator++(int) {GSlice_iter t = *this; ++curr; return t;} //postfix

        T& operator[](size_t i) {return ref(i);}  // C style subscript
        T& operator()(size_t i) {return ref(i);}  // Fortran style subscript
        T& operator*() {return ref(curr);}        // Current element

	friend bool operator==<>(const GSlice_iter& p, const GSlice_iter& q);
	friend bool operator!=<>(const GSlice_iter& p, const GSlice_iter& q);
	friend bool operator< <>(const GSlice_iter& p, const GSlice_iter& q);
    };

//--------------------------------------------------------------
//  Referencing this iterator
//--------------------------------------------------------------
    template<class T> 
    T& GSlice_iter<T>::ref(size_t i) const { 
        size_t loc(0);
        for (size_t ic = 0; ic < gsizes.size()-1; ++ic){
            loc += (i/gsizes[ic+1]) * gs.stride()[ic];
            i %= gsizes[ic+1];
        }
        loc += i * gs.stride()[gsizes.size()-1];
        return (*v)[gs.start()+loc];
    }

//--------------------------------------------------------------
//  Generalized slice iterator constructor
//--------------------------------------------------------------
    template<class T> 
    GSlice_iter<T>::GSlice_iter(valarray<T>* vv,gslice gss) :
       v(vv), gs(gss), curr(0), gsizes(gs.size()){
       for (int ic=1; ic < gsizes.size(); ++ic) 
           gsizes[gss.size().size()-ic-1] *= gsizes[gss.size().size()-ic]; 
    }

//--------------------------------------------------------------
//  Comparison operators for this iterator:
//--------------------------------------------------------------
    template<class T> 
    bool operator==(const GSlice_iter<T>& p, const GSlice_iter<T>& q){
        size_t count(0);
        if (p.curr==q.curr && p.gs.start() == q.gs.start()) 
            while (  (p.gs.stride()[count] == p.gs.stride()[count])
                  && (p.gs.size()[count]   == p.gs.size()[count]  )) 
                if (count++ == q.gs.size().size()-1) return true; 
        return false;
    }

    template<class T> 
    bool operator!=(const GSlice_iter<T>& p, const GSlice_iter<T>& q){
        return !(p==q);
    }
    
    template<class T> 
    bool operator< (const GSlice_iter<T>& p, const GSlice_iter<T>& q){
        size_t count(0);
        if (p.curr<q.curr && p.gs.start() == q.gs.start()) 
            while (  (p.gs.stride()[count] == p.gs.stride()[count])
                  && (p.gs.size()[count]   == p.gs.size()[count]  ))
                if (count++ == q.gs.size().size()-1) return true; 
        return false;
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   2D Matrix with algebra attached
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    template<class T> class Matrix2D {
//--------------------------------------------------------------
//  2D Matrix decleration
//--------------------------------------------------------------
    private:
        valarray<T> *v;
        size_t  d1, d2;            // x, y or ... p, x

    public:
//      Constructors/Destructors
        Matrix2D(size_t x, size_t y);
        Matrix2D(const Matrix2D& other);
        ~Matrix2D();

//      Access
        size_t dim()  const {return d1*d2;}
        size_t dim1() const {return d1;}
        size_t dim2() const {return d2;} 
        valarray<T>& array() const {return *v;}
        T& operator()(size_t i, size_t j); //Fortran-style 
        T& operator()(size_t i);           //1D-style 

        Slice_iter<T>  d2c(size_t i);       // contiguous elements
        Cslice_iter<T> d2c(size_t i) const; 
        Slice_iter<T>  d1c(size_t i);       // non-contiguous (d1 jumps)
        Cslice_iter<T> d1c(size_t i) const; 
        Slice_iter<T>  d2c(size_t b, size_t e);
        GSlice_iter<T> d1c(size_t b, size_t e);

        GSlice_iter<T> SubMatrix2D(size_t st, size_t nx, size_t ny);
        GSlice_iter<T> gsli(size_t st, const valarray<size_t>& sz,const valarray<size_t>& str);

//      Operators
        Matrix2D& operator=(const T& d);
        Matrix2D& operator=(const Matrix2D& other);
        Matrix2D& operator*=(const T& d);
        Matrix2D& operator*=(const Matrix2D& vmulti);
        Matrix2D& operator+=(const T& d);
        Matrix2D& operator+=(const Matrix2D& vadd);
        Matrix2D& operator-=(const T& d);
        Matrix2D& operator-=(const Matrix2D& vmin);

//      Matrix * Vector
        Matrix2D& multid1(const valarray<T>& vmulti); // M*valarray(d1)
        Matrix2D& multid2(const valarray<T>& vmulti); // M*valarray(d2)

//      Central difference
        Matrix2D& Dd1(); // in the direction d1 column-wise
        Matrix2D& Dd2(); // in the direction d2 row-wise

    };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
    template<class T> Matrix2D<T>:: Matrix2D(size_t x,size_t y) : d1(x), d2(y) {
        v = new valarray<T>(d1*d2);

    }
//  Copy constructor
    template<class T> Matrix2D<T>:: Matrix2D(const Matrix2D& other){
        d1  = other.dim1();
        d2  = other.dim2();
        v = new valarray<T>(d1*d2);
        (*v) = other.array();
    }
//  Destructor
    template<class T> Matrix2D<T>:: ~Matrix2D(){
        delete v; // automagically calls v->~valarray<T>; 
    }

//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------

//  constant d2 :: contiguous elements :: access by column 
    template<class T> inline Slice_iter<T> Matrix2D<T>:: d2c(size_t i){
        return Slice_iter<T>(v,slice(i*d1,d1,1)); 
    }
//  constant d2 :: contiguous elements :: access by column 
    template<class T> inline Cslice_iter<T> Matrix2D<T>:: d2c(size_t i) const{
        return Cslice_iter<T>(v,slice(i*d1,d1,1)); 
    }
//  constant d1 :: noncontuguous elements :: access by row 
    template<class T> inline Slice_iter<T> Matrix2D<T>:: d1c(size_t i){
        return Slice_iter<T>(v,slice(i,d2,d1)); 
    }
//  constant d1 :: noncontuguous elements :: access by row 
    template<class T> inline Cslice_iter<T> Matrix2D<T>:: d1c(size_t i) const{
        return Cslice_iter<T>(v,slice(i,d2,d1)); 
    }
//  scan surfaces of multiple constant d2
    template<class T> 
    inline Slice_iter<T> Matrix2D<T>::d2c(size_t b, size_t e){  
         return Slice_iter<T>(v,slice(b*d1,(e-b+1)*d1,1)); 
    }
//  scan surfaces of multiple constant d1
    template<class T> 
    inline GSlice_iter<T> Matrix2D<T>::d1c(size_t b, size_t e){  
                 valarray<size_t> sz(2), str(2);
                 str[1] = d1; str[0] = 1;
                 sz[1]  = d2;  sz[0] = e-b+1;
         return GSlice_iter<T>(v,gslice(b,sz,str));
    }
//  scan Submatrix
    template<class T> 
    inline GSlice_iter<T> Matrix2D<T>::SubMatrix2D(size_t st, size_t nx, size_t ny){ 
                 valarray<size_t> sz(2), str(2);
                 str[1] = 1;  str[0] = dim1();
                 sz[1]  = nx;  sz[0] = ny;
         return GSlice_iter<T>(v,gslice(st,sz,str));
    }
//  scan arbitrary generalized slice
    template<class T> 
    inline GSlice_iter<T> Matrix2D<T>::gsli(size_t st, 
                 const valarray<size_t>& sz,const valarray<size_t>& str){
         return GSlice_iter<T>(v,gslice(st,sz,str));
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


//  Access Fortan-style
    template<class T> inline T& Matrix2D<T>:: operator()(size_t i, size_t j){
        return (*v)[i+j*d1];
    }
//  1D style access
    template<class T> inline T& Matrix2D<T>:: operator() (size_t i){
        return (*v)[i];
    }
//  Access valarray


//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
    template<class T> Matrix2D<T>& Matrix2D<T>::operator=(const T& d){
        (*v) = d;
        return *this;
    }
    template<class T> Matrix2D<T>& Matrix2D<T>::operator=(const Matrix2D& other){
        if (this != &other) {   //self-assignment
           (*v) = other.array();
        }
        return *this;
    }

//  *= 
    template<class T> Matrix2D<T>& Matrix2D<T>::operator*=(const T& d){
        (*v) *=d;
        return *this;
    }
    template<class T> Matrix2D<T>& Matrix2D<T>::operator*=(const Matrix2D& vmulti){
        (*v) *= vmulti.array();
        return *this;
    }

//  +=
    template<class T> Matrix2D<T>& Matrix2D<T>::operator+=(const T& d){
        (*v) +=d;
        return *this;
    }
    template<class T> Matrix2D<T>& Matrix2D<T>::operator+=(const Matrix2D& vadd){
        (*v) += vadd.array();
        return *this;
    }

//  -= 
    template<class T> Matrix2D<T>& Matrix2D<T>::operator-=(const T& d){
        (*v) -=d;
        return *this;
    }
    template<class T> Matrix2D<T>& Matrix2D<T>::operator-=(const Matrix2D& vmin){
        (*v) -= vmin.array();
        return *this;
    }

//--------------------------------------------------------------
//  Matrix2D * Vector
//--------------------------------------------------------------
//  Vector multiplies every column of Matrix2D
    template<class T> Matrix2D<T>& Matrix2D<T>::multid1(const valarray<T>& vmulti){
        for(size_t i=0; i< d1; ++i)
            for(Slice_iter<T> it=d1c(i); it!=it.end(); ++it)
                *it *= vmulti[i];
        return *this;
    }

//  Vector multiplies every row of Matrix2D
    template<class T> Matrix2D<T>& Matrix2D<T>::multid2(const valarray<T>& vmulti){
        for(size_t i=0; i< d2; ++i)
            for (Slice_iter<T> it = d2c(i); it != it.end(); ++it)
               *it *= vmulti[i];
        return *this;
    }

//--------------------------------------------------------------
// Central difference
//--------------------------------------------------------------
// (minus) Central difference for contiguous elements  
// Example: 0 4  8 12 16       -2 -2 -2 -2 -2 
//          1 5  9 13 17  -->  -2 -2 -2 -2 -2 
//          2 6 10 14 18       -2 -2 -2 -2 -2  
//          3 7 11 15 19       -2 -2 -2 -2 19  
    template<class T> Matrix2D<T>& Matrix2D<T>::Dd1(){
        for(long i=0; i< d1*d2-2; ++i)
            (*v)[i] -= (*v)[i+2]; 
        for(long i=d1*d2-3; i>-1; --i)
            (*v)[i+1] = (*v)[i]; 
        return *this;
    }

// (minus) Central difference for elements with distance d1 (row-wise)
// Example: 0 4  8 12 16       -8 -8 -8 -8 16 
//          1 5  9 13 17  -->  -8 -8 -8 -8 17 
//          2 6 10 14 18       -8 -8 -8 -8 18
//          3 7 11 15 19       -8 -8 -8 -8 19
    template<class T> Matrix2D<T>& Matrix2D<T>::Dd2(){
        int twod1 = 2*d1;
        for(long i=0; i< d1*d2-twod1; ++i)
            (*v)[i] -= (*v)[i+twod1]; 
        for(long i=d1*d2-twod1-1; i>-1; --i)
            (*v)[i+d1] = (*v)[i]; 
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   3D Matrix with algebra 
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    template<class T> class Matrix3D {
//--------------------------------------------------------------
//  3D Matrix decleration 
//--------------------------------------------------------------

    private:
        valarray<T> *v;
        size_t  d1, d2, d3;            // x, y, z or ... p, x, y
        size_t d1d2; 

    public:
//      Constructors/Destructors
        Matrix3D(size_t x, size_t y, size_t z);
        Matrix3D(const Matrix3D& other);
        ~Matrix3D();

//      Access
        size_t dim()  const {return d1*d2*d3;}
        size_t dim1() const {return d1;}
        size_t dim2() const {return d2;}
        size_t dim3() const {return d3;}
        valarray<T>&  array() const {return *v;};
        T& operator()(size_t i, size_t j, size_t k); //Fortran-style 
        T& operator()(size_t i);                     //1D-style 
        Slice_iter<T>  d2d3c(size_t j, size_t k);       // contiguous elements
        Cslice_iter<T> d2d3c(size_t j, size_t k) const;  
        Slice_iter<T>  d1d3c(size_t i, size_t k);       // non-contiguous (d1 jumps)
        Cslice_iter<T> d1d3c(size_t i, size_t k) const;  
        Slice_iter<T>  d1d2c(size_t i, size_t j);       // non-contiguous (d1*d2 jumps)
        Cslice_iter<T> d1d2c(size_t i, size_t j) const;  

        GSlice_iter<T> d1c(size_t b, size_t e);
        GSlice_iter<T> d2c(size_t b, size_t e);
        GSlice_iter<T> d3c(size_t b, size_t e);
        GSlice_iter<T> gsli(size_t st, const valarray<size_t>& sz,const valarray<size_t>& str);

        GSlice_iter<T> SubMatrix3D(size_t st, size_t nx, size_t ny, size_t nz);

//      Operators
        Matrix3D& operator=(const T& d);
        Matrix3D& operator=(const Matrix3D& other);
        Matrix3D& operator*=(const T& d);
        Matrix3D& operator*=(const Matrix3D& vmulti);
        Matrix3D& operator+=(const T& d);
        Matrix3D& operator+=(const Matrix3D& vadd);
        Matrix3D& operator-=(const T& d);
        Matrix3D& operator-=(const Matrix3D& vmin);

//      Matrix*Vector
        Matrix3D& multid1(const valarray<T>& vmulti);// M*valarray(d1) 
        Matrix3D& multid2(const valarray<T>& vmulti);// M*valarray(d2) 
        Matrix3D& multid3(const valarray<T>& vmulti);// M*valarray(d3) 

//      Central difference
        Matrix3D& Dd1(); // in the direction d1
        Matrix3D& Dd2(); // in the direction d2
        Matrix3D& Dd3(); // in the direction d3 
    };
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------

//  Constructor
    template<class T> Matrix3D<T>:: 
    Matrix3D(size_t x, size_t y, size_t z) : d1(x), d2(y), d3(z) {
        v = new valarray<T>(d1*d2*d3);
        d1d2 = d1*d2;
    }
//  Copy constructor
    template<class T> Matrix3D<T>:: Matrix3D(const Matrix3D& other){
        d1  = other.dim1();
        d2  = other.dim2();
        d3  = other.dim3();
        d1d2 = d1*d2;
        v = new valarray<T>(d1*d2*d3);
        (*v) = other.array();
    }
//  Destructor
    template<class T> Matrix3D<T>:: ~Matrix3D(){
        delete v; 
    }

//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Slices
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  scan d1 :: contiguous elements
    template<class T> 
    inline Slice_iter<T> Matrix3D<T>:: d2d3c(size_t j, size_t k){
        return Slice_iter<T>(v,slice(j*d1+k*d1d2,d1,1)); 
    }
//  scan d1 :: contiguous elements 
    template<class T> 
    inline Cslice_iter<T> Matrix3D<T>:: d2d3c(size_t j, size_t k) const {
        return Cslice_iter<T>(v,slice(j*d1+k*d1d2,d1,1)); 
    }
//  scan d2 :: noncontuguous elements 
    template<class T> 
    inline Slice_iter<T> Matrix3D<T>:: d1d3c(size_t i,size_t k){
        return Slice_iter<T>(v,slice(i+k*d1d2,d2,d1)); 
    }
//  scan d2 :: noncontuguous elements
    template<class T> 
    inline Cslice_iter<T> Matrix3D<T>:: d1d3c(size_t i,size_t k) const {
        return Slice_iter<T>(v,slice(i+k*d1d2,d2,d1)); 
    }
//  scan d3 :: noncontiguous elements
    template<class T> 
    inline Slice_iter<T> Matrix3D<T>:: d1d2c(size_t i, size_t j){
        return Slice_iter<T>(v,slice(i+j*d1,d3,d1*d2)); 
    }
//  scan d3 :: noncontiguous elements
    template<class T> 
    inline Cslice_iter<T> Matrix3D<T>:: d1d2c(size_t i, size_t j) const {
        return Cslice_iter<T>(v,slice(i+j*d1,d3,d1*d2)); 
    }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Generalized Slices
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  scan surfaces of constant d1
    template<class T> 
    inline GSlice_iter<T> Matrix3D<T>::d1c(size_t b, size_t e){  
                 valarray<size_t> sz(3), str(3);
                 str[2] = d1; str[1] = d1d2; str[0] = 1;
                 sz[2]  = d2; sz[1]  = d3;    sz[0] = e-b+1;
         return GSlice_iter<T>(v,gslice(b,sz,str));
    }
//  scan surfaces of constant d2
    template<class T> 
    inline GSlice_iter<T> Matrix3D<T>::d2c(size_t b, size_t e){  
                 valarray<size_t> sz(3), str(3);
                 str[2] = 1;  str[1] = d1d2;  str[0] = d1;
                 sz[2]  = d1; sz[1]  = d3;    sz[0] = e-b+1;
         return GSlice_iter<T>(v,gslice(d1*b,sz,str));
    }
//  scan surfaces of constant d3
    template<class T> 
    inline GSlice_iter<T> Matrix3D<T>::d3c(size_t b, size_t e){  
                 valarray<size_t> sz(2), str(2);
                 str[1] = 1; str[0] = d1d2;
                 sz[1]  = d1d2; sz[0] = e-b+1;
         return GSlice_iter<T>(v,gslice(d1d2*b,sz,str));
    }
//  scan arbitrary generalized slice
    template<class T> 
    inline GSlice_iter<T> Matrix3D<T>::gsli(size_t st, 
                 const valarray<size_t>& sz,const valarray<size_t>& str){
         return GSlice_iter<T>(v,gslice(st,sz,str));
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  scan Submatrix
    template<class T> 
    inline GSlice_iter<T> Matrix3D<T>::SubMatrix3D(size_t st, size_t nx, size_t ny, size_t nz){ 
                 valarray<size_t> sz(3), str(3);
                 str[2] = 1;  str[1] = dim1(); str[0] = d1d2;
                 sz[2]  = nx;  sz[1] = ny;      sz[0] = nz;
         return GSlice_iter<T>(v,gslice(st,sz,str));
    }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//  Access Fortan-style
    template<class T> 
    inline T& Matrix3D<T>:: operator()(size_t i, size_t j, size_t k){
        return (*v)[i+j*d1+k*d1d2];
    }

//  1D style access
    template<class T> 
    inline T& Matrix3D<T>:: operator()(size_t i){
        return (*v)[i];
    }

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------

//  Copy assignment operator
    template<class T> Matrix3D<T>& Matrix3D<T>::operator=(const T& d){
        (*v) = d;
        return *this;
    }
//  Copy assignment operator
    template<class T> Matrix3D<T>& Matrix3D<T>::operator=(const Matrix3D& other){
        if (this != &other) {   //self-assignment
           (*v) = other.array();
        }
        return *this;
    }

//  *= 
    template<class T> Matrix3D<T>& Matrix3D<T>::operator*=(const T& d){
        (*v) *=d;
        return *this;
    }
    template<class T> Matrix3D<T>& Matrix3D<T>::operator*=(const Matrix3D& vmulti){
        (*v) *= vmulti.array();
        return *this;
    }

//  +=
    template<class T> Matrix3D<T>& Matrix3D<T>::operator+=(const T& d){
        (*v) +=d;
        return *this;
    }
    template<class T> Matrix3D<T>& Matrix3D<T>::operator+=(const Matrix3D& vadd){
        (*v) += vadd.array();
        return *this;
    }

//  -= 
    template<class T> Matrix3D<T>& Matrix3D<T>::operator-=(const T& d){
        (*v) -=d;
        return *this;
    }
    template<class T> Matrix3D<T>& Matrix3D<T>::operator-=(const Matrix3D& vmin){
        (*v) -= vmin.array();
        return *this;
    }
//--------------------------------------------------------------
//  Matrix3D * Vector
//--------------------------------------------------------------
//  Vector multiplies every d1-length of Matrix3D
    template<class T> Matrix3D<T>& Matrix3D<T>::multid1(const valarray<T>& vmulti){
        for(int i=0; i< d1; ++i)
            for(int j=0; j< d2; ++j)
                for(Slice_iter<T> it(d1d2c(i,j)); it!=it.end(); ++it)
                    *it *= vmulti[i];
        return *this;
    }

//  Vector multiplies every d2-length of Matrix3D
    template<class T> Matrix3D<T>& Matrix3D<T>::multid2(const valarray<T>& vmulti){
        for(int j=0; j< d2; ++j)
            for(int k=0; k< d3; ++k)
                for(Slice_iter<T> it(d2d3c(j,k)); it!=it.end(); ++it)
                    *it *= vmulti[j];
        return *this;
    }

//  Vector multiplies every d3-length of Matrix3D
    template<class T> Matrix3D<T>& Matrix3D<T>::multid3(const valarray<T>& vmulti){
        for(int k=0; k< d3; ++k)
            for(int i=0; i< d1; ++i)
                for(Slice_iter<T> it(d1d3c(i,k)); it!=it.end(); ++it)
                    *it *= vmulti[k];
        return *this;
    }

//--------------------------------------------------------------
// Central difference
//--------------------------------------------------------------
// (minus) Central difference for contiguous elements  
    template<class T> Matrix3D<T>& Matrix3D<T>::Dd1(){
        for(long i=0; i< d1*d2*d3-2; ++i)
            (*v)[i] -= (*v)[i+2]; 
        for(long i=d1*d2*d3-3; i>-1; --i)
            (*v)[i+1] = (*v)[i]; 
        return *this;
    }

// (minus) Central difference in the d2 direction
    template<class T> Matrix3D<T>& Matrix3D<T>::Dd2(){
        long twod1 = 2*d1;
        for(long i=0; i< d1*d2*d3-twod1; ++i)
            (*v)[i] -= (*v)[i+twod1]; 
        for(long i=d1*d2*d3-twod1-1; i>-1; --i)
            (*v)[i+d1] = (*v)[i]; 
        return *this;
    }

// (minus) Central difference in the d3 direction
    template<class T> Matrix3D<T>& Matrix3D<T>::Dd3(){
        long twod1d2 = 2*d1d2;
        for(long i=0; i< d1*d2*d3-twod1d2; ++i)
            (*v)[i] -= (*v)[i+twod1d2]; 
        for(long i=d1*d2*d3-twod1d2-1; i>-1; --i)
            (*v)[i+d1d2] = (*v)[i]; 
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//  Axis template definition
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    template<class T> class Axis{
//--------------------------------------------------------------
//  Class used for Axis
//--------------------------------------------------------------
    
    private:
       valarray<T> *ax;

    public:
//      Constructor (for symmetric axis)
        Axis(size_t numx, T xmax){
            ax = new  valarray<T>(numx);
            for (int i = 0; i < numx; ++i)
                (*ax)[i] = i;
            (*ax) *= xmax * (2.0/(static_cast<T>(numx-1)));
        (*ax) -= xmax;
        }
//      Constructor (for nonsymmetric axis)
        Axis(size_t numx, T xmin, T xmax){
            ax = new  valarray<T>(numx);
            for (int i = 0; i < numx; ++i)
                (*ax)[i] = i;
            (*ax) *= (xmax-xmin)/(static_cast<T>(numx-1));
        (*ax) += xmin;
        }
//      Copy Constructor
        Axis(const Axis& other){
            ax = new  valarray<T>(other.dim());
            *ax = other.array();
        }
//      Destructor 
        ~Axis() {delete ax;}

//      Access
        size_t dim()  const {return (*ax).size();}
        valarray<T>& array() const {return *ax;}
        T& operator()(size_t i) {return (*ax)[i];}
        T& operator()(size_t i) const {return (*ax)[i];}
        T dx() {return ((*ax)[1]-(*ax)[0]);}

//      Copy Assignment Operator
        Axis& operator=(const Axis& other){
            if (this != &other){   
               *ax = other.array();
            } 
            return *this;
        }
//      Operators
        Axis& operator*=(const T& d) {*ax *= d; return *this;}
        Axis& operator*=(const Axis& other) {*ax *= other.array(); return *this;}
        Axis& operator+=(const T& d) {*ax += d; return *this;}
        Axis& operator+=(const Axis& other) {*ax += other.array(); return *this;}
        Axis& operator-=(const T& d) {*ax -= d; return *this;}
        Axis& operator-=(const Axis& other) {*ax -= other.array(); return *this;}
    };
//--------------------------------------------------------------
//**************************************************************

    #endif
