///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Aug 31st 2009
///////////////////////////////////////////////////////////

//   
//   This file contains the definitions for the data
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

// Standard Libraries
   #include <iostream>
   #include <vector>
   #include <valarray>
   #include <complex>

// My Libraries
   #include "matrices.h"

// Declerations
   #include "decl-input.h"
   #include "decl-state.h"


//**************************************************************
//**************************************************************
//  Definition of the "Field2D" Class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------

//  Constructor
    Field2D:: Field2D(size_t numx,size_t numy) {
        sh = new  Matrix2D< complex<double> >(numx,numy);
    }
//  Copy constructor
    Field2D:: Field2D(const Field2D& other){
        sh = new  Matrix2D< complex<double> >(other.numx(),other.numy());
        *sh = other.matrix();
    }
//  Destructor
    Field2D:: ~Field2D(){
        delete sh; 
    }

//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------
    GSlice_iter< complex<double> >  Field2D::TrimGuards(size_t cN){
        return (*sh).SubMatrix2D(cN*numx()+cN, numx()-2*cN, numy()-2*cN);
    }
    GSlice_iter< complex<double> >  Field2D::Subarray(size_t st, size_t nx, size_t ny){
        return (*sh).SubMatrix2D(st, nx, ny);
    }

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------

//  Copy assignment operator
    Field2D& Field2D::operator=(const  complex<double>& d){
        *sh = d;
        return *this;
    }
    Field2D& Field2D::operator=(const Field2D& other){
        if (this != &other) {   //self-assignment
           *sh = other.matrix(); 
        }
        return *this;
    }
//  *= 
    Field2D& Field2D::operator*=(const complex<double>& d){
        (*sh) *=d;
        return *this;
    }
    Field2D& Field2D::operator*=(const Field2D& shmulti){
        (*sh) *= shmulti.matrix();
        return *this;
    }
//  +=
    Field2D& Field2D::operator+=(const complex<double>& d){
        (*sh) +=d;
        return *this;
    }
    Field2D& Field2D::operator+=(const Field2D& shadd){
        (*sh) += shadd.matrix(); 
        return *this;
    }
//  -= 
    Field2D& Field2D::operator-=(const complex<double>& d){
        (*sh) -=d;
        return *this;
    }
    Field2D& Field2D::operator-=(const Field2D& shmin){
        (*sh) -= shmin.matrix();
        return *this;
    }
//--------------------------------------------------------------
//   Other Algebra
//--------------------------------------------------------------
    Field2D& Field2D::mxaxis(const valarray< complex<double> >& shmulti){
        (*sh).multid1(shmulti);
        return *this;
    }
    Field2D& Field2D::myaxis(const valarray< complex<double> >& shmulti){
        (*sh).multid2(shmulti);
        return *this;
    }
    Field2D& Field2D::Re(){
        for (int i(0); i < numx()*numy(); ++i) (*sh)(i) = (*sh)(i).real();
        return *this;
    }
    Field2D& Field2D::Im(){
        for (int i(0); i < numx()*numy(); ++i) (*sh)(i) -= (*sh)(i).real();
        return *this;
    }
//--------------------------------------------------------------
    Field2D& Field2D::Dx(){
//--------------------------------------------------------------
        *sh = (*sh).Dd1();
        return *this;
    }
//--------------------------------------------------------------
    Field2D& Field2D::Dy(){
//--------------------------------------------------------------
        *sh = (*sh).Dd2();
        return *this;
    } 
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//  Definition of the Electromagnetic Fields Class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
    EMF2D:: EMF2D(size_t nx, size_t ny) {
        fie = new valarray<Field2D>(Field2D(nx,ny), 6); 

    }
//  Copy constructor
    EMF2D:: EMF2D(const EMF2D& other){
        fie = new valarray<Field2D>(Field2D(other(0).numx(),other(0).numy()),other.dim()); 
        for(int i=0; i < other.dim() ; ++i) (*fie)[i] = other(i); 
    }
//  Destructor
    EMF2D:: ~EMF2D(){
        delete fie;
    }
//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------

//  Copy assignment operator
    EMF2D& EMF2D::operator=(const complex<double>& d){
        for(int i=0; i < dim() ; ++i)  
            (*fie)[i] = d;
        return *this;
    }
    EMF2D& EMF2D::operator=(const Field2D& h){
        for(int i=0; i < dim() ; ++i){  
            if (&((*fie)[i]) != &h) {   //self-assignment
                (*fie)[i] = h;
            }
        }
        return *this;
    }
    EMF2D& EMF2D::operator=(const EMF2D& other){
        if (this != &other) {   //self-assignment
            for(int i=0; i < dim() ; ++i)  
                (*fie)[i] = other(i);
        }
        return *this;
    }
//  *=
    EMF2D& EMF2D::operator*=(const complex<double>& d){
        for(int i=0; i < dim() ; ++i)  
            (*fie)[i] *= d;
        return *this;
    }
    EMF2D& EMF2D::operator*=(const EMF2D& other){
        if (this != &other) {   //self-assignment
            for(int i=0; i < dim() ; ++i)  
                (*fie)[i] *= other(i);
        }
        return *this;
    }
//  +=
    EMF2D& EMF2D::operator+=(const complex<double>& d){
        for(int i=0; i < dim() ; ++i)  
            (*fie)[i] += d;
        return *this;
    }
    EMF2D& EMF2D::operator+=(const EMF2D& other){
        if (this != &other) {   //self-assignment
            for(int i=0; i < dim() ; ++i)  
                (*fie)[i] += other(i);
        }
        return *this;
    }
//  -=
    EMF2D& EMF2D::operator-=(const complex<double>& d){
        for(int i=0; i < dim() ; ++i)  
            (*fie)[i] -= d;
        return *this;
    }
    EMF2D& EMF2D::operator-=(const EMF2D& other){
        if (this != &other) {   //self-assignment
            for(int i=0; i < dim() ; ++i)  
                (*fie)[i] -= other(i);
        }
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//  Definition of the "SHarmonic" Class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------

//  Constructor
    SHarmonic::SHarmonic(size_t nump, size_t numx,size_t numy) {
        sh = new  Matrix3D< complex<double> >(nump,numx,numy);
    }
//  Copy constructor
    SHarmonic::SHarmonic(const SHarmonic& other){
        sh = new  Matrix3D< complex<double> >(other.nump(),other.numx(),other.numy());
        *sh = other.matrix();
    }
//  Destructor
    SHarmonic:: ~SHarmonic(){
        delete sh; 
    }

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------

//  Copy assignment operator
    SHarmonic& SHarmonic::operator=(const  complex<double>& d){
        *sh = d;
        return *this;
    }
    SHarmonic& SHarmonic::operator=(const SHarmonic& other){
        if (this != &other) {   //self-assignment
           *sh = other.matrix(); 
        }
        return *this;
    }
//  *= 
    SHarmonic& SHarmonic::operator*=(const complex<double>& d){
        (*sh) *=d;
        return *this;
    }
    SHarmonic& SHarmonic::operator*=(const SHarmonic& shmulti){
        (*sh) *= shmulti.matrix();
        return *this;
    }
//  +=
    SHarmonic& SHarmonic::operator+=(const complex<double>& d){
        (*sh) +=d;
        return *this;
    }
    SHarmonic& SHarmonic::operator+=(const SHarmonic& shadd){
        (*sh) += shadd.matrix(); 
        return *this;
    }
//  -= 
    SHarmonic& SHarmonic::operator-=(const complex<double>& d){
        (*sh) -=d;
        return *this;
    }
    SHarmonic& SHarmonic::operator-=(const SHarmonic& shmin){
        (*sh) -= shmin.matrix();
        return *this;
    }

//--------------------------------------------------------------
//   Other Algebra
//--------------------------------------------------------------
    SHarmonic& SHarmonic::mpaxis(const valarray< complex<double> >& shmulti){
        (*sh).multid1(shmulti);
        return *this;
    }
    SHarmonic& SHarmonic::mxaxis(const valarray< complex<double> >& shmulti){
        (*sh).multid2(shmulti);
        return *this;
    }
    SHarmonic& SHarmonic::myaxis(const valarray< complex<double> >& shmulti){
        (*sh).multid3(shmulti);
        return *this;
    }
//--------------------------------------------------------------
    SHarmonic& SHarmonic::mxy_matrix(Matrix2D< complex<double> >& mMatrix){
        int st(0), nxt(nump());
        for (int im(0); im < mMatrix.dim(); ++im) {
            for (int ip(st); ip < nxt; ++ip)
                (*sh)(ip) *= mMatrix(im);
            st += nump(); nxt += nump();
        }
        return *this;
    }
    SHarmonic& SHarmonic::Re(){
        for (int i(0); i < dim(); ++i) (*sh)(i) = (*sh)(i).real();
        return *this;
    }


//--------------------------------------------------------------
    SHarmonic& SHarmonic::Dp(){
//--------------------------------------------------------------
//  P-difference with derivative at #0 set to 0, and f at #np equal to #np-1  
//--------------------------------------------------------------
        Matrix2D< complex<double> > ma(numx(), numy());
        ma = 0.0;

        GSlice_iter< complex<double> > it1(p0(nump()-2)), it2(p0(nump()-1));  
        for(int i=0; i< numx()*numy(); ++i){ 
            ma(i)  = *it1; ma(i) -= *it2;
            ++it1; ++it2;
        }

        *sh = (*sh).Dd1();

        it1 = p0(0), it2 = p0(nump()-1);  
        for(int i=0; i< numx()*numy(); ++i){
            *it1 = 0.0; *it2 = ma(i);
             ++it1; ++it2;
        }
        return *this;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    SHarmonic& SHarmonic::Dx(){
//--------------------------------------------------------------
//  X-difference 
//--------------------------------------------------------------
        *sh = (*sh).Dd2();
        return *this;
    }
//--------------------------------------------------------------
    SHarmonic& SHarmonic::Dy(){
//--------------------------------------------------------------
//  Y-difference selection
//--------------------------------------------------------------
        *sh = (*sh).Dd3();
        return *this;
    } 
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//  Definition of the Distribution Function Class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
    DFunc:: DFunc(size_t l, size_t m, size_t np, size_t nx, size_t ny) 
             : lmax(l), mmax(m), ind(l+1,m+1) {

//      Initialize the array of the harmonics
        sz = ((mmax+1)*(2*lmax-mmax+2))/2;
        mf = new valarray<SHarmonic>(SHarmonic(np,nx,ny), sz);

//      Define the index for the triangular array 
        ind = -1;
        for(int il=0; il < lmax+1 ; ++il){ 
            for(int im=0; im < ((mmax < il)? mmax:il)+1; ++im){ 
               ind(il,im) = ((il < mmax+1)?((il*(il+1))/2+im):
               (il*(mmax+1)-(mmax*(mmax+1))/2 + im)); 
             }
        }   
    }

//  Copy constructor
    DFunc:: DFunc(const DFunc& other)
              : lmax(other.l0()), mmax(other.m0()),
                sz(other.dim()), ind(other.l0()+1,other.m0()+1) {

//      Initialize the array of the harmonics
        mf = new valarray<SHarmonic>(SHarmonic(other(0).nump(),other(0).numx(),other(0).numy()),
                                    other.dim());  
        for(int i=0; i < other.dim() ; ++i)  
            (*mf)[i] = other(i);

//      Define the index for the triangular array 
        ind = -1;
        for(int il=0; il < l0()+1 ; ++il){ 
            for(int im=0; im < ((m0() < il)? m0():il)+1; ++im){ 
               ind(il,im) = ((il < m0()+1)?((il*(il+1))/2+im):
               (il*(m0()+1)-(m0()*(m0()+1))/2 + im)); 
             }
        }  
    }
//  Destructor
    DFunc:: ~DFunc(){
        delete mf;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
    DFunc& DFunc::operator=(const complex<double>& d){
        for(int i=0; i < dim() ; ++i)  
            (*mf)[i] = d;
        return *this;
    }
    DFunc& DFunc::operator=(const SHarmonic& h){
        for(int i=0; i < dim() ; ++i){  
            if (&((*mf)[i]) != &h) {   //self-assignment
                (*mf)[i] = h;
            }
        }
        return *this;
    }
    DFunc& DFunc::operator=(const DFunc& other){
        if (this != &other) {   //self-assignment
            for(int i=0; i < dim() ; ++i)  
                (*mf)[i] = other(i);
        }
        return *this;
    }
//  *=
    DFunc& DFunc::operator*=(const complex<double>& d){
        for(int i=0; i < dim() ; ++i)  
            (*mf)[i] *= d;
        return *this;
    }
    DFunc& DFunc::operator*=(const DFunc& other){
        if (this != &other) {   //self-assignment
            for(int i=0; i < dim() ; ++i)  
                (*mf)[i] *= other(i);
        }
        return *this;
    }
//  +=
    DFunc& DFunc::operator+=(const complex<double>& d){
        for(int i=0; i < dim() ; ++i)  
            (*mf)[i] += d;
        return *this;
    }
    DFunc& DFunc::operator+=(const DFunc& other){
        if (this != &other) {   //self-assignment
            for(int i=0; i < dim() ; ++i)  
                (*mf)[i] += other(i);
        }
        return *this;
    }
//  -=
    DFunc& DFunc::operator-=(const complex<double>& d){
        for(int i=0; i < dim() ; ++i)  
            (*mf)[i] -= d;
        return *this;
    }
    DFunc& DFunc::operator-=(const DFunc& other){
        if (this != &other) {   //self-assignment
            for(int i=0; i < dim() ; ++i)  
                (*mf)[i] -= other(i);
        }
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//  Definition of the State Class
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
    Stat:: Stat(){
        using Inputdata::IN;
        ff = new DFunc(IN().inp().l0, IN().inp().m0, IN().inp().pr.dim(), 
                       IN().inp().x.dim(), IN().inp().y.dim());
        em = new EMF2D(IN().inp().x.dim(),IN().inp().y.dim());
    }
    Stat:: Stat( size_t l, size_t m, size_t np, size_t nx, size_t ny){
        ff = new DFunc(l, m, np, nx, ny);
        em = new EMF2D(nx, ny);
    }
//  Copy constructor
    Stat:: Stat(const Stat& other){
        ff = new DFunc(other.DF().l0(), other.DF().m0(),(other.DF()(0)).nump(), 
                         (other.DF()(0)).numx(),(other.DF()(0)).numy()); 
        *ff = other.DF();
        em = new EMF2D(other.FLD(0).numx(),other.FLD(0).numy()); 
        *em = other.EMF();
    }
//  Destructor
    Stat:: ~Stat(){
        delete ff;
        delete em;
    }
//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
    Stat& Stat::operator=(const Stat& other){
        if (this != &other) {   //self-assignment
            *ff = other.DF();
            *em = other.EMF();
        }
        return *this;
    }
//  =
    Stat& Stat::operator=(const complex<double>& d){
        *ff = d;
        *em = d;
        return *this;
    }
//  *=
    Stat& Stat::operator*=(const Stat& other){
        *ff *= other.DF();
        *em *= other.EMF();
        return *this;
    }
//  *=
    Stat& Stat::operator*=(const complex<double>& d){
        *ff *= d;
        *em *= d;
        return *this;
    }
//  +=
    Stat& Stat::operator+=(const Stat& other){
        *ff += other.DF();
        *em += other.EMF();
        return *this;
    }
//  +=
    Stat& Stat::operator+=(const complex<double>& d){
        *ff += d;
        *em += d;
        return *this;
    }
//  -=
    Stat& Stat::operator-=(const Stat& other){
        *ff -= other.DF();
        *em -= other.EMF();
        return *this;
    }
//  -=
    Stat& Stat::operator-=(const complex<double>& d){
        *ff -= d;
        *em -= d;
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************

