///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Jun 4, 2013
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the data
//   structures that characterize the state of the 
//   system:
//
//   1.A. SHarmonic1D :
//       A wrapper class for Array2D<...>. It is used to
//       describe a single spherical harmonic in 1D.
//
//   1.B. Field1D :
//       This is a class decleration for one of the 
//       components of the fields: Ex, Ey, Ez, Bx, By or Bz.
//
//   1.C. Distribution1D:
//       The collection of 1D spherical harmonics. 
//
//   1.D. EMF1D:
//       The collection of 1D Fields. 
//
//   1.E. State1D:
//       The collection of 1D distribution and 1D Fields. 
///////////////////////////////////////////////////////////

// Standard Libraries
   #include <iostream>
   #include <vector>
   #include <valarray>
   #include <complex>

// My Libraries
   #include "lib-array.h"

// Declerations
   #include "state.h"

//**************************************************************
//  Definition of the 1D spherical harmonic
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------

//  Constructor
    SHarmonic1D::SHarmonic1D(size_t nump, size_t numx) {
        sh = new Array2D<double>(nump,numx);
    }
//  Copy constructor
    SHarmonic1D::SHarmonic1D(const SHarmonic1D& other){
        sh = new Array2D<double>(other.nump(),other.numx());
        *sh = other.array();
    }
//  Destructor
    SHarmonic1D:: ~SHarmonic1D(){
        delete sh; 
    }

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------

//  Copy assignment operator
    SHarmonic1D& SHarmonic1D::operator=(const double& d){
        *sh = d;
        return *this;
    }
    SHarmonic1D& SHarmonic1D::operator=(const SHarmonic1D& other){
        if (this != &other) {   //self-assignment
           *sh = other.array(); 
        }
        return *this;
    }
//  *= 
    SHarmonic1D& SHarmonic1D::operator*=(const double& d){
        (*sh) *=d;
        return *this;
    }
    SHarmonic1D& SHarmonic1D::operator*=(const SHarmonic1D& shmulti){
        (*sh) *= shmulti.array();
        return *this;
    }
//  +=
    SHarmonic1D& SHarmonic1D::operator+=(const double& d){
        (*sh) +=d;
        return *this;
    }
    SHarmonic1D& SHarmonic1D::operator+=(const SHarmonic1D& shadd){
        (*sh) += shadd.array(); 
        return *this;
    }
//  -= 
    SHarmonic1D& SHarmonic1D::operator-=(const double& d){
        (*sh) -=d;
        return *this;
    }
    SHarmonic1D& SHarmonic1D::operator-=(const SHarmonic1D& shmin){
        (*sh) -= shmin.array();
        return *this;
    }

//--------------------------------------------------------------
//   Other Algebra
//--------------------------------------------------------------
    SHarmonic1D& SHarmonic1D::mpaxis(const valarray<double>& shmulti){
        (*sh).multid1(shmulti);
        return *this;
    }
    SHarmonic1D& SHarmonic1D::mxaxis(const valarray<double>& shmulti){
        (*sh).multid2(shmulti);
        return *this;
    }
//--------------------------------------------------------------

//  P-difference
    SHarmonic1D& SHarmonic1D::Dp(){

        valarray<double> plast(this->numx());

        for (size_t i(0); i < plast.size(); ++i) {
            plast[i] = (*sh)(nump()-2,i) - (*sh)(nump()-1,i); 
        }
        *sh = (*sh).Dd1();
        for (size_t i(0); i < plast.size(); ++i) {
           // TODO                The Dp at the zeroth cell is taken care off
           // (*sh)(0,i) = 0.0;   separately, both for the E-field and the collisions.                     
            (*sh)(nump()-1,i) = plast[i]; 
        }
             
        return *this;
    }
//--------------------------------------------------------------

//  X-difference 
    SHarmonic1D& SHarmonic1D::Dx(){

        *sh = (*sh).Dd2();           				// Worry about boundaries elsewhere
        return *this;
    }
//--------------------------------------------------------------

//  Filter Pcells
    SHarmonic1D& SHarmonic1D::Filterp(size_t N){
        *sh = (*sh).Filterd1(N);
        return *this;
    }

//**************************************************************

//**************************************************************
//  Definition of the "Field1D" Class
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
    Field1D:: Field1D(size_t numx) {
        fi = new  valarray<double>(numx);
    }
//  Copy constructor
    Field1D:: Field1D(const Field1D& other){
        fi = new valarray<double>(other.numx());
        *fi = other.array();
    }
//  Destructor
    Field1D:: ~Field1D(){
        delete fi; 
    }

//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
    Field1D& Field1D::operator=(const double& d){
        *fi = d;
        return *this;
    }
    Field1D& Field1D::operator=(const valarray<double>& other){
        *fi = other; 
        return *this;
    }
    Field1D& Field1D::operator=(const Field1D& other){
        if (this != &other) {   //self-assignment
           *fi = other.array(); 
        }
        return *this;
    }
//  *= 
    Field1D& Field1D::operator*=(const double& d){
        (*fi) *=d;
        return *this;
    }
    Field1D& Field1D::operator*=(const valarray<double>& fimulti){
        (*fi) *= fimulti;
        return *this;
    }
    Field1D& Field1D::operator*=(const Field1D& fimulti){
        (*fi) *= fimulti.array();
        return *this;
    }
//  +=
    Field1D& Field1D::operator+=(const double& d){
        (*fi) +=d;
        return *this;
    }
    Field1D& Field1D::operator+=(const Field1D& fiadd){
        (*fi) += fiadd.array(); 
        return *this;
    }
//  -= 
    Field1D& Field1D::operator-=(const double& d){
        (*fi) -=d;
        return *this;
    }
    Field1D& Field1D::operator-=(const Field1D& fimin){
        (*fi) -= fimin.array();
        return *this;
    }

//--------------------------------------------------------------
    Field1D& Field1D::Dx(){
//--------------------------------------------------------------
        for(long i(0); i< long(numx())-2; ++i) {
            (*fi)[i] -= (*fi)[i+2]; 
        }
        for(long i(numx()-3); i>-1; --i) {
            (*fi)[i+1] = (*fi)[i]; 
        }
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//  Definition of the 1D distribution function
//--------------------------------------------------------------
//  Constructors and Destructor
//--------------------------------------------------------------
//  Constructor
    DistFunc1D:: DistFunc1D(size_t l, size_t np, size_t nx, double q=1, double _ma=1) 
             : lmax(l), charge(q), ma(_ma) {

//      Initialize the array of the harmonics
        if (lmax < 1) {
            cout << "l0 < 1 is not acceptable.\n";
            exit(1);
        }
 
//      Generate container for the harmonics
        df = new vector<SHarmonic1D>(lmax+1,SHarmonic1D(np,nx)); 

//      list the neighbors of each harmonic
        _Neighbors.push_back(NULL);                
        _Neighbors.push_back(&( (*df)[1] ));
        for(size_t l(1); l < lmax ; ++l){  
            _Neighbors.push_back(&( (*df)[l-1] )); 
            _Neighbors.push_back(&( (*df)[l+1] )); 
        }
        _Neighbors.push_back(&( (*df)[lmax-1] ));
        _Neighbors.push_back(NULL);

    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  Copy constructor
    DistFunc1D:: DistFunc1D(const DistFunc1D& other)
              : lmax(other.l0()), charge(other.q()), ma(other.mass()) {

//      Generate container for the harmonics
        df = new vector<SHarmonic1D>;
        for(size_t l(0); l < other.l0()+1 ; ++l){  
            (*df).push_back(*(other(l)));
        }

//      list the neighbors of each harmonic
        _Neighbors.push_back(NULL);                
        _Neighbors.push_back(&( (*df)[1] ));
        for(size_t l(1); l < lmax ; ++l){  
            _Neighbors.push_back(&( (*df)[l-1] )); 
            _Neighbors.push_back(&( (*df)[l+1] )); 
        }
        _Neighbors.push_back(&( (*df)[lmax-1] ));
        _Neighbors.push_back(NULL);

    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//  Destructor
    DistFunc1D:: ~DistFunc1D(){
        delete df;
    }
//--------------------------------------------------------------
//  Access 
//--------------------------------------------------------------
//  Pointer to the l-th harmonic
    SHarmonic1D* DistFunc1D::operator()(int l) {   
        if ((l < 0) || (l> lmax)) return NULL;
        return &((*df)[size_t(l)]);
    }      

    SHarmonic1D* DistFunc1D::operator()(int l) const {
        if ((l < 0) || (l> lmax)) return NULL;
        return &((*df)[size_t(l)]);
    }

//  Pointer to the "n" neighbor of the l harmonic
    SHarmonic1D* DistFunc1D::Compus(size_t l, _Compus1D n) const {
        return _Neighbors[2*l+n];
    }
//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
    DistFunc1D& DistFunc1D::operator=(const double& d){
        for(size_t i(0); i < dim() ; ++i){  
            (*df)[i] = d;
        }
        return *this;
    }
    DistFunc1D& DistFunc1D::operator=(const SHarmonic1D& h){
        for(size_t i(0); i < dim() ; ++i){  
            if (&((*df)[i]) != &h) {   //self-assignment
                (*df)[i] = h;
            }
        }
        return *this;
    }
    DistFunc1D& DistFunc1D::operator=(const DistFunc1D& other){
        if (this != &other) {   //self-assignment
            for(size_t i(0); i < dim() ; ++i) {  
                (*df)[i] = *(other(i));
            }
        }
        return *this;
    }
//  *=
    DistFunc1D& DistFunc1D::operator*=(const double& d){
        for(size_t i(0); i < dim() ; ++i) { 
            (*df)[i] *= d;
        }
        return *this;
    }
    DistFunc1D& DistFunc1D::operator*=(const DistFunc1D& other){
        if (this != &other) {   //self-assignment
            for(size_t i(0); i < dim() ; ++i) { 
                (*df)[i] *= *(other(i));
            }
        }
        return *this;
    }
//  +=
    DistFunc1D& DistFunc1D::operator+=(const double& d){
        for(size_t i(0); i < dim() ; ++i) { 
            (*df)[i] += d;
        }
        return *this;
    }
    DistFunc1D& DistFunc1D::operator+=(const DistFunc1D& other){
        if (this != &other) {   //self-assignment
            for(size_t i(0); i < dim() ; ++i) {  
                (*df)[i] += *(other(i));
            }
        }
        return *this;
    }
//  -=
    DistFunc1D& DistFunc1D::operator-=(const double& d){
        for(size_t i(0); i < dim() ; ++i) { 
            (*df)[i] -= d;
        }
        return *this;
    }
    DistFunc1D& DistFunc1D::operator-=(const DistFunc1D& other){
        if (this != &other) {   //self-assignment
            for(size_t i(0); i < dim() ; ++i) { 
                (*df)[i] -= *(other(i));
            }
        }
        return *this;
    }

    DistFunc1D& DistFunc1D::Filterp(){
        for(size_t i(0); i < dim() ; ++i) { 
            (*df)[i].Filterp(i);
        }
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//  State for 1D electrostatic code
//--------------------------------------------------------------
//  Constructor and Destructor
//--------------------------------------------------------------
//  Constructor
//    State1D:: State1D(){
//        using Inputdata::IN;
//        ff = new DFunc(IN().inp().l0, IN().inp().m0, IN().inp().pr.dim(), 
//                       IN().inp().x.dim(), IN().inp().y.dim());
//        em = new EMF2D(IN().inp().x.dim(),IN().inp().y.dim());
//    }
    State1D:: State1D( size_t nx, vector<size_t> l0, vector<size_t> np, vector<double> q, vector<double> ma)
	     : ns(l0.size()) {
        if (ns != np.size()) {
            cout << "ERROR:Overdetermined number of species\n";
            exit(1);
        }
        sp = new vector<DistFunc1D>;
        for(size_t s(0); s < ns; ++s){  
            (*sp).push_back(DistFunc1D(l0[s],np[s],nx,q[s],ma[s]));
        }
        ex = new Field1D(nx);
    }
//  Copy constructor
    State1D:: State1D(const State1D& other)
	     : ns(other.Species()) {
        sp = new vector<DistFunc1D>(); 
        for(size_t s(0); s < ns; ++s){  
            (*sp).push_back(DistFunc1D(other.DF(s)));
        }
        ex = new Field1D(other.Ex()); 
    }
//  Destructor
    State1D:: ~State1D(){
        delete sp;
        delete ex;
    }
//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------
//  Copy assignment operator
    State1D& State1D::operator=(const State1D& other){
        if (this != &other) {   //self-assignment
            ns = other.Species();
            for(size_t s(0); s < ns; ++s){  
                (*sp)[s] = other.DF(s);
            }
            *ex = other.Ex();
        }
        return *this;
    }
//  =
    State1D& State1D::operator=(const double& d){
        for(size_t s(0); s < ns; ++s){  
            (*sp)[s] = d;
        }
        *ex = d;
        return *this;
    }
//  *=
    State1D& State1D::operator*=(const State1D& other){
        for(size_t s(0); s < ns; ++s){  
            (*sp)[s] *= other.DF(s);
        }
        *ex *= other.Ex();
        return *this;
    }
//  *=
    State1D& State1D::operator*=(const double& d){
        for(size_t s(0); s < ns; ++s){  
            (*sp)[s] *= d;
        }
        (*ex) *= d;
        return *this;
    }
//  +=
    State1D& State1D::operator+=(const State1D& other){
        for(size_t s(0); s < ns; ++s){  
            (*sp)[s] += other.DF(s);
        }
        *ex += other.Ex();
        return *this;
    }
//  +=
    State1D& State1D::operator+=(const double& d){
        for(size_t s(0); s < ns; ++s){  
            (*sp)[s] += d;
        }
        (*ex) += d;
        return *this;
    }
//  +=
    State1D& State1D::operator-=(const State1D& other){
        for(size_t s(0); s < ns; ++s){  
            (*sp)[s] -= other.DF(s);
        }
        *ex -= other.Ex();
        return *this;
    }
//  +=
    State1D& State1D::operator-=(const double& d){
        for(size_t s(0); s < ns; ++s){  
            (*sp)[s] -= d;
        }
        (*ex) -= d;
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************
