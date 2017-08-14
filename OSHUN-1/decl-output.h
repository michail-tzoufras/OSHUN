///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Aug 20th 2009
///////////////////////////////////////////////////////////
//   
//   This header contains the declerations for the tools
//   necessary to convert the spherical harmonics to a 3D
//   Cartesian phasespace.
///////////////////////////////////////////////////////////
//
// 
//   namespace Savedata::
//
//   1. struct Pout:: 
//      This structure contains the output momentum axis. It also 
//      provides the folowing methods:
//      a) Deposit "sqrt(p1(i)^2+p2(j)^2+p3(k)^2) in a 3D Matrix.
//      b) Deposit costh = pz/pradius in a 3D Matrix (given pradius).
//      c) Deposit arctan2(py/px) in a 3D Matrix.
//
//   2. class PLegendre::
//      A 3D space is generated from p1, p2, p3 axis. The cos8 for
//      this 3D space is calculated and then the Legendre polynomials
//      for each cos8 are calculated. We end up with a triangular l,m
//      array containing 3D matrices of the polynomials for each cos8
//
//   3. class Y_x0_p1p2p3::
//      Generates the 3D output p1 p2 p3 for the sum of the 
//      harmonics at a certain cell. Since it requires a fair
//      number of harmonics it is advisable to use a low number
//      for nump1 nump2 nump3
//
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

    #ifndef DECL_IMPL_OUT_H
    #define DECL_IMPL_OUT_H



//**************************************************************
//**************************************************************
//   Declerations for Savedata
//**************************************************************
//**************************************************************

//**************************************************************
    namespace Savedata{  
//**************************************************************

//--------------------------------------------------------------
//      Evaluation of the Legendre Polynomials for a given x
//--------------------------------------------------------------
        void Legendre(float x, Matrix2D<float>& P_Legendre);
//--------------------------------------------------------------

//--------------------------------------------------------------
        struct Pout1D {
//--------------------------------------------------------------
//      Group of output momentum Axis
//--------------------------------------------------------------
            Axis<float>  pr, px;

//          Struct constructor
            Pout1D(size_t numpx, float pxmax,
                   size_t nump,  float pmin, float pmax) :
                px(numpx,pxmax), pr(nump, pmin, pmax)/*, 
                pr(prin.dim(),prin(0),prin(prin.dim()-1)) */{} 

//          Methods
            void ppolarrad(Matrix2D<float>& pprad);
            void costheta(Matrix2D<float>& costh);
            
        };
//--------------------------------------------------------------

//--------------------------------------------------------------
        struct Pout {
//--------------------------------------------------------------
//      Group of output momentum Axis
//--------------------------------------------------------------
            Axis<float>  p1, p2, p3;

//          Struct constructor
            Pout(size_t nump1, size_t nump2, size_t nump3, float pmax) :
                p1(nump1,pmax), p2(nump2, pmax), p3(nump3,pmax) {} 

//          Methods
            void pradius(Matrix3D<float>& prad);
            void costheta(Matrix3D<float>& costh);
            void atanphi(Matrix3D<float>& atphi);
            
        };
//--------------------------------------------------------------


//--------------------------------------------------------------
        class PLegendre {
//--------------------------------------------------------------
//      A 3D space is generated from p1, p2, p3 axis. The cos8 for
//      this 3D space is calculated and then the Legendre polynomials
//      for each cos8 are calculated. We end up with a triangular l,m
//      array containing 3D matrices of the polynomials for each cos8
//--------------------------------------------------------------
        private:
            valarray< Matrix3D<float> > *plegendre;
            size_t lmax, mmax;

//          Indexing for triangular array
            Matrix2D<int> ind;

        public:
//          Constructors/Destructors
            PLegendre(size_t l, size_t m, Pout& pout);
            PLegendre(const PLegendre& other);
            ~PLegendre();

//          Access
            size_t dim()   const {return (*plegendre).size();}
            size_t l_max() const {return lmax;}
            size_t m_max() const {return mmax;}
            Matrix3D<float>& operator()(size_t i) {return (*plegendre)[i];}
            Matrix3D<float>  operator()(size_t i) const {return (*plegendre)[i];}   
            Matrix3D<float>& operator()(size_t il, size_t im) {return (*plegendre)[ind(il,im)];}

//          Operators
            PLegendre& operator=(const float& d);
            PLegendre& operator=(const Matrix3D<float>& m);
            PLegendre& operator=(const PLegendre& other);
            PLegendre& operator*=(const float& d);
            PLegendre& operator*=(const PLegendre& other);
            PLegendre& operator+=(const float& d);
            PLegendre& operator+=(const PLegendre& other);
            PLegendre& operator-=(const float& d);
            PLegendre& operator-=(const PLegendre& other);

            PLegendre& update(Pout& pout);
        };
//--------------------------------------------------------------

//--------------------------------------------------------------
        class Y_x0_p1p2p3 {
//--------------------------------------------------------------
//      This class converts the state Y at a specific spatial location
//      "x0, y0" into a 3D cartesian grid p1p2p3, which is deposited in
//      a 3D Matrix of floats           
//--------------------------------------------------------------
        public:
//          Constructor/Destructor
            Y_x0_p1p2p3(size_t l0, size_t m0, size_t nump1,  
                        size_t nump2, size_t nump3, float pmax, Axis<double>& pr);
            ~Y_x0_p1p2p3();

//          Methods
            Matrix3D<float>& Convert(Stat& Y, size_t x0, size_t y0);//, Axis<double> pr);

        private:
            size_t lmax, mmax;
            Pout axis;
            PLegendre legendre;
            Matrix3D<float> *pcart;

            Matrix3D<float> pc_data;
            Matrix3D<size_t> pind;
            Axis<float> prf;
        };
//--------------------------------------------------------------



//--------------------------------------------------------------
        class PLegendre1D {
//--------------------------------------------------------------
//      A 2D space is generated from px, |p| axis. The cos8 for
//      this 2D space is calculated and then the Legendre polynomials
//      for each cos8 are calculated. We end up with a 1D array for l
//      containing 2D matrices of the polynomials for each cos8
//--------------------------------------------------------------
        private:
            valarray< Matrix2D<float> > *plegendre;
            size_t lmax;

//          Indexing for triangular array
            //Matrix2D<int> ind;

        public:
//          Constructors/Destructors
            PLegendre1D(size_t l, /*size_t m,*/ Pout1D& pout1D);
            PLegendre1D(const PLegendre1D& other);
            ~PLegendre1D();

//          Access
            size_t dim()   const {return (*plegendre).size();}
            size_t l_max() const {return lmax;}
            //size_t m_max() const {return mmax;}
            Matrix2D<float>& operator()(size_t i) {return (*plegendre)[i];}
            Matrix2D<float>  operator()(size_t i) const {return (*plegendre)[i];}   
            //Matrix2D<float>& operator()(size_t il, size_t im) {return (*plegendre)[ind(il,im)];}

//          Operators
            PLegendre1D& operator=(const float& d);
            PLegendre1D& operator=(const Matrix2D<float>& m);
            PLegendre1D& operator=(const PLegendre1D& other);
            PLegendre1D& operator*=(const float& d);
            PLegendre1D& operator*=(const PLegendre1D& other);
            PLegendre1D& operator+=(const float& d);
            PLegendre1D& operator+=(const PLegendre1D& other);
            PLegendre1D& operator-=(const float& d);
            PLegendre1D& operator-=(const PLegendre1D& other);

            PLegendre1D& update(Pout1D& pout1D);
        };
//--------------------------------------------------------------



//--------------------------------------------------------------
        class P1x1_1D {
//--------------------------------------------------------------
//      This class converts the state Y at a specific spatial location
//      "x0, y0" into a 3D cartesian grid p1p2p3, which is deposited in
//      a 3D Matrix of floats           
//--------------------------------------------------------------
        public:
//          Constructor/Destructor
            P1x1_1D(size_t l0, size_t numpx, float pxmax, Axis<double>& pr);
                        
            ~P1x1_1D();

//          Methods
            valarray<float>& p1x1_out(Stat& Y, size_t x0, size_t y0);//, Axis<double> pr);

        private:
            size_t lmax, mmax;
            Pout1D axis;
            PLegendre1D legendre;
            valarray<float> *p1x1cart;

            Matrix2D<float> pc_costheta;
            Matrix2D<float> pc_polradius;
        };
//--------------------------------------------------------------


   }    

    #endif
