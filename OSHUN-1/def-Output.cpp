///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Aug 19th 2009
///////////////////////////////////////////////////////////

//   Definition of the necessary tools for converting the 
//   spherical harmonics to a 3D cartesian phasespace   
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

//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <fstream>
    #include <iomanip>
    #include <cstdlib>
    #include <sstream>
    #include <string>
    #include <math.h>

//  My libraries
    #include "matrices.h"

//  Declerations
    #include "decl-input.h"
    #include "decl-state.h"
    #include "decl-output.h"



//**************************************************************
//--------------------------------------------------------------
    void Savedata:: Legendre(float x, Matrix2D<float>& P_Legendre){
//--------------------------------------------------------------
//  Evaluate the Legendre polynomials for given x and put them
//  In the P_Legendre Matrix
//--------------------------------------------------------------

//      Local variables
        float r1, sqrtx = sqrt(1.0-x*x), fact = 1.0;
        size_t l0 = P_Legendre.dim1(), m0 = P_Legendre.dim2();
 
//      Initialization 
        P_Legendre = 0.0;
        P_Legendre(0,0) = 1.0;

//      Executable statements
        for (size_t l = 1; l < m0; ++l){
            P_Legendre(l,l) = - P_Legendre(l-1,l-1)*(fact*sqrtx);
            fact += 2.0;
        }

        for (size_t l = 0; l < ((m0 < l0-1) ? m0 : (l0-1)); ++l)
            P_Legendre(l+1,l) = P_Legendre(l,l)*(x*(2.0*l+1.0));

        for (size_t m = 0; m < m0; ++m){
            for (size_t l = m+1; l < l0 - 1; ++l){
                r1 = 1.0 / float(l - m + 1);
                P_Legendre(l+1,m) = P_Legendre(l,m) * (x*(2.0*l+1.0) * r1) -
                                    P_Legendre(l-1,m)*(float(l+m) * r1); 
            }
         }
                
    }               
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//**************************************************************
//   Definition of the output momentum class 
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    void Savedata::Pout1D:: ppolarrad(Matrix2D<float>& pprad){
//--------------------------------------------------------------
//  Calculate the radius for the 3D grid  
//--------------------------------------------------------------
         
//      Local Variables
        float px_sq, pr_sq;

//      Executable Statements
        if ((pprad.dim1() == px.dim()) && 
            (pprad.dim2() == pr.dim())) {
      
            for (size_t i(0); i < px.dim(); ++i){
                px_sq  = px(i) * px(i);
                for (size_t j(0); j < pr.dim(); ++j){
                    pr_sq  = pr(j) * pr(j);
                    pprad(i,j) = pr_sq - px_sq;
                    if (pprad(i,j) < 0.0) {
                        pprad(i,j) = -1.0;
                    }
                    else {
                        pprad(i,j)  = sqrt(abs(pprad(i,j)));
                    }
                }      
            }
        }  
        else { 
            std::cout << "Error Matrix2D and piAxis don't have the same dimensions\n";
    }
}
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Savedata::Pout1D::costheta(Matrix2D<float>& costh){
//--------------------------------------------------------------
//  Calculate the cosine
//--------------------------------------------------------------

//      Executable Statements
        if ((costh.dim1() == px.dim()) && 
            (costh.dim2() == pr.dim())) {
      
            for (size_t i(0); i < px.dim(); ++i){
                for (size_t j(0); j < pr.dim(); ++j){
                    costh(i,j) = px(i)/pr(j);
                }      
            }
               
        }  
        else { 
            std::cout << "Error Matrix2D and piAxis don't have the same dimensions\n";}
    }            
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
    void Savedata::Pout:: pradius(Matrix3D<float>& prad){
//--------------------------------------------------------------
//  Calculate the radius for the 3D grid  
//--------------------------------------------------------------
         
//      Local Variables
        float pz_sq, pzpy_sq;

//      Executable Statements
        if ((prad.dim1() == p1.dim()) && 
           (prad.dim2() == p2.dim()) &&
           (prad.dim3() == p3.dim())) {  
      
            for (size_t k = 0; k < p3.dim(); ++k){
                pz_sq  = p3(k) * p3(k);
                for (size_t j = 0; j < p2.dim(); ++j){
                    pzpy_sq  = p2(j) * p2(j);
                    pzpy_sq += pz_sq;  
                    for (size_t i = 0; i < p1.dim(); ++i)
                        prad(i,j,k) = sqrt(pzpy_sq + p1(i)* p1(i));
                }
            }      
        }  
        else { 
            std::cout << "Error Matrix3D and piAxis don't have the same dimensions\n";}
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Savedata::Pout:: costheta(Matrix3D<float>& costh){
//--------------------------------------------------------------
//  Calculate the cosine
//--------------------------------------------------------------
        if ((costh.dim1() == p1.dim()) && 
           (costh.dim2() == p2.dim()) &&
           (costh.dim3() == p3.dim())) {  
      
            for (size_t i = 0; i < p1.dim(); ++i)
                for (size_t j = 0; j < p2.dim(); ++j){
                    for (size_t k = 0; k < p3.dim(); ++k){
                        costh(i,j,k) = p1(i)/costh(i,j,k);
                
                }
            }      
        }  
        else { 
            std::cout << "Error Matrix3D and piAxis don't have the same dimensions\n";}
    }            
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Savedata::Pout:: atanphi(Matrix3D<float>& atphi){
//--------------------------------------------------------------
//  Calculate the arctan2
//--------------------------------------------------------------
        float atp;

        if ((atphi.dim1() == p1.dim()) && 
        (atphi.dim2() == p2.dim()) &&
        (atphi.dim3() == p3.dim())) {  
            for (size_t j(0); j < p2.dim(); ++j){
                for (size_t k(0); k < p3.dim(); ++k){
                     atp = M_PI+atan2(p3(k),p2(j));
                     for (size_t i(0); i < p1.dim(); ++i) atphi(i,j,k) = atp; 
                }
            }
        }
        else { 
            std::cout << "Error Matrix3D and piAxis don't have the same dimensions\n";}
    }            
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//**************************************************************
//   Definition of the PLegendre1D Class
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
    Savedata:: PLegendre1D:: PLegendre1D(size_t l, Pout1D& pout1D)
        : lmax(l) {
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     
//      Create a pointer to a container of 2D Matrices of float
        plegendre = new valarray< Matrix2D<float> >(Matrix2D<float>(pout1D.px.dim(),pout1D.pr.dim()),lmax+1);

//      Calculate the radius and the radius and the cosine for these p1,p2,p3
        Matrix2D<float> pc_cosines(pout1D.px.dim(),pout1D.pr.dim());
        pout1D.costheta(pc_cosines);

//      Calculate the legendre polynomials for all harmonics for all spatial locations
        Matrix2D<float> legend(lmax+1,2); // where the "2" is a dummy value for m 
        for (size_t i(0); i < pc_cosines.dim(); ++i) {
            Savedata :: Legendre(pc_cosines(i),legend);
            for (size_t l(0); l < lmax+1; ++l){
                ((*plegendre)[l])(i) = legend(l,0); 
                // cout << "L" << l<< "("<< pc_cosines(i) << ")  =  " << ((*plegendre)[l])(i) << "\n";
            }
        }
    }

//--------------------------------------------------------------
//  Copy constructor
//--------------------------------------------------------------
    Savedata:: PLegendre1D:: PLegendre1D(const PLegendre1D& other)
                : lmax(other.l_max()) {
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//      Initialize the container of the 3D Matrices
        plegendre = new valarray< Matrix2D<float> >(Matrix2D<float>(other(0).dim1(),
                          other(0).dim2()), other.dim()); 

        for(size_t i=0; i < other.dim() ; ++i)  
            (*plegendre)[i] = other(i);

     }

//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    Savedata:: PLegendre1D:: ~PLegendre1D(){
        delete plegendre; // automagically calls f->valarray<Harmonic>
    }


//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------

//  Copy assignment operator
    Savedata::PLegendre1D& Savedata::PLegendre1D::operator=(const float& d){
        for(size_t i=0; i < dim() ; ++i)  
            (*plegendre)[i] = d;
        return *this;
    }
    Savedata::PLegendre1D& Savedata::PLegendre1D::operator=(const Matrix2D<float>& m){
        for(size_t i=0; i < dim() ; ++i){  
            if (&((*plegendre)[i]) != &m) {   //self-assignment
                (*plegendre)[i] = m;
            }
        }
        return *this;
    }
    Savedata::PLegendre1D& Savedata::PLegendre1D::operator=(const PLegendre1D& other){
        if (this != &other) {   //self-assignment
            for(size_t i=0; i < dim() ; ++i)  
                (*plegendre)[i] = other(i);
        }
        return *this;
    }


//  *=
    Savedata::PLegendre1D& Savedata::PLegendre1D::operator*=(const float& d){
        for(size_t i=0; i < dim() ; ++i)  
            (*plegendre)[i] *= d;
        return *this;
    }
    Savedata::PLegendre1D& Savedata::PLegendre1D::operator*=(const PLegendre1D& other){
        if (this != &other) {   //self-assignment
            for(size_t i=0; i < dim() ; ++i)  
                (*plegendre)[i] *= other(i);
        }
        return *this;
    }

//  +=
    Savedata::PLegendre1D& Savedata::PLegendre1D::operator+=(const float& d){
        for(size_t i=0; i < dim() ; ++i)  
            (*plegendre)[i] += d;
        return *this;
    }
    Savedata::PLegendre1D& Savedata::PLegendre1D::operator+=(const PLegendre1D& other){
        if (this != &other) {   //self-assignment
            for(size_t i=0; i < dim() ; ++i)  
                (*plegendre)[i] += other(i);
        }
        return *this;
    }

//  -=
    Savedata::PLegendre1D& Savedata::PLegendre1D::operator-=(const float& d){
        for(size_t i=0; i < dim() ; ++i)  
            (*plegendre)[i] -= d;
        return *this;
    }
    Savedata::PLegendre1D& Savedata::PLegendre1D::operator-=(const PLegendre1D& other){
        if (this != &other) {   //self-assignment
            for(size_t i=0; i < dim() ; ++i)  
                (*plegendre)[i] -= other(i);
        }
        return *this;
    }


//--------------------------------------------------------------
    Savedata::PLegendre1D& Savedata::PLegendre1D::update(Pout1D& pout1D){
    
//      Calculate the radius and the radius and the cosine for these p1,p2,p3
        Matrix2D<float> pc_cosines(pout1D.px.dim(),pout1D.pr.dim());
        pout1D.costheta(pc_cosines);

//      Calculate the legendre polynomials for all harmonics for all spatial locations
        Matrix2D<float> legend(lmax+1,2); // where the "2" is a dummy value for m 
        for (size_t i(0); i < pc_cosines.dim(); ++i) {
            Savedata :: Legendre(pc_cosines(i),legend);
            for (size_t l(0); l < lmax+1; ++l){
                ((*plegendre)[l])(i) = legend(l,0); 
            }
        }

        return *this;
    }
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//**************************************************************
//   Definition of the PLegendre Class
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
    Savedata:: PLegendre:: PLegendre(size_t l, size_t m, Pout& pout)
        : lmax(l), mmax(m), ind(l+1,m+1) {
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     
//      Create a pointer to a container of 3D Matrices of float
        plegendre = new valarray< Matrix3D<float> >(Matrix3D<float>(pout.p1.dim(),pout.p2.dim(),
                                    pout.p3.dim()), ((mmax+1)*(2*lmax-mmax+2))/2);  

//      Create the index array for the triangular array
        ind = -1;
        for(size_t il=0; il < lmax+1 ; ++il){ 
            for(size_t im=0; im < ((mmax < il)? mmax:il)+1; ++im){ 
               ind(il,im) = ((il < mmax+1)?((il*(il+1))/2+im):
               (il*(mmax+1)-(mmax*(mmax+1))/2 + im)); 
             }
        } 

//      Calculate the radius and the radius and the cosine for these p1,p2,p3
        Matrix3D<float> pc_data(pout.p1.dim(),pout.p2.dim(),pout.p3.dim());
        pout.pradius(pc_data);
        pout.costheta(pc_data);

//      Calculate the legendre polynomials for all harmonics for all spatial locations
        Matrix2D<float> legend(lmax+1,mmax+1);
        for (size_t i=0; i < pc_data.dim(); ++i) {
            size_t j=0;
            Savedata :: Legendre(pc_data(i),legend);
            for (size_t l=0; l < lmax+1; ++l){
                for (size_t m=0; m < ((mmax < l) ? mmax : l)+1; ++m){
                    ((*plegendre)[j])(i) = legend(l,m); 
                    ++j;
                }
            }
        }
    }

//--------------------------------------------------------------
//  Copy constructor
//--------------------------------------------------------------
    Savedata:: PLegendre:: PLegendre(const PLegendre& other)
                : lmax(other.l_max()), mmax(other.m_max()),
                  ind(other.l_max()+1,other.m_max()+1) {
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//      Initialize the container of the 3D Matrices
        plegendre = new valarray< Matrix3D<float> >(Matrix3D<float>(other(0).dim1(),
                          other(0).dim2(),other(0).dim3()), other.dim()); 

        for(size_t i=0; i < other.dim() ; ++i)  
            (*plegendre)[i] = other(i);

//      Define the index for the triangular array 
        ind = -1;
        for(size_t il=0; il < lmax+1 ; ++il){ 
            for(size_t im=0; im < ((mmax < il)? mmax:il)+1; ++im){ 
               ind(il,im) = ((il < mmax+1)?((il*(il+1))/2+im):
               (il*(mmax+1)-(mmax*(mmax+1))/2 + im)); 
             }
        }
     }

//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    Savedata:: PLegendre:: ~PLegendre(){
        delete plegendre; // automagically calls f->valarray<Harmonic>
    }


//--------------------------------------------------------------
//  Operators
//--------------------------------------------------------------

//  Copy assignment operator
    Savedata::PLegendre& Savedata::PLegendre::operator=(const float& d){
        for(size_t i=0; i < dim() ; ++i)  
            (*plegendre)[i] = d;
        return *this;
    }
    Savedata::PLegendre& Savedata::PLegendre::operator=(const Matrix3D<float>& m){
        for(size_t i=0; i < dim() ; ++i){  
            if (&((*plegendre)[i]) != &m) {   //self-assignment
                (*plegendre)[i] = m;
            }
        }
        return *this;
    }
    Savedata::PLegendre& Savedata::PLegendre::operator=(const PLegendre& other){
        if (this != &other) {   //self-assignment
            for(size_t i=0; i < dim() ; ++i)  
                (*plegendre)[i] = other(i);
        }
        return *this;
    }


//  *=
    Savedata::PLegendre& Savedata::PLegendre::operator*=(const float& d){
        for(size_t i=0; i < dim() ; ++i)  
            (*plegendre)[i] *= d;
        return *this;
    }
    Savedata::PLegendre& Savedata::PLegendre::operator*=(const PLegendre& other){
        if (this != &other) {   //self-assignment
            for(size_t i=0; i < dim() ; ++i)  
                (*plegendre)[i] *= other(i);
        }
        return *this;
    }

//  +=
    Savedata::PLegendre& Savedata::PLegendre::operator+=(const float& d){
        for(size_t i=0; i < dim() ; ++i)  
            (*plegendre)[i] += d;
        return *this;
    }
    Savedata::PLegendre& Savedata::PLegendre::operator+=(const PLegendre& other){
        if (this != &other) {   //self-assignment
            for(size_t i=0; i < dim() ; ++i)  
                (*plegendre)[i] += other(i);
        }
        return *this;
    }

//  -=
    Savedata::PLegendre& Savedata::PLegendre::operator-=(const float& d){
        for(size_t i=0; i < dim() ; ++i)  
            (*plegendre)[i] -= d;
        return *this;
    }
    Savedata::PLegendre& Savedata::PLegendre::operator-=(const PLegendre& other){
        if (this != &other) {   //self-assignment
            for(size_t i=0; i < dim() ; ++i)  
                (*plegendre)[i] -= other(i);
        }
        return *this;
    }


//--------------------------------------------------------------
    Savedata::PLegendre& Savedata::PLegendre::update(Pout& pout){
    
        Matrix3D<float> pc_data(pout.p1.dim(),pout.p2.dim(),pout.p3.dim());
        pout.pradius(pc_data);
        pout.costheta(pc_data);
        
        Matrix2D<float> legend(lmax,mmax);
        for (size_t i=0; i < pc_data.dim(); ++i) {
            Savedata:: Legendre(pc_data(i),legend);
            for (size_t l=0; l < lmax+1; ++l){
                for (size_t m=0; m < ((mmax < l) ? mmax : l)+1; ++m)
                    (*plegendre)[ind(l,m)](i) = legend(l,m); 
            }
        }
        return *this;
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//**************************************************************
//   Definition of the distribution function output class
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
//  Constructor 
//--------------------------------------------------------------
    Savedata::Y_x0_p1p2p3:: Y_x0_p1p2p3( size_t l0, size_t m0, size_t nump1, 
                                         size_t nump2, size_t nump3, float pmax,
                                         Axis<double>& pr)
         : lmax(l0), mmax(m0), 
           axis(nump1, nump2, nump3, pmax),
           legendre(lmax, mmax, axis),
           pind(nump1,nump2,nump3), 
           pc_data(nump1,nump2,nump3),
           prf(pr.dim(),float(pr(0)),float(pr(pr.dim()-1))){

        pcart = new  Matrix3D<float>(nump1,nump2,nump3);

        axis.pradius(pc_data);
        pc_data -= prf(0);
        pc_data *= (1.0/(prf.dx()));
        pc_data += 0.5;
        for (size_t i(0); i < pind.dim(); ++i) pind(i) = static_cast<size_t>(pc_data(i));

        axis.atanphi(pc_data);
    }

//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    Savedata::Y_x0_p1p2p3:: ~Y_x0_p1p2p3(){
        delete pcart; 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Matrix3D<float>&  Savedata::Y_x0_p1p2p3:: Convert(Stat& Y, size_t x0, size_t y0){
//--------------------------------------------------------------
//  Turn the Distribution function at some spatial location (x0,y0) 
//  into a cartesian grid.
//--------------------------------------------------------------

        float YSH_re(0.0), YSH_im(0.0), mphi(0.0);
        float calcos(0.0), calsin(0.0);
        size_t pk(0);
        *pcart = 0; 

        for (size_t k=0; k < pind.dim(); ++k) {
            pk = pind(k);
            if (pk < prf.dim()){ 
                for(size_t l=0; l < lmax+1; ++l){
                    (*pcart)(k) += (static_cast<float>((Y.SH(l,0)(pk,x0,y0)).real())) * legendre(l,0)(k);
                }
                for (size_t m(1); m < mmax+1; ++m){
                    mphi = m*pc_data(k);
                    calcos = cos(mphi);
                    calsin = sin(mphi);
                    for (size_t l(m); l < lmax+1; ++l){
                        YSH_re = static_cast<float>((Y.SH(l,m)(pk,x0,y0)).real());
                        YSH_im = static_cast<float>((Y.SH(l,m)(pk,x0,y0)).imag());
                        YSH_re *= calcos; 
                        YSH_im *= calsin; 
                        YSH_re -= YSH_im;
                        (*pcart)(k) += 2.0*(legendre(l,m)(k))*YSH_re;
                    }
                 }
              }
          }
           /*     for(size_t l=0; l < lmax+1; ++l){
                    YSH_re = static_cast<float>((Y.SH(l,0)(pind(k),x0,y0)).real());
                    (*pcart)(k) += YSH_re * legendre(l,0)(k);
                    for (size_t m=1; m < ((mmax < l)? mmax:l)+1; ++m){            
                        mphi = m*pc_data(k);
                        YSH_re = static_cast<float>((Y.SH(l,m)(pind(k),x0,y0)).real());
                        YSH_im = static_cast<float>((Y.SH(l,m)(pind(k),x0,y0)).imag());
                        YSH_re *= cosf(mphi); YSH_im *= sinf(mphi); YSH_re -= YSH_im;
                        (*pcart)(k) += (legendre(l,m)(k))*YSH_re;
                    }
                }  
            }   
        }    */

        return *pcart;
   } 
//--------------------------------------------------------------



//**************************************************************
//**************************************************************
//   Definition of the p1x1 output class
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
//  Constructor 
//--------------------------------------------------------------
    Savedata::P1x1_1D:: P1x1_1D(size_t l0, size_t numpx, 
                        float pxmax, Axis<double>& pr)  
         : lmax(l0),
           axis(numpx, pxmax,  pr.dim(),float(pr(0)),float(pr(pr.dim()-1))),
           legendre(lmax, axis),
           pc_costheta(numpx,pr.dim()),
           pc_polradius(numpx,pr.dim()){
      //   cout << "Numpx = " << numpx << "\n";
      //   cout << "Pxmax = " << pxmax << "\n";

        axis.ppolarrad(pc_polradius);

      //  cout << "             ";
      //   for (size_t i(0); i < axis.px.dim(); ++i) cout << axis.px(i)<< "     "; cout <<"\n";
      //  for (size_t j(0); j < pc_polradius.dim2(); ++j) {
      //      cout << axis.pr(j) << "   |  ";
      //      for (size_t i(0); i < pc_polradius.dim1(); ++i) {
      //          cout<< pc_polradius(i,j)<< "          ";
      //      }
      //      cout << "\n";
      //  }

      //  cout << "\n\n\n\n\n";

        axis.costheta(pc_costheta);
      //  cout << "             ";
      //  for (size_t i(0); i < axis.px.dim(); ++i) cout << axis.px(i)<< "     "; cout <<"\n";
      //  for (size_t j(0); j < pc_costheta.dim2(); ++j) {
      //      cout << axis.pr(j) << "   |  ";
      //      for (size_t i(0); i < pc_costheta.dim1(); ++i) {
      //          cout<< pc_costheta(i,j)<< "          ";
      //      }
      //      cout << "\n";
      //  }

        p1x1cart = new valarray<float>(numpx);        

        //for (size_t i(0); i < pc_costheta.dim1(); ++i) {
        //   (*p1x1cart)[i] = pc_costheta(i,31); 
        //   cout << pc_costheta(i,31) << ",     " << (*p1x1cart)[i] << "\n";
        //}

        /*axis.pradius(pc_data);
        pc_data -= prf(0);
        pc_data *= (1.0/(prf.dx()));
        pc_data += 0.5;
        for (size_t i(0); i < pind.dim(); ++i) pind(i) = static_cast<size_t>(pc_data(i));

        axis.atanphi(pc_data);*/
    }

//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    Savedata::P1x1_1D:: ~P1x1_1D(){
        delete p1x1cart; 
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    valarray<float>&  Savedata::P1x1_1D:: p1x1_out(Stat& Y, size_t x0, size_t y0){
//--------------------------------------------------------------
//  Turn the Distribution function at some spatial location (x0,y0) 
//  into a cartesian grid.
//--------------------------------------------------------------

	*p1x1cart = 0.0;

        float p_im1(0.0), p_ip1(0.0);
        float integrant_low(0.0);
        float integrant_high(0.0);
       
        for (size_t l(0); l < Y.DF().l0()+1; ++l) { // calculate the integral for each harmonic separately

             valarray<float> InPx(axis.px.dim());
             for (size_t ipx(0); ipx < axis.px.dim()-1; ++ipx) { // at each location in px

                 size_t ip(0);
                 while (pc_polradius(ipx,ip) < 0 ) { ++ip; } // at each location in pr

                 p_ip1 = pc_polradius(ipx,ip);
                 integrant_high =(static_cast<float>((Y.SH(l,0)(ip,x0,y0)).real())) * legendre(l)(ipx, ip);

                 InPx[ipx] += p_ip1 * (0.5*p_ip1) * integrant_high; 
                 // cout << "p("<<ip<<") = " << p_ip1 << ",    Dp = " << p_ip1-p_im1 << ",     " << axis.px(ipx) << " < " << axis.pr(ip) << "\n";
                 // cout << "L_"<< l<< "("<< axis.px(ipx)/axis.pr(ip)<< ") = " << legendre(l)(ipx, ip)<< "\n";
                 ++ip;
                 while ( (ip <  axis.pr.dim()) && (pc_polradius(ipx,ip) >0)  )   {  // at each location in pr
                     p_im1 = p_ip1;
                     integrant_low = integrant_high;

                     p_ip1 = pc_polradius(ipx,ip);
                     integrant_high =(static_cast<float>((Y.SH(l,0)(ip,x0,y0)).real())) * legendre(l)(ipx, ip);
                     InPx[ipx] += p_im1 * (0.5*(p_ip1-p_im1)) * integrant_low; 
                     InPx[ipx] += p_ip1 * (0.5*(p_ip1-p_im1)) * integrant_high; 

                     //cout << "p("<<ip<<") = " << p_ip1 << ",    Dp = " << p_ip1-p_im1 << ",     " << axis.px(ipx) << " < " << axis.pr(ip) << "\n";
                     //cout << "L_"<< l<< "("<< axis.px(ipx)/axis.pr(ip)<< ") = " << legendre(l)(ipx, ip)<< "\n";
                     //if ( pc_polradius(ipx,ip)
                     // cout << "Dp("<< ip <<") = " << p_ip1 - p_im1<< "\n";
                     ++ip;
                 }
                 (*p1x1cart)[ipx] += InPx[ipx];
             }
        }        
        (*p1x1cart) *= 2.0 * M_PI;
   
        /*float YSH_re(0.0), YSH_im(0.0), mphi(0.0);
        float calcos(0.0), calsin(0.0);
        size_t pk(0);
        *pcart = 0; 

        for (size_t k=0; k < pind.dim(); ++k) {
            pk = pind(k);
            if (pk < prf.dim()){ 
                for(size_t l=0; l < lmax+1; ++l){
                    (*pcart)(k) += (static_cast<float>((Y.SH(l,0)(pk,x0,y0)).real())) * legendre(l,0)(k);
                }
                for (size_t m(1); m < mmax+1; ++m){
                    mphi = m*pc_data(k);
                    calcos = cos(mphi);
                    calsin = sin(mphi);
                    for (size_t l(m); l < lmax+1; ++l){
                        YSH_re = static_cast<float>((Y.SH(l,m)(pk,x0,y0)).real());
                        YSH_im = static_cast<float>((Y.SH(l,m)(pk,x0,y0)).imag());
                        YSH_re *= calcos; 
                        YSH_im *= calsin; 
                        YSH_re -= YSH_im;
                        (*pcart)(k) += 2.0*(legendre(l,m)(k))*YSH_re;
                    }
                 }
              }
          }*/
           /*     for(size_t l=0; l < lmax+1; ++l){
                    YSH_re = static_cast<float>((Y.SH(l,0)(pind(k),x0,y0)).real());
                    (*pcart)(k) += YSH_re * legendre(l,0)(k);
                    for (size_t m=1; m < ((mmax < l)? mmax:l)+1; ++m){            
                        mphi = m*pc_data(k);
                        YSH_re = static_cast<float>((Y.SH(l,m)(pind(k),x0,y0)).real());
                        YSH_im = static_cast<float>((Y.SH(l,m)(pind(k),x0,y0)).imag());
                        YSH_re *= cosf(mphi); YSH_im *= sinf(mphi); YSH_re -= YSH_im;
                        (*pcart)(k) += (legendre(l,m)(k))*YSH_re;
                    }
                }  
            }   
        }    */

        return *p1x1cart;
   } 
//--------------------------------------------------------------

//**************************************************************


