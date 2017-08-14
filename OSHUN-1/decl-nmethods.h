///////////////////////////////////////////////////////////
//   Contributing authors :  Michail Tzoufras
//   Modified:               June 14th 2010
//   Last Modified:          July 28th 2011
///////////////////////////////////////////////////////////
//   
//   Iteration loops to invert matrices provided by other
//   parts of the code. Note that the matrices.h header
//   needs to be included ahead of this header so that
//   the matrices in this file are well defined.
//
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H


//*******************************************************************
//-------------------------------------------------------------------
     bool Gauss_Seidel(Matrix2D<double>& A, 
                       valarray< complex<double> >& b,
                       valarray< complex<double> >& xk);
//*******************************************************************

//*******************************************************************
//-------------------------------------------------------------------
     void TridiagonalSolve (const valarray<double>& a, 
                            const valarray<double>& b, 
                                  valarray<double>& c,      
                                  valarray< complex<double> >  d,
                                  valarray< complex<double> >& x);
//-------------------------------------------------------------------


//-------------------------------------------------------------------
     bool Thomas_Tridiagonal(Matrix2D<double>& A, 
                       valarray< complex<double> >& d,
                       valarray< complex<double> >& xk); 
//-------------------------------------------------------------------
//*******************************************************************


//*******************************************************************
//-------------------------------------------------------------------
    complex<double> Det33(/*const valarray< complex<double> >& D, */
                          Matrix2D< complex<double> >& A);
//-------------------------------------------------------------------
    complex<double> Detx33(valarray< complex<double> >& D, 
                           Matrix2D< complex<double> >& A);
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    complex<double> Dety33(valarray< complex<double> >& D, 
                           Matrix2D< complex<double> >& A);
//-------------------------------------------------------------------

//-------------------------------------------------------------------
    complex<double> Detz33(valarray< complex<double> >& D,
                           Matrix2D< complex<double> >& A);
 //-------------------------------------------------------------------
//*******************************************************************


#endif
