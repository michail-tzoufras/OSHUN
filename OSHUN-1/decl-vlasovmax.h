///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//   Modified:  	Sep  12th 2009
//   Last Modified:	July 26th 2011
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the  
//   classes ElectricField, MagneticField, Current, Spatial_Advection 
//   and MaxwellEq.
///////////////////////////////////////////////////////////
//
//   classes :
//
//   1. class ElectricField :
//      The class ElectricField  stores a reference to the slope Yh 
//      and steps it calls Exyz(Yin), where Yin is the current
//      state, in order to calculate the effect of the 
//      electric field on the distribution function.
//
//   2. class MagneticField :
//      The class Magneticfield stores a reference to the slope Yh 
//      and steps it calls Bxyz(Yin), where Yin is the current
//      state, in order to calculate the effect of the 
//      magnetic field on the distribution function.
//
//   3. class Current :
//      The class Current stores a reference to the slope Yh and
//      in subsequent steps it calls Jx(Yin), where Yin is the 
//      state, in order  to calculate effect of the current
//      on the Electric field.
//
//   4. Spatial_Advection :
//      The class Spatial_Advection  stores a reference to the slope Yh and
//      in subsequent steps it calls Axy(Yin), where Yin is the 
//      state, in order to calculate effect of the spatial 
//      advection on the distribution function. 
//
//   5. class MaxwellEq :
//      The class MaxwellEq stores a reference to the slope Yh and
//      in subsequent steps it calls Ampere(Yin) and Faraday(Yin),
//      where Yin is the state, in order to update the electromagnetic
//      fields 
// 
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

    #ifndef DECLERATION_EFIELD_H
    #define DECLERATION_EFIELD_H

//**************************************************************
//**************************************************************
//   Decleration of the Electric Field 
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class ElectricField {
//--------------------------------------------------------------
//      Decleration of the Electric Field
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            ElectricField(Stat& Yslope); 

//          Advance
            Stat& Exyz(Stat& Yin);  
            Stat& Exyz1D(Stat& Yin);  
            Stat& Implicit_Ex(Stat& Yin, Matrix2D< complex<double> >& Ex);  
            Stat& Implicit_Ex1D(Stat& Yin, Matrix2D< complex<double> >& Ex);  
            Stat& Implicit_Eyz(Stat& Yin, Matrix2D< complex<double> >& Em, Matrix2D< complex<double> >& Ep);  

        private:
//          Calculation of G, H
            SHarmonic& G00(SHarmonic& f);
            void MakeGH(size_t el, SHarmonic& f);

//          Variables
            Stat& Yh;
            SHarmonic H, G, TMP; 

//          Sizes
            size_t                      l0, m0;

//          Multiplication parameters
            Matrix2D< complex<double> > A1, A2;
            valarray< complex<double> > B1, B2;
            valarray< complex<double> > C1, C3;
            Matrix2D< complex<double> > C2, C4;
            valarray< complex<double> > Hp0;

//          Axis
            Axis< complex<double> >     pr, invpr;

        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Decleration of the Magnetic Field 
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class MagneticField {
//--------------------------------------------------------------
//      Decleration of the Electric Field
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            MagneticField(Stat& Yslope); 
//          Advance
            Stat& Bxyz(Stat& Yin);  

        private:
//          Variables
            Stat& Yh;
            SHarmonic FLM; 

//          Sizes
            size_t                      l0, m0;

//          Multiplication parameters
            Matrix2D< complex<double> > A2;
            valarray< complex<double> > A1, B1;
            complex<double>             A3;

        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Decleration of the Current
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class Current {
//--------------------------------------------------------------
//      Decleration of the Electric Field
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Current(Stat& Yslope); 
//          Advance
            Stat& Jxyz(Stat& Yin);  

        private:
            Stat& Yh;
            Field2D jayx, jayy, jayz;
            complex<double> Delta_p, small;
            Axis< complex<double> >  p3og;
        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Decleration of the Spatial Advection
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class Spatial_Advection {
//--------------------------------------------------------------
//      Decleration of the Electric Field
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Spatial_Advection(Stat& Yslope); 
//          Advance
            Stat& Axy(Stat& Yin);  

        private:
            Stat&                       Yh;
            size_t                      l0, m0;
            complex<double>             idx, idy;

            Matrix2D< complex<double> > A1, A2, C2, C4;
            valarray< complex<double> > B1, B2, C1, C3;

            SHarmonic                   fd1, fd2;
            Axis< complex<double> >     vr;
        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//   Decleration of Maxwell's Equations
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class MaxwellEq {
//--------------------------------------------------------------
//      Decleration of the Electric Field
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            MaxwellEq(Stat& Yslope); 
//          Faraday's Law
            Stat& Faraday(Stat& Yin);  
//          Ampere's Law (except for the current)
            Stat& Ampere(Stat& Yin);  

        private:
            Stat& Yh;
            Field2D tmpEB;

            complex<double>             idx, idy;
        };
//--------------------------------------------------------------
//**************************************************************


    #endif
