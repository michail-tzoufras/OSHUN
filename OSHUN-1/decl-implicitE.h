///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	First Created:	June  2nd  2011
//	Last Modified:	July 27th  2011
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the 
//   classes required to calculate the effect of the implicit
//   electric field and collisions for the first order harmonic
///////////////////////////////////////////////////////////
//
//   class Implicit_Efield::
//
// 
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

    #ifndef DECLERATION_IMPLICITEFIELD_H
    #define DECLERATION_IMPLICITEFIELD_H

//**************************************************************
//**************************************************************
//   Definition for the Electric_Field_Methods namespace
//**************************************************************
//**************************************************************


    namespace Electric_Field_Methods { 
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>



//**************************************************************
//--------------------------------------------------------------
        class Current_xyz {
//--------------------------------------------------------------
//      Decleration of the Current 
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Current_xyz(Stat& Yin); 

//          Choose component
            Field2D& J(int component);  
            void calculate_J(Stat& Yin);  
            void calculate_J(int component, Stat& Yin);  

//          Components
            Field2D& Jx();  
            void calculate_Jx(Stat& Yin);  
            Field2D& Jy();  
            void calculate_Jy(Stat& Yin);  
            Field2D& Jz();  
            void calculate_Jz(Stat& Yin);  

        private:
            Field2D jayx, jayy, jayz;

            complex<double> Delta_p, small;
            Axis< complex<double> >  p3og;

        };
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
        class Efield_xyz {
//--------------------------------------------------------------
//      Decleration of the Current 
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Efield_xyz(Stat& Yin); 

//          Choose component
            Field2D& E(int component);  

//          Components
            Field2D& Ex();  
            Field2D& Ey();  
            Field2D& Ez();  

        private:
            Field2D efieldx, efieldy, efieldz;

//            complex<double> Delta_p, small;
//            Axis< complex<double> >  p3og;

        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class Efield_Method {
//--------------------------------------------------------------
//      Abstract class for explicit methods
//--------------------------------------------------------------
        public:
            virtual void advance(Explicit_Method* rk, Euler_Backward* eb) = 0; // "covariant" return
            virtual bool implicitE()  const = 0;
            virtual ~Efield_Method() = 0;
        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Explicit_E_Field: public Efield_Method {
//--------------------------------------------------------------
//      class for explicit electric field              
//--------------------------------------------------------------
        public:
//          Constructor
            Explicit_E_Field(Stat& Yin, int tout_start);

//          Main function
            void advance(Explicit_Method* rk, Euler_Backward* eb);    
         
//          Is the electric field implicit?
            bool implicitE() const;

        private:

//          Reference to Y
            Stat& Y;

//          Time
            double t;

//          This class does not do anything at this stage            

        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Implicit_E_Field: public Efield_Method {
//--------------------------------------------------------------
//      class for ==> implicit electric field              
//--------------------------------------------------------------
        public:
//          Constructor
            Implicit_E_Field(Stat& Yin, int tout_start);

//          Main function
            void advance(Explicit_Method* rk, Euler_Backward* eb);    
         
//          Is the electric field implicit?
            bool implicitE() const;

        private:
//          Reference to Y
            Stat& Y;

//          Storage of harmonics
            SHarmonic f00; 
            SHarmonic f10, f11;
            SHarmonic f20, f21, f22;

//          Time
            double t;

//          Current and electric field 
            Current_xyz JN, J0, J_Ex, J_Ey, J_Ez;
            Efield_xyz  EN, E0, DE;

//          Ampere's law JN = rot(B)
            void Ampere();
            void FindDE();

//          Boundary Cells
            int               Nbc, szx, szy;
//          Derivative constants
            complex<double>   idx, idy;
            size_t            l0, m0;

//          Implicit methods
            bool if_implicitES;
            bool if_implicit1D;
        };
//--------------------------------------------------------------
//**************************************************************

    }
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    #endif
