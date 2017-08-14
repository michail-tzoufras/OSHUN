///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Modified:	June 10th 2011
//	Last Modified:	June 26th 2011
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the 
//   actions to be applied to the right hand side of the
//   Runge-Kutta scheme
///////////////////////////////////////////////////////////
//
//   class Actions::
//
//      The Actions class essentially defines the modules
//      that will be used in the calculation of the right
//      hand side of the Runge-Kutta.
//      1. Electric Field
//      2. Current
//      3. Advection
//      4. Magnetic Field
//      5. Faraday's law
//      6. Ampere's law
//      7. A low-pass filter, which is currently deactivated
// 
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

    #ifndef DECLERATION_ACTIONS_H
    #define DECLERATION_ACTIONS_H

//**************************************************************
//**************************************************************
//   Definition for the Output namespace
//**************************************************************
//**************************************************************



//**************************************************************
        void CLEAN_LOW_PC(Stat& Yh);
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class Actions {
//--------------------------------------------------------------
//      Decleration of the Actions
//--------------------------------------------------------------
        private:
//          Define the Actions
            Stat& Yh;
            ElectricField Ef;
            MagneticField Bf;
            Current Jf;
            Spatial_Advection Af;
            MaxwellEq ME;

            size_t l0, m0;

        public:
//          Constructors/Destructors
            Actions(Stat& Yslope); 
         
//          Explicit Advance
            Stat& Adv(Stat& Yin);    

        };
//--------------------------------------------------------------

//--------------------------------------------------------------
        class Implicit_Actions {
//--------------------------------------------------------------
//      Decleration of the Actions
//--------------------------------------------------------------
        private:
//          Define the Actions
            Stat& Yh;
            ElectricField Ef;
            MagneticField Bf;
            Current Jf;
            Spatial_Advection Af;
            MaxwellEq ME;

            bool implicit_E;
            bool if_implicitES;
            bool if_implicit1D;
            size_t l0, m0;

        public:
//          Constructors/Destructors
            Implicit_Actions(Stat& Yslope); 
         
//          Explicit Advance
            Stat& Adv1(Stat& Yin);    
            Stat& Adv2(Stat& Yin);    

//          Implicit Advance
            Stat& ImpEx(Stat& Yin);    
            Stat& ImpEy(Stat& Yin);    
            Stat& ImpEz(Stat& Yin);    
        };
//--------------------------------------------------------------
//**************************************************************

    #endif
