///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Modified:	June 10th 2011
//	Last Modified:	July 26th 2011
///////////////////////////////////////////////////////////

//   
//   This cpp file contains the definition for the 
//   Runge-Kutta methods RK2, RK3, RK4.
//////////////////////////////////////

//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <math.h>
    #include <float.h>

//  My libraries
    #include "matrices.h"

//  Declerations
    #include "decl-input.h"
    #include "decl-state.h"
    #include "decl-vlasovmax.h"
    #include "decl-actions.h"
    #include "decl-RK.h"


//**************************************************************
//--------------------------------------------------------------
//  Definition of the pure virtual destructor for the abstract
//  class Explicit_Method
//--------------------------------------------------------------
    Explicit_Method::~Explicit_Method(){}
//**************************************************************


//**************************************************************
//**************************************************************
//   Definition for the 4th order Runge Kutta Method
//**************************************************************
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
    Runge_Kutta_4::Runge_Kutta_4(Stat& Yin, int tout_start)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       :t(static_cast<double>(tout_start) *
          Inputdata::IN().cont().dt_out),          // Initialize time = t0 *dtout
        Tout(t + Inputdata::IN().cont().dt_out),  // The next output time
        Y0(Yin), Y1(Yin), Yh(Yin),                // 3 local "State" Variables
        Y(Yin),                                   // Initialize a reference to "State" Y
        Act(Yh){                                  // Initialize "Actions"
    
        num_h  = size_t(static_cast<int>(Inputdata::IN().cont().dt_out/Inputdata::IN().inp().CLF))+1; 
        h      = Inputdata::IN().cont().dt_out/static_cast<double>(num_h);
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------
//  Output time
    double& Runge_Kutta_4::tout() {return Tout;}

//  Output time
    size_t  Runge_Kutta_4::numh() const {return num_h;}

//  Real time of the simulation (can only be modified in the RK)
    double Runge_Kutta_4::time() const {return t;}

//  Call Advection Actions 
    Stat& Runge_Kutta_4::F(Stat& Yslope) {

/*        switch (operation) {
            case 1: { // First half of the explicit method
                return Act.Adv1(Yslope); break;
            }
            case 2: { // Second half of the explicit method
                return Act.Adv2(Yslope); break;
            }
            case 3: { // Second half of the explicit method
                return Act.AdvExf1(Yslope); break;
            }
            case 4: { // Second half of the explicit method
                return Act.AdvEyf1(Yslope); break;
            }
            case 5: { // Second half of the explicit method
                return Act.AdvEzf1(Yslope); break;
            }
            default:
                cout<<"This is not valid operation. \n";
                exit(1);
                break;
       }    
       // You should throw an exception here;*/
       return Act.Adv(Yslope);
       
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Runge_Kutta_4& Runge_Kutta_4::advance_p1(){    
//--------------------------------------------------------------
//  Take a step using RK4
//--------------------------------------------------------------

//      Initialization
        Y0 = Y; Y1 = Y; 

//      Step 1
        Yh  = F(Y1);
        Yh *= (0.5*h);   Y1 += Yh;      // Y1 = Y1 + (h/2)*Yh
        Yh *= (1.0/3.0); Y  += Yh;      // Y  = Y  + (h/6)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        Yh  = F(Y1);     Y1  = Y0; 
        Yh *= (0.5*h);   Y1 += Yh;      // Y1 = Y0 + (h/2)*Yh
        Yh *= (2.0/3.0); Y  += Yh;      // Y  = Y  + (h/3)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 3
        Yh  = F(Y1); 
        Yh *= h;          Y0 += Yh;     // Y1 = Y0 + h*Yh
        Yh *= (1.0/3.0);  Y  += Yh;     // Y  = Y  + (h/3)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
//      Step 4
        Yh  = F(Y0);
        Yh *= (h/6.0);    Y += Yh;      // Y  = Y  + (h/6)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

//      Only update the time the first time you use RK
        t += h;

        return *this;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Runge_Kutta_4& Runge_Kutta_4::advance_p2()  {return *this;}    
    Runge_Kutta_4& Runge_Kutta_4::advance_IEx() {return *this;}    
    Runge_Kutta_4& Runge_Kutta_4::advance_IEy() {return *this;}    
    Runge_Kutta_4& Runge_Kutta_4::advance_IEz() {return *this;}    
//--------------------------------------------------------------

//**************************************************************


//**************************************************************
//**************************************************************
//   Definition for the 3rd order Runge Kutta Method
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Runge_Kutta_3::Runge_Kutta_3(Stat& Yin, int tout_start)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       :t(static_cast<double>(tout_start) *
          Inputdata::IN().cont().dt_out),          // Initialize time = t0 *dtout
        Tout(t + Inputdata::IN().cont().dt_out),  // The next output time
        Y0(Yin), Yh(Yin),                         // 3 local "State" Variables
        Y(Yin),                                   // Initialize a reference to "State" Y
        Act(Yh){                                  // Initialize "Actions"
        num_h  = size_t(static_cast<int>(Inputdata::IN().cont().dt_out/Inputdata::IN().inp().CLF))+1; 
        h      = Inputdata::IN().cont().dt_out/static_cast<double>(num_h);
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------
//  Output time
    double& Runge_Kutta_3::tout() {return Tout;}

//  Output time
    size_t  Runge_Kutta_3::numh() const {return num_h;}

//  Real time of the simulation (can only be modified in the RK)
    double Runge_Kutta_3::time() const {return t;}

//  Call Advection Actions 
    Stat& Runge_Kutta_3::F(Stat& Yslope) {

       /* switch (operation) {
            case 1: { // First half of the explicit method
                return Act.Adv1(Yslope); break;
            }
            case 2: { // Second half of the explicit method
                return Act.Adv2(Yslope); break;
            }
            case 3: { // Second half of the explicit method
                return Act.AdvExf1(Yslope); break;
            }
            case 4: { // Second half of the explicit method
                return Act.AdvEyf1(Yslope); break;
            }
            case 5: { // Second half of the explicit method
                return Act.AdvEzf1(Yslope); break;
            }
            default:
                cout<<"This is not valid operation. \n";
                exit(1);
                break;
        }*/
        // You should throw an exception here;
        return Act.Adv(Yslope);
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Runge_Kutta_3& Runge_Kutta_3::advance_p1(){    
//--------------------------------------------------------------
//  Take a step using RK3
//--------------------------------------------------------------

//      Do not proceed for explicit electric field and if the
//      operation # > 1
//	    if ( ( operation > 1 ) ) { 
//            return *this; 
//        } 

//      Initialization
        Y0 = Y; 

//      Step 1
        Yh  = F(Y0);
        Yh *= h;   Y0 += Yh;                                    
//      Y0 = Y0 + h*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        Yh  = F(Y0);  
        Yh *= h; Y0 += Yh; Y *= 3.0; Y0 += Y; Y0 *= 0.25;       
//      Y0 = 1/4 * ( 3*Y + (Y0 + h*Yh) )
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 3
        Yh  = F(Y0); 
        Yh *= h; Y0 += Yh; Y0 *= 2.0; Y *= (1.0/3.0); Y += Y0; Y *= (1.0/3.0); 
//      Y  = 1/3 * ( Y + 2 * (Y0+h*Yh) )
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

//      Only update the time the first time you use RK
//        if (operation == 1) {
            t += h;
//        }

        return *this;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Runge_Kutta_3& Runge_Kutta_3::advance_p2() {return *this;}    
    Runge_Kutta_3& Runge_Kutta_3::advance_IEx() {return *this;}    
    Runge_Kutta_3& Runge_Kutta_3::advance_IEy() {return *this;}    
    Runge_Kutta_3& Runge_Kutta_3::advance_IEz() {return *this;}    
//--------------------------------------------------------------

//**************************************************************


//**************************************************************
//**************************************************************
//   Definition for the 2nd order Runge Kutta Method
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Runge_Kutta_2::Runge_Kutta_2(Stat& Yin, int tout_start)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       :t(static_cast<double>(tout_start) *
          Inputdata::IN().cont().dt_out),          // Initialize time = t0 *dtout
        Tout(t + Inputdata::IN().cont().dt_out),   // The next output time
        Y0(Yin), Yh(Yin),                          // 3 local "State" Variables
        Y(Yin),                                    // Initialize a reference to "State" Y
        Act(Yh){                                   // Initialize "Actions"

        num_h  = size_t(static_cast<int>(Inputdata::IN().cont().dt_out/Inputdata::IN().inp().CLF))+1; 
        h      = Inputdata::IN().cont().dt_out/static_cast<double>(num_h);
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------
//  Output time
    double& Runge_Kutta_2::tout() {return Tout;}

//  Output time
    size_t  Runge_Kutta_2::numh() const {return num_h;}

//  Real time of the simulation (can only be modified in the RK)
    double Runge_Kutta_2::time() const {return t;}

//  Call Advection Actions 
    Stat& Runge_Kutta_2::F(Stat& Yslope) {

        /*switch (operation) {
            case 1: { // First half of the explicit method
                return Act.Adv1(Yslope); break;
            }
            case 2: { // Second half of the explicit method
                return Act.Adv2(Yslope); break;
            }
            case 3: { // Second half of the explicit method
                return Act.AdvExf1(Yslope); break;
            }
            case 4: { // Second half of the explicit method
                return Act.AdvEyf1(Yslope); break;
            }
            case 5: { // Second half of the explicit method
                return Act.AdvEzf1(Yslope); break;
            }
            default:
                cout<<"This is not valid operation. \n";
                exit(1);
                break;
        }*/
        // You should throw an exception here;
        return Act.Adv(Yslope);
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Runge_Kutta_2& Runge_Kutta_2::advance_p1(){    
//--------------------------------------------------------------
//  Take a step using RK2
//--------------------------------------------------------------
         
//      Do not proceed for explicit electric field and if the
//      operation # > 1
	//if ( ( operation > 1 ) ) { 
        //   return *this; 
        //} 

//      Initialization
        Y0 = Y; 

//      Step 1
        Yh  = F(Y0);
        Yh *= h;   Y0 += Yh;      // Y0 = Y0 +   h   * Yh
        Yh *= 0.5; Y  += Yh;      // Y  = Y  + (h/2) * Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        Yh  = F(Y0);    
        Yh *= (0.5*h);   Y += Yh;  // Y = Y + (h/2)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Only update the time the first time you use RK
        //if (operation == 1) {
            t += h;
        //}

        return *this;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Runge_Kutta_2& Runge_Kutta_2::advance_p2() {return *this;}    
    Runge_Kutta_2& Runge_Kutta_2::advance_IEx() {return *this;}    
    Runge_Kutta_2& Runge_Kutta_2::advance_IEy() {return *this;}    
    Runge_Kutta_2& Runge_Kutta_2::advance_IEz() {return *this;}    
//--------------------------------------------------------------

//**************************************************************



//**************************************************************
//**************************************************************
//   Definition for the 2nd order Runge Kutta Method
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Runge_Kutta_2_Imp::Runge_Kutta_2_Imp(Stat& Yin, int tout_start)
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
       :t(static_cast<double>(tout_start) *
          Inputdata::IN().cont().dt_out),          // Initialize time = t0 *dtout
        Tout(t + Inputdata::IN().cont().dt_out),   // The next output time
        Y0(Yin), Yh(Yin),                          // 3 local "State" Variables
        Y(Yin),                                    // Initialize a reference to "State" Y
        Act(Yh){                                   // Initialize "Actions"

        num_h  = size_t(static_cast<int>(Inputdata::IN().cont().dt_out/Inputdata::IN().inp().CLF))+1; 
        h      = Inputdata::IN().cont().dt_out/static_cast<double>(num_h);
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
//  Access
//--------------------------------------------------------------
//  Output time
    double& Runge_Kutta_2_Imp::tout() {return Tout;}

//  Output time
    size_t  Runge_Kutta_2_Imp::numh() const {return num_h;}

//  Real time of the simulation (can only be modified in the RK)
    double Runge_Kutta_2_Imp::time() const {return t;}

//  Call Advection Actions 
    Stat& Runge_Kutta_2_Imp::F1(Stat& Yslope) {
        return Act.Adv1(Yslope);
    }

//  Call Advection Actions 
    Stat& Runge_Kutta_2_Imp::F2(Stat& Yslope) {
        return Act.Adv2(Yslope);
    }

//  Call Advection Actions 
    Stat& Runge_Kutta_2_Imp::ImpEx(Stat& Yslope) {
        return Act.ImpEx(Yslope);
    }

//  Call Advection Actions 
    Stat& Runge_Kutta_2_Imp::ImpEy(Stat& Yslope) {
        return Act.ImpEy(Yslope);
    }

//  Call Advection Actions 
    Stat& Runge_Kutta_2_Imp::ImpEz(Stat& Yslope) {
        return Act.ImpEz(Yslope);
    }
//--------------------------------------------------------------
//--------------------------------------------------------------


//--------------------------------------------------------------
    Runge_Kutta_2_Imp& Runge_Kutta_2_Imp::advance_p1(){    
//--------------------------------------------------------------
//  Take a step using RK2 for the implicit part of the code
//--------------------------------------------------------------
         
//      Do not proceed for explicit electric field and if the
//      operation # > 1
//	if ( ( operation > 1 ) ) { 
//           return *this; 
//        } 

//      Initialization
        Y0 = Y; 

//      Step 1
        Yh  = F1(Y0);
        Yh *= h;   Y0 += Yh;      // Y0 = Y0 +   h   * Yh
        Yh *= 0.5; Y  += Yh;      // Y  = Y  + (h/2) * Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        Yh  = F1(Y0);    
        Yh *= (0.5*h);   Y += Yh;  // Y = Y + (h/2)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Only update the time the first time you use RK
//        if (operation == 1) {
            t += h;
//        }

        return *this;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Runge_Kutta_2_Imp& Runge_Kutta_2_Imp::advance_p2() {    
//--------------------------------------------------------------

//--------------------------------------------------------------

//      Initialization
        Y0 = Y; 

//      Step 1
        Yh  = F2(Y0);
        Yh *= h;   Y0 += Yh;      // Y0 = Y0 +   h   * Yh
        Yh *= 0.5; Y  += Yh;      // Y  = Y  + (h/2) * Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        Yh  = F2(Y0);    
        Yh *= (0.5*h);   Y += Yh;  // Y = Y + (h/2)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return *this;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    Runge_Kutta_2_Imp& Runge_Kutta_2_Imp::advance_IEx() {
//--------------------------------------------------------------

//--------------------------------------------------------------

//      Initialization
        Y0 = Y; 

//      Step 1
        Yh  = ImpEx(Y0);
        Yh *= h;   Y0 += Yh;      // Y0 = Y0 +   h   * Yh
        Yh *= 0.5; Y  += Yh;      // Y  = Y  + (h/2) * Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        Yh  = ImpEx(Y0);    
        Yh *= (0.5*h);   Y += Yh;  // Y = Y + (h/2)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return *this;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Runge_Kutta_2_Imp& Runge_Kutta_2_Imp::advance_IEy() {
//--------------------------------------------------------------

//--------------------------------------------------------------

//      Initialization
        Y0 = Y; 

//      Step 1
        Yh  = ImpEy(Y0);
        Yh *= h;   Y0 += Yh;      // Y0 = Y0 +   h   * Yh
        Yh *= 0.5; Y  += Yh;      // Y  = Y  + (h/2) * Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        Yh  = ImpEy(Y0);    
        Yh *= (0.5*h);   Y += Yh;  // Y = Y + (h/2)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return *this;
    }    
//--------------------------------------------------------------


//--------------------------------------------------------------
    Runge_Kutta_2_Imp& Runge_Kutta_2_Imp::advance_IEz() {
//--------------------------------------------------------------

//--------------------------------------------------------------

//      Initialization
        Y0 = Y; 

//      Step 1
        Yh  = ImpEz(Y0);
        Yh *= h;   Y0 += Yh;      // Y0 = Y0 +   h   * Yh
        Yh *= 0.5; Y  += Yh;      // Y  = Y  + (h/2) * Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//      Step 2
        Yh  = ImpEz(Y0);    
        Yh *= (0.5*h);   Y += Yh;  // Y = Y + (h/2)*Yh
//      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        return *this;
    }    
//--------------------------------------------------------------
//**************************************************************
