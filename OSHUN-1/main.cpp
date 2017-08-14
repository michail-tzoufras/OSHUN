///////////////////////////////////////////////////////////
//   Contributing authors :  Michail Tzoufras
///////////////////////////////////////////////////////////
//   
//   Main program, including the main loop. 
//
///////////////////////////////////////////////////////////

// Standard libraries 
#include <iostream>
#include <vector>
#include <valarray>
#include <complex>

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

// My libraries 
#include "matrices.h"
//#include "it_methods.h"

// Interface 
#include "decl-input.h"
#include "decl-state.h"
#include "decl-setup.h"
#include "decl-export.h"
#include "decl-parallel.h"
#include "decl-vlasovmax.h"
#include "decl-plsource.h"
#include "decl-actions.h"
#include "decl-RK.h"
#include "decl-fokkerplanck.h"
#include "decl-implicitE.h"


//*******************************************************************
//------------------------------------------------------------------- 
    void main_loop(Stat& Y, Parallel_Environment& CX, Explicit_Method* rk, 
                            Euler_Backward* eb, Explicit_EC_FP* ee, 
                            Electric_Field_Methods::Efield_Method* em ){
//--------------------------------------------------------------------------- 
//  Loop for the large output steps and smaller timesteps
//--------------------------------------------------------------------------- 


        // Reading restart files
        if ( CX.READ_RESTART() ) {
            if (CX.RANK() == 0) cout<<"Reading restart files at time = " 
                                    << CX.T_IN()*Inputdata::IN().cont().dt_out<<"\n";
            CX.Read_Restart(Y);
            CX.Output(CX.T_IN(),Y);
            if ( Inputdata::IN().list().only_output ) {  
                CX.Output(CX.T_IN(),Y);
                cout << "This is output \n";
                MPI_Finalize();
                exit(1);
            }
        }
        else {
            // Initial output
            CX.Output(0,Y);
        }

        // Initialize laser heat source
        PLaserSource PLS(Y,CX.T_IN());
        ExternalEfield EFLD(Y,CX.T_IN());
        InverseBremsstrahlung IBs(Y,CX.T_IN());
//        NumericsInfo NMI(CX.T_IN());

        // Loop for output steps 
        for (size_t t_step = CX.T_IN()+1; t_step < Inputdata::IN().cont().n_out+1; ++t_step){
            (*rk).tout() = Inputdata::IN().cont().dt_out * t_step;   // next output step

            // Loop for small timesteps
            for (size_t h_step(0); h_step < (*rk).numh(); ++h_step){ 

                (*rk).advance_p1();                                        // explicit step part 1
                                                                           // Energy history diagnostic
                                                                           //   NMI.Collect_Numerics(Y,(*rk).time(),h_step);  

                if ( Inputdata::IN().list().implicit_E ) {
                    CX.Neighbor_ImplicitE_Communications(Y);               // these are additional communications for the implicit E-field
                }

                (*ee).loop((*rk).time());                                  // explicit e-e collision loop
                                                                           //    NMI.Collect_Numerics(Y,(*rk).time(),h_step);        

                if ( PLS.PLSOURCE() )    PLS.Laser_Source((*rk).time());   // phenomenological laser source
//                if (EFLD.POL() != 0) EFLD.Laser_Efield((*rk).time());    // laser electric field 
                if (IBs.IBSOURCE() == 1) IBs.loop((*rk).time());       // laser electric field 
                                                                           //    NMI.Collect_Numerics(Y,(*rk).time(),h_step);               

                // THIS IS WHERE YOU SHOULD DO IMPLICIT E
                // Collisions for f10,f11
                // Calculate E
                using Electric_Field_Methods::Efield_Method;
                (*em).advance( rk, eb);
                (*rk).advance_p2();                                        // explicit step part 2
                                                                           //    NMI.Collect_Numerics(Y,(*rk).time(),h_step);         

                (*eb).advance_1((*rk).time());                             // implicit collisions for f10, f11
                (*eb).advance_lm((*rk).time());                            // implicit collisions for flm
                                                                           //    NMI.Collect_Numerics(Y,(*rk).time(),h_step);              

                CX.Neighbor_Communications(Y);                             // update guard cells
                                                                           //    NMI.Collect_Numerics(Y,(*rk).time(),h_step);               

                                                                           // Export the buffer with numerics info
                                                                           // if ( NMI.Export_Numerics() ) {
                                                                           //    Export_Facility::Export_numerics(CX.RANK(), NMI.List_Size(), NMI.Energy_Info);
                                                                           // }

                if (CX.RANK() == 0) cout<<"time = " << (*rk).time() <<"\n";
            }

                                                                           // Empty the buffer with numerics information
                                                                           // if ( NMI.Energy_Info.size() > 0 ) {
                                                                           //    Export_Facility::Export_numerics(CX.RANK(), NMI.List_Size(), NMI.Energy_Info);
                                                                           // }

            // Generate Output 
            if (CX.RANK() == 0) cout<<"Generating output...\n";
            CX.Output(t_step,Y);
            // Also export the "numerics buffer here in case there is still info in it.

            // Write restart files
            if ( CX.WRITE_RESTART(t_step) ) {
                if (CX.RANK() == 0) cout<<"Writing restart files...\n";
                CX.Write_Restart(t_step,Y);
            }
        }
    }
//--------------------------------------------------------------------------- 

//--------------------------------------------------------------------------- 
    int main(int argc, char** argv){
//--------------------------------------------------------------------------- 

        MPI_Init(&argc,&argv);

//      Initiate the Paralell Environment / Decompose the Computational Domain  
        Parallel_Environment CX; 

//      Define the State
        Stat Y;

//      Initialize 
        if (CX.RANK() == 0) Setup_Parameters::folders();
        Setup_Y::initialize(Y);

        switch (Inputdata::IN().list().RKLevel) {
            case 2: { // 2nd order Runge-Kutta
                if (CX.RANK() == 0) cout<<"Initializing RK"
                                        << Inputdata::IN().list().RKLevel<<"...\n\n";
                Explicit_EC_FP ee(Y,CX.T_IN());
                Euler_Backward eb(Y,CX.T_IN());
                if ( Inputdata::IN().list().implicit_E ) { 
                    Runge_Kutta_2_Imp rk2(Y,CX.T_IN());
                    Electric_Field_Methods::Implicit_E_Field eim(Y,CX.T_IN());
                    main_loop(Y, CX, &rk2, &eb, &ee, &eim); 
                } 
                else { 
                    Runge_Kutta_2 rk2(Y,CX.T_IN());
                    Electric_Field_Methods::Explicit_E_Field eem(Y,CX.T_IN());
                    main_loop(Y, CX, &rk2, &eb, &ee, &eem);
                }
                break;
            }
            case 3: { // 3nd order Runge-Kutta
                if (CX.RANK() == 0) cout<<"Initializing RK"
                                        << Inputdata::IN().list().RKLevel<<"...\n\n";
                Explicit_EC_FP ee(Y,CX.T_IN());
                Euler_Backward eb(Y,CX.T_IN());
                if ( Inputdata::IN().list().implicit_E ) { 
                    Runge_Kutta_2_Imp rk2(Y,CX.T_IN());
                    Electric_Field_Methods::Implicit_E_Field eim(Y,CX.T_IN());
                    main_loop(Y, CX, &rk2, &eb, &ee, &eim); 
                } 
                else { 
                    Runge_Kutta_3 rk3(Y,CX.T_IN());
                    Electric_Field_Methods::Explicit_E_Field eem(Y,CX.T_IN());
                    main_loop(Y, CX, &rk3, &eb, &ee, &eem);
                }
                break;
            }
            case 4: { // 4th order Runge-Kutta
                if (CX.RANK() == 0) cout<<"Initializing RK"
                                        << Inputdata::IN().list().RKLevel<<"...\n\n";
                Explicit_EC_FP ee(Y,CX.T_IN());
                Euler_Backward eb(Y,CX.T_IN());
                if ( Inputdata::IN().list().implicit_E ) { 
                    Runge_Kutta_2_Imp rk2(Y,CX.T_IN());
                    Electric_Field_Methods::Implicit_E_Field eim(Y,CX.T_IN());
                    main_loop(Y, CX, &rk2, &eb, &ee, &eim); 
                } 
                else { 
                    Runge_Kutta_4 rk4(Y,CX.T_IN());
                    Electric_Field_Methods::Explicit_E_Field eem(Y,CX.T_IN());
                    main_loop(Y, CX, &rk4, &eb, &ee, &eem);
                }
                break;
            }  
            default:
                if (CX.RANK() == 0) cout<<"This RK level is not valid. \n";
                break;
        }
        if (CX.RANK() == 0) std::cout<<"Terminating... \n";

        MPI_Finalize();

    }
//--------------------------------------------------------------------------- 


