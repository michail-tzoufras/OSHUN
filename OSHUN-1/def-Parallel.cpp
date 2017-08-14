///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	May 18 2011
///////////////////////////////////////////////////////////

//   
//   Contains the declerations for the communications
//   between nodes, boundaries and the parallel output
///////////////////////////////////////////////////////////
//
//   This file contains three modules:
//
//   1. class Node_Communications: 
//        Allows the nodes to exchange information in order
//        to update their guard cells. For boundary nodes 
//        it provides the appropriate boundary conditions.
//
//   2. class Parallel_Output: 
//        Collects the output information from all of the
//        nodes and combines them to the final file ready
//        to be exported. For the moments of the distribution
//        function it also performs the integration over 
//        momentum space and for detailed information on 
//        phasespace it calls the appropriate functions from
//        "Output" to convert the spherical harmonics to 
//        cartesian geometry.
// 
//   3. class Parallel_Environment:
//        - It decomposes the computational domain
//        - It controls the node communications
//        - It controls the parallel output
//        - It controls the restart facility
//
///////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

//  Standard libraries
    #include <iostream>
    #include <algorithm>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <cmath>
    #include <stdio.h>
    #include <mpi.h>
    #include <float.h>
    #include <string>

//  My libraries
    #include "matrices.h"

//  Declerations
    #include "decl-input.h"
    #include "decl-state.h"
    #include "decl-output.h"
    #include "decl-export.h"
    #include "decl-parallel.h"



//**************************************************************
//**************************************************************
//  Definition of the Nodes Communications class
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Node_ImplicitE_Communications:: Node_ImplicitE_Communications() : 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        Nbc(Inputdata::IN().list().RKLevel),      // # of boundary cells
        bndX(Inputdata::IN().list().bndX),        // Type of boundary in X
        bndY(Inputdata::IN().list().bndY){        // Type of boundary in Y
         
        using Inputdata::IN; 

        // # of harmonics
        // msg_sizeX = ((IN().inp().m0+1)*(2*IN().inp().l0-IN().inp().m0+2))/2; 
        // (# of harmonics) * (# cells in p)
        // msg_sizeX *= IN().inp().pr.dim();
        // 6 fields: Ex, Ey, Ez, Bx, By, Bz
        //msg_sizeX += 6;  
     
        // 3 components for Bx, By, Bz
        msg_sizeX = 3;  
        msg_sizeY = msg_sizeX; 
        
        msg_sizeX *= (IN().inp().y.dim()*Nbc);   
        msg_sizeY *= (IN().inp().x.dim()*Nbc);  
        msg_bufX = new complex<double>[msg_sizeX];
        msg_bufY = new complex<double>[msg_sizeY];
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Node_ImplicitE_Communications:: ~Node_ImplicitE_Communications(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
        delete[] msg_bufX;
        delete[] msg_bufY;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    int Node_ImplicitE_Communications:: BNDX()  const {return bndX;} 
    int Node_ImplicitE_Communications:: BNDY()  const {return bndY;} 
//--------------------------------------------------------------

//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Send and receive in the X direction
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::Send_right_X(Stat& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the right boundary and send them 
//           to the node on the right
//--------------------------------------------------------------
        // static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().y.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().y.dim()*Nbc);
        size_t bufind(0);

        // Harmonics:x0 "Right-Bound ---> " 
        // for(size_t i = 0; i < Y.DF().dim(); ++i){
        //     copy(Y.SH(i).x0_RB(Nbc), Y.SH(i).x0_RB(Nbc).end(), msg_bufX + bufind); 
        //    bufind += step_h;
        // }
        // Fields:   x0 "Right-Bound --> "
        for(size_t i(3); i < Y.EMF().dim(); ++i){  // "3" as opposed to "0"
            copy(Y.FLD(i).x0_RB(Nbc), Y.FLD(i).x0_RB(Nbc).end(), msg_bufX + bufind);
            bufind += step_f;
        } 

        MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::Recv_from_left_X(Stat& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the left and update
//           the left guard cells
//--------------------------------------------------------------
        // static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().y.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().y.dim()*Nbc);
        size_t bufind(0);
        MPI_Status status; 

        // Receive Data
        MPI_Recv(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &status);

        // Harmonics:x0-"---> Left-Guard" 
        //for(int i = 0; i < Y.DF().dim(); ++i){
        //    copy(msg_bufX + bufind, msg_bufX + bufind + step_h, Y.SH(i).x0_LG(Nbc)); 
        //    bufind += step_h;
        //}
        // Fields:   x0-"---> Left-Guard"
        for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            copy(msg_bufX + bufind, msg_bufX + bufind + step_f, Y.FLD(i).x0_LG(Nbc)); 
            bufind += step_f;
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::Send_left_X(Stat& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the left boundary and send them 
//           to the node on the left 
//--------------------------------------------------------------
        // static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().y.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().y.dim()*Nbc);
        size_t bufind(0); 

        // Harmonics:x0 " <--- Left-Bound "
        //for(int i = 0; i < Y.DF().dim(); ++i){
        //    copy(Y.SH(i).x0_LB(Nbc), Y.SH(i).x0_LB(Nbc).end(), msg_bufX + bufind);
        //    bufind += step_h;
        //} 
        // Fields:   x0 " <--- Left-Bound "
        for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            copy(Y.FLD(i).x0_LB(Nbc), Y.FLD(i).x0_LB(Nbc).end(), msg_bufX + bufind); 
            bufind += step_f;
        }

        MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 1, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::Recv_from_right_X(Stat& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the right and update
//           the right guard cells
//--------------------------------------------------------------
        // static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().y.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().y.dim()*Nbc);
        size_t bufind(0);
        MPI_Status status; 

        // Receive Data
        MPI_Recv(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, origin, 1, MPI_COMM_WORLD, &status);

        // Harmonics:x0-"Right-Guard <--- " 
        // for(int i = 0; i < Y.DF().dim(); ++i){
        //    copy(msg_bufX + bufind, msg_bufX + bufind + step_h, Y.SH(i).x0_RG(Nbc)); 
        //    bufind += step_h;
        // }
        // Fields:   x0-"Right-Guard <--- "
        for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            copy(msg_bufX + bufind, msg_bufX + bufind + step_f, Y.FLD(i).x0_RG(Nbc)); 
            bufind += step_f;
        }
    }
//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Send and receive in the Y direction
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::Send_right_Y(Stat& Y, int dest) {
//--------------------------------------------------------------
//  Y-axis : Read data from the right boundary and send them 
//           to the node on the right
//--------------------------------------------------------------
        // static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().x.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().x.dim()*Nbc);
        size_t bufind(0); 

        // Harmonics:x0 "Right-Bound ---> " 
        // for(int i = 0; i < Y.DF().dim(); ++i){
        //     copy(Y.SH(i).y0_RB(Nbc), Y.SH(i).y0_RB(Nbc).end(), msg_bufY + bufind); 
        //     bufind += step_h;
        // }
        // Fields:   x0 "Right-Bound --> "
        for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            copy(Y.FLD(i).y0_RB(Nbc), Y.FLD(i).y0_RB(Nbc).end(), msg_bufY + bufind); 
            bufind += step_f;
        }

        MPI_Send(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::Recv_from_left_Y(Stat& Y, int origin) {
//--------------------------------------------------------------
//  Y-axis : Receive data from the node on the left and update
//           the left guard cells
//--------------------------------------------------------------
        // static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().x.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().x.dim()*Nbc);
        size_t bufind(0);
        MPI_Status status;
        // Receive Data
        MPI_Recv(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &status);

        // Harmonics:x0-"---> Left-Guard" 
        // for(int i = 0; i < Y.DF().dim(); ++i){
        //    copy(msg_bufY + bufind, msg_bufY + bufind + step_h, Y.SH(i).y0_LG(Nbc)); 
        //    bufind += step_h;
        // }
        // Fields:   x0-"---> Left-Guard"
        for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
            copy(msg_bufY + bufind, msg_bufY + bufind + step_f, Y.FLD(i).y0_LG(Nbc)); 
            bufind += step_f;
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::Send_left_Y(Stat& Y, int dest) {
//--------------------------------------------------------------
//  Y-axis : Read data from the left boundary and send them 
//           to the node on the left 
//--------------------------------------------------------------
        // static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().x.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().x.dim()*Nbc);
        size_t bufind(0); 

        // Harmonics:x0 " <--- Left-Bound "
        // for(int i = 0; i < Y.DF().dim(); ++i) {
        //     copy(Y.SH(i).y0_LB(Nbc),Y.SH(i).y0_LB(Nbc).end(), msg_bufY + bufind); 
        //     bufind += step_h;
        // }
        // Fields:   x0 " <--- Left-Bound "
        for(int i(3); i < Y.EMF().dim(); ++i) { // "3" as opposed to "0"
            copy(Y.FLD(i).y0_LB(Nbc), Y.FLD(i).y0_LB(Nbc).end(), msg_bufY + bufind); 
            bufind += step_f;
        }

        MPI_Send(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, dest, 1, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::Recv_from_right_Y(Stat& Y, int origin) {
//--------------------------------------------------------------
//  Y-axis : Receive data from the node on the right and update
//           the right guard cells
//--------------------------------------------------------------
        // static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().x.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().x.dim()*Nbc);
        size_t bufind(0); 
        MPI_Status status; 
        // Receive Data
        MPI_Recv(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, origin, 1, MPI_COMM_WORLD, &status);

        // Harmonics:x0-"Right-Guard <--- " 
        // for(int i = 0; i < Y.DF().dim(); ++i) {
        //     copy(msg_bufY + bufind, msg_bufY + bufind + step_h, Y.SH(i).y0_RG(Nbc)); 
        //     bufind += step_h;
        // }
        // Fields:   x0-"Right-Guard <--- "
        for(int i(3); i < Y.EMF().dim(); ++i) { // "3" as opposed to "0"
            copy(msg_bufY + bufind, msg_bufY + bufind + step_f, Y.FLD(i).y0_RG(Nbc)); 
            bufind += step_f;
        }
    }
//--------------------------------------------------------------


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  Boundary conditions
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Mirror conditions on one side
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//**************************************************************
//--------------------------------------------------------------
    void Node_ImplicitE_Communications::mirror_bound_Xleft(Stat& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction on the left
//--------------------------------------------------------------
        // int sign(1);
        size_t Nx(Y.SH(0,0).numx());

        // Mirror the harmonics
        // for(int l(0); l < Y.DF().l0(); ++l){
        //     for(int m(0); m < ((Y.DF().m0() < l)? Y.DF().m0():l)+1; ++m){
        //         sign = 1-2*((l+m)%2);          //(-1)^(m+n)

        //         for (int c(0); c < Nbc; ++c)
        //             copy(Y.SH(l,m).x0(2*Nbc-c-1), Y.SH(l,m).x0(2*Nbc-c-1).end(), Y.SH(l,m).x0(c));
        //         for(GSIter_CD it(Y.SH(l,m).x0_LG(Nbc)); it != it.end(); ++it)  *it *= sign; 
        //     }
        // }

        // Mirror the fields 
        for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
           for (int c(0); c < Nbc; ++c) 
              copy(Y.FLD(i).x0(2*Nbc-c-1), Y.FLD(i).x0(2*Nbc-c-1).end(), Y.FLD(i).x0(c));
        }

        //Ey
        // for(GSIter_CD it(Y.EMF().Ey().x0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary
        //Ez
        // for(GSIter_CD it(Y.EMF().Ez().x0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary
        //Bx
        for(GSIter_CD it(Y.EMF().Bx().x0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::mirror_bound_Xright(Stat& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction on the right
//--------------------------------------------------------------
        // int sign(1);
        size_t Nx(Y.SH(0,0).numx());

        // Mirror the harmonics
        // for(int l(0); l < Y.DF().l0(); ++l){
        //     for(int m(0); m < ((Y.DF().m0() < l)? Y.DF().m0():l)+1; ++m){
        //         sign = 1-2*((l+m)%2);          //(-1)^(m+n)

        //         for (int c(0); c < Nbc; ++c) 
        //             copy(Y.SH(l,m).x0(Nx-2*Nbc+c), Y.SH(l,m).x0(Nx-2*Nbc+c).end(), Y.SH(l,m).x0(Nx-c-1));
        //         for(GSIter_CD it(Y.SH(l,m).x0_RG(Nbc)); it != it.end(); ++it)  *it *= sign;
 
        //     }
        // }

        // Mirror the fields 
        for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
           for (int c(0); c < Nbc; ++c) 
               copy(Y.FLD(i).x0(Nx-2*Nbc+c), Y.FLD(i).x0(Nx-2*Nbc+c).end(),  Y.FLD(i).x0(Nx-c-1));
        }

        //Ey
        // for(GSIter_CD it(Y.EMF().Ey().x0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary 
        //Ez
        // for(GSIter_CD it(Y.EMF().Ez().x0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
        //Bx
        for(GSIter_CD it(Y.EMF().Bx().x0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::mirror_bound_Yleft(Stat& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the y direction on the left
//--------------------------------------------------------------
        
        // int sign(1);
        complex< double > cmpl(1,0);
        size_t Ny(Y.SH(0,0).numy());
        // Mirror the harmonics
        // for(int l(1); l < Y.DF().l0(); ++l){
        //     for(int m(0); m < ((Y.DF().m0() < l)? Y.DF().m0():l)+1; ++m){
        //         sign = 1-2*(m%2);          //(-1)^m

                // left boundary
        //         for (int c(0); c < Nbc; ++c) 
        //            copy(Y.SH(l,m).y0(2*Nbc-c-1), Y.SH(l,m).y0(2*Nbc-c-1).end(), Y.SH(l,m).y0(c));
        //        for(GSIter_CD it(Y.SH(l,m).y0_LG(Nbc)); it != it.end(); ++it) {
        //           cmpl  = *it; 
        //           *it  -= 2.0 * cmpl.real();
        //           *it *= ((-1)*sign); 
        //        }
        //    }
        // }

        // Mirror the fields 
        for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
           for (int c(0); c < Nbc; ++c) 
              copy(Y.FLD(i).y0(2*Nbc-c-1), Y.FLD(i).y0(2*Nbc-c-1).end(), Y.FLD(i).y0(c));
        }

        //Ex
        // for(SIter_CD it(Y.EMF().Ex().y0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary
        //Ez
        // for(SIter_CD it(Y.EMF().Ez().y0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary
        //By
        for(SIter_CD it(Y.EMF().By().y0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::mirror_bound_Yright(Stat& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the y direction on the right
//--------------------------------------------------------------
        // int sign(1);
        complex< double > cmpl(1,0);
        size_t Ny(Y.SH(0,0).numy());
        // Mirror the harmonics
        // for(int l(1); l < Y.DF().l0(); ++l){
        //     for(int m(0); m < ((Y.DF().m0() < l)? Y.DF().m0():l)+1; ++m){
        //         sign = 1-2*(m%2);          //(-1)^m

                // right boundary
        //         for (int c(0); c < Nbc; ++c) 
        //             copy(Y.SH(l,m).y0(Ny-2*Nbc+c), Y.SH(l,m).y0(Ny-2*Nbc+c).end(), Y.SH(l,m).y0(Ny-c-1));
        //         for(GSIter_CD it(Y.SH(l,m).y0_RG(Nbc)); it != it.end(); ++it){
        //             cmpl  = *it; 
        //             *it  -= 2.0 * cmpl.real(); 
        //             *it *= ((-1)*sign);
        //         }
        //     }
        // }

        // Mirror the fields 
        for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
           for (int c(0); c < Nbc; ++c) 
               copy(Y.FLD(i).y0(Ny-2*Nbc+c), Y.FLD(i).y0(Ny-2*Nbc+c).end(), Y.FLD(i).y0(Ny-c-1));
        }

        //Ex
        // for(SIter_CD it(Y.EMF().Ex().y0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary 
        //Ez
        // for(SIter_CD it(Y.EMF().Ez().y0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
        //By
        for(SIter_CD it(Y.EMF().By().y0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
    }
//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Boundary conditions on both sides of a node
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::sameNode_bound_X(Stat& Y) {
//--------------------------------------------------------------
//  Choose between boundary conditions in the x direction
//--------------------------------------------------------------
        switch (BNDX()) {
            case 0:                   // periodic
                sameNode_periodic_X(Y);
                break;
            case 1:                   // mirror boundary
                sameNode_mirror_X(Y);
                break;
            default:
                cout<<"Not a valid boundary condition.\n"; 
                break;
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::sameNode_periodic_X(Stat& Y) {
//--------------------------------------------------------------
//  Periodic boundary in the x direction for 1 node
//--------------------------------------------------------------
       // Harmonics:x0 "Right-Bound ---> Left-Guard" 
        // for(int i(0); i < Y.DF().dim(); ++i)
        //     copy(Y.SH(i).x0_RB(Nbc), Y.SH(i).x0_RB(Nbc).end(), Y.SH(i).x0_LG(Nbc)); 
        // Fields:   x0 "Right-Bound ---> Left-Guard"
        for(int i(3); i < Y.EMF().dim(); ++i) // "3" as opposed to "0"
            copy(Y.FLD(i).x0_RB(Nbc), Y.FLD(i).x0_RB(Nbc).end(), Y.FLD(i).x0_LG(Nbc)); 

        // Harmonics:x0 "Left-Bound ---> Right-Guard" 
        // for(int i(0); i < Y.DF().dim(); ++i)
        //     copy(Y.SH(i).x0_LB(Nbc), Y.SH(i).x0_LB(Nbc).end(),  Y.SH(i).x0_RG(Nbc)); 
        // Fields:   x0 "Left-Bound ---> Right-Guard"
        for(int i(3); i < Y.EMF().dim(); ++i) // "3" as opposed to "0"
            copy(Y.FLD(i).x0_LB(Nbc), Y.FLD(i).x0_LB(Nbc).end(), Y.FLD(i).x0_RG(Nbc)); 
       
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::sameNode_mirror_X(Stat& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction for 1 node
//--------------------------------------------------------------
        // int sign(1);
        size_t Nx(Y.SH(0,0).numx());

        // Mirror the harmonics
        // for(int l(0); l < Y.DF().l0(); ++l){
        //     for(int m(0); m < ((Y.DF().m0() < l)? Y.DF().m0():l)+1; ++m){
        //         sign = 1-2*((l+m)%2);          //(-1)^(m+n)

                // right boundary
        //         for (int c(0); c < Nbc; ++c) 
        //             copy(Y.SH(l,m).x0(Nx-2*Nbc+c), Y.SH(l,m).x0(Nx-2*Nbc+c).end(), Y.SH(l,m).x0(Nx-c-1));
        //         for(GSIter_CD it(Y.SH(l,m).x0_RG(Nbc)); it != it.end(); ++it)  *it *= sign;

                // left boundary
        //         for (int c(0); c < Nbc; ++c) 
        //             copy(Y.SH(l,m).x0(2*Nbc-c-1), Y.SH(l,m).x0(2*Nbc-c-1).end(), Y.SH(l,m).x0(c));
        //         for(GSIter_CD it(Y.SH(l,m).x0_LG(Nbc)); it != it.end(); ++it)  *it *= sign; 
        //     }
        // }

        // Mirror the fields 
        for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
           // right boundary
           for (int c(0); c < Nbc; ++c) 
               copy(Y.FLD(i).x0(Nx-2*Nbc+c), Y.FLD(i).x0(Nx-2*Nbc+c).end(), Y.FLD(i).x0(Nx-c-1));
           // left boundary
           for (int c(0); c < Nbc; ++c) 
              copy(Y.FLD(i).x0(2*Nbc-c-1), Y.FLD(i).x0(2*Nbc-c-1).end(), Y.FLD(i).x0(c));
        }

        //Ey
        // for(GSIter_CD it(Y.EMF().Ey().x0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary 
        // for(GSIter_CD it(Y.EMF().Ey().x0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary

        //Ez
        // for(GSIter_CD it(Y.EMF().Ez().x0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
        // for(GSIter_CD it(Y.EMF().Ez().x0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary
       
        //Bx
        for(GSIter_CD it(Y.EMF().Bx().x0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
        for(GSIter_CD it(Y.EMF().Bx().x0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::sameNode_bound_Y(Stat& Y) {
//--------------------------------------------------------------
//  Choose between boundary conditions in the y direction
//--------------------------------------------------------------
        switch (BNDY()) {
            case 0:                   // periodic
                sameNode_periodic_Y(Y);
                break;
            case 1:                   // mirror
                sameNode_mirror_Y(Y);
                break;
            default:
                cout<<"Not a valid boundary condition.\n";
                break;
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::sameNode_periodic_Y(Stat& Y) {
//--------------------------------------------------------------
//  Periodic boundary for 1 node
//--------------------------------------------------------------
      
        // Harmonics:y0 "Right-Bound ---> Left-Guard" 
        // for(int i = 0; i < Y.DF().dim(); ++i)
        //     copy(Y.SH(i).y0_RB(Nbc), Y.SH(i).y0_RB(Nbc).end(), Y.SH(i).y0_LG(Nbc));
        // Fields:   y0 "Right-Bound ---> Left-Guard"
        for(int i(3); i < Y.EMF().dim(); ++i) // "3" as opposed to "0"
            copy(Y.FLD(i).y0_RB(Nbc), Y.FLD(i).y0_RB(Nbc).end(), Y.FLD(i).y0_LG(Nbc)); 

        // Harmonics:y0 "Left-Bound ---> Right-Guard" 
        // for(int i = 0; i < Y.DF().dim(); ++i)
        //     copy(Y.SH(i).y0_LB(Nbc), Y.SH(i).y0_LB(Nbc).end(), Y.SH(i).y0_RG(Nbc));
        // Fields:   y0 "Left-Bound ---> Right-Guard"
        for(int i(3); i < Y.EMF().dim(); ++i) // "3" as opposed to "0"
            copy(Y.FLD(i).y0_LB(Nbc), Y.FLD(i).y0_LB(Nbc).end(), Y.FLD(i).y0_RG(Nbc));

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_ImplicitE_Communications::sameNode_mirror_Y(Stat& Y) {
//--------------------------------------------------------------
//  mirror boundary in the "y" direction for 1 node
//--------------------------------------------------------------
        // int sign(1);
        complex< double > cmpl(1,0);
        size_t Ny(Y.SH(0,0).numy());

        // Mirror the harmonics
        // for(int l(1); l < Y.DF().l0(); ++l){
        //    for(int m(0); m < ((Y.DF().m0() < l)? Y.DF().m0():l)+1; ++m){
        //        sign = 1-2*(m%2);          //(-1)^m

                // right boundary
        //        for (int c(0); c < Nbc; ++c) 
        //            copy(Y.SH(l,m).y0(Ny-2*Nbc+c), Y.SH(l,m).y0(Ny-2*Nbc+c).end(), Y.SH(l,m).y0(Ny-c-1));
        //        for(GSIter_CD it(Y.SH(l,m).y0_RG(Nbc)); it != it.end(); ++it){
        //            cmpl  = *it; 
        //            *it  -= 2.0 * cmpl.real(); 
        //            *it *= ((-1)*sign);
        //        }

                // left boundary
        //        for (int c(0); c < Nbc; ++c) 
        //            copy(Y.SH(l,m).y0(2*Nbc-c-1), Y.SH(l,m).y0(2*Nbc-c-1).end(), Y.SH(l,m).y0(c));
        //        for(GSIter_CD it(Y.SH(l,m).y0_LG(Nbc)); it != it.end(); ++it) {
        //           cmpl  = *it; 
        //           *it  -= 2.0 * cmpl.real();
        //           *it *= ((-1)*sign); 
        //        }
        //    }
        // }

        // Mirror the fields 
        for(int i(3); i < Y.EMF().dim(); ++i){ // "3" as opposed to "0"
           // right boundary
           for (int c(0); c < Nbc; ++c) 
               copy(Y.FLD(i).y0(Ny-2*Nbc+c), Y.FLD(i).y0(Ny-2*Nbc+c).end(), Y.FLD(i).y0(Ny-c-1));
           // left boundary
           for (int c(0); c < Nbc; ++c) 
              copy(Y.FLD(i).y0(2*Nbc-c-1), Y.FLD(i).y0(2*Nbc-c-1).end(), Y.FLD(i).y0(c));
        }

        //Ex
        // for(SIter_CD it(Y.EMF().Ex().y0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary 
        // for(SIter_CD it(Y.EMF().Ex().y0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary

        //Ez
        // for(SIter_CD it(Y.EMF().Ez().y0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
        // for(SIter_CD it(Y.EMF().Ez().y0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary
       
        //By
        for(SIter_CD it(Y.EMF().By().y0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
        for(SIter_CD it(Y.EMF().By().y0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary

    }
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//**************************************************************
//  Definition of the Nodes Communications class
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Node_Communications:: Node_Communications() : 
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        Nbc(Inputdata::IN().list().RKLevel),      // # of boundary cells
        bndX(Inputdata::IN().list().bndX),        // Type of boundary in X
        bndY(Inputdata::IN().list().bndY){        // Type of boundary in Y
         
        using Inputdata::IN; 

        // # of harmonics
        msg_sizeX = ((IN().inp().m0+1)*(2*IN().inp().l0-IN().inp().m0+2))/2; 
        // (# of harmonics) * (# cells in p)
        msg_sizeX *= IN().inp().pr.dim();
        // 6 fields: Ex, Ey, Ez, Bx, By, Bz
        msg_sizeX += 6;  
        msg_sizeY = msg_sizeX; 
        
        msg_sizeX *= (IN().inp().y.dim()*Nbc);   
        msg_sizeY *= (IN().inp().x.dim()*Nbc);  
        msg_bufX = new complex<double>[msg_sizeX];
        msg_bufY = new complex<double>[msg_sizeY];
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    Node_Communications:: ~Node_Communications(){
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
        delete[] msg_bufX;
        delete[] msg_bufY;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    int Node_Communications:: BNDX()  const {return bndX;} 
    int Node_Communications:: BNDY()  const {return bndY;} 
//--------------------------------------------------------------

//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Send and receive in the X direction
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_Communications::Send_right_X(Stat& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the right boundary and send them 
//           to the node on the right
//--------------------------------------------------------------
        static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().y.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().y.dim()*Nbc);
        size_t bufind(0);

        // Harmonics:x0 "Right-Bound ---> " 
        for(size_t i = 0; i < Y.DF().dim(); ++i){
            copy(Y.SH(i).x0_RB(Nbc), Y.SH(i).x0_RB(Nbc).end(), msg_bufX + bufind); 
            bufind += step_h;
        }
        // Fields:   x0 "Right-Bound --> "
        for(size_t i = 0; i < Y.EMF().dim(); ++i){
            copy(Y.FLD(i).x0_RB(Nbc), Y.FLD(i).x0_RB(Nbc).end(), msg_bufX + bufind);
            bufind += step_f;
        } 

        MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::Recv_from_left_X(Stat& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the left and update
//           the left guard cells
//--------------------------------------------------------------
        static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().y.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().y.dim()*Nbc);
        size_t bufind(0);
        MPI_Status status; 

        // Receive Data
        MPI_Recv(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &status);

        // Harmonics:x0-"---> Left-Guard" 
        for(int i = 0; i < Y.DF().dim(); ++i){
            copy(msg_bufX + bufind, msg_bufX + bufind + step_h, Y.SH(i).x0_LG(Nbc)); 
            bufind += step_h;
        }
        // Fields:   x0-"---> Left-Guard"
        for(int i = 0; i < Y.EMF().dim(); ++i){
            copy(msg_bufX + bufind, msg_bufX + bufind + step_f, Y.FLD(i).x0_LG(Nbc)); 
            bufind += step_f;
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::Send_left_X(Stat& Y, int dest) {
//--------------------------------------------------------------
//  X-axis : Read data from the left boundary and send them 
//           to the node on the left 
//--------------------------------------------------------------
        static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().y.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().y.dim()*Nbc);
        size_t bufind(0); 

        // Harmonics:x0 " <--- Left-Bound "
        for(int i = 0; i < Y.DF().dim(); ++i){
            copy(Y.SH(i).x0_LB(Nbc), Y.SH(i).x0_LB(Nbc).end(), msg_bufX + bufind);
            bufind += step_h;
        } 
        // Fields:   x0 " <--- Left-Bound "
        for(int i = 0; i < Y.EMF().dim(); ++i){
            copy(Y.FLD(i).x0_LB(Nbc), Y.FLD(i).x0_LB(Nbc).end(), msg_bufX + bufind); 
            bufind += step_f;
        }

        MPI_Send(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, dest, 1, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::Recv_from_right_X(Stat& Y, int origin) {
//--------------------------------------------------------------
//  X-axis : Receive data from the node on the right and update
//           the right guard cells
//--------------------------------------------------------------
        static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().y.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().y.dim()*Nbc);
        size_t bufind(0);
        MPI_Status status; 

        // Receive Data
        MPI_Recv(msg_bufX, msg_sizeX, MPI_DOUBLE_COMPLEX, origin, 1, MPI_COMM_WORLD, &status);

        // Harmonics:x0-"Right-Guard <--- " 
        for(int i = 0; i < Y.DF().dim(); ++i){
            copy(msg_bufX + bufind, msg_bufX + bufind + step_h, Y.SH(i).x0_RG(Nbc)); 
            bufind += step_h;
        }
        // Fields:   x0-"Right-Guard <--- "
        for(int i = 0; i < Y.EMF().dim(); ++i){
            copy(msg_bufX + bufind, msg_bufX + bufind + step_f, Y.FLD(i).x0_RG(Nbc)); 
            bufind += step_f;
        }
    }
//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Send and receive in the Y direction
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_Communications::Send_right_Y(Stat& Y, int dest) {
//--------------------------------------------------------------
//  Y-axis : Read data from the right boundary and send them 
//           to the node on the right
//--------------------------------------------------------------
        static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().x.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().x.dim()*Nbc);
        size_t bufind(0); 

       // Harmonics:x0 "Right-Bound ---> " 
        for(int i = 0; i < Y.DF().dim(); ++i){
            copy(Y.SH(i).y0_RB(Nbc), Y.SH(i).y0_RB(Nbc).end(), msg_bufY + bufind); 
            bufind += step_h;
        }
        // Fields:   x0 "Right-Bound --> "
        for(int i = 0; i < Y.EMF().dim(); ++i){
            copy(Y.FLD(i).y0_RB(Nbc), Y.FLD(i).y0_RB(Nbc).end(), msg_bufY + bufind); 
            bufind += step_f;
        }

        MPI_Send(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, dest, 0, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::Recv_from_left_Y(Stat& Y, int origin) {
//--------------------------------------------------------------
//  Y-axis : Receive data from the node on the left and update
//           the left guard cells
//--------------------------------------------------------------
        static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().x.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().x.dim()*Nbc);
        size_t bufind(0);
        MPI_Status status;
        // Receive Data
        MPI_Recv(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, origin, 0, MPI_COMM_WORLD, &status);

        // Harmonics:x0-"---> Left-Guard" 
        for(int i = 0; i < Y.DF().dim(); ++i){
            copy(msg_bufY + bufind, msg_bufY + bufind + step_h, Y.SH(i).y0_LG(Nbc)); 
            bufind += step_h;
        }
        // Fields:   x0-"---> Left-Guard"
        for(int i = 0; i < Y.EMF().dim(); ++i){
            copy(msg_bufY + bufind, msg_bufY + bufind + step_f, Y.FLD(i).y0_LG(Nbc)); 
            bufind += step_f;
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::Send_left_Y(Stat& Y, int dest) {
//--------------------------------------------------------------
//  Y-axis : Read data from the left boundary and send them 
//           to the node on the left 
//--------------------------------------------------------------
        static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().x.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().x.dim()*Nbc);
        size_t bufind(0); 

        // Harmonics:x0 " <--- Left-Bound "
        for(int i = 0; i < Y.DF().dim(); ++i) {
            copy(Y.SH(i).y0_LB(Nbc),Y.SH(i).y0_LB(Nbc).end(), msg_bufY + bufind); 
            bufind += step_h;
        }
        // Fields:   x0 " <--- Left-Bound "
        for(int i = 0; i < Y.EMF().dim(); ++i) {
            copy(Y.FLD(i).y0_LB(Nbc), Y.FLD(i).y0_LB(Nbc).end(), msg_bufY + bufind); 
            bufind += step_f;
        }

        MPI_Send(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, dest, 1, MPI_COMM_WORLD);
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::Recv_from_right_Y(Stat& Y, int origin) {
//--------------------------------------------------------------
//  Y-axis : Receive data from the node on the right and update
//           the right guard cells
//--------------------------------------------------------------
        static size_t step_h(Inputdata::IN().inp().pr.dim()*Inputdata::IN().inp().x.dim()*Nbc);
        static size_t step_f(Inputdata::IN().inp().x.dim()*Nbc);
        size_t bufind(0); 
        MPI_Status status; 
        // Receive Data
        MPI_Recv(msg_bufY, msg_sizeY, MPI_DOUBLE_COMPLEX, origin, 1, MPI_COMM_WORLD, &status);

        // Harmonics:x0-"Right-Guard <--- " 
        for(int i = 0; i < Y.DF().dim(); ++i) {
            copy(msg_bufY + bufind, msg_bufY + bufind + step_h, Y.SH(i).y0_RG(Nbc)); 
            bufind += step_h;
        }
        // Fields:   x0-"Right-Guard <--- "
        for(int i = 0; i < Y.EMF().dim(); ++i) {
            copy(msg_bufY + bufind, msg_bufY + bufind + step_f, Y.FLD(i).y0_RG(Nbc)); 
            bufind += step_f;
        }
    }
//--------------------------------------------------------------


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  Boundary conditions
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Mirror conditions on one side
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//**************************************************************
//--------------------------------------------------------------
    void Node_Communications::mirror_bound_Xleft(Stat& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction on the left
//--------------------------------------------------------------
        int sign(1);
        size_t Nx(Y.SH(0,0).numx());

        // Mirror the harmonics
        for(int l(0); l < Y.DF().l0(); ++l){
            for(int m(0); m < ((Y.DF().m0() < l)? Y.DF().m0():l)+1; ++m){
                sign = 1-2*((l+m)%2);          //(-1)^(m+n)

                for (int c(0); c < Nbc; ++c)
                    copy(Y.SH(l,m).x0(2*Nbc-c-1), Y.SH(l,m).x0(2*Nbc-c-1).end(), Y.SH(l,m).x0(c));
                for(GSIter_CD it(Y.SH(l,m).x0_LG(Nbc)); it != it.end(); ++it)  *it *= sign; 
            }
        }

        // Mirror the fields 
        for(int i(0); i < Y.EMF().dim(); ++i){
           for (int c(0); c < Nbc; ++c) 
              copy(Y.FLD(i).x0(2*Nbc-c-1), Y.FLD(i).x0(2*Nbc-c-1).end(), Y.FLD(i).x0(c));
        }

        //Ey
        for(GSIter_CD it(Y.EMF().Ey().x0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary
        //Ez
        for(GSIter_CD it(Y.EMF().Ez().x0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary
        //Bx
        for(GSIter_CD it(Y.EMF().Bx().x0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::mirror_bound_Xright(Stat& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction on the right
//--------------------------------------------------------------
        int sign(1);
        size_t Nx(Y.SH(0,0).numx());

        // Mirror the harmonics
        for(int l(0); l < Y.DF().l0(); ++l){
            for(int m(0); m < ((Y.DF().m0() < l)? Y.DF().m0():l)+1; ++m){
                sign = 1-2*((l+m)%2);          //(-1)^(m+n)

                for (int c(0); c < Nbc; ++c) 
                    copy(Y.SH(l,m).x0(Nx-2*Nbc+c), Y.SH(l,m).x0(Nx-2*Nbc+c).end(), Y.SH(l,m).x0(Nx-c-1));
                for(GSIter_CD it(Y.SH(l,m).x0_RG(Nbc)); it != it.end(); ++it)  *it *= sign;
 
            }
        }

        // Mirror the fields 
        for(int i(0); i < Y.EMF().dim(); ++i){
           for (int c(0); c < Nbc; ++c) 
               copy(Y.FLD(i).x0(Nx-2*Nbc+c), Y.FLD(i).x0(Nx-2*Nbc+c).end(),  Y.FLD(i).x0(Nx-c-1));
        }

        //Ey
        for(GSIter_CD it(Y.EMF().Ey().x0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary 
        //Ez
        for(GSIter_CD it(Y.EMF().Ez().x0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
        //Bx
        for(GSIter_CD it(Y.EMF().Bx().x0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::mirror_bound_Yleft(Stat& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the y direction on the left
//--------------------------------------------------------------
        
        int sign(1);
        complex< double > cmpl(1,0);
        size_t Ny(Y.SH(0,0).numy());
        // Mirror the harmonics
        for(int l(1); l < Y.DF().l0(); ++l){
            for(int m(0); m < ((Y.DF().m0() < l)? Y.DF().m0():l)+1; ++m){
                sign = 1-2*(m%2);          //(-1)^m

                // left boundary
                for (int c(0); c < Nbc; ++c) 
                    copy(Y.SH(l,m).y0(2*Nbc-c-1), Y.SH(l,m).y0(2*Nbc-c-1).end(), Y.SH(l,m).y0(c));
                for(GSIter_CD it(Y.SH(l,m).y0_LG(Nbc)); it != it.end(); ++it) {
                   cmpl  = *it; 
                   *it  -= 2.0 * cmpl.real();
                   *it *= ((-1)*sign); 
                }
            }
        }

        // Mirror the fields 
        for(int i(0); i < Y.EMF().dim(); ++i){
           for (int c(0); c < Nbc; ++c) 
              copy(Y.FLD(i).y0(2*Nbc-c-1), Y.FLD(i).y0(2*Nbc-c-1).end(), Y.FLD(i).y0(c));
        }

        //Ex
        for(SIter_CD it(Y.EMF().Ex().y0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary
        //Ez
        for(SIter_CD it(Y.EMF().Ez().y0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary
        //By
        for(SIter_CD it(Y.EMF().By().y0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::mirror_bound_Yright(Stat& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the y direction on the right
//--------------------------------------------------------------
        int sign(1);
        complex< double > cmpl(1,0);
        size_t Ny(Y.SH(0,0).numy());
        // Mirror the harmonics
        for(int l(1); l < Y.DF().l0(); ++l){
            for(int m(0); m < ((Y.DF().m0() < l)? Y.DF().m0():l)+1; ++m){
                sign = 1-2*(m%2);          //(-1)^m

                // right boundary
                for (int c(0); c < Nbc; ++c) 
                    copy(Y.SH(l,m).y0(Ny-2*Nbc+c), Y.SH(l,m).y0(Ny-2*Nbc+c).end(), Y.SH(l,m).y0(Ny-c-1));
                for(GSIter_CD it(Y.SH(l,m).y0_RG(Nbc)); it != it.end(); ++it){
                    cmpl  = *it; 
                    *it  -= 2.0 * cmpl.real(); 
                    *it *= ((-1)*sign);
                }
            }
        }

        // Mirror the fields 
        for(int i(0); i < Y.EMF().dim(); ++i){
           for (int c(0); c < Nbc; ++c) 
               copy(Y.FLD(i).y0(Ny-2*Nbc+c), Y.FLD(i).y0(Ny-2*Nbc+c).end(), Y.FLD(i).y0(Ny-c-1));
        }

        //Ex
        for(SIter_CD it(Y.EMF().Ex().y0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary 
        //Ez
        for(SIter_CD it(Y.EMF().Ez().y0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
        //By
        for(SIter_CD it(Y.EMF().By().y0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
    }
//--------------------------------------------------------------

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Boundary conditions on both sides of a node
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//--------------------------------------------------------------
    void Node_Communications::sameNode_bound_X(Stat& Y) {
//--------------------------------------------------------------
//  Choose between boundary conditions in the x direction
//--------------------------------------------------------------
        switch (BNDX()) {
            case 0:                   // periodic
                sameNode_periodic_X(Y);
                break;
            case 1:                   // mirror boundary
                sameNode_mirror_X(Y);
                break;
            default:
                cout<<"Not a valid boundary condition.\n"; 
                break;
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::sameNode_periodic_X(Stat& Y) {
//--------------------------------------------------------------
//  Periodic boundary in the x direction for 1 node
//--------------------------------------------------------------
       // Harmonics:x0 "Right-Bound ---> Left-Guard" 
        for(int i(0); i < Y.DF().dim(); ++i)
            copy(Y.SH(i).x0_RB(Nbc), Y.SH(i).x0_RB(Nbc).end(), Y.SH(i).x0_LG(Nbc)); 
        // Fields:   x0 "Right-Bound ---> Left-Guard"
        for(int i(0); i < Y.EMF().dim(); ++i)
            copy(Y.FLD(i).x0_RB(Nbc), Y.FLD(i).x0_RB(Nbc).end(), Y.FLD(i).x0_LG(Nbc)); 

        // Harmonics:x0 "Left-Bound ---> Right-Guard" 
        for(int i(0); i < Y.DF().dim(); ++i)
            copy(Y.SH(i).x0_LB(Nbc), Y.SH(i).x0_LB(Nbc).end(),  Y.SH(i).x0_RG(Nbc)); 
        // Fields:   x0 "Left-Bound ---> Right-Guard"
        for(int i(0); i < Y.EMF().dim(); ++i)
            copy(Y.FLD(i).x0_LB(Nbc), Y.FLD(i).x0_LB(Nbc).end(), Y.FLD(i).x0_RG(Nbc)); 
       
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::sameNode_mirror_X(Stat& Y) {
//--------------------------------------------------------------
//  Mirror boundary in the x direction for 1 node
//--------------------------------------------------------------
        int sign(1);
        size_t Nx(Y.SH(0,0).numx());

        // Mirror the harmonics
        for(int l(0); l < Y.DF().l0(); ++l){
            for(int m(0); m < ((Y.DF().m0() < l)? Y.DF().m0():l)+1; ++m){
                sign = 1-2*((l+m)%2);          //(-1)^(m+n)

                // right boundary
                for (int c(0); c < Nbc; ++c) 
                    copy(Y.SH(l,m).x0(Nx-2*Nbc+c), Y.SH(l,m).x0(Nx-2*Nbc+c).end(), Y.SH(l,m).x0(Nx-c-1));
                for(GSIter_CD it(Y.SH(l,m).x0_RG(Nbc)); it != it.end(); ++it)  *it *= sign;

                // left boundary
                for (int c(0); c < Nbc; ++c) 
                    copy(Y.SH(l,m).x0(2*Nbc-c-1), Y.SH(l,m).x0(2*Nbc-c-1).end(), Y.SH(l,m).x0(c));
                for(GSIter_CD it(Y.SH(l,m).x0_LG(Nbc)); it != it.end(); ++it)  *it *= sign; 
            }
        }

        // Mirror the fields 
        for(int i(0); i < Y.EMF().dim(); ++i){
           // right boundary
           for (int c(0); c < Nbc; ++c) 
               copy(Y.FLD(i).x0(Nx-2*Nbc+c), Y.FLD(i).x0(Nx-2*Nbc+c).end(), Y.FLD(i).x0(Nx-c-1));
           // left boundary
           for (int c(0); c < Nbc; ++c) 
              copy(Y.FLD(i).x0(2*Nbc-c-1), Y.FLD(i).x0(2*Nbc-c-1).end(), Y.FLD(i).x0(c));
        }

        //Ey
        for(GSIter_CD it(Y.EMF().Ey().x0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary 
        for(GSIter_CD it(Y.EMF().Ey().x0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary

        //Ez
        for(GSIter_CD it(Y.EMF().Ez().x0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
        for(GSIter_CD it(Y.EMF().Ez().x0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary
       
        //Bx
        for(GSIter_CD it(Y.EMF().Bx().x0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
        for(GSIter_CD it(Y.EMF().Bx().x0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::sameNode_bound_Y(Stat& Y) {
//--------------------------------------------------------------
//  Choose between boundary conditions in the y direction
//--------------------------------------------------------------
        switch (BNDY()) {
            case 0:                   // periodic
                sameNode_periodic_Y(Y);
                break;
            case 1:                   // mirror
                sameNode_mirror_Y(Y);
                break;
            default:
                cout<<"Not a valid boundary condition.\n";
                break;
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::sameNode_periodic_Y(Stat& Y) {
//--------------------------------------------------------------
//  Periodic boundary for 1 node
//--------------------------------------------------------------
      
        // Harmonics:y0 "Right-Bound ---> Left-Guard" 
        for(int i = 0; i < Y.DF().dim(); ++i)
            copy(Y.SH(i).y0_RB(Nbc), Y.SH(i).y0_RB(Nbc).end(), Y.SH(i).y0_LG(Nbc));
        // Fields:   y0 "Right-Bound ---> Left-Guard"
        for(int i = 0; i < Y.EMF().dim(); ++i)
            copy(Y.FLD(i).y0_RB(Nbc), Y.FLD(i).y0_RB(Nbc).end(), Y.FLD(i).y0_LG(Nbc)); 

        // Harmonics:y0 "Left-Bound ---> Right-Guard" 
        for(int i = 0; i < Y.DF().dim(); ++i)
            copy(Y.SH(i).y0_LB(Nbc), Y.SH(i).y0_LB(Nbc).end(), Y.SH(i).y0_RG(Nbc));
        // Fields:   y0 "Left-Bound ---> Right-Guard"
        for(int i = 0; i < Y.EMF().dim(); ++i)
            copy(Y.FLD(i).y0_LB(Nbc), Y.FLD(i).y0_LB(Nbc).end(), Y.FLD(i).y0_RG(Nbc));

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Node_Communications::sameNode_mirror_Y(Stat& Y) {
//--------------------------------------------------------------
//  mirror boundary in the "y" direction for 1 node
//--------------------------------------------------------------
        int sign(1);
        complex< double > cmpl(1,0);
        size_t Ny(Y.SH(0,0).numy());

        // Mirror the harmonics
        for(int l(1); l < Y.DF().l0(); ++l){
            for(int m(0); m < ((Y.DF().m0() < l)? Y.DF().m0():l)+1; ++m){
                sign = 1-2*(m%2);          //(-1)^m

                // right boundary
                for (int c(0); c < Nbc; ++c) 
                    copy(Y.SH(l,m).y0(Ny-2*Nbc+c), Y.SH(l,m).y0(Ny-2*Nbc+c).end(), Y.SH(l,m).y0(Ny-c-1));
                for(GSIter_CD it(Y.SH(l,m).y0_RG(Nbc)); it != it.end(); ++it){
                    cmpl  = *it; 
                    *it  -= 2.0 * cmpl.real(); 
                    *it *= ((-1)*sign);
                }

                // left boundary
                for (int c(0); c < Nbc; ++c) 
                    copy(Y.SH(l,m).y0(2*Nbc-c-1), Y.SH(l,m).y0(2*Nbc-c-1).end(), Y.SH(l,m).y0(c));
                for(GSIter_CD it(Y.SH(l,m).y0_LG(Nbc)); it != it.end(); ++it) {
                   cmpl  = *it; 
                   *it  -= 2.0 * cmpl.real();
                   *it *= ((-1)*sign); 
                }
            }
        }

        // Mirror the fields 
        for(int i(0); i < Y.EMF().dim(); ++i){
           // right boundary
           for (int c(0); c < Nbc; ++c) 
               copy(Y.FLD(i).y0(Ny-2*Nbc+c), Y.FLD(i).y0(Ny-2*Nbc+c).end(), Y.FLD(i).y0(Ny-c-1));
           // left boundary
           for (int c(0); c < Nbc; ++c) 
              copy(Y.FLD(i).y0(2*Nbc-c-1), Y.FLD(i).y0(2*Nbc-c-1).end(), Y.FLD(i).y0(c));
        }

        //Ex
        for(SIter_CD it(Y.EMF().Ex().y0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary 
        for(SIter_CD it(Y.EMF().Ex().y0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary

        //Ez
        for(SIter_CD it(Y.EMF().Ez().y0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
        for(SIter_CD it(Y.EMF().Ez().y0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary
       
        //By
        for(SIter_CD it(Y.EMF().By().y0_RG(Nbc)); it != it.end(); ++it) *it *= -1.0; // right boundary
        for(SIter_CD it(Y.EMF().By().y0_LG(Nbc)); it != it.end(); ++it) *it *= -1.0; // left  boundary

    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//**************************************************************
//  Definition Parallel Output
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Parallel_Output:: Parallel_Output():
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        Nbc(Inputdata::IN().list().RKLevel),       // # boundary cells                               
        //
        xglob_axis(Inputdata::IN().inp().xglob),   // total x-axis
        yglob_axis(Inputdata::IN().inp().yglob),   // total y-axis
        pr(Inputdata::IN().inp().pr),              // p-radius axis   
        //
        NnodesX(Inputdata::IN().list().NnodesX),   // Number of nodes in X-direction
        //
        numpx(Inputdata::IN().outp().px.dim()),    // numpx cells
        nump1(Inputdata::IN().outp().p1.dim()),    // nump1 cells
        nump2(Inputdata::IN().outp().p2.dim()),    // nump2 cells
        nump3(Inputdata::IN().outp().p3.dim()){    // nump3 cells
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        using Inputdata::IN; 
        szx = IN().inp().x.dim()-2*Nbc;      // size of useful x axis 
        szy = IN().inp().y.dim()-2*Nbc;      // size of useful y axis
        Nnodes = IN().list().NnodesX * IN().list().NnodesY;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    Parallel_Output:: ~Parallel_Output(){ }
//--------------------------------------------------------------

//--------------------------------------------------------------
    int Parallel_Output:: RANK()  const {return rank;} 
    int Parallel_Output:: NODES() const {return Nnodes;} 
    int Parallel_Output:: NODESX() const {return NnodesX;} 
//--------------------------------------------------------------

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Parallel Output for Fields
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//--------------------------------------------------------------
    void Parallel_Output::parallel_Ex(size_t step, Stat& Y) {
//--------------------------------------------------------------
//  Parallel output for Ex
//--------------------------------------------------------------
        MPI_Status status;
        size_t st(0), bi(0);
        int msg_sz(szx*szy);
        float* Exbuf = new float[msg_sz];
        Matrix2D<float> MExglob(xglob_axis.dim(), yglob_axis.dim());

        for(GSIter_CD it(Y.EMF().Ex().TrimGuards(Nbc)); it != it.end(); ++it) 
            Exbuf[bi++] = static_cast<float>((*it).real());

        if (NODES() > 1) { 
            if (RANK()!=0) {
               MPI_Send(Exbuf, msg_sz, MPI_FLOAT, 0, RANK(), MPI_COMM_WORLD);
            } 
            else {            
                // Fill data for rank = 0 
                bi = 0;
                for(GSIter_F it(MExglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                    *it = Exbuf[bi++];
                // Fill data for rank > 0 
                for (int rr = 1; rr < NODES(); ++rr){
                    MPI_Recv(Exbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    st = (rr / NODESX())*szy*xglob_axis.dim()+(rr % NODESX())*szx;
                    bi = 0;
                    for(GSIter_F it(MExglob.SubMatrix2D(st, szx, szy)); it != it.end(); ++it) 
                        *it = Exbuf[bi++]; 
                 }
             }
         }
         // Fill data for Nodes = 0 
         else { 
             bi = 0;
             for(GSIter_F it(MExglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                 *it = Exbuf[bi++];
         }
         if (RANK() == 0) exportdata.Export_Ex(step, MExglob);
         delete[] Exbuf;
     }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Parallel_Output::parallel_Ey(size_t step, Stat& Y) {
//--------------------------------------------------------------
//  Parallel output for Ey
//--------------------------------------------------------------
        MPI_Status status;
        size_t st(0), bi(0);
        int msg_sz(szx*szy);
        float* Eybuf = new float[msg_sz];
        Matrix2D<float> MEyglob(xglob_axis.dim(), yglob_axis.dim());

        for(GSIter_CD it(Y.EMF().Ey().TrimGuards(Nbc)); it != it.end(); ++it) 
            Eybuf[bi++] = static_cast<float>((*it).real());

        if (NODES() > 1) { 
            if (RANK()!=0) {
               MPI_Send(Eybuf, msg_sz, MPI_FLOAT, 0, RANK(), MPI_COMM_WORLD);
            } 
            else {            
                // Fill data for rank = 0 
                bi = 0;
                for(GSIter_F it(MEyglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                    *it = Eybuf[bi++];
                // Fill data for rank > 0 
                for (int rr = 1; rr < NODES(); ++rr){
                    MPI_Recv(Eybuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    st = (rr / NODESX())*szy*xglob_axis.dim()+(rr % NODESX())*szx;
                    bi = 0;
                    for(GSIter_F it(MEyglob.SubMatrix2D(st, szx, szy)); it != it.end(); ++it) 
                        *it = Eybuf[bi++]; 
                 }
             }
         }
         // Fill data for Nodes = 0 
         else { 
             bi = 0;
             for(GSIter_F it(MEyglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                 *it = Eybuf[bi++];
         }
         if (RANK() == 0) exportdata.Export_Ey(step, MEyglob);
         delete[] Eybuf;
     }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Parallel_Output::parallel_Ez(size_t step, Stat& Y) {
//--------------------------------------------------------------
//  Parallel output for Ez
//--------------------------------------------------------------
        MPI_Status status;
        size_t st(0), bi(0);
        int msg_sz(szx*szy);
        float* Ezbuf = new float[msg_sz];
        Matrix2D<float> MEzglob(xglob_axis.dim(), yglob_axis.dim());

        for(GSIter_CD it(Y.EMF().Ez().TrimGuards(Nbc)); it != it.end(); ++it) 
            Ezbuf[bi++] = static_cast<float>((*it).real());

        if (NODES() > 1) { 
            if (RANK()!=0) {
               MPI_Send(Ezbuf, msg_sz, MPI_FLOAT, 0, RANK(), MPI_COMM_WORLD);
            } 
            else {            
                // Fill data for rank = 0 
                bi = 0;
                for(GSIter_F it(MEzglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                    *it = Ezbuf[bi++];
                // Fill data for rank > 0 
                for (int rr = 1; rr < NODES(); ++rr){
                    MPI_Recv(Ezbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    st = (rr / NODESX())*szy*xglob_axis.dim()+(rr % NODESX())*szx;
                    bi = 0;
                    for(GSIter_F it(MEzglob.SubMatrix2D(st, szx, szy)); it != it.end(); ++it) 
                        *it = Ezbuf[bi++]; 
                 }
             }
         }
         // Fill data for Nodes = 0 
         else { 
             bi = 0;
             for(GSIter_F it(MEzglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                 *it = Ezbuf[bi++];
         }
         if (RANK() == 0) exportdata.Export_Ez(step, MEzglob);
         delete[] Ezbuf;
     }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Parallel_Output::parallel_Bx(size_t step, Stat& Y) {
//--------------------------------------------------------------
//  Parallel output for Bx
//--------------------------------------------------------------
        MPI_Status status;
        size_t st(0), bi(0);
        int msg_sz(szx*szy);
        float* Bxbuf = new float[msg_sz];
        Matrix2D<float> MBxglob(xglob_axis.dim(), yglob_axis.dim());

        for(GSIter_CD it(Y.EMF().Bx().TrimGuards(Nbc)); it != it.end(); ++it) 
            Bxbuf[bi++] = static_cast<float>((*it).real());

        if (NODES() > 1) { 
            if (RANK()!=0) {
               MPI_Send(Bxbuf, msg_sz, MPI_FLOAT, 0, RANK(), MPI_COMM_WORLD);
            } 
            else {            
                // Fill data for rank = 0 
                bi = 0;
                for(GSIter_F it(MBxglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                    *it = Bxbuf[bi++];
                // Fill data for rank > 0 
                for (int rr = 1; rr < NODES(); ++rr){
                    MPI_Recv(Bxbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    st = (rr / NODESX())*szy*xglob_axis.dim()+(rr % NODESX())*szx;
                    bi = 0;
                    for(GSIter_F it(MBxglob.SubMatrix2D(st, szx, szy)); it != it.end(); ++it) 
                        *it = Bxbuf[bi++]; 
                 }
             }
         }
         // Fill data for Nodes = 0 
         else { 
             bi = 0;
             for(GSIter_F it(MBxglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                 *it = Bxbuf[bi++];
         }
         if (RANK() == 0) exportdata.Export_Bx(step, MBxglob);
         delete[] Bxbuf;
     }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Parallel_Output::parallel_By(size_t step, Stat& Y) {
//--------------------------------------------------------------
//  Parallel output for By
//--------------------------------------------------------------
        MPI_Status status;
        size_t st(0), bi(0);
        int msg_sz(szx*szy);
        float* Bybuf = new float[msg_sz];
        Matrix2D<float> MByglob(xglob_axis.dim(), yglob_axis.dim());

        for(GSIter_CD it(Y.EMF().By().TrimGuards(Nbc)); it != it.end(); ++it) 
            Bybuf[bi++] = static_cast<float>((*it).real());

        if (NODES() > 1) { 
            if (RANK()!=0) {
               MPI_Send(Bybuf, msg_sz, MPI_FLOAT, 0, RANK(), MPI_COMM_WORLD);
            } 
            else {            
                // Fill data for rank = 0 
                bi = 0;
                for(GSIter_F it(MByglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                    *it = Bybuf[bi++];
                // Fill data for rank > 0 
                for (int rr = 1; rr < NODES(); ++rr){
                    MPI_Recv(Bybuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    st = (rr / NODESX())*szy*xglob_axis.dim()+(rr % NODESX())*szx;
                    bi = 0;
                    for(GSIter_F it(MByglob.SubMatrix2D(st, szx, szy)); it != it.end(); ++it) 
                        *it = Bybuf[bi++]; 
                 }
             }
         }
         // Fill data for Nodes = 0 
         else { 
             bi = 0;
             for(GSIter_F it(MByglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                 *it = Bybuf[bi++];
         }
         if (RANK() == 0) exportdata.Export_By(step, MByglob);
         delete[] Bybuf;
     }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Parallel_Output::parallel_Bz(size_t step, Stat& Y) {
//--------------------------------------------------------------
//  Parallel output for Bz
//--------------------------------------------------------------
        MPI_Status status;
        size_t st(0), bi(0);
        int msg_sz(szx*szy);
        float* Bzbuf = new float[msg_sz];
        Matrix2D<float> MBzglob(xglob_axis.dim(), yglob_axis.dim());

        for(GSIter_CD it(Y.EMF().Bz().TrimGuards(Nbc)); it != it.end(); ++it) 
            Bzbuf[bi++] = static_cast<float>((*it).real());

        if (NODES() > 1) { 
            if (RANK()!=0) {
               MPI_Send(Bzbuf, msg_sz, MPI_FLOAT, 0, RANK(), MPI_COMM_WORLD);
            } 
            else {            
                // Fill data for rank = 0 
                bi = 0;
                for(GSIter_F it(MBzglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                    *it = Bzbuf[bi++];
                // Fill data for rank > 0 
                for (int rr = 1; rr < NODES(); ++rr){
                    MPI_Recv(Bzbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    st = (rr / NODESX())*szy*xglob_axis.dim()+(rr % NODESX())*szx;
                    bi = 0;
                    for(GSIter_F it(MBzglob.SubMatrix2D(st, szx, szy)); it != it.end(); ++it) 
                        *it = Bzbuf[bi++]; 
                 }
             }
         }
         // Fill data for Nodes = 0 
         else { 
             bi = 0;
             for(GSIter_F it(MBzglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                 *it = Bzbuf[bi++];
         }
         if (RANK() == 0) exportdata.Export_Bz(step, MBzglob);
         delete[] Bzbuf;
     }
//--------------------------------------------------------------

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Parallel Output for the moments of the distribution function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//--------------------------------------------------------------
     void Parallel_Output::parallel_mom(  size_t step, Matrix2D<float>& moment, string& filename){
//--------------------------------------------------------------
//   Parallel output for the moments of the distribution function
//--------------------------------------------------------------
         MPI_Status status;
         size_t st(0), bi(0);
         int msg_sz(szx*szy);
         float* Mbuf = new float[msg_sz];
         Matrix2D<float> Mglob(xglob_axis.dim(), yglob_axis.dim());

         for(int i(0); i < moment.dim(); ++i) Mbuf[i] = moment(i);

         if (NODES() > 1) { 
            if (RANK()!=0) {
               MPI_Send(Mbuf, msg_sz, MPI_FLOAT, 0, RANK(), MPI_COMM_WORLD);
            } 
            else {            
                // Fill data for rank = 0 
                bi = 0;
                for(GSIter_F it(Mglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                    *it = Mbuf[bi++];
                // Fill data for rank > 0 
                for (int rr(1); rr < NODES(); ++rr){
                    MPI_Recv(Mbuf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    st = (rr / NODESX())*szy*xglob_axis.dim()+(rr % NODESX())*szx;
                    bi = 0;
                    for(GSIter_F it(Mglob.SubMatrix2D(st, szx, szy)); it != it.end(); ++it) 
                        *it = Mbuf[bi++]; 
                 }
             }
         }
         // Fill data for Nodes = 0 
         else { 
             bi = 0;
             for(GSIter_F it(Mglob.SubMatrix2D(0, szx, szy)); it != it.end(); ++it) 
                 *it = Mbuf[bi++];
         }
         if (RANK() == 0) exportdata.Export_mom(step, Mglob, filename);
         delete[] Mbuf;
         
     }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Parallel_Output::nonrelativistic_moments(size_t step, Stat& Y) {
//--------------------------------------------------------------
//  Calculation of the output matrices related to the moments
//  of the spherical harmonics
//--------------------------------------------------------------

        SHarmonic&  F0(Y.SH(0,0)),
                    F10(Y.SH(1,0)),
                    F11(Y.SH(1,1)),
                    F20(Y.SH(2,0)),
                    F21(Y.SH(2,1)),
                    F22(Y.SH(2,2));

        complex<double> p0p1_sq(pr(0)*pr(0)/(pr(1)*pr(1))),
                        inv_mp0p1_sq(1.0/(1.0-p0p1_sq)),
                        f00;

        Axis< double > gamma(pr), invgamma(pr);
        for (size_t ip(0); ip < pr.dim(); ++ip) gamma(ip)    = sqrt(1.0+pr(ip)*pr(ip)); 
        for (size_t ip(0); ip < pr.dim(); ++ip) invgamma(ip) = 1.0 / gamma(ip); 
        
        // ---------  Density  ------------------------------- 
        Matrix2D<float> density(szx,szy);  

        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){

                density(ix,iy)  = F0(0,ix+Nbc,iy+Nbc).real()*pr(0)*pr(0);
                density(ix,iy) += F0(F0.nump()-1,ix+Nbc,iy+Nbc).real()*pr(pr.dim()-1)*pr(pr.dim()-1);
                density(ix,iy) *= 0.5;
                for (size_t ip(1); ip < F0.nump()-1; ++ip){ 
                    density(ix,iy) += pr(ip)*pr(ip)*F0(ip,ix+Nbc,iy+Nbc).real();
                }

                f00  =  F0(1,ix+Nbc,iy+Nbc)/5.0*p0p1_sq;
                f00 += (F0(0,ix+Nbc,iy+Nbc) - F0(1,ix+Nbc,iy+Nbc)*p0p1_sq)*inv_mp0p1_sq
                          *(1.0/3.0-1.0/5.0*p0p1_sq); 
                f00 *= pr(0)*(pr(0)/pr.dx())*pr(0); 

                density(ix,iy) += f00.real();

            }
        }
        density *= 4.0 * M_PI * pr.dx();

        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                if (density(ix,iy) < (4.0*FLT_EPSILON)) density (ix,iy) = 1.0;//4.0*FLT_EPSILON;
             }
        }
        // --------------------------------------------------- 


        // ---------  Vx  ------------------------------------ 
        Matrix2D<float> Vx(szx,szy);  

        //f0 = Y.SH(1,0); using reference instead
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){

                 Vx(ix,iy)  = F10(0,ix+Nbc,iy+Nbc).real() * pow(pr(0),3)*invgamma(0);
                 Vx(ix,iy) += F10(F10.nump()-1,ix+Nbc,iy+Nbc).real() *
                                  pow(pr(pr.dim()-1),3)*invgamma(pr.dim()-1);
                 Vx(ix,iy) *= 0.5;
                 for (size_t ip(1); ip < F10.nump()-1; ++ip) { 
                     Vx(ix,iy)+= pow(pr(ip),3)*invgamma(ip)*F10(ip,ix+Nbc,iy+Nbc).real();
                 }

                 f00  = F10(1,ix+Nbc,iy+Nbc);
                 f00 *= pow(pr(0),5)/(pr.dx()*pr(1)*5.0); 
                 Vx(ix,iy) += f00.real();

            }
        }
        Vx *= (4.0 * M_PI / 3.0) * pr.dx();

        // Divide by the density to get Vx
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                Vx(ix,iy) /= density(ix,iy); 
            }
        }

        // Decide if you want to output Vx
        if (Inputdata::IN().list().o_Vx) {
            string fname("OUTPUT/MOM/Vx/Vx");
            parallel_mom(step, Vx, fname);
        }
        // --------------------------------------------------------

        // ---------  Vy  ------------------------------------ 
        Matrix2D<float> Vy(szx,szy); 

        //f0 = Y.SH(1,1); using reference instead
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){

                Vy(ix,iy)  = F11(0,ix+Nbc,iy+Nbc).real() * pow(pr(0),3)*invgamma(0);
                Vy(ix,iy) += F11(F11.nump()-1,ix+Nbc,iy+Nbc).real() * 
                                 pow(pr(pr.dim()-1),3)*invgamma(pr.dim()-1);
                Vy(ix,iy) *= 0.5;
                for (size_t ip(1); ip < F11.nump()-1; ++ip){ 
                    Vy(ix,iy)+= pow(pr(ip),3)*invgamma(ip)*F11(ip,ix+Nbc,iy+Nbc).real();
                }

                f00  =  F11(1,ix+Nbc,iy+Nbc);
                f00 *= pow(pr(0),5)/(pr.dx()*pr(1)*5.0); 
                Vy(ix,iy) += f00.real();

            }
        }
        Vy *= (8.0 * M_PI / 3.0) * pr.dx();

        // Divide by the density to get Vy
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                Vy(ix,iy) /= density(ix,iy); 
            }
        }

        // Decide if you want to output Vy
        if (Inputdata::IN().list().o_Vy) {
            string fname("OUTPUT/MOM/Vy/Vy");
            parallel_mom(step, Vy, fname);
        }
        // --------------------------------------------------------


        // ---------  VxVx  -------------------------------------- 
        Matrix2D<float> VxVx(szx,szy);
        Matrix2D<float> moment_1(szx,szy); 
        Matrix2D<float> moment_2(szx,szy); moment_2 = 0.0;

        // integral for "f_0^0"
        //f0 = Y.SH(0,0);
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){

                moment_1(ix,iy)  = F0(0,ix+Nbc,iy+Nbc).real() * pow(pr(0),4) *
                                      pow(invgamma(0),2);
                moment_1(ix,iy) += F0(F0.nump()-1,ix+Nbc,iy+Nbc).real() * pow(pr(pr.dim()-1),4)* 
                                      pow(invgamma(pr.dim()-1),2);
                moment_1(ix,iy) *= 0.5;
                for (size_t ip(1); ip < F0.nump()-1; ++ip){ 
                    moment_1(ix,iy)+= pow(pr(ip),4)*pow(invgamma(ip),2)*
                                      F0(ip,ix+Nbc,iy+Nbc).real();
                }

                f00  =  F0(1,ix+Nbc,iy+Nbc)/7.0*p0p1_sq;
                f00 += (F0(0,ix+Nbc,iy+Nbc) - F0(1,ix+Nbc,iy+Nbc)*p0p1_sq)*inv_mp0p1_sq
                      *(1.0/5.0-1.0/7.0*p0p1_sq); 
                f00 *= pow(pr(0),4)*(pr(0)/pr.dx()); 
                moment_1(ix,iy) += f00.real();

            }
        }
        moment_1 *= 4.0 * M_PI / 3.0 * pr.dx();

        // integral for "f20", notice that f20(p=p0) = 0 
        //f0 = Y.SH(2,0);
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){

                moment_2(ix,iy) += F20(F20.nump()-1,ix+Nbc,iy+Nbc).real()* pow(pr(pr.dim()-1),4)*
                                      pow(invgamma(pr.dim()-1),2);
                moment_2(ix,iy) *= 0.5;
                for (size_t ip(1); ip < F20.nump()-1; ++ip){ 
                    moment_2(ix,iy)+= pow(pr(ip),4)*pow(invgamma(ip),2)*
                                       F20(ip,ix+Nbc,iy+Nbc).real();
                }

            }
        }
        moment_2 *= 8.0 * M_PI / 15.0 * pr.dx();

        VxVx  = moment_1;
        VxVx += moment_2;

        // Divide by the density to get VxVx
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                VxVx(ix,iy) /= density(ix,iy); 
            }
        }

        // Output if necessary
        if (Inputdata::IN().list().o_VxVx) {
            string fname("OUTPUT/MOM/VxVx/VxVx");
            parallel_mom(step, VxVx, fname);
        }
        // --------------------------------------------------------- 

        moment_2 *= 0.5;
        moment_1 -= moment_2; // This is the integral 4pi/3f00 - 4pi/15*f20

        // ---------  VyVy  ----------------------------------------  
        Matrix2D<float> VyVy(szx,szy); VyVy  = 0.0;
        //f0 = Y.SH(2,2);
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){

                VyVy(ix,iy) += F20(F20.nump()-1,ix+Nbc,iy+Nbc).real()* pow(pr(pr.dim()-1),4)*
                                   pow(invgamma(pr.dim()-1),2);
                VyVy(ix,iy) *= 0.5;
                for (size_t ip(1); ip < F20.nump()-1; ++ip){ 
                    VyVy(ix,iy)+= pow(pr(ip),4)*pow(invgamma(ip),2)*
                                  F20(ip,ix+Nbc,iy+Nbc).real();
                }
            }
        }
        VyVy *= 16.0 * M_PI / 5.0 * pr.dx();
        moment_2 = VyVy;  // To use for VzVz
        VyVy += moment_1; // This is the integral 16pi/5f22 + 4pi/3f00 - 4pi/15*f20 


        // Divide by the density to get VyVy
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                VyVy(ix,iy) /= density(ix,iy); 
            }
        }

        // Output if necessary
        if (Inputdata::IN().list().o_VyVy) {
            string fname("OUTPUT/MOM/VyVy/VyVy");
            parallel_mom(step, VyVy, fname);
        }
        // --------------------------------------------------------- 

        // ---------  VzVz  ----------------------------------------  
        Matrix2D<float> VzVz(moment_1); 
        VzVz -= moment_2;

        // Divide by the density to get VzVz
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                VzVz(ix,iy) /= density(ix,iy); 
            }
        }

        // Output if necessary
        if (Inputdata::IN().list().o_VzVz) {
            string fname("OUTPUT/MOM/VzVz/VzVz");
            parallel_mom(step, VzVz, fname);
        }
        // -----------------------------------------------------------  


        // ---------  VxVy  ----------------------------------------- 
        Matrix2D<float> VxVy(szx,szy); VxVy = 0.0;
        //f0 = Y.SH(2,1);
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){

                VxVy(ix,iy) += F21(F21.nump()-1,ix+Nbc,iy+Nbc).real()
                               * pow(pr(pr.dim()-1),4)*pow(invgamma(pr.dim()-1),2);
                VxVy(ix,iy) *= 0.5;
                for (size_t ip(1); ip < F21.nump()-1; ++ip){ 
                    VxVy(ix,iy)+= pow(pr(ip),4)*pow(invgamma(ip),2) * 
                                     F21(ip,ix+Nbc,iy+Nbc).real();
                }
            }
        }
        VxVy *= 8.0 * M_PI / 5.0 * pr.dx();

        // Divide by the density to get VxVy
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                VxVy(ix,iy) /= density(ix,iy); 
            }
        }

        // Output if necessary   
        if (Inputdata::IN().list().o_VxVy) {
            string fname("OUTPUT/MOM/VxVy/VxVy");
            parallel_mom(step, VxVy, fname);
        }
        // -------------------------------------------------------- 


        // ---------  VxVz  ---------------------------------------  
        Matrix2D<float> VxVz(szx,szy); VxVz = 0.0;

        // f0 = Y.SH(2,1) just like before
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){

                VxVz(ix,iy) += F21(F21.nump()-1,ix+Nbc,iy+Nbc).imag()
                                * pow(pr(pr.dim()-1),4)*pow(invgamma(pr.dim()-1),2);
                VxVz(ix,iy) *= 0.5;
                for (size_t ip(1); ip < F21.nump()-1; ++ip){ 
                    VxVz(ix,iy)+= pow(pr(ip),4)*invgamma(ip)*invgamma(ip) *
                                     F21(ip,ix+Nbc,iy+Nbc).imag();
                }
            }
        }
        VxVz *= -(8.0 * M_PI / 5.0) * pr.dx();

        // Divide by the density to get VxVz
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                VxVz(ix,iy) /= density(ix,iy); 
            }
        }

        // Output if necessary 
        if (Inputdata::IN().list().o_VxVz) {
            string fname("OUTPUT/MOM/VxVz/VxVz");
            parallel_mom(step, VxVz, fname);
        }
        // --------------------------------------------------------- 


        // ---------  VyVz  ---------------------------------------  
        Matrix2D<float> VyVz(szx,szy); VyVz = 0.0;

        //f0 = Y.SH(2,2);
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){

                 VyVz(ix,iy) += F22(F22.nump()-1,ix+Nbc,iy+Nbc).imag()
                             * pow(pr(pr.dim()-1),4)*pow(invgamma(pr.dim()-1),2);
                 VyVz(ix,iy) *= 0.5;
                 for (size_t ip(1); ip < F22.nump()-1; ++ip){
                     VyVz(ix,iy)+= pow(pr(ip),4)*pow(invgamma(ip),2)*
                                    F22(ip,ix+Nbc,iy+Nbc).imag();
                 }

            }
        }
        VyVz *= -(16.0 * M_PI / 5.0) * pr.dx();

        // Divide by the density to get VxVz
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                VyVz(ix,iy) /= density(ix,iy); 
            }
        }

        // Output if necessary 
        if (Inputdata::IN().list().o_VyVz) {
            string fname("OUTPUT/MOM/VyVz/VyVz");
            parallel_mom(step, VyVz, fname);
        }
        // -------------------------------------------------------- 


        // ---------  v^2 -----------------------------------------  
        Matrix2D<float> Vsq(szx,szy); Vsq = 0.0;
                    
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                Vsq(ix,iy) += VxVx(ix,iy);
                Vsq(ix,iy) += VyVy(ix,iy);
                Vsq(ix,iy) += VzVz(ix,iy);
                Vsq(ix,iy) -= Vx(ix,iy) * Vx(ix,iy);
                Vsq(ix,iy) -= Vy(ix,iy) * Vy(ix,iy);
            }
        }
            
        if (Inputdata::IN().list().o_Vsq) {
            string fname("OUTPUT/MOM/Vsq/Vsq");
            parallel_mom(step, Vsq, fname);
        }
        if (Inputdata::IN().list().o_Temperature) {
//          Convert the Mean Square Velocity to Temperature in eV
            Vsq *= 170880.776; 
            string fname("OUTPUT/MOM/T_eV/T");
            parallel_mom(step, Vsq, fname);
            Vsq *= 1.0/170880.776;
        }
        if (Inputdata::IN().list().o_Pressure) {
            Matrix2D<float> Pre(Vsq); 
//          Convert the Mean Square Velocity to Pressure in Mbar
            Pre *= 2.73785179e-19 * Inputdata::IN().list().density_np;
            for (int iy(0); iy < szy; ++iy){
                for (int ix(0); ix < szx; ++ix){
                    Pre(ix,iy) *= density(ix,iy); 
                }
            }
            string fname("OUTPUT/MOM/P_Mbar/P");
            parallel_mom(step, Pre, fname);
        }
        // --------------------------------------------------------  

//      ND
        Matrix2D<float> ND(szx,szy); 
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                ND(ix,iy) = 1.72e+9*sqrt(pow(Vsq(ix,iy)*170880.776,3) / 
                           (Inputdata::IN().list().density_np * density(ix,iy)));
            }
        }
        if (Inputdata::IN().list().o_ND) {
            string fname("OUTPUT/MOM/ND/ND");
            parallel_mom(step, ND, fname);
        }

//      Collision frequency
        Matrix2D<float> nu(szx,szy); 
        float lnei;
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                lnei = 24.0 - 0.5*log(Inputdata::IN().list().density_np * density(ix,iy))
                                + log(Vsq(ix,iy)*170880.776);
                lnei *= Inputdata::IN().list().Zeta;
                nu(ix,iy) = sqrt(2.0/M_PI)/9.0 * lnei / ND(ix,iy);
            }
        }
        if (Inputdata::IN().list().o_Nu) {
            string fname("OUTPUT/MOM/Nu/nu");
            parallel_mom(step, nu, fname);
        }
        // --------------------------------------------------------  


        // ---------  Qx -----------------------------------------  
        if (Inputdata::IN().list().o_Qx) {
            Matrix2D<float> Qx(szx,szy); Qx = 0.0;
                    
            for (int iy(0); iy < szy; ++iy){
                for (int ix(0); ix < szx; ++ix){
                    // Qx = -<v^2>*Vx
                    Qx(ix,iy) -= Vx(ix,iy)*(VxVx(ix,iy)+VyVy(ix,iy)+VzVz(ix,iy));  
                    // Qx = -<v^2>*Vx - 2(Vxx*Vx+Vxy*Vy) 
                    Qx(ix,iy) -= 2.0 *( VxVx(ix,iy)*Vx(ix,iy) + VxVy(ix,iy)*Vy(ix,iy)); 
                    // Qx = -<v^2>*Vx - 2( Vxx*Vx+Vxy*Vy - Vx^3-VxVy^2 ) 
                    Qx(ix,iy) += 2.0 *( Vx(ix,iy)*Vx(ix,iy) + Vy(ix,iy)*Vy(ix,iy) ) * Vx(ix,iy); 
                    // Qx = n *[ -<v^2>*Vx - 2( Vxx*Vx+Vxy*Vy - Vx^3-VxVy^2) ]
                    Qx(ix,iy) *= density(ix,iy);
                }
            }

            //f0 = Y.SH(1,0);
            // reuse temporary matrix moment_2
            moment_2 = 0.0;
            for (int iy(0); iy < szy; ++iy){
                for (int ix(0); ix < szx; ++ix){

                    moment_2(ix,iy)  = F10(0,ix+Nbc,iy+Nbc).real() * pow(pr(0),5)*pow(invgamma(0),3);
                    moment_2(ix,iy) += F10(F10.nump()-1,ix+Nbc,iy+Nbc).real() *
                                      pow(pr(pr.dim()-1),5)*pow(invgamma(pr.dim()-1),3);
                    moment_2(ix,iy) *= 0.5;
                    for (size_t ip(1); ip < F10.nump()-1; ++ip) { 
                        moment_2(ix,iy)+= pow(pr(ip),5)*pow(invgamma(ip),3)*F10(ip,ix+Nbc,iy+Nbc).real();
                    }

                    f00  = F10(1,ix+Nbc,iy+Nbc);
                    f00 *= pow(pr(0),7)/(pr.dx()*pr(1)*7.0); 
                    moment_2(ix,iy) += f00.real();

                }
            }
            moment_2 *= (4.0 * M_PI / 3.0) * pr.dx();

            Qx += moment_2;
 
            Qx *= 0.5;

            // Output if necessary 
            string fname("OUTPUT/MOM/Qx/Qx");
            parallel_mom(step, Qx, fname);

            // Evaluate the Heat Conductivity
            Matrix2D<float> Conductivity(Vsq);                    // <v^2>
            Conductivity *= Vsq;                                  // (<v^2>)^2
            Conductivity = Conductivity.Dd1();                    // D(<v^2>)^2= 2 * <v^2> * D<v^2>
            Conductivity *= 0.5/(xglob_axis(0)-xglob_axis(1));    // 2 * <v^2> * D<v^2>/Dx  (this needs to be x2-x0 = 2*dx)
            for (int iy(0); iy < szy; ++iy){
                Conductivity(0,iy) = Conductivity(1,iy);          // Put some value for the boundary...
                Conductivity(szx-1,iy) = Conductivity(szx-2,iy);  
            }
            Conductivity *= density;                              // 2 * n * <v^2> * D<v^2>/Dx
            for (int iy(0); iy < szy; ++iy){
                for (int ix(0); ix < szx; ++ix){
                     Conductivity(ix,iy) /= nu(ix,iy);            // 2 * n/(nu) * <v^2> * D<v^2>/Dx
                 }
            }
            Conductivity *= (-1.0)/18.0;                          // -1/9 * n/(nu) * <v^2> * D<v^2>/Dx

            for (int iy(0); iy < szy; ++iy){
                for (int ix(0); ix < szx; ++ix){
                     if ( (Conductivity(ix,iy) <   4.0 *FLT_MIN) &&
                          (Conductivity(ix,iy) > (-4.0)*FLT_MIN) ) {
                           Conductivity (ix,iy) = 1.0;
                     }                                            //- Qx / [1/9 * n/(nu) * <v^2> * D<v^2>/Dx]
                     Conductivity(ix,iy) = Qx(ix,iy) / Conductivity(ix,iy);
                 }
            }

            string fname1("OUTPUT/MOM/Qx/Kx");
            parallel_mom(step, Conductivity, fname1);
                    
        }
        // -------------------------------------------------------- 

        // ---------  Qy -----------------------------------------  
        if (Inputdata::IN().list().o_Qy) {
            Matrix2D<float> Qy(szx,szy); Qy = 0.0;
                    
            for (int iy(0); iy < szy; ++iy){
                for (int ix(0); ix < szx; ++ix){
                    // Qy = -<v^2>*Vy
                    Qy(ix,iy) -= Vy(ix,iy)*(VxVx(ix,iy)+VyVy(ix,iy)+VzVz(ix,iy));  
                    // Qy = -<v^2>*Vy - 2(Vxy*Vx+Vyy*Vy) 
                    Qy(ix,iy) -= 2.0 *( VxVy(ix,iy)*Vx(ix,iy) + VyVy(ix,iy)*Vy(ix,iy)); 
                    // Qy = -<v^2>*Vy - 2( Vxy*Vx+Vyy*Vy - VxVy^2-Vy^3 ) 
                    Qy(ix,iy) += 2.0 *( Vx(ix,iy)*Vx(ix,iy) + Vy(ix,iy)*Vy(ix,iy) ) * Vy(ix,iy); 
                    // Qy = n *[ -<v^2>*Vy - 2( Vxy*Vx+Vyy*Vy - Vx^2Vy-Vy^3) ]
                    Qy(ix,iy) *= density(ix,iy);
                }
            }

            //f0 = Y.SH(1,0);
            // reuse temporary matrix moment_2
            moment_2 = 0.0;
            for (int iy(0); iy < szy; ++iy){
                for (int ix(0); ix < szx; ++ix){

                    moment_2(ix,iy)  = F11(0,ix+Nbc,iy+Nbc).real() * pow(pr(0),5)*pow(invgamma(0),3);
                    moment_2(ix,iy) += F11(F11.nump()-1,ix+Nbc,iy+Nbc).real() *
                                      pow(pr(pr.dim()-1),5)*pow(invgamma(pr.dim()-1),3);
                    moment_2(ix,iy) *= 0.5;
                    for (size_t ip(1); ip < F11.nump()-1; ++ip) { 
                        moment_2(ix,iy)+= pow(pr(ip),5)*pow(invgamma(ip),3)*F11(ip,ix+Nbc,iy+Nbc).real();
                    }

                    f00  = F11(1,ix+Nbc,iy+Nbc);
                    f00 *= pow(pr(0),7)/(pr.dx()*pr(1)*7.0); 
                    moment_2(ix,iy) += f00.real();

                }
            }
            moment_2 *= (8.0 * M_PI / 3.0) * pr.dx();

            Qy += moment_2;
 
            Qy *= 0.5;

            // Output if necessary 
            string fname("OUTPUT/MOM/Qy/Qy");
            parallel_mom(step, Qy, fname);

            // Evaluate the Heat Conductivity
            Matrix2D<float> Conductivity(Vsq);                    // <v^2>
            Conductivity *= Vsq;                                  // (<v^2>)^2
            Conductivity = Conductivity.Dd2();                    // D(<v^2>)^2= 2 * <v^2> * D<v^2>
            Conductivity *= 0.5/(yglob_axis(0)-yglob_axis(1));    // 2 * <v^2> * D<v^2>/Dy  (this needs to be y2-y0 = 2*dy)
            for (int ix(0); ix < szx; ++ix){
                Conductivity(ix,0) = Conductivity(ix,1);          // Put some value for the boundary...
                Conductivity(ix,szy-1) = Conductivity(ix,szy-2);  
            }
            Conductivity *= density;                              // 2 * n * <v^2> * D<v^2>/Dy
            for (int iy(0); iy < szy; ++iy){
                for (int ix(0); ix < szx; ++ix){
                     Conductivity(ix,iy) /= nu(ix,iy);            // 2 * n/(nu) * <v^2> * D<v^2>/Dy
                 }
            }
            Conductivity *= (-1.0)/18.0;                          // -1/9 * n/(nu) * <v^2> * D<v^2>/Dy

            for (int iy(0); iy < szy; ++iy){
                for (int ix(0); ix < szx; ++ix){
                     if ( (Conductivity(ix,iy) <   4.0 *FLT_MIN) &&
                          (Conductivity(ix,iy) > (-4.0)*FLT_MIN) ) {
                           Conductivity (ix,iy) = 1.0;
                     }                                            //- Qx / [1/9 * n/(nu) * <v^2> * D<v^2>/Dy]
                     Conductivity(ix,iy) = Qy(ix,iy) / Conductivity(ix,iy);
                 }
            }

            string fname1("OUTPUT/MOM/Qy/Ky");
            parallel_mom(step, Conductivity, fname1);
                    
        }
        // -------------------------------------------------------- 
              
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Parallel_Output::out_moments(size_t step, Stat& Y) {
//--------------------------------------------------------------
//  Calculation of the output matrices related to the moments
//  of the spherical harmonics
//--------------------------------------------------------------

        SHarmonic f0(Y.SH(0,0));
        Matrix2D<float> moment(szx,szy), moment_1(szx,szy), moment_2(szx,szy); 
        Matrix2D<float> pth(szx,szy);  pth = 0.0;
        complex<double> p0p1_sq(pr(0)*pr(0)/(pr(1)*pr(1))),
                        inv_mp0p1_sq(1.0/(1.0-p0p1_sq)),
                        f00;

        Axis< double > gamma(pr), invgamma(pr);
        for (size_t ip(0); ip < pr.dim(); ++ip) gamma(ip)    = sqrt(1.0+pr(ip)*pr(ip)); 
        for (size_t ip(0); ip < pr.dim(); ++ip) invgamma(ip) = 1.0 / gamma(ip); 
        
        // ---------  Density  --------- 
        moment = 0.0;
        string fname("OUTPUT/MOM/N/x1x2");
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){

                moment(ix,iy)  = f0(0,ix+Nbc,iy+Nbc).real()*pr(0)*pr(0);
                moment(ix,iy) += f0(f0.nump()-1,ix+Nbc,iy+Nbc).real()*pr(pr.dim()-1)*pr(pr.dim()-1);
                moment(ix,iy) *= 0.5;
                for (size_t ip (1); ip < f0.nump()-1; ++ip) 
                    moment(ix,iy)+= pr(ip)*pr(ip)*f0(ip,ix+Nbc,iy+Nbc).real();

                f00  =  f0(1,ix+Nbc,iy+Nbc)/5.0*p0p1_sq;
                f00 += (f0(0,ix+Nbc,iy+Nbc) - f0(1,ix+Nbc,iy+Nbc)*p0p1_sq)*inv_mp0p1_sq
                          *(1.0/3.0-1.0/5.0*p0p1_sq); 
                f00 *= pr(0)*(pr(0)/pr.dx())*pr(0); 

                moment(ix,iy) += f00.real();

             }
         }
         moment *= 4.0 * M_PI * pr.dx();
         if (Inputdata::IN().list().o_x1x2) parallel_mom(step, moment, fname); 
        // ---------------------------- 
              

        // ---------  pth  --------- 
         if (Inputdata::IN().list().o_pth) {
             fname = "OUTPUT/MOM/Pth/pth";
             for (int iy(0); iy < szy; ++iy){
                 for (int ix(0); ix < szx; ++ix){

                     pth(ix,iy)  = f0(0,ix+Nbc,iy+Nbc).real()* pow(pr(0),4);
                     pth(ix,iy) += f0(f0.nump()-1,ix+Nbc,iy+Nbc).real()*pow(pr(pr.dim()-1),4);
                     pth(ix,iy) *= 0.5;
                     for (size_t ip(1); ip < f0.nump()-1; ++ip) 
                         pth(ix,iy)+= pow(pr(ip),4)*f0(ip,ix+Nbc,iy+Nbc).real();

                     f00  =  f0(1,ix+Nbc,iy+Nbc)/7.0*p0p1_sq;
                     f00 += (f0(0,ix+Nbc,iy+Nbc) - f0(1,ix+Nbc,iy+Nbc)*p0p1_sq)*inv_mp0p1_sq
                          *(1.0/5.0-1.0/7.0*p0p1_sq); 
                     f00 *= (pr(0)/pr.dx())*pow(pr(0),4); 

                     pth(ix,iy) += f00.real();

                 }
             }
             pth *= 4.0 * M_PI * pr.dx() / 3.0;
             for (int iy(0); iy < szy; ++iy){
                 for (int ix(0); ix < szx; ++ix){
                     pth(ix,iy) = pth(ix,iy)/moment(ix,iy); 
                     pth(ix,iy) = sqrt(pth(ix,iy));
                 }
             }
             parallel_mom(step, pth, fname);
         }
         // ---------------------------- 

         // ---------  Gamma  --------- 
         moment = 0.0;
         fname = "OUTPUT/MOM/Gam/Gam";
         if (Inputdata::IN().list().o_G) {
             for (int iy(0); iy < szy; ++iy){
                 for (int ix(0); ix < szx; ++ix){

                     moment(ix,iy)  = f0(0,ix+Nbc,iy+Nbc).real()*pr(0)*pr(0)*gamma(0);
                     moment(ix,iy) += f0(f0.nump()-1,ix+Nbc,iy+Nbc).real()
                                     *pr(pr.dim()-1)*pr(pr.dim()-1)*gamma(pr.dim()-1);
                     moment(ix,iy) *= 0.5;
                     for (size_t ip(1); ip < f0.nump()-1; ++ip) 
                         moment(ix,iy)+= gamma(ip)*pr(ip)*pr(ip)*f0(ip,ix+Nbc,iy+Nbc).real();

                     f00  =  f0(1,ix+Nbc,iy+Nbc)/5.0*p0p1_sq;
                     f00 += (f0(0,ix+Nbc,iy+Nbc) - f0(1,ix+Nbc,iy+Nbc)*p0p1_sq)*inv_mp0p1_sq
                          *(1.0/3.0-1.0/5.0*p0p1_sq); 
                     f00 *= pr(0)*(pr(0)/pr.dx())*pr(0); 

                     moment(ix,iy) += f00.real();

                  }
              }
              moment *= 4.0 * M_PI * pr.dx();
              parallel_mom(step, moment, fname);
          }
          // ---------------------------- 


         // ---------  Px  --------- 
         moment = 0.0;
         if (Inputdata::IN().list().o_Px) {
             fname = "OUTPUT/MOM/Px/GVx";
             f0 = Y.SH(1,0);
             for (int iy(0); iy < szy; ++iy){
                 for (int ix(0); ix < szx; ++ix){

                     moment(ix,iy)  = f0(0,ix+Nbc,iy+Nbc).real() * pow(pr(0),3);
                     moment(ix,iy) += f0(f0.nump()-1,ix+Nbc,iy+Nbc).real() * pow(pr(pr.dim()-1),3);
                     moment(ix,iy) *= 0.5;
                     for (size_t ip(1); ip < f0.nump()-1; ++ip) 
                         moment(ix,iy)+= pow(pr(ip),3)*f0(ip,ix+Nbc,iy+Nbc).real();

                     f00  =  f0(1,ix+Nbc,iy+Nbc);
                     f00 *= pow(pr(0),5)/(pr.dx()*pr(1)*5.0); 
                     moment(ix,iy) += f00.real();

                  }
              }
              moment *= (4.0 * M_PI / 3.0) * pr.dx();
              parallel_mom(step, moment, fname);
         }
         // ---------------------------- 

         // ---------  Py  --------- 
         moment = 0.0;
         if (Inputdata::IN().list().o_Py) {
             fname = "OUTPUT/MOM/Py/GVy";
             f0 = Y.SH(1,1);
             for (int iy(0); iy < szy; ++iy){
                 for (int ix(0); ix < szx; ++ix){

                     moment(ix,iy)  = f0(0,ix+Nbc,iy+Nbc).real() * pow(pr(0),3);
                     moment(ix,iy) += f0(f0.nump()-1,ix+Nbc,iy+Nbc).real() * pow(pr(pr.dim()-1),3);
                     moment(ix,iy) *= 0.5;
                     for (size_t ip(1); ip < f0.nump()-1; ++ip) 
                         moment(ix,iy)+= pow(pr(ip),3)*f0(ip,ix+Nbc,iy+Nbc).real();

                     f00  =  f0(1,ix+Nbc,iy+Nbc);
                     f00 *= pow(pr(0),5)/(pr.dx()*pr(1)*5.0); 
                     moment(ix,iy) += f00.real();

                  }
              }
              moment *= (8.0 * M_PI / 3.0) * pr.dx();
              parallel_mom(step, moment, fname);
         }
         // ---------------------------- 

         // ---------  Pz  --------- 
         moment = 0.0;
         if (Inputdata::IN().list().o_Pz) {
             fname = "OUTPUT/MOM/Pz/GVz";
             f0 = Y.SH(1,1);
             for (int iy(0); iy < szy; ++iy){
                 for (int ix(0); ix < szx; ++ix){

                     moment(ix,iy)  = f0(0,ix+Nbc,iy+Nbc).imag() * pow(pr(0),3);
                     moment(ix,iy) += f0(f0.nump()-1,ix+Nbc,iy+Nbc).imag() * pow(pr(pr.dim()-1),3);
                     moment(ix,iy) *= 0.5;
                     for (size_t ip(1); ip < f0.nump()-1; ++ip) 
                         moment(ix,iy)+= pow(pr(ip),3)*f0(ip,ix+Nbc,iy+Nbc).imag();

                     f00  =  f0(1,ix+Nbc,iy+Nbc);
                     f00 *= pow(pr(0),5)/(pr.dx()*pr(1)*5.0); 
                     moment(ix,iy) += f00.imag();

                  }
              }
              moment *= -(8.0 * M_PI / 3.0) * pr.dx();
              parallel_mom(step, moment, fname);
         }
         // ---------------------------- 

         // ---------  PxVx  --------- 
         moment   = 0.0;
         moment_1 = 0.0;
         moment_2 = 0.0;
         if (Inputdata::IN().list().o_PxPx) {
             fname = "OUTPUT/MOM/PxPx/PxVx";
             // integral for "f00"
             f0 = Y.SH(0,0);
             for (int iy(0); iy < szy; ++iy){
                 for (int ix(0); ix < szx; ++ix){

                     moment_1(ix,iy)  = f0(0,ix+Nbc,iy+Nbc).real() * pow(pr(0),4) * invgamma(0);
                     moment_1(ix,iy) += f0(f0.nump()-1,ix+Nbc,iy+Nbc).real() 
                                      * pow(pr(pr.dim()-1),4)*invgamma(pr.dim()-1);
                     moment_1(ix,iy) *= 0.5;
                     for (size_t ip(1); ip < f0.nump()-1; ++ip) 
                         moment_1(ix,iy)+= pow(pr(ip),4)*invgamma(ip)*f0(ip,ix+Nbc,iy+Nbc).real();

                     f00  =  f0(1,ix+Nbc,iy+Nbc)/7.0*p0p1_sq;
                     f00 += (f0(0,ix+Nbc,iy+Nbc) - f0(1,ix+Nbc,iy+Nbc)*p0p1_sq)*inv_mp0p1_sq
                          *(1.0/5.0-1.0/7.0*p0p1_sq); 
                     f00 *= pow(pr(0),4)*(pr(0)/pr.dx()); 
                     moment_1(ix,iy) += f00.real();

                 }
             }
             moment_1 *= 4.0 * M_PI / 3.0 * pr.dx();

             // integral for "f20"
             f0 = Y.SH(2,0);
             for (int iy(0); iy < szy; ++iy){
                 for (int ix(0); ix < szx; ++ix){

                     moment_2(ix,iy)  = f0(0,ix+Nbc,iy+Nbc).real() * pow(pr(0),4) * invgamma(0);
                     moment_2(ix,iy) += f0(f0.nump()-1,ix+Nbc,iy+Nbc).real()
                                      * pow(pr(pr.dim()-1),4)*invgamma(pr.dim()-1);
                     moment_2(ix,iy) *= 0.5;
                     for (size_t ip(1); ip < f0.nump()-1; ++ip) 
                         moment_2(ix,iy)+= pow(pr(ip),4)*invgamma(ip)*f0(ip,ix+Nbc,iy+Nbc).real();
                                                                               
                     f00  =  f0(1,ix+Nbc,iy+Nbc);
                     f00 *= pow(pr(0),7)/(pr.dx()*pr(1)*pr(1)*7.0); 
                     moment_2(ix,iy) += f00.real();

                  }
              }
              moment_2 *= 8.0 * M_PI / 15.0 * pr.dx();

              moment =  moment_1;
              moment += moment_2;
              parallel_mom(step, moment, fname);
         }
         // ---------------------------- 

         moment_2 *= 0.5;
         moment_1 -= moment_2; // This is the integral 4/3f00 - 4/15f20
         moment_2  = moment_1; // copy of the above

         // ---------  PyVy  --------- 
         moment   = 0.0;
         if (Inputdata::IN().list().o_PyPy) {
             fname = "OUTPUT/MOM/PyPy/PyVy";
             f0 = Y.SH(2,2);
             for (int iy(0); iy < szy; ++iy){
                 for (int ix(0); ix < szx; ++ix){

                     moment(ix,iy)  = f0(0,ix+Nbc,iy+Nbc).real() * pow(pr(0),4) * invgamma(0);
                     moment(ix,iy) += f0(f0.nump()-1,ix+Nbc,iy+Nbc).real()
                                    * pow(pr(pr.dim()-1),4)*invgamma(pr.dim()-1);
                     moment(ix,iy) *= 0.5;
                     for (size_t ip(1); ip < f0.nump()-1; ++ip) 
                         moment(ix,iy)+= pow(pr(ip),4)*invgamma(ip) * f0(ip,ix+Nbc,iy+Nbc).real();

                     f00  =  f0(1,ix+Nbc,iy+Nbc);
                     f00 *= pow(pr(0),7)/(pr.dx()*pr(1)*pr(1)*7.0); 
                     moment(ix,iy) += f00.real();
                 }
             }
             moment   *= 16.0 * M_PI / 5.0 * pr.dx();
             moment_1 += moment; 
             parallel_mom(step, moment_1, fname);
         }
         // ---------------------------- 

         // ---------  PzVz  --------- 
         if (Inputdata::IN().list().o_PzPz) {
             fname = "OUTPUT/MOM/PzPz/PzVz";
             moment_2 -= moment;
             parallel_mom(step, moment_2, fname);
         }
         // ---------------------------- 

         // ---------  PyVx  --------- 
         moment   = 0.0;
         if (Inputdata::IN().list().o_PxPy) {
             fname = "OUTPUT/MOM/PxPy/PyVx";
             f0 = Y.SH(2,1);
             for (int iy(0); iy < szy; ++iy){
                 for (int ix(0); ix < szx; ++ix){

                     moment(ix,iy)  = f0(0,ix+Nbc,iy+Nbc).real() * pow(pr(0),4) * invgamma(0);
                     moment(ix,iy) += f0(f0.nump()-1,ix+Nbc,iy+Nbc).real()
                                    * pow(pr(pr.dim()-1),4)*invgamma(pr.dim()-1);
                     moment(ix,iy) *= 0.5;
                     for (size_t ip(1); ip < f0.nump()-1; ++ip) 
                         moment(ix,iy)+= pow(pr(ip),4)*invgamma(ip) * f0(ip,ix+Nbc,iy+Nbc).real();

                     f00  = f0(1,ix+Nbc,iy+Nbc);
                     f00 *= pow(pr(0),7)/(pr.dx()*pr(1)*pr(1)*7.0); 
                     moment(ix,iy) += f00.real();
                 }
             }
             moment *= 8.0 * M_PI / 5.0 * pr.dx();
             parallel_mom(step, moment, fname);
         }
         // ---------------------------- 

         // ---------  PzVx  --------- 
         moment   = 0.0;
         if (Inputdata::IN().list().o_PxPz) {
             fname = "OUTPUT/MOM/PxPz/PzVx";
             f0 = Y.SH(2,1);
             for (int iy(0); iy < szy; ++iy){
                 for (int ix(0); ix < szx; ++ix){

                     moment(ix,iy)  = f0(0,ix+Nbc,iy+Nbc).imag() * pow(pr(0),4) * invgamma(0);
                     moment(ix,iy) += f0(f0.nump()-1,ix+Nbc,iy+Nbc).imag()
                                    * pow(pr(pr.dim()-1),4)*invgamma(pr.dim()-1);
                     moment(ix,iy) *= 0.5;
                     for (size_t ip(1); ip < f0.nump()-1; ++ip) 
                         moment(ix,iy)+= pow(pr(ip),4)*invgamma(ip) * f0(ip,ix+Nbc,iy+Nbc).imag();

                     f00  = f0(1,ix+Nbc,iy+Nbc);
                     f00 *= pow(pr(0),7)/(pr.dx()*pr(1)*pr(1)*7.0); 
                     moment(ix,iy) += f00.imag();
                 }
             }
             moment *= -(8.0 * M_PI / 5.0) * pr.dx();
             parallel_mom(step, moment, fname);
         }
         // ---------------------------- 


         // ---------  PzVy  --------- 
         moment   = 0.0;
         if (Inputdata::IN().list().o_PyPz) {
             fname = "OUTPUT/MOM/PyPz/PzVy";
             f0 = Y.SH(2,2);
             for (int iy(0); iy < szy; ++iy){
                 for (int ix(0); ix < szx; ++ix){

                     moment(ix,iy)  = f0(0,ix+Nbc,iy+Nbc).imag() * pow(pr(0),4) * invgamma(0);
                     moment(ix,iy) += f0(f0.nump()-1,ix+Nbc,iy+Nbc).imag()
                                    * pow(pr(pr.dim()-1),4)*invgamma(pr.dim()-1);
                     moment(ix,iy) *= 0.5;
                     for (size_t ip(1); ip < f0.nump()-1; ++ip) 
                         moment(ix,iy)+= pow(pr(ip),4)*invgamma(ip) * f0(ip,ix+Nbc,iy+Nbc).imag();

                     f00  = f0(1,ix+Nbc,iy+Nbc);
                     f00 *= pow(pr(0),7)/(pr.dx()*pr(1)*pr(1)*7.0); 
                     moment(ix,iy) += f00.imag();
                 }
             }
             moment *= -(16.0 * M_PI / 5.0) * pr.dx();
             parallel_mom(step, moment, fname);
         }
         // ---------------------------- 

    }
//--------------------------------------------------------------

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Parallel Output for data that requires conversion from 
// spherical harmonics to Cartesian phasespace
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//--------------------------------------------------------------
     void Parallel_Output::parallel_p1p2p3(size_t step, Matrix3D<float>& p1p2p3){
//--------------------------------------------------------------
//   Parallel output for p1p2p3
//--------------------------------------------------------------
         MPI_Status status;
         int msg_sz(nump1*nump2*nump3);
         float* p1p2p3buf = new float[msg_sz];
         Matrix3D<float> Mp1p2p3(p1p2p3);

         for (int i(0); i < p1p2p3.dim(); ++i) p1p2p3buf[i] = p1p2p3(i);

         if (NODES() > 1) { 
            if (RANK()!=0) {
               MPI_Send(p1p2p3buf, msg_sz, MPI_FLOAT, 0, RANK(), MPI_COMM_WORLD);
            } 
            else {            
                // Data is already there for rank = 0 
                // Take care of rank > 0 
                for (int rr(1); rr < NODES(); ++rr){
                    MPI_Recv(p1p2p3buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    for (int i(0); i < Mp1p2p3.dim(); ++i) Mp1p2p3(i) += p1p2p3buf[i]; 
                }
             }
         }
         // Fill data for Nodes = 0 
         if (RANK() == 0) exportdata.Export_p1p2p3(step, Mp1p2p3);
         delete[] p1p2p3buf;
         
     }
//--------------------------------------------------------------


//--------------------------------------------------------------
     void Parallel_Output::parallel_p1x1x2(size_t step, Matrix3D<float>& p1x1x2){
//--------------------------------------------------------------
//   Parallel output for p1x1x2
//--------------------------------------------------------------
         MPI_Status status;


         size_t st(0), bi(0);
         int msg_sz(nump1*szx*szy);
         float* p1x1x2buf = new float[msg_sz];
         Matrix3D<float> Mp1x1x2glob(nump1, xglob_axis.dim(), yglob_axis.dim());

         for(int i(0); i < p1x1x2.dim(); ++i) p1x1x2buf[i] = p1x1x2(i);

         if (NODES() > 1) { 
            if (RANK()!=0) {
               MPI_Send(p1x1x2buf, msg_sz, MPI_FLOAT, 0, RANK(), MPI_COMM_WORLD);
            } 
            else {            
                // Fill data for rank = 0 
                bi = 0;
                for(GSIter_F it(Mp1x1x2glob.SubMatrix3D(0, nump1, szx, szy)); it != it.end(); ++it) 
                    *it = p1x1x2buf[bi++];
                // Fill data for rank > 0 
                for (int rr(1); rr < NODES(); ++rr){
                    MPI_Recv(p1x1x2buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    st = (rr / NODESX())*szy*xglob_axis.dim()*nump1+(rr % NODESX())*szx*nump1;
                    bi = 0;
                    for(GSIter_F it(Mp1x1x2glob.SubMatrix3D(st, nump1, szx, szy)); it != it.end(); ++it) 
                        *it = p1x1x2buf[bi++]; 
                 }
             }
         }
         // Fill data for Nodes = 0 
         else { 
             bi = 0;
             for(GSIter_F it(Mp1x1x2glob.SubMatrix3D(0, nump1, szx, szy)); it != it.end(); ++it) 
                 *it = p1x1x2buf[bi++];
         }
         if (RANK() == 0) exportdata.Export_p1x1x2(step, Mp1x1x2glob);
         delete[] p1x1x2buf;

     }
//--------------------------------------------------------------


//--------------------------------------------------------------
     void Parallel_Output::parallel_p1x1(size_t step, Matrix2D<float>& p1x1){
//--------------------------------------------------------------
//   Parallel output for p1x1x2
//--------------------------------------------------------------
         MPI_Status status;

         size_t st(0), bi(0);
         int msg_sz(numpx*szx);
         float* p1x1buf = new float[msg_sz];
         Matrix2D<float> Mp1x1glob(numpx, xglob_axis.dim());

         for(int i(0); i < p1x1.dim(); ++i) p1x1buf[i] = p1x1(i);

         if (NODES() > 1) { 
            if (RANK()!=0) {
               MPI_Send(p1x1buf, msg_sz, MPI_FLOAT, 0, RANK(), MPI_COMM_WORLD);
            } 
            else {            
                // Fill data for rank = 0 
                bi = 0;
                for(GSIter_F it(Mp1x1glob.SubMatrix2D(0, numpx, szx)); it != it.end(); ++it) 
                    *it = p1x1buf[bi++];
                // Fill data for rank > 0 
                for (int rr(1); rr < NODES(); ++rr){
                    MPI_Recv(p1x1buf, msg_sz, MPI_FLOAT, rr, rr, MPI_COMM_WORLD, &status);
                    st = (rr % NODESX())*szx*numpx;
                    bi = 0;
                    for(GSIter_F it(Mp1x1glob.SubMatrix2D(st, numpx, szx)); it != it.end(); ++it) 
                        *it = p1x1buf[bi++]; 
                 }
             }
         }
         // Fill data for Nodes = 0 
         else { 
             bi = 0;
             for(GSIter_F it(Mp1x1glob.SubMatrix2D(0, numpx, szx)); it != it.end(); ++it) 
                 *it = p1x1buf[bi++];
         }
         if (RANK() == 0)  {
             exportdata.Export_p1x1(step, Mp1x1glob);
         }

         delete[] p1x1buf;

     }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Parallel_Output::parallel_pdistr(size_t step, Stat& Y) {
//--------------------------------------------------------------
//  Conversion to cartesian momentum space and related calculations
//--------------------------------------------------------------

        float  pmaxf(Inputdata::IN().outp().p1(nump1-1));
        Savedata::Y_x0_p1p2p3 cartes(Inputdata::IN().inp().l0, Inputdata::IN().inp().m0,
                                     nump1, nump2, nump3, pmaxf, pr);
        Matrix3D<float> tmp_convert(nump1,nump2,nump3);

//        Matrix2D<float> x1x2(szx,szy); x1x2 = 0;
        Matrix3D<float> p1p2p3(nump1,nump2,nump3); p1p2p3 = 0;
        Matrix3D<float> p1x1x2(nump1,szx,szy); p1x1x2 = 0;
        Axis<float>  p1(nump1,pmaxf), p2(nump2,pmaxf), p3(nump3,pmaxf);
        
        for (int iy(0); iy < szy; ++iy){
            for (int ix(0); ix < szx; ++ix){
                tmp_convert = cartes.Convert(Y,ix+Nbc,iy+Nbc);

//              Old version of x1x2
//                if (o_x1x2) {
//                    for (int i(0); i < tmp_convert.dim(); ++i) 
//                        x1x2(ix,iy) += tmp_convert(i);
//                }

                // Calculate output for p1p2p3
                if (Inputdata::IN().list().o_p1p2p3) { p1p2p3 += tmp_convert; }

                // Calculate output for p1x1x2
                if (Inputdata::IN().list().o_p1x1x2) {
                    for (int ip(0); ip < nump1; ++ip)
                        for (GSIter_F it(tmp_convert.d1c(ip,ip)); it!=it.end(); ++it) 
                            p1x1x2(ip,ix,iy) += *it;
                }
             }
         }

//         Old version of x1x2
//         if (o_x1x2)   { 
//             x1x2 *= p1.dx()*p2.dx()*p3.dx();
//             parallel_x1x2(step, x1x2);     
//         }

         if (Inputdata::IN().list().o_p1p2p3) { 
             p1p2p3 *=static_cast<float>(xglob_axis.dx()*yglob_axis.dx()); 
             parallel_p1p2p3(step, p1p2p3); 
         }
         if (Inputdata::IN().list().o_p1x1x2) { 
             p1x1x2 *= p2.dx()*p3.dx();
             parallel_p1x1x2(step, p1x1x2); 
         }

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Parallel_Output::parallel_pdistr_1D(size_t step, Stat& Y) {
//--------------------------------------------------------------
//  Conversion to cartesian momentum space and related calculations
//--------------------------------------------------------------

        float  pxmax(Inputdata::IN().outp().px(numpx-1));
        Matrix2D<float> p1x1(numpx,szx); p1x1 = 0.0;
        valarray<float> p1_line(numpx); p1_line = 0.0;

        Savedata::P1x1_1D cartesian_p1x1(Inputdata::IN().inp().l0, 
                                         numpx, pxmax, pr);
        
        for (size_t ix(0); ix < szx; ++ix) {
            p1_line = cartesian_p1x1.p1x1_out(Y,ix+Nbc,Nbc);
            for (size_t ip(0); ip < numpx; ++ip) {
//                cout << p1_line[ip] << "\n";
                p1x1(ip,ix) = p1_line[ip] ;
            }
        }
        // cout << "Hi I am passed here\n\n";
        // Y_x0_p1p2p3 cartes(Inputdata::IN().inp().l0, Inputdata::IN().inp().m0,
        //                             nump1, nump2, nump3, pmaxf, pr);

        // for(size_t i(0); i < numpx; ++i) cout << Inputdata::IN().outp().px(i) << "hi \n";

//        Matrix2D<float> x1x2(szx,szy); x1x2 = 0;
        //Axis<float>  px(numpx,pxmax);
        

        if (Inputdata::IN().list().o_p1x1) { 

            parallel_p1x1(step, p1x1); 
        }

    }
//--------------------------------------------------------------


//**************************************************************


//**************************************************************
//**************************************************************
//   Definition of the Parallel Environment
//**************************************************************
//**************************************************************


//**************************************************************

//--------------------------------------------------------------
    Parallel_Environment:: Parallel_Environment() : 
//--------------------------------------------------------------
//  Constructor, domain decomposition
//--------------------------------------------------------------
        bndX(Inputdata::IN().list().bndX),           // Type of boundary
        bndY(Inputdata::IN().list().bndY),           // Type of boundary
        //
        NnodesX(Inputdata::IN().list().NnodesX),   // Number of nodes in X-direction
        NnodesY(Inputdata::IN().list().NnodesY),   // Number of nodes in Y-direction
        //
        restart_time(Inputdata::IN().list().restart_time){

        // Determination of the rank and size of the run
        MPI_Comm_size(MPI_COMM_WORLD, &Nnodes);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        rankx = rank % NnodesX;
        ranky = rank / NnodesX;
        
        if (error_check()) {MPI_Finalize(); exit(1);}

        using Inputdata::IN; 

        // Determination of the local computational domain (i.e. the x-axis and the y-axis) 
        IN().inp().x(0) = IN().inp().xglob(rankx*(IN().inp().x.dim()-2*IN().list().RKLevel))
                         -IN().list().RKLevel*IN().inp().xglob.dx();
        for (int i = 1; i < IN().inp().x.dim(); ++i) 
            IN().inp().x(i) = IN().inp().x(i-1)+IN().inp().xglob.dx();

        IN().inp().y(0) = IN().inp().yglob(ranky*(IN().inp().y.dim()-2*IN().list().RKLevel))
                         -IN().list().RKLevel*IN().inp().yglob.dx();
        for (int i = 1; i < IN().inp().y.dim(); ++i) 
            IN().inp().y(i) = IN().inp().y(i-1)+IN().inp().yglob.dx();


         // Restart files will be generated when output_step % restart_step == 0 
         if ( ( (IN().list().n_outsteps + 1) >  IN().list().n_restarts ) &&
               (IN().list().n_restarts > 0)  ) { 
              restart_step = IN().list().n_outsteps / IN().list().n_restarts ;
         } 
         else restart_step = -1;
         if ( IN().list().restart_time > IN().list().n_restarts ) {
             restart_time = 0; 
         }
         else restart_time = IN().list().restart_time;

    }
//--------------------------------------------------------------
//  Destructor
//--------------------------------------------------------------
    Parallel_Environment:: ~Parallel_Environment(){ }
//--------------------------------------------------------------

//--------------------------------------------------------------
    bool Parallel_Environment:: error_check() {
//--------------------------------------------------------------
//  Temporary error check
//--------------------------------------------------------------
        if ( Nnodes != (NnodesX * NnodesY) ) { 
            if (rank == 0) std:: cout << "the number of nodes in the input deck is "
                                      << (NnodesX * NnodesY) << ", terminating ... \n"; 
            return true; 
        }
        // Test the number of cells
        if (rank == 0){ 
            if (Inputdata::IN().list().numx < 2*Inputdata::IN().list().RKLevel+1){
                std::cout<<"Not enough cells per processor in the x direction \n";
                return true;
            }
            if (Inputdata::IN().list().numy < 2*Inputdata::IN().list().RKLevel+1){ 
                std::cout<<"Not enough cells per processor in the y direction \n";
                return true;
            }
         } 
         return false;
    }
//--------------------------------------------------------------
 
//--------------------------------------------------------------
//  Basic information
//--------------------------------------------------------------
    int Parallel_Environment:: RANK()  const {return rank;} 
    int Parallel_Environment:: RANKX() const {return rankx;} 
    int Parallel_Environment:: RANKY() const {return ranky;} 
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    int Parallel_Environment:: NODES() const {return Nnodes;} 
    int Parallel_Environment:: NODESX() const {return NnodesX;} 
    int Parallel_Environment:: NODESY() const {return NnodesY;} 
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    int Parallel_Environment:: BNDX()  const {return bndX;} 
    int Parallel_Environment:: BNDY()  const {return bndY;} 
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Parallel_Environment::Neighbor_ImplicitE_Communications(Stat& Y){
//--------------------------------------------------------------
//  Information exchange between neighbors 
//--------------------------------------------------------------
        int moduloX(RANKX()%2); 
        int RNx((RANKX()+1)%NODESX() +RANKY()*NODESX()),         // This is the right neighbor 
            LNx((RANKX()-1+NODESX())%NODESX()+RANKY()*NODESX()); // This is the left  neighbor 

        if (NODESX() > 1) {
            //even nodes 
            if (moduloX==0){
                Bfield_Data.Send_right_X(Y,RNx);                  //   (Send) 0 --> 1               
                if ((RANKX() != 0) || (BNDX()==0)){
                    Bfield_Data.Recv_from_left_X(Y,LNx);          //          1 --> 0 (Receive)   
                    Bfield_Data.Send_left_X(Y,LNx);               //          1 <-- 0 (Send)
                }
                else {
                    if (BNDX()==1) {
                        Bfield_Data.mirror_bound_Xleft(Y);        // Update node "0" in the x direction
                    }
                    else {
                        cout<<"Invalid Boundary.\n";
                    }
               } 
               Bfield_Data.Recv_from_right_X(Y,RNx);               // (Receive) 0 <-- 1                            
            }
            //odd nodes 
            else {
                Bfield_Data.Recv_from_left_X(Y,LNx);               //           0 --> 1 (Receive)               
                if ((RANKX()!=(NODESX()-1)) || (BNDX()==0)){
                     Bfield_Data.Send_right_X(Y,RNx);              //   (Send)  1 --> 0 
                     Bfield_Data.Recv_from_right_X(Y,RNx);         // (Receive) 1 <-- 0
                }
                else {
                    if (BNDX()==1) {
                        Bfield_Data.mirror_bound_Xright(Y);        // Update node "N-1" in the x direction
                    }
                    else {
                        cout<<"Invalid Boundary.\n";
                    }
               } 
                Bfield_Data.Send_left_X(Y,LNx);                    //           0 <-- 1 (Send)              
            } 
         }
         else { Bfield_Data.sameNode_bound_X(Y); }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

         int moduloY(RANKY()%2);
         int RNy((RANK()+NODESX())%NODES()),                  // This is the right neighbor 
             LNy((RANK()-NODESX()+NODES())%NODES());          // This is the left  neighbor 

         if (NODESY() > 1) {
            //even nodes 
            if (moduloY==0){
                Bfield_Data.Send_right_Y(Y,RNy);                     //   (Send) 0 --> 1               
                if ((RANKY() != 0) || (BNDY()==0)){
                    Bfield_Data.Recv_from_left_Y(Y,LNy);             //          1 --> 0 (Receive)                
                    Bfield_Data.Send_left_Y(Y,LNy);                  //          1 <-- 0 (Send)
                }
                else {
                    if (BNDY()==1) {
                        Bfield_Data.mirror_bound_Yleft(Y);          //  Update node "0" in the y direction
                    }
                    else {
                        cout<<"Invalid Boundary.\n";
                    }
                }  
                Bfield_Data.Recv_from_right_Y(Y,RNy);                // (Receive) 0 <-- 1                            
            }
            //odd nodes 
            else {
                Bfield_Data.Recv_from_left_Y(Y,LNy);                 //           0 --> 1 (Receive)               
                if ((RANKY()!=(NODESY()-1)) || (BNDY()==0)){
                     Bfield_Data.Send_right_Y(Y,RNy);                //   (Send)  1 --> 0 
                     Bfield_Data.Recv_from_right_Y(Y,RNy);           // (Receive) 1 <-- 0
                }
                else {
                    if (BNDY()==1) {
                        Bfield_Data.mirror_bound_Yright(Y);          // Update node "N-1" in the x direction
                    }
                    else {
                        cout<<"Invalid Boundary.\n";
                    }
                } 
                Bfield_Data.Send_left_Y(Y,LNy);                      //           0 <-- 1 (Send)              
            }
         }
         else { Bfield_Data.sameNode_bound_Y(Y); }
    }
//--------------------------------------------------------------



//--------------------------------------------------------------
    void Parallel_Environment::Neighbor_Communications(Stat& Y){
//--------------------------------------------------------------
//  Information exchange between neighbors 
//--------------------------------------------------------------
        int moduloX(RANKX()%2); 
        int RNx((RANKX()+1)%NODESX() +RANKY()*NODESX()),         // This is the right neighbor 
            LNx((RANKX()-1+NODESX())%NODESX()+RANKY()*NODESX()); // This is the left  neighbor 

        if (NODESX() > 1) {
            //even nodes 
            if (moduloX==0){
                X_Data.Send_right_X(Y,RNx);                  //   (Send) 0 --> 1               
                if ((RANKX() != 0) || (BNDX()==0)){
                    X_Data.Recv_from_left_X(Y,LNx);          //          1 --> 0 (Receive)   
                    X_Data.Send_left_X(Y,LNx);               //          1 <-- 0 (Send)
                }
                else {
                    if (BNDX()==1) {
                        X_Data.mirror_bound_Xleft(Y);        // Update node "0" in the x direction
                    }
                    else {
                        cout<<"Invalid Boundary.\n";
                    }
               } 
               X_Data.Recv_from_right_X(Y,RNx);               // (Receive) 0 <-- 1                            
            }
            //odd nodes 
            else {
                X_Data.Recv_from_left_X(Y,LNx);               //           0 --> 1 (Receive)               
                if ((RANKX()!=(NODESX()-1)) || (BNDX()==0)){
                     X_Data.Send_right_X(Y,RNx);              //   (Send)  1 --> 0 
                     X_Data.Recv_from_right_X(Y,RNx);         // (Receive) 1 <-- 0
                }
                else {
                    if (BNDX()==1) {
                        X_Data.mirror_bound_Xright(Y);        // Update node "N-1" in the x direction
                    }
                    else {
                        cout<<"Invalid Boundary.\n";
                    }
               } 
                X_Data.Send_left_X(Y,LNx);                    //           0 <-- 1 (Send)              
            } 
         }
         else { X_Data.sameNode_bound_X(Y); }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

         int moduloY(RANKY()%2);
         int RNy((RANK()+NODESX())%NODES()),                  // This is the right neighbor 
             LNy((RANK()-NODESX()+NODES())%NODES());          // This is the left  neighbor 

         if (NODESY() > 1) {
            //even nodes 
            if (moduloY==0){
                X_Data.Send_right_Y(Y,RNy);                     //   (Send) 0 --> 1               
                if ((RANKY() != 0) || (BNDY()==0)){
                    X_Data.Recv_from_left_Y(Y,LNy);             //          1 --> 0 (Receive)                
                    X_Data.Send_left_Y(Y,LNy);                  //          1 <-- 0 (Send)
                }
                else {
                    if (BNDY()==1) {
                        X_Data.mirror_bound_Yleft(Y);          //  Update node "0" in the y direction
                    }
                    else {
                        cout<<"Invalid Boundary.\n";
                    }
                }  
                X_Data.Recv_from_right_Y(Y,RNy);                // (Receive) 0 <-- 1                            
            }
            //odd nodes 
            else {
                X_Data.Recv_from_left_Y(Y,LNy);                 //           0 --> 1 (Receive)               
                if ((RANKY()!=(NODESY()-1)) || (BNDY()==0)){
                     X_Data.Send_right_Y(Y,RNy);                //   (Send)  1 --> 0 
                     X_Data.Recv_from_right_Y(Y,RNy);           // (Receive) 1 <-- 0
                }
                else {
                    if (BNDY()==1) {
                        X_Data.mirror_bound_Yright(Y);          // Update node "N-1" in the x direction
                    }
                    else {
                        cout<<"Invalid Boundary.\n";
                    }
                } 
                X_Data.Send_left_Y(Y,LNy);                      //           0 <-- 1 (Send)              
            }
         }
         else { X_Data.sameNode_bound_Y(Y); }
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Parallel_Environment::Output(size_t step, Stat& Y){
//--------------------------------------------------------------
//  Decide on what parallel output needs to be calculated
//--------------------------------------------------------------

        using Inputdata::IN; 

        if (IN().list().o_Ex) {
            Out_Data.parallel_Ex(step,Y);
        } 
        if (IN().list().o_Ey) {
            Out_Data.parallel_Ey(step,Y);
        }
        if (IN().list().o_Ez) {
            Out_Data.parallel_Ez(step,Y);
        }
        if (IN().list().o_Bx) {
            Out_Data.parallel_Bx(step,Y);
        } 
        if (IN().list().o_By) {
            Out_Data.parallel_By(step,Y);
        }
        if (IN().list().o_Bz) {
            Out_Data.parallel_Bz(step,Y);
        }
        if ( IN().list().o_x1x2 || IN().list().o_pth ||
             IN().list().o_G  ||
             IN().list().o_Px || IN().list().o_PxPx ||
             IN().list().o_Py || IN().list().o_PxPy || IN().list().o_PyPy ||
             IN().list().o_Pz || IN().list().o_PxPz || IN().list().o_PyPz || IN().list().o_PzPz )  {
            Out_Data.out_moments(step,Y);
        }     
        if ( IN().list().o_Qx || IN().list().o_Qy ||
             IN().list().o_Vsq  ||
             IN().list().o_Temperature  ||
             IN().list().o_Pressure  ||
             IN().list().o_ND ||
             IN().list().o_Nu ||
             IN().list().o_Vx || IN().list().o_VxVx ||
             IN().list().o_Vy || IN().list().o_VxVy || IN().list().o_VyVy ||
                                 IN().list().o_VxVz || IN().list().o_VyVz || IN().list().o_VzVz )  {
            Out_Data.nonrelativistic_moments(step,Y);
        }

        if ( /*o_x1x2 ||*/ IN().list().o_p1x1x2 || IN().list().o_p1p2p3 )  {
            Out_Data.parallel_pdistr(step,Y);
        }
        if ( IN().list().o_p1x1)  {
            Out_Data.parallel_pdistr_1D(step,Y);
        }
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Restart Files
//--------------------------------------------------------------
    size_t Parallel_Environment::T_IN() const { 
        if ( (restart_time > 0) && (restart_step > 0) ) {
            return restart_time * restart_step;
        }
        else return 0;
    } 
//--------------------------------------------------------------

//--------------------------------------------------------------
    bool Parallel_Environment::READ_RESTART() const { 
        return (restart_time > 0);
    } 
//--------------------------------------------------------------

//--------------------------------------------------------------
    bool Parallel_Environment::WRITE_RESTART(const size_t step) const  {
        return  ( (restart_step > 0) && 
                  (static_cast<int>(step) % restart_step == 0) );
    } 
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Parallel_Environment::Write_Restart(size_t step, Stat& Y) {
//--------------------------------------------------------------
//  This is preparing for writing the restart file
//--------------------------------------------------------------
        Restart_Facility Re;
        int Restart_number(static_cast<int>(step) / restart_step);

        Re.Write(RANK(), Restart_number, Y);
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Parallel_Environment::Read_Restart(Stat& Y) {
//--------------------------------------------------------------
//  This is preparing for read the restart file
//--------------------------------------------------------------
        Restart_Facility Re;

        Re.Read(RANK(), restart_time, Y);
    }
//--------------------------------------------------------------

//**************************************************************




