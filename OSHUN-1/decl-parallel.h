///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Modified:	May 18th 2011
//	Last Modified:	Dec  9th 2011
///////////////////////////////////////////////////////////

//   
//   Contains the declerations for the communications
//   between nodes, boundaries and the parallel output
///////////////////////////////////////////////////////////
//
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

    #ifndef PARALLEL_ENVIRONMENT_H
    #define PARALLEL_ENVIRONMENT_H


    typedef Slice_iter< complex<double> >  SIter_CD;
    typedef GSlice_iter< complex<double> > GSIter_CD;
    typedef GSlice_iter<float> GSIter_F;


//**************************************************************
//--------------------------------------------------------------
        class Node_ImplicitE_Communications {
//--------------------------------------------------------------
//      Declaration of the parallel module
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Node_ImplicitE_Communications(); 
            ~Node_ImplicitE_Communications();
         
//          Boundary conditions
            int BNDX()   const;
            int BNDY()   const;

//          Data exchange in x direction
            void Send_right_X(Stat& Y, int dest);
            void Recv_from_left_X(Stat& Y,int origin);
            void Send_left_X(Stat& Y, int dest);
            void Recv_from_right_X(Stat& Y, int origin); 

//          Data exchange in y direction
            void Send_right_Y(Stat& Y, int dest);
            void Recv_from_left_Y(Stat& Y,int origin);
            void Send_left_Y(Stat& Y, int dest);
            void Recv_from_right_Y(Stat& Y, int origin);

//          Boundaries 
            void mirror_bound_Xleft(Stat& Y);
            void mirror_bound_Xright(Stat& Y);
            void mirror_bound_Yleft(Stat& Y);
            void mirror_bound_Yright(Stat& Y);

//          Boundaries for single-node configurations
            void sameNode_bound_X(Stat& Y);
            void sameNode_bound_Y(Stat& Y);

        private:
//          Domain information
            int Nbc, bndX, bndY;

//          Information exchange
            int  msg_sizeX, msg_sizeY; 
            complex<double> *msg_bufX, *msg_bufY;

//          Boundaries for single-node configurations
            void sameNode_periodic_X(Stat& Y);
            void sameNode_periodic_Y(Stat& Y);
            void sameNode_mirror_X(Stat& Y);
            void sameNode_mirror_Y(Stat& Y);

        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Node_Communications {
//--------------------------------------------------------------
//      Declaration of the parallel module
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Node_Communications(); 
            ~Node_Communications();
         
//          Boundary conditions
            int BNDX()   const;
            int BNDY()   const;

//          Data exchange in x direction
            void Send_right_X(Stat& Y, int dest);
            void Recv_from_left_X(Stat& Y,int origin);
            void Send_left_X(Stat& Y, int dest);
            void Recv_from_right_X(Stat& Y, int origin); 

//          Data exchange in y direction
            void Send_right_Y(Stat& Y, int dest);
            void Recv_from_left_Y(Stat& Y,int origin);
            void Send_left_Y(Stat& Y, int dest);
            void Recv_from_right_Y(Stat& Y, int origin);

//          Boundaries 
            void mirror_bound_Xleft(Stat& Y);
            void mirror_bound_Xright(Stat& Y);
            void mirror_bound_Yleft(Stat& Y);
            void mirror_bound_Yright(Stat& Y);

//          Boundaries for single-node configurations
            void sameNode_bound_X(Stat& Y);
            void sameNode_bound_Y(Stat& Y);

        private:
//          Domain information
            int Nbc, bndX, bndY;

//          Information exchange
            int  msg_sizeX, msg_sizeY; 
            complex<double> *msg_bufX, *msg_bufY;

//          Boundaries for single-node configurations
            void sameNode_periodic_X(Stat& Y);
            void sameNode_periodic_Y(Stat& Y);
            void sameNode_mirror_X(Stat& Y);
            void sameNode_mirror_Y(Stat& Y);

        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Parallel_Output{
//--------------------------------------------------------------
//      Declaration of the parallel module
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Parallel_Output(); 
            ~Parallel_Output(); 
         
//          Parallel Electromagnetic Fields
            void parallel_Ex(size_t step, Stat& Y);
            void parallel_Ey(size_t step, Stat& Y);
            void parallel_Ez(size_t step, Stat& Y);
            void parallel_Bx(size_t step, Stat& Y);
            void parallel_By(size_t step, Stat& Y);
            void parallel_Bz(size_t step, Stat& Y);

//          Parallel Moments
            void parallel_mom(size_t step, Matrix2D<float>& moment, string& filename);
//          Calculation of up to 2nd order moments
            void out_moments(size_t step, Stat& Y);
            void nonrelativistic_moments(size_t step, Stat& Y);

//          Parallel 3D Output 
            void parallel_p1p2p3(size_t step, Matrix3D<float>& p1p2p3);
            void parallel_p1x1x2(size_t step, Matrix3D<float>& p1x1x2);
            void parallel_p1x1(size_t step, Matrix2D<float>& p1x1);
//          Conversion to cartesian grid
            void parallel_pdistr(size_t step, Stat& Y);
            void parallel_pdistr_1D(size_t step, Stat& Y);

        private:
//          Temporary Export module
            Export_Formatted_Data exportdata;

//          Output matrix sizes 
            size_t numpx, nump1, nump2, nump3;
            size_t szx, szy;

//          Output axis
            Axis<double>  xglob_axis; 
            Axis<double>  yglob_axis; 
            Axis<double>  pr;

//          Parallel parameters
            int RANK()   const;
            int NODES()  const;
            int NODESX() const;

            int Nbc;
            int rank, Nnodes, NnodesX;

        };
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
        class Parallel_Environment {
//--------------------------------------------------------------
//      Declaration of the parallel module
//--------------------------------------------------------------
        public:
//          Constructors/Destructors
            Parallel_Environment(); 
            ~Parallel_Environment(); 
         
//          Parallel parameters
            int RANK()  const;
            int NODES() const;
            int BNDX()  const;
            int BNDY()  const;

//          Restart 
            bool READ_RESTART() const;
            void Read_Restart(Stat& Y); 

            bool WRITE_RESTART(const size_t step) const;
            void Write_Restart(size_t step, Stat& Y); 

            size_t T_IN() const;

//          Output
            void Output(size_t step, Stat& Y);
          
//          Information exchange
            void Neighbor_ImplicitE_Communications(Stat& Y);
            void Neighbor_Communications(Stat& Y);

        private:
//          Parallel parameters
            int RANKX()  const;
            int RANKY()  const;
            int NODESX() const;
            int NODESY() const;
            int rank,   rankx,   ranky;
            int Nnodes, NnodesX, NnodesY;

//          Boundaries 
            int bndX, bndY;

//          Information Exchange
            Node_ImplicitE_Communications Bfield_Data;
            Node_Communications X_Data;

//          Parallel Output  
            Parallel_Output Out_Data;

//          Error Checking of the constructor
            bool error_check();

//          Restart files 
            int restart_time;
            int restart_step;
        };
//--------------------------------------------------------------
//**************************************************************




    #endif
