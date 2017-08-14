///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Jun 4, 2013
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the 
//   operators:
//
//   1. Spatial Advection
//
//   2. Electric field
// 
//   3. Current 
// 
//   4. Runge-Kutta 1D Functor
// 
///////////////////////////////////////////////////////////

    #ifndef DECL_VLASOVMAXWELL_H
    #define DECL_VLASOVMAXWELL_H



//**************************************************************

//  Spatial advection
    class Spatial_Advection_1D {
        public:
//      Constructors/Destructors
            Spatial_Advection_1D(size_t Nl, 
                             double pmin, double pmax, size_t Np, 
                             double xmin, double xmax, size_t Nx); 
//          Advance
            void operator()(const DistFunc1D& Din, DistFunc1D& Dh);

        private:
            valarray< double >  A1, A2;
            valarray< double >  vr;
    };
//--------------------------------------------------------------

//  Electric field
    class Electric_Field_1D {
        public:
//      Constructors/Destructors
            Electric_Field_1D(size_t Nl, 
                             double pmin, double pmax, size_t Np, 
                             double xmin, double xmax, size_t Nx); 
//          Advance
            void operator()(const DistFunc1D& Din, const Field1D& FEx, DistFunc1D& Dh);

        private:
            void MakeG00(SHarmonic1D& f);
            void MakeGH( SHarmonic1D& f, size_t l);

            SHarmonic1D H, G;

            valarray< double >  A1, A2;
            valarray< double >  pr, invpr, Hp0;
    };
//--------------------------------------------------------------

//  Current
    class Current_1D {
        public:
//      Constructors/Destructors
            Current_1D( double pmin, double pmax, size_t Np, 
                        size_t Nx ); 
//          Advance
            void operator()(const DistFunc1D& Din, Field1D& FExh);

        private:
            Field1D Jx;
            valarray< double >  pr, invg;
    };
//--------------------------------------------------------------

//  Functor to be used in the Runge-Kutta methods 
    class RKFunctor1D : public Algorithms::AbstFunctor<State1D> {
        public:
//          Constructor
            RKFunctor1D(vector<size_t> Nl, 
                        vector<double> pmax, vector<size_t> Np,  
                        double xmin, double xmax, size_t Nx);
            ~RKFunctor1D(){ };

//          Collect all the operators and apply on Yin
            void operator()(const State1D& Yin, State1D& Yslope);

        private:
            vector<Spatial_Advection_1D> SA;
            vector<Electric_Field_1D>    EF;
            vector<Current_1D>           JX;
    };
//--------------------------------------------------------------
//**************************************************************

    #endif
