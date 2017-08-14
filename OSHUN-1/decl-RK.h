///////////////////////////////////////////////////////////
//   Contributing authors : Michail Tzoufras
//
//   Last Modified: June 2, 2011
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the 
//   Runge-Kutta methods: RK4, RK3, RK2. 
//
//   These methods are derived from an abstract interface
//   called "Explicit_Method". The idea is that standard 
//   Runge-Kutta can be used in explicit solvers. Nevertheless,
//   an implicit solver may employ parts of this algorithm.
//
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

    #ifndef DECLERATION_RUNG_KUT_H
    #define DECLERATION_RUNG_KUT_H


//--------------------------------------------------------------
        class Explicit_Method {
//--------------------------------------------------------------
//      Abstract interface class
//--------------------------------------------------------------
        public:
            virtual Explicit_Method& advance_p1() = 0; // "covariant" return
            virtual Explicit_Method& advance_p2() = 0; // "covariant" return
            virtual Explicit_Method& advance_IEx() = 0; // "covariant" return
            virtual Explicit_Method& advance_IEy() = 0; // "covariant" return
            virtual Explicit_Method& advance_IEz() = 0; // "covariant" return
            virtual double& tout() = 0;
            virtual double  time() const = 0;
            virtual size_t  numh() const = 0;
            virtual ~Explicit_Method() = 0;
        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//   Decleration for the standard RK4 class
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class Runge_Kutta_4 : public Explicit_Method {
//--------------------------------------------------------------
//      Decleration of the 4rth order Runge-Kutta Class
//--------------------------------------------------------------
        public:
//          Constructor
            Runge_Kutta_4(Stat& Yin, int tout_start);

//          Main function
            Runge_Kutta_4& advance_p1();    
            Runge_Kutta_4& advance_p2();    
            Runge_Kutta_4& advance_IEx();    
            Runge_Kutta_4& advance_IEy();    
            Runge_Kutta_4& advance_IEz();    
         
//          Access
            double& tout();
            double  time() const;
            size_t  numh() const;

        private:
//          Helper functions
            Stat& F(Stat& Yslope);

//          Variables
            Stat& Y;
            Stat  Y0, Y1, Yh;

            Actions Act;
            size_t num_h;
            double h, t, Tout;
        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//   Decleration for RK3 (Osher&Shu Method) class
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class Runge_Kutta_3 : public Explicit_Method {
//--------------------------------------------------------------
//      Decleration of the 3rd order Runge-Kutta Class
//--------------------------------------------------------------
        public:
//          Constructor
            Runge_Kutta_3(Stat& Yin, int tout_start);

//          Main function
            Runge_Kutta_3& advance_p1();    
            Runge_Kutta_3& advance_p2();    
            Runge_Kutta_3& advance_IEx();    
            Runge_Kutta_3& advance_IEy();    
            Runge_Kutta_3& advance_IEz();    
         
//          Access
            double& tout();
            double  time() const;
            size_t  numh() const;

        private:
//          Helper functions
            Stat& F(Stat& Yslope);

//          Variables
            Stat& Y;
            Stat  Y0, Yh;
            Actions Act;
            size_t num_h;
            double h, t, Tout;
        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//   Definition for RK2 (Heun's method) class
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class Runge_Kutta_2 : public Explicit_Method {
//--------------------------------------------------------------
//      Decleration of the 2nd order Runge-Kutta Class
//--------------------------------------------------------------
        public:
//          Constructors
            Runge_Kutta_2(Stat& Yin, int tout_start);

//          Main function
            Runge_Kutta_2& advance_p1();    
            Runge_Kutta_2& advance_p2();   
            Runge_Kutta_2& advance_IEx();    
            Runge_Kutta_2& advance_IEy();    
            Runge_Kutta_2& advance_IEz();    
         
//          Access
            double& tout();
            double  time() const;
            size_t  numh() const;

        private:
//          Helper functions
            Stat& F(Stat& Yslope);

//          Variables
            Stat& Y;
            Stat  Y0, Yh;
            Actions Act;
            size_t num_h;
            double h, t, Tout;
        };
//--------------------------------------------------------------
//**************************************************************

//**************************************************************
//   Definition for RK2 (Heun's method) class for implicit E
//**************************************************************

//**************************************************************
//--------------------------------------------------------------
        class Runge_Kutta_2_Imp : public Explicit_Method {
//--------------------------------------------------------------
//      Decleration of the 2nd order Runge-Kutta Class
//--------------------------------------------------------------
        public:
//          Constructors
            Runge_Kutta_2_Imp(Stat& Yin, int tout_start);

//          Main function
            Runge_Kutta_2_Imp& advance_p1();    
            Runge_Kutta_2_Imp& advance_p2();    
            Runge_Kutta_2_Imp& advance_IEx();    
            Runge_Kutta_2_Imp& advance_IEy();    
            Runge_Kutta_2_Imp& advance_IEz();    
         
//          Access
            double& tout();
            double  time() const;
            size_t  numh() const;

        private:
//          Helper functions
            Stat&    F1(Stat& Yslope);
            Stat&    F2(Stat& Yslope);
            Stat& ImpEx(Stat& Yslope);
            Stat& ImpEy(Stat& Yslope);
            Stat& ImpEz(Stat& Yslope);

//          Variables
            Stat& Y;
            Stat  Y0, Yh;
            Implicit_Actions Act;
            size_t num_h;
            double h, t, Tout;
        };
//--------------------------------------------------------------
//**************************************************************

    #endif
