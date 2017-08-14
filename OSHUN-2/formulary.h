///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Apr 3, 2013
///////////////////////////////////////////////////////////

//   
//   Declerations for data structures that are related to
//   simple physics concepts:
//   --> Axis
//   --> AxisBundle: collection of axis
//   --> units: a concrete class that combines a label and a
//       numerical value
//   --> Formulary: A class that enables conversion of units
//       and implements expression from the NRL plasma
//       formulary. 
///////////////////////////////////////////////////////////

    #ifndef DECL_FORMULARY_H
    #define DECL_FORMULARY_H


//**************************************************************
//--------------------------------------------------------------
// Make Gaussian
template<class T> 
valarray<T> Gaussian(const size_t N, const T vmin, const T vmax, const T vth){

//  Make axis first
    valarray<T> G(N);
    for (size_t i(0); i < N; ++i) {
        G[i] = static_cast<T>(i);
    }
    G *= (vmax-vmin)/(static_cast<T>(N-1));
    G += vmin;

    //for (size_t i(0); i < N; ++i) {
    //    cout << G[i]<<"\n";
    //} exit(1);
//  Make Gaussian
    T           C(pow( 1.0/ (sqrt(2.0*M_PI)*vth), 3));
    T           al( (-0.5) / (vth*vth));
    for (size_t i(0); i < N; ++i) {
        G[i]  = exp( al * G[i]*G[i] ); 
    }
    G *= C;

    return G;                              
}


//--------------------------------------------------------------
class units {
    public:
//      Contents
        string label;   // e.g. label = sec
        double d;  // e.g. x[label] = c * x[1/wp] => c = 1/wp 

//      Constructors
        units() : label("default"), d(0.0) {}
        units(string _x, double _d) : label(_x), d(_d){}

//      Copy constructor 
        units(const units& other) { 
           label = other.label; 
           d = other.d;
        } 
        ~units(){ }
};

class Formulary {
    public:
//    Construct an underlying dictionary for units
      Formulary();

//    Access to unit systems
      units Units(string key) { return D[key]; }
      units Units(string key1, string key2) { return D[key1+"_"+key2]; }
      string Label(string key) { return D[key].label; }
      string Label(string key1, string key2) { return D[key1+"_"+key2].label; }
      double Uconv(string key) { return D[key].d; }
      double Uconv(string key1, string key2) { return D[key1+"_"+key2].d; }

//    Coulomb logarithms
      double LOGee(double ne, double Te);
      double LOGee(double ne, string un, double Te, string uT);
      double LOGei(double ne, double Te, double Z);
      double LOGei(double ne, string un, double Te, string uT, double Z);

//    Formulary functions
      double vth(double Te);
      double vth(double Te, string uT);
      double Tau_e(double ne, double Te);
      double Tau_e(double ne, string un, double Te, string uT);
      double MFP(double ne, double Te);
      double MFP(double ne, string un, double Te, string uT);


//    Normalization density 
      static const double n = 1.0e+21;               // cm-3
      const double wp;

//    Physical constants
      static const double c = 2.99792458*1.0e+10;    // cm/sec
      static const double e = 4.80320425*1.0e-10;    // Franklin or statC 
      static const double m = 9.10938215*1.0e-28;    // gr
       
    private:
      map<string,units> D;

      static const double nmin = 1.0e-8;
};
//--------------------------------------------------------------

Formulary& formulary();
//--------------------------------------------------------------
//**************************************************************


    #endif
