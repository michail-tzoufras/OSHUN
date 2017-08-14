///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	Mar 15, 2013
///////////////////////////////////////////////////////////

//   
//   This header file contains the declerations for the 
//   structures that are required for exporting data:
//
//   1. Xport:
//     Output data
//
//   2. Restart Facility:
//     Restart data
///////////////////////////////////////////////////////////

    #ifndef DECL_EXPORT_H
    #define DECL_EXPORT_H



//**************************************************************
//--------------------------------------------------------------
namespace Export_Files{

    namespace ofconventions {
        const int    ofile_digits = 5;
        const string ofile_extension = ".txt";
        const int    ofile_precision = 6; 

        const int    rfile_digits = 3;
        const int    rank_digits = 6;
        const string rfile_extension = ".dat";
    }
//--------------------------------------------------------------

//  Folder and filename functions
    template<typename T> inline std::string stringify(T const& x) {
        std::ostringstream out;
        out << x;
        return out.str();
    }

    void Makefolder(string _name);
//--------------------------------------------------------------

//  Contains all the info you need to construct an output axis
    class DefaultTags {
        public:
//          Contents
        vector<string> time;
        vector<string> space;
        vector<string> momfld;
        vector<string> pvsx;
        
        DefaultTags(size_t species); 
        ~DefaultTags(){}
    };
//--------------------------------------------------------------

//  redefinitions of "<<" operator for data containers
    template <class T> 
    ofstream& operator<<(ofstream& s, const vector<T>& v) {
        s << setprecision(ofconventions::ofile_precision);
        s << 1 <<"\n";
        s << v.size()<<"\n";
        for (size_t i(0); i < v.size(); ++i) {
            s << v[i]<<"\n";
        }
        return s;
    }

    template <class T> 
    ofstream& operator<<(ofstream& s, const valarray<T>& v) {
        s << setprecision(ofconventions::ofile_precision);
        s << 1 <<"\n";
        s << v.size()<<"\n";
        for (size_t i(0); i < v.size(); ++i) {
            s << v[i]<<"\n";
        }
        return s;
    }

    template <class T> 
    ofstream& operator<<(ofstream& s, const Array2D<T>& array2D) {
        s << setprecision(ofconventions::ofile_precision);
        s << 2 <<"\n";
        s << array2D.dim1()<<"\n";
        s << array2D.dim2()<<"\n";
        for (size_t i(0); i < array2D.dim(); ++i) {
            s << array2D(i)<<"\n";
        }
        return s;
    }

    template <class T> 
    ofstream& operator<<(ofstream& s, const Array3D<T>& array3D) {
        s << setprecision(ofconventions::ofile_precision);
        s << 3 <<"\n";
        s << array3D.dim1()<<"\n";
        s << array3D.dim2()<<"\n";
        s << array3D.dim3()<<"\n";
        for (size_t i(0); i < array3D.dim(); ++i) {
            s << array3D(i)<<"\n";
        }
        return s;
    }
//--------------------------------------------------------------

//  Convert data structure to float structure
    valarray<float> vfloat(const valarray<double>& vDouble); 
    vector<float>   vfloat(const vector<double> vDouble); 


//--------------------------------------------------------------

//  Contains all the info you need to construct an output axis
    class oAxis {
        public:
//          Contents
            string label;           // e.g. label = cm 
            float min, max;  
            size_t sz;

//          Constructors
            oAxis();
            oAxis(const float _m, const float _M, const size_t _sz); 
            oAxis(const string _l, const float _m, const float _M, 
                  const size_t _sz);

//          Copy constructor 
            oAxis(const oAxis& other);
            ~oAxis(){}
    };
//--------------------------------------------------------------

//  Facilitates the generation of a header
    class Header {
        public:
//          Constructor
            Header() { };
            Header(oAxis _x,                                        // 1D
                   string _Ql, float _Qc, string _tl, float _tc, string _oD);
            Header(oAxis _x, oAxis _y,                             // 2D 
                   string _Ql, float _Qc,  string _tl, float _tc, string _oD);
            Header(oAxis _x, oAxis _y, oAxis _z,                  // 3D
                   string _Ql, float _Qc, string _tl, float _tc, string _oD);
            Header(vector< oAxis > _xyz,                            // xD
                   string _Ql, float _Qc, string _tl, float _tc, string _oD);
            size_t dim();    

            valarray<float> axis(const size_t i); // this is axis 0, 1, 2 
            string          label(const size_t i);
            float           conv(const size_t i);
            
            string  Title_label(); 
            float   Title_conv(); 
            string  Time_label();
            float   Time_conv();
            string  Directory();
        
        private:
            vector< oAxis >    xyz_axis;         // axis 0, 1, 2  : size 0-2
            string             title,  time;
            float              titleC, timeC;
            string 	       oDir;
    };
//--------------------------------------------------------------

//  Main facility for exporting data 
    class Xport {
        public:
//          Constructor
            Xport(const Algorithms::AxisBundle<double>& _axis, 
                  const vector< string > oTags,
                  string homedir=""); 

//	    export 1D data structure
            void operator() (const string tag,  valarray<float> data,
                             const size_t step, const size_t spec = 0);  
//	           2D data structure
            void operator() (const string tag,  Array2D<float> data, 
                             const size_t step, const size_t spec = 0); 
//	           3D data structure
            void operator() (const string tag,  Array3D<float> data, 
                             const size_t step, const size_t spec = 0);
        
        private:
            map< string, Header > Hdr; // Dictionary of headers
            string oFextension(size_t species, size_t step);
    };
//--------------------------------------------------------------

    class Restart_Facility {

    public:
        Restart_Facility(string homedir="");
        
        void Read(const int rank, const size_t re_step, State1D& Y);
        void Write(const int rank, const size_t re_step, State1D& Y);

    private:
        string hdir;
        string rFextension(const int rank, const size_t rstep);
    }; 
//--------------------------------------------------------------
}
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
namespace Output_Data{

//--------------------------------------------------------------
//  LEGENDRE VALUES
//  A 2D space is generated from px, |p| axis. The cos8 for
//  this 2D space is calculated and then the Legendre polynomials
//  for each cos8 are calculated. We end up with a 1D array for l
//  containing 2D matrices of the polynomials for each cos8
    class PLegendre1D {
    public:
//      Constructors/Destructors
        PLegendre1D(size_t Nl, size_t Np,  float pmin, float pmax,
                               size_t Npx );
        PLegendre1D(const PLegendre1D& other);
        ~PLegendre1D();

//      Access
        size_t dim()   const { return (*plegendre).size(); }
        Array2D<float>& operator()(size_t i)       { return (*plegendre)[i]; }
        Array2D<float>  operator()(size_t i) const { return (*plegendre)[i]; }   
    
private:
        vector< Array2D<float> > *plegendre;
        size_t lmax;
    };
//--------------------------------------------------------------

//  This class converts the state Y at a specific location
//  "x0" and calculates the integral \int\int f*dpy*dpz 
    class p1x1_1D {
    public:
//      Constructor/Destructor
        p1x1_1D(const Grid_Info& _G);
        p1x1_1D(const p1x1_1D& other);
                        
        ~p1x1_1D();

//      Methods
        valarray<float> operator()(DistFunc1D& df, size_t x0, size_t s) ;

//      Access
        size_t Species()         const { return p1x1.size(); }
        size_t Nl(size_t s)      const  { return Pl[s].dim(); }
        size_t Npx(size_t s)     const  { return polarR[s].dim1(); }
        size_t Np(size_t s)      const  { return polarR[s].dim2(); }
        float  Pmin(size_t s)    const  { return pmin[s]; }
        float  Pmax(size_t s)    const  { return pmax[s]; }

        PLegendre1D      PL(size_t s)     const { return Pl[s]; }
        valarray<float>  p1_x1(size_t s)  const { return p1x1[s]; }
        Array2D<float>   PolarR(size_t s) const { return polarR[s]; }

    private:
        vector< PLegendre1D >    Pl;
        vector< Array2D<float> > polarR;

        vector< valarray<float> > p1x1;
        vector<float> pmin, pmax;
    };
//--------------------------------------------------------------

//  Output Functor
    class Output_Preprocessor_1D {
    public:
//      Constructor
        Output_Preprocessor_1D(const Grid_Info& _grid, 
                               const vector< string > _oTags, 
                               string homedir="")  
        : expo( _grid.axis, _oTags, homedir),
          px_x( _grid ),
          oTags(_oTags) { }

//      Functor
        void operator()(const State1D& Y, const size_t tout);

    private:
        Export_Files::Xport expo;
        p1x1_1D             px_x;
        vector< string >    oTags;
    };

}
//**************************************************************

    #endif
