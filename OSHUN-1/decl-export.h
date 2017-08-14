///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Modified:	May 17 2009
//	Last Modifield:	Nov 29 2011
///////////////////////////////////////////////////////////

//   
//   This header contains the declerations for the temporary
//   export facilities.
///////////////////////////////////////////////////////////
//
// 
//   class Export_Formatted_Data::
//
//   This class receives the output matrices and saves the 
//   data in txt files with recognizable names after it attaches
//   a small header with information necessary for plotting. 
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

    #ifndef DECL_EXPORT_H
    #define DECL_EXPORT_H


//**************************************************************
//**************************************************************
//   Declerations for the Output namespace
//**************************************************************

//*************************************************************

// - - - - - - - - - - - - - - - - - - - - - - - - -
    namespace Export_Facility {  
        bool FileExists(const char* _filePath);
        void Export_numerics(int rank, int list_size, vector<float>& Ehist);
        string Numerics_adj_name(string filetype, size_t NodeRank);
    }
// - - - - - - - - - - - - - - - - - - - - - - - - -
//**************************************************************

//--------------------------------------------------------------
    class Export_Formatted_Data{

    public:
        Export_Formatted_Data();
        
//      - - - - - - - - - - - - - - - - - - - - - - - - -
        void Export_Ex(size_t step, Matrix2D<float>& ex);
        void Export_Ey(size_t step, Matrix2D<float>& ey);
        void Export_Ez(size_t step, Matrix2D<float>& ez);
        void Export_Bx(size_t step, Matrix2D<float>& bx);
        void Export_By(size_t step, Matrix2D<float>& by);
        void Export_Bz(size_t step, Matrix2D<float>& bz);
//      - - - - - - - - - - - - - - - - - - - - - - - - -
        void Export_mom(size_t step, Matrix2D<float>& moment, string moment_name);
//      - - - - - - - - - - - - - - - - - - - - - - - - -
        void Export_p1x1(size_t step, Matrix2D<float>& p1x1);
        void Export_p1p2p3(size_t step, Matrix3D<float>& p1p2p3);
        void Export_p1x1x2(size_t step, Matrix3D<float>& p1x1x2);

    private:
        string adj_name(string filetype, size_t step);

        Axis<double> xglob_axis, yglob_axis, pr; 
        Axis<float>  p1, p2, p3, px;
        string x_str, y_str;

        double dtout;

//      de-normalization coefficients
        double B_MG, E_TVpM, x_um, t_fs;
    }; 
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    class Restart_Facility {

    public:
        Restart_Facility();
        
//      - - - - - - - - - - - - - - - - - - - - - - - - -
        void Read(int rank, int re_step, Stat& Y);
        void Write(int rank, int re_step, Stat& Y);
//      - - - - - - - - - - - - - - - - - - - - - - - - -

    private:
        string adj_name(string filetype, int rank, int re_step);
    }; 
//--------------------------------------------------------------
//**************************************************************

    #endif
