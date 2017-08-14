///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Modified:	May 18 2011
//	Last Modified:	Nov 29 2011
///////////////////////////////////////////////////////////

//   
//   This cpp file contains the definitions for the functions
//   required to export the data
///////////////////////////////////////////////////////////
//
// 
//   class Export_Formatted_Data::
//
//   This class receives the output matrices and saves the 
//   data in txt files with recognizable name after it attaches
//   a small header with information necessary for plotting. 
//
// 
//   class Restart_Facility::
//
//   This class writes restart files from each node, and 
//   reads restart files for each node.  
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

//  Standard libraries
    #include <iostream>
    #include <vector>
    #include <valarray>
    #include <complex>
    #include <fstream>
    #include <iomanip>
    #include <cstdlib>
    #include <sstream>
    #include <string>
    #include <cstring>
    #include <math.h>

//  My libraries
    #include "matrices.h"

//  Declerations
    #include "decl-input.h"
    #include "decl-state.h"
    #include "decl-export.h"


//**************************************************************
//**************************************************************
//   Definition for the Export class
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Export_Formatted_Data::Export_Formatted_Data()
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------
        : xglob_axis(Inputdata::IN().inp().xglob),   // total x-axis
          yglob_axis(Inputdata::IN().inp().yglob),   // total y-axis
          pr(Inputdata::IN().inp().pr),              // p-radius axis
          dtout(Inputdata::IN().cont().dt_out),
          px(Inputdata::IN().outp().px),
          p1(Inputdata::IN().outp().p1),
          p2(Inputdata::IN().outp().p2),
          p3(Inputdata::IN().outp().p3){
//--------------------------------------------------------------

          B_MG   = -1.0;
          E_TVpM = -1.0;
          x_um   = -1.0;
          t_fs   = -1.0;

          double wp, cwp;

          wp = 5.64*1.0e+4*sqrt(Inputdata::IN().list().density_np);
          cwp = 3.0*1.0e+10/wp;
       
          // Denormalize length and time
          if (Inputdata::IN().list().denormalize_length) {
              x_um = cwp * 1.0e+4;
              t_fs = 1.0e+15 / wp;
              x_str = "x[um]";
              y_str = "y[um]";
          }
          else {
              x_str = "x[c/wp]";
              y_str = "y[c/wp]";
          }
          // Denormalize Electromagnetic fields, E(TV/m), B(MGauss)
          if (Inputdata::IN().list().denormalize_fields) {
              E_TVpM = (0.511*1.0e-6/cwp)*100;
              B_MG   = (1.602*1.0e-19)*(3.0*1.0e+9)/((2.818*1.0e-13)*cwp)*1.0e-6;
          } 

    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Export_Formatted_Data:: Export_Ex(size_t step, Matrix2D<float>& ex){
//--------------------------------------------------------------
//  Export Ex to OUTPUT/FLD/EX/Ex_*.txt
//--------------------------------------------------------------

//      Local variables
        char *filename; 
        const int   accuracy = 10; 
        string      filetype("OUTPUT/FLD/EX/Ex");
        ofstream    SaveFile;

//      Change the filename to indicate de-normalized variables 
        if ( E_TVpM > 0 ) filetype = "OUTPUT/FLD/EX/Ex_d";

//      Executable statements
        filetype = adj_name(filetype,step);
        filename = new char [filetype.size()+1];
        strcpy(filename, filetype.c_str());
        SaveFile.open(filename); 
        SaveFile.flush();

        SaveFile << setprecision(accuracy);

        // header
        if ( E_TVpM > 0 ) {
            SaveFile << "Ex[TV/m]\n";
        }
        else {
            SaveFile << "Ex[mcwp/e]\n";
        }
        SaveFile <<  x_str << "\n";
        SaveFile <<  y_str << "\n";
        SaveFile << static_cast<float>(dtout*step*abs(t_fs))<<"\n"; 
        SaveFile << ex.dim1()<<"\n" << ex.dim2()<<"\n"; 
        for (size_t i(0); i < ex.dim1(); ++i) 
            SaveFile << static_cast<float>(xglob_axis(i)*abs(x_um))<<"\n"; 
        for (size_t j(0); j < ex.dim2(); ++j) 
            SaveFile << static_cast<float>(yglob_axis(j)*abs(x_um))<<"\n"; 
        // data 
        for (size_t i(0); i< ex.dim1()*ex.dim2(); ++i)
                SaveFile << (ex(i)*abs(E_TVpM)) <<"\n";

        SaveFile.flush(); 
        SaveFile.close();
        delete[] filename;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Export_Formatted_Data:: Export_Ey(size_t step, Matrix2D<float>& ey){
//--------------------------------------------------------------
//  Export Ex to OUTPUT/FLD/Ey/Ey_*.txt
//--------------------------------------------------------------

//      Local variables
        char *filename; 
        const int   accuracy = 10; 
        string      filetype("OUTPUT/FLD/EY/Ey");
        ofstream    SaveFile;

//      Change the filename to indicate de-normalized variables 
        if ( E_TVpM > 0 ) filetype = "OUTPUT/FLD/EY/Ey_d";

//      Executable statements
        filetype = adj_name(filetype,step);
        filename = new char [filetype.size()+1];
        strcpy(filename, filetype.c_str());
        SaveFile.open(filename);
        SaveFile.flush();

        SaveFile << setprecision(accuracy);

        // header
        if ( E_TVpM > 0 ) {
            SaveFile << "Ey[TV/m]\n";
        }
        else {
            SaveFile << "Ey[mcwp/e]\n";
        }
        SaveFile <<  x_str << "\n";
        SaveFile <<  y_str << "\n";
        SaveFile << static_cast<float>(dtout*step*abs(t_fs))<<"\n"; 
        SaveFile << ey.dim1()<<"\n" << ey.dim2()<<"\n"; 
        for (size_t i(0); i < ey.dim1(); ++i) 
            SaveFile << static_cast<float>(xglob_axis(i)*abs(x_um))<<"\n"; 
        for (size_t j(0); j < ey.dim2(); ++j) 
            SaveFile << static_cast<float>(yglob_axis(j)*abs(x_um))<<"\n";
        // data
        for (size_t i(0); i< ey.dim1()*ey.dim2(); ++i)
                SaveFile << (ey(i)*abs(E_TVpM)) <<"\n";

        SaveFile.flush();
        SaveFile.close();
        delete[] filename;
    }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Export_Formatted_Data:: Export_Ez(size_t step, Matrix2D<float>& ez){
//--------------------------------------------------------------
//  Export Ez to OUTPUT/FLD/Ez/Ez_*.txt
//--------------------------------------------------------------

//      Local variables
        char *filename; 
        const int   accuracy = 10; 
        string      filetype("OUTPUT/FLD/EZ/Ez");
        ofstream    SaveFile;

//      Change the filename to indicate de-normalized variables 
        if ( E_TVpM > 0 ) filetype = "OUTPUT/FLD/EZ/Ez_d";

//      Executable statements
        filetype = adj_name(filetype,step);
        filename = new char [filetype.size()+1];
        strcpy(filename, filetype.c_str());
        SaveFile.open(filename);
        SaveFile.flush();

        SaveFile << setprecision(accuracy);

        // header
        if ( E_TVpM > 0 ) {
            SaveFile << "Ez[TV/m]\n";
        }
        else {
            SaveFile << "Ez[mcwp/e]\n";
        }
        SaveFile <<  x_str << "\n";
        SaveFile <<  y_str << "\n";
        SaveFile << static_cast<float>(dtout*step*abs(t_fs))<<"\n"; 
        SaveFile << ez.dim1()<<"\n" << ez.dim2()<<"\n"; 
        for (size_t i(0); i < ez.dim1(); ++i) 
            SaveFile << static_cast<float>(xglob_axis(i)*abs(x_um))<<"\n"; 
        for (size_t j(0); j < ez.dim2(); ++j) 
            SaveFile << static_cast<float>(yglob_axis(j)*abs(x_um))<<"\n";
        // data
        for (size_t i(0); i< ez.dim1()*ez.dim2(); ++i)
                SaveFile << (ez(i)*abs(E_TVpM)) <<"\n";

        SaveFile.flush();
        SaveFile.close();
        delete[] filename;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Export_Formatted_Data:: Export_Bx(size_t step, Matrix2D<float>& bx){
//--------------------------------------------------------------
//  Export Bx to OUTPUT/FLD/Bx/Bx_*.txt
//--------------------------------------------------------------

//      Local variables
        char *filename; 
        const int   accuracy = 10; 
        string      filetype("OUTPUT/FLD/BX/Bx");
        ofstream    SaveFile;

//      Change the filename to indicate de-normalized variables 
        if ( B_MG > 0 ) filetype = "OUTPUT/FLD/BX/Bx_d";

//      Executable statements
        filetype = adj_name(filetype,step);
        filename = new char [filetype.size()+1];
        strcpy(filename, filetype.c_str());
        SaveFile.open(filename);
        SaveFile.flush();

        SaveFile << setprecision(accuracy);

        // header
        if ( B_MG > 0 ) {
            SaveFile << "Bx[MG]\n";
        }
        else {
            SaveFile << "Bx[mcwp/e]\n";
        }
        SaveFile <<  x_str << "\n";
        SaveFile <<  y_str << "\n";
        SaveFile << static_cast<float>(dtout*step*abs(t_fs))<<"\n"; 
        SaveFile << bx.dim1()<<"\n" << bx.dim2()<<"\n"; 
        for (size_t i(0); i < bx.dim1(); ++i) 
            SaveFile << static_cast<float>(xglob_axis(i)*abs(x_um))<<"\n"; 
        for (size_t j(0); j < bx.dim2(); ++j) 
            SaveFile << static_cast<float>(yglob_axis(j)*abs(x_um))<<"\n";
        // data
        for (size_t i(0); i< bx.dim1()*bx.dim2(); ++i)
                SaveFile << (bx(i)*abs(B_MG)) <<"\n";

        SaveFile.flush();
        SaveFile.close();
        delete[] filename;

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Export_Formatted_Data:: Export_By(size_t step, Matrix2D<float>& by){
//--------------------------------------------------------------
//  Export By to OUTPUT/FLD/By/By_*.txt
//--------------------------------------------------------------

//      Local variables
        char *filename; 
        const int   accuracy = 10; 
        string      filetype("OUTPUT/FLD/BY/By");
        ofstream    SaveFile;

//      Change the filename to indicate de-normalized variables 
        if ( B_MG > 0 ) filetype = "OUTPUT/FLD/BY/By_d";

//      Executable statements
        filetype = adj_name(filetype,step);
        filename = new char [filetype.size()+1];
        strcpy(filename, filetype.c_str());
        SaveFile.open(filename);
        SaveFile.flush();

        SaveFile << setprecision(accuracy);
 
        // header
        if ( B_MG > 0 ) {
            SaveFile << "By[MG]\n";
        }
        else {
            SaveFile << "By[mcwp/e]\n";
        }
        SaveFile <<  x_str << "\n";
        SaveFile <<  y_str << "\n";
        SaveFile << static_cast<float>(dtout*step*abs(t_fs))<<"\n"; 
        SaveFile << by.dim1()<<"\n" << by.dim2()<<"\n"; 
        for (size_t i(0); i < by.dim1(); ++i) 
            SaveFile << static_cast<float>(xglob_axis(i)*abs(x_um))<<"\n"; 
        for (size_t j(0); j < by.dim2(); ++j) 
            SaveFile << static_cast<float>(yglob_axis(j)*abs(x_um))<<"\n";
        // data
        for (size_t i(0); i< by.dim1()*by.dim2(); ++i)
                SaveFile << (by(i)*abs(B_MG)) <<"\n";

        SaveFile.flush();
        SaveFile.close();
        delete[] filename;

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Export_Formatted_Data:: Export_Bz(size_t step, Matrix2D<float>& bz){
//--------------------------------------------------------------
//  Export Bz to OUTPUT/FLD/Bz/Bz_*.txt
//--------------------------------------------------------------

//      Local variables
        char *filename; 
        const int   accuracy = 10; 
        string      filetype("OUTPUT/FLD/BZ/Bz");
        ofstream    SaveFile;

//      Change the filename to indicate de-normalized variables 
        if ( B_MG > 0 ) filetype = "OUTPUT/FLD/BZ/Bz_d";

//      Executable statements
        filetype = adj_name(filetype,step);
        filename = new char [filetype.size()+1];
        strcpy(filename, filetype.c_str());
        SaveFile.open(filename);
        SaveFile.flush();

        SaveFile << setprecision(accuracy);

        // header
        if ( B_MG > 0 ) {
            SaveFile << "Bz[MG]\n";
        }
        else {
            SaveFile << "Bz[mcwp/e]\n";
        }
        SaveFile <<  x_str << "\n";
        SaveFile <<  y_str << "\n";
        SaveFile << static_cast<float>(dtout*step*abs(t_fs))<<"\n"; 
        SaveFile << bz.dim1()<<"\n" << bz.dim2()<<"\n"; 
        for (size_t i(0); i < bz.dim1(); ++i) 
            SaveFile << static_cast<float>(xglob_axis(i)*abs(x_um))<<"\n"; 
        for (size_t j(0); j < bz.dim2(); ++j) 
            SaveFile << static_cast<float>(yglob_axis(j)*abs(x_um))<<"\n";
        // data
        for (size_t i(0); i< bz.dim1()*bz.dim2(); ++i)
                SaveFile << (bz(i)*abs(B_MG)) <<"\n";


        SaveFile.flush();
        SaveFile.close();
        delete[] filename;

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Export_Formatted_Data:: Export_mom(size_t step, Matrix2D<float>& moment, string moment_name){
//--------------------------------------------------------------
//  Export a moment of the distribution function
//--------------------------------------------------------------

//      Local variables
        char *filename; 
        const int   accuracy = 10; 
        ofstream    SaveFile;

        string      output_title;

//      Executable statements
        if (moment_name == "OUTPUT/MOM/Vx/Vx")     output_title = "Vx[c]";
        if (moment_name == "OUTPUT/MOM/Vy/Vy")     output_title = "Vy[c]";
        if (moment_name == "OUTPUT/MOM/Vz/Vz")     output_title = "Vy[c]";
        if (moment_name == "OUTPUT/MOM/VxVx/VxVx") output_title = "VxVx[c^2]";
        if (moment_name == "OUTPUT/MOM/VyVy/VyVy") output_title = "VyVy[c^2]";
        if (moment_name == "OUTPUT/MOM/VzVz/VzVz") output_title = "VzVz[c^2]";
        if (moment_name == "OUTPUT/MOM/VxVy/VxVy") output_title = "VxVy[c^2]";
        if (moment_name == "OUTPUT/MOM/VxVz/VxVz") output_title = "VxVz[c^2]";
        if (moment_name == "OUTPUT/MOM/VyVz/VyVz") output_title = "VyVz[c^2]";
        if (moment_name == "OUTPUT/MOM/Vsq/Vsq")   output_title = "V^2[c^2]";
        if (moment_name == "OUTPUT/MOM/T_eV/T")    output_title = "T[eV]";
        if (moment_name == "OUTPUT/MOM/P_Mbar/P")  output_title = "P[Mbar]";
        if (moment_name == "OUTPUT/MOM/ND/ND")     output_title = "ND";
        if (moment_name == "OUTPUT/MOM/Qx/Qx")     output_title = "Qx";
        if (moment_name == "OUTPUT/MOM/Qx/Kx")     output_title = "Kx";
        if (moment_name == "OUTPUT/MOM/Qy/Qy")     output_title = "Qy";
        if (moment_name == "OUTPUT/MOM/Qy/Ky")     output_title = "Ky";
        if (moment_name == "OUTPUT/MOM/N/x1x2")    output_title = "N[np]";
        if (moment_name == "OUTPUT/MOM/Pth/pth")   output_title = "pth[mc]";
        if (moment_name == "OUTPUT/MOM/Gam/Gam")   output_title = "Gamma";
        if (moment_name == "OUTPUT/MOM/Px/GVx")    output_title = "Px";
        if (moment_name == "OUTPUT/MOM/Py/GVy")    output_title = "Py";
        if (moment_name == "OUTPUT/MOM/Pz/GVz")    output_title = "Pz";
        if (moment_name == "OUTPUT/MOM/PxPx/PxVx") output_title = "PxVx";
        if (moment_name == "OUTPUT/MOM/PyPy/PyVy") output_title = "PyVy";
        if (moment_name == "OUTPUT/MOM/PzPz/PzVz") output_title = "PzVz";
        if (moment_name == "OUTPUT/MOM/PxPy/PyVx") output_title = "PyVx";
        if (moment_name == "OUTPUT/MOM/PxPz/PzVx") output_title = "PzVx";
        if (moment_name == "OUTPUT/MOM/PyPz/PzVy") output_title = "PzVy";

        moment_name = adj_name(moment_name,step);
        filename = new char [moment_name.size()+1];
        strcpy(filename, moment_name.c_str());
        SaveFile.open(filename);
        SaveFile.flush();

        SaveFile << setprecision(accuracy);

        // header
         
        SaveFile << output_title <<"\n"; 
        SaveFile <<  x_str << "\n";
        SaveFile <<  y_str << "\n";
        SaveFile << static_cast<float>(dtout*step*abs(t_fs))<<"\n"; 
        SaveFile << moment.dim1()<<"\n" << moment.dim2()<<"\n"; 
        for (size_t i(0); i < moment.dim1(); ++i) 
            SaveFile << static_cast<float>(xglob_axis(i)*abs(x_um))<<"\n"; 
        for (size_t j(0); j < moment.dim2(); ++j) 
            SaveFile << static_cast<float>(yglob_axis(j)*abs(x_um))<<"\n";
        // data
        for (size_t i=0; i< moment.dim1()* moment.dim2(); ++i)
                SaveFile << moment(i)<<"\n";

        SaveFile.flush();
        SaveFile.close();
        delete[] filename;

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Export_Formatted_Data::  Export_p1x1(size_t step, Matrix2D<float>& p1x1){
//--------------------------------------------------------------
//  Export high definition p1x1
//--------------------------------------------------------------
//      Local variables
        char *filename; 
        const int   accuracy = 10; 
        string      filetype = "OUTPUT/DISTR/P1X1/p1x1";
        ofstream    SaveFile;

//      Executable statements
        filetype = adj_name(filetype,step);
        filename = new char [filetype.size()+1];
        strcpy(filename, filetype.c_str());
        SaveFile.open(filename);
        SaveFile.flush();

        SaveFile << setprecision(accuracy);

        // header
        SaveFile << "f(p1x1)\n"; 
        SaveFile <<  "p[mc]\n";
        SaveFile <<  x_str << "\n";
        SaveFile << static_cast<float>(dtout*step*abs(t_fs))<<"\n"; 
                
        
        SaveFile << p1x1.dim1()<<"\n" << p1x1.dim2()<<"\n"; 
        for (size_t i(0); i < px.dim(); ++i) 
            SaveFile << static_cast<float>(px(i))<<"\n"; 
        for (size_t i(0); i < p1x1.dim2(); ++i) 
            SaveFile << static_cast<float>(xglob_axis(i)*abs(x_um))<<"\n"; 
        // data
        for (size_t i(0); i< p1x1.dim(); ++i)
            SaveFile << p1x1(i)<<"\n";

        SaveFile.flush();
        SaveFile.close();
        delete[] filename;
   }
//--------------------------------------------------------------


//--------------------------------------------------------------
    void Export_Formatted_Data:: Export_p1p2p3(size_t step, Matrix3D<float>& p1p2p3){
//--------------------------------------------------------------
//  Export 3D data: p1p2p3
//--------------------------------------------------------------

//      Local variables
        char *filename; 
        const int   accuracy = 10; 
        string      filetype = "OUTPUT/DISTR/P1P2P3/p1p2p3";
        ofstream    SaveFile;

//      Executable statements
        filetype = adj_name(filetype,step);
        filename = new char [filetype.size()+1];
        strcpy(filename, filetype.c_str());
        SaveFile.open(filename);
        SaveFile.flush();

        SaveFile << setprecision(accuracy);

        // header
        SaveFile << static_cast<float>(dtout*step)<<"\n"; 
        SaveFile << p1p2p3.dim1()<<"\n" << p1p2p3.dim2()<<"\n" << p1p2p3.dim3()<<"\n"; 
        for (size_t i(0); i < p1.dim(); ++i) 
            SaveFile << static_cast<float>(p1(i))<<"\n"; 
        for (size_t j(0); j < p2.dim(); ++j) 
            SaveFile << static_cast<float>(p2(j))<<"\n";
        for (size_t k(0); k < p3.dim(); ++k) 
            SaveFile << static_cast<float>(p3(k))<<"\n";
        // data
        for (size_t i(0); i< p1p2p3.dim(); ++i)
            SaveFile << p1p2p3(i)<<"\n";

        SaveFile.flush();
        SaveFile.close();
        delete[] filename;

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Export_Formatted_Data:: Export_p1x1x2(size_t step, Matrix3D<float>& p1x1x2){
//--------------------------------------------------------------
//  Export 3D data: p1x1x2
//--------------------------------------------------------------

//      Local variables
        char *filename; 
        const int   accuracy = 10; 
        string      filetype = "OUTPUT/DISTR/P1X1X2/p1x1x2";
        ofstream    SaveFile;

//      Executable statements
        filetype = adj_name(filetype,step);
        filename = new char [filetype.size()+1];
        strcpy(filename, filetype.c_str());
        SaveFile.open(filename);
        SaveFile.flush();

        SaveFile << setprecision(accuracy);

        // header
        SaveFile << static_cast<float>(dtout*step)<<"\n"; 
        SaveFile << p1x1x2.dim1()<<"\n" << p1x1x2.dim2()<<"\n" << p1x1x2.dim3()<<"\n"; 
        for (size_t i(0); i < p1.dim(); ++i) 
            SaveFile << static_cast<float>(p1(i))<<"\n"; 
        for (size_t i(0); i < xglob_axis.dim(); ++i) 
            SaveFile << static_cast<float>(xglob_axis(i))<<"\n"; 
        for (size_t j(0); j < yglob_axis.dim(); ++j) 
            SaveFile << static_cast<float>(yglob_axis(j))<<"\n";
        // data
        for (size_t i=0; i< p1x1x2.dim(); ++i)
            SaveFile << p1x1x2(i)<<"\n";

        SaveFile.flush();
        SaveFile.close();
        delete[] filename;

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    string Export_Formatted_Data:: adj_name(string filetype, size_t step){
//--------------------------------------------------------------
//   Adjust filename
//--------------------------------------------------------------

//      Local Variables
        stringstream adjfilenum(stringstream::in | stringstream::out | stringstream::trunc);
        stringstream filenum(stringstream::in | stringstream::out | stringstream::trunc);
        const int digits = 5;

//      Executalbe statements
        filenum.flush(); adjfilenum.flush();
        filenum << step;
        adjfilenum << "_";
        for (size_t j = (filenum.str()).length(); j < digits; ++j) 
            adjfilenum << "0";
        adjfilenum << step;

        return filetype.append(adjfilenum.str()).append(".txt");
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    bool Export_Facility::FileExists(const char* _filePath) {
//--------------------------------------------------------------
//  Check if a file exists
//--------------------------------------------------------------
        bool retval = false;
        std::ifstream in(_filePath);
        if (in.good())
            retval = true;
        in.close();
        return retval;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Export_Facility::Export_numerics(int rank, int list_size, vector<float>& Ehist){
//--------------------------------------------------------------
//  Export numerics
//--------------------------------------------------------------
// 
//      Local variables
        char *filename; 
        const int   accuracy = 10; 
        string      filetype = "OUTPUT/NUM/En";
        ofstream    SaveFile;

//      Executable statements
        
//      Create filename
        filetype = Numerics_adj_name(filetype,static_cast<size_t>(rank));
        filename = new char [filetype.size()+1];
        strcpy(filename, filetype.c_str());

//      If the file exists don't create the header
        if ( FileExists(filename) ) { 
            SaveFile.open(filename,ios::app);
            SaveFile.flush();
        } 
        else {
//          Make a header
            SaveFile.open(filename,ios::app);
            SaveFile.flush();
            SaveFile << static_cast<float>(list_size)<<"\n"; 
        }

//      Export data
        for(size_t ind(0); ind < Ehist.size(); ++ind ) {
            SaveFile << Ehist.at(ind) << "\n";
        }

        SaveFile.flush();
        SaveFile.close();

        Ehist.erase(Ehist.begin(), Ehist.end());

        delete[] filename;

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    string Export_Facility::Numerics_adj_name(string filetype, size_t NodeRank){
//--------------------------------------------------------------
//   Adjust filename
//--------------------------------------------------------------

//      Local Variables
        stringstream adjfilenum(stringstream::in | stringstream::out | stringstream::trunc);
        stringstream filenum(stringstream::in | stringstream::out | stringstream::trunc);
        const int digits = 5;

//      Executalbe statements
        filenum.flush(); adjfilenum.flush();
        filenum << NodeRank;
        adjfilenum << "_";
        for (size_t j = (filenum.str()).length(); j < digits; ++j) { 
            adjfilenum << "0";
        }
        adjfilenum << NodeRank;

        return filetype.append(adjfilenum.str()).append(".txt");
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//**************************************************************
//   Definition for the Export class
//**************************************************************
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
    Restart_Facility::Restart_Facility() {
//--------------------------------------------------------------
//  Constructor
//--------------------------------------------------------------

    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Restart_Facility::Read(int rank, int re_step, Stat& Y) {
//--------------------------------------------------------------
//  Read restart file
//--------------------------------------------------------------

        char *filename; 
        string      filetype("RESTART/re");

//      Executable statements
        filetype = adj_name(filetype, rank, re_step);
        filename = new char [filetype.size()+1];
        strcpy(filename, filetype.c_str());

        ifstream  fin(filename, ios::binary);

        // Read restart data
        for(size_t nh(0); nh < Y.DF().dim(); ++nh){
            for(size_t i(0); i < Y.SH(nh).dim(); ++i){
                fin.read((char *)(&Y.SH(nh)(i)), sizeof(Y.SH(nh)(i)));
            }
        }
        for(size_t nf(0); nf < Y.EMF().dim(); ++nf){
            for(size_t i(0); i < Y.FLD(nf).numx()*Y.FLD(nf).numy(); ++i){
                fin.read((char *)(&Y.FLD(nf)(i)), sizeof(Y.FLD(nf)(i)));
            }
        }

        fin.close();

        delete[] filename;
         
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    void Restart_Facility::Write(int rank, int re_step, Stat& Y) {
//--------------------------------------------------------------
//  Write restart file
//--------------------------------------------------------------


        char *filename; 
        string      filetype("RESTART/re");

//      Executable statements
        filetype = adj_name(filetype, rank, re_step);
        filename = new char [filetype.size()+1];
        strcpy(filename, filetype.c_str());

        ofstream  fout(filename, ios::binary);

        // Save restart data
        for(size_t nh(0); nh < Y.DF().dim(); ++nh){
            for(size_t i(0); i < Y.SH(nh).dim(); ++i){
                fout.write((char *)(&Y.SH(nh)(i)), sizeof(Y.SH(nh)(i)));
            }
        }
        for(size_t nf(0); nf < Y.EMF().dim(); ++nf){
            for(size_t i(0); i < Y.FLD(nf).numx()*Y.FLD(nf).numy(); ++i){
                fout.write((char *)(&Y.FLD(nf)(i)), sizeof(Y.FLD(nf)(i)));
            }
        }

        fout.flush();
        fout.close();

        delete[] filename;
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
    string Restart_Facility:: adj_name(string filetype, int rank, int re_step){
//--------------------------------------------------------------
//   Adjust filename
//--------------------------------------------------------------

//      Local Variables
        stringstream adjfilenum_rank(stringstream::in | stringstream::out | stringstream::trunc);
        stringstream filenum_rank(stringstream::in | stringstream::out | stringstream::trunc);
        const int digits_rank = 6;

//      Executalbe statements
        filenum_rank.flush(); adjfilenum_rank.flush();
        filenum_rank << rank;
        adjfilenum_rank << "_";
        for (size_t j = (filenum_rank.str()).length(); j < digits_rank; ++j) 
            adjfilenum_rank << "0";
        adjfilenum_rank << rank;

//      Local Variables
        stringstream adjfilenum_step(stringstream::in | stringstream::out | stringstream::trunc);
        stringstream filenum_step(stringstream::in | stringstream::out | stringstream::trunc);
        const int digits_step = 3;

//      Executalbe statements
        filenum_step.flush(); adjfilenum_step.flush();
        filenum_step << re_step;
        adjfilenum_step << "_";
        for (size_t j = (filenum_step.str()).length(); j < digits_step; ++j) 
            adjfilenum_step << "0";
        adjfilenum_step << re_step;

        return filetype.append(adjfilenum_rank.str()).append(adjfilenum_step.str()).append(".dat");
    }
//--------------------------------------------------------------
//**************************************************************


