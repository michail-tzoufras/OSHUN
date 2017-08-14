///////////////////////////////////////////////////////////
//   Contributing authors :	Michail Tzoufras
//
//	Last Modified:	May 17 2013
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
    #include <algorithm>
    #include <fstream>
    #include <iomanip>
    #include <cstdlib>
    #include <sstream>
    #include <string>
    #include <cstring>

    #include <math.h>
    #include <map>

    #include <sys/stat.h>
    #include <sys/types.h>

//  My libraries
    #include "lib-array.h"
    #include "lib-algorithms.h"

//  Declerations
    #include "state.h"
    #include "formulary.h"
    #include "setup.h"
    #include "export.h"


//**************************************************************
//---------------------------------------------------------------
//   Create a folder
void Export_Files::Makefolder(string _name){

     mode_t _permissions(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
     char*   foldername = new char [_name.size()];
     strcpy(foldername, _name.c_str());
     int   status(mkdir(foldername,_permissions));

     if (status != 0) cout << "Warning: Folder "<< _name<< " exists\n";
         
     delete[] foldername;
        
     return;
}
//--------------------------------------------------------------

//--------------------------------------------------------------
// List of the acceptable Tags 
// Constructor 
Export_Files::DefaultTags::DefaultTags(size_t species){ 

//  Time 
    time.push_back( "Time_cgs");
    time.push_back( "Time_si" );
    time.push_back( "Time_fs" );
    time.push_back( "Time_ps" );
    time.push_back( "Time"    ); // Default tag

//  Space
    space.push_back( "Space_cgs");
    space.push_back( "Space_si" );
    space.push_back( "Space");

//  Fields
    momfld.push_back( "Ex"     ); 
    momfld.push_back( "Ex_cgs" );
    momfld.push_back( "Ex_si"  );
    momfld.push_back( "Ey"     );
    momfld.push_back( "Ey_cgs" );
    momfld.push_back( "Ey_si"  );
    momfld.push_back( "Ez"     );
    momfld.push_back( "Ez_cgs" );
    momfld.push_back( "Ez_si"  );
    momfld.push_back( "Bx"     );
    momfld.push_back( "Bx_cgs" );
    momfld.push_back( "Bx_si"  );
    momfld.push_back( "By"     );
    momfld.push_back( "By_cgs" );
    momfld.push_back( "By_si"  );
    momfld.push_back( "Bz"     );
    momfld.push_back( "Bz_cgs" );
    momfld.push_back( "Bz_si"  );
    momfld.push_back( "Jx"     );
    momfld.push_back( "Jx_cgs" );
    momfld.push_back( "Jx_si"  );
    momfld.push_back( "Jy"     );
    momfld.push_back( "Jy_cgs" );
    momfld.push_back( "Jy_si"  );
    momfld.push_back( "Jz"     );
    momfld.push_back( "Jz_cgs" );
    momfld.push_back( "Jz_si"  );

//  Moments
    momfld.push_back( "P"      );
    momfld.push_back( "P_cgs"  );
    momfld.push_back( "P_si"   );  
    momfld.push_back( "P_Mbar" );
    momfld.push_back( "T"      );
    momfld.push_back( "T_cgs"  );
    momfld.push_back( "T_si"   );  
    momfld.push_back( "T_eV" );
    momfld.push_back( "n"      );
    momfld.push_back( "n_cgs"  );
    momfld.push_back( "n_si"   );  

//  p-x
    for (size_t s(0); s < species; ++s) {
        pvsx.push_back( "px-x_"+stringify(s) ); 
        pvsx.push_back( "px-y_"+stringify(s) ); 
        pvsx.push_back( "px-z_"+stringify(s) ); 
        pvsx.push_back( "py-x_"+stringify(s) ); 
        pvsx.push_back( "py-y_"+stringify(s) ); 
        pvsx.push_back( "py-z_"+stringify(s) ); 
        pvsx.push_back( "pz-x_"+stringify(s) ); 
        pvsx.push_back( "pz-y_"+stringify(s) ); 
        pvsx.push_back( "pz-z_"+stringify(s) ); 
     }

}

// Convert data structure to float structure
valarray<float> Export_Files::vfloat(const valarray<double>& vDouble) {
    valarray<float> vf(vDouble.size());
    for (size_t i(0); i < vf.size(); ++i) {
        vf[i] = static_cast<float>(vDouble[i]);
    }
    return vf;
}
vector<float> Export_Files::vfloat(const vector<double> vDouble) {
    vector<float> vf;
    for (size_t i(0); i < vDouble.size(); ++i) {
        vf.push_back(static_cast<float>(vDouble[i]));
    }
    return vf;
}
//--------------------------------------------------------------

// Definition of the output axis
// Constructor 
Export_Files::oAxis::oAxis() : label(""), min(0.0), max(1.0), sz(3) {}
Export_Files::oAxis::oAxis( const float _m, const float _M, 
                            const size_t _sz) 
                : label(""), min(_m), max(_M), sz(_sz) {}
Export_Files::oAxis::oAxis(const string _l, const float _m, const float _M, 
                   const size_t _sz) 
                : label(_l), min(_m), max(_M), sz(_sz) {}
// Copy constructor 
Export_Files::oAxis::oAxis(const oAxis& other) { 
                label = other.label; 
                min   = other.min; 
                max   = other.max; 
                sz    = other.sz; 
} 
//--------------------------------------------------------------

//--------------------------------------------------------------
// 1D header constructor
Export_Files::Header::Header(oAxis _x, 
                             string _Ql, float _Qc, 
                             string _tl, float _tc,
                             string _oD)  
    : title(_Ql), titleC(_Qc), 
      time(_tl),  timeC(_tc),
      oDir(_oD) {
    xyz_axis.push_back(_x);
}

// 2D header constructor
Export_Files::Header::Header(oAxis _x, oAxis _y, 
                             string _Ql, float _Qc, 
                             string _tl, float _tc,
                             string _oD)  
    : title(_Ql),  time(_tl), 
      titleC(_Qc), timeC(_tc),
      oDir(_oD) {
    xyz_axis.push_back(_x);
    xyz_axis.push_back(_y);
}

// 3D header constructor
Export_Files::Header::Header(oAxis _x, oAxis _y, oAxis _z, 
                             string _Ql, float _Qc, 
                             string _tl, float _tc, 
                             string _oD)  
    : title(_Ql), time(_tl), 
      titleC(_Qc), timeC(_tc), 
      oDir(_oD) {
    xyz_axis.push_back(_x);
    xyz_axis.push_back(_y);
    xyz_axis.push_back(_z);
}

// xD header constructor
Export_Files::Header::Header(vector< oAxis > _xyz,
                             string _Ql, float _Qc, 
                             string _tl, float _tc, 
                             string _oD)  
    : xyz_axis(_xyz), title(_Ql), time(_tl), 
      titleC(_Qc), timeC(_tc), 
      oDir(_oD) {}

// number of header dimensions
size_t Export_Files::Header::dim() { 
    return xyz_axis.size(); 
}     
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

valarray<float> Export_Files::Header::axis(const size_t i) { 
    return Algorithms::MakeAxis(xyz_axis[i].min, xyz_axis[i].max, xyz_axis[i].sz); 
}
string Export_Files::Header::label(const size_t i) { 
    return xyz_axis[i].label; 
}
            
string Export_Files::Header::Title_label() { return title; }
float Export_Files::Header::Title_conv()  { return titleC; }
string Export_Files::Header::Time_label()  { return time; }
float Export_Files::Header::Time_conv()   { return timeC; }
string Export_Files::Header::Directory()   { return oDir; }
//--------------------------------------------------------------


//--------------------------------------------------------------
// Constructor of the export facility for data structures
Export_Files::Xport::Xport(const Algorithms::AxisBundle<double>& _axis, 
                           const vector< string > oTags,
                           string homedir){

    size_t species(_axis.pdim());
    Makefolder(homedir + "OUTPUT");
    DefaultTags dTags(species);

    vector< oAxis > xyz, pxyz;
    xyz.push_back( Export_Files::oAxis(_axis.xmin(0), _axis.xmax(0), _axis.Nx(0)) );
    for (size_t s(0); s < species; ++s) {
        pxyz.push_back( Export_Files::oAxis( (-1.0)*float(_axis.pmax(s)), float(_axis.pmax(s)), _axis.Npx(s)) );
    }

//  Time 
    size_t tloc(0);                  // Find the location of the right tag
    while ( ( tloc < dTags.time.size()-1 ) &&
            ( find(oTags.begin(),oTags.end(), dTags.time[tloc]) == oTags.end() ) ) {
        ++tloc;
    } 
    string tlabel = "t[" +formulary().Label(dTags.time[tloc])+"]";
    float  tconv  =  formulary().Uconv(dTags.time[tloc]);

//  xyz Axis
    size_t xloc(0);                  // Find the location of the right tag
    while ( ( xloc < dTags.space.size()-1 ) &&
            ( find(oTags.begin(),oTags.end(), dTags.space[xloc]) == oTags.end() ) ) {
        ++xloc;
    } 
    xyz[0].label = "x["+ formulary().Label(dTags.space[xloc]) +"]";
    if ( xyz.size() > 1 ) xyz[1].label = "y["+ formulary().Label(dTags.space[xloc]) +"]";
    if ( xyz.size() > 2 ) xyz[2].label = "z["+ formulary().Label(dTags.space[xloc]) +"]";

//  pxyz Axis
    for (size_t i(0); i < species; ++i) {
        pxyz[i].label = "px[mc]";
        if ( pxyz.size()/species > 1 ) pxyz[  species+i].label = "py[mc]";
        if ( pxyz.size()/species > 2 ) pxyz[2*species+i].label = "pz[mc]";
    }

//  Tags for Moments and Fields -->
    for (size_t i(0); i < dTags.momfld.size(); ++i) {

 //     If this tag is an output tag
	if ( find(oTags.begin(),oTags.end(), dTags.momfld[i]) != oTags.end() ) { 

            string nounits = dTags.momfld[i].substr(0, dTags.momfld[i].find("_"));
            string folder = homedir + "OUTPUT/" + nounits + "/";
            Makefolder(folder);
 //         Generate a header file for this tag
            Hdr[dTags.momfld[i]] = Header(xyz, 
	        nounits+"["+formulary().Label(dTags.momfld[i])+"]",
	        formulary().Uconv(dTags.momfld[i]),
                tlabel, tconv, folder);    
        }
    } // <--

//  Tags for p-x -->
    for (size_t i(0); i < dTags.pvsx.size(); ++i) {

 //     If this tag is an output tag
	if ( find(oTags.begin(),oTags.end(), dTags.pvsx[i]) != oTags.end() ) { 

            string folder = homedir + "OUTPUT/" + dTags.pvsx[i] + "/";
            Makefolder(folder);
//          Generate a header file for this tag
//          For each 9 you have a different species 
            Hdr[dTags.pvsx[i]] = Header( pxyz[i/9+((i%9)/3)*species], xyz[i%3], 
                "f"+stringify(i/9), 1.0, tlabel, tconv, folder); 
        }
    } //<--
}
//--------------------------------------------------------------


//--------------------------------------------------------------
//export 1D data structure
void Export_Files::Xport::operator() 
          (const string tag,  valarray<float> data, 
           const size_t step, const size_t spec){

    string      filename(Hdr[tag].Directory());
    ofstream    oFile;

//  Check Header file correctness
    if (Hdr[tag].dim() != 1) {
        cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 1D structure\n"; 
        exit(1);
    }
 
//  Open File
    filename.append(tag).append(oFextension(spec,step));
    oFile.open(filename.c_str()); 
    oFile.flush();

//  Export dimensions
    oFile << Hdr[tag].dim() << "\n";

//  Export time
    oFile << Hdr[tag].Time_label() << "\n";
    oFile << float(step)*Hdr[tag].Time_conv() << "\n";  // 

//  Export all axes
    for (size_t i(0); i< Hdr[tag].dim(); ++i) {
        oFile << Hdr[tag].label(i) << "\n"; 
        oFile << Hdr[tag].axis(i);
    }
//  Renormalize and export data
    oFile << Hdr[tag].Title_label() << "\n"; 
    data *= Hdr[tag].Title_conv();         
    oFile << data; 

    oFile.flush();
    oFile.close();
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//export 2D data structure
void Export_Files::Xport::operator() 
          (const string tag,  Array2D<float> data, 
           const size_t step, const size_t spec){

    string      filename(Hdr[tag].Directory());
    ofstream    oFile;
     
//  Check Header file correctness
    if (Hdr[tag].dim() != 2) {
        cout << "ERROR "<< tag <<" : "  << Hdr[tag].dim() << " dimensions != 2D structure\n"; 
        exit(1);
    }

//  Open File
    filename.append(tag).append(oFextension(spec,step));
    oFile.open(filename.c_str()); 
    oFile.flush();

//  Export dimensions
    oFile << Hdr[tag].dim() << "\n";

//  Export time
    oFile << Hdr[tag].Time_label() << "\n";
    oFile << float(step)*Hdr[tag].Time_conv() << "\n";  // 

//  Export all axes
    for (size_t i(0); i< Hdr[tag].dim(); ++i) {
        oFile << Hdr[tag].label(i) << "\n"; 
        oFile << Hdr[tag].axis(i);
    }
//  Renormalize and export data
    oFile << Hdr[tag].Title_label() << "\n"; 
    data *= Hdr[tag].Title_conv();         
    oFile << data; 

    oFile.flush();
    oFile.close();
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

//export 3D data structure
void Export_Files::Xport::operator() 
          (const string tag,  Array3D<float> data, 
           const size_t step, const size_t spec){

    string      filename(Hdr[tag].Directory());
    ofstream    oFile;
     
//  Check Header file correctness
    if (Hdr[tag].dim() != 3) {
        cout << "ERROR "<<tag<<" : "  << Hdr[tag].dim() << " dimensions != 3D structure\n"; 
        exit(1);
    }

//  Open File
    filename.append(tag).append(oFextension(spec,step));
    oFile.open(filename.c_str()); 
    oFile.flush();

//  Export dimensions
    oFile << Hdr[tag].dim() << "\n";

//  Export time
    oFile << Hdr[tag].Time_label() << "\n";
    oFile << float(step)*Hdr[tag].Time_conv() << "\n";  // 

//  Export all axes
    for (size_t i(0); i< Hdr[tag].dim(); ++i) {
        oFile << Hdr[tag].label(i) << "\n"; 
        oFile << Hdr[tag].axis(i);
    }
//  Renormalize and export data
    oFile << Hdr[tag].Title_label() << "\n"; 
    data *= Hdr[tag].Title_conv();         
    oFile << data; 

    oFile.flush();
    oFile.close();
}
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Adjust filenames with zeros to reach some prescribed length. 
//  Add the filename extension.  
string Export_Files::Xport::oFextension(size_t species, size_t step){
     stringstream sFilename;
     sFilename << "_s" << species <<"_";
        
     // Number of zeros to add to the filename
     int Nzeros(ofconventions::ofile_digits - stringify(step).length()); 
     while (Nzeros-- > 0) { 
         sFilename << "0"; 
     }
        
     sFilename << step << ofconventions::ofile_extension;

     return sFilename.str();
}
//--------------------------------------------------------------
//**************************************************************



//**************************************************************
//--------------------------------------------------------------
    Export_Files::Restart_Facility::Restart_Facility(string homedir) {
        hdir = homedir;
        Makefolder(hdir+"RESTART/");
//        if (makefolder(hdir+"RESTART/") != 0) cout<<"Warning: Folder "<< hdir+"RESTART/"<<" exists\n";
    }

//--------------------------------------------------------------
//  Read restart file
    void Export_Files::Restart_Facility::Read(const int rank, const size_t re_step, State1D& Y) {

//      Generate filename 
        string   filename(hdir+"RESTART/re_1D_");
        filename.append(rFextension(rank,re_step));

//      Open file
        ifstream  fin(filename.c_str(), ios::binary);

//      Read distribution functions
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t nh(0); nh < Y.DF(s).dim(); ++nh) {
                for(size_t i(0); i < Y.SH(s,nh).dim(); ++i) {
                    fin.read((char *)(&Y.SH(s,nh)(i)), sizeof(Y.SH(s,nh)(i)));
                }
            }
        }

//      Read Ex
        for(size_t i(0); i < Y.Ex().numx(); ++i){
            fin.read((char *)(&Y.Ex()(i)), sizeof(Y.Ex()(i)));
        }

        fin.close();
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Write restart file
    void Export_Files::Restart_Facility::Write(const int rank, const size_t re_step, State1D& Y) {

//      Generate filename 
        string   filename(hdir+"RESTART/re_1D_");
        filename.append(rFextension(rank,re_step));

//      Open file
        ofstream  fout(filename.c_str(), ios::binary);

//      Write distribution functions
        for(size_t s(0); s < Y.Species(); ++s) {
            for(size_t nh(0); nh < Y.DF(s).dim(); ++nh) {
                for(size_t i(0); i < Y.SH(s,nh).dim(); ++i) {
                    fout.write((char *)(&Y.SH(s,nh)(i)), sizeof(Y.SH(s,nh)(i)));
                }
            }
        }

//      Write Ex
        for(size_t i(0); i < Y.Ex().numx(); ++i) {
            fout.write((char *)(&Y.Ex()(i)), sizeof(Y.Ex()(i)));
        }

        fout.flush();
        fout.close();
    }
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Adjust filenames with zeros to reach some prescribed length. 
//  Add the filename extension. 
string Export_Files::Restart_Facility::rFextension(const int rank, const size_t rstep){
        stringstream sFilename;

        // Number of zeros to add to the filename
        int Nzeros(ofconventions::rank_digits - stringify(rank).length()); 
        while (Nzeros-- > 0) { 
            sFilename << "0"; 
        } 
        
        sFilename << rank << "_";

        // Number of zeros to add to the filename
        Nzeros = ofconventions::rfile_digits - stringify(rstep).length(); 
        while (Nzeros-- > 0) { 
            sFilename << "0";
        } 
        
        sFilename << rstep << ofconventions::rfile_extension;

        return sFilename.str();
    }
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Output_Data::PLegendre1D::PLegendre1D( size_t Nl, size_t Np,  
              float pmin, float pmax, size_t Npx ) {

    valarray<float> p( Algorithms::MakeAxis(pmin, pmax, Np)),
                    px(Algorithms::MakeAxis(float(-1.0)*pmax, pmax, Npx));
     
//  Generate the structure to save the polynomials
    plegendre = new vector< Array2D<float> > ; 
    for (size_t l(0); l < Nl; ++l) { 
        (*plegendre).push_back( Array2D<float>(Npx,Np) );
    }

//  Generate the polynomial values for each cos(theta) = px/p 
    for (size_t j(0); j < p.size(); ++j) {
        float invp(1.0/p[j]);
        for (size_t i(0); i < px.size(); ++i) {

//          For given px/p generate all the polynomial values ...
            valarray<float> vL( Algorithms::Legendre( px[i]*invp, Nl ) );
//          ... and save them
            for (size_t l(0); l < Nl; ++l) {
                (*plegendre)[l](i,j) = vL[l];
            }
        }      
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::PLegendre1D::PLegendre1D( const PLegendre1D& other ) {

//  Generate the structure to save the polynomials
    plegendre = new vector< Array2D<float> > ; 
    for (size_t l(0); l < other.dim(); ++l) { 
        (*plegendre).push_back( other(l) );
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::PLegendre1D::~PLegendre1D(){
    delete plegendre;
}
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------
Output_Data::p1x1_1D::p1x1_1D( const Grid_Info& _G) {

//  Generate the required structures
    for (size_t s(0); s < _G.axis.pdim(); ++s) {
        pmax.push_back( static_cast<float>(_G.axis.pmax(s)) );
        pmin.push_back( static_cast<float>(_G.axis.pmin(s)) ); 
        Pl.push_back( PLegendre1D( _G.Nl[s], _G.axis.Np(s), static_cast<float>(_G.axis.pmin(s)), 
                                             static_cast<float>(_G.axis.pmax(s)), _G.axis.Npx(s) ) );
        polarR.push_back( Array2D<float>( _G.axis.Npx(s), _G.axis.Np(s) ) );
        p1x1.push_back( valarray<float>( _G.axis.Npx(s)) );
    }

//  Calculate a 2D array for the polar radius for each px and pr   
    for (size_t s(0); s < _G.axis.pdim(); ++s) {
        valarray<float> p( Algorithms::MakeAxis( pmin[s], pmax[s], _G.axis.Np(s) )),
                        px(Algorithms::MakeAxis( float(-1.0)* pmax[s], pmax[s], _G.axis.Npx(s) ));
//      --->
        for (size_t i(0); i < _G.axis.Npx(s); ++i){
            for (size_t j(0); j < _G.axis.Np(s); ++j){
                float polarR_sq  = p[j] * p[j] - px[i] * px[i];
                if (polarR_sq < 0.0) {
                    polarR[s](i,j) = -1.0;
                }
                else {
                    polarR[s](i,j)  = sqrt(abs(polarR_sq));
                }
            }      
        }
//      <---
    }
//    p1x1 = new valarray<float>(Npx);
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::p1x1_1D::p1x1_1D( const p1x1_1D& other) {

    for (size_t s(0); s < other.Species(); ++s) {
        pmin.push_back( other.Pmin(s) );
        pmax.push_back( other.Pmax(s) );
        Pl.push_back( other.PL(s) );
        polarR.push_back( other.PolarR(s) );
        p1x1.push_back( other.p1_x1(s) );
    }
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Output_Data::p1x1_1D::~p1x1_1D(){
   // delete p1x1;
}
//--------------------------------------------------------------

//--------------------------------------------------------------
//  Turn the Distribution function at some spatial location x0 
//  into a cartesian grid.
valarray<float>  Output_Data::p1x1_1D::operator()(DistFunc1D& df, size_t x0, size_t s) {

    p1x1[s] = 0.0; // this is a valarray<float>

    float p_im1(0.0), p_ip1(0.0);
    float integrant_low(0.0),
          integrant_high(0.0);

//  Calculate the integral for each harmonic separately
    for (size_t l(0); l < Nl(s); ++l) { 

        valarray<float> InPx( Npx(s) ); 
        for (size_t ipx(0); ipx < Npx(s)-1; ++ipx) { // at each location in px

//          At each location in |p|
            size_t ip(0);
            while (polarR[s](ipx,ip) < 0 ) { ++ip; }          

            p_ip1 = polarR[s](ipx,ip);
            integrant_high =(float( (*df(l))(ip,x0) )) * Pl[s](l)(ipx, ip);

            InPx[ipx] += p_ip1 * (0.5*p_ip1) * integrant_high; 
            ++ip;
            while ( (ip < Np(s) ) && (polarR[s](ipx,ip) > 0)  )   {  
                p_im1 = p_ip1;
                integrant_low = integrant_high;

                p_ip1 = polarR[s](ipx,ip);
                integrant_high =(float( (*df(l))(ip,x0) )) * Pl[s](l)(ipx, ip);
                InPx[ipx] += p_im1 * (0.5*(p_ip1-p_im1)) * integrant_low; 
                InPx[ipx] += p_ip1 * (0.5*(p_ip1-p_im1)) * integrant_high; 

                ++ip;
            }
            (p1x1[s])[ipx] += InPx[ipx];
        }
    }        
    p1x1[s] *= 2.0 * M_PI;

    return p1x1[s];
} 
//--------------------------------------------------------------
//**************************************************************


//**************************************************************
//--------------------------------------------------------------

void Output_Data::Output_Preprocessor_1D::operator()(const State1D& Y, const size_t tout) { 

    size_t Nx( Y.SH(0,0).numx() );
    valarray<float> dens( Nx );

    for (size_t s(0); s < Y.Species(); ++s){

//      Distribution function
        Array2D<float> px_x_data( px_x.Npx(s), Nx);
        for (size_t j(0); j < Nx; ++j) {
            valarray<float> data1D = px_x( Y.DF(s), j, s);
            for (size_t i(0); i < px_x.Npx(s); ++i) {
                px_x_data(i,j) = data1D[i];
            }
        }
        expo("px-x_"+ Export_Files::stringify(s), px_x_data, tout, s);

//      Density (moments)
        valarray<float> pra( Algorithms::MakeAxis( px_x.Pmin(s), px_x.Pmax(s), px_x.Np(s) ) );
        for (size_t i(0); i < Y.SH(0,0).numx(); ++i) {
            dens[i] = 4.0*M_PI*Algorithms::moment(  Export_Files::vfloat( (Y.SH(s,0)).xVec(i) ), pra, 2);
        }
        expo("n",dens, tout, s);

    }
//  Fields
    expo("Ex",Export_Files::vfloat(Y.Ex().array()), tout);

}
//--------------------------------------------------------------
//**************************************************************

