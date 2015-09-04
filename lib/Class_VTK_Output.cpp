

#include "Class_VTK.hpp"

using namespace std;


// =================================================================================== //
void VTK::Write( ){

    CalcAppendedOffsets() ;

    WriteMetaData() ;
    WriteData() ;

    if( nr_procs != 0  && my_proc == 0)  WriteCollection() ;

    fh.IncrementCounter() ;

    return ;
};

// =================================================================================== //
void VTK::WriteDataHeader( fstream &str, bool parallel ){

  string        location, line ;
  stringstream  scalars, vectors ;

  for( int j=0; j<2; j++){

    if( j==0 ) location = "Point" ;
    if( j==1 ) location = "Cell" ;

//    scalars.clear();
//    vectors.clear();

    scalars.str("");
    vectors.str("");

    //Creating Scalar and Vector Lists
    scalars << "\"" ;
    vectors << "\"" ;

    for( int i=0; i< nr_data; i++ ){

      if( data[i].GetLocation() == location){
        if(      data[i].GetComponents() == 1 ) scalars <<  data[i].GetName() << " " ;
        else if( data[i].GetComponents() == 3 ) vectors <<  data[i].GetName() << " " ;
      };

    };

    scalars << "\"" ;
    vectors << "\"" ;

    if(      location == "Point") {
      str << "      <" ;
      if( parallel )  str << "P" ;
      str << "PointData " ;
    }

    else if( location == "Cell")  {
      str << "      <" ;
      if( parallel )  str << "P" ;
      str << "CellData " ;
    };

    str << " Scalars=" << scalars.str()
        << " Vectors=" << vectors.str()
        << ">" << endl;

    //Writing DataArray
    for( int i=0; i< nr_data; i++){
      if( data[i].GetLocation() == location && !parallel) WriteDataArray( str, data[i] ) ;
      if( data[i].GetLocation() == location &&  parallel) WritePDataArray( str, data[i] ); 
    };

    str << "      </" ;
    if( parallel )  str << "P" ;

    if( location == "Point") str << "PointData> " << endl;
    if( location == "Cell")  str << "CellData> "  << endl;

  };

  return ;

};

// =================================================================================== //
void VTK::WriteDataArray( fstream &str, VTK::Field_C &field_ ){

    string          line;

    DataArrayToString( line, field_ ) ;

    str << line << endl ;
    
    str << "        </DataArray>" << endl ;

    return ;

};

// =================================================================================== //
void VTK::WritePDataArray( fstream &str, VTK::Field_C &field_ ){

    string          line;

    PDataArrayToString( line, field_ ) ;

    str << line << endl ;
    str << "        </DataArray>" << endl ;

    return ;

};

// =================================================================================== //
void VTK::WriteData( ){

    fstream             str ;
    fstream::pos_type   position_insert, position_eof ;

    int                 length;
    char*               buffer ;

    str.open( fh.GetName( ), ios::in | ios::out ) ;

    { // Write Ascii

        position_insert = str.tellg();
        VTK::Field_C    temp ;

        //Writing first point data then cell data
        for( int i=0; i< nr_data; i++){
            if( data[i].GetCodification() == "ascii" && data[i].GetLocation() == "Point") {
            str.seekg( position_insert);
            ReadDataArray( str, data[i] ) ;

            str.seekg( data[i].GetPosition() ) ;
            CopyUntilEOFInString( str, buffer, length );

            Flush( str, "ascii", data[i].GetName() ) ;

            position_insert = str.tellg();
            str << endl ;
            flush_binary( str, buffer, length) ;

            delete [] buffer ;

            };
        }; 

        for( int i=0; i< nr_data; i++){
            if( data[i].GetCodification() == "ascii" && data[i].GetLocation() == "Cell") {
            str.seekg( position_insert);
            ReadDataArray( str, data[i] ) ;

            str.seekg( data[i].GetPosition() ) ;
            CopyUntilEOFInString( str, buffer, length );

            Flush( str, "ascii", data[i].GetName() ) ;

            position_insert = str.tellg();
            str << endl ;
            flush_binary( str, buffer, length) ;

            delete [] buffer ;
            };
        }; 

        for( int i=0; i< geometry.size(); i++){
            if( geometry[i].GetCodification() == "ascii" ) {
            str.seekg( position_insert);
            ReadDataArray( str, geometry[i] ) ;

            str.seekg( geometry[i].GetPosition() ) ;
            CopyUntilEOFInString( str, buffer, length );

            Flush( str, "ascii", geometry[i].GetName() ) ;

            position_insert = str.tellg();
            str << endl ;
            flush_binary( str, buffer, length) ;

            delete [] buffer ;
            };
        }; 

        str.seekg( temp.GetPosition() ) ;


    }
    
    { // Write Appended

        char                c_;
        string              line ;
        fstream::pos_type   position_appended ;

        //Go to the initial position of the appended section
        while( getline(str, line) && (! Keyword_In_String( line, "<AppendedData")) ){} ;
        
        str >> c_;
        while( c_ != '_') str >> c_;
        
        position_insert = str.tellg();
        CopyUntilEOFInString( str, buffer, length );

        str.close();
        str.clear();
        
        
        //Reopening in binary mode
        str.open( fh.GetName( ), ios::out | ios::in | ios::binary);
        str.seekg( position_insert) ;

        //str.open( "data.dat", ios::out | ios::binary);

        //Writing first point data then cell data
        for( int i=0; i< nr_data; i++){
            if( data[i].GetCodification() == "appended" && data[i].GetLocation() == "Point") {
                if( GetHeaderType() == "UInt32"){
                    uint32_t    nbytes = data[i].GetNbytes() ;
                    flush_binary(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = data[i].GetNbytes() ;
                    flush_binary(str, nbytes) ;
                };
                Flush( str, "binary", data[i].GetName() ) ;
            };
        } 
        
        for( int i=0; i< nr_data; i++){
            if( data[i].GetCodification() == "appended" && data[i].GetLocation() == "Cell") {
                if( GetHeaderType() == "UInt32"){
                    uint32_t    nbytes = data[i].GetNbytes() ;
                    flush_binary(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = data[i].GetNbytes() ;
                    flush_binary(str, nbytes) ;
                };
                Flush( str, "binary", data[i].GetName() ) ;
            };
        } 
        
        //Writing Geometry Data
        for(int i=0; i<geometry.size(); i++){
        //for(int i=0; i<4; i++){
            if( geometry[i].GetCodification() == "appended" ) {
                if( GetHeaderType() == "UInt32"){
                    uint32_t    nbytes = geometry[i].GetNbytes() ;
                    flush_binary(str, nbytes) ;
                }

                else{
                    uint64_t    nbytes = geometry[i].GetNbytes() ;
                    flush_binary(str, nbytes) ;
                };
                Flush( str, "binary", geometry[i].GetName() ) ;           
            };
        };
   
        // { 
        // fstream             str2 ;
        // str2.open( "test2.dat", ios::out | ios::binary);
        // flush_binary( str2, buffer, length) ;
        // str2.close();
        // }

        flush_binary( str, buffer, length) ;

        delete [] buffer ;
    };
    
    // Closing Appended Secyion
    str.close();
    
    return ;

};
