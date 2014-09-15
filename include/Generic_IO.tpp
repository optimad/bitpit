// =================================================================================== //
// For Writing

template< class data_T >
void flush_ascii( fstream &str, const data_T data ){

  str << setprecision(8) << scientific ;
  str << data << " ";

  return ;
};

// =================================================================================== //
template< class data_T >
void flush_ascii( fstream &str, int elements_per_line, const vector<data_T> &data ){

  int i(0), j(0), k(0) ;
  int nr ;

  int lines, this_line ;

  bool next(true) ;

  nr = data.size() ;
  lines = (nr-1) /elements_per_line + 1;

  str << setprecision(8) << scientific ;

  
  while( next ) {

    this_line = min( elements_per_line, nr - i * elements_per_line ) ;
    for( j=0; j<this_line; j++){
      str << data[k] << " ";
      k++ ;
    };

    i++ ;
    if( i<lines){
      str << endl;
    }
    else{
      next = false;
    };

  };

  return ;
};

// =================================================================================== //
template< class data_T, size_t d >
void flush_ascii( fstream &str, int elements_per_line, const array<data_T,d> &data ){

  int i(0), j(0), k(0) ;
  int nr ;

  int lines, this_line ;

  bool next(true) ;

  nr = d ;
  lines = nr /elements_per_line ;

  str << setprecision(8) << scientific ;

  
  while( next ) {

    this_line = min( elements_per_line, nr - i * elements_per_line ) ;
    for( j=0; j<this_line; j++){
      str << data[k] << " ";
      k++ ;
    };

    i++ ;
    if( i<lines){
      str << endl;
    }
    else{
      next = false;
    };

  };

  return ;
};

// =================================================================================== //
template< class data_T >
void flush_ascii( fstream &str, int elements_per_line, const data_T *data, int nr ){

  int i(0), j(0), k(0) ;

  int lines, this_line ;

  bool next(true) ;

  lines = nr /elements_per_line ;

  str << setprecision(8) << scientific ;

  
  while( next ) {

    this_line = min( elements_per_line, nr - i * elements_per_line ) ;
    for( j=0; j<this_line; j++){
      str << data[k] << " ";
      k++ ;
    };

    i++ ;
    if( i<lines){
      str << endl;
    }
    else{
      next = false;
    };

  };

  return ;
};

// =================================================================================== //
template< class data_T >
void flush_binary( fstream &str, const data_T data ){

  int nbytes;
  nbytes = sizeof(data_T) ;

  str.write( reinterpret_cast<const char*>(&data), nbytes ) ;

  return ;
};

// =================================================================================== //
template< class data_T >
void flush_binary( fstream &str, const vector<data_T> &data ){

  int i, nbytes, nr;
  nr = data.size() ;
  nbytes = sizeof(data_T) ;


  for(i=0; i<nr; i++){
    str.write( reinterpret_cast<const char*>(&data[i]), nbytes ) ;
  };

  return ;
};

// =================================================================================== //
template< class data_T, size_t d >
void flush_binary( fstream &str, const array<data_T,d> &data ){

  int nbytes;
  array<int,8> data1 ;
  nbytes = sizeof(data_T)*d ;

  str.write( reinterpret_cast<const char*>(&data[0]), nbytes ) ;

  return ;
};

// =================================================================================== //
template< class data_T >
void flush_binary( fstream &str, int nr, const data_T *data ){

  int i, nbytes;
  nbytes = sizeof(data_T) *nr ;

  str.write( reinterpret_cast<const char*>(data), nbytes ) ;

  return ;
};

// =================================================================================== //
// FOR READING

template< class data_T >
void  line_stream( fstream &str, data_T &data){

  vector<data_T> temp;
  data_T         x_;
  string         line;
  int            expected, read(0) ;

  expected = 1;

  getline( str, line );
  trim( line );

  stringstream ss( line );

  while( ss.good() ){
    ss >> x_ ;
    temp.push_back(x_);
    read++ ;
  };

  if( read != expected){
    cout << " Not expected nr of element in line" << endl;
    cout << " Expected number: "<< expected << endl ; 
    cout << " Actual number: "<< read << endl ; 
  }

  else{
    data=temp[0];
  };

  return;

};

// =================================================================================== //
template< class data_T >
void  line_stream( fstream &str, vector<data_T> &data){

  vector<data_T> temp;
  data_T         x_;
  string         line;
  int            expected, read(0) ;

  expected = data.size() ;

  getline( str, line );
  trim( line );

  stringstream ss( line );
  
  while( ss.good() ){
    ss >> x_ ;
    temp.push_back(x_);
    read++ ;
  };

  if( expected == 0) {
    data = temp ;
  }

  else{
    if( expected == read){
      data = temp ;
    }
    else{
      cout << " Not expected nr of element in line" << endl;
      cout << " Expected number: "<< expected << endl ; 
      cout << " Actual number: "<< read << endl ; 
    };
  };

  return;
  
};

// =================================================================================== //
template< class data_T, size_t d >
void  line_stream( fstream &str, array<data_T,d> &data){

  vector<data_T> temp;
  data_T         x_;
  string         line;
  int            expected, read(0), i ;

  expected = d ;

  getline( str, line );
  trim( line );

  stringstream ss( line );
  
  while( ss.good() ){
    ss >> x_ ;
    temp.push_back(x_);
    read++ ;
  };

  if( expected == read){
    for(i=0; i<read; i++) data[i] = temp[i] ;
  }
  else{
    cout << " Not expected nr of element in line" << endl;
    cout << " Expected number: "<< expected << endl ; 
    cout << " Actual number: "<< read << endl ; 
  };

  return;
  
};

// =================================================================================== //
template< class data_T >
void  line_stream( fstream &str, data_T *data, int nr ){

  vector<data_T> temp;
  data_T         x_;
  string         line;
  int            expected, read(0), i ;

  expected = nr ;

  getline( str, line );
  trim( line );

  stringstream ss( line );

  while( ss.good() ){
    ss >> x_ ;
    temp.push_back(x_);
    read++ ;
  };

  if( expected == read){
    for(i=0; i<read; i++) data[i] = temp[i] ;
  }
  else{
    cout << " Not expected nr of element in line" << endl;
    cout << " Expected number: "<< expected << endl ; 
    cout << " Actual number: "<< read << endl ; 
  };
  
  return;

};


// =================================================================================== //
template< class data_T >
void absorb_ascii( fstream &str, data_T &data ){

  str >> data ;

  return ;
};

// =================================================================================== //
template< class data_T >
void absorb_ascii( fstream &str, vector<data_T> &data ){

  vector<data_T>   temp;
  int              i ;
  int              nr, read, new_ ;

  read = 0 ;
  nr = data.size() ;


  while( str.good() ) {

    line_stream( str, temp) ;
    new_ = temp.size() ; 

    for( i=0; i<new_; i++){
      if( read < nr ){
        data[read] = temp[i] ;
        read++ ;
      };
    };

  };

  if( read != nr) {
    cout << "Not enough elements found to fill vector" << endl ;
    cout << "Size of vector: " << nr << endl ;
    cout << "Elements within file: " << read << endl ;
  };


  return ;
};

// =================================================================================== //
template< class data_T, size_t d >
void absorb_ascii( fstream &str, array<data_T,d> &data ){

  vector<data_T>   temp;
  int              i ;
  int              nr, read, new_ ;

  read = 0 ;
  nr = d ;


  while( str.good() && read<nr ) {

    line_stream( str, temp) ;
    new_ = temp.size() ; 

    for( i=0; i<new_; i++){
      if( read < nr ){
        data[read] = temp[i] ;
        read++ ;
      };
    };

  };

  if( read != nr) {
    cout << "absorb_ascii<> array overloading" << endl ;
    cout << "Not enough elements found to fill array" << endl ;
    cout << "Size of array: " << nr << endl ;
    cout << "Elements within file: " << read << endl ;
  };



  return ;
};

// =================================================================================== //
template< class data_T >
void absorb_ascii( fstream &str, data_T *data, int nr ){

  vector<data_T>   temp;
  int              i ;
  int              read, new_ ;

  read = 0 ;

  while( str.good() ) {

    line_stream( str, temp) ;
    new_ = temp.size() ; 

    for( i=0; i<new_; i++){
      if( read < nr ){
        data[read] = temp[i] ;
        read++ ;
      };
    };

  };

  if( read != nr) {
    cout << "Not enough elements found to fill c array" << endl ;
    cout << "Size of c array: " << nr << endl ;
    cout << "Elements within file: " << read << endl ;
  };


  return ;
};

// =================================================================================== //
template< class data_T >
void absorb_binary( fstream &str, data_T &data ){

  int nbytes ;
  nbytes = sizeof(data_T) ;

  str.read( reinterpret_cast<char*>(&data), nbytes ) ;

  return ;
};


// =================================================================================== //
template< class data_T >
void absorb_binary( fstream &str, vector<data_T> &data ){

  int i, nbytes, nr;
  nr = data.size() ;
  nbytes = sizeof(data_T) ;

  for(i=0; i<nr; i++){
    str.read( reinterpret_cast<char*>(&data[i]), nbytes ) ;
  };

  return ;
};

// =================================================================================== //
template< class data_T, size_t d >
void absorb_binary( fstream &str, array<data_T,d> &data ){

  int nbytes;
  nbytes = sizeof(data_T) *d ;

  str.read( reinterpret_cast<char*>(&data[0]), nbytes ) ;
  
  return ;
};

// =================================================================================== //
template< class data_T >
void absorb_binary( fstream &str, data_T *data, int nr ){

  int i, nbytes;
  nbytes = sizeof(data_T) *nr ;

  str.read( reinterpret_cast<char*>(data), nbytes ) ;

  return ;
};
