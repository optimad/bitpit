// ================================================================================== //
//                               STRING OPERATORS                                     //
//                                                                                    //
// Basic operators for string.                                                        //
// ================================================================================== //
// INFO                                                                               //
// ================================================================================== //
// Author     : Alessandro Alaia                                                      //
// Date       : Jul 16, 2013                                                          //
// Version    : v1.0                                                                  //
//                                                                                    //
// All rights reserved.                                                               //
// ================================================================================== //

// Trimming operators =============================================================== //

// ---------------------------------------------------------------------------------- //
static inline string &ltrim(
    string &s
) {

// ================================================================================== //
// static inline string &ltrim(                                                       //
//     string &s)                                                                     //
//                                                                                    //
// String left trimming. Remove left trailing spaces from string.                     //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - s           : string, input string                                               //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - s           : string, trimmed string                                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
// none

// Counters
// none

// ================================================================================== //
// TRIM STRING                                                                        //
// ================================================================================== //
s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));

return(s); };

// ---------------------------------------------------------------------------------- //
static inline string &rtrim(
    string &s
) {

// ================================================================================== //
// static inline string &rtrim(                                                       //
//     string &s)                                                                     //
//                                                                                    //
// String right trimming. Remove right trailing spaces from string.                   //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - s           : string, input string                                               //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - s           : string, trimmed string                                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
// none

// Counters
// none

// ================================================================================== //
// TRIM STRING                                                                        //
// ================================================================================== //
s.erase(find_if(s.rbegin(), s.rend(), not1(std::ptr_fun<int, int>(isspace))).base(), s.end());

return(s); };

// ---------------------------------------------------------------------------------- //
static inline string &trim(
    string &s
) {

// ================================================================================== //
// static inline string &trim(                                                        //
//     string &s)                                                                     //
//                                                                                    //
// String trimming. Remove left & right trailing spaces from string.                  //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - s           : string, input string                                               //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - s           : string, trimmed string                                             //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

// Local variables
// none

// Counters
// none

// ================================================================================== //
// TRIM STRING                                                                        //
// ================================================================================== //
return(ltrim(rtrim(s))); };

// Padding operators ================================================================ //

// ---------------------------------------------------------------------------------- //
static inline string ZeroPadNumber(
    int nchar,
    int num
) {

// ================================================================================== //
// static inline string ZeroPadNumber(                                                //
//     int nchar,                                                                     //
//     int num)                                                                       //
//                                                                                    //
// Covert integer into string with format 000xxx.                                     //
// ================================================================================== //
// INPUT                                                                              //
// ================================================================================== //
// - nchar       : number of character in the final string                            //
// - num         : int, integer to be padded                                          //
// ================================================================================== //
// OUTPUT                                                                             //
// ================================================================================== //
// - s           : string, with padded padded integer                                 //
// ================================================================================== //

// ================================================================================== //
// VARIABLES DECLARATION                                                              //
// ================================================================================== //

//  Local variables
ostringstream          ss;

// Counters
// none

// ================================================================================== //
// PADDING                                                                            //
// ================================================================================== //
ss << setw(nchar) << setfill('0') << num;

return (ss.str()); };

// ---------------------------------------------------------------------------------- //
static inline bool Keyword_In_String( string line_, string key_){
  bool exist= ( line_.find( key_) != string::npos ) ;
  return exist  ;
};
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //
// converts a string to fundamental data types and vectors or arrays of them

//---------------------------------------------------------------------------------
template <class T>
void  convert_string( string input_, T &output_ ){

  vector<T>    temp_;
  T            x_;

  temp_.clear() ;

  trim( input_ );
  stringstream ss_( input_ );

  while( ss_.good() ){
    ss_ >> x_ ;
    temp_.push_back(x_);
  };

 
  if( temp_.size() == 0 ){
    cout << " no useful information in string " << input_   << endl;
    cout << " casting zero                   " <<  endl;

    x_ = static_cast<T> (0) ; 
  } 

  else if( temp_.size() == 1){
    x_ = temp_[0] ;
  }

  else if( temp_.size() > 1){
    cout << " more than one element in string " << input_   << endl;
    cout << " assigning first element             "  << endl;

    x_ = temp_[0] ; 
  };


  output_ = x_ ;
  return ;

};

//---------------------------------------------------------------------------------
template <class T>
void  convert_string( string input_, vector<T> &output_){

  vector<T>    temp_;
  T            x_;

  temp_.clear() ;

  trim( input_ );
  stringstream ss_( input_ );

  while( ss_.good() ){
    ss_ >> x_ ;
    temp_.push_back(x_);
  };

 
  if( temp_.size() == 0 ){
    cout << " no useful information in string " << input_   << endl;
    cout << " returning void vector          " <<  endl;
  } ;

  output_= temp_ ;
  return ;
};

//---------------------------------------------------------------------------------
template <class T, size_t n>
void  convert_string( string input_, array<T,n> &output_) {

  vector<T>    temp_;
  T            x_;

  temp_.clear() ;

  trim( input_ );
  stringstream ss_( input_ );

  while( ss_.good() ){
    ss_ >> x_ ;
    temp_.push_back(x_);
  };


  if( temp_.size() < n ){
    cout << " not enough useful information in string " << input_   << endl;
    cout << " casting zero into missing elements      " <<  endl;

    x_ = static_cast<T> (0) ;
    output_.fill( x_ ) ;

    for( int i=0; i<temp_.size(); i++){
      output_[i] = temp_[i] ;
    };
  }

  else if( temp_.size() == n){
    for( int i=0; i<n; i++){
      output_[i] = temp_[i] ;
    };
  }

  else if( temp_.size() > n){
    cout << " more than " << n << " elements in string " << input_   << endl;
    cout << " assigning first element " << n << " elements "   << endl;

    for( int i=0; i<n; i++){
      output_[i] = temp_[i] ;
    };
  };


  return ;
};
