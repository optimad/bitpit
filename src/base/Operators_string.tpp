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
/*!
    Left-trim operator for std::string.
    Remove left trailing spaces from string. For instance, if the input string is
    "  test_string  ", on output this function returns "test_string  "

    \param[in] s input string

    \result reference to input string
*/
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
/*!
    Right-trim operator for std::string.
    Remove right blank spaces from string. For instance, if the input string is
    "  test_string  ", on output this function returns "  test_string"

    \param[in] s input string

    \result reference to input string
*/
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
/*!
    Trim operator for std::string.
    Remove left/right blank spaces from string. For instance, if the input string is
    "  test_string  ", on output this function returns "test_string"

    \param[in] s input string

    \result reference to input string
*/
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
/*!
    Given an integer, returns a string of length nchar, composed of the input number
    and nchar - ndigits '0' characters (where ndigits is the number of digits of the input integer)
    in the following format "000xxx".
    If ndigits > nchar, the output string will contaiens ndigits characters storing the
    digits of the input number.
    For instance, if nchar = 4 and num = 3, this function returns the string "0003".
    If nchar = 4, and num = 12345, this function returns "12345".

    \param[in] nchar string length
    \param[in] num input integer

    \result returns a string storing the input number in the format 000xxx.
*/
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
/*!
    Check whether a string contains the kwyword or not.

    \param[in] line_ input string
    \param[in] key_ search key
    
    \result boolean flag (true) if the keyword has been found, (false) otherwise.
*/
static inline bool Keyword_In_String( string line_, string key_){
  bool exist= ( line_.find( key_) != string::npos ) ;
  return exist  ;
};

// ---------------------------------------------------------------------------------- //
/*!
    Convertes a string into fundamental data type.

    If no data of type T can be extracted from the input string a 0 value,
    will be stored in output_.
    If multiple values can be extracted from the input string, only the first
    value will be saved in output_.

    \param[in] input_ input string
    \param[in,out] output_ variable storing the value read from string.
*/
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

// ---------------------------------------------------------------------------------- //
/*!
    Convertes a string into a vector of fundamental data type.

    If no data of type T can be extracted from the input string a void vector is returned.
    Values extracted from string are pushed at the end of the vector.

    \param[in] input_ input string
    \param[in,out] output_ vector storing the value extracted from string.
*/
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

// ---------------------------------------------------------------------------------- //
/*!
    Convertes a string into a arrayof fundamental data type.

    If no data of type T can be extracted from the input string a void array with null
    elements is returned.
    If the number of elements which can be extracted from the input string is larger
    than the array size, only the first n elements are saved in the array.

    \param[in] input_ input string
    \param[in,out] output_ array storing the value extracted from string.
*/
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
