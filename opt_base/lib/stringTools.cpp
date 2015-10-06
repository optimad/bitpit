#include"Operators.hpp"
using namespace std;      

bool Get_After_Keyword( string line_, string key_, char del_, string& result_){

  std::size_t c1, c2, pos ;

  string::iterator it;

  result_.clear() ;

  pos = line_.find( key_);

  if( pos == string::npos ){
      return false;
  }

  else{

    it= line_.begin();
    advance(it, pos);
    advance(it, key_.size() ) ;

    while( (*it) != del_){
      it++;
    };
    c1= it- line_.begin() +1;

    it++;

    while( (*it) != del_){
      it++;
    };
    c2= it- line_.begin()-1;

    pos= c2 -c1 +1;

    result_= line_.substr( c1, pos) ;
    trim( result_ ) ;

    return true;
  };

};
