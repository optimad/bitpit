#include"Operators.hpp"
using namespace std;      

string Get_After_Keyword( string line_, string key_, char del_){

  unsigned c1, c2, pos ;
  string val;

  string::iterator it;

  pos = line_.find( key_);

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

  val= line_.substr( c1, pos) ;
  trim( val ) ;

  return val;

};
