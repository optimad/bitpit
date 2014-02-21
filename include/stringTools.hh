// =================================================================================== //
// STRING TOOLS                                                                        //
//                                                                                     //
// Functions for string handling.                                                      //
//                                                                                     //
// LIST OF FUNCTIONS                                                                   //
// - ltrim     : string left trimming                                                  //
// - rtrim     : string right trimming                                                 //
// - trim      : string trimming                                                       //
// =================================================================================== //
// INFO                                                                                //
// =================================================================================== //
// Author   : Marco Cisternino                                                         //
// Company  : Optimad Engineering srl                                                  //
// Date     : Jul 2, 2013                                                              //
// Version  : v 1.0                                                                    //
//                                                                                     //
// All rights reserved                                                                 //
// =================================================================================== //

// =================================================================================== //
// PRE-COMPILATION INSTRUCTIONS                                                        //
// =================================================================================== //
#ifndef STRINGTOOLS_HH
#define STRINGTOOLS_HH

// includes
#include<algorithm>
#include<string>
#include<functional>
#include<cctype>

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //
// Left Trimming
static inline std::string &ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //
// Right trimming
static inline std::string &rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}

// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //
// String trimming (left & right)
static inline std::string &trim(std::string &s) {
	return ltrim(rtrim(s));
}


#endif
