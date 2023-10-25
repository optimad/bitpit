/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include"bitpit_common.hpp"

#include"fileHandler.hpp"

namespace bitpit{

/*!
 * Default constructor. 
 * Refers to a serial file in the current directory called "file.dat"
 * Default values of counters and blocks are set to 0.
 */
FileHandler::FileHandler()
    : directory("./"), name("file"), appendix("dat"),
      series(false), parallel(false),
      counter(0), block(0) {
}

/*!
 * Constructor. 
 * Refers to a serial file named according to input parameters.
 * Default values of counters and blocks are set to 0.
 * @param[in]   dir_    directory of file
 * @param[in]   name_   name of the file
 * @param[in]   app_    appendix of the file
 */
FileHandler::FileHandler(const std::string &dir_, const std::string &name_, const std::string &app_)
    : directory(dir_), name(name_), appendix(app_),
      series(false), parallel(false),
      counter(0), block(0) {
}

/*!
 * Copy Constructor. 
 * @param[in]   other    FileHandler to be copied.
 */
FileHandler::FileHandler(const FileHandler& other){

    (*this) = other ;

}

/*!
 * Assignment operator.
 */
FileHandler&    FileHandler::operator=(const FileHandler& other){

    directory   = other.directory ;
    name        = other.name ;
    appendix    = other.appendix ;
    series      = other.series ;
    parallel    = other.parallel ;
    counter     = other.counter ;
    block       = other.block ;


    return *this ;
}

/*!
 * Checks if file exists.
 * @return  [true/false] if file exists
 */
bool   FileHandler::exists() {
    std::ifstream f( getPath() );
    return f.good() ;
}

/*!
 * Composes the filename.
 * @return  complete filename
 */
std::string  FileHandler::getPath(){
  std::stringstream filename ;

  if (!directory.empty()) filename << directory << "/";
  filename << name;
  if(parallel) filename <<".b"<< bitpit::utils::string::zeroPadNumber(4, block)    ;
  if(series)   filename <<"." << bitpit::utils::string::zeroPadNumber(4, counter)  ;
  filename <<"."<< appendix  ;

  return filename.str() ;

}

/*!
 * Get the directory where the file will be saved
 * @return The directory where the file will be saved.
 */
std::string  FileHandler::getDirectory() const {

  return directory;
}

/*!
 * Get the name of the file
 * @return The name of the file.
 */
std::string  FileHandler::getName() const {

  return name;
}

/*!
 * Get the appendix of the file
 * @return The appendix of the file.
 */
std::string  FileHandler::getAppendix() const {

  return appendix;
}

/*!
 * sets the directory
 * @param[in]   d_      name of directory
 */
void    FileHandler::setDirectory(const std::string &d_){
    directory=d_;
    return;
}

/*!
 * sets the name
 * @param[in]   n_      name of file
 */
void    FileHandler::setName(const std::string &n_){
    name=n_;
    return;
}

/*!
 * sets the appendix
 * @param[in]   a_      appendix of file
 */
void    FileHandler::setAppendix(const std::string &a_){
    appendix=a_;
    return;
}

/*!
 * Checks if file belongs to a time series
 * @return Returns true if the file belongs to a time series, false otherwise
 */
bool    FileHandler::isSeries() const {
    return series;
}

/*!
 * sets if file belongs to a time series
 * @param[in]   s_      [true/false] if series
 */
void    FileHandler::setSeries( bool s_){ 
    series=s_; 
    return;
}

/*!
 * Get the time index of the following file.
 * If the file doen't belong to a time series, a negative number will
 * be returned.
 * @return counter
 */
int  FileHandler::getCounter() const {
  return counter;
}

/*!
 * sets the counter of the series
 * @param[in]   c_      index of series
 */
void    FileHandler::setCounter(int c_){ 
    counter=c_; 
    return; 
}

/*!
 * Increments the counter used for series.
 */
void    FileHandler::incrementCounter(){ 
    counter++; 
    return; 
}

/*!
 * Checks if file belongs to a parallel output
 * @return Returns true if the file belongs to a parallel output, false otherwise
 */
bool    FileHandler::isParallel() const {
    return parallel;
}

/*!
 * sets if file belongs to a parallel output
 * @param[in]   p_      [true/false] if parallel
 */
void    FileHandler::setParallel( bool p_){
    parallel=p_;
    return;
}

/*!
 * gets the index of the parallel block
 * @return the index of the parallel block
 */
int    FileHandler::getBlock( ) const {
    return block;
}

/*!
 * sets the counter of the series
 * @param[in]   b_      index of block of parallel set
 */
void    FileHandler::setBlock( int b_){
    block=b_;
    return;
}

}
