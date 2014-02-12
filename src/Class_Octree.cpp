/*
 * Class_Octree.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: Edoardo Lombardi
 */

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "Class_Octree.hpp"

// =================================================================================== //
// CONSTRUCTORS                                                                        //
// =================================================================================== //
Class_Octant::Class_Octant(){
	x = y = z = 0;
	level = 0;
	marker = 0;
	bool dummy[16] = {false};
	memcpy(info, dummy, 16);
};

Class_Octant::Class_Octant(int8_t level, int32_t x, int32_t y, int32_t z){
	this->x = x;
	this->y = y;
	this->z = z;
	this->level = level;
	marker = 0;
	bool dummy[16] = {false};
	memcpy(info, dummy, 16);
};

Class_Octant::Class_Octant(Class_Octant & octant){
	x = octant.x;
	y = octant.y;
	z = octant.z;
	level = octant.level;
	marker = octant.marker;
	memcpy(info,octant.info,16);
};

// =================================================================================== //
// METHODS                                                                        //
// =================================================================================== //

// =================================================================================== //
// Get/Set methods
// =================================================================================== //

uint32_t Class_Octant::getx() const {
	return x;
}

uint32_t Class_Octant::gety() const {
	return y;
}

uint32_t Class_Octant::getz() const {
	return z;
}

uint8_t Class_Octant::getlevel() const {
	return level;
}

int8_t Class_Octant::getmarker() const {
	return marker;
}

bool Class_Octant::getbound(uint8_t face) const{
	return info[face];
}

bool Class_Octant::getpbound(uint8_t face) const{
	return info[6+face];
}

bool Class_Octant::getisnewR() const{
	return info[12];
}

bool Class_Octant::getisnewC() const{
	return info[13];
}

bool Class_Octant::getbalance() const{
	return info[14];
}

bool Class_Octant::getisghost() const{
	return info[15];
}

void Class_Octant::setmarker(int8_t marker) {
	this->marker = marker;
}

void Class_Octant::setbalance(bool balance) {
	info[14] = balance;
}

void Class_Octant::setlevel(uint8_t level) {
	this->level = level;
}

// =================================================================================== //
