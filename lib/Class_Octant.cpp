/*
 * Class_Octree.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: Edoardo Lombardi
 */

// =================================================================================== //
// INCLUDES                                                                            //
// =================================================================================== //
#include "Class_Octant.hpp"

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

Class_Octant::Class_Octant(const Class_Octant &octant){
	x = octant.x;
	y = octant.y;
	z = octant.z;
	level = octant.level;
	marker = octant.marker;
	memcpy(info,octant.info,16);
};

// =================================================================================== //
// METHODS                                                                             //
// =================================================================================== //

// =================================================================================== //
// Basic Get/Set methods															   //
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
// Other Get/Set methods															   //
// =================================================================================== //

uint32_t Class_Octant::getsize() const {
	uint32_t size = uint32_t(pow(double(2),double(MAX_LEVEL-level)));
	return size;
}

uint32_t Class_Octant::getvolume() const {
	uint64_t volume = uint32_t(pow(double(this->getsize()),2.0));
	return volume;
}

void Class_Octant::buildchildren(vector<Class_Octant>& children) {
	children.clear();
	if (this->level < MAX_LEVEL){
		for (int i=0; i<nchildren; i++){
			switch (i) {
			case 0 :
			{
				Class_Octant oct(*this);
				oct.setmarker(max(0,oct.marker-1));
				oct.setlevel(oct.level+1);
				children.push_back(oct);
			}
			break;
			case 1 :
			{
				Class_Octant oct(*this);
				oct.setmarker(max(0,oct.marker-1));
				oct.setlevel(oct.level+1);
				uint32_t dh = oct.getsize();
				oct.x += dh;
				children.push_back(oct);
			}
			break;
			case 2 :
			{
				Class_Octant oct(*this);
				oct.setmarker(max(0,oct.marker-1));
				oct.setlevel(oct.level+1);
				uint32_t dh = oct.getsize();
				oct.y += dh;
				children.push_back(oct);
			}
			break;
			case 3 :
			{
				Class_Octant oct(*this);
				oct.setmarker(max(0,oct.marker-1));
				oct.setlevel(oct.level+1);
				uint32_t dh = oct.getsize();
				oct.x += dh;
				oct.y += dh;
				children.push_back(oct);
			}
			break;
			case 4 :
			{
				Class_Octant oct(*this);
				oct.setmarker(max(0,oct.marker-1));
				oct.setlevel(oct.level+1);
				uint32_t dh = oct.getsize();
				oct.z += dh;
				children.push_back(oct);
			}
			break;
			case 5 :
			{
				Class_Octant oct(*this);
				oct.setmarker(max(0,oct.marker-1));
				oct.setlevel(oct.level+1);
				uint32_t dh = oct.getsize();
				oct.x += dh;
				oct.z += dh;
				children.push_back(oct);
			}
			break;
			case 6 :
			{
				Class_Octant oct(*this);
				oct.setmarker(max(0,oct.marker-1));
				oct.setlevel(oct.level+1);
				uint32_t dh = oct.getsize();
				oct.y += dh;
				oct.z += dh;
				children.push_back(oct);
			}
			break;
			case 7 :
			{
				Class_Octant oct(*this);
				oct.setmarker(max(0,oct.marker-1));
				oct.setlevel(oct.level+1);
				uint32_t dh = oct.getsize();
				oct.x += dh;
				oct.y += dh;
				oct.z += dh;
				children.push_back(oct);
			}
			break;
			}
		}
	}
	else{
		Write_Log("Max level reached ---> No Children Built");
	}
}

// =================================================================================== //

