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
#include <iostream>
// =================================================================================== //
// CONSTRUCTORS AND OPERATORS                                                          //
// =================================================================================== //
Class_Octant::Class_Octant(){
	x = y = z = 0;
	level = 0;
	marker = 0;
	bool dummy[16] = {false};
	memcpy(info, dummy, 16);
	for (int i=0; i<nface; i++){
		info[i] = true;
	}
};

Class_Octant::Class_Octant(int8_t level, int32_t x, int32_t y, int32_t z){
	this->x = x;
	this->y = y;
	this->z = z;
	this->level = level;
	marker = 0;
	bool dummy[16] = {false};
	memcpy(info, dummy, 16);
	if (level==0){
		for (int i=0; i<nface; i++){
			info[i] = true;
		}
	}
};

Class_Octant::Class_Octant(const Class_Octant &octant){
	x = octant.x;
	y = octant.y;
	z = octant.z;
	level = octant.level;
	marker = octant.marker;
	memcpy(info,octant.info,16);
};

bool Class_Octant::operator ==(const Class_Octant & oct2){
	bool check = true;
	check = check && (x == oct2.x);
	check = check && (y == oct2.y);
	check = check && (z == oct2.z);
	check = check && (level == oct2.level);
	return check;
}

// =================================================================================== //
// METHODS                                                                             //
// =================================================================================== //

// =================================================================================== //
// Basic Get/Set methods															   //
// =================================================================================== //

uint32_t Class_Octant::getX() const {
	return x;
}

uint32_t Class_Octant::getY() const {
	return y;
}

uint32_t Class_Octant::getZ() const {
	return z;
}

uint8_t Class_Octant::getLevel() const {
	return level;
}

int8_t Class_Octant::getMarker() const {
	return marker;
}

bool Class_Octant::getBound(uint8_t face) const{
	return info[face];
}

bool Class_Octant::getPbound(uint8_t face) const{
	return info[6+face];
}

bool Class_Octant::getIsNewR() const{
	return info[12];
}

bool Class_Octant::getIsNewC() const{
	return info[13];
}

bool Class_Octant::getNotBalance() const{
	return info[14];
}

bool Class_Octant::getIsGhost() const{
	return info[15];
}

void Class_Octant::setMarker(int8_t marker) {
	this->marker = marker;
}

void Class_Octant::setBalance(bool balance) {
	info[14] = balance;
}

void Class_Octant::setLevel(uint8_t level) {
	this->level = level;
}

// =================================================================================== //
// Other Get/Set methods															   //
// =================================================================================== //

uint32_t Class_Octant::getSize() const {
	uint32_t size = uint32_t(pow(double(2),double(MAX_LEVEL-level)));
	return size;
}

uint32_t Class_Octant::getArea() const {
	uint32_t area = uint32_t(pow(double(getSize()),2.0));
	return area;
}

uint32_t Class_Octant::getVolume() const {
	uint64_t volume = uint32_t(pow(double(getSize()),3.0));
	return volume;
}


uint32_t (*Class_Octant::getNodes())[DIM] {
	uint8_t		i, cx, cy, cz;
	uint32_t	dh;

	dh = getSize();
	uint32_t (*nodes)[DIM] = new uint32_t[nnodes][DIM];

	for (i = 0; i < nnodes; i++){
		cx = uint8_t(i%2);
		cy = uint8_t((i-4*(i/4))/2);
		cz = uint8_t(i/4);
		nodes[i][0] = x + cx*dh;
		nodes[i][1] = y + cy*dh;
		if (DIM ==3){
			nodes[i][2] = z + cz*dh;
		}
	}
	return nodes;
}

double*	Class_Octant::getCenter() {
	uint8_t		i;
	double	dh;

	dh = double(getSize())/2.0;
	double *center = new double[DIM];

	center[0] = (double)x + dh;
	center[1] = (double)y + dh;
#if DIM == 3
	center[2] = (double)z + dh;
#endif
	return center;
}

// =================================================================================== //
// Other methods													    			   //
// =================================================================================== //

Class_Octant Class_Octant::buildLastDesc() {
	uint32_t delta = (uint32_t)pow(2.0,(double)((uint8_t)MAX_LEVEL - level)) - 1;
	Class_Octant last_desc(MAX_LEVEL,x+delta,y+delta,z+delta);
	return last_desc;
}

// =================================================================================== //

Class_Octant Class_Octant::buildFather() {
	uint32_t deltax = x%(uint32_t(pow(2.0,(double)((uint8_t)MAX_LEVEL - (level-1)))));
	uint32_t deltay = y%(uint32_t(pow(2.0,(double)((uint8_t)MAX_LEVEL - (level-1)))));
#if DIM == 3
	uint32_t deltaz = z%(uint32_t(pow(2.0,(double)((uint8_t)MAX_LEVEL - (level-1)))));
#endif
	Class_Octant father(level-1, x-deltax, y-deltay, z-deltaz);
	return father;
}

// =================================================================================== //

Class_Octant* Class_Octant::buildChildren() {
	uint8_t xf,yf,zf;


	if (this->level < MAX_LEVEL){
		Class_Octant* children = new Class_Octant[nchildren];
		for (int i=0; i<nchildren; i++){
			switch (i) {
			case 0 :
			{
				Class_Octant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				// Update interior face bound and pbound
				xf=1; yf=3; zf=5;
				oct.info[xf] = oct.info[xf+nface] = false;
				oct.info[yf] = oct.info[yf+nface] = false;
				oct.info[zf] = oct.info[zf+nface] = false;
				children[0] = oct;
			}
			break;
			case 1 :
			{
				Class_Octant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize();
				oct.x += dh;
				// Update interior face bound and pbound
				xf=0; yf=3; zf=5;
				oct.info[xf] = oct.info[xf+nface] = false;
				oct.info[yf] = oct.info[yf+nface] = false;
				oct.info[zf] = oct.info[zf+nface] = false;
				children[1] = oct;
			}
			break;
			case 2 :
			{
				Class_Octant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize();
				oct.y += dh;
				// Update interior face bound and pbound
				xf=1; yf=2; zf=5;
				oct.info[xf] = oct.info[xf+nface] = false;
				oct.info[yf] = oct.info[yf+nface] = false;
				oct.info[zf] = oct.info[zf+nface] = false;
				children[2] = oct;
			}
			break;
			case 3 :
			{
				Class_Octant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize();
				oct.x += dh;
				oct.y += dh;
				// Update interior face bound and pbound
				xf=0; yf=2; zf=5;
				oct.info[xf] = oct.info[xf+nface] = false;
				oct.info[yf] = oct.info[yf+nface] = false;
				oct.info[zf] = oct.info[zf+nface] = false;
				children[3] = oct;
			}
			break;
			case 4 :
			{
				Class_Octant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize();
				oct.z += dh;
				// Update interior face bound and pbound
				xf=1; yf=3; zf=4;
				oct.info[xf] = oct.info[xf+nface] = false;
				oct.info[yf] = oct.info[yf+nface] = false;
				oct.info[zf] = oct.info[zf+nface] = false;
				children[4] = oct;
			}
			break;
			case 5 :
			{
				Class_Octant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize();
				oct.x += dh;
				oct.z += dh;
				// Update interior face bound and pbound
				xf=0; yf=3; zf=4;
				oct.info[xf] = oct.info[xf+nface] = false;
				oct.info[yf] = oct.info[yf+nface] = false;
				oct.info[zf] = oct.info[zf+nface] = false;
				children[5] = oct;
			}
			break;
			case 6 :
			{
				Class_Octant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize();
				oct.y += dh;
				oct.z += dh;
				// Update interior face bound and pbound
				xf=1; yf=2; zf=4;
				oct.info[xf] = oct.info[xf+nface] = false;
				oct.info[yf] = oct.info[yf+nface] = false;
				oct.info[zf] = oct.info[zf+nface] = false;
				children[6] = oct;
			}
			break;
			case 7 :
			{
				Class_Octant oct(*this);
				oct.setMarker(max(0,oct.marker-1));
				oct.setLevel(oct.level+1);
				oct.info[12]=true;
				uint32_t dh = oct.getSize();
				oct.x += dh;
				oct.y += dh;
				oct.z += dh;
				// Update interior face bound and pbound
				xf=0; yf=2; zf=4;
				oct.info[xf] = oct.info[xf+nface] = false;
				oct.info[yf] = oct.info[yf+nface] = false;
				oct.info[zf] = oct.info[zf+nface] = false;
				children[7] = oct;
			}
			break;
			}
		}
		return children;
	}
	else{
		Class_Octant* children = new Class_Octant[0];
		writeLog("Max level reached ---> No Children Built");
		return children;
	}
}

// =================================================================================== //

uint64_t Class_Octant::computeMorton() const {
	uint64_t morton = 0;
	morton = mortonEncode_magicbits(this->x,this->y,this->z);
	return morton;
}
uint64_t Class_Octant::computeMorton(){
	uint64_t morton = 0;
	morton = mortonEncode_magicbits(this->x,this->y,this->z);
	return morton;
}

// =================================================================================== //

uint64_t* Class_Octant::computeHalfSizeMorton(uint8_t iface,
											 uint32_t & sizehf) {
	uint32_t dh,dh2;
	uint64_t morton;
	uint32_t nneigh;
	uint8_t i,cx,cy,cz;

	nneigh = (level < MAX_LEVEL) ? nchildren/2 : 1;
	dh = (level < MAX_LEVEL) ? getSize()/2 : getSize();
	dh2 = getSize();

	if (info[iface]){
//		uint64_t* Morton = new uint64_t[0];
		sizehf = 0;
//		return Morton;
		return NULL;
	}
	else{
		uint64_t* Morton = new uint64_t[nneigh];
		switch (iface) {
		case 0 :
		{
			for (i=0; i<nneigh; i++){
				cy = (i==1)||(i==3);
				cz = (i==2)||(i==3);
				Morton[i] = mortonEncode_magicbits(this->x-dh,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 1 :
		{
			for (i=0; i<nneigh; i++){
				cy = (i==1)||(i==3);
				cz = (i==2)||(i==3);
				Morton[i] = mortonEncode_magicbits(this->x+dh2,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 2 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1)||(i==3);
				cz = (i==2)||(i==3);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y-dh,this->z+dh*cz);
			}
		}
		break;
		case 3 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1)||(i==3);
				cz = (i==2)||(i==3);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2,this->z+dh*cz);
			}
		}
		break;
		case 4 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1)||(i==3);
				cy = (i==2)||(i==3);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z-dh);
			}
		}
		break;
		case 5 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i==1)||(i==3);
				cy = (i==2)||(i==3);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh2);
			}
		}
		break;
		}
		sizehf = nneigh;
		return Morton;
	}

}

// =================================================================================== //

uint64_t* Class_Octant::computeMinSizeMorton(uint8_t iface,
											const uint8_t & maxdepth,
											uint32_t & sizem) {
	uint32_t dh,dh2;
	uint64_t morton;
	uint32_t nneigh, nline;
	uint32_t i,cx,cy,cz;

	nneigh = (level < MAX_LEVEL) ? uint32_t(pow(2.0,double((DIM-1)*(maxdepth-level)))) : 1;
	dh = (level < MAX_LEVEL) ? uint32_t(pow(2.0,double(MAX_LEVEL - maxdepth))) : getSize();
	dh2 = getSize();
	nline = uint32_t(pow(2.0,double((maxdepth-level))));

	if (info[iface]){
//		uint64_t* Morton = new uint64_t[0];
		sizem = 0;
//		return Morton;
		return NULL;
	}
	else{
		uint64_t* Morton = new uint64_t[nneigh];
		switch (iface) {
		case 0 :
		{
			for (i=0; i<nneigh; i++){
				cy = (i/nline);
				cz = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->x-dh,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 1 :
		{
			for (i=0; i<nneigh; i++){
				cy = (i/nline);
				cz = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->x+dh2,this->y+dh*cy,this->z+dh*cz);
			}
		}
		break;
		case 2 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i/nline);
				cz = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y-dh,this->z+dh*cz);
			}
		}
		break;
		case 3 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i/nline);
				cz = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh2,this->z+dh*cz);
			}
		}
		break;
		case 4 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i/nline);
				cy = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z-dh);
			}
		}
		break;
		case 5 :
		{
			for (i=0; i<nneigh; i++){
				cx = (i/nline);
				cy = (i%nline);
				Morton[i] = mortonEncode_magicbits(this->x+dh*cx,this->y+dh*cy,this->z+dh2);
			}
		}
		break;
		}
		sizem = nneigh;
		sort(Morton,Morton+nneigh);
		return Morton;
	}
}

void Class_Octant::setPbound(uint8_t face, bool flag) {
	info[nface+face] = flag;
}

uint64_t* Class_Octant::computeVirtualMorton(uint8_t iface,
											const uint8_t& maxdepth,
											uint32_t& sizeneigh) {

	if (getNotBalance()){
		return computeMinSizeMorton(iface,
									maxdepth,
									sizeneigh);
	}
	else{
		return computeHalfSizeMorton(iface,
									sizeneigh);
	}
		}

// =================================================================================== //
