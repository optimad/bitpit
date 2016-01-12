#ifndef INLINEDFUNCT_HPP_
#define INLINEDFUNCT_HPP_

#include <stdint.h>
#include <limits.h>

// method to seperate bits from a given integer 3 positions apart
inline uint64_t splitBy3(unsigned int a){
	uint64_t x = a & 0x1fffff; // we only look at the first 21 bits
	x = (x | x << 32) & 0x1f00000000ffff;  // shift left 32 bits, OR with self, and 00011111000000000000000000000000000000001111111111111111
	x = (x | x << 16) & 0x1f0000ff0000ff;  // shift left 32 bits, OR with self, and 00011111000000000000000011111111000000000000000011111111
	x = (x | x << 8) & 0x100f00f00f00f00f; // shift left 32 bits, OR with self, and 0001000000001111000000001111000000001111000000001111000000000000
	x = (x | x << 4) & 0x10c30c30c30c30c3; // shift left 32 bits, OR with self, and 0001000011000011000011000011000011000011000011000011000100000000
	x = (x | x << 2) & 0x1249249249249249;
	return x;
}

inline uint64_t mortonEncode_magicbits(unsigned int x, unsigned int y, unsigned int z){
	uint64_t answer = 0;
	answer |= splitBy3(x) | splitBy3(y) << 1 | splitBy3(z) << 2;
	return answer;
}

inline uint64_t splitBy2(unsigned int a){
	uint64_t x = a;
	x = (x | x << 16) & 0xFFFF0000FFFF;  // shift left 16 bits, OR with self, and 0000000000000000111111111111111100000000000000001111111111111111
	x = (x | x << 8) & 0xFF00FF00FF00FF;  // shift left 8 bits, OR with self, and 0000000011111111000000001111111100000000111111110000000011111111
	x = (x | x << 4) & 0xF0F0F0F0F0F0F0F; // shift left 4 bits, OR with self, and 0000111100001111000011110000111100001111000011110000111100001111
	x = (x | x << 2) & 0x3333333333333333; // shift left 2 bits, OR with self, and 0011001100110011001100110011001100110011001100110011001100110011
	x = (x | x << 1) & 0x5555555555555555; // shift left 1 bits, OR with self, and 0101010101010101010101010101010101010101010101010101010101010101
	return x;
}

inline uint64_t mortonEncode_magicbits(unsigned int x, unsigned int y){
	uint64_t answer = 0;
	answer |= splitBy2(x) | splitBy2(y) << 1;
	return answer;
}



inline uint64_t keyXY(uint64_t x, uint64_t y, int8_t max_level){
	uint64_t answer = 0;
	answer |= x | (y << max_level);
	return answer;
}

inline uint64_t keyXYZ(uint64_t x, uint64_t y, uint64_t z, int8_t max_level){
	uint64_t answer = 0;
	answer |= x | (y << max_level) | (z << 2*max_level);
	return answer;
}


#endif /* INLINEDFUNCT_HPP_ */
