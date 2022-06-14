/*

	C++ implementation of a hashing function based on new operations in the symmetric group.
	Author: Adrian Perez Keilty
	
	Other files:	hash_sn.cpp
					main.cpp

*/

#include <iostream>
#include <cstring>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/format.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <stdio.h>
#include <chrono>
#include <fstream>
#include <sstream>

#define current_stamp(a) asm volatile("rdtsc" : "=a"(((unsigned int *)(a))[0]),"=d"(((unsigned int *)a)[1]))

using namespace boost::multiprecision;

typedef number<cpp_int_backend<132, 132, unsigned_magnitude, unchecked, void> >   uint132;
typedef number<cpp_int_backend<264, 264, unsigned_magnitude, unchecked, void> >   uint264;

//typedef uint256_t uint256;
typedef unsigned int uint;
typedef boost::basic_format<char> boost_string;

using namespace std;

typedef basic_string<unsigned char> ustring;

class HASH_SN
{
	public:
	
	HASH_SN(string str);
	boost_string output_hex;
	
	private:
	
	uint132 n=0;
	uint132 k,k_aux;
	uint264 k_last;
	
	uint8_t fact_j,j;
	uint i;
	int8_t h,l,val;
	
	static inline uint132 p1=1;
	static inline uint132 p2=(p1<<131)-347; 
	static inline uint132 p=p2+(p1<<131); // p -> prime
	
	inline const unsigned char *uc_str(const char *s){
	return reinterpret_cast<const unsigned char*>(s);
	}
	
	const static uint8_t block_size=16;
	const static uint8_t half_block_size=block_size/2;
	const static uint8_t count_block_size=block_size*8-8;
	const static uint8_t count_half_block_size=block_size*4-8;
	const static uint8_t t=35; // t! smallest factorial larger than 2^s
	vector<uint8_t> e={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
	18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};
	
	uint len_pad,len_new,len_orig,len_min,count;
	uint8_t count_index_vector,count_index_array;
	char M_char;
	
	uint132 fact[t];
	uint8_t a[t],b[t],b_aux[t],buf_a_rotate_b[t],a_index[t];
	
	void hash_blocks(string str);
	void a_rotate_b(uint8_t*,uint8_t*,uint8_t*,uint8_t);
	void a_rotate_b_inv(uint8_t*,uint8_t*,uint8_t*,uint8_t);
	uint8_t index_array(uint8_t*,uint8_t,uint8_t);
	uint8_t index_vector(vector<uint8_t> , uint8_t);
	uint132 per_to_int(uint8_t*,vector<uint8_t>, uint8_t);
	void int_to_per(uint132*,uint8_t*,uint132,vector<uint8_t>,uint8_t);
	void print_arr(uint8_t*,uint8_t); 
	
	uint132 two_pow_128=boost::multiprecision::pow(uint132(2),128);
	
	uint132 fact_20=2432902008176640000; // Hignest factorial lower than 2^64;
	
	uint132 fact_34=fact_20*21*22*23*24*25*26*27*28*29*30*31*32*33*34;
};

boost_string hash_sn(const string str);

