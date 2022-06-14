/*

	C++ implementation of a hashing function based on new operations in the symmetric group.
	Author: Adrian Perez Keilty
	
	Other files:	hash_sn.h
					main.cpp

*/

#include "hash_sn.h"

using namespace std;

boost_string hash_sn(const string str)
{
	long long unsigned int start, finish, diff; 


	// Uncomment for performance stats
	
	//chrono::steady_clock::time_point begin = chrono::steady_clock::now();
	//current_stamp(&start);
	
    HASH_SN hash=HASH_SN(str);
	
	//current_stamp(&finish);
	//diff = finish - start;
	
	//printf("Number of cycles: %lld cycles.\n", diff);
    //chrono::steady_clock::time_point end = chrono::steady_clock::now();

    //cout << "Ellapsed time: " << chrono::duration_cast<chrono::microseconds>(end - begin).count()/1000000.0 << "\n";

    return hash.output_hex;
}

// Constructor
HASH_SN::HASH_SN(string str)
{
	hash_blocks(str);
}

void HASH_SN::hash_blocks(string M){
	
	// Padding
	len_orig=M.length();
	len_pad=half_block_size-(len_orig%half_block_size);
	
	if (!(len_pad%block_size)){
		len_pad+=half_block_size;
	};
	
	// Append bits 10000000
	M.push_back(0x80);
	
	// Append zeros until we reach the next half block
	for(i=1;i<len_pad;i++){
		M.push_back(0);
	};
	
	k,k_aux=0;
	count=count_half_block_size;
	len_min=(len_orig>half_block_size) ? half_block_size : len_orig;
	
	// Append reverse last half block 
	for(i=0;i<len_min;i++){
		
		M_char=M[len_orig-1-i]+0;
		M.push_back(M_char);
		k_aux+=M_char;
		k+=k_aux<<count;
		k_aux=0;
		count-=8;
		
	};
	
	// Reverse last block*2^64
	k=k<<64;
	count=count_block_size;
		
	// Initial value -> (Reverse last block*2^64)<<len_orig* 2^64
	int_to_per(fact,b,static_cast<uint132>(len_orig)<<64,e,t);
	int_to_per(fact,a,k,e,t);
	
	for(h=0;h<t;h++){
		a_index[a[h]]=h;
	};
	
	a_rotate_b_inv(a,b,buf_a_rotate_b,t);
	
	// First character, loop case i=0
	k=M[0];
	n+=(k<<count);
	count-=8;
	
	// Initialize accumulative block product variable
	k_last=1;
	
	len_new=len_orig+len_pad+half_block_size;
	
	// Process message blocks
	for(i=1;i!=len_new;i++){
		
		if(!(i%block_size)){
			
			k_last*=n;
			k_last%=p;
			
			int_to_per(fact,b,n,e,t);
			
			if(buf_a_rotate_b[(i-1)%t]%2){
				a_rotate_b_inv(buf_a_rotate_b,b,buf_a_rotate_b,t);
			}
			else {
				a_rotate_b(buf_a_rotate_b,b,buf_a_rotate_b,t);
			}
			
			n=0;
			count=count_block_size;
		}
		k=M[i];	
		n+=(k<<count);
		count-=8;
	};
	
	k_last*=n;
	k_last%=p;
	int_to_per(fact,b,n,e,t);
	if(buf_a_rotate_b[(i-1)%t]%2){
		a_rotate_b_inv(buf_a_rotate_b,b,buf_a_rotate_b,t);
	}
	else {
		a_rotate_b(buf_a_rotate_b,b,buf_a_rotate_b,t);
	}
	
	// Last state: k_last << buf_a_rotate_b
	int_to_per(fact,a,static_cast<uint132>(k_last),e,t);
	a_rotate_b_inv(a,buf_a_rotate_b,buf_a_rotate_b,t);
	
	//Hash resulting from computing last state of a>>b modulo 2^128
	output_hex=boost::format("%1$032x") %(per_to_int(buf_a_rotate_b,e,t)%two_pow_128);
	
}

// a>>b for a and b in symmetric group
// Stores value in buf_a_rotate_b
void HASH_SN::a_rotate_b(uint8_t a[],uint8_t b[],uint8_t buf_a_rotate_b[],uint8_t t){
	
    uint8_t n,r;
	
    for(h=0;h<t;h++){
		
        n=a[h];
        r=b[h];
		
        if(n!=0){
			
            for(l=h;l<(n+h);l++){
                b[l%t]=b[(l+1)%t];
            }
            b[(h+n)%t]=r;
        }
        else {
			
            for(l=0;l<t;l++){
                b_aux[l]=b[l];
            };
            for(l=0;l<t;l++){
				b[l]=b_aux[a_index[t-1-a[l]]];
            }
        }

    };
	for(h=0;h<t;h++){
		val=b[h];
		buf_a_rotate_b[h]=val;
		a_index[val]=h;
	};
};

// b<<a for a and b in symmetric group
// Stores value in buf_a_rotate_b
void HASH_SN::a_rotate_b_inv(uint8_t a[],uint8_t b[],uint8_t buf_a_rotate_b[],uint8_t t){
	
    uint8_t n,r,h1;
	
    for(h=t-1;h>-1;h--){
		
        n=a[h];
        r=b[(h+n)%t];
		
        if(n!=0){
			
            for(l=n+h;l>h;l--){
				
                b[l%t]=b[(l-1+t)%t];
            
			}
            b[h]=r;
			
        }
		// Assuming t is odd
        else {
			
            for(l=0;l<t;l++){
				
                b_aux[l]=b[l];
            };
            for(l=0;l<t;l++){
				b[l]=b_aux[a_index[t-1-a[l]]];
            };
        };
    };
	for(h=0;h<t;h++){
		val=b[h];
		buf_a_rotate_b[h]=val;
		a_index[val]=h;
	};
};

void HASH_SN::int_to_per(uint132 fact[],uint8_t a[],uint132 b, vector<uint8_t> e, uint8_t t){
	
	// To factoradic
	fact[t-1]=0;
	for(j=1;j!=t;j++){
		
		fact[t-1-j]=b%(j+1);
		b/=j+1;
		
	};
	
	// To permutation
	for(j=0;j!=t;j++){
		
		fact_j=static_cast<uint8_t>(fact[j]);
		a[j]=e.at(fact_j);
		e.erase(e.begin()+fact_j);
		
	};
};

// Convert permutation to integer
uint132 HASH_SN::per_to_int(uint8_t a[],vector<uint8_t> e, uint8_t t){
	
	uint8_t index,e_index,a0;
	uint132 fact_count=fact_34;
	
    // To factorial system
	a0=a[0];
	e.erase(e.begin()+a0);
	uint132 k_per_to_int=a0*fact_count;
	
	for(uint iter=1;iter<t;iter++){
		fact_count/=t-iter;
		index=index_vector(e,a[iter]);
		e.erase(e.begin()+index);
		k_per_to_int+=index*fact_count;
	}
    return k_per_to_int;
};

// Return index of element in array
uint8_t HASH_SN::index_array(uint8_t a[], uint8_t n,uint8_t size){
	
	count_index_array=0;
    uint8_t* p=a;
	
    for ( auto it = p; count_index_array != size; ++it){
		
        if(*it == n) {
            break;
        };
		
        count_index_array+=1;
        p++;
    }
    return count_index_array;
};


// Return index of element in vector
uint8_t HASH_SN::index_vector(vector<uint8_t> v, uint8_t n){
	
    vector<uint8_t>::iterator it;
    count_index_vector=0;
	
    for(it = v.begin(); it < v.end(); it++,count_index_vector++){
		
    if(*it == n) {
        break;
        };
    };
    return count_index_vector;
};


// Print content of array
void HASH_SN::print_arr(uint8_t a[],uint8_t size){
	
    for(uint iter=0;iter<size;iter++){
        std::cout << unsigned(a[iter]) << " " ;
    }
    std::cout << "\n";
};
