# A hashing algorithm in the symmetric group

Implementations based on the MSc thesis "A hashing algorithm based on a one-way function in the symmetric group" :      
https://lnu.diva-portal.org/smash/record.jsf?pid=diva2%3A1662325&dswid=5727

Two new operations $(\gg, \ll)$ work as mirrored one way functions in $S_n$ , in the sense that:  
  - Given $b$ and $a\gg b$ or $b\ll a$ it is hard to retrieve $a$.  
  - $(a\gg b)\ll a = a\gg(b\ll a) = b$.

We develop a hashing algorithm based on these operations by encoding blocks into permutations and exploting the algebraic incompatibility of the operations with the usual product of elements in the field $\mathbb Z_p^*$ .

## C++ implementation
  
  main.cpp  
  hash_sn.cpp  
  hash_sn.h  
  
  The Boost library is needed for the compilation (https://www.boost.org/).
  
  Example on windows command line:  
  
  Compilation: 
  
  ```
  C:\path>g++ main.cpp hash_sn.cpp -I"include" -o hash_sn.exe  
  ```
  Execution:  
  
  ```
  C:\path> hash_sn.exe digest_message.txt  
  Reading file...  
  Digesting message...  
  Hash of digest_message.txt: 528a81a97baa7d4f4a08465852ce7369  
  ``` 
  
## Python implementation + empirical tests.
  
  hash_sn.py                  -> Hashing algorithm    
  hash_aux.py                 -> Auxiliary functions for hash_sn  
  hash_testing.py             -> Empirical tests on the S_n constructions (>>, <<)  
  two_block_attack_hash_sn.py -> Simulate 2-block attack for small sizes (16-bit, 24-bit etc),  
  finite_field_arithmetic.py  -> Compute inverses modulo a prime using Fermat's little theorem.  
  
### Obtain hashing value of string  
  Parameters:  
  m -> String to hash  
  s -> Block size (128, 256, 512, 1024...)  
  t -> Factorial for embedding in the symmetric group. If s=128, then t>34.  
  p -> Prime between s and t! (https://bigprimes.org/).  

  Small string digest:  
  ```
  >>> p = 9272585787985760943894005456578885141087
  >>> hash_sn.hash_sn('The quick brown fox jumps over the lazy dog',128 ,35, p).hash()
  'f308af47709e72a537c1545eb91d0e67'
  ```
  File digest:
  ```
  >>> p=1288079068764670493881163748072332651218703668359555082283935082653270651535749
  >>> hash_sn.hash_sn(open('digest_message.txt').read(), 256, 58, p).hash()
  '94d6f93bfc4b68cc941ad6ba2a6209c5f5479df07a64fc8c093674eb760dc363'
  ```
### Empirical tests  
  
