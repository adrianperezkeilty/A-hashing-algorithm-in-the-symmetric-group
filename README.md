# Hash_sn

- C++ implementation of a hashing function based on a one way function in the symmetric group.
  
  The Boost library is needed for the compilation (https://www.boost.org/).
  
  Example on windows command line:
  
  Compilation: g++ main.cpp hash_sn.cpp -I"include" -o hash_sn.exe
  
  Execution: hash_sn.exe digest_message_1MB.txt

- Python implementation of hash_sn + empirical tests.
  
  hash_sn.py                 -> hashing algorithm + tests,
  
  two_block_attack_hash_sn   -> Simulate 2-block attack for small sizes (16-bit, 24-bit etc),
  
  finite_field_arithmetic.py -> compute inverses modulo a prime using Fermat's little theorem.

