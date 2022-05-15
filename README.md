# hash_sn

C++ implementation of a hashing function based on new operations in the symmetric group.

The Boost library is needed for the compilation. 

For windows command line:
Example of compilation using minGW and storing the folder "boost" in "include":
g++ main.cpp hash_sn.cpp -I"include" -o hash_sn.exe

Then execute 
hash_sn.exe digest_message_1MB.txt


