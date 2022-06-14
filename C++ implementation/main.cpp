/*

	C++ implementation of a hashing function based on new operations in the symmetric group.
	Author: Adrian Perez Keilty
	
	Other files:	hash_sn.cpp
					hash_sn.h

*/

#include "hash_sn.h"

using namespace std;


// Convert File into string
string readFileIntoString(const string& path) {
    
    ifstream input_file(path);
	
    if (!input_file.is_open()) {
        cerr << "Could not open the file - '"
             << path << "'" << endl;
        exit(EXIT_FAILURE);
    }
	
    return string((istreambuf_iterator<char>(input_file)), istreambuf_iterator<char>());
}

int main(int argc, char** argv){

    string file_contents;
    string file_name=argv[1];
	
    cout << "Reading file..." << endl;
    file_contents = readFileIntoString(file_name);

    cout << "Digesting message..." << endl;
	
    cout << "Hash of " << file_name << ":" << hash_sn(file_contents) << endl;
	    
    return 0;	
}
