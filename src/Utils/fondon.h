#ifndef fondon_h
#define fondon_h

#include <exception>
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
using namespace std;

constexpr uint8_t k_void_nuc = 127;

std::vector<std::string> split(const std::string &s, char delim) {
        std::vector<std::string> result;
        std::stringstream ss(s);
        std::string item;
        while (getline(ss, item, delim)) 
            result.push_back(item);
        return result;
}

map<char, string> map_fasta = {
	{'F',"Phe"}, 
	{'L',"Leu"}, 
	{'S',"Ser"}, 
	{'Y',"Tyr"}, 
	{'*',"STOP"}, 
	{'C',"Cys"}, 
	{'W',"Trp"}, 
	{'P',"Pro"}, 
	{'H',"His"}, 
	{'Q',"Gln"}, 
	{'R',"Arg"}, 
	{'I',"Ile"}, 
	{'M',"Met"}, 
	{'T',"Thr"}, 
	{'N',"Asn"}, 
	{'K',"Lys"}, 
	{'V',"Val"}, 
	{'D',"Asp"}, 
	{'E',"Glu"}, 
	{'G',"Gly"}, 
	{'A',"Ala"}
};

std::map<std::string, char> k_map_3_1 = {
    {"Phe", 'F'},
    {"Leu", 'L'},
    {"Ser", 'S'},
    {"Tyr", 'Y'},
    {"STOP", '*'},
    {"Cys", 'C'},
    {"Trp", 'W'},
    {"Pro", 'P'},
    {"His", 'H'},
    {"Gln", 'Q'},
    {"Arg", 'R'},
    {"Ile", 'I'},
    {"Met", 'M'},
    {"Thr", 'T'},
    {"Asn", 'N'},
    {"Lys", 'K'},
    {"Val", 'V'},
    {"Asp", 'D'},
    {"Glu", 'E'},
    {"Gly", 'G'},
    {"Ala", 'A'}
};

static bool cvt_to_seq(const string& fasta, string& nucs) {
	    nucs.reserve(4 * fasta.length());
	    for(auto aa: fasta) {
	        if (map_fasta.count(aa)) {
	            nucs.append(map_fasta[aa] + " ");
	        } else {
	            cout << "invalid protein sequence!\n" << endl;
	            return false;
	        }
	    }
	    nucs.pop_back();
	    return true;
}

#endif
