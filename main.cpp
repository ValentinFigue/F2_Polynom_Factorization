#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <F2_Polynom.hpp>
#include <Arithmetic.hpp>
#include <utils.hpp>

int main()
{

    // Inputs reading
    int size;
    std::cin >> size;

    unsigned int* encoded = new unsigned int[size / 16];
    for (int i = 0; i < size / 16; i++) {   // Read size / 16 integers to a
            std::cin >> std::hex >> encoded[i];
    }

    // Transformation to polynoms
    F2_polynom encoded_polynom = convert_integers_to_polynoms(encoded,size);

    // Factorisation of encoded polynoms
    std::vector<F2_polynom> factorization = F2_Arithmetic::polynom_factorization(encoded_polynom);
    
    // Transformation to integers
    std::vector< std::vector<unsigned int> > factorization_integers = convert_polynoms_to_integers(factorization, size);

    // Hexadecimal transformation
    std::set<std::string> final_result;

    for (int i = 0;i < factorization_integers.size(); i++){
        std::stringstream sstream;
        for (int j = 0;j < factorization_integers[i].size();j++){
            if (j < (factorization_integers[i].size()-1)){
                sstream << std::setfill('0') << std::setw(8) << std::hex << factorization_integers[i][j]<<" ";
            }
            else{
                sstream << std::setfill('0') << std::setw(8) << std::hex << factorization_integers[i][j];
            }
        }
        final_result.insert(sstream.str());
    }

    // Alphabetical sort
    std::set<std::string>::iterator it;
    for (it = final_result.begin(); it != final_result.end(); ++it){
        std::cout<<*it<<std::endl;
    }

    return 0;
}
