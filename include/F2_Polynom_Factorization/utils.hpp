//
//  utils.hpp
//  F2_Polynom_Factorization
//
//  Created by Valentin Figué  on 25/12/2019.
//

#ifndef utils_hpp
#define utils_hpp

#include <stdio.h>
#include <vector>
#include "F2_Polynom.hpp"

F2_polynom convert_integers_to_polynoms(unsigned int  *integers_array, const int size);
std::vector < std::vector <unsigned int> > convert_polynoms_to_integers(const std::vector < F2_polynom > &polynoms_vector, int size);

#endif /* utils_hpp */
