//
//  Arithmetic.hpp
//  F2_Polynom_Factorization
//
//  Created by Valentin Figué  on 25/12/2019.
//

#ifndef Arithmetic_hpp
#define Arithmetic_hpp

#include <stdio.h>
#include "F2_Polynom.hpp"

class F2_Arithmetic{
  
public:
    static F2_polynom compute_greatest_common_divisor(const F2_polynom &i_polynom_left, const F2_polynom &i_polynom_right);
    static std::vector<F2_polynom> compute_square_free_factorization(const F2_polynom &i_polynom);
    static std::vector<F2_polynom> compute_Berlekamp_matrix(const std::vector<F2_polynom> &i_vectors);
    static std::vector<F2_polynom> resolve_Berlekamp_matrix(const std::vector<F2_polynom> &i_matrix);
    static std::vector<F2_polynom> find_null_vectors(const std::vector<F2_polynom> &i_matrix);
    static std::vector<F2_polynom> compute_Berlekamp_factorization(const F2_polynom &i_polynom);
    static std::vector<F2_polynom> polynom_factorization(const F2_polynom &i_polynom);
};


#endif /* Arithmetic_hpp */
