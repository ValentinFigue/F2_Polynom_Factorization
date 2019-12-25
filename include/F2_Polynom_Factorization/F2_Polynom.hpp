//
//  F2_Polynom.hpp
//  F2_Polynom_Factorization
//
//  Created by Valentin Figué  on 25/12/2019.
//

#ifndef F2_Polynom_hpp
#define F2_Polynom_hpp

#include <stdio.h>
#include <vector>

class F2_polynom{
public:


    // Data Members
    int max_order; // max order of the polynom, must be a multiple of 32
    std::vector<int> coefficients;  // vector which represents the coefficients of the polynom

    F2_polynom();
    F2_polynom(const std::vector<int> &i_coefficients);


    const void print();
    
    static F2_polynom addition(const F2_polynom &i_polynom_left,const F2_polynom &i_polynom_right);
    static F2_polynom substraction(const F2_polynom &i_polynom_left,const F2_polynom &i_polynom_right);
    static F2_polynom multiplication(const F2_polynom &i_polynom_left,const F2_polynom &i_polynom_right);
    static std::vector<F2_polynom> division(const F2_polynom &i_dividende_polynom, const F2_polynom &i_diviseur_polynom);
    static F2_polynom derivation(const F2_polynom &i_polynom);
    static F2_polynom square_root(const F2_polynom &i_polynom);
    
};
    
#endif /* F2_Polynom_hpp */
