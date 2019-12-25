//
//  F2_Polynom.cpp
//  F2_Polynom_Factorization
//
//  Created by Valentin Figué  on 25/12/2019.
//

#include <iostream>
#include <iomanip>

#include "F2_Polynom.hpp"

F2_polynom::F2_polynom(){};

F2_polynom::F2_polynom(const std::vector<int> &i_coefficients){
    /*
     Initialization function of the polynom
     */

    this->coefficients = i_coefficients;
    this->max_order = -1;
    for (int i = 0;i<i_coefficients.size();i++){
        if (i_coefficients[i] == 1){
            this->max_order = i;
        }
    }
}

void const F2_polynom::print(){
    /*
     Display function
     */

    for (int i=max_order; i>-1;i--){
        std::cout<<this->coefficients[i];
    }
    std::cout<<std::endl;
}

F2_polynom F2_polynom::addition(const F2_polynom &i_polynom_left,const F2_polynom &i_polynom_right){
    /*
     Definition of the addition between two polynoms
     */

    // Computation of the new max order
    int max_order = max(i_polynom_left.max_order,i_polynom_right.max_order);

    // Initialization of the resulting vector
    std::vector<int> coefficients;
    for (int i = 0; i<max_order+1;i++){
       coefficients.push_back(0);
    }

    if (i_polynom_left.max_order<max_order){
        for (int i = 0; i<i_polynom_left.max_order+1;i++){
            coefficients[i] = i_polynom_left.coefficients[i]^i_polynom_right.coefficients[i];
        }
        for (int i = i_polynom_left.max_order+1; i<max_order+1;i++){
            coefficients[i] = i_polynom_right.coefficients[i];
        }
    }
    else{
        for (int i = 0; i<i_polynom_right.max_order+1;i++){
            coefficients[i] = i_polynom_left.coefficients[i]^i_polynom_right.coefficients[i];
        }
        for (int i = i_polynom_right.max_order+1; i<max_order+1;i++){
            coefficients[i] = i_polynom_left.coefficients[i];
        }
    }

    return F2_polynom(coefficients);
}
