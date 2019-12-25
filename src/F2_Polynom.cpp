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
    int max_order = std::max(i_polynom_left.max_order,i_polynom_right.max_order);

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

F2_polynom F2_polynom::substraction(const F2_polynom &i_polynom_left,const F2_polynom &i_polynom_right){
    /*
     Definition of the substraction between two polynoms
     */


    // Computation of the new max order
    int max_order = std::max(i_polynom_left.max_order,i_polynom_right.max_order);

    // Initialization of the resulting vector
    std::vector<int> coefficients;
    for (int i = 0; i<max_order+1;i++){
       coefficients.push_back(0);
    }

    if (i_polynom_left.max_order<max_order){
        for (int i = 0; i<i_polynom_left.max_order+1;i++){
            coefficients[i] = i_polynom_left.coefficients[i]-i_polynom_right.coefficients[i];
        }
        for (int i = i_polynom_left.max_order+1; i<max_order+1;i++){
            coefficients[i] = -i_polynom_right.coefficients[i];
        }
    }
    else{
        for (int i = 0; i<i_polynom_right.max_order+1;i++){
            coefficients[i] = i_polynom_left.coefficients[i]-i_polynom_right.coefficients[i];
        }
        for (int i = i_polynom_right.max_order+1; i<max_order+1;i++){
            coefficients[i] = i_polynom_left.coefficients[i];
        }
    }

    for (int i = 0; i<max_order+1;i++){
        if (coefficients[i] == -1){
            coefficients[i] = 1;
        }
    }

    return F2_polynom(coefficients);

}

F2_polynom F2_polynom::multiplication(const F2_polynom &i_polynom_left,const F2_polynom &i_polynom_right){
    /*
     Definition of the multiplication between two polynoms
     */

    // Computation of the new max order
    int max_order = i_polynom_left.max_order+i_polynom_right.max_order;

    // Initialization of the resulting vector
    std::vector<int> coefficients;
    for (int i = 0; i<max_order+1;i++){
       coefficients.push_back(0);
    }

    for (int i = 0; i<i_polynom_left.max_order+1;i++){
        for (int j =0; j<i_polynom_right.max_order+1;j++){
           coefficients[i+j] ^= i_polynom_left.coefficients[i]&i_polynom_right.coefficients[j];
        }
    }

    return F2_polynom(coefficients);
}

std::vector<F2_polynom> F2_polynom::division(const F2_polynom &i_dividende_polynom, const F2_polynom &i_diviseur_polynom){
    /*
     Definition of the division between two polynoms
     */

    // Initialization of the quotient polynom

    // Computation of the quotient max order
    int quotient_order = i_dividende_polynom.max_order-i_diviseur_polynom.max_order;

    // Initialization of the resulting vector
    std::vector<int> quotient_coefficients;
    for (int i = 0; i<quotient_order+1;i++){
       quotient_coefficients.push_back(0);
    }

    F2_polynom quotient_polynom = F2_polynom(quotient_coefficients);

    // Initialization of the rest polynom
    F2_polynom rest_polynom = i_dividende_polynom;

    F2_polynom s;

    // Division
    while((rest_polynom.max_order>=i_diviseur_polynom.max_order) && (i_diviseur_polynom.max_order > -1)){
        quotient_coefficients[rest_polynom.max_order-i_diviseur_polynom.max_order] = 1;
        s = F2_polynom(quotient_coefficients);
        quotient_coefficients[rest_polynom.max_order-i_diviseur_polynom.max_order] = 0;

        quotient_polynom = addition(quotient_polynom,s);
        rest_polynom = substraction(rest_polynom,multiplication(i_diviseur_polynom,s));
    }


    std::vector<F2_polynom> result;
    result.push_back(quotient_polynom);
    result.push_back(rest_polynom);

    return result;
}

F2_polynom F2_polynom::compute_derivation(const F2_polynom &i_polynom){
    /*
     Function which computes the derivate of a polynom.
     */

    // Initialization of the derivativate of the input polynom
    std::vector<int> derivative_coeffs;

    // Computation of the derivation in F2
    for (int i = 0; i<i_polynom.max_order;i++){
        if ((i%2)==0){
            derivative_coeffs.push_back(i_polynom.coefficients[i+1]);
        }
        else{
            derivative_coeffs.push_back(0);
        }
    }

    return F2_polynom(derivative_coeffs);

}

F2_polynom F2_polynom::compute_square_root(const F2_polynom &i_polynom){
    /*
     Function which computes the square root of a polynom.
     */

    // Initialization of the square root the input polynom
    std::vector<int> square_coeffs;
    for (int i = 0; i<(i_polynom.max_order/2)+1;i++){
        square_coeffs.push_back(0);
    }

    // Computation of the derivation in F2
    for (int i = 0; i<(i_polynom.max_order/2)+1;i++){
        square_coeffs[i] =  i_polynom.coefficients[2*i];
    }

    return F2_polynom(square_coeffs);

}
