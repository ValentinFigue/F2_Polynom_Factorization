//
//  Arithmetic.cpp
//  F2_Polynom_Factorization
//
//  Created by Valentin Figué  on 25/12/2019.
//

#include "Arithmetic.hpp"

F2_polynom F2_Arithmetic::compute_greatest_common_divisor(const F2_polynom &i_polynom_left, const F2_polynom &i_polynom_right){
    /*
     Function which computes the greatest common divisor between two polynoms using Euclidean algorithm.
     */

    // Initialization of variables
    F2_polynom dividende_polynom;
    F2_polynom diviseur_polynom;
    F2_polynom remainder;
    if (i_polynom_left.max_order>i_polynom_right.max_order){
        dividende_polynom = i_polynom_left;
        diviseur_polynom = i_polynom_right;
    }
    else {
        dividende_polynom = i_polynom_right;
        diviseur_polynom = i_polynom_left;
    }
    remainder = diviseur_polynom;

    std::vector<F2_polynom> division_result;

    bool start = false;

    while(remainder.max_order>-1){
        if (start){
            dividende_polynom = diviseur_polynom;
            diviseur_polynom = remainder;

        }
        start = true;
        division_result = F2_polynom::division(dividende_polynom,diviseur_polynom);
        remainder = division_result[1];
    }

    if (diviseur_polynom.max_order == -1){
        std::vector<int> coeff;
        coeff.push_back(1);
        diviseur_polynom = F2_polynom(coeff);
    }
    return diviseur_polynom;

}


std::vector<F2_polynom> F2_Arithmetic::compute_square_free_factorization(const F2_polynom &i_polynom){
    /*
     Function which computes the square free factorization of a polynom.
     */

    // Initialization of variables
    std::vector<F2_polynom> R, recursive_R;
    F2_polynom common_divisor_derivation, square_root, factor;
    std::vector<F2_polynom> division_result;

    // Derivation
    F2_polynom derivation = F2_polynom::derivation(i_polynom);

    // If the derivation is not null
    if (derivation.max_order>-1){

        // Computation of the greatest common divisor between the polynom and the derivation
        common_divisor_derivation = compute_greatest_common_divisor(i_polynom,derivation);

        // Check that there the derivation divide the polynom
        if (common_divisor_derivation.max_order>0){
            division_result = F2_polynom::division(i_polynom,common_divisor_derivation);
            factor = division_result[0];

            // Recursive call over the two factors
            recursive_R = compute_square_free_factorization(common_divisor_derivation);
            for (int i = 0; i <recursive_R.size();i++){
                R.push_back(recursive_R[i]);
            }
            recursive_R = compute_square_free_factorization(factor);
            for (int i = 0; i <recursive_R.size();i++){
                R.push_back(recursive_R[i]);
            }
        }

        else{
            R.push_back(i_polynom);
        }

    }

    // If the derivation is null
    else {
        // Computation of the square root
        square_root = F2_polynom::square_root(i_polynom);
        // Recursion over each element of the square root
        recursive_R = compute_square_free_factorization(square_root);
        for (int i = 0; i <recursive_R.size();i++){
            for (int j = 0; j<2;j++){
                R.push_back(recursive_R[i]);
            }
        }
    }

    return R;
}




std::vector<F2_polynom> F2_Arithmetic::compute_Berlekamp_matrix(const std::vector<F2_polynom> &i_vectors){
    /*
     Function which computes the Berlekamp matrix which can be used to compute polynom factorization.
     */

    std::vector<F2_polynom> matrix;

    // Loop over input vectors
    for (int i = 0; i<i_vectors.size();i++){

        std::vector<int> coeff_vector;
        for (int j = 0;j<i_vectors.size();j++){
            coeff_vector.push_back(0);
        }
        for (int j = 0;j<i_vectors[i].max_order+1;j++){
            coeff_vector[j] =i_vectors[i].coefficients[j] ;
        }
        if (coeff_vector[i]==1){
            coeff_vector[i] = 0;
        }
        else{
            coeff_vector[i] = 1;
        }

        std::vector<int> coeff_identity;
        for (int j = 0;j<i_vectors.size();j++){
            coeff_identity.push_back(0);
        }
        coeff_identity[i]= 1;

        std::vector<int> new_coeff;
        for (int j = 0;j<2*i_vectors.size();j++){
            new_coeff.push_back(0);
            if (j<i_vectors.size()){
                new_coeff[j] = coeff_identity[j];
            }
            else{
                new_coeff[j] = coeff_vector[j-i_vectors.size()];
            }
        }

        matrix.push_back(F2_polynom(new_coeff));
    }
    return matrix;
}

std::vector<F2_polynom> F2_Arithmetic::resolve_Berlekamp_matrix(const std::vector<F2_polynom> &i_matrix){
    /*
     Function which resolve the Berlekamp matrix which can be used to compute polynom factorization.
     */

    std::vector<F2_polynom> matrix;
    F2_polynom temp_row;
    matrix = i_matrix;
    int row_indicator = 0;

    while(row_indicator<matrix.size()){

        // Find the left most row
        int left_most_row = row_indicator;
        for (int i = row_indicator;i<matrix.size();i++){
            if (matrix[i].max_order>matrix[left_most_row].max_order){
                left_most_row  = i;
            }
        }

        // Swap left most row and current row
        temp_row = matrix[row_indicator];
        matrix[row_indicator] = matrix[left_most_row];
        matrix[left_most_row] = temp_row;

        // Add current row to each row containaing the actual left_bits_raw
        for (int i = 0;i<matrix.size();i++){
            if (i != row_indicator){
                if (matrix[i].coefficients[matrix[row_indicator].max_order] == 1){
                    matrix[i] = F2_polynom::addition(matrix[i],matrix[row_indicator]);
                }
            }
        }

        row_indicator++;
    }

    return matrix;
}


std::vector<F2_polynom> F2_Arithmetic::find_null_vectors(const std::vector<F2_polynom> &i_matrix){
    /*
     Function which finds the null vectors of the matrix.
     */

    std::vector<F2_polynom> null_vectors;

    for (int i = 0; i<i_matrix.size(); i++){
        std::vector<int> left_coeff, right_coeff;
        F2_polynom left_polynom, right_polynom;
        for (int j = 0; j<i_matrix.size(); j++){
            left_coeff.push_back(0);
            if ((i_matrix.size()+j)<i_matrix[i].max_order+1){
                left_coeff[j] = i_matrix[i].coefficients[i_matrix.size()+j];
            }
            right_coeff.push_back(0);
            if (j<i_matrix[i].max_order+1){
                right_coeff[j] = i_matrix[i].coefficients[j];
            }
        }
        left_polynom = F2_polynom(left_coeff);
        right_polynom = F2_polynom(right_coeff);
        if (left_polynom.max_order == -1){
            if (right_polynom.max_order >0){
                null_vectors.push_back(right_polynom);
            }
        }
    }

    return null_vectors;
}

std::vector<F2_polynom> F2_Arithmetic::compute_Berlekamp_factorization(const F2_polynom &i_polynom){
    /*
     Function which computes the Berlekamp factorization of a square free polynom.
     */

    // Initialization of variables
    std::vector<F2_polynom> division_result, vector_result, matrix, null_vectors, factorization, recursive_factorization;
    F2_polynom vector, factor_1, factor_2;

    // Computation of new vectors
    for (int i = 0; i<i_polynom.max_order;i++){
        std::vector<int> coeff;
        for (int j = 0; j<2*i;j++){
            coeff.push_back(0);
        }
        coeff.push_back(1);
        vector = F2_polynom(coeff);
        division_result = F2_polynom::division(vector,i_polynom);
        vector_result.push_back(division_result[1]);
    }

    // Computation of matrix
    matrix = compute_Berlekamp_matrix(vector_result);

    // Resolution of the matrix
    matrix = resolve_Berlekamp_matrix(matrix);

    // Resolution of the matrix
    null_vectors = find_null_vectors(matrix);
    
    std::vector<int> coeff;
    coeff.push_back(1);
    F2_polynom one_polynom = F2_polynom(coeff);

    if (null_vectors.size()>0){
        factor_1 = compute_greatest_common_divisor(i_polynom,null_vectors[0]);
        recursive_factorization = compute_Berlekamp_factorization(factor_1);
        for (int i = 0; i<recursive_factorization.size();i++){
            factorization.push_back(recursive_factorization[i]);
        }
        division_result = F2_polynom::division(i_polynom,factor_1);
        factor_2 = division_result[0];
        recursive_factorization = compute_Berlekamp_factorization(factor_2);
        for (int i = 0; i<recursive_factorization.size();i++){
            factorization.push_back(recursive_factorization[i]);
        }
    }
    else {
        factorization.push_back(i_polynom);
    }

    return factorization;
}


std::vector<F2_polynom> F2_Arithmetic::polynom_factorization(const F2_polynom &i_polynom){
    /*
     Function which computes the factorization of a polynom.
     */

    std::vector<F2_polynom> factorization, square_free_factorization, recursive_factorization;

    square_free_factorization = compute_square_free_factorization(i_polynom);

    for (int i = 0; i <square_free_factorization.size();i++){
        recursive_factorization = compute_Berlekamp_factorization(square_free_factorization[i]);
        for (int j = 0; j<recursive_factorization.size();j++){
            factorization.push_back(recursive_factorization[j]);
        }
    }

    return factorization;

}
