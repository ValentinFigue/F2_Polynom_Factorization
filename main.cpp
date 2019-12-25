#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <set>

using namespace std;

F2_polynom compute_greatest_common_divisor(const F2_polynom &i_polynom_left, const F2_polynom &i_polynom_right){
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
        division_result = division(dividende_polynom,diviseur_polynom);
        remainder = division_result[1];
    }

    if (diviseur_polynom.max_order == -1){
        std::vector<int> coeff = {1};
        diviseur_polynom = F2_polynom(coeff);
    }
    return diviseur_polynom;

}


std::vector<F2_polynom> compute_square_free_factorization(const F2_polynom &i_polynom){
    /*
     Function which computes the square free factorization of a polynom.
     */

    // Initialization of variables
    std::vector<F2_polynom> R, recursive_R;
    F2_polynom common_divisor_derivation, square_root, factor;
    std::vector<F2_polynom> division_result;

    // Derivation
    F2_polynom derivation = compute_derivation(i_polynom);

    // If the derivation is not null
    if (derivation.max_order>-1){

        // Computation of the greatest common divisor between the polynom and the derivation
        common_divisor_derivation = compute_greatest_common_divisor(i_polynom,derivation);

        // Check that there the derivation divide the polynom
        if (common_divisor_derivation.max_order>0){
            division_result = division(i_polynom,common_divisor_derivation);
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
        square_root = compute_square_root(i_polynom);
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


std::vector<F2_polynom> compute_Berlekamp_matrix(const std::vector<F2_polynom> &i_vectors){
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


std::vector<F2_polynom> resolve_Berlekamp_matrix(const std::vector<F2_polynom> &i_matrix){
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
                    matrix[i] = addition(matrix[i],matrix[row_indicator]);
                }
            }
        }

        row_indicator++;
    }

    return matrix;
}


std::vector<F2_polynom> find_null_vectors(const std::vector<F2_polynom> &i_matrix){
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


std::vector<F2_polynom> compute_Berlekamp_factorization(const F2_polynom &i_polynom){
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
        division_result = division(vector,i_polynom);
        vector_result.push_back(division_result[1]);
    }

    // Computation of matrix
    matrix = compute_Berlekamp_matrix(vector_result);

    // Resolution of the matrix
    matrix = resolve_Berlekamp_matrix(matrix);

    // Resolution of the matrix
    null_vectors = find_null_vectors(matrix);

    F2_polynom one_polynom = F2_polynom(std::vector<int>({1}));

    if (null_vectors.size()>0){
        factor_1 = compute_greatest_common_divisor(i_polynom,null_vectors[0]);
        recursive_factorization = compute_Berlekamp_factorization(factor_1);
        for (int i = 0; i<recursive_factorization.size();i++){
            factorization.push_back(recursive_factorization[i]);
        }
        division_result = division(i_polynom,factor_1);
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


std::vector<F2_polynom> polynom_factorization(const F2_polynom &i_polynom){
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


std::vector<std::vector<F2_polynom>> compute_solution_from_factorization(const std::vector<F2_polynom> &factorization, int size){
    /*
     Function which computes the different polynom solutions according to the factorization of a polynom.
     */

    std::vector<std::vector<F2_polynom>> result;

    std::vector<std::vector<int>> possible_solutions;
    possible_solutions.push_back(std::vector<int>());

    std::vector<int> possible_vector;
    std::vector<int> temp_vector;

    // Compute all the different solutions
    for (int i = 0; i< factorization.size();i++){
        std::vector<std::vector<int>> temp_solutions;
        for (int j = 0; j <possible_solutions.size();j++){
            possible_vector = possible_solutions[j];
            temp_vector = possible_vector;
            temp_vector.push_back(0);
            temp_solutions.push_back(temp_vector);
            temp_vector = possible_vector;
            temp_vector.push_back(1);
            temp_solutions.push_back(temp_vector);
        }
        possible_solutions = temp_solutions;
    }

    // Compute for each solution the polynoms
    for (int i = 0; i< possible_solutions.size();i++){
        F2_polynom left = F2_polynom(std::vector<int>({1}));
        F2_polynom right = F2_polynom(std::vector<int>({1}));

        for (int j = 0; j<factorization.size();j++){
            if (possible_solutions[i][j] == 0){
                left = multiplication(left,factorization[j]);
            }
            else {
                right = multiplication(right,factorization[j]);
            }
        }

        if ((left.max_order<size)&&(right.max_order<size)){
            std::vector<F2_polynom> polynom_vector;

            polynom_vector.push_back(left);
            polynom_vector.push_back(right);

            result.push_back(polynom_vector);
        }
    }

    return result;
}


int main()
{

     /*
     Test of creation of F2 polynoms,addition, substraction, multiplication and division
     */


//    std::vector<int> coefficient_left = {1,1};
//    std::vector<int> coefficient_right = {0,1,1,0,1,0};
//    std::vector<int> coefficient = {1,0,0,1};

//    F2_polynom polynom_left = F2_polynom(coefficient_left);
//    F2_polynom polynom_right = F2_polynom(coefficient_right);
//    F2_polynom polynom_middle = F2_polynom(coefficient);


//    polynom_left.print();
//    polynom_right.print();

//    F2_polynom polynom_a, polynom_b, polynom;
//    polynom = addition(polynom_left,polynom_right);
//    polynom.print();

//    polynom = substraction(polynom_left,polynom_right);
//    polynom.print();

//    polynom_a = multiplication(polynom_left,polynom_right);
//    polynom_b = multiplication(polynom_right,polynom_middle);

//    polynom.print();

//    std::vector<F2_polynom> result;

//    result = division(polynom_a,polynom_right);
//    result[0].print();
//    result[1].print();

//    result;
//    result = division(polynom_b,polynom_left);
//    result[0].print();
//    result[1].print();

//    polynom = compute_greatest_common_divisor(polynom_a,polynom_b);
//    polynom.print();

//    polynom = compute_greatest_common_divisor(polynom_right,polynom_left);
//    polynom.print();

//    polynom_a.print();
//    polynom = multiplication(polynom_left,polynom_middle);
//    polynom.print();
//    cout<<endl;

//    F2_polynom res;

//    result = compute_square_free_factorization(polynom);
//    for (int i = 0;i<result.size();i++){
//        if (i == 0){
//            res = result[i];
//        }
//        else{
//            res  = multiplication(res,result[i]);
//        }
//        result[i].print();
//    }
//    res.print();
//    cout<<endl;

//    coefficient = {1,1,1,0,0,1,1};
//    coefficient = {1,1,1,0,0,1,1};
//    polynom_a = F2_polynom(coefficient);
//    polynom_b = F2_polynom(coefficient);
//    polynom = multiplication(polynom_a,polynom_b);
//    result = polynom_factorization(polynom);
//    for (int i = 0;i<result.size();i++){
//        result[i].print();
//    }

    // 32
    //46508fb7 6677e201

    // Inputs reading
    int size;
    cin >> size;

    unsigned int* encoded = new unsigned int[size / 16];
    for (int i = 0; i < size / 16; i++) {   // Read size / 16 integers to a
            cin >> hex >> encoded[i];
    }

    // Transformation to polynoms
    F2_polynom encoded_polynom = convert_integers_to_polynoms(encoded,size);

    // Factorisation of encoded polynoms
    std::vector<F2_polynom> factorization = polynom_factorization(encoded_polynom);

    // Definition of the different function solutions and filtering of the different polynoms solutions under the problem constraint
    std::vector<std::vector<F2_polynom>> polynoms_solutions = compute_solution_from_factorization(factorization,size);

    // Conversion to integers for the couples of polynoms found
    std::vector<std::vector<unsigned int>> decoded_integers = convert_polynoms_to_integers(polynoms_solutions,size);

    // Hexadecimal transformation
    std::set<std::string> final_result;

    for (int i = 0;i < decoded_integers.size(); i++){
        std::stringstream sstream;
        for (int j = 0;j < decoded_integers[i].size();j++){
            if (j < (decoded_integers[i].size()-1)){
                sstream << setfill('0') << setw(8) << hex << decoded_integers[i][j]<<" ";
            }
            else{
                sstream << setfill('0') << setw(8) << hex << decoded_integers[i][j];
            }
        }
        final_result.insert(sstream.str());
    }

    // Alphabetical sort
    std::set<std::string>::iterator it;
    for (it = final_result.begin(); it != final_result.end(); ++it){
        cout<<*it<<endl;
    }

    return 0;
}
