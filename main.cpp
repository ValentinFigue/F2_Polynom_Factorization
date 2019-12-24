#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace std;


// Definition of class which represents all the polynomials whith coefficients in F2

class F2_polynom{
public:


   // Data Members
   int max_order; // max order of the polynom, must be a multiple of 32
   unsigned int* coefficients;  // vector which represents the coefficients of the polynom


   F2_polynom(){


   }


   F2_polynom(unsigned int* i_coefficients, int i_order){
       /*
        Initialization function of the polynom
        */


       this->coefficients = new unsigned int[(i_order/32)+1];
       for (int i = 0; i < (i_order/32)+1; i++){
           this->coefficients[i] = i_coefficients[i];
       }


       int indice = i_order;
       while((indice>-1)&&(((coefficients[indice/32]>>(indice%32))&1)==0)){
           indice--;
       }


       this->max_order = indice;


   }


   void print() const{
       /*
        Display function
        */


       for (int i=(max_order/32); i>-1;i--){
           cout<<this->coefficients[i];
       }
       cout<<endl;
   }






};



F2_polynom multiplication(const F2_polynom &i_polynom_left,const F2_polynom &i_polynom_right, unsigned int* temp_coeff){
   /*
    Definition of the multiplication between two polynoms
    */


   // Computation of the new max order
   int max_order = i_polynom_left.max_order+i_polynom_right.max_order;

   for (int i = 0; i<i_polynom_left.max_order+1;i++){
       for (int j =0; j<i_polynom_right.max_order+1;j++){
           temp_coeff[(i+j)/32] ^= ((i_polynom_left.coefficients[i/32]>>(i%32))&(i_polynom_right.coefficients[j/32]>>(j%32))&1)<<((i+j));
       }
   }


   F2_polynom result =  F2_polynom(temp_coeff, max_order);

   for (int i = 0; i<(max_order/32)+1;i++){
       temp_coeff[i] = 0;
   }

   return result;

}

std::vector<F2_polynom> division(const F2_polynom &i_dividende_polynom, const F2_polynom &i_diviseur_polynom, unsigned int *quotient_coefficients){
   /*
    Definition of the division between two polynoms
    */

   // Initialization of the quotient polynom

   // Computation of the quotient max order
   int quotient_order = i_dividende_polynom.max_order-i_diviseur_polynom.max_order;
   std::vector<F2_polynom> result;
   F2_polynom rest_polynom;
   F2_polynom quotient_polynom;

   if (quotient_order<0){
       quotient_order = 0;
   }
   quotient_polynom = F2_polynom(quotient_coefficients,quotient_order);


   // Initialization of the rest polynom
   rest_polynom = i_dividende_polynom;

   F2_polynom s;

   // Division
   while((rest_polynom.max_order>=i_diviseur_polynom.max_order) && (i_diviseur_polynom.max_order > -1)){
       quotient_coefficients[(rest_polynom.max_order-i_diviseur_polynom.max_order)/32] = (1<<((rest_polynom.max_order-i_diviseur_polynom.max_order)%32));
       s = F2_polynom(quotient_coefficients,(rest_polynom.max_order-i_diviseur_polynom.max_order));
       quotient_coefficients[(rest_polynom.max_order-i_diviseur_polynom.max_order)/32] = 0;

       quotient_polynom = addition(quotient_polynom,s, quotient_coefficients);
       rest_polynom = substraction(rest_polynom,multiplication(i_diviseur_polynom,s,quotient_coefficients),quotient_coefficients);

   }
   result.push_back(quotient_polynom);
   result.push_back(rest_polynom);

   return result;
}


F2_polynom convert_integers_to_polynoms(unsigned int  *integers_array, const int size){
   /*
    Function which transforms integers to polynom
    */




   int max_order = 2*size-1;


   return F2_polynom(integers_array,max_order);
}


std::vector<std::vector<unsigned int>> convert_polynoms_to_integers(const std::vector<std::vector<F2_polynom>> &polynoms_vector, int size){
   /*
    Function which transforms polynom vectors to integers
    */




   std::vector<std::vector<unsigned int>> integers_vectors;
   for (int i =0; i<polynoms_vector.size();i++){
       integers_vectors.push_back(std::vector<unsigned int>());
       for (int k = 0; k<(size/32);k++){
           integers_vectors[i].push_back(polynoms_vector[i][0].coefficients[k]);
       }


       for (int k = 0; k<(size/32);k++){
           integers_vectors[i].push_back(polynoms_vector[i][1].coefficients[k]);
       }
   }


   return integers_vectors;
}



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

   unsigned int* temp_coeff = new unsigned int[(remainder.max_order/32)+1];

   std::vector<F2_polynom> division_result;


   bool start = false;


   while(remainder.max_order>-1){
       if (start){
           dividende_polynom = diviseur_polynom;
           diviseur_polynom = remainder;

       }
       start = true;
       division_result = division(dividende_polynom,diviseur_polynom,temp_coeff);
       remainder = division_result[1];
   }


   if (diviseur_polynom.max_order == -1){
       unsigned int *coeff = new unsigned int [1];
       coeff[0] = 1;
       diviseur_polynom = F2_polynom(coeff,0);
   }
   return diviseur_polynom;


}





std::vector<F2_polynom> compute_Berlekamp_matrix(const std::vector<F2_polynom> &i_vectors, int size){
   /*
    Function which computes the Berlekamp matrix which can be used to compute polynom factorization.
    */


   std::vector<F2_polynom> matrix;


   // Loop over input vectors
   for (int i = 0; i<i_vectors.size();i++){


       unsigned int *coeff_vector = new unsigned int[(size/32)+1];
       for (int j = 0; j<(size/32)+1;j++){
           coeff_vector[j] = 0;
       }
       for (int j = 0;j<(i_vectors[i].max_order/32)+1;j++){
           coeff_vector[j] = i_vectors[i].coefficients[j] ;
       }
       coeff_vector[(i/32)] ^= (1<<(i%32));


       unsigned int *coeff_identity = new unsigned int[(size/32)+1];
       for (int j = 0; j<(size/32)+1;j++){
           coeff_identity[j] = 0;
       }
       coeff_identity[i/32] ^= (1<<(i%32));


       unsigned int *new_coeff = new unsigned int[2*((size/32)+1)];
       for (int j = 0; j<(size/32)+1;j++){
           new_coeff[j] = coeff_identity[j];
       }
       for (int j = (size/32)+1; j<2*((size/32)+1);j++){
           new_coeff[j] = coeff_vector[j-((size/32)+1)];
       }


       matrix.push_back(F2_polynom(new_coeff,64*((size/32)+1)));
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

   unsigned int* temp_coeff = new unsigned int[2*((matrix.size()/32)+1)];

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
               if (((matrix[i].coefficients[(matrix[row_indicator].max_order)/32]>>((matrix[row_indicator].max_order)%32))&1) == 1){
                   matrix[i] = addition(matrix[i],matrix[row_indicator],temp_coeff);
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

       if (i_matrix[i].max_order<(32*(i_matrix.size()/32+1))){
           unsigned int *coeff = new unsigned int[(i_matrix.size()/32)+1];
           for (int j = 0; j < (i_matrix.size()/32)+1; j++){
               coeff[j] = i_matrix[i].coefficients[j];
           }

           F2_polynom polynom = F2_polynom(coeff,32*(i_matrix.size()/32)+31);

           if (polynom.max_order>0){
               null_vectors.push_back(polynom);
           }
       }
   }


   return null_vectors;
}


std::vector<F2_polynom> compute_Berlekamp_factorization(const F2_polynom &i_polynom, unsigned int* temp_coeff){
   /*
    Function which computes the Berlekamp factorization of a square free polynom.
    */


   // Initialization of variables
   std::vector<F2_polynom> division_result, vector_result, matrix, null_vectors, factorization, recursive_factorization;
   F2_polynom vector, factor_1, factor_2;

   unsigned int *coeff = new unsigned int [((2*i_polynom.max_order)/32)+1];
   for (int j = 0; j < ((2*i_polynom.max_order)/32)+1; j++){
       coeff[j] =0;
   }

   // Computation of new vectors
   for (int i = 0; i<i_polynom.max_order/2;i++){
       coeff[(2*i)/32] ^= (1<<((2*i)%32));
       vector_result.push_back(F2_polynom(coeff,2*i));
       coeff[(2*i)/32] = 0;
   }
   for (int i = i_polynom.max_order/2; i<i_polynom.max_order;i++){
       coeff[(2*i)/32] ^= (1<<((2*i)%32));
       vector = F2_polynom(coeff,2*i);
       division_result = division(vector,i_polynom,temp_coeff);
       vector_result.push_back(division_result[1]);
       coeff[(2*i)/32] = 0;
   }


   // Computation of matrix
   matrix = compute_Berlekamp_matrix(vector_result, i_polynom.max_order-1);


   // Resolution of the matrix
   matrix = resolve_Berlekamp_matrix(matrix);


   // Resolution of the matrix
   null_vectors = find_null_vectors(matrix);


   if (null_vectors.size()>0){


       std::vector<F2_polynom> sorted_factor;
       for (int i = 0; i<null_vectors.size();i++){
           factor_1 = compute_greatest_common_divisor(i_polynom,null_vectors[0]);
           division_result = division(i_polynom,factor_1,temp_coeff);
           factor_2 = division_result[0];


           if (factor_1.max_order < factor_2.max_order){
               int num = 0;
               for (int j = 0; j<sorted_factor.size();j++){
                   if (factor_1.max_order>sorted_factor[num].max_order){
                       num ++;
                   }
               }
               sorted_factor.insert(sorted_factor.begin()+num,factor_1);
               sorted_factor.insert(sorted_factor.end()-num,factor_2);
           }
           else{
               int num = 0;
               for (int j = 0; j<sorted_factor.size();j++){
                   if (factor_2.max_order>sorted_factor[num].max_order){
                       num ++;
                   }
               }
               sorted_factor.insert(sorted_factor.begin()+num,factor_2);
               sorted_factor.insert(sorted_factor.end()-num,factor_1);
           }

       }

       int cumul_order = 0;
       int indice = 0;
       while (cumul_order<i_polynom.max_order){
           cumul_order+=sorted_factor[indice].max_order;
           factorization.push_back(sorted_factor[indice]);
           indice++;
       }

   }
   else {
       factorization.push_back(i_polynom);
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


       unsigned int *coeff = new unsigned int [size/16];
       coeff[0] = 1;
       F2_polynom left = F2_polynom(coeff,0);
       F2_polynom right = F2_polynom(coeff,0);
       coeff[0] = 0;



       for (int j = 0; j<factorization.size();j++){
           if (possible_solutions[i][j] == 0){
               left = multiplication(left,factorization[j],coeff);
           }
           else {
               right = multiplication(right,factorization[j],coeff);
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


/////////////
/// \brief compute_degree_polynom
/// \param i_polynom
/// \param maximum_size
/// \param o_degree
///
///
///
///
///
///
///
///
///
///
///
/// /////////


void slide_left_polynom( uint32_t* polynom, const int &polynom_size, int slide_value )
{

  while (slide_value > 0)
  {
    int top_num_uint = polynom_size-1;
    int fractional_slide_value      = slide_value & 0x1f;

    if (fractional_slide_value == 0)
    {
      while (top_num_uint > 0)
      {
        polynom[ top_num_uint ] = polynom[top_num_uint-1];
        top_num_uint--;
      }
      polynom[0] = 0;
      fractional_slide_value = 32;
    }
    else
    {
      while (top_num_uint > 0)
      {
        uint32_t top_num = polynom[ top_num_uint ] << fractional_slide_value;
        uint32_t low_num = polynom[ top_num_uint-1 ] >> (32-fractional_slide_value);

        polynom[ top_num_uint ] = top_num | low_num;
        top_num_uint--;
      }
      polynom[0] <<= fractional_slide_value;
    }

    slide_value -= fractional_slide_value;
  }
}


void slide_right_polynom( uint32_t* polynom, const int &polynom_size, int slide_value )
{

  while (slide_value > 0)
  {
    int top_num_uint = 0;
    int fractional_slide_value      = slide_value & 0x1f;

    if (fractional_slide_value == 0)
    {
      while (top_num_uint < polynom_size-1)
      {
        polynom[ top_num_uint ] = polynom[top_num_uint+1];
        top_num_uint++;
      }
      polynom[polynom_size-1] = 0;
      fractional_slide_value = 32;
    }
    else
    {
      while (top_num_uint < polynom_size-1)
      {
        uint32_t low_num = polynom[ top_num_uint ] >> fractional_slide_value;
        uint32_t top_num = polynom[ top_num_uint+1 ] << (32-fractional_slide_value);

        polynom[ top_num_uint ] = top_num | low_num;
        top_num_uint++;
      }
      polynom[polynom_size-1] >>= fractional_slide_value;
    }

    slide_value -= fractional_slide_value;
  }
}

int compute_degree_polynom(uint32_t* i_polynom,const int &polynom_size){
    /*
     Compute the highest degree of a polynom recursively to optimize the research
     */

    int upper_decompositon = polynom_size-1;

    do {
        if (i_polynom[upper_decompositon]){
            break;
        }
        upper_decompositon --;
    }
    while(upper_decompositon>=0);

    int degree = 32*upper_decompositon+31-__builtin_clz(i_polynom[upper_decompositon]);

    return degree;
}


bool is_zero_polynom(uint32_t* i_polynom,const int &polynom_size){
    /*
     Checks if the polynom is null
     */

    int upper_decompositon = polynom_size-1;

    do {
        if (i_polynom[upper_decompositon]){
            return 0;
        }
        upper_decompositon --;
    }
    while(upper_decompositon>=0);

    return 1;
}

void addition(uint32_t *io_polynom_left,uint32_t *i_polynom_right, const int& polynom_size){
   /*
    Definition of the addition between two polynoms
    */

   int upper_uint_num = polynom_size-1;

   do{
       io_polynom_left[upper_uint_num] ^= i_polynom_right[upper_uint_num];
       upper_uint_num--;
   }
   while(upper_uint_num>=0);

}

void multiplication(uint32_t *io_polynom_left, uint32_t *i_polynom_right, const int& polynom_size){
   /*
    Definition of the multiplication between two polynoms
    */

    // Initialization of variables

    int num_byte = 0;
    int max_bytes = 32*polynom_size;
    uint32_t temp_polynom[polynom_size];
    memcpy(temp_polynom,io_polynom_left,sizeof (uint32_t)*polynom_size);

    int top_uint = polynom_size-1;
    do{io_polynom_left[top_uint]=0;}
    while(top_uint>=0);

    // Multiplication

    do{
        if ((i_polynom_right[(num_byte>>5)]>>(num_byte&0x1f))&1){
            addition(io_polynom_left,temp_polynom,polynom_size);
        }
        slide_left_polynom(temp_polynom,polynom_size,1);
        num_byte++;
    }
    while(num_byte<(32*polynom_size));

}

void division(uint32_t *i_polynom_dividend, uint32_t *i_polynom_divisor, uint32_t* polynom_quotient, uint32_t* polynom_remainder, const int& polynom_size){
   /*
    Definition of the multiplication between two polynoms
    */


    // Initialization of the remainder and quotient polynoms
    memset(polynom_quotient,0,sizeof(uint32_t)*polynom_size);
    memcpy(polynom_remainder, i_polynom_dividend, sizeof(uint32_t)*polynom_size);

    uint32_t temp_polynom[polynom_size];

    int degree_dividend = compute_degree_polynom(i_polynom_dividend, polynom_size);
    int degree_divisor = compute_degree_polynom(i_polynom_divisor, polynom_size);


    int degree_quotient = degree_dividend-degree_divisor;
    int current_degree = 0;

    // Division
    do {
        slide_left_polynom(polynom_quotient,polynom_size,1);
        int max_degree_remainder = degree_dividend-current_degree;

        if ((polynom_remainder[(max_degree_remainder)>>5]<<((max_degree_remainder)&0x1f))&1){
            polynom_quotient[0]^=1;
            memcpy(temp_polynom, i_polynom_divisor, sizeof(uint32_t)*polynom_size);
            slide_left_polynom(temp_polynom,polynom_size,degree_quotient-current_degree);
            addition(polynom_remainder,temp_polynom,polynom_size);
        }
        current_degree++;
    }
    while(current_degree<=degree_quotient);
}


void compute_derivation(uint32_t* encoded_polynom, const int &polynom_size,  uint32_t* derivation){
   /*
    Function which computes the derivative of a polynom.
    */

   // First set all the odd elements to 0
   memset( derivation, 0xAA, sizeof(uint32_t) * polynom_size);

   // Multiplication with encoded polynom to keep only odd elements
   int top_num_uint = polynom_size-1;
   do
   {
     derivation[ top_num_uint ] = derivation[ top_num_uint ] & encoded_polynom[ top_num_uint ];
     top_num_uint--;
   }
   while (top_num_uint >= 0);

   // Slide all the elements to the right from one byte
    slide_right_polynom(derivation,polynom_size,1);

}


void compute_square_root(uint32_t* encoded_polynom, const int &polynom_size,  uint32_t* square_polynom){
   /*
    Function which computes the square root of a polynom.
    */


   memset(square_polynom,0, sizeof (uint32_t)*polynom_size);
   int byte = 0;
   do {
       square_polynom[byte>>5] ^= ((encoded_polynom[byte>>4]>>((byte<<1)&0x1f))&1)<<(byte&0x1f);
   }
   while((byte<<1)<((polynom_size<<5)-1));

}


void compute_greatest_common_divisor(uint32_t* polynom_a, uint32_t* polynom_b, uint32_t* common_divisor, const int &polynom_size){
   /*
    Function which computes the greatest common divisor between two polynoms using Euclidean algorithm.
    */


   // Initialization of variables
   uint32_t dividende_polynom[polynom_size];
   uint32_t diviseur_polynom[polynom_size];
   uint32_t remainder_polynom[polynom_size];
   uint32_t quotient_polynom[polynom_size];

   int degree_a = compute_degree_polynom(polynom_a,polynom_size);;
   int degree_b = compute_degree_polynom(polynom_b,polynom_size);
;

   if (degree_a>degree_b){
       memcpy(dividende_polynom,polynom_a,sizeof(uint32_t)*polynom_size);
       memcpy(diviseur_polynom,polynom_b,sizeof(uint32_t)*polynom_size);
   }
   else {
       memcpy(dividende_polynom,polynom_b,sizeof(uint32_t)*polynom_size);
       memcpy(diviseur_polynom,polynom_a,sizeof(uint32_t)*polynom_size);
   }


   division(dividende_polynom, diviseur_polynom,quotient_polynom,remainder_polynom,polynom_size);

   do {
        memcpy(dividende_polynom,diviseur_polynom,sizeof (uint32_t)*polynom_size);
        memcpy(diviseur_polynom,remainder_polynom,sizeof (uint32_t)*polynom_size);
        division(dividende_polynom, diviseur_polynom,quotient_polynom,remainder_polynom,polynom_size);
   }
   while(~is_zero_polynom(remainder_polynom,polynom_size));

   if(is_zero_polynom(diviseur_polynom,polynom_size)){
        common_divisor[0] = 1;
        return;
   }

   memcpy(common_divisor,diviseur_polynom,sizeof(uint32_t)*polynom_size);

}


void compute_square_free_factorization(uint32_t* polynom, const int &polynom_size, uint32_t* factorization, int &num_factor){
   /*
    Function which computes the square free factorization of a polynom.
    */


   // Derivation
   uint32_t derivation[polynom_size];
   compute_derivation(polynom,polynom_size,derivation);

   // If the derivation is not null
   if (~is_zero_polynom(derivation,polynom_size)){

       uint32_t common_divisor_derivation[polynom_size];
       // Computation of the greatest common divisor between the polynom and the derivation
       compute_greatest_common_divisor(polynom,derivation,common_divisor_derivation,polynom_size);

       // Check that the derivation divide the polynom
       int degree = compute_degree_polynom(common_divisor_derivation,polynom_size);
       if (degree>0){
           uint32_t quotient[polynom_size];
           uint32_t remainder[polynom_size];
           division(polynom,common_divisor_derivation,quotient,remainder,polynom_size);

           // Recursive call over the two factors
           compute_square_free_factorization(common_divisor_derivation,polynom_size,factorization,num_factor);
           compute_square_free_factorization(quotient,polynom_size,factorization,num_factor);
       }
       else{
           memcpy(factorization+polynom_size*num_factor,polynom,sizeof(uint32_t)*polynom_size);
       }
   }


   // If the derivation is null
   else {
       // Computation of the square root
       uint32_t square_root[polynom_size];
       compute_square_root(polynom,polynom_size,square_root);

       // Recursion over an element of the square root
       int current_num_factor = num_factor;
       compute_square_free_factorization(square_root,polynom_size,factorization,num_factor);

       // Addition of the factors corresponding to the other part of the factor
       memcpy(factorization+polynom_size*num_factor,factorization+polynom_size*current_num_factor,sizeof(uint32_t)*polynom_size*(num_factor-current_num_factor));
       num_factor += (num_factor-current_num_factor);
   }
}


void polynom_factorization(uint32_t* encoded_polynom, const int &polynom_size, const int &encoded_polynom_degree,  uint32_t* factorization){
   /*
    Function which computes the factorization of a polynom.
    */


   uint32_t square_free_factorization [(polynom_size>>5)*encoded_polynom_degree];

   square_free_factorization = compute_square_free_factorization(i_polynom,temp_coeff);

   for (int i = 0; i <square_free_factorization.size();i++){
       recursive_factorization = compute_Berlekamp_factorization(square_free_factorization[i], temp_coeff);
       for (int j = 0; j<recursive_factorization.size();j++){
           factorization.push_back(recursive_factorization[j]);
       }
   }



}



int main()
{

   // Inputs reading
   int size;
   cin >> size;

   int polynom_size = (size>>4);

   uint32_t encoded [polynom_size];
   for (int i = 0; i < size / 16; i++) {   // Read size / 16 integers to a
       cin >> hex >> encoded[i];
   }

   // Computation of the maximum degree
   int encoded_polynom_degree = compute_degree_polynom(encoded,polynom_size);

   // Initialization of the factors
   uint32_t  factorization[(polynom_size>>1)*encoded_polynom_degree];

   // Factorisation of encoded polynoms
   polynom_factorization(encoded,polynom_size,encoded_polynom_degree,factorization);

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
