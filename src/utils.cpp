//
//  utils.cpp
//  F2_Polynom_Factorization
//
//  Created by Valentin Figué  on 25/12/2019.
//

#include "utils.hpp"

F2_polynom convert_integers_to_polynoms(unsigned int  *integers_array, const int size){
   /*
    Function which transforms integers to polynom
    */


   std::vector<int> polynom_coefficient;
   for (int i =0; i<2*size;i++){
       int coeff_value = (integers_array[i/32]>>(i%32))&1;
       polynom_coefficient.push_back(coeff_value);
   }

   return F2_polynom(polynom_coefficient);
};

std::vector < std::vector <unsigned int> > convert_polynoms_to_integers(const std::vector < std::vector < F2_polynom > > &polynoms_vector, int size){
   /*
    Function which transforms polynom vectors to integers
    */


   std::vector < std::vector < unsigned int > > integers_vectors;
   for (int i =0; i<polynoms_vector.size();i++){
       integers_vectors.push_back(std::vector<unsigned int>());
       for (int k = 0; k<(size/32);k++){
           unsigned int first = 0;
           unsigned int factor = 1;
           for (int j = 32*(k); j<32*(k+1);j++){
               if (j<polynoms_vector[i][0].max_order+1){
                   first += polynoms_vector[i][0].coefficients[j]*factor;
                   factor = 2*factor;
               }
           }
           integers_vectors[i].push_back(first);
       }

       for (int k = 0; k<(size/32);k++){
           unsigned int second = 0;
           unsigned int factor = 1;
           for (int j = 32*(k); j<32*(k+1);j++){
               if (j<polynoms_vector[i][1].max_order+1){
                   second += polynoms_vector[i][1].coefficients[j]*factor;
                   factor = 2*factor;
               }
           }
           integers_vectors[i].push_back(second);
       }
   }

   return integers_vectors;
}
