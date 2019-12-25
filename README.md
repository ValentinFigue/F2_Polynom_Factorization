# F2 Polynom Factorization

## Goal

This repository implements a way to factorize a F2 polynom in irreducible factors. This was inspired by a coding challenge from Nintendo (https://www.codingame.com/training/expert/nintendo-sponsored-contest).

The factorization used to reduce F2 polynoms is an implementation of famous Berlekamp algorithm. For more information, you can read this article :

https://en.wikipedia.org/wiki/Berlekamp%27s_algorithm

Feel free to use the different functions for your own purposes.

## Structure

The code is decomposed in three different source files :

- **F2_Polynom** regroups the object which represents a F2 polynom and the different basic operations (addition, multiplication, division, ...).
- **Utils** gathers different conversion function between integers vectors and polynoms.
- **Arithmetic** explicits the factorization itself and the different sub algorithms required.

## How to use it

This repository implements a sample which factorize polynom. To use it, you just need to run the following commands in the code folder of the  :

```terminal
cmake
make
```

Then you just have to specify the degree of the polynom to factorize and specify the different coefficient in hexadecimal, when you run the build **F2_Polynom_Factorization**.
