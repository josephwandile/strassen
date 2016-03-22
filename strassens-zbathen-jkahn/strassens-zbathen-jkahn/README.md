# Programming Assnt. 2

> Optimizing Strassen's Algorithm for Matrix Multiplication

## TODOs

* Create a class for matrices
 * `[[],[],[]]` where each subarray is a row
 * build from text file row by row
 * print diagonal
 * all calculations up to 32 bits
 * Randomly populate for testing with {0,1}, {0,1,2}, {-1,0,1}

* Writeup
 * What optimizations in the code
 * How low was cross-over point
 * Describe bugs and difficulties
 * What types of matrices were tested? Does it matter?

* Normal Mult
 * What does "looping through in right order" refer to?

* Strassen
 * Padding
 * Avoid copying data (i.e. perhaps build smaller matrices by simply referring to parts of the larger one)
