# Programming Assnt. 2

> Optimizing Strassen's Algorithm for Matrix Multiplication

## TODOs

* All calculations up to 32 bits
* Randomly populate for testing with {0,1}, {0,1,2}, {-1,0,1}

* Writeup
 * left_matrix v. right_matrix to promote cache efficiency
 * inline strassen for space
 * could have used DP to optimize standard matrix mult further
 * How low was cross-over point
 * Describe bugs and difficulties
 * What types of matrices were tested? Does it matter?

* Strassen
 * Padding
 * Avoid copying data (i.e. perhaps build smaller matrices by simply referring to parts of the larger one)
 * When to use auxiliary matrices
