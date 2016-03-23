# Programming Assnt. 2

> Optimizing Strassen's Algorithm for Matrix Multiplication

## TODOs

* All calculations up to 32 bits
* Randomly populate for testing with {0,1}, {0,1,2}, {-1,0,1}
* Inline strassen implementation
* Increase modularity
* Excel analysis / i.e. printing to CSV
* Do mathy analysis 
* Padding (lazy data structures?)
* Decide how to do memory management for auxiliary matrices

* Writeup
 * left_matrix v. right_matrix to promote cache efficiency
 * inline strassen for space
 * could have used DP to optimize standard matrix mult further
 * How low was cross-over point
 * Describe bugs and difficulties
 * What types of matrices were tested? Does it matter?

* Analysis
 * Building the Ms involves 10 additions/subtractions of matrices each of dim n/2. 7 multiplications.
 * Building the Cs involves 8 additions/subtractions each of dim n/2.
 * Conventionally we consider the Ms to be auxiliary. 
 * "Both initial matrices must have their dimensions expanded to the next power of 2, which results in storing up to four times as many elements, and the seven auxiliary matrices each contain a quarter of the elements in the expanded ones" ~ Wikipedia
 * Change matrix passed in as param inline with the C matrices? 
