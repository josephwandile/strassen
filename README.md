# Programming Assnt. 2

> Optimizing Strassen's Algorithm for Matrix Multiplication

## Analysis

* All timing data is generated with format `[strassen|traditional]_[smallest_dim]_[largest_dim]_[trials_per_dim]` and is found in the root directory.

## Style Rules

* `thisIsAFunction`
* `this_is_a_variable`
* Separate functions by two lines
* Block off logical sections of code with descriptive comments
* `pointer*`, rather than `pointer *`
* When writing code which has different behavior for Strassen's vs. traditional multiplication use the variable `bool using_strassen`

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
 * inline strassen for space
 * could have used DP to optimize standard matrix mult further
 * How low was cross-over point
 * Describe bugs and difficulties
 * What types of matrices were tested? Does it matter?
 * m matrices generation
 * addition occurs using references, not new matrices
 * there is a trade off between space and speed for modifying submatrices. it can be done with 4 passes but then requires that more Pis are stored. Or you can store only one Pi at a time and do it in 12 passes. Note it might be possible to temporarily store some of the Pis in C.

* Analysis
 * Building the Ms involves 10 additions/subtractions of matrices each of dim n/2. 7 multiplications.
 * Building the Cs involves 8 additions/subtractions each of dim n/2.
 * Conventionally we consider the Ms to be auxiliary.
 * "Both initial matrices must have their dimensions expanded to the next power of 2, which results in storing up to four times as many elements, and the seven auxiliary matrices each contain a quarter of the elements in the expanded ones" ~ Wikipedia
 * Change matrix passed in as param inline with the C matrices?
