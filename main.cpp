#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <assert.h>
#include <ctime>
#include <random>
#include <math.h>

using namespace std;
int CUTOFF = 115;
const bool IN_DEV = false; // Runs a couple simple tests before executing main commands as a sanity check
const string OUTPUT_SEPERATOR = "----------------------------------------\n\n";

default_random_engine generator;
uniform_int_distribution<int> distribution(-1,2);


/*

 Definition of Matrix as well as utility functions for generating, populating, and printing them.

 */
typedef struct Matrix {
    int dimension;
    vector<vector<int>> entries;
} Matrix;

Matrix* instantiateMatrix(int dimension) {

    vector<vector<int>> entries(dimension, vector<int>(dimension));

    Matrix* matrix = (Matrix*) malloc(1000000);

    matrix->dimension = dimension;
    matrix->entries = entries;

    return matrix;
}


// Does the actual IO work for populating matrices
void populateMatrix(Matrix* matrix, string infile, int position, int dimension) {

    ifstream inputfile(infile);
    string element = "";

    if (inputfile.is_open()) {

        // Skips to the appropriate line
        for (int i = 0; i < position; i++){
            getline(inputfile, element);
        }

        for (int row = 0; row < dimension; row++) {
            for (int col = 0; col < dimension; col++) {
                if (getline(inputfile, element)) {
                    matrix->entries[row][col] = stoi(element);
                    element = "";
                }
                else {
                    cout << "File does not contain enough data" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
}


// Instantiates and populates a matrix from a specific point in the input file
Matrix* buildMatrix(string infile, int position, int dimension) {

    Matrix* matrix = instantiateMatrix(dimension);
    populateMatrix(matrix, infile, position, dimension);
    return matrix;
}


// Required output format for assignment
void printMatrix(Matrix* matrix, bool printing_diagonal=true) {

    if (printing_diagonal) {

        if (IN_DEV) {
            cout << "Matrix Diagonal" << endl;
        }
        for (int i = 0; i < matrix->dimension; i++) {
            cout << matrix->entries[i][i] << endl;
        }
        cout << "\n";

    } else {

        cout << "Formatted Matrix" << endl;
        int row, col, entry;
        string entry_as_str;
        for (row = 0; row < matrix->dimension; row++) {
            for (col = 0; col < matrix->dimension; col++) {
                entry = matrix->entries[row][col];
                cout << setw(5) << entry;
            }
            cout << endl;
        }
    }
}


/*

 STANDARD MULTIPLICATION

 Used by Strassen when the dimension is below an experimentally defined
 cutoff point.

 */
Matrix* tradMult(Matrix* l_matrix, Matrix* r_matrix) {

    if (l_matrix->dimension != r_matrix->dimension) {
        throw invalid_argument("tradMult expects square matrices.");
    }

    int dimension = l_matrix->dimension;

    // Always outputs a matrix in standard form (i.e. Matrix[rows][cols])
    Matrix* ans_mat = instantiateMatrix(dimension);

    int i, j, k;

    /*

     Indices carefully chosen to improve cache performance. Note that the innermost loop only jumps
     between contiguous blocks of memory (from vector to vector) the minimum number of times possible.

     Performance difference on PC?

     73.4703 seconds to multiply matrices of dimension = 2048 without cache optization.

     vs.

     6.0186 seconds with optimization. Total insanity.

     */
    for (i = 0; i < dimension; i++) { // Rows
        for (j = 0; j < dimension; j++) { // Cols
            for (k = 0; k < dimension; k++) { // Iterator
                ans_mat->entries[i][k] += l_matrix->entries[i][j] * r_matrix->entries[j][k];
            }
        }
    }
	return ans_mat;
}


/*

 STRASSEN'S METHOD

 A variety of helper functions as well as the actual recursive implementation.

 Instantiate new matric (which will be passed as a parameter to Strassen's
 by adding or subtracting two arbitrary submatrices.
 */
Matrix* combineSubmatrices(Matrix* refmat, int j_1, int i_1, int j_2, int i_2, bool add) {

    int dimension = refmat->dimension / 2;
	Matrix* ans_mat = instantiateMatrix(dimension);

	if (add) {
		for (int row = 0; row < dimension; row++) {
			for (int col = 0; col < dimension; col++) {
				ans_mat->entries[row][col] = refmat->entries[i_1 + row][j_1 + col] + refmat->entries[i_2 + row][j_2 + col];
			}
		}
	} else
	{
		for (int row = 0; row < dimension; row++) {
			for (int col = 0; col < dimension; col++) {
				ans_mat->entries[row][col] = refmat->entries[i_1 + row][j_1 + col] - refmat->entries[i_2 + row][j_2 + col];
			}
		}
	}
	return ans_mat;
}


/*

 As with combineSubmatrices(), except we instantiate new matrix from
 single arbitrary submatrix.

 */
Matrix* extractSubmatrix(Matrix* refmat, int j, int i) {

    int dimension = refmat->dimension / 2;
	Matrix* ans_mat = instantiateMatrix(dimension);

	for (int row = 0; row < dimension; row++) {
		for (int col = 0; col < dimension; col++) {
			ans_mat->entries[row][col] = refmat->entries[i + row][j + col];
		}
	}
	return ans_mat;
}


// Add a row and a col of zeroes to matrix
void pad(Matrix* mat) {

    int dimension = mat->dimension;

	// Pad a row at the bottom
	vector<int> zeroes(dimension + 1, 0);
	mat->entries.push_back(zeroes);

	// Pad a column on the right
	for (int i = 0; i < dimension; i++) {
		mat->entries[i].push_back(0);
	}

    mat->dimension++;
}


// Remove the last row and last col of a matrix
void unpad(Matrix* matrix) {

    int dimension = matrix->dimension;

	// Remove last row from bottom
	matrix->entries.pop_back();

	// Remove last col from right
	for (int i = 0; i < dimension; i++) {
		matrix->entries[i].pop_back();
	}
	matrix->dimension = dimension - 1;
}


Matrix* strassenMult(Matrix* A, Matrix* B) {

    if (A->dimension != B->dimension) {
        throw invalid_argument("strassenMult expects square matrices.");
    }

    int dimension = A->dimension;

	if (dimension <= CUTOFF){
		return tradMult(A, B); // TODO what is the implication of strassen referencing P with the cutoff point and normal mult?
	} else {

		bool padding = false;

		// If matrix is of odd dimension, pad a row and col of zeroes
		if ((dimension & 1) != 0) {
			padding = true;
			pad(A);
			pad(B);
			dimension++; // Keep local representation of dim up to date
		}

        // Matrix to be filled and returned
        Matrix* C = instantiateMatrix(dimension);

        int half_dim = dimension / 2;

        /*

         Using an auxiliary matrix to store the currently needed product increases space efficiency significantly,
         but (at least in our implementation) makes it harder to minimize the number of algebraic operations.

         There's a clear trade-off between time and space here (and in some sense simplicity, since making time and space usage
         both optimal requires difficult-to-reason-about code with the majority of operations happening in place. (Though this can
         be made less burdensome for a 'client' if the proper abstractions are used.)

         */
		Matrix* m1a = combineSubmatrices(A, 0, 0, half_dim, half_dim, true);
		Matrix* m1b = combineSubmatrices(B, 0, 0, half_dim, half_dim, true);
        Matrix* m1 = strassenMult(m1a, m1b);
        free( m1a);
        free( m1b);

		Matrix* m2a = combineSubmatrices(A, 0, half_dim, half_dim, half_dim, true);
		Matrix* m2b = extractSubmatrix(B, 0, 0);
        Matrix* m2 = strassenMult(m2a, m2b);
        free( m2a);
        free( m2b);

		Matrix* m3a = extractSubmatrix(A, 0, 0);
		Matrix* m3b = combineSubmatrices(B, half_dim, 0, half_dim, half_dim, false);
        Matrix* m3 = strassenMult(m3a, m3b);
        free( m3a);
        free( m3b);

        Matrix* m6a = combineSubmatrices(A, 0, half_dim, 0, 0, false);
        Matrix* m6b = combineSubmatrices(B, 0, 0, half_dim, 0, true);
        Matrix* m6 = strassenMult(m6a, m6b);
        free( m6a);
        free( m6b);

        // c22 = m1 - m2 + m3 + m6
        for (int row = 0; row < half_dim; row++) {
            for (int col = 0; col < half_dim; col++) {
                C->entries[row + half_dim][col + half_dim] =
                    m1->entries[row][col] - m2->entries[row][col]
                    + m3->entries[row][col] + m6->entries[row][col];
            }
        }

        free(m6);

		Matrix* m4a = extractSubmatrix(A, half_dim, half_dim);
		Matrix* m4b = combineSubmatrices(B, 0, half_dim, 0, 0, false);
        Matrix* m4 = strassenMult(m4a, m4b);
        free( m4a);
        free( m4b);

        // c21 = m2 + m4
        for (int row = 0; row < half_dim; row++) {
            for (int col = 0; col < half_dim; col++) {
                C->entries[row + half_dim][col] = m2->entries[row][col] + m4->entries[row][col];
            }
        }

        free(m2);

		Matrix* m5a = combineSubmatrices(A, 0, 0, half_dim, 0, true);
		Matrix* m5b = extractSubmatrix(B, half_dim, half_dim);
        Matrix* m5 = strassenMult(m5a, m5b);
        free( m5a);
        free( m5b);

        // c12 = m3 + m5
        for (int row = 0; row < half_dim; row++) {
            for (int col = 0; col < half_dim; col++) {
                C->entries[row][col + half_dim] = m3->entries[row][col] + m5->entries[row][col];
            }
        }

        free(m3);

		Matrix* m7a = combineSubmatrices(A, half_dim, 0, half_dim, half_dim, false);
		Matrix* m7b = combineSubmatrices(B, 0, half_dim, half_dim, half_dim, true);
        Matrix* m7 = strassenMult(m7a, m7b);
        free( m7a);
        free( m7b);

        // c11 = m1 + m4 - m5 + m7
        for (int row = 0; row < half_dim; row++) {
            for (int col = 0; col < half_dim; col++) {
				C->entries[row][col] = m1->entries[row][col] + m4->entries[row][col]
					- m5->entries[row][col] + m7->entries[row][col];
			}
		}

        free(m1);
        free(m4);
        free(m5);
        free(m7);

		// Remove padding if necessary
		if (padding) {
			unpad(C);
		}

		return C;
	}
}

/*

 TESTING

 */
Matrix* genRandMatrix(int dimension) {

    int new_entry, i, j;
    Matrix* matrix = instantiateMatrix(dimension);

    for (i = 0; i < dimension; i++) {
        for (j=0; j < dimension; j++) {

            new_entry = distribution(generator);
            matrix->entries[i][j] = new_entry;
        }
    }
    return matrix;
}

bool matricesAreEqual(Matrix* A, Matrix* B) {

    bool are_equal = true;
    if (A->dimension != B->dimension) {
        return false;
    }
    int dimension = A->dimension;

    int i, j;
    for (i = 0; i < dimension; i++) {
        for (j = 0; j < dimension; j++) {

            if (A->entries[i][j] != B->entries[i][j]) {
                return false;
            }
        }
    }
    return are_equal;
}

void testingUtility(string infile, int dimension, bool use_random_matrices=true, bool using_strassen=true, bool printing_matrix=false) {

    if (use_random_matrices) {

        Matrix* A;
        Matrix* B;

        cout << "Testing randomly generated matrices of size: " << dimension << endl;

        /*
         Once upon time, there was a terrible, terrible bug.

         Which taught two young lads of yore, never to share global variables
         between functions.

         Amen
         */
        A = genRandMatrix(dimension);
        B = genRandMatrix(dimension);

        Matrix* C_trad = tradMult(A, B);
        Matrix* C_strass = strassenMult(A, B);

        free( A);
        free( B);

        /*
         Can't "know" the correct result a priori, so we assume that if each type of
         matrix multiplication is working correcly that they will return the same result
         */
        if (printing_matrix) {
            printMatrix(C_trad, false);
            printMatrix(C_strass, false);
        }

        assert(matricesAreEqual(C_strass, C_trad));
        cout << "Assertion True" << endl;

        free( C_strass);
        free( C_trad);

    } else {

        Matrix* A;
        Matrix* B;

        A = buildMatrix(infile, 0, dimension);
        B = buildMatrix(infile, dimension*dimension, dimension);

        Matrix* C;

        // Deterministic test so we can test Strassen and Trad independently of each other
        if (using_strassen) {
            cout << "Testing Strassen" << endl;
            C = strassenMult(A,B);
        } else {
            cout << "Testing Traditional Mult" << endl;
            C = tradMult(A,B);
        }

        free( A);
        free( B);

        // Left matrix built from test files
        Matrix* correct_C = buildMatrix(infile, dimension*dimension*2, dimension);

        if (printing_matrix) {
            printMatrix(C, false);
            printMatrix(correct_C, false);
        }

        assert(matricesAreEqual(correct_C, C));
        cout << "Assertion True" << endl;


        free( C);
        free( correct_C);

    }
}

/*

 TIMING UTILITY FUNCTIONS

 timeMatrixFromFile exists in case TFs/we want to time the execution of Strassen's on a known matrix pair

 findCrossover is a utility function to find the optimal cutoff point for matrix multiplication of arbitrary dimension

 timingUtility generates time data over a range of dimensions for either Strassen's or Traditional Multiplication

 */
void timeMatrixFromFile(string infile, int dimension) {

    clock_t construct_start = clock();
    Matrix* A = buildMatrix(infile, 0, dimension);
    Matrix* B = buildMatrix(infile, dimension*dimension, dimension);
    double construct_total = (clock() - construct_start) / (double)(CLOCKS_PER_SEC);

    // Assuming we'll only want to time Strassen's method direcly from input file
    cout << "Timing Strassen from input file." << endl << OUTPUT_SEPERATOR;
    clock_t mult_start = clock();
    Matrix* C = strassenMult(A,B);
    double mult_total = (clock() - mult_start) / (double)(CLOCKS_PER_SEC);

    free(A);
    free(B);
    free(C);

    cout << "Construction for matrix of size " << dimension << " took " << construct_total << "s" << endl;
    cout << "Multiplication took " << mult_total << "s" << endl;
}


int findCrossover(int dimension, int trials) {

    cout << "Finding optimal crossover point for matrix of size " << dimension << endl << OUTPUT_SEPERATOR;

    string output_file_name = "strassen_crossover_" + to_string(dimension) + "_" + to_string(trials) + ".txt";
    ofstream output_file;
    output_file.open(output_file_name);

    double min_time = 10000;
    int min_cutoff = 10000;

    int increasing_function_counter = 0;

    for (int cutoff = 16; cutoff < dimension; cutoff+=10) { // TODO Try to get ballpark

        CUTOFF = cutoff; // Mutate global variable

        double total_time = 0;

        for (int trial = 0; trial < trials; trial++) {

            Matrix* A = genRandMatrix(dimension);
            Matrix* B = genRandMatrix(dimension);

            clock_t start_mult = clock();
            Matrix* C = strassenMult(A,B);
            double mult_time = (clock() - start_mult) / (double)(CLOCKS_PER_SEC);

            free(A);
            free(B);
            free(C);

            total_time += mult_time;
        }

        double avg_time = total_time / trials;
        output_file << CUTOFF << "\t" << avg_time << endl;

        if (avg_time < min_time) {

            cout << "Found better cutoff of " << CUTOFF << endl;

            min_time = avg_time;
            min_cutoff = cutoff;
            increasing_function_counter = 0;

        } else {

            increasing_function_counter++;
        }
    }

    cout << "Optimal (lowest) cutoff of " << min_cutoff << " results in an average running time of " << min_time << endl;

    output_file.close();

    return min_cutoff;
}


void timingUtility(int lower_bound, int upper_bound, int trials, int interval, bool using_strassen=true) {

    int cur_matrix_dimension;
    string output_file_name;

    if (using_strassen) {
        cout << "Timing Strassen" << endl << OUTPUT_SEPERATOR;
        output_file_name = "strassen_";
    } else {
        cout << "Timing Traditional" << endl << OUTPUT_SEPERATOR;
        output_file_name = "traditional_";
    }

    // Build ouput file name
    output_file_name += to_string(lower_bound) + "_" + to_string(upper_bound) + "_" + to_string(interval) + "_" + to_string(trials) + ".txt";
    ofstream output_file;
    output_file.open(output_file_name);

    cur_matrix_dimension = lower_bound;
    while(cur_matrix_dimension <= upper_bound) {

//        CUTOFF = ceil(cur_matrix_dimension / 2);
        
        // For matrix of cur_matrix_dimension = n
        double total_mult_time = 0;
        double avg_mult_time = 0;

        cout << "Matrix of Size: " << cur_matrix_dimension << endl << OUTPUT_SEPERATOR;

        clock_t construct_start = clock();
        Matrix* A = genRandMatrix(cur_matrix_dimension);
        Matrix* B = genRandMatrix(cur_matrix_dimension);
        double construct_time = (clock() - construct_start) / (double)(CLOCKS_PER_SEC);
        cout << "Time taken to construct matrices: " << construct_time << "s" << endl;

        if (using_strassen) {
            cout << "CUTOFF of " << CUTOFF << endl;
        }
        
        for (int trial = 0; trial < trials; trial++) {

            clock_t mult_start = clock();
            Matrix* C;
            if (using_strassen) {
                C = strassenMult(A, B);
            } else {
                C = tradMult(A, B);
            }

            double mult_total = (clock() - mult_start) / (double)(CLOCKS_PER_SEC);
            total_mult_time += mult_total;

            free( C);
        }

        avg_mult_time = total_mult_time / trials;

        cout << "Average Time for Mult over " << trials << " trial/s: " << avg_mult_time << endl;
//        output_file << cur_matrix_dimension << "\t" << avg_mult_time << endl;

        cur_matrix_dimension += interval;

        free(A);
        free(B);
    }

    output_file.close();
}


/*

 PROGRAM INTERFACE

 */

int main(int argc, char* argv[]) {
	// check to ensure correct number of arguments
	if (argc != 4) {
        cout << "Please enter the correct number of arguments" << endl;
		return -1;
	}

    int flag = stoi(argv[1]);
    int dimension = stoi(argv[2]);
    string infile = argv[3];

    if (IN_DEV) {

        // Simple test cases to make sure nothing has gone totally wrong.
        testingUtility("test33.txt", 3, false, true); // Strassen
        testingUtility("test4d.txt", 4, false, true); // Strassen
        testingUtility("test33.txt", 3, false, false); // Traditional
        testingUtility("", 39, true, true); // Random Matrices
        testingUtility("", 16, true, true); // Random Matrices
        cout << "Basic Tests Pass. Executing instructions from command line." << endl << OUTPUT_SEPERATOR;
    }

    if (flag == 0) { // Production behavior

        if (!IN_DEV) {

            // buildMatrix(filename, read_from_position, buffer_length)
            Matrix* A = buildMatrix(infile, 0, dimension);
            Matrix* B = buildMatrix(infile, dimension*dimension, dimension);
            
            // Different CUTOFF for even vs. Odd
            if ((dimension & 1) != 0) {
                CUTOFF = 215; // Odd Cutoff
            } else {
                CUTOFF = 115; // Even Cutoff
            }

            
        	Matrix* C = strassenMult(A, B);
            printMatrix(C);
            free(A);
            free(B);
            free(C);
        }
        return 0;
    }

	if (flag == 1) { // Testing on randomly generated matrices of arbitrary size

        testingUtility(infile, dimension, true , true, true);
        return 0;
    }

    if (flag == 2) { // Run deterministic tests on Strassen and Trad

        // Normal
        testingUtility(infile, dimension, false, true, true);
        testingUtility(infile, dimension, false, false, true);

        return 0;
    }

    // lower bound - upper bound - trials - intervals - bool using_strassen
    if (flag == 3) { // Generate time data for Strassen

        // Strassen
        timingUtility(dimension, dimension, 1, 1, true);
    }

    if (flag == 4) { // Generate time data for Trad

        // Traditional
        timingUtility(1000, dimension, 3, 100, false);
    }

    if (flag == 5) { // Finding the crossover!

        string output_file_name = "strassen_crossover_even_dims_v_cutoff.txt";
        ofstream output_file;
        output_file.open(output_file_name);
        
        string next_output_file_name = "strassen_crossover_odd_dims_v_cutoff.txt";
        ofstream next_output_file;
        next_output_file.open(next_output_file_name);

        for (int i = 16; i <= 512; i*=2) {
            output_file << i << "\t" << findCrossover(i, 3) << endl;
            next_output_file << i+1 << "\t" << findCrossover(i+1, 3) << endl;
        }

        output_file.close();
        next_output_file.close();
    }
}
