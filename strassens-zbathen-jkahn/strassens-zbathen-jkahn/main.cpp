#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <assert.h>
#include <ctime>
#include <random>


// TODO func to return two matrices to be multiplied


/*

 PROGRAM SETUP

 */
using namespace std;
const int CUTOFF = 1;
const bool IN_DEV = false;
const string OUTPUT_SEPERATOR = "-----------------------------\n\n";

default_random_engine generator;
uniform_int_distribution<int> distribution(-1,2);

/*

 Definition of Matrix as well as utility functions for generating, populating, and printing them.

 */
// Basic Matrix struct
typedef struct Matrix {
    bool left_matrix;
    int dimension;
    vector<vector<int>> entries;
} Matrix;

Matrix* instantiateMatrix(int dimension, bool left_matrix=true) {

    Matrix* matrix = new Matrix();

    matrix->left_matrix = left_matrix;
    matrix->dimension = dimension;

    vector<vector<int>> entries(dimension, vector<int>(dimension));
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
                    if (matrix->left_matrix) {
                        matrix->entries[row][col] = stoi(element);
                        element = "";
                    } else {
                        matrix->entries[col][row] = stoi(element);
                        element = "";
                    }
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
Matrix* buildMatrix(string infile, int position, int dimension, bool left_matrix=true) {

    Matrix* matrix = instantiateMatrix(dimension, left_matrix);
    populateMatrix(matrix, infile, position, dimension);
    return matrix;
}

// Required format for assignment output
void printMatrix(Matrix* matrix, bool printing_diagonal=true) {

    if (printing_diagonal) {

        cout << "Matrix Diagonal" << endl;
        for (int i = 0; i < matrix->dimension; i++) {
            cout << matrix->entries[i][i] << endl;
        }
        // Trailing newline
        cout << "\n";

    } else {

        bool is_left_matrix = matrix->left_matrix;

        if (is_left_matrix){
            cout << "Formatted Left Matrix" << endl;
        } else {
            cout << "Formatted Right Matrix" << endl;
        }
        int row, col, entry;
        string entry_as_str;
        for (row = 0; row < matrix->dimension; row++) {
            for (col = 0; col < matrix->dimension; col++) {
                if (is_left_matrix) {
                    entry = matrix->entries[row][col];
                    cout << setw(5) << entry;
                } else {
                    entry = matrix->entries[col][row];
                    cout << setw(5) <<  entry;
                }
            }
            cout << endl;
        }
    }
}

/*

 TRADITIONAL MATRIX MULT
 With option for optimized for cache performance
 (not used to avoid unecessary complexity)

 */

// mutiplies matrices mat1 and mat2 traditionally
Matrix* tradMult(Matrix* l_matrix, Matrix* r_matrix) {

    if (l_matrix->dimension != r_matrix->dimension) {
        throw invalid_argument("Function only multiplies square matrices.");
    }

    int dimension = l_matrix->dimension;

    // Always outputs a matrix in standard form (i.e. Matrix[rows][cols])
    Matrix* ans_mat = instantiateMatrix(dimension);

    int i, j, r, cur_dot_prod;

    // If matrix on RHS is in transposed form, optimize cache use
    if (!r_matrix->left_matrix) {
        for (r = 0; r < dimension; r++) {
            for (i = 0; i < dimension; i++) {
                cur_dot_prod = 0;
                for (j = 0; j < dimension; j++) {
                    cur_dot_prod += l_matrix->entries[i][j] * r_matrix->entries[r][j];
                }
                ans_mat->entries[i][r] = cur_dot_prod;
            }
        }

    // Standard matrix multiplication if both are in standard form
    } else {
        for (r = 0; r < dimension; r++) {
            for (i = 0; i < dimension; i++) {
                cur_dot_prod = 0;
                for (j = 0; j < dimension; j++) {
                    cur_dot_prod += l_matrix->entries[i][j] * r_matrix->entries[j][r];
                }
                ans_mat->entries[i][r] = cur_dot_prod;
            }
        }

    }

	return ans_mat;
}

/*

 STRASSEN'S METHOD

 */

// given matrix, distance from left and distance from top of both parts, and dimension of output matrix, returns sum
Matrix* refArith(Matrix* refmat, int leftref1, int downref1, int leftref2, int downref2, int dimension, bool add) {
	Matrix* ans_mat = instantiateMatrix(dimension);
	if (add) {
		for (int row = 0; row < dimension; row++) {
			for (int col = 0; col < dimension; col++) {
				ans_mat->entries[row][col] = refmat->entries[downref1 + row][leftref1 + col] + refmat->entries[downref2 + row][leftref2 + col];
			}
		}
	} else
	{
		for (int row = 0; row < dimension; row++) {
			for (int col = 0; col < dimension; col++) {
				ans_mat->entries[row][col] = refmat->entries[downref1 + row][leftref1 + col] - refmat->entries[downref2 + row][leftref2 + col];
			}
		}
	}
	return ans_mat;
}


Matrix* refReturn(Matrix* refmat, int leftref, int downref, int dimension) {
	Matrix* ans_mat = instantiateMatrix(dimension);
	for (int row = 0; row < dimension; row++) {
		for (int col = 0; col < dimension; col++) {
			ans_mat->entries[row][col] = refmat->entries[downref + row][leftref + col];
		}
	}
	return ans_mat;
}

// add a row and a col of zeroes to matrix
void pad(Matrix* mat, int dimension) {

	// pad a row at the bottom
	vector<int> zeroes(dimension + 1, 0);
	mat->entries.push_back(zeroes);

	// pad a column on the right
	for (int i = 0; i < dimension; i++) {
		mat->entries[i].push_back(0);
	}
}

// remove the last row and last col of a matrix
void unpad(Matrix* matrix, int dimension) {

	// remove last row from bottom
	matrix->entries.pop_back();

	// remove last col from right
	for (int i = 0; i < dimension; i++) {
		matrix->entries[i].pop_back();
	}
	matrix->dimension = dimension - 1;
}

Matrix* strassenMult(Matrix* mata, Matrix* matb, int dimension) {

	if (dimension <= CUTOFF){
		return tradMult(mata, matb);
	} else {

        // TODO make even more modular
        // TODO need to pass in references to initial matrices for inline strass

		bool padding = false;
		
		// if matrix is of odd dimension, pad a row and col of zeroes
		if (dimension%2 == 1) {
			padding = true;
			pad(mata, dimension);
			pad(matb, dimension);
			dimension += 1;
		} 

		Matrix* m1a = refArith(mata, 0, 0, dimension / 2, dimension / 2, dimension / 2, true);
		Matrix* m1b = refArith(matb, 0, 0, dimension / 2, dimension / 2, dimension / 2, true);
		Matrix* m1 = strassenMult(m1a, m1b, dimension / 2);

		Matrix* m2a = refArith(mata, 0, dimension / 2, dimension / 2, dimension / 2, dimension / 2, true);
		Matrix* m2b = refReturn(matb, 0, 0, dimension / 2);
		Matrix* m2 = strassenMult(m2a, m2b, dimension / 2);

		Matrix* m3a = refReturn(mata, 0, 0, dimension / 2);
		Matrix* m3b = refArith(matb, dimension / 2, 0, dimension / 2, dimension / 2, dimension / 2, false);
		Matrix* m3 = strassenMult(m3a, m3b, dimension / 2);

		Matrix* m4a = refReturn(mata, dimension / 2, dimension / 2, dimension / 2);
		Matrix* m4b = refArith(matb, 0, dimension / 2, 0, 0, dimension / 2, false);
		Matrix* m4 = strassenMult(m4a, m4b, dimension / 2);

		Matrix* m5a = refArith(mata, 0, 0, dimension / 2, 0, dimension / 2, true);
		Matrix* m5b = refReturn(matb, dimension / 2, dimension / 2, dimension / 2);
		Matrix* m5 = strassenMult(m5a, m5b, dimension / 2);

		Matrix* m6a = refArith(mata, 0, dimension / 2, 0, 0, dimension / 2, false);
		Matrix* m6b = refArith(matb, 0, 0, dimension / 2, 0, dimension / 2, true);
		Matrix* m6 = strassenMult(m6a, m6b, dimension / 2);

		Matrix* m7a = refArith(mata, dimension / 2, 0, dimension / 2, dimension / 2, dimension / 2, false);
		Matrix* m7b = refArith(matb, 0, dimension / 2, dimension / 2, dimension / 2, dimension / 2, true);
		Matrix* m7 = strassenMult(m7a, m7b, dimension / 2);

		// fill answer matrix using found m matrices
		Matrix* ans_mat = instantiateMatrix(dimension);

		// c11 = m1 + m4 - m5 + m7
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				ans_mat->entries[row][col] = m1->entries[row][col] + m4->entries[row][col]
					- m5->entries[row][col] + m7->entries[row][col];
			}
		}

		// c12 = m3 + m5
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				ans_mat->entries[row][col + dimension/2] = m3->entries[row][col] + m5->entries[row][col];
			}
		}

		// c21 = m2 + m4
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				ans_mat->entries[row + dimension/2][col] = m2->entries[row][col] + m4->entries[row][col];
			}
		}

		// c22 = m1 - m2 + m3 + m6
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				ans_mat->entries[row + dimension/2][col + dimension/2] =
					m1->entries[row][col] - m2->entries[row][col]
					+ m3->entries[row][col] + m6->entries[row][col];
			}
		}
		
		// if we padded at the beginning, remove the extra zeroes
		if (padding) {
			unpad(ans_mat, dimension);
		}

		return ans_mat;
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

// TODO extend to allow for random matrix testing
void testingUtility(string infile, int dimension, bool use_random_matrices=true, bool using_strassen=true, bool printing_matrix=false) {

    Matrix* A;
    Matrix* B;

    if (use_random_matrices) {

        cout << "Testing randomly generated matrices of size: " << dimension << endl << OUTPUT_SEPERATOR;

        A = genRandMatrix(dimension);
        B = genRandMatrix(dimension);

        Matrix* C_strass = strassenMult(A, B, dimension);
        Matrix* C_trad = tradMult(A, B);

        /*
         Can't "know" the correct result a priori, so we assume that if each type of
         matrix multiplication is working correcly that they will return the same result
         */
        assert(matricesAreEqual(C_strass, C_trad));

        if (printing_matrix) {
            printMatrix(C_strass);
            printMatrix(C_trad);
        }

    } else {

        A = buildMatrix(infile, 0, dimension);
        B = buildMatrix(infile, dimension*dimension, dimension);
        Matrix* C;

        // Deterministic test so we can test Strassen and Trad independently of each other
        if (using_strassen) {
            cout << "Testing Strassen" << endl;
            C = strassenMult(A,B, dimension);
        } else {
            cout << "Testing Traditional Mult" << endl;
            C = tradMult(A,B);
        }

        // Left matrix built from test files
        Matrix* correct_C = buildMatrix(infile, dimension*dimension*2, dimension);

        assert(matricesAreEqual(correct_C, C));

        if (printing_matrix) {
            printMatrix(C);
            printMatrix(correct_C);
        }

    }
}

/*

 TIMING UTILITY FUNCTIONS
 
 timeMatrixFromFile exists in case TFs/we want to time the execution of Strassen's on a known matrix pair
 
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
    Matrix* C = strassenMult(A,B, dimension);
    double mult_total = (clock() - mult_start) / (double)(CLOCKS_PER_SEC);

    cout << "Construction for matrix of size " << dimension << " took " << construct_total << "s" << endl;
    cout << "Multiplication took " << mult_total << "s" << endl;
}


void timingUtility(int lower_bound, int upper_bound, int trials, bool using_strassen=true) {

    int cur_matrix_dimension;
    string output_file_name;

    if (using_strassen) {
        cout << "Timing Strassen" << endl;
        output_file_name = "strassen_";
    } else {
        cout << "Timing Traditional" << endl;
        output_file_name = "traditional_";
    }

    // Build ouput file name
    output_file_name += to_string(lower_bound) + "_" + to_string(upper_bound) + "_" + to_string(trials) + ".txt";
    ofstream output_file;
    output_file.open(output_file_name);

    for (cur_matrix_dimension = lower_bound; cur_matrix_dimension <= upper_bound; cur_matrix_dimension++) {

        // For matrix of cur_matrix_dimension = n
        double total_mult_time = 0;
        double avg_mult_time = 0;

        cout << "Matrix of Size: " << cur_matrix_dimension << endl << OUTPUT_SEPERATOR;

        // TODO Currently not doing anything with this construction data. Might be useful later.
//        clock_t construct_start = clock();
        Matrix* A = genRandMatrix(cur_matrix_dimension); // TODO Free this? Memory leaks potentially from recreating matrices again and again?
        Matrix* B = genRandMatrix(cur_matrix_dimension);
//        double construct_time = (clock() - construct_start) / (double)(CLOCKS_PER_SEC);

        for (int trial = 0; trial < trials; trial++) {

            clock_t mult_start = clock();
            Matrix* C; // TODO Free?
            if (using_strassen) {
                C = strassenMult(A,B, cur_matrix_dimension);
            } else {
                C = tradMult(A,B);
            }

            double mult_total = (clock() - mult_start) / (double)(CLOCKS_PER_SEC);
            total_mult_time += mult_total;
        }

        avg_mult_time = total_mult_time / trials;

        cout << "Average Time for Mult:    " << avg_mult_time << endl;
        output_file << cur_matrix_dimension << "\t" << avg_mult_time << endl;
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
        testingUtility(infile, dimension, false, true);
        testingUtility(infile, dimension, false, false);
        cout << "Basic Tests Pass. Executing instructions from command line." << endl << OUTPUT_SEPERATOR;
    }

    if (flag == 0) { // Production settings

        Matrix* A = buildMatrix(infile, 0, dimension, true);

        /*
         The second parameter determines where to start reading in from the text file.

         The last parameter refers to the fact that this is a "right_matrix"-- basically that it is
         an array of columns instead of rows, which yields improved caching performance during matrix multiplication.
         */
        Matrix* B = buildMatrix(infile, dimension*dimension, dimension, false);

		Matrix* C = tradMult(A,B);
        printMatrix(C,false);
    }

	if (flag == 1) { // Testing on randomly generated matrices of arbitrary size

        testingUtility(infile, dimension, true,true,true);
        return 0;
    }

    if (flag == 2) { // Run deterministic tests on Strassen and Trad

        // Normal
        testingUtility(infile, dimension, false, true,true);
        testingUtility(infile, dimension, false, false,true);
        return 0;
    }

    if (flag == 3) { // Generate time data for Strassen

        // Strassen
        timingUtility(dimension, dimension, 5);
    }

    if (flag == 4) { // Generate time data for Trad

        // Traditional
        timingUtility(dimension, dimension, 5, false);
    }
}
