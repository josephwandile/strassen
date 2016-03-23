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
const bool IN_DEV = true;
const string OUPUT_SEPERATOR = "-----------------------------\n\n";

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


Matrix* strassenMult(Matrix* mata, Matrix* matb, int dimension) {

	if (dimension <= CUTOFF){
		return tradMult(mata, matb);
	} else {

        // TODO make even more modular
        // TODO need to pass in references to initial matrices for inline strass

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
		return ans_mat;
	}
}

/*

 TESTING

 */

// Assumes left_matrix
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
void testMultiplication(string infile, int dimension, bool using_strassen=true, bool use_random_matrices=true, bool printing_matrix=false) {

    Matrix* A = buildMatrix(infile, 0, dimension);
    Matrix* B = buildMatrix(infile, dimension*dimension, dimension);
    Matrix* C;

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

// TODO Should be able to more precisely define range for testing
// TODO Do this with random matrix
void timingUtility(string infile, int lower_bound, int upper_bound, int trials, bool using_strassen=true) {

    int cur_matrix_dimension;
    string file_name;

    if (using_strassen) {
        cout << "Strassen" << endl;
        file_name = "strassen_";
    } else {
        cout << "Traditional" << endl;
        file_name = "traditional_";
    }

    // Build ouput file name
    file_name += to_string(lower_bound) + "_" + to_string(upper_bound) + "_" + to_string(trials) + ".txt";
    ofstream output_file;
    output_file.open(file_name);

    for (cur_matrix_dimension = lower_bound; cur_matrix_dimension <= upper_bound; cur_matrix_dimension++) {

        double total_construct_time = 0;
        double total_mult_time = 0;
        double avg_construct_time = 0;
        double avg_mult_time = 0;

        cout << "Matrix of Size: " << cur_matrix_dimension << endl << OUPUT_SEPERATOR;

        for (int trial = 0; trial < trials; trial++) {

            clock_t construct_start = clock();
            Matrix* A = buildMatrix(infile, 0, cur_matrix_dimension, true);
            Matrix* B = buildMatrix(infile, cur_matrix_dimension*cur_matrix_dimension, cur_matrix_dimension, false);
            double construct_total = (clock() - construct_start) / (double)(CLOCKS_PER_SEC);
            total_construct_time += construct_total;

            clock_t mult_start = clock();
            Matrix* C;
            if (using_strassen) {
                C = strassenMult(A,B, cur_matrix_dimension);
            } else {
                C = tradMult(A,B);
            }

            double mult_total = (clock() - mult_start) / (double)(CLOCKS_PER_SEC);
            total_mult_time += mult_total;

//            cout << construct_total << "s" << " for construction during trial " << trial << endl;
//            cout << mult_total << "s" << " for multiplication during trial " << trial << endl;
        }

        avg_mult_time = total_mult_time / trials;
        avg_construct_time = total_construct_time / trials;

        cout << "Average Time for Construction:    " << avg_construct_time << endl << "Average Time for Mult:    " << avg_mult_time << endl;
        output_file << cur_matrix_dimension << "\t" << avg_mult_time << endl;
    }

    output_file.close();
}

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
        testMultiplication(infile, dimension, true);
        testMultiplication(infile, dimension, false);
        cout << "Basic Tests Pass. Executing instructions from command line." << endl << OUPUT_SEPERATOR;
    }
    
//    Matrix* T = genRandMatrix(5);
//    printMatrix(T, false);

    if (flag == 0) {

        Matrix* A = buildMatrix(infile, 0, dimension, true);

        /*
         The second parameter determines where to start reading in from the text file.

         The last parameter refers to the fact that this is a "right_matrix"-- basically that it is
         an array of columns instead of rows, which yields improved caching performance during matrix multiplication.
         */
        Matrix* B = buildMatrix(infile, dimension*dimension, dimension, false);

        Matrix* C = tradMult(A,B);
        printMatrix(C);
    }

	if (flag == 1) {

        // Strassen
        testMultiplication(infile, dimension, true);
        return 0;
    }

    if (flag == 2) {

        // Normal
        testMultiplication(infile, dimension, false);
        return 0;
    }

    if (flag == 3) {

        cout << "Testing Cache-optimized Traditional Mult" << endl;
    }

    if (flag == 4) {

        // Strassen
        timingUtility(infile, dimension, dimension, 5);
    }

    if (flag == 5) {

        // Traditional
        timingUtility(infile, dimension, dimension, 5, false);
    }
}
