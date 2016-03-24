#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <assert.h>
#include <ctime>
#include <random>


/*

 PROGRAM SETUP

 */
using namespace std;
const int CUTOFF = 2;
const bool IN_DEV = true; // Runs a couple simple tests before executing main commands as a sanity check
const string OUTPUT_SEPERATOR = "-----------------------------\n\n";

default_random_engine generator;
uniform_int_distribution<int> distribution(-1,2);


/*

 Definition of Matrix as well as utility functions for generating, populating, and printing them.

 */
// Basic Matrix struct
typedef struct Matrix {
    int dimension;
    vector<vector<int>> entries;
} Matrix;

Matrix* instantiateMatrix(int dimension) {

    Matrix* matrix = new Matrix();

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

 */

// Given matrix, distance from left and distance from top of both parts, and dimension of output matrix, returns sum
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


// Inserts the return value of strassenMult into auxiliary matrix
void updateAuxMatrix(Matrix* P_aux, Matrix* P_new) {

    /*
     We only ever need to look at some square matrix in the
     upper right corner of P_aux, so adjust dimension accordingly
     */
    P_aux->dimension = P_new->dimension;

    for (int i = 0; i < P_new->dimension; i++) {
        for (int j = 0; j < P_new->dimension; i++) {
            P_aux->entries[i][j] = P_new->entries[i][j];
        }
    }
}


void modifySubmatrix(Matrix* C, Matrix* P, int i, int j, bool adding=true) {

    int row, col;

    if (adding) {

        for (row = 0; row < P->dimension; row++) {
            for (col = 0; col < P->dimension; col++) {
                C->entries[i + row][j + col] += P->entries[row][col];
            }
        }

    } else {

        for (row = 0; row < P->dimension; row++) {
            for (col = 0; col < P->dimension; col++) {
                C->entries[i + row][j + col] -= P->entries[row][col];
            }
        }
    }
}


// Hardcoded to deal with Strassen's messiness
void modifyC(Matrix* C, int P_i) {}


// add a row and a col of zeroes to matrix
void pad(Matrix* mat, int dimension) {

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


// TODO dimension shouldn't be passed around-- each matrix should always know it's dimension
Matrix* strassenMult(Matrix* A, Matrix* B) {

    if (A->dimension != B->dimension) {
        throw invalid_argument("strassenMult expects square matrices.");
    }

    int dimension = A->dimension;

	if (dimension <= CUTOFF){
		return tradMult(A, B);
	} else {

        // TODO make even more modular
        // TODO need to pass in references to initial matrices for inline strass

		bool padding = false;

		// If matrix is of odd dimension, pad a row and col of zeroes
		if ((dimension & 1) != 0) {
			padding = true;
			pad(A, dimension);
			pad(B, dimension);
			dimension++; // Keep local representation of dim up to date
		}

        // Matrix to be filled and returned
        Matrix* C = instantiateMatrix(dimension);

        int half_dim = dimension / 2;

		Matrix* m1a = combineSubmatrices(A, 0, 0, half_dim, half_dim, true);
		Matrix* m1b = combineSubmatrices(B, 0, 0, half_dim, half_dim, true);
		Matrix* m1 = strassenMult(m1a, m1b);

        modifySubmatrix(C, m1, 0, half_dim); // Add tr
        modifySubmatrix(C, m1, half_dim, half_dim); // Add br

        delete m1a;
        delete m1b;
        delete m1;

		Matrix* m2a = combineSubmatrices(A, 0, half_dim, half_dim, half_dim, true);
		Matrix* m2b = extractSubmatrix(B, 0, 0);
		Matrix* m2 = strassenMult(m2a, m2b);

        modifySubmatrix(C, m2, 0, 0, false); // Subtract tl
        modifySubmatrix(C, m2, 0, half_dim); // Add tr

        delete m2a;
        delete m2b;
        delete m2;

		Matrix* m3a = extractSubmatrix(A, 0, 0);
		Matrix* m3b = combineSubmatrices(B, half_dim, 0, half_dim, half_dim, false);
		Matrix* m3 = strassenMult(m3a, m3b);

        modifySubmatrix(C, m3, half_dim, 0); // Add bl
        modifySubmatrix(C, m3, half_dim, half_dim, false); // Subtract br

        delete m3a;
        delete m3b;
        delete m3;

		Matrix* m4a = extractSubmatrix(A, half_dim, half_dim);
		Matrix* m4b = combineSubmatrices(B, 0, half_dim, 0, 0, false);
		Matrix* m4 = strassenMult(m4a, m4b);

        modifySubmatrix(C, m4, 0, 0); // Add tl
        modifySubmatrix(C, m4, half_dim, 0); // Add bl

        delete m4a;
        delete m4b;
        delete m4;

		Matrix* m5a = combineSubmatrices(A, 0, 0, half_dim, 0, true);
		Matrix* m5b = extractSubmatrix(B, half_dim, half_dim);
		Matrix* m5 = strassenMult(m5a, m5b);

        modifySubmatrix(C, m5, 0, 0); // Add tl
        modifySubmatrix(C, m5, half_dim, half_dim); // Add br

        delete m5a;
        delete m5b;
        delete m5;

		Matrix* m6a = combineSubmatrices(A, 0, half_dim, 0, 0, false);
		Matrix* m6b = combineSubmatrices(B, 0, 0, half_dim, 0, true);
		Matrix* m6 = strassenMult(m6a, m6b);

        modifySubmatrix(C, m6, 0, 0); // Add tl

        delete m6a;
        delete m6b;
        delete m6;

		Matrix* m7a = combineSubmatrices(A, half_dim, 0, half_dim, half_dim, false);
		Matrix* m7b = combineSubmatrices(B, 0, half_dim, half_dim, half_dim, true);
		Matrix* m7 = strassenMult(m7a, m7b);

        modifySubmatrix(C, m7, half_dim, half_dim, false); // Subtract br

        delete m7a;
        delete m7b;
        delete m7;

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

    Matrix* A;
    Matrix* B;

    if (use_random_matrices) {

        cout << "Testing randomly generated matrices of size: " << dimension << endl;

        A = genRandMatrix(dimension);
        B = genRandMatrix(dimension);

        Matrix* C_strass = strassenMult(A, B);
        Matrix* C_trad = tradMult(A, B);

        delete A;
        delete B;

        /*
         Can't "know" the correct result a priori, so we assume that if each type of
         matrix multiplication is working correcly that they will return the same result
         */
        if (printing_matrix) {
            printMatrix(C_strass, false);
            printMatrix(C_trad, false);
        }
        
        assert(matricesAreEqual(C_strass, C_trad));

        delete C_strass;
        delete C_trad;

    } else {

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

        delete A;
        delete B;

        // Left matrix built from test files
        Matrix* correct_C = buildMatrix(infile, dimension*dimension*2, dimension);

        assert(matricesAreEqual(correct_C, C));

        if (printing_matrix) {
            printMatrix(C);
            printMatrix(correct_C);
        }

        delete C;
        delete correct_C;

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
    Matrix* C = strassenMult(A,B);
    double mult_total = (clock() - mult_start) / (double)(CLOCKS_PER_SEC);
    free(C);

    cout << "Construction for matrix of size " << dimension << " took " << construct_total << "s" << endl;
    cout << "Multiplication took " << mult_total << "s" << endl;
}


void timingUtility(int lower_bound, int upper_bound, int trials, int interval, bool using_strassen=true) {

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
    output_file_name += to_string(lower_bound) + "_" + to_string(upper_bound) + "_" + to_string(interval) + "_" + to_string(trials) + ".txt";
    ofstream output_file;
    output_file.open(output_file_name);

    cur_matrix_dimension = lower_bound;
    while(cur_matrix_dimension <= upper_bound) {

        // For matrix of cur_matrix_dimension = n
        double total_mult_time = 0;
        double avg_mult_time = 0;

        cout << "Matrix of Size: " << cur_matrix_dimension << endl << OUTPUT_SEPERATOR;

        clock_t construct_start = clock();
        Matrix* A = genRandMatrix(cur_matrix_dimension);
        Matrix* B = genRandMatrix(cur_matrix_dimension);
        double construct_time = (clock() - construct_start) / (double)(CLOCKS_PER_SEC);
        construct_time = 0; // No need for this measurement, but might use it later.


        for (int trial = 0; trial < trials; trial++) {

            clock_t mult_start = clock();
            Matrix* C;
            if (using_strassen) {
                C = strassenMult(A,B);
            } else {
                C = tradMult(A,B);
            }

            double mult_total = (clock() - mult_start) / (double)(CLOCKS_PER_SEC);
            total_mult_time += mult_total;

            delete C;
        }

        avg_mult_time = total_mult_time / trials;

        cout << "Average Time for Mult:    " << avg_mult_time << endl;
        output_file << cur_matrix_dimension << "\t" << avg_mult_time << endl;

        cur_matrix_dimension += interval;

        delete A;
        delete B;
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
//        testingUtility("test33.txt", 3, false, true); // Strassen
//        testingUtility("test33.txt", 3, false, false); // Traditional
        testingUtility("", 4, true, true, true); // Random Matrices

//        testingUtility("", 16, true, true); // Random Matrices
//        testingUtility("", 39, true, true); // Random Matrices
        cout << "Basic Tests Pass. Executing instructions from command line." << endl << OUTPUT_SEPERATOR;
    }

    if (flag == 0) { // Production behavior

        if (!IN_DEV) {

            // buildMatrix(filename, read_from_position, buffer_length)
            Matrix* A = buildMatrix(infile, 0, dimension);
            Matrix* B = buildMatrix(infile, dimension*dimension, dimension);
        	Matrix* C = strassenMult(A, B);
            printMatrix(C);
            delete A;
            delete B;
            delete C;
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

    if (flag == 3) { // Generate time data for Strassen

        // Strassen
        timingUtility(dimension, dimension, 1, 1, true);
    }

    if (flag == 4) { // Generate time data for Trad

        // Traditional
        timingUtility(dimension, dimension, 1, 1, false);
    }
}
