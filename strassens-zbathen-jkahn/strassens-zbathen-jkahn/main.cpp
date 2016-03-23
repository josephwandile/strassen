#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;
const int CUTOFF = 1;

/*
 
 Definition of Matrix as well as utility functions for generating, populating, and printing them.
 
 */

// Basic Matrix struct
typedef struct Matrix {
    bool left_matrix;
    int dimension;
    vector<vector<int>> entries;
} Matrix;

// Memory allocation
Matrix* instantiate_matrix(int dimension, bool left_matrix=true) {

    Matrix* matrix = new Matrix();
    matrix->left_matrix = left_matrix;
    matrix->dimension = dimension;
    vector<vector<int>> entries(dimension, vector<int>(dimension));
    matrix->entries = entries;
    
    return matrix;
}

// Does the actual IO work for populating matrices
void populate_matrix(Matrix* matrix, string infile, int position, int dimension) {
    
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
Matrix* build_matrix(string infile, int position, int dimension, bool left_matrix=true) {
    
    Matrix* matrix = instantiate_matrix(dimension, left_matrix);
    populate_matrix(matrix, infile, position, dimension);
    return matrix;
}

// Required format for assignment output
void print_diagonal(Matrix* matrix) {
    cout << "Printing Matrix Diagonal" << endl;
    for (int i = 0; i < matrix->dimension; i++) {
        cout << matrix->entries[i][i] << endl;
    }
}

// Print in grid form
void print_formatted_matrix(Matrix* matrix) {
    
    bool is_left_matrix = matrix->left_matrix;
    
    if (is_left_matrix){
        cout << "Printing Left Matrix" << endl;
    } else {
        cout << "Printing Right Matrix" << endl;
    }
    int row, col;
    for (row = 0; row < matrix->dimension; row++) {
        for (col = 0; col < matrix->dimension; col++) {
            if (is_left_matrix) {
                cout << matrix->entries[row][col] << "  ";
            } else {
                cout << matrix->entries[col][row] << "  ";
            }
        }
        cout << endl;
    }
}

/* 
 
 TRADITIONAL MATRIX MULT
 
 Optimized for cache performance.
 
 */

// mutiplies matrices mat1 and mat2 traditionally
Matrix* trad_mult(Matrix* l_matrix, Matrix* r_matrix) {
    
    if (l_matrix->left_matrix && r_matrix->left_matrix) {
        throw invalid_argument("Can't multiply two left matrices.");
    }
    
    if (l_matrix->dimension != r_matrix->dimension) {
        throw invalid_argument("Function only multiplies square matrices.");
    }
    
    int dimension = l_matrix->dimension;

    // TODO More efficient to index into a right_matrix, but using left_matrix for now.
    Matrix* ansmat = instantiate_matrix(dimension, true);
    
    int i, j, r, cur_dot_prod;
    
    for (r = 0; r < dimension; r++) {
        
        for (i = 0; i < dimension; i++) {
            
            cur_dot_prod = 0;
            
            for (j = 0; j < dimension; j++) {
                
                cur_dot_prod += l_matrix->entries[i][j] * r_matrix->entries[r][j];
            }
        
            // TODO Assumes left_matrix to be ouput. Invert if right_matrix
            ansmat->entries[i][r] = cur_dot_prod;
        }
    }
	return ansmat;
}

/*
 
 PROGRAM INTERFACE 
 
 */

int main(int argc, char* argv[])
{
	// check to ensure correct number of arguments
	if (argc != 4) {
        cout << "Please enter the correct number of arguments" << endl;
		return -1;
	}
	
    int flag = stoi(argv[1]);
    int dimension = stoi(argv[2]);
    string infile = argv[3];
    
	// TEST CODE
	if (flag == 1) {
        cout << "Run test suite 1" << endl;
        
        // TODO Write test function for matrix equivalency
	}
    
    Matrix* A = build_matrix(infile, 0, dimension, true);
    
    /*
     The second parameter determines where to start reading in from the text file.
     
     The last parameter refers to the fact that this is a "right_matrix"-- basically that it is
     an array of columns instead of rows, which yields improved caching performance during matrix multiplication.
     */
    Matrix* B = build_matrix(infile, dimension*dimension, dimension, false);

    Matrix* C = trad_mult(A,B);
    print_formatted_matrix(C);
}






