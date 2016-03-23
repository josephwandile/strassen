#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iomanip>

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
 
 STRASSEN'S METHOD
 
 */

// multiplies matrices using strassen's method
Matrix * strassenmult(Matrix* mata, Matrix* matb, int dimension) {
	
	if (dimension <= CUTOFF){
		return trad_mult(mata, matb);
	} else {
		Matrix * ansmat = instantiate_matrix(dimension);

		// create matrces for "clever multiplication"
		Matrix * m1 = instantiate_matrix(dimension/2);
		Matrix * m2 = instantiate_matrix(dimension/2);
		Matrix * m3 = instantiate_matrix(dimension/2);
		Matrix * m4 = instantiate_matrix(dimension/2);
		Matrix * m5 = instantiate_matrix(dimension/2);
		Matrix * m6 = instantiate_matrix(dimension/2);
		Matrix * m7 = instantiate_matrix(dimension/2);

		// fill m1 = (a11 + a22)(b11+b22)
		Matrix * m1a = instantiate_matrix(dimension/2); // all of these can be reduced to just 
		Matrix * m1b = instantiate_matrix(dimension/2, false); // 2 "working matrices" later
		for (int col = 0; col < dimension/2; col++) {
			for (int row = 0; row < dimension/2; row++) {
				m1a->entries[row][col] = mata->entries[row][col] + 
					mata->entries[row + dimension/2][col + dimension/2];
				m1b->entries[col][row] = matb->entries[col][row] + 
					matb->entries[(col + dimension/2)][row + dimension/2];
			}
		}
		m1 = strassenmult(m1a, m1b, dimension/2);

		// fill m2 = (a21 + a22)b11
		Matrix * m2a = instantiate_matrix(dimension/2); 
		Matrix * m2b = instantiate_matrix(dimension/2, false);
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				m2a->entries[row][col] = mata->entries[row+dimension/2][col] + 
					mata->entries[row + dimension/2][col + dimension/2];
				m2b->entries[col][row] = matb->entries[col][row];
			}
		}
		m2 = strassenmult(m2a, m2b, dimension/2);

		// fill m3 = a11(b12-b22)
		Matrix * m3a = instantiate_matrix(dimension/2);
		Matrix * m3b = instantiate_matrix(dimension/2, false);
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				m3a->entries[row][col] = mata->entries[row][col];
				m3b->entries[col][row] = matb->entries[col + dimension/2][row] -
					matb->entries[(col + dimension / 2)][row + dimension / 2];
			}
		}
		m3 = strassenmult(m3a, m3b, dimension/2);

		// fill m4 = a22(b21-b11)
		Matrix * m4a = instantiate_matrix(dimension/2);
		Matrix * m4b = instantiate_matrix(dimension/2, false);
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				m4a->entries[row][col] = mata->entries[row + dimension / 2][col + dimension / 2];
				m4b->entries[col][row] = matb->entries[col+dimension / 2][row] -
					matb->entries[col][row];
			}
		}
		m4 = strassenmult(m4a, m4b, dimension/2);

		// fill m5 = (a11+a12)b22
		Matrix * m5a = instantiate_matrix(dimension/2);
		Matrix * m5b = instantiate_matrix(dimension/2, false);
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				m5a->entries[row][col] = mata->entries[row][col] +
					mata->entries[row][col + dimension/2];
				m5b->entries[col][row] = matb->entries[col + dimension/2][row + dimension/2];
			}
		}
		m5 = strassenmult(m5a, m5b, dimension / 2);

		// fill m6 = (a21-a11)(b11+b12)
		Matrix * m6a = instantiate_matrix(dimension/2);
		Matrix * m6b = instantiate_matrix(dimension/2, false);
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				m6a->entries[row][col] = mata->entries[row+dimension/2][col] -
					mata->entries[row][col];
				m6b->entries[col][row] = matb->entries[col][row] +
					matb->entries[col][row + dimension/2];
			}
		}
		m6 = strassenmult(m6a, m6b, dimension / 2);

		// fill m7 = (a12-a22)(b21+b22)
		Matrix * m7a = instantiate_matrix(dimension/2);
		Matrix * m7b = instantiate_matrix(dimension/2, false);
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				m7a->entries[row][col] = mata->entries[row][col + dimension/2] -
					mata->entries[row + dimension/2][col + dimension/2];
				m7b->entries[col][row] = matb->entries[col + dimension/2][row] +
					matb->entries[col + dimension/2][row + dimension / 2];
			}
		}
		m7 = strassenmult(m7a, m7b, dimension / 2);
	
		
		// fill answer matrix using found m matrices

		// c11 = m1 + m4 - m5 + m7
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				ansmat->entries[row][col] = m1->entries[row][col] + m4->entries[row][col]
					- m5->entries[row][col] + m7->entries[row][col];
			}
		}

		// c12 = m3 + m5
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				ansmat->entries[row][col + dimension/2] = m3->entries[row][col] + m5->entries[row][col];
			}
		}

		// c21 = m2 + m4
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				ansmat->entries[row + dimension/2][col] = m2->entries[row][col] + m4->entries[row][col];
			}
		}

		// c22 = m1 - m2 + m3 + m6
		for (int col = 0; col < dimension / 2; col++) {
			for (int row = 0; row < dimension / 2; row++) {
				ansmat->entries[row + dimension/2][col + dimension/2] =
					m1->entries[row][col] - m2->entries[row][col]
					+ m3->entries[row][col] + m6->entries[row][col];
			}
		}
		return ansmat;
	}
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
	print_formatted_matrix(strassenmult(A, B, 2));
}
