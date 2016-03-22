#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;
const int CUTOFF = 1;

/*
 
 Definition of Matrix as well as utility functions for generating, populating, and printing them.
 
 */

// Basic Matrix struct
typedef struct Matrix {
    int dimension;
    vector<vector<int>> entries;
} Matrix;

// Memory allocation
Matrix* instantiate_matrix(int dimension) {
    
    Matrix* matrix = new Matrix();
    matrix->dimension = dimension;
    vector<vector<int>> entries(dimension, vector<int>(dimension));
    matrix->entries = entries;
    
    return matrix;
}

// Does the actual IO work for populating matrices
void populate_matrix(Matrix* matrix, string infile, int position, int dimension, bool left_matrix=true) {
    
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
                    if (left_matrix) {
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
    
    Matrix* matrix = instantiate_matrix(dimension);
    populate_matrix(matrix, infile, position, dimension, left_matrix);
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
void print_formatted_matrix(Matrix* matrix, bool left_matrix=true) {
    cout << "Printing Matrix" << endl;
    int row, col;
    for (row = 0; row < matrix->dimension; row++) {
        for (col = 0; col < matrix->dimension; col++) {
            if (left_matrix) {
                cout << matrix->entries[row][col] << "  ";
            } else {
                cout << matrix->entries[col][row] << "  ";
            }
        }
        cout << endl;
    }
}

// Print each matrix element, a(0,0), a(0,1), ..., a(dimension,dimension) on new line
void print_whole_matrix(Matrix* matrix) {
    cout << "Printing matrix" << endl;
	int row,col;
	for (row = 0; row < matrix->dimension; row++) {
		for (col = 0; col < matrix->dimension; col++) {
			cout << matrix->entries[row][col] << endl;
		}
	}
}

/*
 
 UTILITY FUNCTIONS FOR MATRIX MULT
 
 */

// calculates the dot product of row r of mat1 and col c of mat2
int dotprod(int* mat1, int*mat2, int r, int c, int dimension) {
	int answer = 0;
	for (int i = 0; i < dimension; i++) {
		answer += (mat1[r*dimension + i] * mat2[i*dimension + c]);
	}
	return answer;
}

/* 
 
 TRADITIONAL MATRIX MULT
 
 */

// mutiplies matrices mat1 and mat2 traditionally
int * tradmult(int* mat1, int*mat2, int dimension) {
	int * ansmat = new int[dimension*dimension];
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			ansmat[i*dimension + j] = dotprod(mat1, mat2,i,j, dimension);
		}
	}
	return ansmat;
}

/* 
 
 STRASSEN'S METHOD
 
 */

// multiplies matrices using strassen's method
int * strassenmult(int* mata, int*matb, int dimension) {
	
	if (dimension <= CUTOFF){
		return tradmult(mata, matb, dimension);
	}
	else {
		int * strassans = new int[dimension*dimension];

		// create matrces for "clever multiplication"
		int * m1 = new int[dimension*dimension / 4];
		int * m2 = new int[dimension*dimension / 4];
		int * m3 = new int[dimension*dimension / 4];
		int * m4 = new int[dimension*dimension / 4];
		int * m5 = new int[dimension*dimension / 4];
		int * m6 = new int[dimension*dimension / 4];
		int * m7 = new int[dimension*dimension / 4];

		// fill m1 = (a11 + a22)(b11+b22)
		int * m1a = new int[dimension*dimension/4]; // all of these can be reduced to just 
		int * m1b = new int[dimension*dimension/4]; // 2 "working matrices" later
		for (int c = 0; c < dimension/2; c++) {
			for (int r = 0; r < dimension/2; r++) {
				m1a[r*dimension/2 + c] = mata[r*dimension + c] + 
					mata[(r + dimension/2)*dimension + c + dimension/2];
				m1b[r*dimension/2 + c] = matb[r*dimension + c] + 
					matb[(r + dimension/2)*dimension + c + dimension/2];
			}
		}
		m1 = strassenmult(m1a, m1b, dimension/2);

		// fill m2 = (a21+a22)b11
		int * m2a = new int[dimension*dimension/4];
		int * m2b = new int[dimension*dimension/4];
		for (int c = 0; c < dimension/2; c++) {
			for (int r = 0; r < dimension/2; r++) {
				m2a[r*dimension/2 + c] = mata[(r+dimension/2)*dimension + c] + 
					mata[(r + dimension / 2)*dimension + c + dimension / 2];
				m2b[r*dimension / 2 + c] = matb[r*dimension + c];
			}
		}
		m2 = strassenmult(m2a, m2b, dimension/2);

		// fill m3 = a11(b12-b22)
		int * m3a = new int[dimension*dimension/4];
		int * m3b = new int[dimension*dimension/4];
		for (int c = 0; c < dimension/2; c++) {
			for (int r = 0; r < dimension / 2; r++) {
				m3a[r*dimension/2 + c] = mata[r*dimension + c];
				m3b[r*dimension/2 + c] = matb[r*dimension + c + dimension/2] -
					matb[(r + dimension / 2)*dimension + c + dimension / 2];
			}
		}
		m3 = strassenmult(m3a, m3b, dimension/2);

		// fill m4 = a22(b21-b11)
		int * m4a = new int[dimension*dimension/4];
		int * m4b = new int[dimension*dimension/4];
		for (int c = 0; c < dimension/2; c++) {
			for (int r = 0; r < dimension / 2; r++) {
				m4a[r*dimension/2 + c] = mata[(r + dimension / 2)*dimension + c + dimension / 2];
				m4b[r*dimension/2 + c] = matb[(r+dimension/2)*dimension + c] -
					matb[r*dimension + c];
			}
		}
		m4 = strassenmult(m4a, m4b, dimension/2);

		// fill m5 = (a11+a12)b22
		int * m5a = new int[dimension*dimension/4];
		int * m5b = new int[dimension*dimension/4];
		for (int c = 0; c < dimension/2; c++) {
			for (int r = 0; r < dimension/2; r++) {
				m5a[r*dimension/2 + c] = mata[r*dimension + c] + 
					mata[r*dimension + c + dimension/2];
				m5b[r*dimension/2 + c] = matb[(r + dimension/2)*dimension + c + dimension/2];
			}
		}
		m5 = strassenmult(m5a, m5b, dimension / 2);

		// fill m6 = (a21-a11)(b11+b12)
		int * m6a = new int[dimension*dimension/4];
		int * m6b = new int[dimension*dimension/4];
		for (int c = 0; c < dimension/2; c++) {
			for (int r = 0; r < dimension/2; r++) {
				m6a[r*dimension/2 + c] = mata[(r+dimension/2)*dimension + c] - 
					mata[r*dimension + c];
				m6b[r*dimension/2 + c] = matb[r*dimension + c] + 
					matb[r*dimension + c + dimension/2];
			}
		}
		m6 = strassenmult(m6a, m6b, dimension / 2);

		// fill m7 = (a12-a22)(b21+b22)
		int * m7a = new int[dimension*dimension/4];
		int * m7b = new int[dimension*dimension/4];
		for (int c = 0; c < dimension/2; c++) {
			for (int r = 0; r < dimension/2; r++) {
				m7a[r*dimension/2 + c] = mata[r*dimension + c + dimension/2] - 
					mata[(r + dimension/2)*dimension + c + dimension/2];
				m7b[r*dimension/2 + c] = matb[(r + dimension/2)*dimension + c] + 
					matb[(r + dimension/2)*dimension + c + dimension / 2];
			}
		}
		m7 = strassenmult(m7a, m7b, dimension / 2);
	
		
		// fill answer matrix using found m matrices

		// c11 = m1 + m4 - m5 + m7
		for (int c = 0; c < dimension / 2; c++) {
			for (int r = 0; r < dimension / 2; r++) {
				strassans[r*dimension + c] = m1[r*dimension/2 + c] + m4[r*dimension / 2 + c] 
					- m5[r*dimension / 2 + c] + m7[r*dimension / 2 + c];
			}
		}

		// c12 = m3 + m5
		for (int c = 0; c < dimension / 2; c++) {
			for (int r = 0; r < dimension / 2; r++) {
				strassans[r*dimension + c + dimension/2] = m3[r*dimension / 2 + c] + m5[r*dimension / 2 + c];
			}
		}

		// c21 = m2 + m4
		for (int c = 0; c < dimension / 2; c++) {
			for (int r = 0; r < dimension / 2; r++) {
				strassans[(r+dimension/2)*dimension + c] = m2[r*dimension / 2 + c] + m4[r*dimension / 2 + c];
			}
		}

		// c22 = m1 - m2 + m3 + m6
		for (int c = 0; c < dimension / 2; c++) {
			for (int r = 0; r < dimension / 2; r++) {
				strassans[(r+dimension/2)*dimension + c + dimension/2] = 
					m1[r*dimension / 2 + c] - m2[r*dimension / 2 + c]
					+ m3[r*dimension / 2 + c] + m6[r*dimension / 2 + c];
			}
		}
		return strassans;
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
	}
    
    Matrix* A = build_matrix(infile, 0, dimension, true);
    
    /*
     The second parameter determines where to start reading in from the text file.
     
     The last parameter refers to the fact that this is a "right_matrix"-- basically that it is
     an array of columns instead of rows, which yields improved caching performance during matrix multiplication.
     */
    Matrix* B = build_matrix(infile, dimension*dimension, dimension, false);
    
    print_formatted_matrix(A);
    print_formatted_matrix(B);
}