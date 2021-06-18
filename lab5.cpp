
#define _CRT_SECURE_NO_WARNINGS
#include "mpi.h"
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>

#define SIZE 16

using namespace std;

void InputMatrix(complex<double>, int);
void OutputMatrix(complex<double>, int);

void MatrixAdd(complex<double> A, complex<double> B, complex<double> Result, int N);
void MatrixSubtrac(complex<double> A, complex<double> B, complex<double> Result, int N);

void StrassenAlgorithm(complex<double> A, complex<double> B, complex<double> C, int N);

void InputMatrix(complex<double> A[][SIZE], int N)
{
	int i, j;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			cin >> A[i][j];
		}
	}
}

void OutputMatrix(complex<double> A[][SIZE], int N)
{
	int i, j;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			cout << A[i][j].real() << " + " << A[i][j].imag() << "i   ";
		}
		cout << endl;
	}
}
void MatrixAdd(complex<double> A[][SIZE], complex<double> B[][SIZE], complex<double> Result[][SIZE], int N)
{
	int i, j;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			Result[i][j] = A[i][j] + B[i][j];
		}
	}

}

void MatrixSubtrac(complex<double> A[][SIZE], complex<double> B[][SIZE], complex<double> Result[][SIZE], int N)
{
	int i, j;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			Result[i][j] = A[i][j] - B[i][j];
		}
	}
}
void StrassenAlgorithm(complex<double> A[][SIZE], complex<double> B[][SIZE], complex<double> C[][SIZE], int N)
{
	if (N == 1)
	{
		C[0][0] = A[0][0] * B[0][0];
		return;
	}
	else
	{
		int Divide = (int)(N / 2);

		complex<double> A11[SIZE][SIZE], A12[SIZE][SIZE], A21[SIZE][SIZE], A22[SIZE][SIZE];
		complex<double> B11[SIZE][SIZE], B12[SIZE][SIZE], B21[SIZE][SIZE], B22[SIZE][SIZE];
		complex<double> C11[SIZE][SIZE], C12[SIZE][SIZE], C21[SIZE][SIZE], C22[SIZE][SIZE];
		complex<double> P1[SIZE][SIZE], P2[SIZE][SIZE], P3[SIZE][SIZE], P4[SIZE][SIZE], P5[SIZE][SIZE], P6[SIZE][SIZE], P7[SIZE][SIZE];
		complex<double> AResult[SIZE][SIZE], BResult[SIZE][SIZE];

		int i, j;

		for (i = 0; i < Divide; i++)
		{
			for (j = 0; j < Divide; j++)
			{
				A11[i][j] = A[i][j];
				A12[i][j] = A[i][j + Divide];
				A21[i][j] = A[i + Divide][j];
				A22[i][j] = A[i + Divide][j + Divide];

				B11[i][j] = B[i][j];
				B12[i][j] = B[i][j + Divide];
				B21[i][j] = B[i + Divide][j];
				B22[i][j] = B[i + Divide][j + Divide];
			}
		}

		MatrixAdd(A11, A22, AResult, Divide);
		MatrixAdd(B11, B22, BResult, Divide);
		StrassenAlgorithm(AResult, BResult, P1, Divide);

		MatrixAdd(A21, A22, AResult, Divide);
		StrassenAlgorithm(AResult, B11, P2, Divide);

		MatrixSubtrac(B12, B22, BResult, Divide);
		StrassenAlgorithm(A11, BResult, P3, Divide);

		MatrixSubtrac(B21, B11, BResult, Divide);
		StrassenAlgorithm(A22, BResult, P4, Divide);

		MatrixAdd(A11, A12, AResult, Divide);
		StrassenAlgorithm(AResult, B22, P5, Divide);

		MatrixSubtrac(A21, A11, AResult, Divide);
		MatrixAdd(B11, B12, BResult, Divide);
		StrassenAlgorithm(AResult, BResult, P6, Divide);

		MatrixSubtrac(A12, A22, AResult, Divide);
		MatrixAdd(B21, B22, BResult, Divide);
		StrassenAlgorithm(AResult, BResult, P7, Divide);

		MatrixAdd(P3, P5, C12, Divide);
		MatrixAdd(P2, P4, C21, Divide);

		MatrixAdd(P1, P4, AResult, Divide);
		MatrixAdd(AResult, P7, BResult, Divide);
		MatrixSubtrac(BResult, P5, C11, Divide);

		MatrixAdd(P1, P3, AResult, Divide);
		MatrixAdd(AResult, P6, BResult, Divide);
		MatrixSubtrac(BResult, P2, C22, Divide);
		for (i = 0; i < Divide; i++)
		{
			for (j = 0; j < Divide; j++)
			{
				C[i][j] = C11[i][j];
				C[i][j + Divide] = C12[i][j];
				C[i + Divide][j] = C21[i][j];
				C[i + Divide][j + Divide] = C22[i][j];
			}
		}

	}

}


int main(int* argc, char** argv)
{
	MPI_Init(argc, &argv);
	int i, j;
	int N, M, Count = 0;
	complex<double> A[SIZE][SIZE], B[SIZE][SIZE], C[SIZE][SIZE];
	int rank, size, newrank, newsize;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Group oldgroup, newgroup1, newgroup2;
	MPI_Comm_group(MPI_COMM_WORLD, &oldgroup);
	int ranks[4] = { 0,1,2,3};
	MPI_Group_incl(oldgroup, 4, ranks, &newgroup1);
	MPI_Comm New1, New2;
	MPI_Comm_create(MPI_COMM_WORLD, newgroup1, &New1);
	if (rank == 0) {
		cout << "Enter the size of the matrices ";
		cin >> N;
	}
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Datatype NewType;
	MPI_Type_vector(1, SIZE, SIZE, MPI_DOUBLE_COMPLEX, &NewType);
	MPI_Type_commit(&NewType);
	M = N;
	if (rank == 0) {
		cout << "Matrix 1:\n";
		InputMatrix(A, M);
		cout << "Matrix 2:\n";
		InputMatrix(B, M);
	}
	MPI_Bcast(&A, 1, NewType, 0, MPI_COMM_WORLD);
	MPI_Bcast(&B, 1, NewType, 0, MPI_COMM_WORLD);
	if (M > 1)
	{
		while (M >= 2)
		{
			M /= 2;
			Count++;
		}

		M = N;

		if (M != (pow(2.0, Count)))
		{
			N = pow(2.0, Count + 1);

			for (i = 0; i < N; i++)
			{
				for (j = 0; j < N; j++)
				{
					if ((i >= M) || (j >= M))
					{
						A[i][j] = 0.0;
						B[i][j] = 0.0;
					}
				}
			}
		}
	}
	if (New1 != MPI_COMM_NULL) {
		if (N == 1)
		{
			C[0][0] = A[0][0] * B[0][0];
		}
		else
		{


			MPI_Comm_size(New1, &newsize);
			MPI_Comm_rank(New1, &newrank);
			int Divide = (int)(N / 2);

			complex<double> A11[SIZE][SIZE], A12[SIZE][SIZE], A21[SIZE][SIZE], A22[SIZE][SIZE];
			complex<double> B11[SIZE][SIZE], B12[SIZE][SIZE], B21[SIZE][SIZE], B22[SIZE][SIZE];
			complex<double> C11[SIZE][SIZE], C12[SIZE][SIZE], C21[SIZE][SIZE], C22[SIZE][SIZE];
			complex<double> P1[SIZE][SIZE], P2[SIZE][SIZE], P3[SIZE][SIZE], P4[SIZE][SIZE], P5[SIZE][SIZE], P6[SIZE][SIZE], P7[SIZE][SIZE];
			complex<double> AResult[SIZE][SIZE], BResult[SIZE][SIZE];

			int i, j;

			for (i = 0; i < 1; i++)
			{
				for (j = 0; j < Divide; j++)
				{
					A11[i][j] = A[i][j];
					A12[i][j] = A[i][j + Divide];
					A21[i][j] = A[i + Divide][j];
					A22[i][j] = A[i + Divide][j + Divide];

					B11[i][j] = B[i][j];
					B12[i][j] = B[i][j + Divide];
					B21[i][j] = B[i + Divide][j];
					B22[i][j] = B[i + Divide][j + Divide];
				}
			}
			MatrixAdd(A11, A22, AResult, Divide);
			MatrixAdd(B11, B22, BResult, Divide);
			StrassenAlgorithm(AResult, BResult, P1, Divide);

			MatrixAdd(A21, A22, AResult, Divide);
			StrassenAlgorithm(AResult, B11, P2, Divide);

			MatrixSubtrac(B12, B22, BResult, Divide);
			StrassenAlgorithm(A11, BResult, P3, Divide);

			MatrixSubtrac(B21, B11, BResult, Divide);
			StrassenAlgorithm(A22, BResult, P4, Divide);

			MatrixAdd(A11, A12, AResult, Divide);
			StrassenAlgorithm(AResult, B22, P5, Divide);

			MatrixSubtrac(A21, A11, AResult, Divide);
			MatrixAdd(B11, B12, BResult, Divide);
			StrassenAlgorithm(AResult, BResult, P6, Divide);

			MatrixSubtrac(A12, A22, AResult, Divide);
			MatrixAdd(B21, B22, BResult, Divide);
			StrassenAlgorithm(AResult, BResult, P7, Divide);
			MatrixAdd(P3, P5, C12, Divide);
			MatrixAdd(P2, P4, C21, Divide);

			MatrixAdd(P1, P4, AResult, Divide);
			MatrixAdd(AResult, P7, BResult, Divide);
			MatrixSubtrac(BResult, P5, C11, Divide);

			MatrixAdd(P1, P3, AResult, Divide);
			MatrixAdd(AResult, P6, BResult, Divide);
			MatrixSubtrac(BResult, P2, C22, Divide);
			for (i = 0; i < Divide; i++)
			{
				for (j = 0; j < Divide; j++)
				{
					C[i][j] = C11[i][j];
					C[i][j + Divide] = C12[i][j];
					C[i + Divide][j] = C21[i][j];
					C[i + Divide][j + Divide] = C22[i][j];
				}
			}
			MPI_Comm_free(&New1);
		}
		if (rank == 0) {
			cout << "Matrix 1:\n\n";
			OutputMatrix(A, M);
			cout << "Matrix 2:\n\n";
			OutputMatrix(B, M);
			cout << "Product:\n\n";
			OutputMatrix(C, M);
		}
	}
	MPI_Type_free(&NewType);
	MPI_Group_free(&newgroup1);
	MPI_Group_free(&oldgroup);
	MPI_Finalize();
	return 0;
}

