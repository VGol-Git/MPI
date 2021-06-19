#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "mpi.h"
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <vector>
using namespace std;
int noOfDigit(long x) {
	int n = 0;
	while (x > 0) {
		x /= 10;
		n++;
	}
	return n;
}
void schonhageStrassenMultiplication(long a, long b, int n, int m) {
	int* linearConvolution = new int[n + m - 1];
	for (int i = 0; i < (n + m - 1); i++)
		linearConvolution[i] = 0;
	long int p = a;
	for (int i = 0; i < m; i++) {
		a = p;
		for (int j = 0; j < n; j++) {
			linearConvolution[i + j] += (b % 10) * (a % 10);
			a /= 10;
		}
		b /= 10;
	}
	cout << "Convolutions: \n";
	for (int i = (n + m - 2); i >= 0; i--) {
		cout << linearConvolution[i] << " ";
	}


	long int product = 0;
	vector<long int> prod(n + m,0);
	int nextCarry = 0, base = 1;
	for (int i = 0; i < n + m - 1; i++) {
		linearConvolution[i] += nextCarry;
		product = product + (base * (linearConvolution[i] % 10));
		prod[i] = prod[i] + (base * (linearConvolution[i] % 10));
		nextCarry = linearConvolution[i] / 10;
		//base *= 10;
	}
	if (nextCarry)
		prod[n + m - 1] = nextCarry;

	cout << "\nresult: \n" << endl;
	for (int i = 0; i < n + m; i++)
		cout  << prod[i];

	cout << "\nresult: \n" << endl;
	for (int i = n + m-1 ; i >= 0; i--)
		cout  << prod[i];

	//cout << "\nresult: \n" << product << endl;
}
int main(int* argc, char** argv)
{
	MPI_Init(argc, &argv);
	int index[2] = { 1, 2 };
	int edges[2] = { 1, 0 };
	int size, rank, reorder = 1;
	long int a, b;
	MPI_Comm graphcomm;
	MPI_Graph_create(MPI_COMM_WORLD, 2, index, edges, reorder, &graphcomm);
	int nneighbors, neighbors;
	if (graphcomm != MPI_COMM_NULL) {
		MPI_Comm_size(graphcomm, &size);
		MPI_Comm_rank(graphcomm, &rank);
		MPI_Graph_neighbors_count(graphcomm, rank, &nneighbors);
		cout << "Count of Process " << rank << " neighbors: " << nneighbors << endl;
		MPI_Graph_neighbors(graphcomm, rank, 1, &neighbors);
		cout << "Rank of Process " << rank << " neighbors: " << neighbors << endl;
		if (rank == 0) {
			cout << "Enter the numbers:\n";
			cin >> a;
			cin >> b;
		}
		MPI_Bcast(&a, 1, MPI_LONG_INT, 0, graphcomm);
		MPI_Bcast(&b, 1, MPI_LONG_INT, 0, graphcomm);
		if (rank == 1) {
			int n = noOfDigit(a);
			int m = noOfDigit(b);
			
			schonhageStrassenMultiplication(a, b, n, m);
		}
	}
	MPI_Finalize();
}
