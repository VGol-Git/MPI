#include "mpi.h"
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

int Program(int argc, char* argv[]);
int* ConventionalMultiplication(int* first_number, int first_number_length,
	int* second_number, int second_number_length, int ProcRank, int ProcNum);
int* StrassenMultiplication(int* first_number, int first_number_length,
	int* second_number, int second_number_length, int ProcRank, int ProcNum);

int main(int argc, char* argv[])
{
	Program(argc, argv); 

	return 0;
}

int Program(int argc, char* argv[])
{
	int ProcNum, ProcRank;
	string num1 = "", num2 = "";

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Datatype MPI_My_Long1;
	MPI_Datatype MPI_My_Long2;
	if (ProcRank == 0)
	{
	//	printf("%d", ProcNum);
		getline(cin, num1);
		getline(cin, num2);

		MPI_Type_contiguous(num1.length(), MPI_CHAR, &MPI_My_Long1);
		MPI_Type_contiguous(num2.length(), MPI_CHAR, &MPI_My_Long2);
		MPI_Type_commit(&MPI_My_Long1);
		MPI_Type_commit(&MPI_My_Long2);


		for (int i = 1; i < ProcNum; i++)
			MPI_Send(num1.c_str(), 1 , MPI_My_Long1, i, 0, MPI_COMM_WORLD);
		//	MPI_Send(num1.c_str(), num1.length(), MPI_CHAR, i, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		for (int i = 1; i < ProcNum; i++)
			MPI_Send(num2.c_str(), 1 , MPI_My_Long2, i, 0, MPI_COMM_WORLD);
		//MPI_Send(num2.c_str(), num2.length(), MPI_CHAR, i, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		//Отправляем числа всем процессорам
	}
	else
	{
		int number_amount;
		MPI_Status status;
	
		MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_CHAR, &number_amount);
		char* number_buf1 = (char*)malloc(sizeof(char) * number_amount);

		MPI_Recv(number_buf1, number_amount, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Barrier(MPI_COMM_WORLD);
		num1 = number_buf1;
		free(number_buf1);
		num1.erase(number_amount, num1.length());


		MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_CHAR, &number_amount);

		char* number_buf2 = (char*)malloc(sizeof(char) * number_amount);
		MPI_Recv(number_buf2, number_amount, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Barrier(MPI_COMM_WORLD);
		num2 = number_buf2;
		free(number_buf2);
		num2.erase(number_amount, num2.length());
		// принимаем от процессоров num1 and num2
	}

	int* first_number = (int*)malloc(sizeof(int) * num1.length());
	int* second_number = (int*)malloc(sizeof(int) * num2.length());
	for (int i = 0; i < num1.length(); i++)
		first_number[i] = atoi(num1.substr(i, 1).c_str());
	for (int i = 0; i < num2.length(); i++)
		second_number[i] = atoi(num2.substr(i, 1).c_str());
	// делаем из числа массив из его цифр

	int* conventionalResult = ConventionalMultiplication(first_number, num1.length(),
		second_number, num2.length(), ProcRank, ProcNum);


	if (ProcRank == 0)
	{
		if (conventionalResult != NULL)
		{
			cout << endl;
			for (int i = 0; i < (num1.length() + num2.length()); i++)
				cout << conventionalResult[i];
			cout << endl;
		}
		else
		{
			cout << "Количество процессов и цифр во втором числе должно быть одинаковым" << endl;
		}
	}

	free(conventionalResult);
	free(first_number);
	free(second_number);

	MPI_Finalize();

	return 0;
}



int* ConventionalMultiplication(int* first_number, int first_number_length,
	int* second_number, int second_number_length, int ProcRank, int ProcNum)
{
	int* result = (int*)malloc(sizeof(int) * (first_number_length + second_number_length));
	for (int i = 0; i < first_number_length + second_number_length; i++)
		result[i] = 0;

	if (second_number_length == ProcNum)
	{
		int i = ProcRank;
		int rest = 0;
		for (int j = first_number_length - 1; j >= 0; j--)
		{
			int intTemp = second_number[i] * first_number[j] + rest;
			if (intTemp < 10)
			{
				if (result[(i + j + 1)] + intTemp < 10)
				{
					result[(i + j + 1)] += intTemp;
				}
				else
				{
					div_t yetAnotherRest = div(result[(i + j + 1)] + intTemp, 10);
					result[(i + j + 1)] = yetAnotherRest.rem;
					result[(i + j)] += yetAnotherRest.quot;
				}
			}
			else
			{
				if (result[(i + j + 1)] + (intTemp % 10) < 10)
				{
					result[(i + j + 1)] += intTemp % 10;
					result[(i + j)] += div(intTemp, 10).quot;
				}
				else
				{
					div_t yetAnotherRest = div(result[(i + j + 1)] + (intTemp % 10), 10);
					result[(i + j + 1)] = yetAnotherRest.rem;
					result[(i + j)] += div(intTemp, 10).quot + yetAnotherRest.quot;
				}
			}
		}

		int* resultsFromAllProcesses = NULL;
		if (ProcRank == 0)
			resultsFromAllProcesses = (int*)malloc(sizeof(int) *
			(first_number_length + second_number_length) * ProcNum);
		MPI_Gather(result, (first_number_length + second_number_length), MPI_INT,
			resultsFromAllProcesses, (first_number_length + second_number_length), MPI_INT,
			0, MPI_COMM_WORLD);

		if (ProcRank == 0)
		{
			int tail;
			for (int i = 0; i < (first_number_length + second_number_length); i++)
			{
				result[i] = 0;
			}
			for (int i = 0; i < (first_number_length + second_number_length); i++)
			{
				tail = 0;
				for (int j = 0; j < (first_number_length + second_number_length) *
					ProcNum; j = j + (first_number_length + second_number_length))
				{
					tail += resultsFromAllProcesses[i + j];
				}
				if (tail < 10)
				{
					result[i] += tail;
				}
				else
				{
					div_t rest = div(tail, 10);
					result[i] = rest.rem;
					int k = i;
					while (result[k - 1] + rest.quot >= 10)
					{
						k--;
						rest = div(result[k] + rest.quot, 10);
						result[k] = rest.rem;
					}
					result[k - 1] += rest.quot;
					tail = 0;
				}
			}
		}

		free(resultsFromAllProcesses);

		return result;
	}
	else
	{
		return NULL;
	}
}
