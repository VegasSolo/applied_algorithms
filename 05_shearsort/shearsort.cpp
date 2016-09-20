/*
Author:		Dillon VanBuskirk & Connor Nesbitt
Date:		7-01-2016
Program:	Shear Sort [Parallel]
Purpose:	Sort data. This program generates data at random and then sorts it and 
			outputs it to a file named sorted.txt. Error checking is coded in this
			but is commented out for cluster runs.
*/

//#include "stdafx.h" // Visual Studio
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
//#include <random>
//#include <vector>
using namespace std;


int compare(const void * elem1, const void * elem2) {
	int f = *((int*)elem1);
	int s = *((int*)elem2);
	if (f > s) return  1;
	if (f < s) return -1;
	return 0;
}

int compare_reversed(const void * elem1, const void * elem2) {
	int f = *((int*)elem1);
	int s = *((int*)elem2);
	if (f < s) return  1;
	if (f > s) return -1;
	return 0;
}

void quicksort_column(int **a, int top, int bottom, int column) {
	int i = top, j = bottom;

	int pivot = bottom;
	swap(a[(i + j) / 2][column], a[bottom][column]);

	while (i < j) {
		while (i < j && a[i][column] < a[pivot][column])
			i++;
		while (i < j && a[j][column] >= a[pivot][column])
			j--;
		swap(a[i][column], a[j][column]);
	}

	swap(a[i][column], a[pivot][column]);
	pivot = i;

	if (top < pivot - 1)
		quicksort_column(a, top, pivot - 1, column);
	if (pivot + 1 < bottom)
		quicksort_column(a, pivot + 1, bottom, column);
}

void shearsort(int **a, const int NUM_DATA, const int NUM_ROWS, const int NUM_COLUMNS, const int THREAD_COUNT) {
	const int SHEARSORT_RUN_COUNT = (int)(ceil(log2(NUM_ROWS)) + 1);

	for (int i = 0; i < SHEARSORT_RUN_COUNT*2; i++) {
		if (i % 2 == 0) {
#pragma omp parallel for num_threads(THREAD_COUNT)
			for (int j = 0; j < NUM_ROWS; j++) {
				if (j % 2 == 0) {
					qsort(a[j], NUM_COLUMNS, sizeof(*a[j]),compare);
				}
				else {
					qsort(a[j], NUM_COLUMNS, sizeof(*a[j]),compare_reversed);
				}
			}
		}
		else {
#pragma omp parallel for num_threads(THREAD_COUNT)
			for (int j = 0; j < NUM_COLUMNS; j++) {
				quicksort_column(a, 0, NUM_ROWS - 1, j);
			}
		}
	}
}

int main(int argc, char* argv[])
{
	const int MAX_DATA = 1000000;
	const int MIN_DATA = 0;
	const int NUM_ROWS = 128; // a multiple of the number of processors used
	const int NUM_COLUMNS = 2500;
	const int THREAD_COUNT = argc > 1 ? strtol(argv[1], NULL, 10) : 16;
	const int NUM_DATA = NUM_ROWS * NUM_COLUMNS;

	int **a;
	a = new int *[NUM_ROWS];
	int seed = 73;
	srand(seed);


	//default_random_engine generator;
	//uniform_int_distribution<unsigned long int> distribution(0, 1000000);

	for (unsigned long int r = 0; r < NUM_ROWS; r++) {
		a[r] = new int[NUM_COLUMNS];
	}

	for (unsigned long int r = 0; r < NUM_ROWS; r++) {
		for (unsigned long int c = 0; c < NUM_COLUMNS; c++) {
			//a[r][c] = distribution(generator);
			a[r][c] = rand() % MAX_DATA;
		}
	}

	/* Unsorted Output */
	/*
	for (unsigned long int r = 0; r < NUM_ROWS - 1; r++) {
	for (unsigned long int c = 0; c < NUM_COLUMNS - 1; c++) {
	cout << a[r][c] << " ";
	}
	cout << endl;
	}
	*/
	/* End Unsorted Output */


	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	shearsort(a, NUM_DATA, NUM_ROWS, NUM_COLUMNS, THREAD_COUNT);
	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);


	/* Output */
	ofstream outfile;
	outfile.open("sorted.txt");
	bool sorted = true;
	int prev = a[0][0];
	vector<vector<int>> failed;

	for (int r = 0; r < NUM_ROWS; r++) {
		if (r % 2 == 0) {
			for (int c = 0; c < NUM_COLUMNS; c++) {
				outfile << a[r][c] << " ";
				if (prev > a[r][c]) {
					sorted = false;
					vector<int> f;
					f.push_back(prev);
					f.push_back(a[r][c]);
					failed.push_back(f);
				}

				prev = a[r][c];
			}
		}
		else {
			for (int c = NUM_COLUMNS - 1; c >= 0; c--) {
				outfile << a[r][c] << " ";
				if (prev > a[r][c]) {
					sorted = false;
					vector<int> f;
					f.push_back(prev);
					f.push_back(a[r][c]);
					failed.push_back(f);
				}
				prev = a[r][c];
			}
		}
		outfile << endl;
	}
	outfile << "\n" << endl;
	outfile << "Shearsort took " << time_span.count() * 1000 << " milliseconds to sort " << NUM_DATA << " data ranging from "
		<< MIN_DATA << "-" << MAX_DATA << " with " << THREAD_COUNT << " threads." << endl;
	if (sorted)
		outfile << "Shearsort was a success." << endl;
	else
		outfile << "Shearsort FAILED." << endl;
	for (unsigned int i = 0; i < failed.size(); i++) {
		for (unsigned int j = 0; j < failed.at(i).size(); j++) {
			outfile << failed.at(i).at(j) << " ";
		}
		outfile << endl;
	}
	outfile.close();
	/* End Output */
	delete[] a;

	return 0;
}
