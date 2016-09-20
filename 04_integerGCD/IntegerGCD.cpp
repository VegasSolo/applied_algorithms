/*
Author:		Dillon VanBuskirk
Date:		6-26-2016
Program:	Integer GCD
Purpose:	Compute the greatest common denominator of two given integers (up to 2^32).
			This will be done using basic Euclidean algorithm, an extended Euclidean
			algorithm, binary gcd algorithm, and k-ary algorithm from Sorenson's paper.
			The goal will be to compare the time between each algorithm for many u,v integer
			pairs as well as find the relation between time and the size of k (which is used in k-ary).
*/

#include "stdafx.h"
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <time.h>
#include <chrono>
#include <fstream>
using namespace std;


// Returns the GCD of two integers, u,v. Typically done through recursion
int gcd_euclidean(int u, int v) {
	unsigned int x;
	while (v)
	{
		x = u % v;
		u = v;
		v = x;
	}
	return u;
}

// Returns the GCD of two integers, u,v. Gets a,b reference parameters which are the coefficients to the integers gcd(u,v)
int gcd_extended_euclidean(int u, int v, int &a, int &b) {
	a = 1;
	int g = u;
	int x = 0;
	int y = v;
	while (y > 0) {
		int q = g / y;
		int t = g % y;
		int s = a - q * x;
		a = x;
		g = y;
		x = s;
		y = t;
	}
	b = (g - u * a) / v;
	return g;
}

// Returns the GCD of two integers, u,v. 
int gcd_binary(int u, int v) {
	int g = 1;
	while (u % 2 == 0 && v % 2 == 0) {
		g = 2 * g;
		u = u / 2;
		v = v / 2;
	}
	while (u != 0 && v != 0) {
		if (u % 2 == 0)
			u = u / 2;
		else if (v % 2 == 0)
			v = v / 2;
		else {
			if (u > v)
				u = abs(u - v) / 2;
			else
				v = abs(u - v) / 2;
		}
	}
	return (g * (u + v));
}

// Returns GCD of two integers by reducing by a factor of k each step. 
int gcd_sorenson_kary(int u, int v, int k) {
	int g = 1;
	int g_u, g_v, a, b;

	while (u != 0 && v != 0) {
		g_u = gcd_binary(u, k);
		g_v = gcd_binary(v, k);
		if (u == v)
			return u;
		else if (g_u > 1) {
			u = u / g_u;
		}
		else if (g_v > 1) {
			v = v / g_v;
		}
		else {			/* Jebelean-Weber Algorithm */			int c = (u / v) % k;
			int f1[2] = { k, 0 };
			int f2[2] = { c, 1 };
			while (f2[0] >= sqrt(k)) {
				f1[0] = f1[0] - floor((f1[0] / f2[0])) * f2[0];
				f1[1] = f1[1] - floor((f1[0] / f2[0])) * f2[1];
				swap(f1, f2);
			}
			/* End Jebelean-Weber Algorithm */
			a = f2[0];
			b = f2[0];
			
			if (u > v) {
				u = abs(a * u + b * v) / k;
			}
			else {
				v = abs(a * u + b * v) / k;
			}
		}
	}
	if (u == 0)
		return v;
	else
		return u;
}


int main()
{
	int u, v; // Container for two integers that are being tested.
	int g; // Container for GCD
	int a, b; // Container for the two coefficients on (u,v)
	int gcd; // Container for the solution
	
	ifstream infile;
	infile.open("numbers.txt");
	ofstream outfile;
	outfile.open("gcd.txt");
	outfile.precision(4);
	outfile.setf(ios::fixed);
	outfile.setf(ios::showpoint);

	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

	while (infile >> u) {
		infile >> v;
		infile >> gcd;
		outfile << "For GCD(" << u << "," << v << "). It should be " << gcd << ":" << endl;
		if (u == 0 || v == 0) {
			outfile << "One of the integers was zero. The GCD is 0. This took 0 microseconds." << endl;
		}
		else {
			if (abs(u + v) < pow(2, 16)) {
				std::chrono::high_resolution_clock::time_point e_start = std::chrono::high_resolution_clock::now();
				g = gcd_euclidean(u, v);
				std::chrono::high_resolution_clock::time_point e_end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> e_time_span = std::chrono::duration_cast<std::chrono::duration<double>>(e_end - e_start);
				outfile << "Euclid returned a gcd of " << g << " and took " << e_time_span.count()*1000000 << " microseconds." << endl;
			}

			std::chrono::high_resolution_clock::time_point ee_start = std::chrono::high_resolution_clock::now();
			g = gcd_extended_euclidean(u, v, a, b);
			std::chrono::high_resolution_clock::time_point ee_end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> ee_time_span = std::chrono::duration_cast<std::chrono::duration<double>>(ee_end - ee_start);
			outfile << "Extended Euclid returned a gcd of " << g << " and (a,b) of (" << a << " ," << b << ") and took " << ee_time_span.count()*1000000 << " microseconds." << endl;

			std::chrono::high_resolution_clock::time_point b_start = std::chrono::high_resolution_clock::now();
			g = gcd_binary(u, v);
			std::chrono::high_resolution_clock::time_point b_end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> b_time_span = std::chrono::duration_cast<std::chrono::duration<double>>(b_end - b_start);
			outfile << "Binary returned a gcd of " << g << " and took " << b_time_span.count()*1000000 << " microseconds." << endl;

			int k[21] = { 2, 3, 4, 6, 11, 19, 28, 31, 42, 73, 74, 200, 209, 1951, 1954, 7919, 64000, 65537, 8420742, 2147483647, sqrt(abs(u + v)) };
			outfile << "Running Sorenson's k-ary with " << size(k) << " values of k:" << endl;
			for (int i = 0; i < size(k); i++) {
				std::chrono::high_resolution_clock::time_point k_start = std::chrono::high_resolution_clock::now();
				g = gcd_sorenson_kary(u, v, k[i]);
				std::chrono::high_resolution_clock::time_point k_end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> k_time_span = std::chrono::duration_cast<std::chrono::duration<double>>(k_end - k_start);
				outfile << "     k value of " << k[i] << " resulted in gcd of " << g << " and took " << k_time_span.count()*1000000 << " microseconds." << endl;
			}
		}
		outfile << "Timing tests are complete for this pair. The GCD should have been " << gcd << "\n" << endl;
	}

	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
	outfile << "-- Complete. -- Total time: " << time_span.count()*1000 << " milliseconds. -- Complete. --" << endl;

	infile.close();
	outfile.close();

    return 0;
}

