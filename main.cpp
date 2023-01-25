#include <iostream>
#include <ctime>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex>
#include <fstream> 


using namespace std;

bool is_exp_of_2(int n)
{
	return (n & (n - 1)) == 0;
}

int is_file_correct(const char* name, int n) {
	if (is_exp_of_2(n) == 1) {
		ifstream file1(name);
		if (!file1.is_open()) cout << "Error to open file\n";
		file1.close();
		ifstream file2(name);
		double k = 0.0; int count = 0;
		while (file2 >> k)
			count++;
		file2.close();
		if (count == 2 * n)
			return 1;
		else return -1;
	}
	else {
		cout << "N is not power of 2" << endl;
		return -1;
	}
}

complex<double>* FillArray(complex<double> arr[], const int size)
{
	ifstream file("your_file");
	double* x = new double[size];
	double* y = new double[size];
	for (int i = 0; i < size; i++) 
	{
		file >> x[i];
		file >> y[i];
		arr[i] = complex<double>(x[i], y[i]);
		cout << arr[i] << endl;
	}
	return arr;
}

void PrintArray(complex<double> arr[], const int size)
{
	for (int i = 0; i < size; i++)
	{
		cout << arr[i] << " ";
	}
}

complex<double>* DPF(complex<double> const X[], const int size)
{
	complex <double> im(0.0, 1.0); // im - мнимая единица
	complex <double> *Y = new complex<double>(size);
	for (int i = 0; i < size; i++) 
	{
		Y[i] = 0;
	}
	for (int k = 0; k < size; k++)
	{
		for (int j = 0; j < size; j++)
		{
			Y[k] += 1/sqrt(size) * X[j] * exp( (-2 * k / size * j) * M_PI * im);
		}
	}
	return Y;
	//delete[] Y;
}


int main()
{
	srand(time(0));

	ifstream info("your_file");
	double k = 0.0;
	int count = 0;
	while (info >> k)
		count++; //кол-во элементов в файле
	info.close();
	int N = count / 2; //кол-во пар(x+iy)
	complex<double>* X = new complex<double>[N];
	cout << "Array received from the file: " << endl;
	FillArray(X,N);
	complex<double> *Y = DPF(X, N);
	cout << "Converted array: " << endl;
	PrintArray(Y, N);
	delete[] X;
}
