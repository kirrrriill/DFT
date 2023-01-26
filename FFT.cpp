#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <complex>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream> 
#include <string>
#include <chrono>
using namespace std;

bool is_pow_of_2(int n)//функция для проверки, является ли количество элементов точной степенью двойки
{
	return (n & (n - 1)) == 0;
}

string dec2bin(int num, int n) //функция для перевода числа из десятичного в двоичное в виде строки
{
	char res1[20];
	string str;
	str = _itoa(num,res1,2);
	for (size_t i = str.length(); i < n ;i++) {
		str = "0"+str;
	}
	return str;
}

int bin2dec(int n)  //функция для преобразования двоичного числа в десятичное
{
	int num = n;
	int dec_value = 0; 
	int base = 1;
	int temp = num;
	while (temp) {
		int last_digit = temp % 10;
		temp = temp / 10;
		dec_value += last_digit * base;
		base = base * 2;
	}
	return dec_value;
}

int is_file_correct(const char* name1) //функция для проверки того, что с файлом все хорошо
{
	ifstream file(name1);
	if (!file.is_open()) cout << "File isn`t available or doesn`t exist\n"; // если не открылся
	else if (file.peek() == EOF) cout << "File is empty\n"; // если первый символ конец файла
	file.close();
	ifstream fin(name1);
	double k = 0.0; int count = 0;
	while (fin >> k)
		count++;
	fin.close();
	return count;
}

complex<double>* inv(int N, int n, complex<double>* Z1) //функция для выполнения двоично-инверсных перестановок
{
	int* A = new int[N];
	for (int k = 0; k < N; k++) {
		A[k] = 0;
	}
	int i = 0;
	while (i < N) {
		string j = dec2bin(i, n);
		size_t length1 = j.length();
		string j1(length1, '0');
		for (size_t k = 0; k < j.size(); k++) {
			j1[k] = j[j.size() - k - 1];
		}
		int n1 = atoi(j.c_str());
		int n2 = atoi(j1.c_str());
		int n_1 = bin2dec(n1);
		int n_2 = bin2dec(n2);
		A[i] = n_1; int h = 0;
		for (int k = 0; k < i; k++) {
			if (n_2 == A[k])
				h++;
		}
		if (h == 0) {
			swap(Z1[n_1], Z1[n_2]);
			i++;
		}
		else i++;
	}
	return Z1;
	delete[]A;
}

void print(int N, complex<double>* Z1) //функция для печати полученного вектора, наш вектор нормирован на выходе БПФ
{
	for (int i = 0; i < N; i++) {
		cout << Z1[i] << endl;
	}
	cout << "_______________________________" << endl;
}

complex<double>* my_fft(int N, complex<double>* Z1, complex<double>* Z, complex<double> i1) //функция для выполнения БПФ прореживанием по частоте
{
	int n = static_cast<int>(log2(N));
	for (int k = 1; k < n + 1; k++) {
		for (int j = 0; j < (int)pow(2, (k - 1)); j++) {
			for (int l = 0; l < (int)pow(2, (n - k)); l++) {
				int p = (j * (int)pow(2, (n + 1 - k)) + l);
				int p1 = (int)pow(2, (n - k));
				Z1[p] = Z[p] + Z[p + p1];
				Z1[p + p1] = (Z[p] - Z[p + p1]) * exp((-2.0 * l / (int)pow(2, (n + 1 - k))) * M_PI * i1);
				Z[p] = Z1[p]; Z[p + p1] = Z1[p + p1];
			}
		}
	}
	Z1 = inv(N, n, Z1);
	for (int i = 0; i < N; i++) {
		Z1[i] = 1 / sqrt(N) * Z1[i];
	}
	return Z1;	
}

complex<double>* inv_my_fft(int N, complex<double>* Z2, complex<double>* Z1, complex<double> i1) //функция для вычисления обратного БПФ аналогично прямому БПФ
{
	for (int i = 0; i < N; i++) {
		Z1[i] = conj(Z1[i]);
	}
	Z2 = my_fft(N, Z2, Z1, i1); 
	for (int i = 0; i < N; i++) {
		Z2[i] = conj(Z2[i]);
	}
	return Z2;
}

int main() {
	setlocale(0, "RUS");
	complex<double> i1(0.0, 1.0);
	char R[20] = "";
	cout << "Введите название файла "; cin >> R;
	int N = is_file_correct(R) / 2;
	cout << "Всего элементов в векторе " << N << endl;
	complex<double>* Z = new complex<double>[N];
	double* x = new double[N];
	double* y = new double[N];
	const char* X1 = R;
	complex<double>* Z1 = new complex<double>[N];
	for (int i = 0; i < N; i++) {
		Z1[i] = 0;
	}
	complex<double>* Z2 = new complex<double>[N];
	for (int i = 0; i < N; i++) {
		Z2[i] = 0;
	}
	if (is_pow_of_2(N) == 1) {
		ifstream fin(X1);
		for (int i = 0; i < N; i++) {
			fin >> x[i];
			fin >> y[i];
			Z[i] = complex<double>(x[i], y[i]);
		}
		cout << "Входной вектор" << endl;
		print(N, Z);
		auto begin = chrono::steady_clock::now();
		Z1 = my_fft(N, Z1, Z, i1);//прямое БПФ
		auto end = chrono::steady_clock::now();
		auto elapsed_ms = chrono::duration_cast<chrono::microseconds>(end - begin);	
		cout << "Прямое быстрое преобразование Фурье прореживанием по частоте с перестановками" << endl;
		print(N, Z1);
		cout << "Время работы БПФ при N = " << N << ": " << elapsed_ms.count() << " мкс\n";
		cout << "Обратное быстрое преобразование Фурье прореживанием по частоте с перестановками" << endl;
		Z2 = inv_my_fft(N, Z2, Z1, i1);//обратное БПФ
		print(N, Z2);
	}
	else cout << "Ошибка в заполнении файла" << endl;
	delete[]Z;
	delete[]Z1;
	delete[]Z2;
	delete[]x;
	delete[]y;
	return 0;
}
