﻿#include<iostream>
#include<string>
#include <fstream>
#include<cstdlib>
#include<cmath>
using namespace std;

std::ofstream fout;

double NORM(double **x, double n) {
	double normiA = 0;
	for (int j = 0; j < n; j++) {
		normiA += abs(x[0][j]);
	}
	for (int i = 1; i < n; i++) {
		double normA = 0;
		for (int j = 0; j < n; j++) {
			normA += abs(x[i][j]);

		}
		if (normA > normiA) {
			normiA = normA;
		}
	}
	cout << "NORMI:" << normiA << endl << endl;
	return normiA;
}


bool shodimost(double *xk, double *xkp)
{
	double norm = 0;
	for (int i = 0; i < 10; i++)
		norm += (xk[i] - xkp[i])*(xk[i] - xkp[i]);
	return (sqrt(norm) < 0.0000001);
}

void inversion(double **A, int N)
{
	double temp;

	double **E = new double *[N];

	for (int i = 0; i < N; i++)
		E[i] = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			E[i][j] = 0.0;

			if (i == j)
				E[i][j] = 1.0;
		}

	for (int k = 0; k < N; k++)
	{
		temp = A[k][k];

		for (int j = 0; j < N; j++)
		{
			A[k][j] /= temp;
			E[k][j] /= temp;
		}

		for (int i = k + 1; i < N; i++)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = N - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = E[i][j];

	for (int i = 0; i < N; i++)
		delete[] E[i];

	delete[] E;
}

double * gauss_method() {
double **x, *f, *t;
	int n = 10;
	x = new double*[n];
	f = new double[n];
	double a[9], b[9], d[9], c[9],q[10];
	 b[0]=1,c[0]=0,f[0]=1,q[0]=1,q[9]=1;

for (int i = 0; i < n; i++) {
			a[i] = 1;
		}
for (int i = 1; i < n; i++) {
			c[i] = 1;
			b[i] = -2;
			q[i] = 2;
		}
for (int i = 1; i < n; i++) {
			f[i] = 2/pow(i,2);
			
		}
f[n-1] = -n/3;
	for (int i = 0; i < n; i++) {
		x[i] = new double[n];
		for (int j = 0; j < n; j++) {
			x[i][j] = 0;
		}
	}


		for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j)
				x[i][j] = b[j];
			if (i > 0) {
				if (j == i - 1)
					x[i][j] = a[i];}
			if (i < n - 1) {
				if (j == i + 1)
					x[i][j] = c[i];
			}
			if (i==9)
				{x[i][j]=q[j];}
				x[0][0]=b[0];
				x[0][1]=c[0];
				x[0][2]=0;
	}}
	for (int i = 0; i < n; i++) {
		x[i][n] = f[i]; 
	}
	double *r;
	r = new double[n];
	double m = 0;
	for (int k = 1; k < n; k = k + 1) {
		for (int j = k; j < n; j = j + 1) {
			m = x[j][k - 1] / x[k - 1][k - 1];
			for (int i = 0; i < n + 1; i = i + 1) {
				x[j][i] = x[j][i] - m * x[k - 1][i];
			}
		}
	}
	for (int i = n - 1; i >= 0; i = i - 1) {
		r[i] = x[i][n] / x[i][i];
		for (int c = n - 1; c > i; c = c - 1) {
			r[i] = r[i] - x[i][c] * r[c] / x[i][i];
		}
	}
	return r;
}

int main(int argc, char *argv[])
{

		
double **x, *f, *t;
	int n = 10;
	x = new double*[n];
	f = new double[n];
	double a[9], b[9], d[9], c[9],q[10];
	 b[0]=1,c[0]=0,f[0]=1,q[0]=1,q[9]=1;

for (int i = 0; i < n; i++) {
			a[i] = 1;
		}
for (int i = 1; i < n; i++) {
			c[i] = 1;
			b[i] = -2;
			q[i] = 2;
		}
for (int i = 1; i < n; i++) {
			f[i] = 2/pow(i,2);
			
		}
f[n-1] = -n/3;
for (int i = 0; i < n-1; i++) {
	if (i==0)
		cout<< b[i]<<"x"<<i+1<<"+"<<c[i]<<"x"<<i+2<<"="<<f[i]<<endl;
	else{
		cout<< a[i]<<"x"<<i<< b[i]<<"x"<<i+1<<"+"<<c[i]<<"x"<<i+2<<"="<<f[i]<<endl;
		}
	}
for (int i = 0; i < n; i++) {
	if (i==9)
	cout<< q[i]<<"x"<<i+1<<"="<<f[i]<<endl;
	else{
			cout<< q[i]<<"x"<<i+1<<"+";
		}
	}

	//Задание матрицы
	cout << "Matrix A:" << endl;
	for (int i = 0; i < n; i++) {
		x[i] = new double[n];
		for (int j = 0; j < n; j++) {
			x[i][j] = 0;
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j)
				x[i][j] = b[j];
			if (i > 0) {
				if (j == i - 1)
					x[i][j] = a[i];}
			if (i < n - 1) {
				if (j == i + 1)
					x[i][j] = c[i]; 
			}
			if (i==9)
				{x[i][j]=q[j];}
				x[0][0]=b[0];
				x[0][1]=c[0];
				x[0][2]=0;
	}}


	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << x[i][j] << ", ";
		}
		cout << endl;
	}
	cout << endl;
	//Высчитываем норму для начальной матрицы (задание 3)
	double norm1;
	norm1 = NORM(x, n);

	//Гаус
	cout << "Reshenie po Gausu:" << endl;
	std::ofstream fout;
	t = gauss_method();
	fout.open("ans1.dat");
	for (int i = 0; i < n; i++) {
		cout << "x[" << i << "]=" << t[i] << endl;
		fout << t[i] << endl;
	}
	cout << endl;

	//Вектор невязки для метода Гауса
	double podstav, vector;

	double normVec = 0;
	cout << "Norm vector neviazki Gaus V=";
	for (int i = 0; i < 1; i++) {
		podstav = 0;
		for (int j = 0; j < n; j++) {
			if (x[i][j] != 0)
				podstav += t[j] * x[i][j];
		}
		vector = abs(f[i] - podstav);
		normVec = vector;
		
	}
	for (int i = 1; i < n; i++) {
		podstav = 0;
		for (int j = 0; j < n; j++) {
			if (x[i][j] != 0)
				podstav += t[j] * x[i][j];
		}
		vector = abs(f[i] - podstav);
		if (vector >= normVec) {
			normVec = vector;
		}
		
	}
	
	cout << normVec << endl << endl;
	//////////////////////////////////////////////////////////////////////

	//Гаус-Зейдель
	cout << "Reshenie po Gausu-Zeidely:" << endl;
	std::ofstream fout1;
	fout1.open("ans2.dat");
	double *p, *r;
	p = new double[n];
	r = new double[n];
	int iteration = 0;
	do
	{
		for (int i = 0; i < n; i++)
			p[i] = r[i];
		for (int i = 0; i < n; i++)
		{
			double var = 0;
			for (int j = 0; j < i; j++)
				var += (x[i][j] * r[j]);
			for (int j = i + 1; j < n; j++)
				var += (x[i][j] * p[j]);
			r[i] = (f[i] - var) / x[i][i];
		}
		iteration++;
	} while (!shodimost(r, p));
	for (int i = 0; i < n; i++) {
		cout << "x[" << i << "]=" << r[i] << endl;
		fout1 <<  r[i] << endl;
	}
	cout << "Iteration:" << iteration << endl;
	cout << endl;

	//Вектор невязки для метода Гауса-Зейделя
	double podstav1, vector1;
	double normVec1 = 0;
	cout << "Vector neviazki Gaus-Zeidel V=";
	for (int i = 0; i < 1; i++) {
		podstav1 = 0;
		for (int j = 0; j < n; j++) {
			if (x[i][j] != 0)
				podstav1 += r[j] * x[i][j];
		}
		vector1 = abs(f[i] - podstav1);
		normVec1 = vector1;
		
	}
	for (int i = 1; i < n; i++) {
		podstav1 = 0;
		for (int j = 0; j < n; j++) {
			if (x[i][j] != 0)
				podstav1 += r[j] * x[i][j];
		}
		vector1 = abs(f[i] - podstav1);
		if (vector1 >= normVec1) {
			normVec1 = vector1;
		}
		
	}
	
	cout << normVec1 << endl << endl;




	//Обратная матрица, её норма и вычисление числа обусловленности (делается последним)
	cout << "Inversion matrix A:" << endl;
	inversion(x, n);

	cout << "Normi inversion matrix: ";
	double norm2 = NORM(x, n);
	cout << "Chislo obyslovlennosti u=" << norm1 * norm2 << endl;



	return 0;
}
