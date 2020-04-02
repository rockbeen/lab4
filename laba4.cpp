#include<iostream>
#include<string>
#include<cstdlib>
#include<cmath>
using namespace std;

double NORM(double **x, double n) {
double normiA = 0;
for (int j = 0; j < n; j++) {
	normiA += abs(x[0][j]);
}
for (int i = 1; i < n; i++) {
	double normA = 0;
	for (int j = 0; j < n; j++) {
		normA +=abs(x[i][j]);

	}
	if (normA > normiA) {
		normiA = normA;
	}
}
cout << "NORM:" << normiA << endl << endl;
return normiA;
}
bool shodimost(double *xk, double *xkp)
{
	double norm = 0;
	for (int i = 0; i < 20; i++)
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

double * gauss(double **a, double *y, int n)
{
	double *x, max;
	int k, index;
	const double eps = 0.00001;  // точность
	x = new double[n];
	k = 0;
	while (k < n)
	{
		
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			cout << "Reshenie nevosmoghno polychit` iz-za null stolbca";
			cout << index << " matrix A" << endl;
			return 0;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < n; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}
	return x;
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

	cout << "Matrix A:" << endl;
	for (int i = 0; i < n; i++) {
		x[i] = new double[n];
		for (int j = 0; j < n; j++) {
			x[i][j] = 0;
		}
	}


	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
				x[i][j] = 0;
			if (i >0 and j ==0)
				x[i][j] = a[j];
			if (i >0 and j ==1)
				x[i][j] = b[j];
			if (i >0 and j ==2)
				x[i][j] = c[j];
			if (i==9)
				{x[i][j]=q[j];}
				x[0][0]=b[0];
				x[0][1]=c[0];
				x[0][2]=0;
	}}



	

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << x[i][j] << " ";
		}
		cout <<f[i]<< endl;
	}
	cout << endl;
	//(задание 3)
	double norm1;
	norm1 = NORM(x, n);

	//Гаус
	cout << "Reshenie po Gausu:" << endl;
	t = gauss(x, f, n);
	for (int i = 0; i < n; i++) {
		cout << "x[" << i << "]=" << t[i] << endl;
	}
	cout << endl;

	//Вектор невязки для метода Гауса
	double podstav, *vector;
	vector = new double[n];
	cout << "Vector neviazki Gaus V=(";
	for (int i = 0; i < n; i++) {
		podstav = 0;
		for (int j = 0; j < n; j++) {
			if(x[i][j]!=0)
			podstav += t[j] * x[i][j];
		}
		vector[i] = abs(f[i]-podstav);
		cout << vector[i] << ",";
	}
	cout << ")" << endl<<endl;

	//Гаус-Зейдель
	cout << "Reshenie po Gausu-Zeidely:" << endl;
	double *p, *r;
	p = new double[n];
	r = new double[n];
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
	} 
	while (!shodimost(r, p));
	for (int i = 0; i < n; i++) {
		cout << "x[" << i << "]=" << r[i] << endl;
	}
	cout << endl;

	//Вектор невязки для метода Гауса-Зейделя
	double podstav1, *vector1;
	vector1 = new double[n];
	cout << "Vector neviazki Gaus-Zeidel V=(";
	for (int i = 0; i < n; i++) {
		podstav1 = 0;
		for (int j = 0; j < n; j++) {
			if (x[i][j] != 0)
				podstav1 += r[j] * x[i][j];
		}
		vector1[i] = abs(f[i] - podstav1);
		cout << vector1[i] << ",";
	}
	cout << ")" << endl<<endl;

	
	//Обратная матрица, её норма и вычисление числа обусловленности (делается последним)
	cout << "Inversion matrix A:" << endl;
	inversion(x, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			
			if (fabs(x[i][j]) < 0.00001)
				x[i][j] = 0;
			cout << x[i][j] << ", ";
		}
		cout << endl;
	}
	double norm2 = NORM(x, n);
	cout << "Chislo obyslovlennosti u=" << norm1 * norm2 <<endl;


	return 0;
}
