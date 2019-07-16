
#include <stdio.h>
#include <cstring>
#include <cmath>
#include <iostream>

int n;
using namespace std;

//*************************************************************************/
//Загрузка данных
//*************************************************************************/
/*int Load_Data()
{
	printf("Input filename with data\n");
	char FileName[20];
	int i = 0;
	FILE* File;
	scanf_s("%s", &FileName, 20);
	if (!fopen_s(&File, FileName, "r"))
	{
		printf("file is opened\n");
	}
	else
	{
		printf("%s: file doesn't exist\n", FileName);
		return -1;
	}
	double x, f, f2;
	fscanf_s(File, "%d", &n);
	KnotArray = new knot[n + 2];
	while (!feof(File))
	{
		fscanf_s(File, "%lf%lf%lf", &x, &f, &f2);
		KnotArray[i].Add(x, f, f2);
		i++;
		if (i == n)
			return 1;
	}
	fclose(File);
	return 1;
}*/

/*
int Load_Data_to_mass()//считывает в массив
{
	printf("Input filename with data\n");
	int i = 0;
	FILE* File;
	if (!fopen_s(&File, "D:\data_SIN.txt", "r"))
	{
		printf("file is opened\n");
	}
	else
	{
		printf("data_SIN.txt doesn't exist\n");
		return -1;
	}
	//double x, f, f2;
	fscanf_s(File, "%d", &n);
	
	double* x = new double[n];
	double* y = new double[n];
	double* z = new double[n];

	while (!feof(File))
	{
		fscanf_s(File, "%lf%lf%lf", &x[i], &y[i], &z[i]);
		i++;
		if (i == n)
		{
			fclose(File);
			return 1;
		}
	}
	fclose(File);
	return 1;
}*/
void Spline(double* x, double* y, // input
	double* A, double* B,     // output
	double* C, double* D)     // output
{
	int N = n - 1;
	double* w = new double[N];
	double* h = new double[N];
	double* ftt = new double[N+1];

	for (int i = 0; i < N; i++)
	{
		w[i] = (x[i + 1] - x[i]);
		h[i] = (y[i + 1] - y[i]) / w[i];
	}

	ftt[0] = 0;
	for (int i = 0; i < N - 1; i++)
		ftt[i + 1] = 3 * (h[i + 1] - h[i]) / (w[i + 1] + w[i]);
	ftt[N] = 0;

	for (int i = 0; i < N; i++)
	{
		A[i] = (ftt[i + 1] - ftt[i]) / (6 * w[i]);
		B[i] = ftt[i] / 2;
		C[i] = h[i] - w[i] * (ftt[i + 1] + 2 * ftt[i]) / 6;
		D[i] = y[i];
	}
	delete[] w, h, ftt;

}
void PrintSpline(double* x,            // input
	double* A, double* B, // input
	double* C, double* D) // input
{
	int N = n - 1;
	for (int i = 0; i < N; i++)
	{
		cout << x[i] << " <= x <= " << x[i + 1] << " : f(x) = ";
		cout << A[i] << "(x-" << x[i] << ")^3 + ";
		cout << B[i] << "(x-" << x[i] << ")^2 + ";
		cout << C[i] << "(x-" << x[i] << ")^1 + ";
		cout << D[i] << "(x-" << x[i] << ")^0";
		cout << endl;
	}
}
double Interpolate(double x_interp, double* x, double* A, double* B, double* C, double* D)
{
	//double result;
	int i = 0;

	while ((x[i] < x_interp) && (i < n - 1))
		i++;
	if (i > 0) i--;
	return D[i] + C[i] * (x_interp - x[i]) + B[i] * powf((x_interp - x[i]), 2) +A[i] * powf((x_interp - x[i]), 3);

}

int main()
{
	FILE* File_Out;
	FILE* File_In;

	int N;//количество значений функции-интерполянты
	printf("Input the number of interpolant values:\n");
	scanf_s("%d", &N);

	//------открываем файл с данными------
	if (!fopen_s(&File_In, "D:\data_SIN.txt", "r"))
		printf("file is opened\n");
	else
		printf("data_SIN.txt doesn't exist\n");

	//------считываем количество точек-----
	fscanf_s(File_In, "%d", &n);

	//------массивы под исходные данные-----
	double* x = new double[n];
	double* y = new double[n];
	double* z = new double[n];
	//------массивы под коэффициенты сплайнов-----
	double* A = new double[n-1];
	double* B = new double[n-1];
	double* C = new double[n-1];
	double* D = new double[n-1];
	//------массивы под значения функции-интеполянты-----
	double* init = new double[N];
	double* fun = new double[N];
	double* fun_z = new double[N];
	//------считываем данные-----
	int i = 0;
	while (!feof(File_In))
	{
		fscanf_s(File_In, "%lf%lf%lf", &x[i], &y[i], &z[i]);
		i++;
		if (i == n)
			break;
	}

	for (int i = 0; i < n; i++)
	{
		printf("x:  %lf\ty:   %lf\tz:   %lf\n", x[i], y[i], z[i]);
	}
	//------строим сплайн для y=f(x)----
	Spline(x, y, A, B, C, D);

	for (int i = 0; i < n-1; i++)
	{
		printf("A:  %lf\tB:   %lf\tC:   %lf\tD:   %lf\n", A[i], B[i], C[i], D[i]);
	}

	PrintSpline(x, A, B, C, D);

	//------интерполируем y=f(x)----
	for (int i = 0; i < N; i++)
	{
		init[i] = x[0] + (x[n - 1] - x[0]) / (N - 1) * i;
		fun[i] = Interpolate(init[i], x, A, B, C, D);
		//printf("i=%d init = %lf\tfun = %lf\n", i, init[i], fun[i]);
	}
	//------строим сплайн для z=f(x)----
	Spline(x, z, A, B, C, D);
	
	//------интерполируем y=f(x)----
	for (int i = 0; i < N; i++)
	{
		fun_z[i] = Interpolate(init[i], x, A, B, C, D);
		//printf("i=%d init = %lf\tfun = %lf\n", i, init[i], fun_z[i]);
	}
	//-----записываем данные в файл----
	if (fopen_s(&File_Out, "D:\interp.dat", "wb"))
		printf("File could not be opened\n");
	else
	{
		for (int i = 0; i < N; i++)
		{
			printf("x:  %lf\ty_interp:   %lf\tz_interp:   %lf\n", init[i], fun[i], fun_z[i]);
			fwrite(&init[i], sizeof(double), 1, File_Out);
			fwrite(&fun[i], sizeof(double), 1, File_Out);
			fwrite(&fun_z[i], sizeof(double), 1, File_Out);
		}
	}

/*	if (Load_Data() != -1)
	{
		for (int i = 0; i <= N; i++)
		{
			init[i] = KnotArray[0].x + (KnotArray[n - 1].x - KnotArray[0].x) / N * i;

		}
		BuildSpline();
		fun[0] = KnotArray[0].f;
		fun_z[0] = KnotArray[0].f2;
		fun[N] = KnotArray[n - 1].f;
		fun_z[N] = KnotArray[n - 1].f2;

		double x_add_r = KnotArray[n - 1].x + 0.5236;
		double x_add_l = KnotArray[0].x - 0.5236;
		double f_add_r = Interpolate(x_add_r);
		double f_add_l = Interpolate(x_add_l);
		printf("x_add=%lf\tf_add=%lf\n", x_add_r, f_add_r);
		printf("x_add=%lf\tf_add=%lf\n", x_add_l, f_add_l);


		for (int i = 1; i < N; i++)
		{
			fun[i] = Interpolate(init[i]);
			//		printf("i=%d f=%lf\n", i, fun[i]);
		}
		BuildSpline_z();
		for (int i = 1; i < N; i++)
		{
			fun_z[i] = Interpolate(init[i]);
		}
		if (fopen_s(&file_out, "D:\interp.dat", "wb"))
			printf("File could not be opened\n");
		else
		{
			for (int i = 0; i <= N; i++)
			{
				printf("x:  %lf\ty_interp:   %lf\tz_interp:   %lf\n", init[i], fun[i], fun_z[i]);
				fwrite(&init[i], sizeof(double), 1, file_out);
				fwrite(&fun[i], sizeof(double), 1, file_out);
				fwrite(&fun_z[i], sizeof(double), 1, file_out);
			}
		}
	}
	delete[] init;
	delete[] fun;
	delete[] fun_z;
	for (int count = 0; count < n - 1; count++)
		delete[] Coef[count];*/
	
	delete[] x,y,z;
	delete[] A, B, C, D;
	delete[] init, fun, fun_z;
	fclose(File_Out);
	fclose(File_In);
	system("pause");
	return 0;
}
