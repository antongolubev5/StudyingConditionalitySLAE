// Lab2slau.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <malloc.h>

#include <cmath>
using namespace std;
typedef double mytype;

mytype** readfromfile(string filename, int &n)
{
	char buff[500];
	mytype element;
	int count = 0;
	ifstream fin(filename); // открыли файл для чтения		

	while (!fin.eof())//подсчитываем число строк
	{
		fin.getline(buff, 500);
		count++;
	}
	//выделяем память под массив
	int i;
	mytype **extMatrix;
	extMatrix = new mytype*[count];
	for (i = 0; i<count; i++)
		//a[i], a[i] адресует М элементов типа double
		extMatrix[i] = new mytype[count + 1];

	fin.close(); // закрываем файл
	fin.open(filename); // открыли файл для чтения
	int elcnt = 0;
	while (!fin.eof())
	{
		fin >> element; // считали очередной элемент
		div_t ij = div(elcnt, count + 1);
		extMatrix[ij.quot][ij.rem] = element;
		elcnt++;
	}
	n = count;
	fin.close(); // закрываем файл
	return extMatrix;
}

void matrix_destroyer(mytype** ary, int n)
{
	if (ary != nullptr)
	{
		for (int i = 0; i < n; i++) {
			delete[] ary[i];
		}
		delete[] ary;
		ary = nullptr;
	}
}

void vector_destroyer(mytype* vec, int n)
{
	if (vec != nullptr)
	{
		
		delete[] vec;
		vec = nullptr;
	}
}

mytype** readfromscreen(int &k)
{
	int i, j, raz;
	cout << "Введите размерность матрицы: ";
	cin >> raz;
	k = raz;
	mytype **extMatrix;
	extMatrix = new mytype*[raz];

	for (i = 0; i<raz; i++)
		extMatrix[i] = new mytype[raz + 1];

	for (i = 0; i < raz; i++)
	{
		for (j = 0; j < raz + 1; j++)
		{
			cin >> extMatrix[i][j];
		}
	}
	return extMatrix;
}

mytype** readfromscreensquare(int &k)
{
	int i, j, raz;
	cout << "Введите размерность матрицы A: ";
	cin >> raz;
	k = raz;
	mytype **extMatrix;
	extMatrix = new mytype*[raz];

	for (i = 0; i<raz; i++)
		//a[i], a[i] адресует М элементов типа double
		extMatrix[i] = new mytype[raz];

	for (i = 0; i < raz; i++)
	{
		for (j = 0; j < raz; j++)
		{
			cin >> extMatrix[i][j];
		}
	}
	return extMatrix;
}

void printMatrix(mytype** extMatrix, int k, int m, bool extended)
{
	cout << endl;
	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < m; j++)
		{
			char simb;
			if ((j == k - 1) && extended) { simb = '='; }
			else { simb = ' '; }
			cout << extMatrix[i][j] << simb;

		}
		cout << endl;
	}

}

void printVector(mytype* Vector, int k)
{
	cout << endl;
	for (int i = 0; i < k; i++)
	{
		cout << Vector[i] /*<< " " << i+1 */<< endl;
	}

}

mytype** multiplyMatrix(mytype** Matrix1, mytype** Matrix2, int n)
{
	mytype **result;
	result = new mytype*[n];
	mytype s = 0;
	for (int ii = 0; ii<n; ii++)
		result[ii] = new mytype[n];

	for (int i = 0; i < n; i++)
	{
		for (int l = 0; l < n; l++)
		{
			s = 0;

			for (int j = 0; j < n; j++)
			{
				s += Matrix1[i][j] * Matrix2[j][l];
			}
			result[i][l] = s;
		}
	}
	return result;
}

mytype* multiplyMatrixVector(mytype** Matrix, mytype* Vector, int n)
{
	mytype *result;
	result = new mytype[n];
	mytype s = 0;

	for (int i = 0; i < n; i++)
	{
		s = 0;

		for (int j = 0; j < n; j++)
		{
			s += Matrix[i][j] * Vector[j];
		}
		result[i] = s;

	}
	return result;
}

mytype** multiplyMatrixNumber(mytype** Matrix, mytype Number, int n)
{
	mytype **result;
	result = new mytype*[n];
	for (int ii = 0; ii<n; ii++)
		result[ii] = new mytype[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			result[i][j] = Matrix[i][j] * Number;
		}
	}
	return result;
}

mytype* multiplyVectorNumber(mytype* Vector, mytype Number, int n)
{
	mytype *result;
	result = new mytype[n];

	for (int i = 0; i < n; i++)
	{
		result[i] = Vector[i] * Number;
	}
	return result;
}

mytype* substractVector(mytype* Vector1, mytype* Vector2, int n)
{
	mytype *result;
	result = new mytype[n];

	for (int i = 0; i < n; i++)
	{
		result[i] = Vector1[i] - Vector2[i];
	}

	return result;
}

mytype** substractMatrix(mytype** Matrix1, mytype** Matrix2, int n)
{
	mytype **result;
	result = new mytype*[n];

	for (int ii = 0; ii<n; ii++)
		result[ii] = new mytype[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			result[i][j] = Matrix1[i][j] - Matrix2[i][j];
		}
	}
	return result;
}

mytype** sumMatrix(mytype** Matrix1, mytype** Matrix2, int n)
{
	mytype **result;
	result = new mytype*[n];

	for (int ii = 0; ii<n; ii++)
		result[ii] = new mytype[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			result[i][j] = Matrix1[i][j] + Matrix2[i][j];
		}
	}
	return result;
}

mytype** transposeMatrix(mytype** Matrix, int n)
{
	mytype **result;
	result = new mytype*[n];

	for (int ii = 0; ii<n; ii++)
		result[ii] = new mytype[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			result[i][j] = Matrix[j][i];
			result[j][i] = Matrix[i][j];
		}
	}

	return result;
}

void unitMatrix(mytype** Matrix, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) { Matrix[i][j] = 1; }
			else  Matrix[i][j] = 0;
		}
	}
}

mytype discrepancy(mytype** extMatrix, mytype* solution, int n)
{
	mytype *b1;
	b1 = new mytype[n];
	mytype result = 0;

	b1 = multiplyMatrixVector(extMatrix, solution, n);

	for (int i = 0; i < n - 1; i++)
	{
		result += pow((b1[i] - extMatrix[i][n]), 2);
	}
	result = sqrt(result);
	return result;
}

mytype normVectorUnit(mytype* Vector, int n)
{
	mytype max = 0;
	for (int i = 0; i < n; i++)
	{
		max += abs(Vector[i]);
	}
	return max;
}

mytype normVectorInfinity(mytype* Vector, int n)
{
	mytype max = 0;
	for (int i = 0; i < n; i++)
	{
		if (abs(Vector[i]) > max) max = abs(Vector[i]);
	}
	return max;
}

mytype normMatrixUnit(mytype** Matrix, int n)//макс столбец
{
	mytype max = 0;
	mytype *Vector;
	Vector = new mytype[n];
	for (int i = 0; i < n; i++)
	{
		Vector[i] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Vector[j] += Matrix[i][j];
		}
	}
	max = normVectorInfinity(Vector, n);
	return max;
}

mytype normMatrixInfinity(mytype** Matrix, int n)//макс строка
{
	mytype max = 0;
	mytype *Vector;
	Vector = new mytype[n];
	for (int i = 0; i < n; i++)
	{
		Vector[i] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Vector[i] += Matrix[i][j];
		}
	}
	max = normVectorInfinity(Vector, n);
	return max;
}

mytype* gaussmethod(mytype** extMatrix, int n)
{
	mytype *solution;
	solution = new mytype[n];
	mytype maxvalue = 0;
	int imax;

	for (int cnt = 0; cnt < n; cnt++)
	{
		solution[cnt] = 0;
	}
	//double det = determinant(extMatrix, n);
	//if (abs(det) < 1e-30)
	//{
	//cout << "Определитель равен 0. Не существует единственного решения." << endl;
	//solution = nullptr;
	//}
	//else
	//{

	for (int i = 0; i < n - 1; i++)//цикл по строкам, которые вычитаются из нижележащих
	{
		//выбор макс элемента из i-го столбца
		maxvalue = 0;
		for (int il = i; il < n; il++)
		{
			if (maxvalue < abs(extMatrix[il][i]))
			{
				maxvalue = abs(extMatrix[il][i]);
				imax = il;
			}
		}

		if (maxvalue < 1e-10)
		{
			cout << "Не существует единственного решения." << endl;
			return nullptr;
		}

		if (imax != i)
		{
			mytype* buf = extMatrix[imax];
			extMatrix[imax] = extMatrix[i];
			extMatrix[i] = buf;
		}

		//extMatrix[i][n] = extMatrix[i][n] / extMatrix[i][i];
		mytype aii = extMatrix[i][i];

		if (abs(aii) < 1e-10)
		{
			cout << "Не существует единственного решения. Последняя строка диагонализированной матрицы - нулевая" << endl;
			return nullptr;
		}

		for (int j = i; j <= n; j++)//цикл по элементам строками, которая вычитается из нижележащих  от i+1???
		{
			extMatrix[i][j] = extMatrix[i][j] / aii;
		}

		for (int ii = i + 1; ii < n; ii++)//вычитание из низлежащих строк i-ой строки
		{
			mytype a_ii_i = extMatrix[ii][i];
			for (int jj = i; jj <= n; jj++)
			{
				extMatrix[ii][jj] -= a_ii_i * extMatrix[i][jj];
			}
		}
	}
	//нормируем нижнюю строку
	mytype	 aii = extMatrix[n - 1][n - 1];
	if (abs(aii) < 1e-10)
	{
		cout << "Не существует единственного решения. Последняя строка диагонализированной матрицы - нулевая" << endl;
		return nullptr;
	}
	for (int j = n - 1; j <= n; j++)//цикл по элементам строками, которая вычитается из нижележащих  от i+1???
	{
		extMatrix[n - 1][j] = extMatrix[n - 1][j] / aii;
	}
	//printMatrix(extMatrix, n, n + 1, true);
	//обратный ход

	mytype sum = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		sum = 0;
		for (int j = i + 1; j < n; j++) //суммируем все более старшие переменные  взвешенные на коэффициенты текущей строки
		{
			sum += solution[j] * extMatrix[i][j];
		}
		solution[i] = extMatrix[i][n] - sum;//вычитаем из правой части 
	}

	//printMatrix(extMatrix, n);//печать диагонализированной (для проверки)
	return solution;
}

mytype iterationnumbertheor(mytype eps, mytype* solution, mytype* x0, mytype normC, int n)
{
	mytype k;
	mytype normsubstract = normVectorUnit((substractVector(solution, x0, n)), n);
	k = (log(eps / normsubstract)) / (log(normC));
	return k;
}

mytype estimateepssimpleiterationaprior(mytype* solution, mytype* x0, mytype* xk, mytype normC, int n, int k)
{
	mytype est;
	mytype normsubstract = normVectorUnit((substractVector(x0, solution, n)), n);
	est = pow(normC, k) * normsubstract;
	return est;
}

mytype estimateepssimpleiterationaposter( mytype* solution, mytype* x0, mytype* xk, mytype normC, int n)
{
	mytype est;
	mytype normsubstract = normVectorUnit((substractVector(solution, xk, n)), n);
	est = (normC / (1 - normC)) *normsubstract;
	return est;
}

mytype estimateepsrelaxationaprior(mytype* solution, mytype* x0, mytype* xk, mytype normCu, mytype normCl, int n, int k)
{
	mytype est;
	mytype normsubstract = normVectorUnit((substractVector(x0, xk, n)), n);
	mytype q = normCu/(1-normCl);
	est = pow(q, k) * normsubstract;
	return est;
}

mytype estimateepsrelaxationaposter(mytype* solution, mytype* x0, mytype* xk, mytype normC, mytype normCu, int n)
{
	mytype est;
	mytype normsubstract = normVectorUnit((substractVector(solution, xk, n)), n);
	est = (normCu / (1 - normC)) *normsubstract;
	return est;
}

mytype* simpleiteration(mytype** extMatrix, mytype* x0, mytype eps, int n)

{// для 1-й системы eps=1е-10 tau= 0.000512 -с 1е-6 умножали на 2 в цикле 200 раз (норма 1)
	//mytype *x;//true solution
	/*x = new mytype[n];
	x[0] = 5;
	x[1] = -7;
	x[2] = 12;
	x[3] = 4;*/
	mytype *vectokill;
	mytype *solution;
	solution = new mytype[n];
	mytype *b;
	b = new mytype[n];
	mytype *y;
	y = new mytype[n];
	mytype **A;
	A = new mytype*[n];
	mytype **boofA = nullptr;
	boofA = new mytype*[n];
	mytype **C;
	C = new mytype*[n];
	mytype **E;
	E = new mytype*[n];
	mytype **memtokill;
	mytype iterparam = 	1e-9;
	mytype normCmin = 1e+60;
	mytype normC, iterparammin;

	//diag norm
	mytype diagii = 0;
	for (int str = 0; str < n; str++)
	{
		diagii = extMatrix[str][str];
		for (int col = 0; col < n+1; col++)
		{
			extMatrix[str][col] = extMatrix[str][col] / diagii;
		}
	}
	cout << "Преобразованная матрица:" << endl;
	printMatrix(extMatrix, n, n+1, true);

	for (int ii = 0; ii < n; ii++)
	{
		A[ii] = new mytype[n];
		C[ii] = new mytype[n];
		E[ii] = new mytype[n];
	}

	for (int i = 0; i < n; i++)//заполняем матрицу Е нулями
	{
		for (int j = 0; j < n; j++)
		{
			E[i][j] = 0;
		}
	}

	unitMatrix(E, n);//инициализируем единичную матрицу Е

	for (int i = 0; i < n; i++)//вытаскиваем свободный член из матрицы  в вектор 
	{
		b[i] = extMatrix[i][n];
	}

	for (int i = 0; i < n; i++)//переносим элементы из матрицы  в квадратную матрицу 
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j] = extMatrix[i][j];
		}
	}
	boofA = nullptr;
	//находим такой итерационный параметр, при котором норма матрицы С будет наименьшей
	for (int i = 1; i < 200; i++)
	{
		//вычисляем матрицу С по формуле С=-(tau*A-E)
		//matrix_destroyer(mytype** ary, int n)
		memtokill = boofA;
		boofA = multiplyMatrixNumber(A, iterparam, n);
		matrix_destroyer(memtokill, n);
		memtokill = boofA;
		boofA = substractMatrix(boofA, E, n);
		matrix_destroyer(memtokill, n);
		memtokill = C;
		C = multiplyMatrixNumber(boofA, -1, n);
		matrix_destroyer(memtokill, n);
		/*normC = normMatrixInfinity(C, n);*/
		normC = normMatrixUnit(C, n);
		if (normC < normCmin)
		{
			normCmin = normC;
			iterparammin = iterparam;
		}

		/*printMatrix(C, n, n, false);
		cout << endl;*/
		iterparam *= 2.0;
	}

	//iterparammin = 0.1;
	//находим матрицу С и вектор у, исходя из найденного итерационного параметра
	memtokill = boofA;
	boofA = multiplyMatrixNumber(A, iterparammin, n);

	boofA = substractMatrix(boofA, E, n);
	C = multiplyMatrixNumber(boofA, -1, n);
	matrix_destroyer(memtokill, n);
	y = multiplyVectorNumber(b, iterparammin, n);

	cout << endl;
	cout << "Норма inf матрицы С:" << normCmin << endl;
	cout << "Норма unit матрицы С:" << normMatrixUnit(C, n) << endl;
	eps *= abs((1 - normCmin) / normCmin);
	cout << "tau = " << iterparammin << endl;
	cout << endl;
	cout <<  "Матрица С: " << endl;
	printMatrix(C, n, n, false);
	cout << endl;
	cout << "Вектор у: " << endl;
	printVector(y, n);

	for (int i = 0; i < n; i++)// solution = x0
	{
		solution[i] = x0[i];
	}

	mytype* xk = multiplyVectorNumber(solution, 10, n);
	xk[0] += 10;
	mytype* minusy = multiplyVectorNumber(y, -1, n);// y = -y
	int cnt = 0;
	while (normVectorUnit(substractVector(solution, xk, n), n) > eps&&cnt<1e7)
	{
		for (int w = 0; w < n; w++)
		{
			xk[w] = solution[w];
		}

		vectokill= solution;
		solution = multiplyMatrixVector(C, solution, n);
		vector_destroyer(vectokill, n);
		vectokill = solution;
		solution = substractVector(solution, minusy, n);
		vector_destroyer(vectokill, n);
		/*cout << "iteration = " << cnt << endl;
		cout <<"xk+1 = " << endl;
		printVector(solution, n);
		cout << "xk = " << endl;
		printVector(xk, n);
		cout << "xk+1-xk " << endl;
		cout << normVectorUnit(substractVector(solution, xk, n), n) << endl;
		cout << "xk+1-xtrue " << endl;
		//cout << normVectorUnit(substractVector(solution, x, n), n) << endl;*/
		cnt++;
	}
	cout << endl;
	cout << "Количество итераций: " << cnt << endl;
	cout << endl;
	mytype kol = round(iterationnumbertheor(eps, solution, x0, normCmin, n));
	cout << "Теоретическая оценка количества итераций: " << kol << endl;
	cout << endl;
	mytype est = estimateepssimpleiterationaprior(solution, x0, xk, normCmin, n, cnt);
	cout << "Априорная погрешность решения = " << est << endl;
	cout << endl;

	est = estimateepssimpleiterationaposter(solution, x0, xk, normCmin, n);
	cout << "Апостериорная погрешность решения = " << est << endl;
	return solution;
}

mytype* jacobimethod(mytype** extMatrix, mytype* x0, mytype eps, int n)
{
	mytype normC;
	mytype *vectokill;
	mytype *solution;
	solution = new mytype[n];
	mytype *b;
	b = new mytype[n];
	mytype *y;
	y = new mytype[n];

	mytype **A;
	A = new mytype*[n];
	mytype **boofA = nullptr;
	boofA = new mytype*[n];
	mytype **C;
	C = new mytype*[n];
	mytype **E;
	E = new mytype*[n];
	mytype **L;
	L = new mytype*[n];
	mytype **D;
	D = new mytype*[n];
	mytype **U;
	U = new mytype*[n];

	for (int ii = 0; ii < n; ii++)
	{
		A[ii] = new mytype[n];
		C[ii] = new mytype[n];
		E[ii] = new mytype[n];
		L[ii] = new mytype[n];
		D[ii] = new mytype[n];
		U[ii] = new mytype[n];
	}
	//подготавливаем матрицу Е, достаем из расширенной матрицы матрицу А и столбец b соответственно
	for (int i = 0; i < n; i++)//заполняем новые матрицы нулями
	{
		for (int j = 0; j < n; j++)
		{
			E[i][j] = 0;
			L[i][j] = 0;
			D[i][j] = 0;
			U[i][j] = 0;
		}
	}

	unitMatrix(E, n);//инициализируем единичную матрицу Е

	for (int i = 0; i < n; i++)//вытаскиваем свободный член из матрицы  в вектор 
	{
		b[i] = extMatrix[i][n];
	}

	for (int i = 0; i < n; i++)//переносим элементы из матрицы  в квадратную матрицу 
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j] = extMatrix[i][j];
		}
	}
	boofA = nullptr;
	//представляем матрицу А в виде трех матриц: верхнетреугольная L, диагональная D  и нижнетреугольная U

	for (int i = 0; i < n; i++)//нижнетреугольная L
	{
		for (int j = 0; j < n; j++)
		{
			if (i > j) { L[i][j] = A[i][j]; }
		}
	}

	/*cout << "Матрица L: " << endl;
	printMatrix(L, n, n, false);
	cout << endl;*/

	for (int i = 0; i < n; i++)//диагональная D
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) { D[i][j] = A[i][j]; }
		}
	}

	/*cout << "Матрица D: " << endl;
	printMatrix(D, n, n, false);
	cout << endl;*/

	for (int i = 0; i < n; i++)//верхнетреугольная U
	{
		for (int j = 0; j < n; j++)
		{
			if (i < j) { U[i][j] = A[i][j]; }
		}
	}

	/*cout << "Матрица U: " << endl;
	printMatrix(U, n, n, false);
	cout << endl;*/

	//заполняем матрицу С элементами
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i != j) { C[i][j] = -(A[i][j]/A[i][i]) ; }
			else (C[i][j] = 0);
		}
	}

	cout << "Матрица С: "<< endl;
	printMatrix(C, n, n, false);
	cout << endl;
	normC = normMatrixUnit(C, n);
	cout << "Норма unit матрицы С:" << normC << endl;
	cout << "Норма inf матрицы С:" << normMatrixUnit(C, n) << endl;
	cout << endl;
	eps *= abs((1 - normC) / normC);

	//заполняем столбец y элементами
	for (int i = 0; i < n; i++)
	{
		y[i] = b[i] / A[i][i];
	}

	cout << "Вектор у: " << endl;
	printVector(y, n);

	for (int i = 0; i < n; i++)// solution = x0
	{
		solution[i] = x0[i];
	}

	mytype* xk = multiplyVectorNumber(solution, 10, n);
	xk[0] += 10;
	mytype* minusy = multiplyVectorNumber(y, -1, n);// y = -y
	int cnt = 0;

	while (normVectorUnit(substractVector(solution, xk, n), n) > eps&&cnt<1e7)
	{
		for (int w = 0; w < n; w++)
		{
			xk[w] = solution[w];
		}

		vectokill = solution;
		solution = multiplyMatrixVector(C, solution, n);
		vector_destroyer(vectokill, n);
		vectokill = solution;
		solution = substractVector(solution, minusy, n);
		vector_destroyer(vectokill, n);
		/*cout << "iteration = " << cnt << endl;
		cout <<"xk+1 = " << endl;
		printVector(solution, n);
		cout << "xk = " << endl;
		printVector(xk, n);
		cout << "xk+1-xk " << endl;
		cout << normVectorUnit(substractVector(solution, xk, n), n) << endl;
		cout << "xk+1-xtrue " << endl;
		//cout << normVectorUnit(substractVector(solution, x, n), n) << endl;*/
		cnt++;
	}
	cout << endl;
	cout << "Количество итераций: " << cnt << endl;
	cout << endl;
	mytype kol = round(iterationnumbertheor(eps, solution, x0, normC, n));
	cout << "Теоретическая оценка количества итераций: " << kol << endl;
	cout << endl;

	mytype est = estimateepssimpleiterationaprior(solution, x0, xk, normC, n, cnt);
	cout << "Априорная погрешность решения = " << est << endl;
	cout << endl;

	est = estimateepssimpleiterationaposter(solution, x0, xk, normC, n);
	cout << "Апостериорная погрешность решения = " << est << endl;

	return solution;
}

mytype* relaxationmethoddiag(mytype* vecta, mytype* vectb, mytype* vectc, mytype* vectd, mytype* x0, mytype eps, int n)
{
	mytype *solution;
	solution = new mytype[n];
	mytype relaxparam = 1;
    
	for (int i = 0; i < n; i++)// solution = x0
	{
		solution[i] = x0[i];
	}

	mytype* xk = multiplyVectorNumber(solution, 10, n);
	xk[0] += 10;//xk!=x0
	int cnt = 0;

	while (normVectorUnit(substractVector(solution, xk, n), n) > eps&&cnt<1e7)
	{
		for (int w = 0; w < n; w++)
		{
			xk[w] = solution[w];
		}

		solution[0] = (1 - relaxparam)*xk[0] - relaxparam*vectc[0] / vectb[0] * xk[1] + relaxparam*vectd[0] / vectb[0];

		for (int i = 1; i < n-1; i++)
		{
			solution[i] = -relaxparam*vecta[i] / vectb[i] * solution[i - 1] + (1 - relaxparam)*xk[i]
				- relaxparam*vectc[i] / vectb[i] * xk[i + 1] + relaxparam*vectd[i] / vectb[i];
		}

		solution[n-1] = -relaxparam*vecta[n - 2] / vectb[n-1] * solution[n - 2] + (1 - relaxparam)*xk[n-1]
			 + relaxparam*vectd[n-1] / vectb[n-1];

		cnt++;
		//printVector(solution, n);
	}

    cout << endl;
	cout << "Количество итераций = " << cnt << endl;
	
	return solution;
}

mytype** reverseMatrix(mytype** Matrix, int n)
{
	mytype **extMatrix;
	extMatrix = new mytype*[n];
	mytype **result;
	result = new mytype*[n];
	mytype *ordinary;
	ordinary = new mytype[n];

	for (int ii = 0; ii<n; ii++)
		result[ii] = new mytype[n];
	for (int ii = 0; ii<n; ii++)
		extMatrix[ii] = new mytype[n + 1];

	for (int j = 0; j < n; j++)//заполняем правый столбец extMatrix нулями
	{
		extMatrix[j][n] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		for (int ii = 0; ii < n; ii++)//переносим элементы из матрицы в расширенную матрицу
		{
			for (int jj = 0; jj < n; jj++)
			{
				extMatrix[ii][jj] = Matrix[ii][jj];
			}
		}
		for (int jj = 0; jj < n; jj++)//заполняем правый столбец extMatrix нулями
		{
			extMatrix[jj][n] = 0;
		}

		extMatrix[i][n] = 1;

		ordinary = gaussmethod(extMatrix, n);
		if (ordinary == nullptr)
		{
			cout << "Не существует обратной матрицы." << endl;
			return nullptr;
		}
		for (int j = 0; j < n; j++)
		{
			result[j][i] = ordinary[j];
		}

	}

	return result;
}

mytype* relaxationmethodgeneral(mytype** extMatrix, mytype* x0, mytype eps, int n)
{
	mytype relaxparam = 1;
	mytype normCu, normC, normCl;
	mytype *solution;
	solution = new mytype[n];
	mytype **A;
	A = new mytype*[n];
	mytype **E;
	E = new mytype*[n];
	mytype *b;
	b = new mytype[n];
	mytype **L;
	L = new mytype*[n];
	mytype **D;
	D = new mytype*[n];
	mytype **U;
	U = new mytype*[n];
	mytype **C;
	C = new mytype*[n]; 
	mytype **Cu;
	Cu = new mytype*[n];
	mytype **Cl;
	Cl = new mytype*[n];
	mytype **boofMatrix1;
	boofMatrix1 = new mytype*[n];
	mytype **boofMatrix2;
	boofMatrix2 = new mytype*[n];
	mytype **boofMatrix3;
	boofMatrix3 = new mytype*[n];

	for (int ii = 0; ii < n; ii++)
	{
		A[ii] = new mytype[n];
		E[ii] = new mytype[n];
		C[ii] = new mytype[n];
		Cu[ii] = new mytype[n];
		Cl[ii] = new mytype[n];
		L[ii] = new mytype[n];
		D[ii] = new mytype[n];
		U[ii] = new mytype[n];
		boofMatrix1[ii] = new mytype[n];
		boofMatrix2[ii] = new mytype[n];
		boofMatrix3[ii] = new mytype[n];
	}

	unitMatrix(E, n);//инициализируем единичную матрицу Е

	//разделим расширенную матрицу на квадратную и столбец
	for (int i = 0; i < n; i++)//переносим элементы из матрицы  в квадратную матрицу 
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j] = extMatrix[i][j];
		}
	}

	for (int i = 0; i < n; i++)//вытаскиваем свободный член из матрицы  в вектор 
	{
		b[i] = extMatrix[i][n];
	}

	for (int i = 0; i < n; i++)//заполняем новые матрицы нулями
	{
		for (int j = 0; j < n; j++)
		{
			L[i][j] = 0;
			D[i][j] = 0;
			U[i][j] = 0;
			Cu[i][j] = 0;
			Cl[i][j] = 0;
			boofMatrix1[i][j] = 0;
			boofMatrix2[i][j] = 0;
			boofMatrix3[i][j] = 0;
		}
	}

	//представляем матрицу А в виде трех матриц: нижнетреугольная L, диагональная D  и верхнетреугольная U

	for (int i = 0; i < n; i++)//нижнетреугольная L
	{
		for (int j = 0; j < n; j++)
		{
			if (i > j) { L[i][j] = A[i][j]; }
		}
	}

	cout << "Матрица L: " << endl;
	printMatrix(L, n, n, false);
	cout << endl;

	for (int i = 0; i < n; i++)//диагональная D
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) { D[i][j] = A[i][j]; }
		}
	}

	cout << "Матрица D: " << endl;
	printMatrix(D, n, n, false);
	cout << endl;

	for (int i = 0; i < n; i++)//верхнетреугольная U
	{
		for (int j = 0; j < n; j++)
		{
			if (i < j) { U[i][j] = A[i][j]; }
		}
	}

	cout << "Матрица U: " << endl;
	printMatrix(U, n, n, false);
	cout << endl;
	//находим матрицу С из матричного уравнения
	boofMatrix1 = reverseMatrix(D, n);
	boofMatrix1 = multiplyMatrixNumber(boofMatrix1, relaxparam, n);
	boofMatrix1 = multiplyMatrix(boofMatrix1, L, n);
	boofMatrix1 = sumMatrix(E, boofMatrix1, n);
	boofMatrix2 = multiplyMatrixNumber(boofMatrix2, 1-relaxparam, n);
	boofMatrix2 = multiplyMatrix(boofMatrix2, E, n);
	boofMatrix3 = reverseMatrix(D, n);
	boofMatrix3 = multiplyMatrixNumber(boofMatrix3, relaxparam, n);
	boofMatrix3 = multiplyMatrix(boofMatrix3, U, n);
	boofMatrix2 = substractMatrix(boofMatrix2, boofMatrix3, n);
	boofMatrix1 = reverseMatrix(boofMatrix1, n);
	boofMatrix3 = multiplyMatrix(boofMatrix3, boofMatrix1, n);

	normC = normMatrixUnit(boofMatrix3, n);
	//выделяем верхнетреугольную матрицу Сu
	for (int i = 0; i < n; i++)//верхнетреугольная U
	{
		for (int j = 0; j < n; j++)
		{
			if (i < j) { Cu[i][j] = boofMatrix3[i][j]; }
		}
	}

	normCu = normMatrixUnit(Cu, n);

	//выделяем нижнетреугольную матрицу Сl
	for (int i = 0; i < n; i++)//нижнетреугольная L
	{
		for (int j = 0; j < n; j++)
		{
			if (i > j) { Cl[i][j] = boofMatrix3[i][j]; }
		}
	}

	normCl = normMatrixUnit(Cu, n);

	for (int i = 0; i < n; i++)// solution = x0
	{
		solution[i] = x0[i];
	}

	mytype* xk = multiplyVectorNumber(solution, 10, n);
	xk[0] += 10;//xk!=x0
	int cnt = 0;
	mytype sum1, sum2;
	eps *= abs((1 - normC) / normCu);

	while (normVectorUnit(substractVector(solution, xk, n), n) > eps&&cnt<1e7)
	{
		for (int w = 0; w < n; w++)
		{
			xk[w] = solution[w];
		}

		//solution[0] = (1 - relaxparam)*xk[0] - relaxparam*vectc[0] / vectb[0] * xk[1] + relaxparam*vectd[0] / vectb[0];

		for (int i = 0; i < n; i++)
		{
			sum1 = sum2 = 0;
			for (int j = 0; j < n; j++)
			{
				if (j < i) { sum1 += A[i][j] / A[i][i] * solution[j]; }
				if (j > i) { sum2 += A[i][j] / A[i][i] * xk[j]; }
			}
			solution[i] = -relaxparam*sum1 + (1 - relaxparam)*xk[i] - relaxparam*sum2 + relaxparam*b[i] / A[i][i];
		}

		//solution[n - 1] = -relaxparam*vecta[n - 2] / vectb[n - 1] * solution[n - 2] + (1 - relaxparam)*xk[n - 1]
			//+ relaxparam*vectd[n - 1] / vectb[n - 1];

		cnt++;
		//printVector(solution, n);
	}

	cout << "Количество итераций = " << cnt << endl;
	cout << endl;
	mytype kol = round(iterationnumbertheor(eps, solution, x0, normC, n));
	cout << "Теоретическая оценка количества итераций: " << kol << endl;
	cout << endl;

	mytype est = estimateepsrelaxationaprior(solution, x0, xk, normCu, normCl, n, cnt);
	cout << "Априорная погрешность решения = " << est << endl;
	cout << endl;

	est = estimateepsrelaxationaposter(solution, x0, xk, normC, normCu, n);
	cout << "Апостериорная погрешность решения = " << est << endl;
	return solution;
}

int main()
{
	setlocale(LC_ALL, "rus");
	int k;
	mytype eps = 1e-10;
	
	mytype** extMatrix = readfromfile("test1.txt", k);
	mytype *x0;
	x0 = new mytype[k];
	for (int i = 0; i < k; i++)
	{
		x0[i] = 1;
	}

	/*x0[0] = 5;
	x0[1] = -7;
	x0[2] = 12;
	x0[3] = 4;*/
    
	x0[0] = 1;
	x0[1] = 1;
	x0[2] = 1;
	x0[3] = 1;
	
	cout << "Вектор начального приближения: "<< endl;
	printVector(x0, k);;
	cout << endl;

	cout << "Решение СЛАУ методом простой итерации: " << endl;
	printMatrix(extMatrix, k, k + 1, true);
	cout << endl;
	mytype *solution = simpleiteration(extMatrix, x0, eps, k);

	cout << endl;
	cout << "Решение: " << endl;
	printVector(solution, k);
	cout << endl;

	mytype disсrep = discrepancy(extMatrix, solution, k);
	cout << "Норма вектора невязки = "<< disсrep << endl;
	cout << endl;

	/*cout << "Решение СЛАУ методом Якоби: " << endl;
	printMatrix(extMatrix, k, k + 1, true);
	cout << endl;
	mytype *solution = jacobimethod(extMatrix, x0, eps, k);

	cout << endl;
	cout << "Решение: " << endl;
	printVector(solution, k);
	cout << endl;

	mytype disсrep = discrepancy(extMatrix, solution, k);
	cout << "Норма вектора невязки = " << disсrep << endl;
	cout << endl;*/

	//определим 4 вектора для трехдиагональной матрицы (для методов Зейделя и релаксации)
	k = 205;//размерность по условию

	mytype *vecta;//поддиагональ
	vecta = new mytype[k-1];
	mytype *vectb;//диагональ
	vectb = new mytype[k];
	mytype *vectc;//наддиагональ
	vectc = new mytype[k-1];
	mytype *vectd;//правая часть
	vectd = new mytype[k];

	mytype  *x0large;//вектор начального приближение
	x0large = new mytype[k];
	for (int i = 0; i < k; i++)
	{
		x0large[i] = 1;
		//x0large[i] = 2 - ((i + 1) % 2);
	}
	
	for (int i = 0; i < k-1; i++)
	{
		vecta[i] = vectc[i] = 1;
	}


	for (int i = 0; i < k; i++)
	{
		vectb[i] = 4;

		if (i == 0)
		{
			vectd[i] = 6;
		}
		else  if (i == k - 1)
		{
			vectd[i] = 9 - 3 * ((k % 2));
		}
		else 
		{
			vectd[i] = 10 - 2 * ((i+1) % 2);
		}
	}

	/* printVector(vecta, k-1);
	 cout << endl;
	 printVector(vectb, k);
	 cout << endl;
	 printVector(vectc, k - 1);
	 cout << endl;
	 printVector(vectd, k);
	 cout << endl;*/

	/*cout << endl;
	cout << "Решение трехдиагональной СЛАУ большой размерности методом релаксации: " << endl;
	mytype *solution = relaxationmethoddiag(vecta, vectb, vectc, vectd, x0large, eps, k);

	cout << endl;
	cout << "Решение: " << endl;
	printVector(solution, k);
	cout << endl;*/

	/*cout << "Решение СЛАУ методом релаксации: " << endl;
	extMatrix = readfromfile("test1.txt", k);
	printMatrix(extMatrix, k, k + 1, true);
	cout << endl;
	mytype *solution = relaxationmethodgeneral(extMatrix, x0, eps, k);

	cout << "Решение: " << endl;
	printVector(solution, k);
	cout << endl;*/
	system("pause");
	return 0;
}




