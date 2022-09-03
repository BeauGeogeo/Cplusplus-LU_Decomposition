#ifndef MATRICE_H
#define MATRICE_H
#include <cmath>
#include <cassert>
#include <iostream>

class Matrice {

private:

	int m_rows;
	int m_columns; 
	float** m_values;

public:

	Matrice(int, int);
	Matrice(const Matrice&);
	~Matrice();
	Matrice& operator=(const Matrice&);
	Matrice& operator+(const Matrice&);
	Matrice operator+=(const Matrice&) const;
	Matrice& operator*(const Matrice&);
	Matrice operator*=(const Matrice&) const; 
	Matrice& operator-(const Matrice&);
	Matrice operator-=(const Matrice&) const;
	void switchRows(int, int);
	Matrice& addEye(float, int);
	Matrice& makeZeros();
	auto decompoLU() const; 
	float determinant() const;
	void setValue(int, int, float);
	float getValue(int, int) const;
	Matrice power(int) const;
	Matrice makeId(int) const;
	Matrice subMatrix(int, int, int, int) const;
	void printValues() const;
	void getSize(int*) const;
};

template<typename T>
T power(const T& var, int n)
{
	T res(1);
	if (!n)
	{
		return res;
	}
	for (int i = 0; i < n; i++)
	{
		res *= var;
	}
	return res;
}

template <typename T>
T polyCarTest(const T& var)
{
	//return ((var * power2(var)) * -1) + power2(var) + ((var * 3) - (power0(var) * 3));
	return -power(var, 3) + power(var, 2) + 3 * var - 3 * power(var, 0);
}


/*
template <typename T>
Matrice& operator*(const T& scalar, Matrice& matrice)
{	
	int size[2];
	matrice.getSize(size);
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			matrice.setValue(i, j, matrice.getValue(i, j) * scalar);
		}
	}
	return *this;
}

template <typename T>
Matrice& operator*(Matrice& matrice, const T& scalar)
{	
	int size[2];
	matrice.getSize(size);
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			matrice.setValue(i, j, matrice.getValue(i, j) * scalar);
		}
	}
	return *this;
}
Pb pour celles-ci, car du coup j'ai mat * matDet * matDet qui deconne et me mets les 3 surcharges pour * sont similaires
https://stackoverflow.com/questions/1726740/c-error-operator-2-overloads-have-similar-conversions elts de rep ici
AH ok pb arrive car le type T même si s'appelle scalar bah ça peut être une matrice et dans ce cas fait même chose que autres opérateurs *...
*/

// on met juste un argument T sinon le template s'attend à + et pb
// si je fais des *= de deux matrices je risque d'avoir le mm pb. Comment corriger, en enlevant le scalar ou plutot en fixant un type et conversion 
// implicite ? ou en appelant l'operateur avec les <> ?


// DEMANDER qd je renvoie res là ça marche pas. Je dois le stocker qque part. Ok pas de sens en fait de faire *= et de renvoyer qqch faut bosser sur la ref
// et eventuellement renvoyer la matrice elle-même. mais si je renvoie et stocke pas, useless
template <typename T>
Matrice& operator*=(const T& scalar, Matrice& matrice)
{
	int size[2];
	matrice.getSize(size);
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			matrice.setValue(i, j, matrice.getValue(i, j) * scalar);
		}
	}
	return matrice;
}

template <typename T>
Matrice& operator*=(Matrice& matrice, const T& scalar)
{	
	int size[2];
	matrice.getSize(size);
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			matrice.setValue(i, j, matrice.getValue(i, j) * scalar);
		}
	}
	return matrice;
}

#endif