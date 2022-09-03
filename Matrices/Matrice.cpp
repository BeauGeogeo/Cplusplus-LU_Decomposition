#include "Matrice.h"

#include <iostream>
#include <cassert>

Matrice::Matrice(int rows, int columns)
{
	m_rows = rows;
	m_columns = columns;
	m_values = new float* [m_rows]; 
	for (int i = 0; i < m_rows; i++)
	{	
		m_values[i] = new float[m_columns];
	}
}

Matrice::Matrice(const Matrice& matrice)
{
	m_rows = matrice.m_rows;
	m_columns = matrice.m_columns;
	m_values = new float* [m_rows];

	for (int i = 0; i < m_rows; i++)
	{	
		m_values[i] = new float[m_columns];
		for (int j = 0; j < m_columns; j++)
		{	
			setValue(i, j, matrice.getValue(i, j));
		}
	}
}

Matrice::~Matrice()
{
	for (int i = 0; i < m_rows; i++)
	{
		delete[] m_values[i];
	}
	delete[] m_values;

}

Matrice Matrice::subMatrix(int firstRow, int lastRow, int firstCol, int lastCol) const
{	
	assert(0 <= firstRow <= lastRow < m_rows);
	assert(0 <= firstCol <= lastCol < m_columns);
	Matrice newMatrix(lastRow - firstRow + 1, lastCol - firstCol + 1); 
	for (int rowIndex = 0; rowIndex < newMatrix.m_rows; rowIndex++)
	{
		for (int colIndex = 0; colIndex < newMatrix.m_columns; colIndex++)
		{
			newMatrix.setValue(rowIndex, colIndex, getValue(firstRow + rowIndex, firstCol + colIndex));
		}
	}
	return newMatrix;
} 

Matrice& Matrice::operator=(const Matrice& matrice) 
{
	for (int i = 0; i < m_rows; i++)
	{
		delete[] m_values[i];
	}
	delete[] m_values;
	m_rows = matrice.m_rows;
	m_columns = matrice.m_columns;
	m_values = new float* [m_rows]; 
	for (int i = 0; i < m_rows; i++)
	{
		m_values[i] = new float[m_columns];
		for (int j = 0; j < m_columns; j++)
		{
			setValue(i, j, matrice.getValue(i, j));
		}
	}
	return *this;
}

Matrice& Matrice::operator+(const Matrice& matrice)
{	
	assert(m_rows == matrice.m_rows);
	assert(m_columns == matrice.m_columns);
	for (int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_columns; j++)
		{
			setValue(i, j, getValue(i, j) + matrice.getValue(i, j));
		}
	}
	return *this;
} 

Matrice Matrice::operator+=(const Matrice& matrice) const
{
	assert(m_rows == matrice.m_rows);
	assert(m_columns == matrice.m_columns);
	Matrice res(m_rows, m_columns);
	for (int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_columns; j++)
		{
			res.setValue(i, j, getValue(i, j) + matrice.getValue(i, j));
		}
	}
	return res;
}

Matrice& Matrice::operator*(const Matrice& matrice)
{
	assert(m_columns == matrice.m_rows);
	for (int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < matrice.m_columns; j++)
		{
			float coeff_ij(0);
			for (int k = 0; k < m_columns; k++)
			{
				coeff_ij += getValue(i, k) * matrice.getValue(k, j);
			}
			setValue(i, j, coeff_ij);
		}
	}
	return *this;
}

Matrice Matrice::operator*=(const Matrice& matrice) const
{
	assert(m_columns == matrice.m_rows);
	Matrice res(m_rows, matrice.m_columns);
	for (int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < matrice.m_columns; j++)
		{
			float coeff_ij(0);
			for (int k = 0; k < m_columns; k++)
			{
				coeff_ij += getValue(i, k) * matrice.getValue(k, j);
			}
			res.setValue(i, j, coeff_ij);
		}
	}
	return res;
}

Matrice& Matrice::operator-(const Matrice& matrice)
{
	assert(m_rows == matrice.m_rows);
	assert(m_columns == matrice.m_columns);
	for (int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_columns; j++)
		{
			setValue(i, j, getValue(i, j) - matrice.getValue(i, j));
		}
	}
	return *this;
}

Matrice Matrice::operator-=(const Matrice& matrice) const
{
	assert(m_rows == matrice.m_rows);
	assert(m_columns == matrice.m_columns);
	Matrice res(m_rows, m_columns);
	for (int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_columns; j++)
		{
			res.setValue(i, j, getValue(i, j) - matrice.getValue(i, j));
		}
	}
	return res;
}

Matrice Matrice::power(int n) const
{
		if (!n)
		{
			Matrice res(m_columns, m_columns);
			for (int i = 0; i < m_columns; i++)
			{
				for (int j = 0; j < m_columns; j++)
				{
					i == j ? res.setValue(i, j, 1) : res.setValue(i, j, 0);
				}
			}
			return res;
		}
		Matrice res(m_columns, m_columns);
		for (int k = 0; k < n; k++)
		{
			res = res * *this;
		}
		return res;
}


void Matrice::setValue(int i, int j, float x)
{	
	*(m_values[i] + j) = x;
}

float Matrice::getValue(int i, int j) const
{
	return *(m_values[i] + j);
}

Matrice& Matrice::addEye(float value, int diag)
{	
	assert(m_rows == m_columns);
	assert(-m_rows < diag < m_rows);
	if (diag > 0)
	{
		for (int i = 0; i < m_columns - diag; i++)
		{
			setValue(i, i + diag, value);
		}
	}
	else if (diag < 0)
	{
		for (int i = 0; i < m_columns + diag; i++)
		{
			setValue(i - diag, i, value);
		}
	}
	else
	{
		for (int i = 0; i < m_columns; i++)
		{
			setValue(i, i, value);
		}
	}
	return *this;
}

Matrice& Matrice::makeZeros()
{
	for (int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_columns; j++)
		{
			setValue(i, j, 0);
		}
	}
	return *this;
}

void Matrice::switchRows(int row1, int row2)
{	
	assert(0 <= row1 < m_rows); 
	assert(0 <= row2 < m_rows);
	assert(row1 != row2);
	float* tempPtr = m_values[row1]; 
	m_values[row1] = m_values[row2];
	m_values[row2] = tempPtr;
}

void Matrice::printValues() const
{
	std::cout << std::endl;
	for (int i = 0; i < m_rows; i++)
	{	
		std::cout << "[" << getValue(i, 0);
		for (int j = 1; j < m_columns; j++)
		{
			std::cout << " " << getValue(i, j);
		}
		std::cout << "]" << std::endl;
	}
	std::cout << std::endl;
}

auto Matrice::decompoLU() const
{
	assert(m_rows == m_columns);
	assert(m_rows >= 2);
	
	struct finalValues 
	{
		Matrice R, P, L; 
		int nbPermutations;
	};

	Matrice R(*this);
	Matrice P(m_rows, m_columns);
	P = makeId(m_rows);
	Matrice L(m_rows, m_columns);
	L.makeZeros();
	int nbPermutations(0);
	
	for (int k = 0; k < m_rows - 1; k++) 
	{	
		int indexPivot(k);
		float pivot(R.getValue(indexPivot, indexPivot));

		if (!pivot)
		{	
			indexPivot++;
			while ((indexPivot < m_rows) && (!pivot)) 
			{	
				if (R.getValue(indexPivot, k))
				{
					pivot = R.getValue(indexPivot, k);
					P.switchRows(indexPivot, k);
					L.switchRows(indexPivot, k);
					R.switchRows(indexPivot, k);
					nbPermutations++;
				}
				else 
				{
					indexPivot++; 
				}
			}
		}
		
		for (int i = indexPivot + 1; i < m_rows; i++)
		{	
			float mulCoeff(-R.getValue(i, k) / pivot);
			L.setValue(i, k, -mulCoeff);

			for (int j = k; j < m_columns; j++)
			{	
				float Rij(R.getValue(i, j));
				float newValue(R.getValue(k, j) * mulCoeff + Rij);
				R.setValue(i, j, newValue);
			}
		}
	}
	L.addEye(1, 0);
	return finalValues {R, P, L, nbPermutations};
}

float Matrice::determinant() const
{	
	assert(m_rows == m_columns);
	if (!m_rows)
	{
		return 1;
	}
	else if (m_rows == 1)
	{
		return getValue(0, 0);
	}
	else
	{	
		float res(1);
		auto values = decompoLU();
		for (int i = 0; i < m_rows; i++)
		{
			res *= values.R.getValue(i, i);
		}
		return res * pow(-1, values.nbPermutations);
	}
}

Matrice Matrice::makeId(int n) const
{
	Matrice res(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			i == j ? res.setValue(i, j, 1) : res.setValue(i, j, 0);
		}
	}
	return res;
}

void Matrice::getSize(int* size) const
{
	size[0] = m_rows;
	size[1] = m_columns;
}

int main()
{	
	int totalSizes[4];
	Matrice mat(3, 3);
	mat.getSize(totalSizes); 
	for (int i = 0; i < 3; i++) 
	{
		for (int j = 0; j < 3; j++)
		{
			mat.setValue(i, j, i);
		}
	}
	Matrice matDet(mat);
	mat.printValues();
	Matrice mat2(5, 3);
	int size[2];
	mat2.getSize(size);
	mat2.getSize(totalSizes + 2);
	for (int i = 0; i < 5; i++) 
	{
		for (int j = 0; j < 3; j++)
		{
			mat2.setValue(i, j, i);
		}
	}
	Matrice mat3(mat2);
	std::cout << "mat3 avant changement" << std::endl;
	mat3.printValues();
	mat3.setValue(4, 2, 0);
	std::cout << "mat3 APRES changement, puis mat2" << std::endl;
	mat3.printValues();
	mat2.printValues();
	mat = mat3 = mat2;
	mat.printValues();
	Matrice mat4(5,3);
	mat4 = mat3 - mat2; 
	mat4.printValues();
	Matrice mat4_bis(5, 3);
	mat4_bis = mat4 - mat2 - mat2;
	mat4_bis.printValues();
	mat4_bis -= mat2;
	mat4_bis.printValues();
	mat * matDet * matDet;
	mat.printValues();
	int nb(0);
	mat *= nb;
	mat.printValues();
	mat.setValue(0, 0, 1);
	mat.printValues();
	mat.switchRows(0, 1);
	mat.printValues();
	std::cout << nb << std::endl;
	Matrice matDet2(matDet);
	matDet2.makeZeros();
	matDet2.setValue(0, 1, 1);
	matDet2.setValue(0, 2, 1);
	matDet2.setValue(1, 0, 1);
	matDet2.setValue(1, 1, 2);
	matDet2.setValue(1, 2, 1);
	matDet2.setValue(2, 0, 2);
	matDet2.setValue(2, 1, 7);
	matDet2.setValue(2, 2, 9);
	
	const auto val = matDet2.decompoLU(); 

	std::cout << std::endl;
	std::cout << "DECOMPOSITION PA = LU (AVEC P LA MATRICE DES PERMUTATIONS)" << std::endl;
	std::cout << std::endl;
	std::cout << "LA MATRICE DE DEPART" << std::endl;

	matDet2.printValues();

	std::cout << "LA MATRICE DE U" << std::endl;

	val.R.printValues();

	std::cout << "LA MATRICE L" << std::endl;

	val.L.printValues();

	std::cout << "LA MATRICE DE PERMUTATIONS (NECESSAIRE POUR CHANGER LE PIVOT QUAND IL EST NUL)" << std::endl;

	val.P.printValues();

	std::cout << "LE NOMBRE DE PERMUTATIONS FAITES PENDANT LA DECOMPOSITION PA = LU " << std::endl;
	std::cout << std::endl;

	std::cout << val.nbPermutations << std::endl;
	std::cout << std::endl;

	std::cout << "LE DETERMINANT DE LA MATRICE INITIALE QU ON CALCULE GRACE A LA FACTORISATION PA = LU" << std::endl;
	std::cout << std::endl;

	std::cout << matDet2.determinant() << std::endl;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			matDet2.setValue(i, j, 1);
		}
	}

	std::cout << std::endl;
	std::cout << "DECOMPOSITION PA = LU (AVEC P LA MATRICE DES PERMUTATIONS)" << std::endl;
	std::cout << std::endl;
	std::cout << "LA MATRICE DE DEPART" << std::endl;

	matDet2.printValues();


	const auto val2 = matDet2.decompoLU();

	std::cout << "LA MATRICE DE U" << std::endl;

	val2.R.printValues();

	std::cout << "LA MATRICE L" << std::endl;

	val2.L.printValues();

	std::cout << "LA MATRICE DE PERMUTATIONS (NECESSAIRE POUR CHANGER LE PIVOT QUAND IL EST NUL)" << std::endl;

	val2.P.printValues();

	std::cout << "LE NOMBRE DE PERMUTATIONS FAITES PENDANT LA DECOMPOSITION PA = LU " << std::endl;
	std::cout << std::endl;

	std::cout << val2.nbPermutations << std::endl;
	std::cout << std::endl;

	std::cout << "LE DETERMINANT DE LA MATRICE INITIALE QU ON CALCULE GRACE A LA FACTORISATION PA = LU" << std::endl;
	std::cout << std::endl;

	std::cout << matDet2.determinant() << std::endl;
	
	matDet2.makeZeros();
	const auto val3 = matDet2.decompoLU();

	std::cout << std::endl;
	std::cout << "DECOMPOSITION PA = LU (AVEC P LA MATRICE DES PERMUTATIONS)" << std::endl;
	std::cout << std::endl;
	std::cout << "LA MATRICE DE DEPART" << std::endl;

	matDet2.printValues();

	std::cout << "LA MATRICE DE U" << std::endl;

	val3.R.printValues();

	std::cout << "LA MATRICE L" << std::endl;

	val3.L.printValues();

	std::cout << "LA MATRICE DE PERMUTATIONS (NECESSAIRE POUR CHANGER LE PIVOT QUAND IL EST NUL)" << std::endl;

	val3.P.printValues();

	std::cout << "LE NOMBRE DE PERMUTATIONS FAITES PENDANT LA DECOMPOSITION PA = LU " << std::endl;
	std::cout << std::endl;

	std::cout << val3.nbPermutations << std::endl;
	std::cout << std::endl;

	std::cout << "LE DETERMINANT DE LA MATRICE INITIALE QU ON CALCULE GRACE A LA FACTORISATION PA = LU" << std::endl;
	std::cout << std::endl;

	std::cout << matDet2.determinant() << std::endl;
	return 0;
}
