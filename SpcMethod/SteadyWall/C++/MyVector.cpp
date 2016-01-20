#include "MyVector.h"
//In-class functions
template <typename T>
vector<T>::vector(int s)
{
	size = s;
	_elem = new T[s];
}

template <typename T>
vector<T>::~vector()
{
	delete [] _elem;
}

template <typename T>
T& vector<T>::operator[](int rank)
{
	return _elem[rank];
}

template <typename T>
vector<T>& vector<T>::operator=(const vector<T> &v2)
{
	for (int i = 0; i < size; i++)
	{
		_elem[i] = v2[i];
	}
	return *this;
}

template <typename T>
vector<T>& vector<T>::operator+=(const vector<T> &v2)
{
	for (int i = 0; i < size; i++)
	{
		_elem[i] += v2[i];
	}
	return *this;
}

template <typename T>
vector<T>& vector<T>::operator-=(const vector<T> &v2)
{
	for (int i = 0; i < size; i++)
	{
		_elem[i] -= v2[i];
	}
	return *this;
}

template <typename T>
vector<T>& vector<T>::operator*=(const vector<T> &v2)
{
	for (int i = 0; i < size; i++)
	{
		_elem[i] *= v2[i];
	}
	return *this;
}

template <typename T>
vector<T>& vector<T>::operator/=(const vector<T> &v2)
{
	for (int i = 0; i < size; i++)
	{
		_elem[i] /= v2[i];
	}
	return *this;
}

template <typename T>
vector<T>& vector<T>::operator=(double scaler)
{
	for (int i = 0; i < size; i++)
	{
		_elem[i] = scaler;
	}
	return *this;
}

template <typename T>
vector<T>& vector<T>::operator+=(double scaler)
{
	for (int i = 0; i < size; i++)
	{
		_elem[i] += scaler;
	}
	return *this;
}

template <typename T>
vector<T>& vector<T>::operator-=(double scaler)
{
	for (int i = 0; i < size; i++)
	{
		_elem[i] -= scaler;
	}
	return *this;
}

template <typename T>
vector<T>& vector<T>::operator*=(double scaler)
{
	for (int i = 0; i < size; i++)
	{
		_elem[i] *= scaler;
	}
	return *this;
}

template <typename T>
vector<T>& vector<T>::operator/=(double scaler)
{
	for (int i = 0; i < size; i++)
	{
		_elem[i] /= scaler;
	}
	return *this;
}

//Out-class functions
template <typename T>
vector<T> operator+(const vector<T> &v1, const vector<T> &v2)
{
	vector<T> temp(v1.size);
	for (int i = 0; i < v1.size; i++)
	{
		temp[i] = v1[i] + v2[i];
	}
	return temp;
}

template <typename T>
vector<T> operator-(const vector<T> &v1, const vector<T> &v2)
{
	vector<T> temp(v1.size);
	for (int i = 0; i < v1.size; i++)
	{
		temp[i] = v1[i] - v2[i];
	}
	return temp;
}

template <typename T>
vector<T> operator*(const vector<T> &v1, const vector<T> &v2)
{
	vector<T> temp(v1.size);
	for (int i = 0; i < v1.size; i++)
	{
		temp[i] = v1[i] * v2[i];
	}
	return temp;
}

template <typename T>
vector<T> operator/(const vector<T> &v1, const vector<T> &v2)
{
	vector<T> temp(v1.size);
	for (int i = 0; i < v1.size; i++)
	{
		temp[i] = v1[i] / v2[i];
	}
	return temp;
}

template <typename T>
vector<T> operator+(const vector<T> &v, double scaler)
{
	vector<T> temp(v.size);
	for (int i = 0; i < v.size; i++)
	{
		temp[i] = v[i] + scaler;
	}
	return temp;
}

template <typename T>
vector<T> operator-(const vector<T> &v, double scaler)
{
	vector<T> temp(v.size);
	for (int i = 0; i < v.size; i++)
	{
		temp[i] = v[i] - scaler;
	}
	return temp;
}

template <typename T>
vector<T> operator*(const vector<T> &v, double scaler)
{
	vector<T> temp(v.size);
	for (int i = 0; i < v.size; i++)
	{
		temp[i] = v[i] * scaler;
	}
	return temp;
}

template <typename T>
vector<T> operator/(const vector<T> &v, double scaler)
{
	vector<T> temp(v.size);
	for (int i = 0; i < v.size; i++)
	{
		temp[i] = v[i] / scaler;
	}
	return temp;
}

template <typename T>
vector<T> operator+(double scaler, const vector<T> &v)
{
	vector<T> temp(v.size);
	for (int i = 0; i < v.size; i++)
	{
		temp[i] = scaler + v[i];
	}
	return temp;
}

template <typename T>
vector<T> operator-(double scaler, const vector<T> &v)
{
	vector<T> temp(v.size);
	for (int i = 0; i < v.size; i++)
	{
		temp[i] = scaler - v[i];
	}
	return temp;
}

template <typename T>
vector<T> operator*(double scaler, const vector<T> &v)
{
	vector<T> temp(v.size);
	for (int i = 0; i < v.size; i++)
	{
		temp[i] = scaler * v[i];
	}
	return temp;
}

template <typename T>
vector<T> operator/(double scaler, const vector<T> &v)
{
	vector<T> temp(v.size);
	for (int i = 0; i < v.size; i++)
	{
		temp[i] = scaler / v[i];
	}
	return temp;
}
