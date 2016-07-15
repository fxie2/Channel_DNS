#pragma once
//#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <memory>

template<typename T>
class Field3d;

template<typename T>
std::ostream& operator<<(std::ostream& out, const Field3d<T>& field);

template<typename T>
class Field3d
{
public:
	// Constructors and Dstructor
	Field3d();
	Field3d(const Field3d<T>& field);
	Field3d(Field3d<T>&& field);
	Field3d(size_t up1, size_t up2, size_t up3, size_t dn1 = 1, size_t dn2 = 1, size_t dn3 = 1);
	~Field3d();
	// Memory Allocation
	void reAlloc(size_t up1, size_t up2, size_t up3, size_t dn1 = 1, size_t dn2 = 1, size_t dn3 = 1);
	void deAlloc();
	// Copy, ValueSet and Math operators
	Field3d<T>& operator=(const Field3d<T>& field2);
	Field3d<T>& operator=(Field3d<T>&& temp_field2);
	Field3d<T>& operator=(const T& c);
	T& operator()(size_t i, size_t j, size_t k);
	T& operator()(size_t i, size_t j, size_t k) const;
	Field3d<T>& operator+=(const Field3d<T>& field2);
	Field3d<T>& operator*=(const Field3d<T>& field2);
	Field3d<T>& operator-=(const Field3d<T>& field2);
	Field3d<T>& operator/=(const Field3d<T>& field2);
	Field3d<T>& operator+=(const T& c);
	Field3d<T>& operator*=(const T& c);
	Field3d<T>& operator-=(const T& c);
	Field3d<T>& operator/=(const T& c);
	int dim1() { return _up1 - _dn1 + 1; };
	int dim2() { return _up2 - _dn2 + 1; };
	int dim3() { return _up3 - _dn3 + 1; };

	//friend Field3d<T>&& operator+(const Field3d<T>& f1, const Field3d<T>& f2);
	//friend Field3d<T>&& operator-(const Field3d<T>& f1, const Field3d<T>& f2);
	//friend Field3d<T>&& operator*(const Field3d<T>& f1, const Field3d<T>& f2);
	//friend Field3d<T>&& operator/(const Field3d<T>& f1, const Field3d<T>& f2);
	//friend Field3d<T>&& operator+(const Field3d<T>& f, const T& c);
	//friend Field3d<T>&& operator-(const Field3d<T>& f, const T& c);
	//friend Field3d<T>&& operator*(const Field3d<T>& f, const T& c);
	//friend Field3d<T>&& operator/(const Field3d<T>& f, const T& c);
	// Initialize functions
	void setRand();
	void setZero();
	T maxval()
	{
		int mi, mj, mk;
		mi = 1;
		mj = 1;
		mk = 1;
		T mx = _data[0];
		for (size_t i = _dn1; i <= _up1; i++)
		{
			for (size_t j = _dn2; j <= _up2; j++)
			{
				for (size_t k = _dn3; k <= _up3; k++)
				{
					if (this->operator()(i, j, k) > mx)
					{
						mx = this->operator()(i, j, k);
						mi = i; mj = j; mk = k;
					}
				}
			}
		}
		std::cout << mi << '\t' << mj << '\t' << mk << '\t' << mx << endl;
		return mx;
	}
	// IO functions
	//friend std::ostream& operator<<(std::ostream& out, const Field3d<T>& field);
	friend Field3d<T> operator+(const Field3d<T>& f1, const Field3d<T>& f2)
	{
		Field3d<T> temp(f1);
		temp += f2;
		return temp;
	}

	friend Field3d<T> operator-(const Field3d<T>& f1, const Field3d<T>& f2)
	{
		Field3d<T> temp(f1);
		temp -= f2;
		return temp;
	}

	friend Field3d<T> operator*(const Field3d<T>& f1, const Field3d<T>& f2)
	{
		Field3d<T> temp(f1);
		temp *= f2;
		return temp;
	}

	friend Field3d<T> operator/(const Field3d<T>& f1, const Field3d<T>& f2)
	{
		Field3d<T> temp(f1);
		temp /= f2;
		return temp;
	}

	friend Field3d<T> operator+(const Field3d<T>& f, const T& c)
	{
		Field3d<T> temp(f);
		temp += c;
		return temp;
	}

	friend Field3d<T> operator-(const Field3d<T>& f, const T& c)
	{
		Field3d<T> temp(f);
		temp -= c;
		return temp;
	}

	friend Field3d<T> operator*(const Field3d<T>& f, const T& c)
	{
		Field3d<T> temp(f);
		temp *= c;
		return temp;
	}

	friend Field3d<T> operator/(const Field3d<T>& f, const T& c)
	{
		Field3d<T> temp(f);
		temp /= c;
		return temp;
	}

	friend std::ostream & operator<<(std::ostream & out, const Field3d<T>& field)
	{
		out << std::scientific << std::setprecision(8) << std::endl;
		out << "\n";
		out << "dimx = " << field._dn1 << ":" << field._up1 << "\t";
		out << "dimy = " << field._dn2 << ":" << field._up2 << "\t";
		out << "dimz = " << field._dn1 << ":" << field._up3 << "\n";
		for (size_t k = field._dn3; k <= field._up3; k++)
		{
			out << "\nz = " << k << "\n\n";
			out << "From x = " << field._dn1 << " to x = " << field._up1 << " , y = " << field._dn2 << " to y = " << field._up2 << " :\n";
			for (size_t i = field._dn1; i <= field._up1; i++)
			{
				for (size_t j = field._dn1; j <= field._up2; j++)
				{
					out << field(i, j, k) << " ";
				}
				out << "\n";
			}
		}
		return out;
	}
private:
	size_t _up1, _up2, _up3;
	size_t _dn1, _dn2, _dn3;
	T* _data;
	bool _alloced;
};

template<typename T>
Field3d<T>::Field3d() :_up1(0), _up2(0), _up3(0), _dn1(0), _dn2(0), _dn3(0), _data(nullptr), _alloced(false) {}

template<typename T>
inline Field3d<T>::Field3d(const Field3d<T>& field)
	:_up1(field._up1), _up2(field._up2), _up3(field._up3), _dn1(field._dn1), _dn2(field._dn2), _dn3(field._dn3), _alloced(field._alloced)
{
	if (_alloced)
	{
		_data = new T[(_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn2 + 1)];
		memcpy(_data, field._data, sizeof(T) * (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1));
	}
	else
	{
		_data = nullptr;
	}
}

template<typename T>
inline Field3d<T>::Field3d(Field3d<T>&& field)
{
	_data = field._data;
	_alloced = field._alloced;
	_up1 = field._up1;
	_up2 = field._up2;
	_up3 = field._up3;
	_dn1 = field._dn1;
	_dn2 = field._dn2;
	_dn3 = field._dn3;
	field._data = nullptr;
	field._alloced = false;
}

template<typename T>
inline Field3d<T>::Field3d(size_t up1, size_t up2, size_t up3, size_t dn1, size_t dn2, size_t dn3)
{
	_up1 = up1;
	_up2 = up2;
	_up3 = up3;
	_dn1 = dn1;
	_dn2 = dn2;
	_dn3 = dn3;
	_data = new T[(_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1)];
	_alloced = true;
}

template<typename T>
Field3d<T>::~Field3d()
{
	if (_alloced) delete[] _data;
}

template<typename T>
inline void Field3d<T>::reAlloc(size_t up1, size_t up2, size_t up3, size_t dn1, size_t dn2, size_t dn3)
{
	if (!_alloced)
	{
		_data = new T[(up1 - dn1 + 1) * (up2 - dn2 + 1) * (up3 - dn3 + 1)];
		_alloced = true;
	}
	else if((up1 - dn1) != (_up1 - _dn1) || (up2 - dn2) != (_up2 - _dn2) || (up3 - dn3) != (_up3 - _dn3))
	{
		delete[] _data;
		_data = new T[(up1 - dn1 + 1) * (up2 - dn2 + 1) * (up3 - dn3 + 1)];
	}
	_up1 = up1;
	_dn1 = dn1;
	_up2 = up2;
	_dn2 = dn2;
	_up3 = up3;
	_dn3 = dn3;
}

template<typename T>
inline void Field3d<T>::deAlloc()
{
	if (_alloced) delete[] _data;
	_alloced = false;
}

template<typename T>
inline T& Field3d<T>::operator()(size_t i, size_t j, size_t k)
{
	return _data[(i - _dn1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1) + (j - _dn2) * (_up3 - _dn3 + 1) + k - _dn3];
}

template<typename T>
inline T & Field3d<T>::operator()(size_t i, size_t j, size_t k) const
{
	return _data[(i - _dn1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1) + (j - _dn2) * (_up3 - _dn3 + 1) + k - _dn3];
}

template<typename T>
inline Field3d<T>& Field3d<T>::operator=(const Field3d<T>& field2)
{
	for (size_t i = 0; i < (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1); i++)
	{
		_data[i] = field2._data[i];
	}
	return *this;
}

template<typename T>
inline Field3d<T>& Field3d<T>::operator=(Field3d<T>&& temp_field2)
{
	_data = temp_field2._data;
	temp_field2._data = nullptr;
	return *this;
}

template<typename T>
inline Field3d<T>& Field3d<T>::operator=(const T& c)
{
	for (size_t i = 0; i < (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1); i++)
	{
		_data[i] = c;
	}
	return *this;
}

template<typename T>
inline Field3d<T>& Field3d<T>::operator+=(const Field3d<T>& field2)
{
	for (size_t i = 0; i < (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1); i++)
	{
		_data[i] += field2._data[i];
	}
	return *this;
}

template<typename T>
inline Field3d<T>& Field3d<T>::operator+=(const T& c)
{
	for (size_t i = 0; i < (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1); i++)
	{
		_data[i] += c;
	}
	return *this;
}

template<typename T>
inline Field3d<T>& Field3d<T>::operator*=(const Field3d<T>& field2)
{
	for (size_t i = 0; i < (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1); i++)
	{
		_data[i] *= field2._data[i];
	}
	return *this;
}

template<typename T>
inline Field3d<T>& Field3d<T>::operator*=(const T& c)
{
	for (size_t i = 0; i < (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1); i++)
	{
		_data[i] *= c;
	}
	return *this;
}

template<typename T>
inline Field3d<T>& Field3d<T>::operator-=(const Field3d<T>& field2)
{
	for (size_t i = 0; i < (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1); i++)
	{
		_data[i] -= field2._data[i];
	}
	return *this;
}

template<typename T>
inline Field3d<T>& Field3d<T>::operator-=(const T& c)
{
	for (size_t i = 0; i < (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1); i++)
	{
		_data[i] -= c;
	}
	return *this;
}

template<typename T>
inline Field3d<T>& Field3d<T>::operator/=(const Field3d<T>& field2)
{
	for (size_t i = 0; i < (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1); i++)
	{
		_data[i] /= field2._data[i];
	}
	return *this;
}

template<typename T>
inline Field3d<T>& Field3d<T>::operator/=(const T& c)
{
	for (size_t i = 0; i < (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1); i++)
	{
		_data[i] /= c;
	}
	return *this;
}

template<typename T>
inline void Field3d<T>::setRand()
{
	srand(time(nullptr));
	for (size_t i = 0; i < (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1); i++)
	{
		_data[i] = (T)(rand() / (double)(RAND_MAX));
	}
}

template<typename T>
inline void Field3d<T>::setZero()
{
	for (size_t i = 0; i < (_up1 - _dn1 + 1) * (_up2 - _dn2 + 1) * (_up3 - _dn3 + 1); i++)
	{
		_data[i] = (T)(0);
	}
}