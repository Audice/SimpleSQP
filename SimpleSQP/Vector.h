#pragma once
#include <set>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>


namespace QuadProg
{

	enum MType { DIAG };

	template <typename T>
	class Vector
	{
	public:
		Vector();
		Vector(const unsigned int n);
		Vector(const T& a, const unsigned int n); //initialize to constant value 
		Vector(const T* a, const unsigned int n); // Initialize to array 
		Vector(const Vector& rhs); // copy constructor 
		~Vector(); // destructor

		inline void set(const T* a, const unsigned int n);
		Vector<T> extract(const std::set<unsigned int>& indexes) const;
		inline T& operator[](const unsigned int& i); //i-th element 
		inline const T& operator[](const unsigned int& i) const;

		inline unsigned int size() const;
		inline void resize(const unsigned int n);
		inline void resize(const T& a, const unsigned int n);

		inline T* GetData() const;

		Vector<T>& operator=(const Vector<T>& rhs); //assignment 
		Vector<T>& operator=(const T& a); //assign a to every element 
		inline Vector<T>& operator+=(const Vector<T>& rhs);
		inline Vector<T>& operator-=(const Vector<T>& rhs);
		inline Vector<T>& operator*=(const Vector<T>& rhs);
		inline Vector<T>& operator/=(const Vector<T>& rhs);
		inline Vector<T>& operator^=(const Vector<T>& rhs);
		inline Vector<T>& operator+=(const T& a);
		inline Vector<T>& operator-=(const T& a);
		inline Vector<T>& operator*=(const T& a);
		inline Vector<T>& operator/=(const T& a);
		inline Vector<T>& operator^=(const T& a);
	private:
		unsigned int n; // size of array. upper index is n-1 
		T* v;           // storage for data
	};

	template <typename T>
	Vector<T>::Vector() : n(0), v(0)
	{}

	template <typename T>
	Vector<T>::Vector(const unsigned int n) : v(new T[n])
	{
		this->n = n;
	}

	template <typename T>
	Vector<T>::Vector(const T& a, const unsigned int n) : v(new T[n])
	{
		this->n = n;
		for (unsigned int i = 0; i < n; i++)
			v[i] = a;
	}

	template <typename T>
	Vector<T>::Vector(const T* a, const unsigned int n) : v(new T[n])
	{
		this->n = n;
		for (unsigned int i = 0; i < n; i++)
			v[i] = *a++;
	}

	template <typename T>
	Vector<T>::Vector(const Vector<T>& rhs) : v(new T[rhs.n])
	{
		this->n = rhs.n;
		for (unsigned int i = 0; i < n; i++)
			v[i] = rhs[i];
	}

	template <typename T>
	Vector<T>::~Vector()
	{
		if (v != 0)
			delete[](v);
	}

	template <typename T>
	void Vector<T>::resize(const unsigned int n)
	{
		if (n == this->n)
			return;
		if (v != 0)
			delete[](v);
		v = new T[n];
		this->n = n;
	}

	template <typename T>
	void Vector<T>::resize(const T& a, const unsigned int n)
	{
		resize(n);
		for (unsigned int i = 0; i < n; i++)
			v[i] = a;
	}


	template <typename T>
	inline Vector<T>& Vector<T>::operator=(const Vector<T>& rhs)
		// postcondition: normal assignment via copying has been performed; 
		// if vector and rhs were different sizes, vector 
		// has been resized to match the size of rhs 
	{
		if (this != &rhs)
		{
			resize(rhs.n);
			for (unsigned int i = 0; i < n; i++)
				v[i] = rhs[i];
		}
		return *this;
	}

	template <typename T>
	inline Vector<T>& Vector<T>::operator=(const T& a) //assign a to every element 
	{
		for (unsigned int i = 0; i < n; i++)
			v[i] = a;
		return *this;
	}

	template <typename T>
	inline T& Vector<T>::operator[](const unsigned int& i) //subscripting 
	{
		return v[i];
	}

	template <typename T>
	inline const T& Vector<T>::operator[](const unsigned int& i) const //subscripting 
	{
		return v[i];
	}

	template <typename T>
	inline unsigned int Vector<T>::size() const
	{
		return n;
	}

	template <typename T>
	inline void Vector<T>::set(const T* a, unsigned int n)
	{
		resize(n);
		for (unsigned int i = 0; i < n; i++)
			v[i] = a[i];
	}

	template <typename T>
	inline T* Vector<T>::GetData() const
	{
		T* dataCopy = new T[n];
		for (size_t i = 0; i < n; i++) {
			dataCopy[i] = v[i];
		}
		return dataCopy;
	}

	template <typename T>
	inline Vector<T> Vector<T>::extract(const std::set<unsigned int>& indexes) const
	{
		Vector<T> tmp(indexes.size());
		unsigned int i = 0;

		for (std::set<unsigned int>::const_iterator el = indexes.begin(); el != indexes.end(); el++)
		{
			if (*el >= n)
				throw std::logic_error("Error extracting subvector: the indexes are out of vector bounds");
			tmp[i++] = v[*el];
		}

		return tmp;
	}

	template <typename T>
	inline Vector<T>& Vector<T>::operator+=(const Vector<T>& rhs)
	{
		if (this->size() != rhs.size())
			throw std::logic_error("Operator+=: vectors have different sizes");
		for (unsigned int i = 0; i < n; i++)
			v[i] += rhs[i];

		return *this;
	}


	template <typename T>
	inline Vector<T>& Vector<T>::operator+=(const T& a)
	{
		for (unsigned int i = 0; i < n; i++)
			v[i] += a;

		return *this;
	}

	template <typename T>
	inline Vector<T> operator+(const Vector<T>& rhs)
	{
		return rhs;
	}

	template <typename T>
	inline Vector<T> operator+(const Vector<T>& lhs, const Vector<T>& rhs)
	{
		if (lhs.size() != rhs.size())
			throw std::logic_error("Operator+: vectors have different sizes");
		Vector<T> tmp(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++)
			tmp[i] = lhs[i] + rhs[i];

		return tmp;
	}

	template <typename T>
	inline Vector<T> operator+(const Vector<T>& lhs, const T& a)
	{
		Vector<T> tmp(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++)
			tmp[i] = lhs[i] + a;

		return tmp;
	}

	template <typename T>
	inline Vector<T> operator+(const T& a, const Vector<T>& rhs)
	{
		Vector<T> tmp(rhs.size());
		for (unsigned int i = 0; i < rhs.size(); i++)
			tmp[i] = a + rhs[i];

		return tmp;
	}

	template <typename T>
	inline Vector<T>& Vector<T>::operator-=(const Vector<T>& rhs)
	{
		if (this->size() != rhs.size())
			throw std::logic_error("Operator-=: vectors have different sizes");
		for (unsigned int i = 0; i < n; i++)
			v[i] -= rhs[i];

		return *this;
	}


	template <typename T>
	inline Vector<T>& Vector<T>::operator-=(const T& a)
	{
		for (unsigned int i = 0; i < n; i++)
			v[i] -= a;

		return *this;
	}

	template <typename T>
	inline Vector<T> operator-(const Vector<T>& rhs)
	{
		return (T)(-1) * rhs;
	}

	template <typename T>
	inline Vector<T> operator-(const Vector<T>& lhs, const Vector<T>& rhs)
	{
		if (lhs.size() != rhs.size())
			throw std::logic_error("Operator-: vectors have different sizes");
		Vector<T> tmp(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++)
			tmp[i] = lhs[i] - rhs[i];

		return tmp;
	}

	template <typename T>
	inline Vector<T> operator-(const Vector<T>& lhs, const T& a)
	{
		Vector<T> tmp(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++)
			tmp[i] = lhs[i] - a;

		return tmp;
	}

	template <typename T>
	inline Vector<T> operator-(const T& a, const Vector<T>& rhs)
	{
		Vector<T> tmp(rhs.size());
		for (unsigned int i = 0; i < rhs.size(); i++)
			tmp[i] = a - rhs[i];

		return tmp;
	}

	template <typename T>
	inline Vector<T>& Vector<T>::operator*=(const Vector<T>& rhs)
	{
		if (this->size() != rhs.size())
			throw std::logic_error("Operator*=: vectors have different sizes");
		for (unsigned int i = 0; i < n; i++)
			v[i] *= rhs[i];

		return *this;
	}


	template <typename T>
	inline Vector<T>& Vector<T>::operator*=(const T& a)
	{
		for (unsigned int i = 0; i < n; i++)
			v[i] *= a;

		return *this;
	}

	template <typename T>
	inline Vector<T> operator*(const Vector<T>& lhs, const Vector<T>& rhs)
	{
		if (lhs.size() != rhs.size())
			throw std::logic_error("Operator*: vectors have different sizes");
		Vector<T> tmp(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++)
			tmp[i] = lhs[i] * rhs[i];

		return tmp;
	}

	template <typename T>
	inline Vector<T> operator*(const Vector<T>& lhs, const T& a)
	{
		Vector<T> tmp(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++)
			tmp[i] = lhs[i] * a;

		return tmp;
	}

	template <typename T>
	inline Vector<T> operator*(const T& a, const Vector<T>& rhs)
	{
		Vector<T> tmp(rhs.size());
		for (unsigned int i = 0; i < rhs.size(); i++)
			tmp[i] = a * rhs[i];

		return tmp;
	}

	template <typename T>
	inline Vector<T>& Vector<T>::operator/=(const Vector<T>& rhs)
	{
		if (this->size() != rhs.size())
			throw std::logic_error("Operator/=: vectors have different sizes");
		for (unsigned int i = 0; i < n; i++)
			v[i] /= rhs[i];

		return *this;
	}


	template <typename T>
	inline Vector<T>& Vector<T>::operator/=(const T& a)
	{
		for (unsigned int i = 0; i < n; i++)
			v[i] /= a;

		return *this;
	}

	template <typename T>
	inline Vector<T> operator/(const Vector<T>& lhs, const Vector<T>& rhs)
	{
		if (lhs.size() != rhs.size())
			throw std::logic_error("Operator/: vectors have different sizes");
		Vector<T> tmp(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++)
			tmp[i] = lhs[i] / rhs[i];

		return tmp;
	}

	template <typename T>
	inline Vector<T> operator/(const Vector<T>& lhs, const T& a)
	{
		Vector<T> tmp(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++)
			tmp[i] = lhs[i] / a;

		return tmp;
	}

	template <typename T>
	inline Vector<T> operator/(const T& a, const Vector<T>& rhs)
	{
		Vector<T> tmp(rhs.size());
		for (unsigned int i = 0; i < rhs.size(); i++)
			tmp[i] = a / rhs[i];

		return tmp;
	}

	template <typename T>
	inline Vector<T> operator^(const Vector<T>& lhs, const Vector<T>& rhs)
	{
		if (lhs.size() != rhs.size())
			throw std::logic_error("Operator^: vectors have different sizes");
		Vector<T> tmp(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++)
			tmp[i] = pow(lhs[i], rhs[i]);

		return tmp;
	}

	template <typename T>
	inline Vector<T> operator^(const Vector<T>& lhs, const T& a)
	{
		Vector<T> tmp(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++)
			tmp[i] = pow(lhs[i], a);

		return tmp;
	}

	template <typename T>
	inline Vector<T> operator^(const T& a, const Vector<T>& rhs)
	{
		Vector<T> tmp(rhs.size());
		for (unsigned int i = 0; i < rhs.size(); i++)
			tmp[i] = pow(a, rhs[i]);

		return tmp;
	}

	template <typename T>
	inline Vector<T>& Vector<T>::operator^=(const Vector<T>& rhs)
	{
		if (this->size() != rhs.size())
			throw std::logic_error("Operator^=: vectors have different sizes");
		for (unsigned int i = 0; i < n; i++)
			v[i] = pow(v[i], rhs[i]);

		return *this;
	}

	template <typename T>
	inline Vector<T>& Vector<T>::operator^=(const T& a)
	{
		for (unsigned int i = 0; i < n; i++)
			v[i] = pow(v[i], a);

		return *this;
	}

	template <typename T>
	inline bool operator==(const Vector<T>& v, const Vector<T>& w)
	{
		if (v.size() != w.size())
			throw std::logic_error("Vectors of different size are not confrontable");
		for (unsigned i = 0; i < v.size(); i++)
			if (v[i] != w[i])
				return false;
		return true;
	}

	template <typename T>
	inline bool operator!=(const Vector<T>& v, const Vector<T>& w)
	{
		if (v.size() != w.size())
			throw std::logic_error("Vectors of different size are not confrontable");
		for (unsigned i = 0; i < v.size(); i++)
			if (v[i] != w[i])
				return true;
		return false;
	}

	template <typename T>
	inline bool operator<(const Vector<T>& v, const Vector<T>& w)
	{
		if (v.size() != w.size())
			throw std::logic_error("Vectors of different size are not confrontable");
		for (unsigned i = 0; i < v.size(); i++)
			if (v[i] >= w[i])
				return false;
		return true;
	}

	template <typename T>
	inline bool operator<=(const Vector<T>& v, const Vector<T>& w)
	{
		if (v.size() != w.size())
			throw std::logic_error("Vectors of different size are not confrontable");
		for (unsigned i = 0; i < v.size(); i++)
			if (v[i] > w[i])
				return false;
		return true;
	}

	template <typename T>
	inline bool operator>(const Vector<T>& v, const Vector<T>& w)
	{
		if (v.size() != w.size())
			throw std::logic_error("Vectors of different size are not confrontable");
		for (unsigned i = 0; i < v.size(); i++)
			if (v[i] <= w[i])
				return false;
		return true;
	}

	template <typename T>
	inline bool operator>=(const Vector<T>& v, const Vector<T>& w)
	{
		if (v.size() != w.size())
			throw std::logic_error("Vectors of different size are not confrontable");
		for (unsigned i = 0; i < v.size(); i++)
			if (v[i] < w[i])
				return false;
		return true;
	}

	/**
	Input/Output
	*/
	template <typename T>
	inline std::ostream& operator<<(std::ostream& os, const Vector<T>& v)
	{
		os << std::endl << v.size() << std::endl;
		for (unsigned int i = 0; i < v.size() - 1; i++)
			os << std::setw(20) << std::setprecision(16) << v[i] << ", ";
		os << std::setw(20) << std::setprecision(16) << v[v.size() - 1] << std::endl;

		return os;
	}

	template <typename T>
	std::istream& operator>>(std::istream& is, Vector<T>& v)
	{
		int elements;
		char comma;
		is >> elements;
		v.resize(elements);
		for (unsigned int i = 0; i < elements; i++)
			is >> v[i] >> comma;

		return is;
	}

	/**
	Index utilities
	*/

	std::set<unsigned int> seq(unsigned int s, unsigned int e);

	std::set<unsigned int> singleton(unsigned int i);

	template <typename T>
	class CanonicalBaseVector : public Vector<T>
	{
	public:
		CanonicalBaseVector(unsigned int i, unsigned int n);
		inline void reset(unsigned int i);
	private:
		unsigned int e;
	};

	template <typename T>
	CanonicalBaseVector<T>::CanonicalBaseVector(unsigned int i, unsigned int n)
		: Vector<T>((T)0, n), e(i)
	{
		(*this)[e] = (T)1;
	}

	template <typename T>
	inline void CanonicalBaseVector<T>::reset(unsigned int i)
	{
		(*this)[e] = (T)0;
		e = i;
		(*this)[e] = (T)1;
	}

#include <stdexcept>

	template <typename T>
	inline T sum(const Vector<T>& v)
	{
		T tmp = (T)0;
		for (unsigned int i = 0; i < v.size(); i++)
			tmp += v[i];

		return tmp;
	}

	template <typename T>
	inline T prod(const Vector<T>& v)
	{
		T tmp = (T)1;
		for (unsigned int i = 0; i < v.size(); i++)
			tmp *= v[i];

		return tmp;
	}

	template <typename T>
	inline T mean(const Vector<T>& v)
	{
		T sum = (T)0;
		for (unsigned int i = 0; i < v.size(); i++)
			sum += v[i];
		return sum / v.size();
	}

	template <typename T>
	inline T median(const Vector<T>& v)
	{
		Vector<T> tmp = sort(v);
		if (v.size() % 2 == 1) // it is an odd-sized vector
			return tmp[v.size() / 2];
		else
			return 0.5 * (tmp[v.size() / 2 - 1] + tmp[v.size() / 2]);
	}

	template <typename T>
	inline T stdev(const Vector<T>& v, bool sample_correction = false)
	{
		return sqrt(var(v, sample_correction));
	}

	template <typename T>
	inline T var(const Vector<T>& v, bool sample_correction = false)
	{
		T sum = (T)0, ssum = (T)0;
		unsigned int n = v.size();
		for (unsigned int i = 0; i < n; i++)
		{
			sum += v[i];
			ssum += (v[i] * v[i]);
		}
		if (!sample_correction)
			return (ssum / n) - (sum / n) * (sum / n);
		else
			return n * ((ssum / n) - (sum / n) * (sum / n)) / (n - 1);
	}

	template <typename T>
	inline T max(const Vector<T>& v)
	{
		T value = v[0];
		for (unsigned int i = 1; i < v.size(); i++)
			value = std::max(v[i], value);

		return value;
	}

	template <typename T>
	inline T min(const Vector<T>& v)
	{
		T value = v[0];
		for (unsigned int i = 1; i < v.size(); i++)
			value = std::min(v[i], value);

		return value;
	}

	template <typename T>
	inline unsigned int index_max(const Vector<T>& v)
	{
		unsigned int max = 0;
		for (unsigned int i = 1; i < v.size(); i++)
			if (v[i] > v[max])
				max = i;

		return max;
	}

	template <typename T>
	inline unsigned int index_min(const Vector<T>& v)
	{
		unsigned int min = 0;
		for (unsigned int i = 1; i < v.size(); i++)
			if (v[i] < v[min])
				min = i;

		return min;
	}


	template <typename T>
	inline T dot_prod(const Vector<T>& a, const Vector<T>& b)
	{
		T sum = (T)0;
		if (a.size() != b.size())
			throw std::logic_error("Dotprod error: the vectors are not the same size");
		for (unsigned int i = 0; i < a.size(); i++)
			sum += a[i] * b[i];

		return sum;
	}

	/**
	Single element mathematical functions
	*/

	template <typename T>
	inline Vector<T> exp(const Vector<T>& v)
	{
		Vector<T> tmp(v.size());
		for (unsigned int i = 0; i < v.size(); i++)
			tmp[i] = exp(v[i]);

		return tmp;
	}

	template <typename T>
	inline Vector<T> log(const Vector<T>& v)
	{
		Vector<T> tmp(v.size());
		for (unsigned int i = 0; i < v.size(); i++)
			tmp[i] = log(v[i]);

		return tmp;
	}

	template <typename T>
	inline Vector<T> sqrt(const Vector<T>& v)
	{
		Vector<T> tmp(v.size());
		for (unsigned int i = 0; i < v.size(); i++)
			tmp[i] = sqrt(v[i]);

		return tmp;
	}

	template <typename T>
	inline Vector<T> pow(const Vector<T>& v, double a)
	{
		Vector<T> tmp(v.size());
		for (unsigned int i = 0; i < v.size(); i++)
			tmp[i] = pow(v[i], a);

		return tmp;
	}

	template <typename T>
	inline Vector<T> abs(const Vector<T>& v)
	{
		Vector<T> tmp(v.size());
		for (unsigned int i = 0; i < v.size(); i++)
			tmp[i] = (T)fabs(v[i]);

		return tmp;
	}

	template <typename T>
	inline Vector<T> sign(const Vector<T>& v)
	{
		Vector<T> tmp(v.size());
		for (unsigned int i = 0; i < v.size(); i++)
			tmp[i] = v[i] > 0 ? +1 : v[i] == 0 ? 0 : -1;

		return tmp;
	}

	template <typename T>
	inline unsigned int partition(Vector<T>& v, unsigned int begin, unsigned int end)
	{
		unsigned int i = begin + 1, j = begin + 1;
		T pivot = v[begin];
		while (j <= end)
		{
			if (v[j] < pivot) {
				std::swap(v[i], v[j]);
				i++;
			}
			j++;
		}
		v[begin] = v[i - 1];
		v[i - 1] = pivot;
		return i - 2;
	}


	template <typename T>
	inline void quicksort(Vector<T>& v, unsigned int begin, unsigned int end)
	{
		if (end > begin)
		{
			unsigned int index = partition(v, begin, end);
			quicksort(v, begin, index);
			quicksort(v, index + 2, end);
		}
	}

	template <typename T>
	inline Vector<T> sort(const Vector<T>& v)
	{
		Vector<T> tmp(v);

		quicksort<T>(tmp, 0, tmp.size() - 1);

		return tmp;
	}

	template <typename T>
	inline Vector<double> rank(const Vector<T>& v)
	{
		Vector<T> tmp(v);
		Vector<double> tmp_rank(0.0, v.size());

		for (unsigned int i = 0; i < tmp.size(); i++)
		{
			unsigned int smaller = 0, equal = 0;
			for (unsigned int j = 0; j < tmp.size(); j++)
				if (i == j)
					continue;
				else
					if (tmp[j] < tmp[i])
						smaller++;
					else if (tmp[j] == tmp[i])
						equal++;
			tmp_rank[i] = smaller + 1;
			if (equal > 0)
			{
				for (unsigned int j = 1; j <= equal; j++)
					tmp_rank[i] += smaller + 1 + j;
				tmp_rank[i] /= (double)(equal + 1);
			}
		}

		return tmp_rank;
	}
}

