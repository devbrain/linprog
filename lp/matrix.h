#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <assert.h>
#include <algorithm>
#include <memory>
#include <vector>
#include <iostream>

template<class T, class A = std::allocator<T> >
struct Matrix {
	typedef T value_type;
	typedef std::vector<value_type, A> Container;

	Matrix() : _b(0) {}
	Matrix(int rows, int cols, value_type const& initial = value_type())
		: _b(0)
	{
		resize(rows, cols, initial);
	}
	Matrix(Matrix const& other)
		: _data(other._data), _b(other._b)
	{}
	
	Matrix& operator=(Matrix copy) {
		swap(*this, copy);
		return *this;
	}

	bool empty() const { return _data.empty(); }
	void clear() { _data.clear(); _b = 0; }

	int numRows() const { return _b ? _data.size() / _b : 0; }
	int numCols() const { return _b; }

	bool rowIsNull(int row) {
		for (int i = 0; i<numCols(); i++)
			if (_data[row * _b + i] != 0)
				return false;
		return true;
	}

	bool colIsNull(int col) {
		for (int i = 0; i<numRows(); i++)
			if (_data[i * _b + col] != 0)
				return false;
		return true;
	}

	value_type& operator()(int i) {
		return _data[i];
	}

	const value_type& operator()(int i) const {
		return _data[i];
	}

	value_type& operator()(int row, int col) {
		return _data[row * _b + col];
	}

	const value_type& operator()(int row, int col) const {
		return _data[row * _b + col];
	}

	void resize(int a, int b, value_type const& initial = value_type()) {
		if (a == 0) {
			b = 0;
		}
		_data.resize(a * b, initial);
		_b = b;
	}

	friend void swap(Matrix& a, Matrix& b) {
		using std::swap;
		swap(a._data, b._data);
		swap(a._b, b._b);
	}

	std::ostream& print(std::ostream& s)
	{
		s << "<Matrix" << numRows() << 'x' << numCols();
		if (!empty())
		{
			bool first = true;
			typedef typename Container::const_iterator Iter;
			Iter i = _data.begin(), end = _data.end();
			while (i != end)
			{
				s << (first ? " [[" : "], [");
				first = false;
				s << *i;
				++i;
				for (int b = _b - 1; b; --b) {
					s << ", " << *i;
					++i;
				}
			}
			s << "]]";
		}
		s << '>' << std::endl;
		return s;
	}

private:
	Container _data;
	int _b;
};


template <class M>
struct MatrixBuilder
{
	typedef typename M::value_type value_type;
	typedef M return_type;

	MatrixBuilder(int a) :_m(a, a) {
	}

	MatrixBuilder(int a, int b) :_m(a, b) {
	}

	MatrixBuilder(const M& m) :_m(m) {
	}

	MatrixBuilder& zero() {
		_m = M(_m.numRows(), _m.numCols(), value_type(0));
		return *this;
	}


	MatrixBuilder& random() {
		_m = M(_m.numRows(), _m.numCols(), value_type(0));
		for (int i = 0; i<_m.numRows(); i++)
			for (int j = 0; j<_m.numCols(); j++)
				_m(i, j) = rand();
		return *this;
	}

	MatrixBuilder& identity() {
		assert(_m.numRows() == _m.numCols());
		_m = M(_m.numRows(), _m.numCols(), value_type(0));
		for (int i = 0; i<_m.numRows(); i++) {
			_m(i, i) = 1;
		}
		return *this;
	}

	MatrixBuilder& addToRow(int row, value_type v) {
		assert(!_m.empty());
		assert(row <= _m.numRows());
		for (int i = 0; i<_m.numCols(); i++) {
			_m(row, i) += v;
		}
		return *this;
	}

	MatrixBuilder& addToColumn(int col, value_type v) {
		assert(!_m.empty());
		assert(col <= _m.numCols());
		for (int i = 0; i<_m.numRows(); i++) {
			_m(i, col) += v;
		}
		return *this;
	}

	MatrixBuilder& multiplyRowBy(int row, value_type v) {
		assert(!_m.empty());
		assert(row <= _m.numRows());
		for (int i = 0; i<_m.numCols(); i++) {
			_m(row, i) *= v;
		}
		return *this;
	}

	MatrixBuilder& multiplyColumnBy(int col, value_type v) {
		assert(!_m.empty());
		assert(col <= _m.numCols());
		for (int i = 0; i<_m.numRows(); i++) {
			_m(i, col) *= v;
		}
		return *this;
	}

	MatrixBuilder& divideRowBy(int row, value_type v) {
		assert(!_m.empty());
		assert(row <= _m.numRows());
		for (int i = 0; i<_m.numCols(); i++) {
			_m(row, i) /= v;
		}
		return *this;
	}

	MatrixBuilder& divideColumnBy(int col, value_type v) {
		assert(!_m.empty());
		assert(col <= _m.numCols());
		for (int i = 0; i<_m.numRows(); i++) {
			_m(i, col) /= v;
		}
		return *this;
	}

	MatrixBuilder& swapRows(int row1, int row2) {
		assert(!_m.empty());
		assert(row1 <= _m.numRows());
		assert(row2 <= _m.numRows());

		for (int i = 0; i<_m.numCols(); i++) {
			std::swap(_m(row1, i), _m(row2, i));
		}
		return *this;
	}

	MatrixBuilder& swapColumns(int col1, int col2) {
		assert(!_m.empty());
		assert(col1 <= _m.numCols());
		assert(col2 <= _m.numCols());

		for (int i = 0; i<_m.numRows(); i++) {
			std::swap(_m(i, col1), _m(i, col2));
		}
		return *this;
	}

	MatrixBuilder& addRowElements(int row1, int row2, value_type multiplier) {
		assert(!_m.empty());
		assert(row1 <= _m.numRows());
		assert(row2 <= _m.numRows());

		for (int i = 0; i<_m.numCols(); i++) {
			_m(row1, i) += _m(row2, i) * multiplier;
		}
		return *this;
	}

	MatrixBuilder& subRowElements(int row1, int row2, value_type multiplier) {
		assert(!_m.empty());
		assert(row1 <= _m.numRows());
		assert(row2 <= _m.numRows());

		for (int i = 0; i<_m.numCols(); i++) {
			_m(row1, i) -= _m(row2, i) * multiplier;
		}
		return *this;
	}

	M result() const  {
		return _m;
	}

private:
	M _m;
};


template <class M>
struct GuassJordanSolver
{
	typedef typename M::value_type value_type;
	typedef M return_type;

	GuassJordanSolver(const M& m) :_m(m), _b(m), _i(m.numRows()) {
		assert(!_m.empty());
		assert(_m.numRows() == _m.numCols());
		_i.identity();
	}

	bool solve()
	{
		for (int i = 0; i<_m.numCols(); i++)
		{
			value_type v = _b.result()(i, i);

			// Interchange rows to bring a non-zero element, only look for rows greater than current
			if (v == 0 && !fixNullPivot(i, i))
			{
				return false;
			}

			v = _b.result()(i, i);
			_b.divideRowBy(i, v);
			_i.divideRowBy(i, v);

			//_b.result().print(std::cout);

			for (int j = 0; j<_m.numRows(); j++)
			{
				if (i == j) continue;

				v = _b.result()(j, i);

				if (v == 0) continue;

				_b.subRowElements(j, i, v);
				_i.subRowElements(j, i, v);

				//_b.result().print(std::cout);
			}

		}

		return true;
	}

	M original() const  {
		return _m;
	}

	M transformed() const  {
		return _b.result();
	}

	M inverse() const  {
		return _i.result();
	}

private:
	bool fixNullPivot(int r, int c)
	{
		for (int i = r + 1; i<_m.numRows(); i++)
		{
			if (_b.result()(i, c) != 0)
			{
				_b.swapRows(r, i);
				_i.swapRows(r, i);
			//	_b.result().print(std::cout);
				return true;
			}
		}
		return false;
	}

private:
	M _m;
	MatrixBuilder<return_type> _b;
	MatrixBuilder<return_type> _i;
};

template<class T, class A = std::allocator<T> >
Matrix <T, A> invert(const Matrix <T, A>& m)
{
	GuassJordanSolver <Matrix <T, A> > s(m);
	bool rc = s.solve();
	assert(rc);
	return s.inverse();
}


#endif

