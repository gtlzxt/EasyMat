#include <iostream>
#include <math.h>
#include "EasyMat.h"
#include "EasyVec.h"
using namespace em;

const double EasyMat::mPrecision = 0.000001;

EasyMat::EasyMat(): mCols(0),mRows(0), mCapacity(0), mDim(EM_BY_ROW),mData(nullptr)
{

}

EasyMat::EasyMat(unsigned long cols): mCols(cols),
	mRows(0), mDim(EM_BY_ROW)
{
	mCapacity = cols > INIT_CAPACITY ? cols : INIT_CAPACITY;
	initMem();
}

EasyMat::EasyMat(unsigned long rows, unsigned long cols) : mCols(cols),
	mRows(rows), mDim(EM_BY_ROW)
{
	unsigned long totalCount = rows * cols;
	mCapacity = totalCount > INIT_CAPACITY ? totalCount : INIT_CAPACITY;
	initMem();
}

EasyMat::EasyMat(const EasyMat & rhs)
{
	mData = nullptr;
	operator=(rhs);
}

EasyMat::EasyMat(EasyMat&& rhs):mData(rhs.mData),mRows(rhs.mRows),
	mCols(rhs.mCols),mCapacity(rhs.mCapacity), mDim(rhs.mDim)
{
	rhs.mData = nullptr;
}

const EasyMat & EasyMat::operator=(const EasyMat & rhs)
{
	if (this == &rhs)
		return *this;
	mCols = rhs.mCols;
	mRows = rhs.mRows;
	mCapacity = rhs.mCapacity;
	mDim = rhs.mDim;

	if(mData != nullptr)
		delete[] mData;
	mData = new double[mCapacity];
	memcpy(mData, rhs.mData, rhs.mCapacity * sizeof(double));
	return *this;
}

const EasyMat& EasyMat::operator=(EasyMat&& rhs)
{
	if (this == &rhs)
		return *this;
	mCols = rhs.mCols;
	mRows = rhs.mRows;
	mCapacity = rhs.mCapacity;
	mDim = rhs.mDim;
	mData = rhs.mData;
	rhs.mData = nullptr;
	return *this;
}

EasyMat EasyMat::operator+(const EasyMat & rhs) const
{
	if (mCols != rhs.mCols || mRows != rhs.mRows)
		throw "add can only be applied to Mat with same cols and rows";
	EasyMat result(*this);
	result += rhs;
	return result;
}

EasyMat EasyMat::operator+(double rhs) const
{
	EasyMat result(*this);
	result += rhs;
	return result;
}

EasyMat EasyMat::operator-(const EasyMat & rhs) const
{
	if (mCols != rhs.mCols || mRows != rhs.mRows)
		throw "minus can only be applied to Mat with same cols and rows";
	EasyMat result(*this);
	result -= rhs;
	return result;
}

EasyMat EasyMat::operator-(double rhs) const
{
	EasyMat result(*this);
	result -= rhs;
	return result;
}

EasyMat EasyMat::operator*(const EasyMat & rhs) const
{
	if (mCols != rhs.mRows)
		throw "multi can not be applied.";
	EasyMat result(*this);
	result *= rhs;
	return result;
}

EasyMat EasyMat::operator*(double rhs) const
{
	EasyMat result(*this);
	result *= rhs;
	return result;
}

EasyMat EasyMat::operator/(double rhs) const
{
	EasyMat result(*this);
	result /= rhs;
	return result;
}

EasyMat& EasyMat::operator+=(const EasyMat & rhs)
{
	if (mCols != rhs.mCols || mRows != rhs.mRows)
		throw "add can only be applied to Mat with same cols and rows";
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < mCols; j++)
		{
			operator()(i, j) = operator()(i, j) + rhs(i, j);
		}
	}
	return *this;
}

EasyMat& EasyMat::operator+=(double rhs)
{
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < mCols; j++)
		{
			operator()(i, j) = operator()(i, j) + rhs;
		}
	}
	return *this;
}

EasyMat& EasyMat::operator-=(const EasyMat & rhs)
{
	if (mCols != rhs.mCols || mRows != rhs.mRows)
		throw "minus can only be applied to Mat with same cols and rows";
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < mCols; j++)
		{
			operator()(i, j) = operator()(i, j) - rhs(i, j);
		}
	}
	return *this;
}

EasyMat& EasyMat::operator-=(double rhs)
{
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < mCols; j++)
		{
			operator()(i, j) = operator()(i, j) - rhs;
		}
	}
	return *this;
}

EasyMat& EasyMat::operator*=(const EasyMat & rhs)
{
	if (mCols != rhs.mRows)
		throw "multi can not be applied.";
	EasyMat tmp(mRows, rhs.mCols);
	double sum;
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < rhs.mCols; j++)
		{
			sum = 0;
			for (unsigned long k = 0; k < mCols; k++)
			{
				sum += operator()(i, k) * rhs(k, j);
			}
			tmp(i, j) = sum;
		}
	}

	*this = std::move(tmp);
	return *this;
}

EasyMat& EasyMat::operator*=(double rhs)
{
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < mCols; j++)
		{
			operator()(i, j) = operator()(i, j) * rhs;
		}
	}
	return *this;
}

EasyMat& EasyMat::operator/=(double rhs)
{
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < mCols; j++)
		{
			operator()(i, j) = operator()(i, j) / rhs;
		}
	}
	return *this;
}

bool EasyMat::operator==(const EasyMat& rhs) const
{
	if (mRows != rhs.rows() || mCols != rhs.cols() || mDim != rhs.mDim)
		return false;
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < mCols; j++)
		{
			if (std::abs(mData[i * mCols + j] - rhs(i, j)) > mPrecision)
				return false;
		}
	}	
	return true;
}

bool EasyMat::operator!=(const EasyMat& rhs) const
{
	return !operator==(rhs);
}

EasyMat EasyMat::transpose()
{
	EasyMat result(mCols, mRows);
	for (unsigned long i = 0; i < mCols; i++)
	{
		for (unsigned long j = 0; j < mRows; j++)
		{
			result(i, j) = operator()(j, i);
		}
	}
	return result;
}

double & EasyMat::operator()(unsigned long i, unsigned long j)
{
	if (i >= mRows || j >= mCols)
		throw "index outof range";
	return mData[i * mCols + j];
}

const double & EasyMat::operator()(unsigned long i, unsigned long j) const
{
	if (i >= mRows || j >= mCols)
		throw "index outof range";
	return mData[i * mCols + j];
}

void EasyMat::append(const EasyVec & vec, EM_DIMENSION dim)
{
	//add a row
	if (dim == EM_BY_ROW)
	{
		insert(mRows, vec, dim);
	}
	else
	{
		insert(mCols, vec, dim);
	}
}

void EasyMat::insert(unsigned long index, const EasyVec & vec, EM_DIMENSION dim)
{
	//insert a row
	if (dim == EM_BY_ROW)
	{
		if (vec.getDimension() == EM_BY_COL)
			throw "vec and mat not fit, you may need to transpose the vec";
		if (mCols != vec.size())
			throw "vec size and mat Cols not fit";
		if (index > mRows)
			throw "index outof range";

		unsigned long totalCount = mCols * mRows;
		if ((totalCount + mCols) > mCapacity)
			resize((totalCount + mCols) * 2);

		unsigned long dataIndex = index * mCols;
		if (index != mRows)
		{			
			memcpy(mData + dataIndex + mCols,mData + dataIndex, 
				(totalCount - dataIndex) * sizeof(double));
		}		
		memcpy(mData + dataIndex, vec.mData, mCols * sizeof(double));
		++mRows;
	}
	else
	{
		if (vec.getDimension() == EM_BY_ROW)
			throw "vec and mat not fit, you may need to transpose the vec";
		if (mRows != vec.size())
			throw "vec size and mat Rows not fit";
		if (index > mCols)
			throw "index outof range";

		unsigned long newCols = mCols + 1;
		if (newCols * mRows > mCapacity)
			resize(newCols * mRows * 2);
		moveCol(index, 1);
		for (unsigned long i = 0; i < mRows; i++)
		{
			mData[i * newCols + index] = vec[i];
		}
		mCols = newCols;
	}
}

void EasyMat::remove(unsigned long index, EM_DIMENSION dim)
{
	//remove a row
	if (dim == EM_BY_ROW)
	{
		if (index > mRows)
			throw "index outof range";

		if (index != mRows)
		{
			unsigned long dataIndex = index * mCols;
			memcpy(mData + dataIndex,mData + dataIndex + mCols, 
				(mCols * mRows - dataIndex - mCols)*sizeof(double));
		}
		--mRows;
	}
	else
	{
		if (index > mCols)
			throw "index outof range";

		moveCol(index+1, -1);
		--mCols;
	}
}

void EasyMat::set(unsigned long index, const EasyVec& vec, EM_DIMENSION dim)
{
	//set a row
	if (dim == EM_BY_ROW)
	{
		if (vec.getDimension() == EM_BY_COL)
			throw "vec and mat not fit, you may need to transpose the vec";
		if (mCols != vec.size())
			throw "vec size and mat Cols not fit";
		if (index > mRows)
			throw "index outof range";

		memcpy(mData + index * mCols, vec.mData, mCols * sizeof(double));
	}
	else
	{
		if (vec.getDimension() == EM_BY_ROW)
			throw "vec and mat not fit, you may need to transpose the vec";
		if (mRows != vec.size())
			throw "vec size and mat Rows not fit";
		if (index > mCols)
			throw "index outof range";

		for (unsigned long i = 0; i < mRows; i++)
		{
			operator()(i, index) = vec[i];
		}
	}
}

void EasyMat::set(unsigned long index, double value, EM_DIMENSION dim)
{
	//set a row
	if (dim == EM_BY_ROW)
	{
		if (index > mRows)
			throw "index outof range";

		for (unsigned long i = 0; i < mCols; i++)
		{
			operator()(index, i) = value;
		}
	}
	else
	{
		if (index > mCols)
			throw "index outof range";

		for (unsigned long i = 0; i < mRows; i++)
		{
			operator()(i, index) = value;
		}
	}
}

void EasyMat::add(const EasyVec & vec, EM_DIMENSION dim)
{
	//add a value to a row
	if (dim == EM_BY_ROW)
	{
		if (vec.getDimension() == EM_BY_COL)
			throw "vec and mat not fit, you may need to transpose the vec";
		if (mCols != vec.size())
			throw "vec size and mat Cols not fit";

		for (unsigned long i = 0; i < mRows; i++)
		{
			for (unsigned long j = 0; j < mCols; j++)
			{
				operator()(i, j) += vec[j];
			}
		}
	}
	else
	{
		if (vec.getDimension() == EM_BY_ROW)
			throw "vec and mat not fit, you may need to transpose the vec";
		if (mRows != vec.size())
			throw "vec size and mat Rows not fit";

		for (unsigned long i = 0; i < mCols; i++)
		{
			for (unsigned long j = 0; j < mRows; j++)
			{
				operator()(j, i) += vec[j];
			}
		}
	}
}

void EasyMat::sub(const EasyVec& vec, EM_DIMENSION dim)
{
	//add a value to a row
	if (dim == EM_BY_ROW)
	{
		if (vec.getDimension() == EM_BY_COL)
			throw "vec and mat not fit, you may need to transpose the vec";
		if (mCols != vec.size())
			throw "vec size and mat Cols not fit";

		for (unsigned long i = 0; i < mRows; i++)
		{
			for (unsigned long j = 0; j < mCols; j++)
			{
				operator()(i, j) -= vec[j];
			}
		}
	}
	else
	{
		if (vec.getDimension() == EM_BY_ROW)
			throw "vec and mat not fit, you may need to transpose the vec";
		if (mRows != vec.size())
			throw "vec size and mat Rows not fit";

		for (unsigned long i = 0; i < mCols; i++)
		{
			for (unsigned long j = 0; j < mRows; j++)
			{
				operator()(j, i) -= vec[j];
			}
		}
	}
}

void EasyMat::mul(const EasyVec& vec, EM_DIMENSION dim)
{
	//add a value to a row
	if (dim == EM_BY_ROW)
	{
		if (vec.getDimension() == EM_BY_COL)
			throw "vec and mat not fit, you may need to transpose the vec";
		if (mCols != vec.size())
			throw "vec size and mat Cols not fit";

		for (unsigned long i = 0; i < mRows; i++)
		{
			for (unsigned long j = 0; j < mCols; j++)
			{
				operator()(i, j) *= vec[j];
			}
		}
	}
	else
	{
		if (vec.getDimension() == EM_BY_ROW)
			throw "vec and mat not fit, you may need to transpose the vec";
		if (mRows != vec.size())
			throw "vec size and mat Rows not fit";

		for (unsigned long i = 0; i < mCols; i++)
		{
			for (unsigned long j = 0; j < mRows; j++)
			{
				operator()(j, i) *= vec[j];
			}
		}
	}
}

void EasyMat::div(const EasyVec& vec, EM_DIMENSION dim)
{
	//add a value to a row
	if (dim == EM_BY_ROW)
	{
		if (vec.getDimension() == EM_BY_COL)
			throw "vec and mat not fit, you may need to transpose the vec";
		if (mCols != vec.size())
			throw "vec size and mat Cols not fit";

		for (unsigned long i = 0; i < mRows; i++)
		{
			for (unsigned long j = 0; j < mCols; j++)
			{
				operator()(i, j) /= vec[j];
			}
		}
	}
	else
	{
		if (vec.getDimension() == EM_BY_ROW)
			throw "vec and mat not fit, you may need to transpose the vec";
		if (mRows != vec.size())
			throw "vec size and mat Rows not fit";

		for (unsigned long i = 0; i < mCols; i++)
		{
			for (unsigned long j = 0; j < mRows; j++)
			{
				operator()(j, i) /= vec[j];
			}
		}
	}
}

void EasyMat::swap(unsigned long index1, unsigned long index2, EM_DIMENSION dim)
{
	//swap two row
	if (dim == EM_BY_ROW)
	{
		unsigned long dataSize = mCols * sizeof(double);
		double* buf = new double[dataSize];
		memcpy(buf, mData + index1 * mCols, dataSize);
		memcpy(mData + index1 * mCols, mData + index2 * mCols, dataSize);
		memcpy(mData + index2 * mCols, buf, dataSize);
		delete[] buf;
	}
	else
	{
		double tmp;
		for (unsigned long i = 0; i < mRows; i++)
		{
			tmp = operator()(i, index1);
			operator()(i, index1) = operator()(i, index2);
			operator()(i, index2) = tmp;
		}
	}
}


void EasyMat::sortByRowAsc(unsigned long index)
{
	//use choose sort, avoid to much swap
	for (unsigned long i = 0; i < mRows; i++)
	{
		unsigned long smallestIndex = i;
		for (unsigned long j = i; j < mRows; j++)
		{
			if (operator()(smallestIndex,index) > operator()(j, index))
				smallestIndex = j;
		}
		if (i != smallestIndex)
		{
			swap(i, smallestIndex, EM_BY_ROW);
		}
	}
}

void EasyMat::sortByRowDesc(unsigned long index)
{
	//use choose sort, avoid to much swap
	for (unsigned long i = 0; i < mRows; i++)
	{
		unsigned long BiggestIndex = i;
		for (unsigned long j = i; j < mRows; j++)
		{
			if (operator()(BiggestIndex, index) < operator()(j, index))
				BiggestIndex = j;
		}
		if (i != BiggestIndex)
		{
			swap(i, BiggestIndex, EM_BY_ROW);
		}
	}
}

EasyMat EasyMat::col(const EasyVec & vec) const
{
	EasyMat result(mRows, vec.size());
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < vec.size(); j++)
		{
			result(i, j) = operator()(i, (unsigned long)vec[j]);
		}
	}
	return result;
}

EasyMat EasyMat::row(const EasyVec & vec) const
{
	EasyMat result(vec.size(), mCols);
	for (unsigned long i = 0; i < vec.size(); i++)
	{
		memcpy(result.mData + i * mCols, 
			mData + (unsigned long)vec[i] * mCols, mCols * sizeof(double));
	}
	return result;
}

EasyVec EasyMat::col(unsigned long index) const
{
	EasyVec result(mRows, EM_BY_COL);
	for (unsigned long i = 0; i < mRows;i++)
	{
		result[i] = operator()(i, index);
	}
	return result;
}

EasyVec EasyMat::row(unsigned long index) const
{
	EasyVec result(mCols, EM_BY_ROW);
	memcpy(result.mData, mData + index * mCols, mCols * sizeof(double));
	return result;
}

double EasyMat::max(unsigned long index, EM_DIMENSION dim) const
{
	double result;
	if(dim == EM_BY_COL)
	{
		if (index >= mCols)
			throw "index out of range";
		if (mCols == 0 || mRows == 0)
			throw "empty mat";
		result = operator()(0, index);
		for (unsigned long i = 1; i < mRows; i++)
		{
			if (result < operator()(i, index))
				result = operator()(i, index);
		}
	}
	else
	{
		if (index >= mRows)
			throw "index out of index";
		if (mCols == 0 || mRows == 0)
			throw "empty mat";
		result = operator()(index, 0);
		for (unsigned long i = 1; i < mCols; i++)
		{
			if (result < operator()(index, i))
				result = operator()(index, i);
		}
	}
	return result;
}

double EasyMat::min(unsigned long index, EM_DIMENSION dim) const
{
	double result;
	if (dim == EM_BY_COL)
	{
		if (index >= mCols)
			throw "index out of range";
		if (mCols == 0 || mRows == 0)
			throw "empty mat";
		result = operator()(0, index);
		for (unsigned long i = 1; i < mRows; i++)
		{
			if (result > operator()(i, index))
				result = operator()(i, index);
		}
	}
	else
	{
		if (index >= mRows)
			throw "index out of index";
		if (mCols == 0 || mRows == 0)
			throw "empty mat";
		result = operator()(index, 0);
		for (unsigned long i = 1; i < mCols; i++)
		{
			if (result > operator()(index, i))
				result = operator()(index, i);
		}
	}
	return result;
}

double EasyMat::avg(unsigned long index, EM_DIMENSION dim) const
{
	double total = sum(index, dim);
	double result;
	if (dim == EM_BY_COL)
	{
		result = total / mRows;
	}
	else
	{
		result = total / mCols;
	}
	return result;
}

double EasyMat::sum(unsigned long index, EM_DIMENSION dim) const
{
	double result = 0.0;
	if (dim == EM_BY_COL)
	{
		if (index >= mCols)
			throw "index out of range";
		if (mCols == 0 || mRows == 0)
			throw "empty mat";
		for (unsigned long i = 0; i < mRows; i++)
		{
			result += operator()(i, index);
		}
	}
	else
	{
		if (index >= mRows)
			throw "index out of range";
		if (mCols == 0 || mRows == 0)
			throw "empty mat";
		for (unsigned long i = 0; i < mCols; i++)
		{
			result += operator()(index, i);
		}
	}
	return result;
}

EasyVec EasyMat::max(EM_DIMENSION dim) const
{
	if (dim == EM_BY_COL)
	{
		if (mCols == 0 || mRows == 0)
			throw "empty mat";
		EasyVec result = row(0);
		for (unsigned long i = 1; i < mRows; i++)
		{
			for (unsigned long j = 0; j < mCols; j++)
			{
				if (result[j] < operator()(i, j))
				{
					result[j] = operator()(i, j);
				}
			}
		}
		return result;
	}
	else
	{
		if (mCols == 0 || mRows == 0)
			throw "empty mat";
		EasyVec result = col(0);
		for (unsigned long i = 1; i < mCols; i++)
		{
			for (unsigned long j = 0; j < mRows; j++)
			{
				if (result[j] < operator()(j, i))
				{
					result[j] = operator()(j, i);
				}
			}
		}
		return result;
	}
}

EasyVec EasyMat::min(EM_DIMENSION dim) const
{
	if (dim == EM_BY_COL)
	{
		if (mCols == 0 || mRows == 0)
			throw "empty mat";
		EasyVec result = row(0);
		for (unsigned long i = 1; i < mRows; i++)
		{
			for (unsigned long j = 0; j < mCols; j++)
			{
				if (result[j] > operator()(i, j))
				{
					result[j] = operator()(i, j);
				}
			}
		}
		return result;
	}
	else
	{
		if (mCols == 0 || mRows == 0)
			throw "empty mat";
		EasyVec result = col(0);
		for (unsigned long i = 1; i < mCols; i++)
		{
			for (unsigned long j = 0; j < mRows; j++)
			{
				if (result[j] > operator()(j, i))
				{
					result[j] = operator()(j, i);
				}
			}
		}
		return result;
	}
}

EasyVec EasyMat::avg(EM_DIMENSION dim) const
{
	EasyVec result = sum(dim);
	if (dim == EM_BY_COL)
	{
		result /= mRows;
	}
	else
	{
		result /= mCols;
	}
	return result;
}

EasyVec EasyMat::sum(EM_DIMENSION dim) const
{
	if (dim == EM_BY_COL)
	{
		if (mCols == 0 || mRows == 0)
			throw "empty mat";
		EasyVec result = row(0);
		for (unsigned long i = 1; i < mRows; i++)
		{
			for (unsigned long j = 0; j < mCols; j++)
			{
				result[j] += operator()(i, j);
			}
		}
		return result;
	}
	else
	{
		if (mCols == 0 || mRows == 0)
			throw "empty mat";
		EasyVec result = col(0);
		for (unsigned long i = 1; i < mCols; i++)
		{
			for (unsigned long j = 0; j < mRows; j++)
			{
				result[j] += operator()(j, i);
			}
		}
		return result;
	}
}

void EasyMat::abs()
{
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < mCols; j++)
		{
			operator()(i, j) = std::abs(operator()(i, j));
		}
	}
}

void EasyMat::sqrt()
{
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < mCols; j++)
		{
			operator()(i, j) = std::sqrt(operator()(i, j));
		}
	}
}

void EasyMat::pow(double p)
{
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < mCols; j++)
		{
			operator()(i, j) = std::pow(operator()(i, j), p);
		}
	}
}

void EasyMat::reserve(unsigned long newCapacity)
{
	if (newCapacity < mRows * mCols)
		throw "new capacity is too small";
	if (newCapacity == 0)
	{
		if (mData != nullptr)
		{
			delete[] mData;
			mData = nullptr;
			return;
		}
	}
	double* buf = new double[newCapacity];
	if (mCols != 0 && mRows != 0)
	{
		memcpy(buf, mData, mRows * mCols * sizeof(double));
	}
	if (mData != nullptr)
	{
		delete[] mData;
	}
	mData = buf;
	mCapacity = newCapacity;
}

void EasyMat::shrink()
{
	unsigned long dataSize = mRows * mCols;
	reserve(dataSize);
}

void EasyMat::show() const
{
	for (unsigned long i = 0; i < mRows; i++)
	{
		for (unsigned long j = 0; j < mCols; j++)
		{
			std::cout << mData[i * mCols + j] << " ";
		}
		std::cout << std::endl;
	}
}

EasyMat EasyMat::ones(unsigned long rows, unsigned long cols)
{
	EasyMat result(rows, cols);
	if (rows != 0 && cols != 0)
	{
		for (unsigned long i = 0; i < cols; i++)
		{
			result.mData[i] = 1;
		}
		for (unsigned long i = 1; i < rows; i++)
		{
			memcpy(result.mData + i * cols, result.mData, cols * sizeof(double));
		}
	}
	return result;
}

EasyMat EasyMat::zeros(unsigned long rows, unsigned long cols)
{
	return EasyMat(rows,cols);
}

EasyMat EasyMat::eyes(unsigned long rows, unsigned long cols)
{
	EasyMat result(rows, cols);
	unsigned long min = rows > cols ?cols :rows ;
	for (unsigned long i = 0; i < min; i++)
	{
		result(i, i) = 1;
	}
	return result;
}

EasyMat::~EasyMat()
{
	if (mData != nullptr)
	{
		delete[] mData;
		mData = nullptr;
	}
}

void EasyMat::resize(unsigned long newSize)
{
	if (newSize <= mCapacity)
		return;
	if (newSize < INIT_CAPACITY)
		newSize = INIT_CAPACITY;
	double* buf = new double[newSize];
	if (mRows != 0 && mCols != 0 && mData != nullptr)
	{
		memcpy(buf, mData, mRows*mCols * sizeof(double));
	}
	if (mData != nullptr)
	{
		delete[] mData;
	}
	mData = buf;
	mCapacity = newSize;
}

void EasyMat::initMem()
{
	mData = new double[mCapacity];
	memset(mData, 0, sizeof(double) * mCapacity);
}

void EasyMat::moveCol(unsigned long index, long steps)
{
	long tmpCols = mCols + steps;
	if (tmpCols < 0 || tmpCols * mRows > mCapacity)
		throw "move too far...";
	unsigned long totalCount = mCols * mRows;
	for (unsigned long i = 0; i < mRows; i++)
	{
		memcpy(mData + i * tmpCols + index + steps, mData + i * tmpCols + index, 
			(totalCount - (i * mCols + index)) * sizeof(double));
	}
}