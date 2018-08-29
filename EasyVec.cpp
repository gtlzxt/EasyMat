#include "EasyVec.h"
#include <iostream>
using namespace em;

EasyVec::EasyVec(EM_DIMENSION dim)
{
	mDim = dim;
	if (dim == EM_BY_ROW)
		mRows = 1;
	else
		mCols = 1;
}

EasyVec::EasyVec(unsigned long size, EM_DIMENSION dim):EasyMat(1, size)
{
	mDim = dim;
	if (dim == EM_BY_COL)
	{
		mRows = size;
		mCols = 1;
	}
}

EasyVec::EasyVec(const EasyVec& rhs):EasyMat(rhs)
{

}
EasyVec::EasyVec(EasyVec&& rhs)
{
	*this = std::move(rhs);
}

EasyVec::EasyVec(const EasyMat& rhs):EasyMat(rhs)
{

}

EasyVec::EasyVec(EasyMat&& rhs)
{
	*this = std::move(EasyVec(rhs));
}

const EasyVec& EasyVec::operator=(const EasyVec& rhs)
{
	EasyMat::operator=(rhs);
	return *this;
}

const EasyVec& EasyVec::operator=(EasyVec&& rhs)
{
	EasyMat::operator=(std::move(rhs));
	return *this;
}

EasyVec::~EasyVec()
{
}

void EasyVec::push_back(double value)
{
	if (mDim == EM_BY_ROW)
	{
		++mCols;
	}
	else
	{
		++mRows;
	}
	unsigned long size = this->size();
	if (size > mCapacity)
		resize(size * 2);
	mData[size - 1] = value;
}

void EasyVec::pop_back()
{
	unsigned long size = this->size();
	if (size == 0)
		throw "empty vector, can not pop";
	if (mDim == EM_BY_ROW)
	{
		--mCols;
	}
	else
	{
		--mRows;
	}
}

void EasyVec::insert(unsigned long index, double value)
{
	unsigned long size = this->size();
	if (index > size)
		throw "index is bigger than the size of vector";
	if (mDim == EM_BY_ROW)
	{
		++mCols;
	}
	else
	{
		++mRows;
	}
	++size;
	if (size > mCapacity)
		resize(size * 2);
	memcpy(mData + index + 1, mData + index, (size - index) * sizeof(double));
	mData[index] = value;
}

void EasyVec::erase(unsigned long index)
{
	unsigned long size = this->size();
	if (index > size)
		throw "index is bigger than the size of vector";
	if (mDim == EM_BY_ROW)
	{
		--mCols;
	}
	else
	{
		--mRows;
	}
	if(index != size)
		memcpy(mData + index, mData + index + 1, (size - index) * sizeof(double));
}

void EasyVec::clear()
{
	mCols = 0;
	mRows = 0;
}

double EasyVec::max() const
{
	unsigned long size = this->size();
	if (size == 0)
		throw "the vector is empty";
	double result = mData[0];
	for (unsigned long i = 1; i < size; i++)
	{
		if (mData[i] > result)
			result = mData[i];
	}
	return result;
}

double EasyVec::min() const
{
	unsigned long size = this->size();
	if (size == 0)
		throw "the vector is empty";
	double result = mData[0];
	for (unsigned long i = 1; i < size; i++)
	{
		if (mData[i] < result)
			result = mData[i];
	}
	return result;
}

double EasyVec::avg() const
{
	return sum() / size();
}

double EasyVec::sum() const
{
	unsigned long size = this->size();
	if (size == 0)
		throw "the vector is empty";
	double result = 0;
	for (unsigned long i = 0; i < size; i++)
	{
		result += mData[i];
	}
	return result;
}

unsigned long EasyVec::size() const
{
	if (mDim == EM_BY_ROW)
	{
		return mCols;
	}
	else
	{
		return mRows;
	}
}

EasyVec EasyVec::transpose()
{
	EasyVec result(*this);
	if (mDim == EM_BY_ROW)
	{
		result.mDim = EM_BY_COL;
	}
	else
	{
		result.mDim = EM_BY_ROW;
	}
	result.mRows = mCols;
	result.mCols = mRows;
	return result;
}

double & EasyVec::operator[](unsigned long index)
{
	unsigned long size = this->size();
	if (index >= size)
		throw "index outof range";
	return mData[index];
}

const double & EasyVec::operator[](unsigned long index) const
{
	unsigned long size = this->size();
	if (index >= size)
		throw "index outof range";
	return mData[index];
}
