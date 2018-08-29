#pragma once
#include "EasyMat.h"
namespace em
{
	class EasyVec : public EasyMat
	{
	public:
		explicit EasyVec(EM_DIMENSION dim = EM_BY_ROW);
		explicit EasyVec(unsigned long size, EM_DIMENSION dim = EM_BY_ROW);
		EasyVec(const EasyVec& rhs);
		EasyVec(EasyVec&& rhs);
		~EasyVec();
		void push_back(double value);
		void pop_back();
		void insert(unsigned long index, double value);
		void erase(unsigned long index);
		void clear();
		unsigned long size() const;
		double max() const;
		double min() const;
		double avg() const;
		double sum() const;
		EasyVec transpose();
		double& operator[](unsigned long index);
		const double& operator[](unsigned long index) const;
		EM_DIMENSION getDimension() const { return mDim; }
	};
}

