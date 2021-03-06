#include "EasyMat.h"
#include "EasyVec.h"
#include "gtest/gtest.h"
#include <random>

using namespace em;

TEST(EasyMatTest, DefaultConstruction)
{
	EasyMat em;
	ASSERT_EQ(0, em.cols());
	ASSERT_EQ(0, em.rows());
}

TEST(EasyMatTest, ColConstruction)
{
	unsigned long col = 10;
	EasyMat em(col);
	ASSERT_EQ(col, em.cols());
	ASSERT_EQ(0, em.rows());
}

TEST(EasyMatTest, RowColConstruction)
{
	unsigned long row = 10;
	unsigned long col = 20;
	EasyMat em(row, col);
	ASSERT_EQ(row, em.rows());
	ASSERT_EQ(col, em.cols());	
}

TEST(EasyMatTest, AppendRowData)
{
	unsigned long col = 800;
	unsigned long row = 30;
	EasyMat em(col);
	EasyVec ev;
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ev.push_back(i * col + j);
		}
		em.append(ev);
		ev.clear();
	}

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(i * col + j, em(i, j));
		}
	}
}

TEST(EasyMatTest, AppendColData)
{
	unsigned long col = 300;
	unsigned long row = 80;
	EasyMat em(row,0);
	EasyVec ev(EM_BY_COL);
	for (unsigned long i = 0; i < col; i++)
	{
		for (unsigned long j = 0; j < row; j++)
		{
			ev.push_back(i * row + j);
		}
		em.append(ev,EM_BY_COL);
		ev.clear();
	}

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(j * row + i, em(i, j));
		}
	}
}

TEST(EasyMatTest, InsertRowData)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em(i, j) = i * col + j;
		}
	}

	EasyVec insertRowBegin;
	EasyVec insertRowMiddle;
	EasyVec insertRowEnd;
	for (unsigned long i = 0; i < col; i++)
	{
		insertRowBegin.push_back(100 + i);
		insertRowMiddle.push_back(200 + i);
		insertRowEnd.push_back(300 + i);
	}

	em.insert(0, insertRowBegin);
	em.insert(200, insertRowMiddle);
	em.insert(em.rows(), insertRowEnd);

	ASSERT_EQ(col, em.cols());
	ASSERT_EQ(row + 3, em.rows());

	unsigned long index = 0;
	for (unsigned long i = 0; i < row; i++)
	{
		if (i == 0)
		{
			for (unsigned long j = 0; j < col; j++)
			{
				ASSERT_EQ(100 + j, em(i, j));
			}
		}
		else if ( i == 200)
		{
			for (unsigned long j = 0; j < col; j++)
			{
				ASSERT_EQ(200 + j, em(i, j));
			}
		}
		else if (i == em.rows() - 1)
		{
			for (unsigned long j = 0; j < col; j++)
			{
				ASSERT_EQ(300 + j, em(i, j));
			}
		}
		else
		{
			for (unsigned long j = 0; j < col; j++)
			{
				if (index * col + j != em(i, j))
					double result = em(i, j);
				ASSERT_EQ(index * col + j, em(i, j));
			}
			++index;
		}
	}
}

TEST(EasyMatTest, InsertColData)
{
	unsigned long col = 300;
	unsigned long row = 80;
	EasyMat em(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em(i, j) = i * col + j;
		}
	}
	EasyVec insertRowBegin(EM_BY_COL);
	EasyVec insertRowMiddle(EM_BY_COL);
	EasyVec insertRowEnd(EM_BY_COL);
	for (unsigned long i = 0; i < row; i++)
	{
		insertRowBegin.push_back(100 + i);
		insertRowMiddle.push_back(200 + i);
		insertRowEnd.push_back(300 + i);
	}

	em.insert(0, insertRowBegin, EM_BY_COL);
	em.insert(15, insertRowMiddle, EM_BY_COL);
	em.insert(em.cols(), insertRowEnd, EM_BY_COL );


	ASSERT_EQ(col + 3, em.cols());
	ASSERT_EQ(row, em.rows());

	unsigned long index = 0;
	for (unsigned long i = 0; i < em.cols(); i++)
	{
		if (i == 0)
		{
			for (unsigned long j = 0; j < row; j++)
			{
				ASSERT_EQ(100 + j, em(j, i));
			}
		}
		else if (i == 15)
		{
			for (unsigned long j = 0; j < row; j++)
			{
				ASSERT_EQ(200 + j, em(j, i));
			}
		}
		else if (i == em.cols() - 1)
		{
			for (unsigned long j = 0; j < row; j++)
			{
				ASSERT_EQ(300 + j, em(j, i));
			}
		}
		else
		{
			for (unsigned long j = 0; j < row; j++)
			{
				if (j * col + index != em(j, i))
					double result = em(j, i);
				ASSERT_EQ(j * col + index, em(j, i));
			}
			++index;
		}
	}
}

TEST(EasyMatTest, RemoveData)
{
	unsigned long col = 80;
	unsigned long row = 800;
	EasyMat em(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em(i,j) = i * col + j;
		}
	}

	em.remove(0);
	em.remove(30);
	em.remove(em.rows() - 1);
	em.remove(0, EM_BY_COL);
	em.remove(40, EM_BY_COL);
	em.remove(em.cols() - 1, EM_BY_COL);

	ASSERT_EQ(row - 3, em.rows());
	ASSERT_EQ(col - 3, em.cols());
	unsigned long indexRow = 0;
	unsigned long indexCol = 0;
	for (unsigned long i = 0; i < em.rows(); i++)
	{
		if (i == 0 || i == 30)
			indexRow++;
		for (unsigned long j = 0; j < em.cols(); j++)
		{
			if (j == 0 || j == 40)
				indexCol++;
			ASSERT_EQ((i + indexRow) * col + j + indexCol, em(i, j));
		}
		indexCol = 0;
	}
}

TEST(EasyMatTest, SetVecTest)
{
	unsigned long col = 800;
	unsigned long row = 80;
	EasyMat em1(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = i * col + j;
		}
	}
	EasyMat em2(em1);

	EasyVec vr1(col);
	EasyVec vr2(col);
	EasyVec vr3(col);

	for (unsigned long i = 0; i < col; i++)
	{
		vr1[i] = 100 + i;
		vr2[i] = 200 + i;
		vr3[i] = 300 + i;
	}
	em1.set(0, vr1);
	em1.set(30, vr2);
	em1.set(row - 1, vr3);
	for (unsigned long i = 0; i < em1.rows(); i++)
	{
		if (i == 0)
		{
			for (unsigned long j = 0; j < em1.cols(); j++)
			{
				ASSERT_EQ(100 + j, em1(i, j));
			}
		}
		else if (i == 30)
		{
			for (unsigned long j = 0; j < em1.cols(); j++)
			{
				ASSERT_EQ(200 + j, em1(i, j));
			}
		}
		else if (i == row - 1)
		{
			for (unsigned long j = 0; j < em1.cols(); j++)
			{
				ASSERT_EQ(300 + j, em1(i, j));
			}
		}
		else
		{
			for (unsigned long j = 0; j < em1.cols(); j++)
			{
				ASSERT_EQ(i * col + j, em1(i, j));
			}
		}
	}

	EasyVec vc1(row,EM_BY_COL);
	EasyVec vc2(row, EM_BY_COL);
	EasyVec vc3(row, EM_BY_COL);

	for (unsigned long i = 0; i < row; i++)
	{
		vc1[i] = 400 + i;
		vc2[i] = 500 + i;
		vc3[i] = 600 + i;
	}
	em2.set(0, vc1, EM_BY_COL);
	em2.set(40, vc2, EM_BY_COL);
	em2.set(col - 1, vc3, EM_BY_COL);
	for (unsigned long i = 0; i < em2.cols(); i++)
	{
		if (i == 0)
		{
			for (unsigned long j = 0; j < em2.rows(); j++)
			{
				ASSERT_EQ(400 + j, em2(j, i));
			}
		}
		else if (i == 40)
		{
			for (unsigned long j = 0; j < em2.rows(); j++)
			{
				ASSERT_EQ(500 + j, em2(j, i));
			}
		}
		else if (i == col - 1)
		{
			for (unsigned long j = 0; j < em2.rows(); j++)
			{
				ASSERT_EQ(600 + j, em2(j, i));
			}
		}
		else
		{
			for (unsigned long j = 0; j < em2.rows(); j++)
			{
				ASSERT_EQ(j * col + i, em2(j, i));
			}
		}
	}
}

TEST(EasyMatTest, SetDataTest)
{
	unsigned long col = 800;
	unsigned long row = 80;
	EasyMat em1(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = i * col + j;
		}
	}
	EasyMat em2(em1);

	em1.set(0, 100);
	em1.set(30, 200);
	em1.set(row - 1, 300);
	for (unsigned long i = 0; i < em1.rows(); i++)
	{
		if (i == 0)
		{
			for (unsigned long j = 0; j < em1.cols(); j++)
			{
				ASSERT_EQ(100, em1(i, j));
			}
		}
		else if (i == 30)
		{
			for (unsigned long j = 0; j < em1.cols(); j++)
			{
				ASSERT_EQ(200, em1(i, j));
			}
		}
		else if (i == row - 1)
		{
			for (unsigned long j = 0; j < em1.cols(); j++)
			{
				ASSERT_EQ(300, em1(i, j));
			}
		}
		else
		{
			for (unsigned long j = 0; j < em1.cols(); j++)
			{
				ASSERT_EQ(i * col + j, em1(i, j));
			}
		}
	}

	em2.set(0, 400, EM_BY_COL);
	em2.set(40, 500, EM_BY_COL);
	em2.set(col - 1, 600, EM_BY_COL);
	for (unsigned long i = 0; i < em2.cols(); i++)
	{
		if (i == 0)
		{
			for (unsigned long j = 0; j < em2.rows(); j++)
			{
				ASSERT_EQ(400, em2(j, i));
			}
		}
		else if (i == 40)
		{
			for (unsigned long j = 0; j < em2.rows(); j++)
			{
				ASSERT_EQ(500, em2(j, i));
			}
		}
		else if (i == col - 1)
		{
			for (unsigned long j = 0; j < em2.rows(); j++)
			{
				ASSERT_EQ(600, em2(j, i));
			}
		}
		else
		{
			for (unsigned long j = 0; j < em2.rows(); j++)
			{
				ASSERT_EQ(j * col + i, em2(j, i));
			}
		}
	}
}

TEST(EasyMatTest, OperatorAdd)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = i * col + j;
		}
	}

	EasyMat em2(em1);
	EasyMat result(row, col);
	result = em1 + em2;
	
	for (unsigned long i = 0; i < row;i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ((i * col + j) * 2, result(i, j));
			ASSERT_EQ(i * col + j, em1(i, j));
			ASSERT_EQ(i * col + j, em2(i, j));
		}
	}

	double tmp = 1.0;
	result = em1 + tmp;
	for (unsigned long i = 0; i < row;i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(i * col + j + tmp, result(i, j));
			ASSERT_EQ(i * col + j, em1(i, j));
		}
	}

	em1 += em2;
	for (unsigned long i = 0; i < row;i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ((i * col + j) * 2, em1(i, j));
		}
	}

	em2 += tmp;
	for (unsigned long i = 0; i < row;i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(i * col + j + tmp, em2(i, j));
		}
	}
}

TEST(EasyMatTest, OperatorSub)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = i * col + j;
		}
	}

	EasyMat em2(em1);
	EasyMat result(row, col);
	result = em1 - em2;
	for (unsigned long i = 0; i < row;i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(0, result(i, j));
			ASSERT_EQ(i * col + j, em1(i, j));
			ASSERT_EQ(i * col + j, em2(i, j));
		}
	}

	double tmp = 1.0;
	result = em1 - tmp;
	for (unsigned long i = 0; i < row;i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(i * col + j - tmp, result(i, j));
			ASSERT_EQ(i * col + j, em1(i, j));
		}
	}

	em1 -= em2;
	for (unsigned long i = 0; i < row;i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(0, em1(i, j));
		}
	}

	em2 -= tmp;
	for (unsigned long i = 0; i < row;i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(i * col + j - tmp, em2(i, j));
		}
	}
}

TEST(EasyMatTest, OperatorMulti)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = i * col + j;
		}
	}

	EasyMat em2 = em1.transpose();
	EasyMat result(row, row);
	result = em1 * em2;
	for (unsigned long i = 0; i < row;i++)
	{
		for (unsigned long j = 0; j < row; j++)
		{
			double element = 0;
			for (unsigned long k = 0; k < col; k++)
			{
				element += em1(i, k) * em2(k, j);
			}
			ASSERT_EQ(element, result(i, j));
		}
	}

	double tmp = 2.0;
	result = em1 * tmp;
	for (unsigned long i = 0; i < row;i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ((i * col + j) * tmp, result(i, j));
			ASSERT_EQ(i * col + j, em1(i, j));
		}
	}

	em1 *= EasyMat::ones(col,col);
	for (unsigned long i = 0; i < row;i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			double element = 0; 
			for (unsigned long k = 0; k < col; k++)
			{
				element += i * col + k;
			}
			ASSERT_EQ(element, em1(i, j));
		}
	}

	em2 *= tmp;
	for (unsigned long i = 0; i < row;i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ((i * col + j) * tmp, em2(j, i));
		}
	}
}

TEST(EasyMatTest, OperatorDivision)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em(i, j) = i * col + j;
		}
	}

	double tmp = 2.0;
	EasyMat result(row, col);
	result = em / tmp;
	for (unsigned long i = 0; i < em.rows(); i++)
	{
		for (unsigned long j = 0; j < em.cols(); j++)
		{
			ASSERT_EQ((i*col + j) / tmp, result(i, j));
			ASSERT_EQ(i*col + j, em(i, j));
		}
	}

	em /= tmp;
	for (unsigned long i = 0; i < em.rows(); i++)
	{
		for (unsigned long j = 0; j < em.cols(); j++)
		{
			ASSERT_EQ((i*col + j) / tmp, em(i, j));
		}
	}
}

TEST(EasyMatTest, OperatorEqualandNotEqual)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = i * col + j;
		}
	}

	EasyMat em2(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em2(i, j) = i * col + j;
		}
	}

	ASSERT_TRUE(em1 == em2);
	em1(0, 0) = 1;
	ASSERT_TRUE(em1 != em2);
}

TEST(EasyMatTest, TransposeTest)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em(i, j) = i * col + j;
		}
	}

	EasyMat t = em.transpose();
	ASSERT_EQ(row, t.cols());
	ASSERT_EQ(col, t.rows());
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(i * col + j, em(i, j));
			ASSERT_EQ(i * col + j, t(j, i));
		}
	}
}

TEST(EasyMatTest, SwapTest)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = i * col + j;
		}
	}
	EasyMat em2(em1);

	em1.swap(2,4);
	for (unsigned long i = 0; i < row; i++)
	{
		unsigned long index = i;
		if (i == 2)
			index = 4;
		else if (i == 4)
			index = 2;
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(index * col + j, em1(i, j));
		}
	}

	em2.swap(0, 1,EM_BY_COL);
	for (unsigned long i = 0; i < col; i++)
	{
		unsigned long index = i;
		if (i == 0)
			index = 1;
		else if (i == 1)
			index = 0;
		for (unsigned long j = 0; j < row; j++)
		{
			ASSERT_EQ(j * col + index, em2(j, i));
		}
	}
}

TEST(EasyMatTest, SortTest)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	std::default_random_engine e((unsigned int)time(nullptr));
	std::uniform_real_distribution<double> u(0, 10);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = u(e);
		}
	}
	EasyMat em2(em1);
	em2.sortByRowAsc(0);

	double sum1 = 0;
	double sum2 = 0;
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			sum1 += em1(i, j);
			sum2 += em2(i, j);
		}
	}
	ASSERT_GE(EasyMat::precision(), abs(sum1 - sum2));

	for (unsigned long i = 1; i < row; i++)
	{
		ASSERT_GE(em2(i, 0), em2(i - 1, 0));
	}

	em2.sortByRowDesc(1);
	sum1 = 0;
	sum2 = 0;
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			sum1 += em1(i, j);
			sum2 += em2(i, j);
		}
	}
	ASSERT_GE(EasyMat::precision(), abs(sum1 - sum2));

	for (unsigned long i = 1; i < row; i++)
	{
		ASSERT_LE(em2(i, 1), em2(i - 1, 1));
	}
}

TEST(EasyMatTest, RowAndColTest)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em(i, j) = i * col + j;
		}
	}

	EasyVec ev1 = em.row(1);
	ASSERT_EQ(ev1.getDimension(), EM_BY_ROW);
	for (unsigned long i = 0; i < col; i++)
	{
		ASSERT_EQ(ev1[i], em(1, i));
	}

	EasyVec ev2 = em.col(2);
	ASSERT_EQ(ev2.getDimension(), EM_BY_COL);
	for (unsigned long i = 0; i < row; i++)
	{
		ASSERT_EQ(ev2[i], em(i, 2));
	}

	EasyVec tmp(2);
	tmp[0] = 1;
	tmp[1] = 4;
	EasyMat em1 = em.row(tmp);
	ASSERT_EQ(2, em1.rows());
	ASSERT_EQ(col, em1.cols());
	for (unsigned long i = 0; i < tmp.size(); i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(em1(i, j), em((unsigned long)tmp[i], j));
		}
	}

	tmp[0] = 0;
	tmp[1] = 2;
	EasyMat em2 = em.col(tmp);
	ASSERT_EQ(row, em2.rows());
	ASSERT_EQ(2, em2.cols());
	for (unsigned long i = 0; i < tmp.size(); i++)
	{
		for (unsigned long j = 0; j < row; j++)
		{
			ASSERT_EQ(em2(j, i), em(j, (unsigned long)tmp[i]));
		}
	}
}

TEST(EasyMatTest, MaxMinSumAvgTest)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em(row, col);
	std::default_random_engine e((unsigned int)time(nullptr));
	std::uniform_real_distribution<double> u(0, 10);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em(i, j) = u(e);
		}
	}

	unsigned long rowIndex = 4;
	double rowMax = em.max(rowIndex, EM_BY_ROW);
	double rowMin = em.min(rowIndex, EM_BY_ROW);
	double rowSum = em.sum(rowIndex, EM_BY_ROW);
	double rowAvg = em.avg(rowIndex, EM_BY_ROW);

	unsigned long colIndex = 1;
	double colMax = em.max(colIndex);
	double colMin = em.min(colIndex);
	double colSum = em.sum(colIndex);
	double colAvg = em.avg(colIndex);

	double rowRealMax = em(rowIndex, 0);
	double rowRealMin = em(rowIndex, 0);
	double rowRealSum = 0.0;
	double rowRealAvg = 0.0;
	double colRealMax = em(0, colIndex);
	double colRealMin = em(0, colIndex);
	double colRealSum = 0.0;
	double colRealAvg = 0.0;

	for (unsigned long i = 0; i < col; i++)
	{
		if (em(rowIndex, i) > rowRealMax)
			rowRealMax = em(rowIndex, i);
		if (em(rowIndex, i) < rowRealMin)
			rowRealMin = em(rowIndex, i);
		rowRealSum += em(rowIndex, i);
	}
	rowRealAvg = rowRealSum / col;

	for (unsigned long i = 0; i < row; i++)
	{
		if (em(i, colIndex) > colRealMax)
			colRealMax = em(i, colIndex);
		if (em(i, colIndex) < colRealMin)
			colRealMin = em(i, colIndex);
		colRealSum += em(i, colIndex);
	}
	colRealAvg = colRealSum / row;

	ASSERT_EQ(rowRealMax, rowMax);
	ASSERT_EQ(rowRealMin, rowMin);
	ASSERT_EQ(rowRealSum, rowSum);
	ASSERT_EQ(rowRealAvg, rowAvg);
	ASSERT_EQ(colRealMax, colMax);
	ASSERT_EQ(colRealMin, colMin);
	ASSERT_GE(EasyMat::precision(), std::abs(colRealSum - colSum));	
	ASSERT_GE(EasyMat::precision(), std::abs(colRealAvg - colAvg));
}

TEST(EasyMatTest, MaxMinSumAvgVecTest)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em(row, col);
	std::default_random_engine e((unsigned int)time(nullptr));
	std::uniform_real_distribution<double> u(0, 10);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em(i, j) = u(e);
		}
	}

	EasyVec rowMaxVec = em.max(EM_BY_ROW);
	EasyVec rowMinVec = em.min(EM_BY_ROW);
	EasyVec rowSumVec = em.sum(EM_BY_ROW);
	EasyVec rowAvgVec = em.avg(EM_BY_ROW);
	EasyVec colMaxVec;
	colMaxVec = em.max();
	EasyVec colMinVec = em.min();
	EasyVec colSumVec = em.sum();
	EasyVec colAvgVec = em.avg();

	EasyVec rowRealMaxVec(row, EM_BY_COL);
	EasyVec rowRealMinVec(row, EM_BY_COL);
	EasyVec rowRealSumVec(row, EM_BY_COL);
	EasyVec rowRealAvgVec(row, EM_BY_COL);
	EasyVec colRealMaxVec(col);
	EasyVec colRealMinVec(col);
	EasyVec colRealSumVec(col);
	EasyVec colRealAvgVec(col);

	for (unsigned long i = 0;i < row; i++)
	{
		rowRealMaxVec[i] = em.max(i, EM_BY_ROW);
		rowRealMinVec[i] = em.min(i, EM_BY_ROW);
		rowRealSumVec[i] = em.sum(i, EM_BY_ROW);
		rowRealAvgVec[i] = em.avg(i, EM_BY_ROW);
	}

	for (unsigned long i = 0; i < col; i++)
	{
		colRealMaxVec[i] = em.max(i);
		colRealMinVec[i] = em.min(i);
		colRealSumVec[i] = em.sum(i);
		colRealAvgVec[i] = em.avg(i);
	}

	if (rowRealMaxVec != rowMaxVec)
	{
		rowRealMaxVec.show();
		rowMaxVec.show();
	}
	ASSERT_EQ(rowRealMaxVec, rowMaxVec);
	ASSERT_EQ(rowRealMinVec, rowMinVec);
	ASSERT_EQ(rowRealSumVec, rowSumVec);
	ASSERT_EQ(rowRealAvgVec, rowAvgVec);
	ASSERT_EQ(colRealMaxVec, colMaxVec);
	ASSERT_EQ(colRealMinVec, colMinVec);
	ASSERT_EQ(colRealSumVec, colSumVec);
	ASSERT_EQ(colRealAvgVec, colAvgVec);
}

TEST(EasyMatTest, OnesTest)
{
	unsigned long row = 800;
	unsigned long col = 30;
	EasyMat em = EasyMat::ones(row, col);

	ASSERT_EQ(row, em.rows());
	ASSERT_EQ(col, em.cols());
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(1, em(i, j));
		}
	}
}

TEST(EasyMatTest, ZerosTest)
{
	unsigned long row = 800;
	unsigned long col = 30;
	EasyMat em = EasyMat::zeros(row, col);

	ASSERT_EQ(row, em.rows());
	ASSERT_EQ(col, em.cols());
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(0, em(i, j));
		}
	}
}

TEST(EasyMatTest, eyesTest)
{
	unsigned long row = 800;
	unsigned long col = 30;
	EasyMat em = EasyMat::eyes(row, col);

	ASSERT_EQ(row, em.rows());
	ASSERT_EQ(col, em.cols());
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			if (i != j)
			{
				ASSERT_EQ(0, em(i, j));
			}
			else
			{
				ASSERT_EQ(1, em(i, j));
			}
		}
	}
}

TEST(EasyMatTest, ReserveTest)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em(i, j) = i * col + j;
		}
	}

	unsigned long capacity = em.capacity();
	ASSERT_ANY_THROW(em.reserve(0));
	ASSERT_ANY_THROW(em.reserve(capacity - 1));
	em.reserve(capacity +1);
	ASSERT_EQ(capacity + 1, em.capacity());
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(i * col + j, em(i, j));
		}
	}
}

TEST(EasyMatTest, ShrinkTest)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em(i, j) = i * col + j;
		}
	}
	unsigned long capacity = em.capacity();
	em.shrink();
	ASSERT_EQ(em.rows() * em.cols(), em.capacity());
	ASSERT_GE(capacity, em.capacity());
}

TEST(EasyMatTest, AddVector)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = i * col + j;
		}
	}
	EasyMat em2(em1);

	EasyVec vec1;
	for (unsigned long i = 0; i < col;i++)
	{
		vec1.push_back(i);
	}
	em1.add(vec1);

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(i * col + j + j, em1(i, j));
		}
	}

	EasyVec vec2(EM_BY_COL);
	for (unsigned long i = 0; i < row; i++)
	{
		vec2.push_back(i);
	}
	em2.add(vec2, EM_BY_COL);

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(i * col + j + i, em2(i, j));
		}
	}
}

TEST(EasyMatTest, SubVector)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = i * col + j;
		}
	}
	EasyMat em2(em1);

	EasyVec vec1;
	for (unsigned long i = 0; i < col;i++)
	{
		vec1.push_back(i);
	}
	em1.sub(vec1);

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(i * col + j - j, em1(i, j));
		}
	}

	EasyVec vec2(EM_BY_COL);
	for (unsigned long i = 0; i < row; i++)
	{
		vec2.push_back(i);
	}
	em2.sub(vec2, EM_BY_COL);

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ(i * col + j - i, em2(i, j));
		}
	}
}

TEST(EasyMatTest, MulVector)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = i * col + j;
		}
	}
	EasyMat em2(em1);

	EasyVec vec1;
	for (unsigned long i = 0; i < col;i++)
	{
		vec1.push_back(i);
	}
	em1.mul(vec1);

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ((i * col + j) * j, em1(i, j));
		}
	}

	EasyVec vec2(EM_BY_COL);
	for (unsigned long i = 0; i < row; i++)
	{
		vec2.push_back(i);
	}
	em2.mul(vec2, EM_BY_COL);

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ((i * col + j) * i, em2(i, j));
		}
	}
}

TEST(EasyMatTest, DivVector)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = i * col + j;
		}
	}
	EasyMat em2(em1);

	EasyVec vec1;
	for (unsigned long i = 0; i < col;i++)
	{
		vec1.push_back(i+1);
	}
	em1.div(vec1);

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ((double)(i * col + j) / (j + 1), em1(i, j));
		}
	}

	EasyVec vec2(EM_BY_COL);
	for (unsigned long i = 0; i < row; i++)
	{
		vec2.push_back(i + 1);
	}
	em2.div(vec2, EM_BY_COL);

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_EQ((double)(i * col + j) / (i + 1), em2(i, j));
		}
	}
}

TEST(EasyMatTest, AbsTest)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	std::default_random_engine e((unsigned int)time(nullptr));
	std::uniform_real_distribution<double> u(0, 10);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = u(e) - 5;
		}
	}
	EasyMat em2(em1);
	em1.abs();

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			if (em2(i, j) < 0)
			{
				ASSERT_EQ(em2(i, j), -em1(i, j));
			}
			else
			{
				ASSERT_EQ(em2(i, j), em1(i, j));
			}
		}
	}
}

TEST(EasyMatTest, SqrtTest)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	std::default_random_engine e((unsigned int)time(nullptr));
	std::uniform_real_distribution<double> u(0, 10);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = u(e);
		}
	}
	EasyMat em2(em1);
	em1.sqrt();

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_GE(EasyMat::precision(), std::abs(em1(i, j) * em1(i, j) - em2(i, j)));
		}
	}
}

TEST(EasyMatTest, PowTest)
{
	unsigned long col = 30;
	unsigned long row = 800;
	EasyMat em1(row, col);
	std::default_random_engine e((unsigned int)time(nullptr));
	std::uniform_real_distribution<double> u(0, 10);
	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			em1(i, j) = u(e);
		}
	}
	EasyMat em2(em1);
	em1.pow(2);

	for (unsigned long i = 0; i < row; i++)
	{
		for (unsigned long j = 0; j < col; j++)
		{
			ASSERT_GE(EasyMat::precision(), std::abs(em1(i, j) - em2(i, j) * em2(i, j)));
		}
	}
}