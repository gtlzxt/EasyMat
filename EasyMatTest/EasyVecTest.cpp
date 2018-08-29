#include "EasyMat.h"
#include "EasyVec.h"
#include "gtest/gtest.h"

using namespace em;

TEST(EasyVecTest, push_backTest)
{
	unsigned long size = 100;
	EasyVec vec1;
	for (unsigned long i = 0; i < size; i++)
	{
		vec1.push_back(i);
	}

	ASSERT_EQ(size, vec1.size());
	ASSERT_EQ(EM_BY_ROW, vec1.getDimension());
	for (unsigned long i = 0; i < size; i++)
	{
		ASSERT_EQ(i, vec1[i]);
	}

	EasyVec vec2(EM_BY_COL);
	for (unsigned long i = 0; i < size; i++)
	{
		vec2.push_back(i);
	}

	ASSERT_EQ(size, vec2.size());
	ASSERT_EQ(EM_BY_COL, vec2.getDimension());
	for (unsigned long i = 0; i < size; i++)
	{
		ASSERT_EQ(i, vec2[i]);
	}
}