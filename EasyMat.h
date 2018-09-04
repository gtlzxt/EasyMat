#pragma once
/*
*To make it simple and avaiable to most kind of arithmetical opertion,
*the raw data is of kind double,
*
*/
namespace std
{
	class thread;
}
namespace em
{
	constexpr auto THREAD_NUMBER = 4;
	class EasyVec;
	enum EM_DIMENSION
	{
		EM_BY_COL,
		EM_BY_ROW
	};
	class EasyMat
	{
	public:
		EasyMat();
		explicit EasyMat(unsigned long cols);
		EasyMat(unsigned long rows, unsigned long cols);
		EasyMat(const EasyMat& rhs);
		EasyMat(EasyMat&& rhs);
		virtual ~EasyMat();
		const EasyMat& operator=(const EasyMat& rhs);
		const EasyMat& operator=(EasyMat&& rhs);
		EasyMat operator+(const EasyMat& rhs) const;
		EasyMat operator+(double rhs) const;
		EasyMat operator-(const EasyMat& rhs) const;
		EasyMat operator-(double rhs) const;
		EasyMat operator*(const EasyMat& rhs) const;
		EasyMat operator*(double rhs) const;
		EasyMat operator/(double rhs) const;
		EasyMat& operator+=(const EasyMat& rhs);
		EasyMat& operator+=(double rhs);
		EasyMat& operator-=(const EasyMat& rhs);
		EasyMat& operator-=(double rhs);
		EasyMat& operator*=(const EasyMat& rhs);
		EasyMat& operator*=(double rhs);
		EasyMat& operator/=(double rhs);
		bool operator==(const EasyMat& rhs) const;
		bool operator!=(const EasyMat& rhs) const;
		EasyMat transpose();
		double& operator()(unsigned long i, unsigned long j);
		const double& operator()(unsigned long i, unsigned long j) const;
		void append(const EasyVec& vec, EM_DIMENSION dim = EM_BY_ROW);
		void insert(unsigned long index,const EasyVec& vec, EM_DIMENSION dim = EM_BY_ROW);
		void remove(unsigned long index, EM_DIMENSION dim = EM_BY_ROW);
		void set(unsigned long index, const EasyVec& vec, EM_DIMENSION dim = EM_BY_ROW);
		void set(unsigned long index, double value, EM_DIMENSION dim = EM_BY_ROW);
		void add(const EasyVec& vec, EM_DIMENSION dim = EM_BY_ROW);
		void sub(const EasyVec& vec, EM_DIMENSION dim = EM_BY_ROW);
		void mul(const EasyVec& vec, EM_DIMENSION dim = EM_BY_ROW);
		void div(const EasyVec& vec, EM_DIMENSION dim = EM_BY_ROW);
		void swap(unsigned long index1, unsigned long index2, EM_DIMENSION dim = EM_BY_ROW);
		void sortByRowAsc(unsigned long index);
		void sortByRowDesc(unsigned long index);
		EasyMat col(const EasyVec& vec) const;
		EasyMat row(const EasyVec& vec) const;
		EasyVec col(unsigned long index) const;
		EasyVec row(unsigned long index) const;
		unsigned long cols() const{ return mCols; }
		unsigned long rows() const{ return mRows; }
		unsigned long capacity() const { return mCapacity; }
		static double precision() { return mPrecision; }
		double max(unsigned long index, EM_DIMENSION dim = EM_BY_COL) const;
		double min(unsigned long index, EM_DIMENSION dim = EM_BY_COL) const;
		double avg(unsigned long index, EM_DIMENSION dim = EM_BY_COL) const;
		double sum(unsigned long index, EM_DIMENSION dim = EM_BY_COL) const;
		EasyVec max(EM_DIMENSION dim = EM_BY_COL) const;
		EasyVec min(EM_DIMENSION dim = EM_BY_COL) const;
		EasyVec avg(EM_DIMENSION dim = EM_BY_COL) const;
		EasyVec sum(EM_DIMENSION dim = EM_BY_COL) const;
		void abs();
		void sqrt();
		void pow(double p);
		void reserve(unsigned long newCapacity);
		void shrink();
		void show() const;
		static EasyMat ones(unsigned long rows, unsigned long cols);
		static EasyMat zeros(unsigned long rows, unsigned long cols);
		static EasyMat eyes(unsigned long rows, unsigned long cols);
		
	protected:
		void resize(unsigned long newSize);
		unsigned long mCols;
		unsigned long mRows;
		double* mData;
		unsigned long mCapacity;
		EM_DIMENSION mDim;
		static const double mPrecision;
	private:
		class Thread_Guard
		{
		public:
			explicit Thread_Guard(std::thread& t);
			~Thread_Guard();
			Thread_Guard(const Thread_Guard&) = delete;
			Thread_Guard& operator=(const Thread_Guard&) = delete;
		private:
			std::thread& mThread;
		};
		void initMem();
		//process data with multi thread
		template<typename fun>
		void processData(fun f) const
		{
			if (mRows > NEED_MULTITHREAD && THREAD_NUMBER > 1)
			{
				std::thread t[THREAD_NUMBER];
				unsigned long batch = mRows / THREAD_NUMBER;
				for (int i = 0; i < THREAD_NUMBER - 1;i++)
				{
					t[i] = std::thread(f, i*batch, (i + 1)*batch);
				}
				t[THREAD_NUMBER - 1] = std::thread(f, (THREAD_NUMBER - 1) * batch, mRows);

				for (int i = 0; i < THREAD_NUMBER;i++)
				{
					t[i].join();
				}
			}
			else
			{
				f(0, mRows);
			}
		}
		/*
		move (index to mCols) cols data, used for add, insert or remove cols
		*/
		void moveCol(unsigned long index, long steps);
		const static unsigned long INIT_CAPACITY = 64;
		const static unsigned long NEED_MULTITHREAD = 64;
	};
}

