#ifndef PARALUTIONINTERFACE_H_
#define PARALUTIONINTERFACE_H_

#include <string>

#include "paralution.hpp"

enum class PRECOND {ILU_SIMPLE, ILU_SERIOUS, ILUT, ILU_GMRES};

class ParSolver
{
	enum class RETURN_TYPE { NO_CRITERIA, ABS_CRITERION, REL_CRITERION, DIV_CRITERIA, MAX_ITER };
public:
	typedef paralution::LocalMatrix<double> Matrix;
	typedef paralution::LocalVector<double> Vector;
protected:
	Vector x, Rhs;
	Matrix Mat;
	paralution::BiCGStab<Matrix,Vector,double> bicgstab;
	void SolveBiCGStab();
	void SolveBiCGStab_Simple();
	paralution::GMRES<Matrix,Vector,double> gmres;
	void SolveGMRES();
	paralution::ILU<Matrix,Vector,double> p;
	paralution::ILUT<Matrix, Vector, double> p_ilut;

	bool isAssembled;
	bool isPrecondBuilt;
	bool isTheSameMatrix;
	bool isCleared;
	int matSize;
	RETURN_TYPE status;

	inline void writeSystem()
	{
		Mat.WriteFileMTX("snaps/mat.mtx");
		Rhs.WriteFileASCII("snaps/rhs.dat");
		x.WriteFileASCII("snaps/x.dat");
	};

	double initRes, finalRes;
	int iterNum;
	const std::string resHistoryFile;
	void getResiduals();
public:
	void Init(const int vecSize, const double relTol, const double dropTol);
	void Assemble(const int* ind_i, const int* ind_j, const double* a, const int counter, const int* ind_rhs, const double* rhs);
	void Solve();
	void Solve(const PRECOND key);
	void SetSameMatrix();
	void Clear();

	const Vector& getSolution() { return x; };
	void getInvert(const int* ind_i, const int* ind_j, const double* a, const int counter, int*& offset, int*& col, double*& dmat)
	{
		Mat.Zeros();
		Mat.Assemble(ind_i, ind_j, a, counter, "A", matSize, matSize);
		Mat.Invert();
		/*Mat.LeaveDataPtrDENSE(&dmat);
		double tmp;
		for (int i = 0; i < matSize; i++)
			for (int j = i + 1; j < matSize; j++)
			{
				tmp = dmat[i * matSize + j];
				dmat[i * matSize + j] = dmat[j * matSize + i];
				dmat[j * matSize + i] = tmp;
			}*/
		Mat.LeaveDataPtrCSR(&offset, &col, &dmat);
	};

	ParSolver();
	~ParSolver();
};

#endif /* PARALUTIONINTERFACE_H_ */