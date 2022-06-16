#pragma once
#include <functional>
#include <vector>
#include <array>
#include "Vector.h"
#include "Matrix.h"
#include "problem_interface.h"
#include "stronginc3_problem.h"

typedef std::function<double(QuadProg::Vector<double>)> PFunc;

class TaskPrep
{
private:
	/// <summary>
	/// Указатель на целевую функцию
	/// </summary>
	PFunc _problem;
	/// <summary>
	/// Вектор указателей на ограничения, заданные равенством
	/// </summary>
	std::vector<PFunc> _equalConstraint;
	/// <summary>
	/// Вектор указателей на ограничения, заданные неравенством
	/// </summary>
	std::vector<PFunc> _inequalConstraint;

	const std::array<std::vector<int>, 4> coeff =
	{ { {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} } };
	const std::array<std::vector<int>, 4> coeff2 =
	{ { {1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} } };
	const std::array<int, 4> dd = { 2, 12, 60, 840 };

	QuadProg::Vector<double> _startPoint;

public:
	TaskPrep(QuadProg::Vector<double> startPoint, PFunc problem) : _problem(problem), _startPoint(startPoint){

	}

	void SetEqConstr(std::vector<PFunc> eqContr);
	void SetIqConstr(std::vector<PFunc> iqContr);

private:
	void FiniteGradient(const QuadProg::Vector<double>& x, QuadProg::Vector<double>& grad, PFunc func, int accuracy);

	void FiniteHessian(const QuadProg::Vector<double>& x, QuadProg::Matrix<double>& hessian, int accuracy);

	/// <summary>
	/// Расчёт градиента целевой функции
	/// </summary>
	/// <param name="x">Точка расчёта</param>
	/// <param name="c">Полученный градиент</param>
	void Gradient(const QuadProg::Vector<double> x, QuadProg::Vector<double>& grad);

	//The Hessian Matrix of the objective function
	void Hessian(const QuadProg::Vector<double>& x, QuadProg::Matrix<double>& Q);

	//The gradient of equation constraint function
	void GradientCE(const QuadProg::Vector<double>& x, QuadProg::Matrix<double>& CE);

	//The constant term of equation constraint function
	void CE0(const QuadProg::Vector<double>& x, QuadProg::Vector<double>& ce0);

	//The gradient of inequation constraint function
	void GradientCI(const QuadProg::Vector<double>& x, QuadProg::Matrix<double>& CI);

	//The constant term of inequation constraint function
	void CI0(const QuadProg::Vector<double>& x, QuadProg::Vector<double>& ci0);

	// The value of the Quadratic function
	double Quadfun(QuadProg::Matrix<double> Q, QuadProg::Vector<double> c, QuadProg::Vector<double> d);


public:
	void Solve();



};

