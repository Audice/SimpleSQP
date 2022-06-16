#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#include "TaskPrep.h"
#include "stronginc3_problem.h"

using namespace QuadProg;

//double FuncValue(const Vector<double>& x) {
//	double f = pow(x[1] - x[0] * x[0], 2) + pow(1 - x[0], 2);
//	return f;
//}
//
//double InConstr1(const Vector<double>& x) {
//	double f = x[0] * x[1] - 4;
//	return f;
//}
//
//double InConstr2(const Vector<double>& x) {
//	double f = -0.25 * x[0] * x[0] + x[1] - 2;
//	return f;
//}

double FuncValue(const Vector<double>& x) {
	TStronginC3Problem stPr = TStronginC3Problem();
	double* point = x.GetData();
	std::cout << x << std::endl;
	return stPr.CalculateFunctionals(point, 3);
}

double InConstr1(const Vector<double>& x) {
	TStronginC3Problem stPr = TStronginC3Problem();
	double* point = x.GetData();
	return stPr.CalculateFunctionals(point, 0);
}

double InConstr2(const Vector<double>& x) {
	TStronginC3Problem stPr = TStronginC3Problem();
	double* point = x.GetData();
	return stPr.CalculateFunctionals(point, 1);
}

double InConstr3(const Vector<double>& x) {
	TStronginC3Problem stPr = TStronginC3Problem();
	double* point = x.GetData();
	return stPr.CalculateFunctionals(point, 2);
}

int main() {
	Vector<double> x(2);
	x[0] = 2;
	x[1] = 2;

	TaskPrep tp = TaskPrep(x, FuncValue);
	tp.SetIqConstr({ InConstr1 , InConstr2, InConstr3 });

	tp.Solve();

	system("pause");  // This line is valid only in Visual Studio IDE. Delete it, if in Linux environment.

	return 0;
}