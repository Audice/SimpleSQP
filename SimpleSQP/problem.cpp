#include "problem.h"
#include <vector>
#include <array>

double objfun(double x1, double x2)
{
	double f = pow(x2-x1*x1,2)+pow(1-x1,2);
	return f;
}

//double FuncValue(const Vector<double>& x) {
//	double f = pow(x[1] - x[0] * x[0], 2) + pow(1 - x[0], 2);
//	return f;
//}

void IC_value(const Vector<double>& x) {

}

void Gradient(Vector<double> x, Vector<double>& c)
{
	//Добавил численный расчёт градиента под вектора разного размера
	finiteGradient(x, c, 1);
}

static const std::array<std::vector<int>, 4> coeff =
{ { {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} } };
static const std::array<std::vector<int>, 4> coeff2 =
{ { {1, -1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} } };
static const std::array<int, 4> dd = { 2, 12, 60, 840 };


void finiteGradient(const Vector<double>& x, Vector<double>& grad, int accuracy) {
	// accuracy can be 0, 1, 2, 3
	const double eps = 2.2204e-6;

	grad.resize(x.size());
	Vector<double>& xx = const_cast<Vector<double>&>(x);

	const int innerSteps = 2 * (accuracy + 1);
	const double ddVal = dd[accuracy] * eps;

	for (size_t d = 0; d < x.size(); d++) {
		grad[d] = 0;
		for (int s = 0; s < innerSteps; ++s)
		{
			double tmp = xx[d];
			xx[d] += coeff2[accuracy][s] * eps;
			grad[d] += coeff[accuracy][s] * 1;// FuncValue(xx);
			xx[d] = tmp;
		}
		grad[d] /= ddVal;
	}
}





// The Hessian Matrix of the object function
void Hessian(double x1, double x2, Matrix<double>& Q)
{
	double Hess[2][2] = {0};

	Hess[0][0] = 12*pow(x1,2)-4*x2+2;  
	Hess[0][1] = -4*x1;
	Hess[1][0] = -4*x1;
	Hess[1][1] = 2;

	for(int i = 0; i < Q.nrows(); i++)
	{
		for(int j = 0; j < Q.ncols(); j++)
		{
			Q[i][j] = Hess[i][j];
		}
	}
}

//The gradient of equality constraint function
void GradientCE(double x1, double x2, Matrix<double>& CE)
{

}

//The constant term of equality constraint function
void CE0(double x1, double x2, Vector<double>& ce0)
{
	
}


//for inequality constraint І»µИКЅФјКшµДПµКэ
void GradientCI(double x1, double x2, Matrix<double>& CI)  
{
	double gradie[2][2]={0};

	gradie[0][0] = x2;
	gradie[0][1] = -0.5*x1;
	gradie[1][0] = x1;
	gradie[1][1] = 1;

	for(int i = 0; i < CI.nrows(); i++)
	{
		for(int j = 0; j < CI.ncols(); j++)
		{
			CI[i][j] = gradie[i][j];
		}
	}

}



//void GradientCI(const Vector<Vector<double>>& constr, Matrix<double>& CI)
//{
//	Vector<Vector<double>> grads(constr.size());
//	for (int i = 0; i < grads.size(); i++) {
//		Vector<double> grad(grads[i].size());
//		finiteGradient(constr[i], grad, 1);
//	}
//	
//
//
//	double gradie[2][2] = { 0 };
//
//	gradie[0][0] = x2;
//	gradie[0][1] = -0.5 * x1;
//	gradie[1][0] = x1;
//	gradie[1][1] = 1;
//
//	for (int i = 0; i < CI.nrows(); i++)
//	{
//		for (int j = 0; j < CI.ncols(); j++)
//		{
//			CI[i][j] = gradie[i][j];
//		}
//	}
//
//}




//І»µИКЅФјКшµДіЈКэПо
void CI0(double x1, double x2, Vector<double>& ci0)
{
	double ci0x[2] = {0};
	
	// The inequality constraint 1 І»µИКЅФјКш 1
	double ine1 = x1*x2-4;
	
	// The inequality constraint 2 І»µИКЅФјКш 2
	double ine2 = -0.25*pow(x1,2)+x2-2;

	ci0x[0] = ine1;
    ci0x[1] = ine2;

	for(int i = 0; i< ci0.size(); i++)
	{
		ci0[i] = ci0x[i];
	}
}

double Quadfun(Matrix<double> Q, Vector<double> c, Vector<double> d)
{
	//f(x)=1/2*(d^T)*Q*d+(c^T)*d
	double sum = 0.0;
	for (int i = 0; i < Q.nrows(); i++)
	{
		for (int j = 0; j < Q.ncols(); j++)
		{
			sum += d[i] * Q[i][j] * d[j];
		}
	}

	sum *= 0.5;	

	for (int i = 0; i < c.size(); i++)
	{
		sum += c[i] * d[i];
	}

	return sum;
}
