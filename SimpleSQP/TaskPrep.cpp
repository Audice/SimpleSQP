#include "TaskPrep.h"
#include "QuadProg.h"



void TaskPrep::SetEqConstr(std::vector<PFunc> qeConstr) {
	this->_equalConstraint = qeConstr;
}
void TaskPrep::SetIqConstr(std::vector<PFunc> iqConstr) {
	this->_inequalConstraint = iqConstr;
}


void TaskPrep::FiniteGradient(const QuadProg::Vector<double>& x, QuadProg::Vector<double>& grad, PFunc func, int accuracy) {
	// accuracy can be 0, 1, 2, 3
	const double eps = 2.2204e-6;

	grad.resize(x.size());
	QuadProg::Vector<double>& xx = const_cast<QuadProg::Vector<double>&>(x);

	const int innerSteps = 2 * (accuracy + 1);
	const double ddVal = dd[accuracy] * eps;

	for (size_t d = 0; d < x.size(); d++) {
		grad[d] = 0;
		for (int s = 0; s < innerSteps; ++s)
		{
			double tmp = xx[d];
			xx[d] += coeff2[accuracy][s] * eps;
			grad[d] += coeff[accuracy][s] * func(xx);
			xx[d] = tmp;
		}
		grad[d] /= ddVal;
	}
}

void TaskPrep::Gradient(const QuadProg::Vector<double> x, QuadProg::Vector<double>& grad)
{
	//Добавил численный расчёт градиента под вектора разного размера
	this->FiniteGradient(x, grad, this->_problem, 1);
}



//void TaskPrep::IC_value(const QuadProg::Vector<double>& x) {
//
//}


// The Hessian Matrix of the object function
void TaskPrep::Hessian(const QuadProg::Vector<double>& x, QuadProg::Matrix<double>& Q)
{
	/*
	double Hess[2][2] = { 0 };

	Hess[0][0] = 12 * pow(x[0], 2) - 4 * x[1] + 2;
	Hess[0][1] = -4 * x[0];
	Hess[1][0] = -4 * x[0];
	Hess[1][1] = 2;

	for (int i = 0; i < Q.nrows(); i++)
	{
		for (int j = 0; j < Q.ncols(); j++)
		{
			Q[i][j] = Hess[i][j];
		}
	}
	*/
	this->FiniteHessian(x, Q, 0);

	std::cout << Q[0][0] << " " << Q[0][1] << std::endl;
	std::cout << Q[1][0] << " " << Q[1][1] << std::endl;
	int s = 0;

}

void TaskPrep::FiniteHessian(const QuadProg::Vector<double>& x, QuadProg::Matrix<double>& hessian, int accuracy = 0) {
	
	const double eps = std::numeric_limits<double>::epsilon() * 10e7;

	hessian.resize(x.size(), x.size()); //???
	QuadProg::Vector<double>& xx = const_cast<QuadProg::Vector<double>&>(x);

	if (accuracy == 0) {
		for (size_t i = 0; i < x.size(); i++) {
			for (size_t j = 0; j < x.size(); j++) {
				double tmpi = xx[i];
				double tmpj = xx[j];

				double f4 = this->_problem(xx);
				xx[i] += eps;
				xx[j] += eps;
				double f1 = this->_problem(xx);
				xx[j] -= eps;
				double f2 = this->_problem(xx);
				xx[j] += eps;
				xx[i] -= eps;
				double f3 = this->_problem(xx);
				hessian[i][j] = (f1 - f2 - f3 + f4) / (eps * eps);
				xx[i] = tmpi;
				xx[j] = tmpj;
			}
		}
	}
	else {
		for (size_t i = 0; i < x.size(); i++) {
			for (size_t j = 0; j < x.size(); j++) {
				double tmpi = xx[i];
				double tmpj = xx[j];

				double term_1 = 0.0;
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1 * eps;  xx[j] += -2 * eps;  term_1 += this->_problem(xx);
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2 * eps;  xx[j] += -1 * eps;  term_1 += this->_problem(xx);
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2 * eps; xx[j] += 1 * eps;   term_1 += this->_problem(xx);
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1 * eps; xx[j] += 2 * eps;   term_1 += this->_problem(xx);

				double term_2 = 0.0;
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1 * eps; xx[j] += -2 * eps;  term_2 += this->_problem(xx);
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2 * eps; xx[j] += -1 * eps;  term_2 += this->_problem(xx);
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1 * eps;  xx[j] += 2 * eps;   term_2 += this->_problem(xx);
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2 * eps;  xx[j] += 1 * eps;   term_2 += this->_problem(xx);

				double term_3 = 0.0;
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2 * eps;  xx[j] += -2 * eps;  term_3 += this->_problem(xx);
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2 * eps; xx[j] += 2 * eps;   term_3 += this->_problem(xx);
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += -2 * eps; xx[j] += -2 * eps;  term_3 -= this->_problem(xx);
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += 2 * eps;  xx[j] += 2 * eps;   term_3 -= this->_problem(xx);

				double term_4 = 0.0;
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1 * eps; xx[j] += -1 * eps;  term_4 += this->_problem(xx);
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1 * eps;  xx[j] += 1 * eps;   term_4 += this->_problem(xx);
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += 1 * eps;  xx[j] += -1 * eps;  term_4 -= this->_problem(xx);
				xx[i] = tmpi; xx[j] = tmpj; xx[i] += -1 * eps; xx[j] += 1 * eps;   term_4 -= this->_problem(xx);

				xx[i] = tmpi;
				xx[j] = tmpj;

				hessian[i][j] = (-63 * term_1 + 63 * term_2 + 44 * term_3 + 74 * term_4) / (600.0 * eps * eps);
			}
		}
	}
}




//The gradient of equality constraint function
void TaskPrep::GradientCE(const QuadProg::Vector<double>& x, QuadProg::Matrix<double>& CE)
{
	if (this->_equalConstraint.size() > 0) {
		for (int i = 0; i < this->_equalConstraint.size(); i++) {
			QuadProg::Vector<double> localConstrGrad(x.size());
			this->FiniteGradient(x, localConstrGrad, this->_equalConstraint[i], 1);
			for (int j = 0; j < CE.nrows(); j++)
				CE[j][i] = localConstrGrad[j];
		}
	}
	else {
		//Нет ограничений
	}
}

//The constant term of equality constraint function
void TaskPrep::CE0(const QuadProg::Vector<double>& x, QuadProg::Vector<double>& ce0)
{
	if (this->_equalConstraint.size() > 0) {
		if (this->_equalConstraint.size() == ce0.size())
			for (int i = 0; i < ce0.size(); i++)
				ce0[i] = _equalConstraint[i](x);
		else
			throw std::invalid_argument("Не совпадает количество ограничений!");
	}


}


//for inequality constraint І»µИКЅФјКшµДПµКэ
void TaskPrep::GradientCI(const QuadProg::Vector<double>& x, QuadProg::Matrix<double>& CI)
{	
	if (this->_inequalConstraint.size() > 0) {
		for (int i = 0; i < this->_inequalConstraint.size(); i++) {
			QuadProg::Vector<double> localConstrGrad(x.size());
			this->FiniteGradient(x, localConstrGrad, this->_inequalConstraint[i], 1);
			for (int j = 0; j < CI.nrows(); j++)
				CI[j][i] = localConstrGrad[j];
		}
	}
	else {
		//Нет ограничений
	}
}


void TaskPrep::CI0(const QuadProg::Vector<double>& x, QuadProg::Vector<double>& ci0)
{
	if (this->_inequalConstraint.size() > 0) {
		if (this->_inequalConstraint.size() == ci0.size())
			for (int i = 0; i < ci0.size(); i++)
				ci0[i] = _inequalConstraint[i](x);
		else
			throw std::invalid_argument("Не совпадает количество ограничений!");
	}
}

double TaskPrep::Quadfun(QuadProg::Matrix<double> Q, QuadProg::Vector<double> c, QuadProg::Vector<double> d)
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

void TaskPrep::Solve()
{
	QuadProg::Matrix<double> Q, CE, CI;

	QuadProg::Vector<double> c, ce0, ci0;

	QuadProg::Vector<double> d;  //The optimal solution of a Quadratic Programming(QP)

	//The initial start point:
	QuadProg::Vector<double> x = QuadProg::Vector<double>(this->_startPoint);


	int n;   // The number of variables in objective function
	int m;   // The number of equality constraints
	int p;   // The number of inequality constraints

	//The value of quadratic function.
	double sum = 0.0;

	n = 2;
	d.resize(n);

	// The Hessian Matrix of Objective function
	Q.resize(n, n);
	c.resize(n);

	//The number of equality constraints
	m = this->_equalConstraint.size();
	//coefficient matrix of equality constraints 
	CE.resize(n, m);
	ce0.resize(m);

	//The number of inequality constraint
	p = this->_inequalConstraint.size();
	//coefficient matrix of inequality constraints 
	CI.resize(n, p);
	ci0.resize(p);

	double epsilon = 1e-6;
	double err = 10;

	// The loop of the SQP
	while (err > epsilon)
	{
		Hessian(x, Q);
		Gradient(x, c);



		//The equality constraint
		GradientCE(x, CE);
		CE0(x, ce0);

		//The inequality constraint
		GradientCI(x, CI);
		CI0(x, ci0);

		// Call a QP-solver in the SQP loop.
		solve_quadprog(Q, c, CE, ce0, CI, ci0, d);

		std::cout << "x = [" << x[0] << " " << x[1] << "]" << std::endl;
		std::cout << "d = [" << d[0] << " " << d[1] << "]" << std::endl;

		/* FOR DOUBLE CHECKING COST since in the solve_quadprog routine the matrix Q is modified */
		Hessian(x, Q);

		QuadProg::Vector<double> tmp(x.size());
		tmp[0] = x[0] + d[0];
		tmp[1] = x[1] + d[1];


		/*Calculate the value of the Quadratic Function in the approximating QP problem.*/
		sum = Quadfun(Q, c, d);
		std::cout << "Quadratic Fun(" << d[0] << "," << d[1] << ") =" << sum << std::endl;

		double temp = 0.0;
		for (int i = 0; i < d.size(); i++)
		{
			temp = +d[i] * d[i];
		}

		err = pow(temp, 0.5);

		std::cout << "x_new = [" << tmp[0] << " " << tmp[1] << "]" << std::endl;

		double objf = this->_problem(tmp); //objfun(xa, xb);
		std::cout << "f(" << tmp[0] << "," << tmp[1] << ") = " << objf << std::endl;

		std::cout << "err= " << err << std::endl;
		std::cout << std::endl;

		//system("pause");

		x[0] = tmp[0];
		x[1] = tmp[1];

	}
}
