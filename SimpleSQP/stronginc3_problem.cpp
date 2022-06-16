#include "stronginc3_problem.h"

#include <math.h>

// ------------------------------------------------------------------------------------------------
TStronginC3Problem::TStronginC3Problem()
{
  mIsInitialized = false;
  mDimension = mSupportedDimension;
}

// ------------------------------------------------------------------------------------------------
int TStronginC3Problem::SetConfigPath(const std::string& configPath)
{
  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TStronginC3Problem::SetDimension(int dimension)
{
  if(dimension == mSupportedDimension)
  {
    mDimension = dimension;
    return IProblem::OK;
  }
  else
    return IProblem::ERROR;
}

// ------------------------------------------------------------------------------------------------
int TStronginC3Problem::GetDimension() const
{
  return mDimension;
}

// ------------------------------------------------------------------------------------------------
int TStronginC3Problem::Initialize()
{
  if (mDimension == mSupportedDimension)
  {
    mIsInitialized = true;
    return IProblem::OK;
  }
  else
    return IProblem::ERROR;
}

// ------------------------------------------------------------------------------------------------
void TStronginC3Problem::GetBounds(double* lower, double *upper)
{
  if (mIsInitialized)
  {
    lower[0] =  0.;
    upper[0] =  4.;
    lower[1] = -1.;
    upper[1] =  3.;
  }
}

// ------------------------------------------------------------------------------------------------
int TStronginC3Problem::GetOptimumValue(double& value) const
{
  if (!mIsInitialized)
    return IProblem::UNDEFINED;

  value = -1.489444;

  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TStronginC3Problem::GetOptimumPoint(double* point) const
{
  if (!mIsInitialized)
    return IProblem::UNDEFINED;

  point[0] = 0.941176;
  point[1] = 0.941176;

  return IProblem::OK;
}

// ------------------------------------------------------------------------------------------------
int TStronginC3Problem::GetNumberOfFunctions() const
{
  return 4;
}

// ------------------------------------------------------------------------------------------------
int TStronginC3Problem::GetNumberOfConstraints() const
{
  return 3;
}

// ------------------------------------------------------------------------------------------------
int TStronginC3Problem::GetNumberOfCriterions() const
{
  return 1;
}

// ------------------------------------------------------------------------------------------------
double TStronginC3Problem::CalculateFunctionals(const double* y, int fNumber)
{
  double res = 0.0;
  double x1 = y[0], x2 = y[1];
  switch (fNumber)
  {
  case 0: // constraint 1
    res = 0.01 * ((x1 - 2.2) * (x1 - 2.2) + (x2 - 1.2) * (x2 - 1.2) - 2.25);
    break;
  case 1: // constraint 2
    res = 100.0 * (1.0 - ((x1 - 2.0) / 1.2) * ((x1 - 2.0) / 1.2) -
      (x2 / 2.0) * (x2 / 2.0));
    break;
  case 2: // constraint 3
    res = 10.0 * (x2 - 1.5 - 1.5 * sin(6.283 * (x1 - 1.75)));
    break;
  case 3: // criterion
  {
    double t1 = pow(0.5 * x1 - 0.5, 4.0);
    double t2 = pow(x2 - 1.0, 4.0);
    res = 1.5 * x1 * x1 * exp(1.0 - x1 * x1 - 20.25 * (x1 - x2) * (x1 - x2));
    res = res + t1 * t2 * exp(2.0 - t1 - t2);
    res = -res;
  }
    break;
  }

  return res;
}

// ------------------------------------------------------------------------------------------------
TStronginC3Problem::~TStronginC3Problem()
{

}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API IProblem* create()
{
  return new TStronginC3Problem();
}

// ------------------------------------------------------------------------------------------------
LIB_EXPORT_API void destroy(IProblem* ptr)
{
  delete ptr;
}
// - end of file ----------------------------------------------------------------------------------
