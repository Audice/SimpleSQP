#ifndef __STRONGINC3PROBLEM_H__
#define __STRONGINC3PROBLEM_H__

#include "problem_interface.h"

class TStronginC3Problem : public IProblem
{
protected:

  int mDimension;
  bool mIsInitialized;
  static const int mSupportedDimension = 2;

public:

  TStronginC3Problem();

  virtual int SetConfigPath(const std::string& configPath);
  virtual int SetDimension(int dimension);
  virtual int GetDimension() const;
  virtual int Initialize();

  virtual void GetBounds(double* upper, double *lower);
  virtual int GetOptimumValue(double& value) const;
  virtual int GetOptimumPoint(double* x) const;

  virtual int GetNumberOfFunctions() const;
  virtual int GetNumberOfConstraints() const;
  virtual int GetNumberOfCriterions() const;

  virtual double CalculateFunctionals(const double* x, int fNumber);

  ~TStronginC3Problem();
};

extern "C" LIB_EXPORT_API IProblem* create();
extern "C" LIB_EXPORT_API void destroy(IProblem* ptr);

#endif
// - end of file ----------------------------------------------------------------------------------
