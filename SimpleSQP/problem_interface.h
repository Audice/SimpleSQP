/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2016 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      problem_interface.h                                         //
//                                                                         //
//  Purpose:   Header file for ExaMin problem interface                    //
//                                                                         //
//                                                                         //
//  Author(s): Sovrasov V.                                                 //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
\file problem_interface.h

\authors Соврасов В.
\date 2016
\copyright ННГУ им. Н.И. Лобачевского

\brief Объявление абстрактного класса #TIProblem

\details Объявление абстрактного класса #TIProblem и сопутствующих типов данных
*/

#ifndef __PROBLEMINTERFACE_H__
#define __PROBLEMINTERFACE_H__

#include <vector>
#include <string>
#include <stdexcept>

/**
Базовый класс-интерфейс, от которого наследуются классы, описывающие задачи оптимизации.

В классе #TIProblem описаны прототипы методов, которые должны быть реализованы в подключамых модулях с задачами.
*/
class IProblem
{
public:

  /// Код ошибки, возвращаемый, если операция завершена успешно
  static const int OK = 0;
  /** Код ошибки, возвращаемый методами #GetOptimumValue и #GetOptimumPoint,
  если соответствующие параметры задачи не определены,
  */
  static const int UNDEFINED = -1;
  /// Код ошибки, возвращаемый, если операция не выполнена
  static const int ERROR = -2;

  /** Задание пути до конфигурационного файла

  Данный метод должн вызываться перед #Initialize
  \param[in] configPath строка, содержащая путь к конфигурационному файлу задачи
  \return Код ошибки
  */
  virtual int SetConfigPath(const std::string& configPath) = 0;

  /** Метод задаёт размерность задачи

  Данный метод должен вызываться перед #Initialize. Размерность должна быть в
  списке поддерживаемых.
  \param[in] dimension размерность задачи
  \return Код ошибки
  */
  virtual int SetDimension(int dimension) = 0;
  ///Возвращает размерность задачи, можно вызывать после #Initialize
  virtual int GetDimension() const = 0;
  ///Инициализация задачи
  virtual int Initialize() = 0;

  /** Метод возвращает границы области поиска
  */
  virtual void GetBounds(double* lower, double *upper) = 0;
  /** Метод возвращает значение целевой функции в точке глобального минимума
  \param[out] value оптимальное значение
  \return Код ошибки (#OK или #UNDEFINED)
  */
  virtual int GetOptimumValue(double& value) const = 0;
  /** Метод возвращает значение функции с номером index в точке глобального минимума
  \param[out] value оптимальное значение
  \return Код ошибки (#OK или #UNDEFINED)
  */
  virtual int GetOptimumValue(double& value, int index) const
  {
    return IProblem::UNDEFINED;
  }
  /** Метод возвращает координаты точки глобального минимума целевой функции
  \param[out] y точка, в которой достигается оптимальное значение
  \return Код ошибки (#OK или #UNDEFINED)
  */
  virtual int GetOptimumPoint(double* y) const = 0;
  /** Метод возвращает координаты всех точек глобального минимума целевой функции
  и их количество
  \param[out] y координаты точек, в которых достигается оптимальное значение
  \param[out] n количество точек, в которых достигается оптимальное значение
  \return Код ошибки (#OK или #UNDEFINED)
  */
  virtual int GetAllOptimumPoint(double* y, int& n) const
  {
    return IProblem::UNDEFINED;
  }

  /** Метод возвращает число общее функций в задаче (оно равно число ограничений + число критериев)
  \return Число функций
  */
  virtual int GetNumberOfFunctions() const = 0;
  /** Метод возвращает число ограничений в задаче
  \return Число ограничений
  */
  virtual int GetNumberOfConstraints() const = 0;
  /** Метод возвращает число критериев в задаче
  \return Число критериев
  */
  virtual int GetNumberOfCriterions() const = 0;

  /** Метод, вычисляющий функции задачи

  \param[in] y Точка, в которой необходимо вычислить значение
  \param[in] fNumber Номер вычисляемой функции. 0 соответствует первому ограничению,
  #GetNumberOfFunctions() - 1 -- последнему критерию
  \return Значение функции с указанным номером
  */
  virtual double CalculateFunctionals(const double* y, int fNumber);



  virtual bool isOptimal(const double* y, double *minv, double *maxv) { return true; }

  ///Деструктор
  virtual ~IProblem() = 0;
};

//// ------------------------------------------------------------------------------------------------

/**
Базовый класс-интерфейс, от которого наследуются классы, описывающие задачи оптимизации
использующие для вычислений GPU.

В классе #IGPUProblem описаны прототипы методов, которые должны быть реализованы в подключамых
модулях с решением на GPU.
*/
class IGPUProblem : public IProblem
{
public:

  /** Метод, вычисляющий функции задачи в нескольких точках одновременно

\param[in] y массив, содержащий последовательно записанные многомерные точки, в которых необходимо
вычислить функционалы задачи
\param[in] Номер вычисляемой функции
\param[in] numPoints количество передаваемых точек
\param[out] values массив, в который будут записаны вычисленные значения функционалов
*/
  virtual void CalculateFunctionals(double* y, int fNumber, int& numPoints, double* values)
  {
    throw std::runtime_error(std::string("Required overload of the following method is not implemented: ")
      + std::string(__FUNCTION__));
  }
};


//// ------------------------------------------------------------------------------------------------

/**
Базовый класс-интерфейс, от которого наследуются классы, описывающие задачи оптимизации
с дискретными параметрами.

В классе #TIProblem описаны прототипы методов, которые должны быть реализованы в
подключамых модулях в которых задачи используют дискретный параметр.

Для дискретных параметра задаются все допустимые значения.
Дискретные параметры являются последними в векторе параметров y.
*/

class IIntegerProgrammingProblem : public IGPUProblem
{
public:
  /// Код ошибки, возвращаемый, если попытались получить значения для недискретного параметра
  static const int ERROR_DISCRETE_VALUE = -201;
  /// Возвращает число дискретных параметров, дискретные параметры всегда последние в векторе y
  virtual int GetNumberOfDiscreteVariable() = 0;
  /**
  Возвращает число значений дискретного параметра discreteVariable.
  GetDimension() возвращает общее число параметров.
  (GetDimension() - GetNumberOfDiscreteVariable()) - номер начальной дискретной переменной
  Для не дискретных переменных == -1
  */
  virtual int GetNumberOfValues(int discreteVariable) = 0;
  /**
  Определяет значения дискретного параметра с номером discreteVariable
  Возвращает код ошибки.
  \param[out] values массив, в который будут сохранены значения дискретного параметра
  нулевой элемент это левая граница, поледний элемент это правая граница.
  */
  virtual int GetAllDiscreteValues(int discreteVariable, double* values) = 0;
  /**
  Определяет значения дискретного параметра с номером discreteVariable после номера previousNumber
  Возвращает код ошибки.
  \param[in] previousNumber - номер значения после которого возвращается значение
  -2 - значение по умолчанию, возвращает следующее значение
  -1 - возвращает после -1, т.е. левую границу области
  \param[out] value переменная в которую сохраняется значение дискретного параметра
  */
  virtual int GetNextDiscreteValues(int* mCurrentDiscreteValueIndex, double& value, int discreteVariable, int previousNumber = -2) = 0;
  /// Проверяет является ли value допустимым значением для параметра с номером discreteVariable
  virtual bool IsPermissibleValue(double value, int discreteVariable) = 0;
};

////
//// ------------------------------------------------------------------------------------------------
//void IGPUProblem::CalculateFunctionals(double* y, int fNumber, int& numPoints, double* values)
//{
//  throw std::runtime_error(std::string("Required overload of the following method is not implemented: ")
//    + std::string(__FUNCTION__));
//}

// ------------------------------------------------------------------------------------------------
inline double IProblem::CalculateFunctionals(const double* y, int fNumber)
{
  throw std::runtime_error(std::string("Required overload of the following method is not implemented: ")
    + std::string(__FUNCTION__));
}

// ------------------------------------------------------------------------------------------------
inline IProblem::~IProblem() {}

///Тип функции-фабрики, которая экспортируется подключаемой библиотекой с задачей
typedef IProblem* create_t();
///Тип функции-деструктора, которая экспортируется подключаемой библиотекой с задачей
typedef void destroy_t(IProblem*);

///Префикс для фуккций, экспортируемых подключаемой библиотекой с задачей
#ifdef WIN32
#define LIB_EXPORT_API __declspec(dllexport)
#else
#define LIB_EXPORT_API
#endif

#endif
// - end of file ----------------------------------------------------------------------------------
