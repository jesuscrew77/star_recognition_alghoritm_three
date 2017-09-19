#ifndef MATHFUNC_H
#define MATHFUNC_H
#include <utility>
#include <iterator>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <QtGlobal>
#include <QDateTime>
#include <QDebug>
#include <QPointF>



typedef int (*CorrectBoardQuat)(wchar_t *patch_cat, /*Имя папки c бортовыми каталогами*/
                                wchar_t *patch_dtmi, /*Имя папки для создания файлов с тестовой информацией, если NULL - файлы не создаются*/
                                double ArrBokz[][8], /*Массив кватернионов*/
unsigned int &nBokz, /*Число измерений кватернионов*/
unsigned int nBokzMax, /*Макс. размер массива кватернионов*/
float ArrLoc[][77], /*Массив локализованных объектов*/
unsigned int nLoc, /*Размер массива локализованных объектов*/
const double ArrVect[][18], /*Массив ПДСМ*/
unsigned int nVect, /*Размер массива ПДСМ*/
unsigned int &CurrentIndex, /*Текущий индекс массива кватернионов*/
unsigned char &Stop); /*Признак останова программы*/



enum AXIS_TYPE
{
    X_AXIS,
    Y_AXIS,
    Z_AXIS
};


struct RecognizedInfo
{
    double distance = 0;
    int numInCat = 0;
    bool checked = false;
    bool recognized = false;
};




namespace BOKZMath {

constexpr const double MSecsInSec = 1000;
constexpr const double transToAngularMinutes = 3437.7467747131374;
constexpr const double transToGradus = 57.29577957855229;
constexpr const double transToRadians = 0.0174532925;
constexpr const double julianZeroDate = 2451545.;
constexpr const double countOfDaysInCentury = 36525.;




float correctDTMIFloat(float element) noexcept;

void correctLocArray(float locArray[16][4]) noexcept;

void correctAngularVArray(float Wo[3]) noexcept;

quint32 correctDTMITimePr(quint32 incrTimePr) noexcept;

double calculateJulianDate(const QDateTime& dateTime) noexcept;

double calculateJulianCenturies(const QDateTime& dateTime) noexcept;

void multMatrixVector(const double M[3][3], const double V[3], double VR[3]) noexcept;

double calculateDatePr(const QDateTime& dateTime, int timePr) noexcept;

double toTDateTimeValue(const QDateTime& datetime) noexcept;

QDateTime fromTDateTime(double value) noexcept;

/*смотри http://law.rufox.ru/view/standarti/1672.htm*/
void calculateMoonDirection(double JD_DEV, double* pMoonI) noexcept;

void calculateSunDirection(double JD_DEV, double* pSunI) noexcept;

double starsTime2(double JD) noexcept;

void getGSKISKMatrix(double JD, double MG[3][3]) noexcept;

void  transSpeedGSKtoISK(double JD, double CG[3], double VG[3], double VI[3]) noexcept;

void  transCoordsGSKtoISK(double JD, double CG[3], double CI[3]) noexcept;

void quatToMatr(const double Quat[], double M_ornt[3][3]) noexcept;

double atan2m(double yf,double xf) noexcept;

void matrToEkvA(double M_ornt[3][3], double& al, double& dl, double& Az) noexcept;

void quatToEkvA(const double Quat[4], double EkvA[3]) noexcept;

double acosm(double xf) noexcept;

void multMatrix(const double Matr1[3][3],const double Matr2[3][3], double Matr[3][3]) noexcept;

void getAngularDisplacementFromOrientMatr(const double M_ornt_pr[3][3],const double M_ornt[3][3], double Wop [3]) noexcept;

double calculateAngleAxis(const double quat1[], const double quat2[], AXIS_TYPE axis) noexcept;

QVector <QPointF> createHistogramm(size_t histSize, float shiftRange, const QVector <float>& firstX,const QVector <float>& firstY,const QVector <float>& secondX,const QVector <float>& secondY);

QVector <QPointF> createHistogramm(size_t histSize, float shiftRange, const QVector<QPointF>& firstCoords, const QVector<QPointF>& secondCoords);

qint64 roundUp(qint64 numToRound, qint64 multiple);

QDateTime timePrToDateTime(const double timePr, const QDateTime& timePrStartData) noexcept;

QDateTime timePrToDateTime(const quint32 timePr, const QDateTime& timePrStartData) noexcept;

QDateTime timePrToDateTime(const quint64 timePr, const QDateTime& timePrStartData) noexcept;

QVector< QVector<float> > calcTransitionMatrix(double pointAlpha, double pointBeta, double pointAzimut) noexcept;

double calcScalarProduct(double l_oz, double l_st, double m_oz, double m_st, double n_oz, double n_st);

QVector <quint32> firstMinDistanceTable(RecognizedInfo** distMatrix, quint32 countOfMins, quint32 objectIndex, quint32 size);


template <class InputIterator, class T>
std::pair<T, T> calculateMeanStDv (InputIterator first, InputIterator last, T init) noexcept
{
    if (first == last) return std::pair<T,T>(*first, T());

    T dispersio = 0;
    for (InputIterator i = first;i < last;i ++)
    {
        init += *i;
        dispersio += pow(*i, 2);
    }
    auto count = std::distance(first,last);
    T mean = init/count;
    dispersio = (dispersio / count) - pow(mean, 2);

    return std::pair <T,T> (mean, sqrt(dispersio));

}



template <class InputIterator, class T, class UnaryOperation>
std::pair<T,T> calculateMeanStDv (InputIterator first, InputIterator last, T init, UnaryOperation extractWtC) noexcept
{
    if (first == last) return std::pair<T,T>(extractWtC(first), T());

    T dispersio = 0;
    for (InputIterator i = first;i < last;i ++)
    {
        init += extractWtC(i);
        dispersio += pow(extractWtC(i), 2);
    }
    auto count = std::distance(first,last);
    T mean = init/count;
    dispersio = (dispersio / count) - pow(mean, 2);

    return std::pair <T,T> (mean, sqrt(dispersio));

}
template <class InputIterator, class T, class UnaryOperation>
T calculateMean (InputIterator first, InputIterator last, T init, UnaryOperation extractWtC) noexcept
{
    if (first == last) return extractWtC(first);

    for (InputIterator i = first;i < last;i ++)
    {
        init += extractWtC(i);
    }
    auto count = std::distance(first,last);
    return init/count;
}

template <class InputIterator, class T>
T calculateMean (InputIterator first, InputIterator last, T init) noexcept
{
    if (first == last) return *first;

    for (InputIterator i = first;i < last;i ++)
    {
        init += *i;
    }
    auto count = std::distance(first,last);
    return init/count;
}


/*функция считающая пересечения в векторах, и оставляющая только общие значения по времени для каждого*/
template <typename Data, typename Functor>
void setIntersectionData(QVector <Data> &intersData, Functor&& functor)
{

    /*если данные только для одного прибора, пересечения считать не с чем, выходим*/
    if (intersData.size() == 1)
    {
        return;
    }

    /*формируем первый (нулевой) вектор, как вектор без непересекающимеся значениями относительно остальных*/
    for (qint32 i = 1;i < intersData.size();i ++)
    {
        Data intersectedVector;// вектор, который будет содержать непересекающиеся значения, его будем присваивать первому вектору
        std::set_intersection (intersData[0].begin(), intersData[0].end(), intersData[i].begin(), intersData[i].end(),std::back_inserter(intersectedVector),
                functor);
        intersData[0] = std::move(intersectedVector);
    }

    for (qint32 i = 1;i < intersData.size();i ++)
    {
        Data intersectedVector;// вектор, который будет содержать непересекающиеся значения, его будем присваивать оставшимся векторам
        std::set_intersection (intersData[i].begin(), intersData[i].end(), intersData[0].begin(), intersData[0].end(), std::back_inserter(intersectedVector),
                functor);
        intersData[i] = std::move(intersectedVector);
    }
}
template <typename Matrix>
std::pair <quint32, quint32> matrixMax(const Matrix& matrix)
{
    if (matrix.size() < 1)
    {
        return std::pair <quint32, quint32> (0, 0);
    }
    quint32 row = 0;
    quint32 column = 0;
    auto maxMatrix = decltype(matrix[0][0]){0};
    for (int i = 0 ; i < matrix.size(); i ++)
    {
        for (int j = 0 ; j < matrix[i].size(); j ++)
        {
            if (maxMatrix < matrix[i][j])
            {
                maxMatrix = matrix[i][j];
                row = i;
                column = j;
            }
        }
    }
    return std::pair <quint32, quint32> (row, column);
}

template <typename Matrix>
std::pair <quint32, quint32> matrixMin(const Matrix& matrix)
{
    if (matrix.size() < 1)
    {
        return std::pair <quint32, quint32> (0, 0);
    }
    quint32 row = 0;
    quint32 column = 0;
    auto minMatrix = matrix[0][0];

    for (int i = 0 ; i < matrix.size(); i ++)
    {
        for (int j = 0 ; j < matrix[i].size(); j ++)
        {
            if (minMatrix < matrix[i][j])
            {
                minMatrix = matrix[i][j];
                row = i;
                column = j;
            }
        }
    }
    return std::pair <quint32, quint32> (row, column);
}


template <class InputIterator>
QVector <quint32> firstMinElements (InputIterator first, InputIterator last, quint32 countOfMins)
{
    QVector <float> vec(countOfMins, std::numeric_limits <float>::max());
    QVector <quint32> minIndexes(countOfMins, std::numeric_limits <quint32>::max());

    for (auto it = first; it != last; it++)
    {
        for (int i = 0; i < countOfMins; i++)
        {
            if (*it < vec[i])
            {
                for (int j = countOfMins - 1; j >= i; j--)
                {
                    if (j == i)
                    {
                        vec[j] = *it;
                        minIndexes[j] = it - first;
                        break;
                    }
                    else
                    {
                        vec[j] = vec[j - 1];
                        minIndexes[j] = minIndexes[j - 1];
                    }
                }
                break;
            }
        }
    }
    return minIndexes;
}


template <class InputIterator, typename Functor>
QVector <quint32> firstMinElements (InputIterator first, InputIterator last, quint32 countOfMins, Functor functor)
{
    QVector <float> vec(countOfMins, std::numeric_limits <decltype(functor(*first))>::max());
    QVector <quint32> minIndexes(countOfMins, std::numeric_limits <quint32>::max());

    for (auto it = first; it != last; it++)
    {
        for (int i = 0; i < countOfMins; i++)
        {
            if (functor(*it) < vec[i])
            {
                for (int j = countOfMins - 1; j >= i; j--)
                {
                    if (j == i)
                    {
                        vec[j] = functor(*it);
                        minIndexes[j] = it - first;
                        break;
                    }
                    else
                    {
                        vec[j] = vec[j - 1];
                        minIndexes[j] = minIndexes[j - 1];
                    }
                }
                break;
            }
        }
    }
    return minIndexes;
}
//        if (*it < std::get <0> (tuple))
//        {
//            std::get <3> (tuple) = std::get <2> (tuple);
//            std::get <2> (tuple) = std::get <1> (tuple);
//            std::get <1> (tuple) = std::get <0> (tuple);
//            std::get <0> (tuple) = *it;
//        }
//        else if (*it < std::get <1> (tuple))
//        {
//            std::get <3> (tuple) = std::get <2> (tuple);
//            std::get <2> (tuple) = std::get <1> (tuple);
//            std::get <1> (tuple) = *it;
//        }
//        else if (*it < std::get <2> (tuple))
//        {
//            std::get <3> (tuple) = std::get <2> (tuple);
//            std::get <2> (tuple) = *it;
//        }
//        else if (*it < std::get <3> (tuple))
//        {
//            std::get <3> (tuple) = *it;
//        }



//template<class T, T v>
//struct ttrue_false {
//    static constexpr T value = v;
//};

//template<>
//struct ttrue_false<bool,true> {
//    static constexpr bool value = true;
//};
//typedef ttrue_false<bool,true> true_typep;

//template<>
//struct ttrue_false<bool,false> {
//    static constexpr bool value = false;
//};
//typedef ttrue_false<bool,false> false_typep;


//template <typename T, typename U>
//struct is_integral  : false_typep {};

//template <>
//struct is_integral<int,int>  : true_typep {};

//template <typename T, typename U>
//struct is_same : false_typep{};

//template <typename T>
//struct is_same<T,T> : true_typep{};

//template<class T ,
//         class = typename std::enable_if<std::is_integral<T>::value>::type >
//T foo3(T t) // обратите внимание, сигнатура функции не меняется
//{
//    return t;
//}


}
template<class ForwardIt, class T>
ForwardIt lower_bound(ForwardIt first, ForwardIt last, const T& value)
{
    ForwardIt it;
    typename std::iterator_traits<ForwardIt>::difference_type count, step;
    count = std::distance(first, last);

    while (count > 0) {
        it = first;
        step = count / 2;
        std::advance(it, step);
        if (*it < value) {
            first = ++it;
            count -= step + 1;
        }
        else
            count = step;
    }
    return first;
}
#endif // MATHFUNC_H
