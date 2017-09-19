#ifndef STARRECALGH_H
#define STARRECALGH_H
#include <mathfunc.h>
#include <QDebug>
#include <QPointF>
#include <QTextStream>
#include <QFile>
#include <QIODevice>
#include <QDataStream>
#include <QElapsedTimer>
#include <QDir>
#include <QVector2D>
#include <catalog.h>
//#define _TEST_

#pragma pack(push,1)

struct StarInfoTwo
{
    quint16 id;
    float angle; // град
    quint16 id1;
    float distance1; // град
    quint16 id2;
    float distance2; // град
};
#pragma pack(pop)

class StarRecognizeThree
{
public:
    StarRecognizeThree(const QString& filename) {
        bool status;
        QString error;
        catalog.openCatalog(filename ,status, error);
        if (!status) {
            qDebug() << error;
        }
    }
    void makeTest();

private:
    Catalog catalog;
    QVector <double> lf;
    QVector <double> mf;
    QVector <double> nf;

    QVector <StarInfoTwo> prepareStarInfoTwoData(float focus, float pixelSize);
    bool checkAngleInRange (double angle, double compAngle, double range);
    bool recognizeStar(const QVector <QPointF>& starsCoords, const QVector <StarInfoTwo>& vec, RecognizedInfo** recMatrix, quint32& i , quint32 f, quint32 s, quint32& starNumber, quint32& nPos);
    void testRecognitionCalibrationFile(const QString& fileName, QVector <StarInfoTwo>& vec, float focus, quint32 martixSize, float pixelSize, QVector <quint32>& recCount, QVector<quint32>& times);
};
#endif // STARRECALGH_H
