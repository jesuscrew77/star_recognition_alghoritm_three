#include "starrecalgh.h"
#include <QDebug>


using namespace BOKZMath;
using namespace std;

QVector <StarInfoTwo> StarRecognizeThree::prepareStarInfoTwoData(float focus, float pixelSize)
{

    QVector <StarInfoTwo> starInfoVector;

    const auto countOfStars = catalog.alphaVec().size();

    for(int i = 0;i < countOfStars; i++)
    {
        double cos_b = cos(catalog.betaVec()[i] * transToRadians);
        double cos_a = cos(catalog.alphaVec()[i] * transToRadians);
        double sin_b = sin(catalog.betaVec()[i] * transToRadians);
        double sin_a = sin(catalog.alphaVec()[i] * transToRadians);
        lf.push_back(cos_b * cos_a);
        mf.push_back(cos_b * sin_a);
        nf.push_back(sin_b);
    }


    for (int cStar = 0; cStar < countOfStars; cStar ++) {
        StarInfoTwo starInfo;
        QVector <QPair <double, int>> distToStar;
        for( int j = 0; j < countOfStars; j ++) {
            auto scalar_product = calcScalarProduct(lf[cStar], lf[j], mf[cStar], mf[j] ,nf[cStar], nf[j]);
            distToStar.append(QPair <double, int> (acos(scalar_product), j + 1));
        }

        std::sort(distToStar.begin(), distToStar.end(), [](const auto& f, const auto& s) {return f.first < s.first;});
        starInfo.id = distToStar[0].second;

        int fID = 0;
        int sID = 0;
        if  (distToStar[1].first * transToGradus < 0.02)
        {
            starInfo.id1 = distToStar[2].second;
            starInfo.distance1 = distToStar[2].first * transToGradus;
            starInfo.id2 = distToStar[3].second;
            starInfo.distance2 = distToStar[3].first * transToGradus;
            fID = 2;
            sID = 3;

        }
        // ТУТ ИСПРАВИТЬ
        else if (abs(distToStar[1].first * transToGradus - distToStar[2].first * transToGradus) < 0.01666667)
        {

            starInfo.id1 = distToStar[2].second;
            starInfo.distance1 = distToStar[2].first * transToGradus;
            starInfo.id2 = distToStar[3].second;
            starInfo.distance2 = distToStar[3].first * transToGradus;
            fID = 2;
            sID = 3;
        }
        else
        {
            starInfo.id1 = distToStar[1].second;
            starInfo.distance1 = distToStar[1].first * transToGradus;
            starInfo.id2 = distToStar[2].second;
            starInfo.distance2 = distToStar[2].first * transToGradus;
            fID = 1;
            sID = 2;
        }


        auto trMatrix = calcTransitionMatrix(catalog.alphaVec()[cStar], catalog.betaVec()[cStar], 0);
        float CC;
        QVector <double> filtered_l_st {lf[distToStar[0].second - 1], lf[distToStar[fID].second - 1], lf[distToStar[sID].second - 1]};
        QVector <double> filtered_m_st {mf[distToStar[0].second - 1], mf[distToStar[fID].second - 1], mf[distToStar[sID].second - 1]};
        QVector <double> filtered_n_st {nf[distToStar[0].second - 1], nf[distToStar[fID].second - 1], nf[distToStar[sID].second - 1]};

        QVector <float> xCoords;
        QVector <float> yCoords;

        for (int i = 0; i < 3; i ++) {

            float x_coord;
            float y_coord;
            float x_coord_mm;
            float y_coord_mm;

            CC = trMatrix[2][0] * filtered_l_st[i]
                    +
                    trMatrix[2][1] * filtered_m_st[i]
                    +
                    trMatrix[2][2] * filtered_n_st[i];

            x_coord_mm = (- focus * (trMatrix[0][0] * filtered_l_st[i]
                    + trMatrix[0][1] * filtered_m_st[i]
                    + trMatrix[0][2] * filtered_n_st[i]) / CC);

            x_coord = x_coord_mm / pixelSize + 0.5;

            y_coord_mm = ( - focus * (trMatrix[1][0] * filtered_l_st[i]
                    + trMatrix[1][1] * filtered_m_st[i]
                    + trMatrix[1][2] * filtered_n_st[i]) / CC);

            y_coord = y_coord_mm / pixelSize + 0.5;

            double length = sqrt(x_coord_mm * x_coord_mm + y_coord_mm * y_coord_mm + focus * focus);

            filtered_l_st[i] = - x_coord_mm / length;
            filtered_m_st[i] = - y_coord_mm / length;
            filtered_n_st[i] =  focus / length;

            xCoords.append(x_coord);
            yCoords.append(y_coord);

            // x_coord += slideData.slideSizeX / 2; y_coord += slideData.slideSizeY / 2;// x_coord,slideSizeX- ось X, y_coord,slideSizeY - ось Y
        }

        QVector2D vecFirst (xCoords[1], yCoords[1]);
        QVector2D vecSecond (xCoords[2], yCoords[2]);
        vecFirst.normalize();
        vecSecond.normalize();
        starInfo.angle = acos(QVector2D::dotProduct(vecFirst, vecSecond)) * transToGradus;
        // if(flag) qDebug() << starInfo.angle;
        starInfoVector.append(std::move(starInfo));

    }
    std::sort(starInfoVector.begin(), starInfoVector.end(), [](const auto& f, const auto& s){return f.angle < s.angle;});

    for (auto it = starInfoVector.begin(); it != starInfoVector.end(); it ++)
    {
        auto itf = starInfoVector.begin();
        while (it->id1 != itf->id) itf++;
        it->id1 = std::distance(starInfoVector.begin(), itf);
        itf = starInfoVector.begin();
        while (it->id2 != itf->id) itf++;
        it->id2 = std::distance(starInfoVector.begin(), itf);
    }


    QFile text ("star_info_text.txt");
    if (text.open(QFile::WriteOnly | QFile::Truncate)) {
        QTextStream outText(&text);
        for (const auto& i : starInfoVector) {
            QString row = QString("%1%2%3%4%5%6\n")
                    .arg(i.angle, 12,'g', 7)
                    .arg(i.id, 12)
                    .arg(i.id1, 12)
                    .arg(i.distance1, 12 ,'g', 7)
                    .arg(i.id2, 12)
                    .arg(i.distance2, 12, 'g', 7);
            outText << row;

        }
        text.close();
    }
    // qDebug() << "size" << sizeof(StarInfoTwo) << starInfoVector.size();
    QFile bin ("star_info_bin.bin");
    if (bin.open(QFile::WriteOnly | QFile::Truncate)) {
        QDataStream outBin(&bin);
        for (const auto& i : starInfoVector) {
            outBin.writeRawData((const char*)&i, sizeof(StarInfoTwo));
        }
        bin.close();
    }
    return starInfoVector;
}


bool StarRecognizeThree::checkAngleInRange (double angle, double compAngle, double range)
{
    return (angle >= compAngle - range) && (angle <= compAngle + range);
}


// nPos - если 1, то 1 и 2 сосед , если 2, то 2 и 1.
bool StarRecognizeThree::recognizeStar(const QVector <QPointF>& starsCoords, const QVector <StarInfoTwo>& vec, RecognizedInfo** recMatrix, quint32& i , quint32 f, quint32 s, quint32& starNumber, quint32& nPos)
{

    QVector2D vecFirst (starsCoords[f].x() - starsCoords[i].x()
            , starsCoords[f].y() - starsCoords[i].y());
    vecFirst.normalize();

    QVector2D vecSecond (starsCoords[s].x() - starsCoords[i].x()
            , starsCoords[s].y() - starsCoords[i].y());
    vecSecond.normalize();

    double angleToFind = acos(QVector2D::dotProduct(vecFirst, vecSecond)) * transToGradus;
    constexpr const double error = 0.8;

    auto it = std::lower_bound(vec.begin(), vec.end(), angleToFind - error, [](const auto& fs, const auto& ss){return fs.angle < ss;});
    float fAngle;
    if (f > i) fAngle = recMatrix[f][i].distance;
    else fAngle = recMatrix[i][f].distance;
    //qDebug() << "Углы";
    //qDebug() << fAngle;

    float sAngle;
    if (s > i) sAngle = recMatrix[s][i].distance;
    else sAngle = recMatrix[i][s].distance;
    //qDebug() << sAngle;

    while (it->angle < angleToFind + error && it->angle < 180) {

        if (checkAngleInRange(it->distance1, fAngle, 0.0167) && checkAngleInRange(it->distance2, sAngle, 0.0167)) {
            nPos = 1;
            starNumber = std::distance(vec.begin(), it);
            return true;
        }

        if (checkAngleInRange(it->distance1, sAngle, 0.0167) && checkAngleInRange(it->distance2, fAngle, 0.0167)) {
            nPos = 2;
            starNumber = std::distance(vec.begin(), it);
            return true;
        }

        else {
            ++it;
        }
    }
    return false;
}




void StarRecognizeThree::testRecognitionCalibrationFile(const QString& fileName, QVector <StarInfoTwo>& vec, float focus, quint32 martixSize, float pixelSize, QVector <quint32>& recCount, QVector<quint32>& times)
{
    focus = focus / pixelSize;
    QFile file (fileName);
    if (file.open(QFile::ReadOnly)) {
        QTextStream in(&file);
        QString currentLine;
        constexpr const quint32 countOfStars = 20;
        QVector <QPointF> starsCoords;
        double alpha = 0;
        double delta = 0;
        double azimut = 0;
        bool coordsReaded = false;

#ifdef _TEST_
        QVector <int> withoutNoise;
        QVector <int> withNoise;
#endif

        while  (in.readLineInto(&currentLine)) {
            if (currentLine.contains("№      Br       X        Y       Ne     Max")) {
                for (int i = 0; i < countOfStars; i ++) {
                    in.readLineInto(&currentLine);
                    QStringList list = currentLine.split("  ", QString::SkipEmptyParts);
                    starsCoords.append(QPointF(list[2].toFloat() - martixSize / 2, list[3].toFloat() - martixSize / 2));
                }
                coordsReaded = true;
            }
            if (currentLine.contains("Матрица ориентации:") && coordsReaded) {
                double orientMatrix[3][3];
                for (int i = 0; i < 3; i ++) {
                    for (int j = 0; j < 3; j ++) {
                        in >> orientMatrix[i][j];
                    }
                }
                matrToEkvA(orientMatrix, alpha, delta, azimut);
                break;
            }
        }

        if  (!coordsReaded || !alpha) {
            qDebug() << starsCoords << " - неполный протокол";
        }

        qDebug() << "Следующий кадр";
#ifdef _TEST_
        withoutNoise.clear();

        for (int c = 0; c < 20; c++) {
            withNoise.clear();
            if (c > 0) {
                starsCoords.erase(starsCoords.begin() + countOfStars, starsCoords.end());
                for (int i = 0; i < 45; i++) {
                    starsCoords.append(QPointF(static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 2048)), static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 2048))));
                }
            }
#endif
            QVector <float> l;
            QVector <float> m;
            QVector <float> n;

            QElapsedTimer timer;
            timer.start();


            for (int i = 0; i < starsCoords.size(); i ++)
            {
                double length = sqrt(starsCoords[i].x() * starsCoords[i].x() + starsCoords[i].y() * starsCoords[i].y() + focus * focus);
                l.append(- starsCoords[i].x() / length);
                m.append(- starsCoords[i].y() / length);
                n.append(focus / length);
            }


            int objectsCount = starsCoords.size();
            RecognizedInfo** recMatrix = new RecognizedInfo* [objectsCount];
            for (int i = 0; i < starsCoords.size(); i++) {
                recMatrix[i] = new RecognizedInfo [objectsCount];
            }

            int tempObjectsCount = objectsCount - 2;
            for (int i = objectsCount - 1; i >= 0; i--) {
                for (int j = tempObjectsCount; j >= 0; j--) {
                    recMatrix[i][j].distance = acos(calcScalarProduct(l[i], l[j], m[i], m[j], n[i], n[j])) * transToGradus;
                }
                --tempObjectsCount;
            }

#ifdef _TEST_
            for (int i = 0; i < objectsCount; i++) {
                QString row;
                for (int j = 0; j < objectsCount; j ++) {
                    row.append(QString::number(recMatrix[i][j].distance) + " ");

                }
                qDebug() << row;
            }
#endif

            for (quint32 i = 0; i < starsCoords.size(); i ++) {
                auto mins = BOKZMath::firstMinDistanceTable(recMatrix, 4, i, objectsCount);
                quint32 starNumber = 0;
                quint32 nPos = 0;
                for (int fn = 0; fn < mins.size() - 1 ; fn++) {
                    for (int sn = fn + 1; sn < mins.size(); sn++) {
                        if (recognizeStar(starsCoords, vec, recMatrix, i, mins[fn], mins[sn], starNumber, nPos)) {
                            recMatrix[i][mins[fn]].checked = true;
                            recMatrix[i][mins[sn]].checked = true;
                            recMatrix[mins[fn]][i].checked = true;
                            recMatrix[mins[sn]][i].checked = true;

                            recMatrix[i][i].recognized = true;
                            recMatrix[mins[sn]][mins[sn]].recognized = true;
                            recMatrix[mins[fn]][mins[fn]].recognized = true;


                            recMatrix[i][i].numInCat = vec[starNumber].id - 1;
                            if (nPos == 1) {
                                recMatrix[mins[fn]][mins[fn]].numInCat = vec[vec[starNumber].id1].id - 1;
                                recMatrix[mins[sn]][mins[sn]].numInCat = vec[vec[starNumber].id2].id - 1;

                            }
                            else if (nPos == 2) {
                                recMatrix[mins[fn]][mins[fn]].numInCat = vec[vec[starNumber].id2].id - 1;
                                recMatrix[mins[sn]][mins[sn]].numInCat = vec[vec[starNumber].id1].id - 1;

                            }

                        }
                    }
                }
            }

            quint32 recognizedCount = 0;
            for (int i = 0; i < objectsCount; i++) {
                quint32 checked = 0;
                quint32 counter = 0;
                if (!recMatrix[i][i].recognized) continue;
                for (int j = 0; j < objectsCount; j++) {
                    if (j == i || !recMatrix[j][j].recognized) continue;
                    if (recMatrix[i][j].checked) {
                        checked++;
                        continue;
                    }
                    double cos_bs = cos(catalog.betaVec()[recMatrix[j][j].numInCat] * transToRadians);
                    double cos_as = cos(catalog.alphaVec()[recMatrix[j][j].numInCat] * transToRadians);
                    double sin_bs = sin(catalog.betaVec()[recMatrix[j][j].numInCat] * transToRadians);
                    double sin_as = sin(catalog.alphaVec()[recMatrix[j][j].numInCat] * transToRadians);
                    double ls = cos_bs * cos_as;
                    double ms = cos_bs * sin_as;
                    double ns = sin_bs;
                    double angle = acos(
                                calcScalarProduct
                                (lf[recMatrix[i][i].numInCat], ls, mf[recMatrix[i][i].numInCat], ms, nf[recMatrix[i][i].numInCat], ns)) * transToGradus;

                    double rDistance = 0;
                    if (recMatrix[j][i].distance == 0) {
                        rDistance = recMatrix[i][j].distance;
                    }
                    else {
                        rDistance = recMatrix[j][i].distance;
                    }
                    double diff = rDistance - angle;
                    if (abs(diff) < 0.0187) {
                        ++counter;
                    }


                }
               // qDebug() << counter << checked;
                if(counter >= 3  || checked >= 5) ++recognizedCount;
            }
            //qDebug() <<recognizedCount << timer.nsecsElapsed();
            times.append(timer.nsecsElapsed());
            recCount.append(recognizedCount);
            //  qDebug() << recognizedCount;
            //        std::sort(recognizedNumbers.begin(), recognizedNumbers.end());
            //        recognizedNumbers.erase(std::unique(recognizedNumbers.begin(),recognizedNumbers.end()), recognizedNumbers.end());
            //qDebug() <<  recognizedCount <<"\t" << timer.nsecsElapsed();
            //if (c == 0) withoutNoise.append(recognizedNumbers.size());
            // else withNoise.append(recognizedNumbers.size());

            // if (c > 0 && withNoise.last() - withoutNoise.last()  != 0) {
            //   counter++;
            //     qDebug() << "Влияние помех" << withNoise.last() - withoutNoise.last();
            // }

#ifdef _TEST_
        }
#endif
    }
    file.close();
}

void StarRecognizeThree::makeTest()
{
    QVector <StarInfoTwo> vec = prepareStarInfoTwoData(31.8894, 0.0055);
    QDir testDir("C:/Users/Public/Documents/catalog/test");
    auto files = testDir.entryList(QDir::Files);
    QVector <quint32> recCount;
    QVector <quint32> times;
    for (const auto& file : files) {
        testRecognitionCalibrationFile(testDir.path() + "/" + file, vec, 31.8894, 2048, 0.0055, recCount, times);
    }
    std::sort(recCount.begin(), recCount.end());
    auto meanstd = calculateMeanStDv(recCount.begin(), recCount.end(), 0);
    qDebug() << meanstd.first;
    auto min = min_element(recCount.begin(), recCount.end());
    qDebug() << *min;
    auto max = max_element(recCount.begin(), recCount.end());
    qDebug() << *max;
    if(*min == 0)
    {
        auto zerosCount = std::count(recCount.begin(), recCount.end(), 0);
        qDebug() << "Число нулей" <<zerosCount;
    }
    std::remove(recCount.begin(), recCount.end(), 0);
    auto minWithoutZero = min_element(recCount.begin(), recCount.end());
    qDebug() << *minWithoutZero;

    std::sort(times.begin(), times.end());
    auto meanstdt = calculateMeanStDv(times.begin(), times.end(), 0);
    qDebug() << meanstdt.first;
    auto mint = min_element(times.begin(), times.end());
    qDebug() << *mint;
    auto maxt = max_element(times.begin(), times.end());
    qDebug() << *maxt;
}
