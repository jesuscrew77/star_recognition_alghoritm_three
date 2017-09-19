#include <mathfunc.h>
#include <qDebug>


namespace BOKZMath {


double calcScalarProduct(double l_oz, double l_st, double m_oz, double m_st, double n_oz, double n_st)
{
    double scalar_product = l_oz * l_st + m_oz * m_st + n_oz * n_st;

    if (scalar_product > 1.0)// проверка на выход значения косинуса за диапазон [-1,1]
    {
        scalar_product = 1;
    }
    else if (scalar_product < -1.0)
    {
        scalar_product = -1;
    }
    return scalar_product;
}


void calculateSunDirection(double JD_DEV, double *pSunI) noexcept
{
    double Ts, Ls, Ms, eS, Es, Vs, Qs, Ekl, AlphaS, DeltaS, QsVid, Om;

    Ts = JD_DEV+1;
    Ls = 279.69668+36000.76892*Ts+0.0003025*Ts*Ts;
    Ls /= 360.;
    Ls -= (int)Ls;
    Ls *= 360.;

    Ms = 358.47583+35999.04975*Ts-0.000150*Ts*Ts-0.0000033*Ts*Ts*Ts;
    Ms /= 360.;
    Ms -= (int)Ms;
    Ms *= 360.;

    Ms *= M_PI/180.;

    eS = 0.01675104-0.0000418*Ts-0.000000126 * Ts * Ts;
    Es = Ms;
    for (int i = 0; i < 10; i++) Es = Ms+eS*sin(Es);
    Vs=2*atan(sqrt((1+eS)/(1-eS))*tan(Es/2));
    Qs=Ls*M_PI/180.+Vs-Ms;

    Om = 259.18-1934.142*Ts;
    Om /= 360.;
    Om -= (int)Om;
    Om *= 2*M_PI;

    QsVid = Qs*180/M_PI-0.00569-0.00479*sin(Om);
    QsVid = QsVid*M_PI/180.;

    Ekl = 23.452294-0.0130125*Ts-0.00000164*Ts*Ts+0.000000503*Ts*Ts*Ts;
    Ekl += 0.00256*cos(Om);
    Ekl *= M_PI/180.;

    AlphaS = atan2(cos(Ekl)*sin(QsVid),cos(QsVid));
    DeltaS = asin(sin(Ekl)*sin(QsVid));

    if (AlphaS < 0) AlphaS += 2*M_PI;

    pSunI[0] = cos(AlphaS)*cos(DeltaS);
    pSunI[1] = sin(AlphaS)*cos(DeltaS);
    pSunI[2] = sin(DeltaS);

}

void calculateMoonDirection(double JD_DEV, double *pMoonI) noexcept
{
    double Ts, L1, Ms, M1, D, F, Om, e, Lam, B, Ekl, AlphaM, DeltaM;

    Ts=JD_DEV+1;
    //знак при Ts^2 должен быть минус, а не плюс
    L1=270.434164+481267.8831*Ts-0.001133*Ts*Ts+0.0000019*Ts*Ts*Ts;
    L1/=360.; L1-=(int)L1; L1*=360.;

    Ms=358.47583+35999.04975*Ts-0.000150*Ts*Ts-0.0000033*Ts*Ts*Ts;
    Ms/=360.; Ms-=(int)Ms; Ms*=2*M_PI;

    M1=296.104608+477198.8492*Ts+0.009192*Ts*Ts+0.0000144*Ts*Ts*Ts;
    M1/=360.; M1-=(int)M1; M1*=2*M_PI;

    D=350.737486+445267.1142*Ts-0.001436*Ts*Ts+0.0000019*Ts*Ts*Ts;
    D/=360.; D-=(int)D; D*=2*M_PI;

    F=11.250889+483202.0251*Ts-0.003211*Ts*Ts-0.000003*Ts*Ts*Ts;
    F/=360.; F-=(int)F; F*=2*M_PI;

    Om=259.183275-1934.1420*Ts+0.002078*Ts*Ts+0.0000022*Ts*Ts*Ts;
    Om/=360.; Om-=(int)Om; Om*=2*M_PI;

    e=1-0.002495*Ts-0.00000752*Ts*Ts;

    Lam=L1+6.288750*sin(M1)+1.274018*sin(2*D-M1)+0.658309*sin(2*D)
            +0.213616*sin(2*M1)-e*0.185596*sin(Ms)-0.114336*sin(2*F)
            +0.058793*sin(2*D-2*M1)+e*0.057212*sin(2*D-Ms-M1)+0.053320*sin(2*D+M1)
            +e*0.045874*sin(2*D-Ms)+e*0.041024*sin(M1-Ms)-0.034718*sin(D)-e*0.030465*sin(Ms+M1)
            +0.015326*sin(2*D-2*F)-0.012528*sin(2*F+M1)-0.01098*sin(2*F-M1)+0.010674*sin(4*D-M1)+0.010034*sin(3*M1);

    Lam/=360; Lam-=(int)Lam; Lam*=2*M_PI;

    B=5.128189*sin(F)+0.280606*sin(M1+F)+0.277693*sin(M1-F)+0.173238*sin(2*D-F)+0.055413*sin(2*D+F-M1)+0.046272*sin(2*D-F-M1)
            +0.032573*sin(2*D+F)+0.017198*sin(2*M1+F)+0.009267*sin(2*D+M1-F)+0.008823*sin(2*M1-F);
    B/=360; B-=(int)B; B*=2*M_PI;

    Ekl=23.452294-0.0130125*Ts-0.00000164*Ts*Ts+0.000000503*Ts*Ts*Ts;
    Ekl*=M_PI/180.;

    AlphaM=atan2((sin(Lam)*cos(Ekl)-tan(B)*sin(Ekl)),cos(Lam));
    DeltaM= asin(sin(B)*cos(Ekl)+cos(B)*sin(Ekl)*sin(Lam));

    pMoonI[0] = cos(AlphaM) * cos(DeltaM);
    pMoonI[1] = sin(AlphaM) * cos(DeltaM);
    pMoonI[2] = sin(DeltaM);
}




void multMatrixVector(const double M[3][3], const double V[3], double VR[3]) noexcept
{
    VR[0] = M[0][0]*V[0]+M[0][1]*V[1]+M[0][2]*V[2];
    VR[1] = M[1][0]*V[0]+M[1][1]*V[1]+M[1][2]*V[2];
    VR[2] = M[2][0]*V[0]+M[2][1]*V[1]+M[2][2]*V[2];
}


/*меняем местами 1,2 и 3,4 байты qint32 местами*/
quint32 correctDTMITimePr(quint32 incrTimePr) noexcept
{
    quint16 parts2Bytes[2];
    parts2Bytes[0] = *((quint16*)&incrTimePr);
    parts2Bytes[1] = *((quint16*)&incrTimePr+1);
    *((quint16*)&incrTimePr) = parts2Bytes[1];
    *((quint16*)&incrTimePr+1) = parts2Bytes[0];

    return incrTimePr;
}





/*меняем местами 1,2 и 3,4 байты float местами*/
float correctDTMIFloat(float element) noexcept
{
    quint16 parts2Bytes[2];
    parts2Bytes[0] = *((quint16*)&element);
    parts2Bytes[1] = *((quint16*)&element+1);
    *((quint16*)&element) = parts2Bytes[1];
    *((quint16*)&element+1) = parts2Bytes[0];
    return element;
}





/*корректируем весь массив локализованных объектов*/
void correctLocArray(float locArray[16][4]) noexcept
{
    quint16 locColumns = 16;
    quint16 locRows = 4;

    /* конвертируем Loc в формат postgres массива (хранить будем его как одномерный, иначе долго)*/
    for(qint32 locColumn = 0;locColumn < locColumns;locColumn ++)
    {
        for(qint32 locRow=0;locRow<locRows;locRow++)
        {
            locArray[locColumn][locRow] = correctDTMIFloat(locArray[locColumn][locRow]);
        }
    }
}


void correctAngularVArray(float Wo[3]) noexcept
{
    quint16 WoLenght = 3;
    for(qint32 i = 0;i < WoLenght;i ++)
    {
        Wo[i] = correctDTMIFloat(Wo[i]);
    }

}


// расчёт юлианской даты
double calculateJulianDate(const QDateTime & dateTime) noexcept
{
    double julianDate =
            static_cast <double> (dateTime.date().toJulianDay()) - 0.5;
    julianDate += static_cast <double> (dateTime.time().hour()) / 24
            + static_cast <double> (dateTime.time().minute()) / 1440
            + static_cast <double> (dateTime.time().second()) / 86400;
    return julianDate;
}

double calculateJulianCenturies(const QDateTime& dateTime) noexcept
{
    double jd = calculateJulianDate(dateTime);
    return (jd - BOKZMath::julianZeroDate) / BOKZMath::countOfDaysInCentury;
}


inline double mSecsToSecs(quint64 mSecsSinceEpoch) noexcept
{
    return (double)mSecsSinceEpoch / 1000;
}


// для вычисления дат времени привязки (юлианская дата на полночь)
double calculateDatePr(const QDateTime& dateTime, int timePr) noexcept
{
    int julianDate  = dateTime.addSecs(-timePr).date().toJulianDay();
    return static_cast <double> (julianDate) - 0.5;
}


// вычисляет дату как TDateTime::value, НЕ УЧТЕНЫ миллисекунды.
double toTDateTimeValue(const QDateTime& datetime) noexcept
{
    auto daysToTDateTimeValue = abs(datetime.daysTo(QDateTime(QDate(1899,12,30))));
    double partOfDay = static_cast <double> (abs(datetime.time().secsTo(QTime(0,0,0))))/86400;
    return static_cast <double> (daysToTDateTimeValue) + partOfDay;
}

QDateTime fromTDateTime(double value) noexcept
{
    return QDateTime(QDate(1899,12,30)).addDays(static_cast<quint32>(value)).addSecs(round((value - static_cast<quint32>(value)) * 86400));
}


/*смотри http://law.rufox.ru/view/standarti/1672.htm*/
double starsTime2(double JD) noexcept
{
    double     l1;      // средняя аномалия Солнца
    double     F;       // средний аргумент широты Луны
    double     D;       // разность средних долгот Луны и Солнца
    double     W;       // средняя долгота восходящего узла орбиты Луны на эклиптике
    double     Nf;      // значения нутации в долготе
    double     Ne;      //значения нутации в наклоне
    double     e0;      // средний наклон эклиптики к экватору, выраженный в радианах
    double     M;       //всемирное время  рассматриваемой даты, выраженное в долях суток
    double     tau;     //интервал времени от эпохи T0  до эпохи t в юлианских столетиях по 36525 средних солнечных суток
    double     tau_2;   // тоже самое в квадрате
    double     tau_3;   // тоже самоее в кубе
    double     Ssr_rad; // гринвичское среднее звездное время


    tau=JD/36525.;
    tau_2=tau*tau;
    tau_3=tau_2*tau;

    M=JD+0.5-(qint32)JD;
    Ssr_rad=1.7533685592+0.0172027918051*JD+6.2831853072*M+6.7707139e-6*tau_2
            -4.50876e-10*tau_2*tau;

    l1=6.24003594 + 628.30195602 *tau - 2.7974e-6 *tau_2-5.82e-8*tau_3;
    F =1.62790193 + 8433.46615831 *tau - 6.42717e-5*tau_2 + 5.33e-8*tau_3;
    D =5.19846951 + 7771.37714617 *tau - 3.34085e-5*tau_2 + 9.21e-8*tau_3;
    W =2.182438624 - 33.757045936*tau + 3.61429e-5*tau_2 + 3.88e-8*tau_3;
    Nf=-0.83386e-4*sin(W)+0.9997e-6*sin(2*W)-0.63932e-5*sin(2*(F-D+W))+
            0.6913e-6*sin(l1)-0.11024e-5*sin(2*(F+W));
    Ne=0.44615e-4*cos(W)+0.27809e-5*cos(2*(F-D+W))+0.474e-6*cos(2*(F+W));
    e0=0.4090928042-0.2269655e-3*tau-0.29e-8*tau_2+0.88e-8*tau_3;

    //учитываем нутацию в прямом восхождении
    Ssr_rad+=Nf*cos(e0+Ne);

    //приводим звездное время к "нормальному" виду
    while (Ssr_rad>2*M_PI) Ssr_rad=Ssr_rad-2*(M_PI);

    return Ssr_rad;
}


void getGSKISKMatrix(double JD, double MG[3][3]) noexcept
{
    double ST = starsTime2(JD); // истинное звездное время

    MG[0][0]=cos(ST);
    MG[0][1]=-sin(ST);
    MG[0][2]=0;

    MG[1][0]=sin(ST);
    MG[1][1]=cos(ST);
    MG[1][2]=0;

    MG[2][0]=0;
    MG[2][1]=0;
    MG[2][2]=1;
}


void  transSpeedGSKtoISK(double JD, double CG[3], double VG[3], double VI[3]) noexcept
{
    double MG[3][3]; // матрица перехода от ГСК к ИСК
    getGSKISKMatrix(JD, MG);

    constexpr const double radianToSec=180*3600/M_PI;
    constexpr const double W=15.041/radianToSec;

    double CI [3];
    multMatrixVector(MG, CG, CI);
    multMatrixVector(MG, VG, VI);

    for (qint32 i=0; i<3; i++)
    {
        CI[i]/=1000.; // в километры
        VI[i]/=1000.;
    }

    VI[0]-=W*CI[1];  VI[1]+=W*CI[0];
}

void  transCoordsGSKtoISK(double JD, double CG[3], double CI[3]) noexcept
{

    double MG[3][3]; // матрица перехода от ГСК к ИСК
    getGSKISKMatrix(JD,MG);
    multMatrixVector(MG, CG, CI);

    for (qint32 i=0; i<3; i++)
    {
        CI[i]/=1000.; // в километры
    }

}


void quatToMatr(const double Quat[], double M_ornt[3][3]) noexcept
{

    M_ornt[0][0]=Quat[0]*Quat[0]+Quat[1]*Quat[1]-Quat[2]*Quat[2]-Quat[3]*Quat[3];
    M_ornt[0][1]=2.0*(Quat[1]*Quat[2]+Quat[0]*Quat[3]);
    M_ornt[0][2]=2.0*(Quat[1]*Quat[3]-Quat[0]*Quat[2]);
    M_ornt[1][0]=2.0*(Quat[1]*Quat[2]-Quat[0]*Quat[3]);
    M_ornt[1][1]=Quat[0]*Quat[0]-Quat[1]*Quat[1]+Quat[2]*Quat[2]-Quat[3]*Quat[3];
    M_ornt[1][2]=2.0*(Quat[2]*Quat[3]+Quat[0]*Quat[1]);
    M_ornt[2][0]=2.0*(Quat[1]*Quat[3]+Quat[0]*Quat[2]);
    M_ornt[2][1]=2.0*(Quat[2]*Quat[3]-Quat[0]*Quat[1]);
    M_ornt[2][2]=Quat[0]*Quat[0]-Quat[1]*Quat[1]-Quat[2]*Quat[2]+Quat[3]*Quat[3];
}


double atan2m(double yf, double xf) noexcept
{
    double ang;
    if (fabs(xf)>1e-10)
    {
        ang=atan2(yf,xf);
    }
    else
    {
        if (yf>0)
        {
            ang=M_PI/2.;
        }
        else
        {
            ang=-M_PI/2.;
        }
    }
    return ang;
}

void matrToEkvA(double M_ornt[3][3], double& al, double& dl, double& Az) noexcept
{
    dl=asin(M_ornt[2][2]);
    al=atan2m(M_ornt[2][1],M_ornt[2][0]);
    if (al<0)
    {
        al+=2*M_PI;
    }
    Az=atan2m(M_ornt[0][2],M_ornt[1][2]);
    if (Az<0)
    {
        Az+=2*M_PI;
    }
}

void quatToEkvA(double Quat[], double EkvA[3]) noexcept
{
    double Matr[3][3],Al, Dl, Az;
    quatToMatr(Quat,Matr);
    matrToEkvA(Matr,Al,Dl,Az);

    constexpr double trans_to_grad=57.29577957855229;
    EkvA[0]=Al*trans_to_grad; EkvA[1]=Dl*trans_to_grad; EkvA[2]=Az*trans_to_grad;
}


double acosm(double xf) noexcept
{
    if (xf>1.) xf=1.;
    else if (xf<-1.) xf=-1.;

    return acos(xf);
}

void multMatrix(const double Matr1[3][3],const double Matr2[3][3], double Matr[3][3]) noexcept
{
    double buf;
    int i,j,k;

    for(i = 0;i < 3;i ++)
    {
        for(j = 0;j < 3;j ++)
        {
            buf = 0;
            for (k = 0;k < 3;k ++)
            {
                buf = buf + Matr1[i][k] * Matr2[k][j];
            }
            Matr[i][j] = buf;
        }
    }
}

void getAngularDisplacementFromOrientMatr(const double M_ornt_pr[3][3],const double M_ornt[3][3], double Wop [3]) noexcept

{
    double MT_ornt_pr[3][3], dMB[3][3];
    double delta, sdt;
    // считаем матрицу, дающую противоположный поворот
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            MT_ornt_pr[i][j] = M_ornt_pr[j][i];
        }
    }
    // перемножаем, чтобы получить разницу
    multMatrix(M_ornt,MT_ornt_pr,dMB);


    delta=(dMB[0][0] + dMB[1][1] + dMB[2][2] - 1.)/2.;

    // если угол очень маленький, то считаем дельта равной 0
    if (delta>1) delta=1;

    if (fabs(delta-1.0)< 1E-20)
    {
        Wop[0]=0;
        Wop[1]=0;
        Wop[2]=0;
    }
    else
    {
        delta = acosm(delta);
        sdt = sin(delta);
        Wop[0]= -delta*(dMB[2][1] - dMB[1][2]) / (2.0 * sdt);
        Wop[1]= -delta*(dMB[0][2] - dMB[2][0]) / (2.0 * sdt);
        Wop[2]= -delta*(dMB[1][0] - dMB[0][1]) / (2.0 * sdt);
    }
}

double calculateAngleAxis(const double quat1[], const double quat2[], AXIS_TYPE axis) noexcept
{
    double orientMatrix1[3][3];
    BOKZMath::quatToMatr(quat1,orientMatrix1);

    /*считаем углы между заданной осью*/
    qint32 axisBetween = static_cast <qint32> (axis);

    double lfirst = orientMatrix1[axisBetween][0];
    double mfirst = orientMatrix1[axisBetween][1];
    double nfirst = orientMatrix1[axisBetween][2];

    double orientMatrix2[3][3];
    BOKZMath::quatToMatr(quat2,orientMatrix2);

    double lsecond = orientMatrix2[axisBetween][0];
    double msecond = orientMatrix2[axisBetween][1];
    double nsecond = orientMatrix2[axisBetween][2];


    return std::acos(lfirst * lsecond + mfirst * msecond + nfirst * nsecond) * BOKZMath::transToGradus;
}


QVector <QPointF> createHistogramm(size_t histSize, float shiftRange, QVector <float> firstX, QVector <float> firstY, QVector <float> secondX, QVector <float> secondY)
{
    QVector <QVector <quint32>> histogramm;
    histogramm.resize(histSize);
    for(auto& i : histogramm)
    {
        i.resize(histSize);
    }
    quint32 centerX, centerY, range;
    centerX = centerY = (histSize / 2) + 1;
    range = histSize / 2;

    quint32 countObjectsFirst = firstX.size();
    quint32 countObjectsSecond = secondX.size();

    qint32 dx, dy;
    for(quint32 i = 0; i < countObjectsFirst; i ++)
    {
        if (firstX[i] == 0 && firstY[i] == 0)
        {
            break;
        }
        for (quint32 j = 0; j < countObjectsSecond; j ++)
        {
            dx = round(firstX[i] - secondX[j]);
            dy = round(firstY[i] - secondY[j]);

            if (abs(dx) < range && abs(dy) < range )
            {
                ++histogramm[centerX + dx][centerY + dy];
            }
        }
    }
    auto maxHist = matrixMax(histogramm);
    qint32 shiftX = maxHist.first - centerX;
    qint32 shiftY = maxHist.second - centerY;

    QVector <QPointF> shifts;
    float currShiftX = 0.0;
    float currShiftY = 0.0;

    // считаем смещения для всех совпадающих по порядку объектов
    for(quint32 i = 0; i < countObjectsFirst; i ++)
    {
        currShiftX = firstX[i] - secondX[i];
        currShiftY = firstY[i] - secondY[i];

        if (currShiftX <= (shiftX + shiftRange) && currShiftX >= (shiftX - shiftRange)
                && currShiftY <= (shiftY + shiftRange) && currShiftY >= (shiftY - shiftRange))
        {
            shifts.append(QPointF(currShiftX,currShiftY));
        }
        else
        {
            shifts.append(QPointF(0.0,0.0));
        }
    }

    // считаем смещения для тех объектов, которые могут находиться в разных порядках
    // а возможно он уже просто пропал, тогда ничего не изменится
    for(quint32 i = 0; i < shifts.size(); i ++)
    {
        if (shifts[i].x() == 0 && shifts[i].y() == 0)
        {
            for(quint32 j = 0; j < countObjectsSecond; j++)
            {
                currShiftX = firstX[i] - secondX[j];
                currShiftY = firstY[i] - secondY[j];

                if (currShiftX <= (shiftX + shiftRange) && currShiftX >= (shiftX - shiftRange)
                        && currShiftY <= (shiftY + shiftRange) && currShiftY >= (shiftY - shiftRange))
                {
                    shifts[i] = QPointF(currShiftX,currShiftY);
                }
            }
        }
    }
    return shifts;
}

// потом отрефакторить
QVector<QPointF> createHistogramm(size_t histSize, float shiftRange, const QVector<QPointF>& firstCoords, const QVector<QPointF>& secondCoords)
{
    QVector <QVector <quint32>> histogramm;
    histogramm.resize(histSize);
    for(auto& i : histogramm)
    {
        i.resize(histSize);
    }
    quint32 centerX, centerY, range;
    centerX = centerY = (histSize / 2) + 1;
    range = histSize / 2;

    quint32 countObjectsFirst = firstCoords.size();
    quint32 countObjectsSecond = secondCoords.size();
    qint32 dx, dy;
    for(quint32 i = 0; i < countObjectsFirst; i ++)
    {
        if (firstCoords[i].x() == 0 && firstCoords[i].y() == 0)
        {
            break;
        }
        for (quint32 j = 0; j < countObjectsSecond; j ++)
        {
            dx = round(firstCoords[i].x() - secondCoords[j].x());
            dy = round(firstCoords[i].y() - secondCoords[j].y());

            if (abs(dx) < range && abs(dy) < range)
            {
                ++histogramm[centerX + dx][centerY + dy];
            }
        }
    }
    auto maxHist = matrixMax(histogramm);
    qint32 shiftX = maxHist.first - centerX;
    qint32 shiftY = maxHist.second - centerY;

    QVector <QPointF> shifts;
    float currShiftX = 0.0;
    float currShiftY = 0.0;

    // считаем смещения для всех совпадающих по порядку объектов
    for(quint32 i = 0; i < countObjectsFirst; i ++)
    {
        currShiftX = firstCoords[i].x() - secondCoords[i].x();
        currShiftY = firstCoords[i].y() - secondCoords[i].y();

        if (currShiftX <= (shiftX + shiftRange) && currShiftX >= (shiftX - shiftRange)
                && currShiftY <= (shiftY + shiftRange) && currShiftY >= (shiftY - shiftRange))
        {
            shifts.append(QPointF(currShiftX,currShiftY));
        }
        else
        {
            shifts.append(QPointF(0.0,0.0));
        }
    }
    // считаем смещения для тех объектов, которые могут находиться в разных порядках
    // а возможно он уже просто пропал, тогда ничего не изменится
    for(int i = 0; i < shifts.size(); i ++)
    {
        if (shifts[i].x() == 0 && shifts[i].y() == 0)
        {
            for(quint32 j = 0; j < countObjectsSecond; j++)
            {
                currShiftX = firstCoords[i].x() - secondCoords[j].x();
                currShiftY = firstCoords[i].y() - secondCoords[j].y();

                if (currShiftX <= (shiftX + shiftRange) && currShiftX >= (shiftX - shiftRange)
                        && currShiftY <= (shiftY + shiftRange) && currShiftY >= (shiftY - shiftRange))
                {
                    shifts[i] = QPointF(currShiftX,currShiftY);
                }
            }
        }
    }
    return shifts;
}

qint64 roundUp(qint64 numToRound, qint64 multiple)
{
    if (multiple == 0)
        return numToRound;

    auto remainder = std::remainder(numToRound, multiple);
    if (remainder == 0)
        return numToRound;

    if (numToRound < 0)
        return -(abs(numToRound) - remainder);
    else
        return numToRound + multiple - remainder;
}

QDateTime timePrToDateTime(const double timePr, const QDateTime& timePrStartData) noexcept
{
    qint64 dtInMsecs = timePr * BOKZMath::MSecsInSec;
    return timePrStartData.addMSecs(dtInMsecs);
}

QDateTime timePrToDateTime(const quint32 timePr, const QDateTime& timePrStartData) noexcept
{
    qint64 dtInMsecs = static_cast<qint64> (timePr) * BOKZMath::MSecsInSec;
    return timePrStartData.addMSecs(dtInMsecs);
}

QDateTime timePrToDateTime(const quint64 timePr, const QDateTime& timePrStartData) noexcept
{
    return timePrStartData.addMSecs(timePr * BOKZMath::MSecsInSec);
}

QVector< QVector<float> > calcTransitionMatrix(double pointAlpha, double pointBeta, double pointAzimut) noexcept
{

    QVector< QVector<float> > trMat(3);
    float PS,PC,QS,QC,RS,RC;
    for(int i = 0;i < trMat.size();i ++) trMat[i].resize(3);

    PS = sin(pointAzimut * BOKZMath::transToRadians); PC= cos(pointAzimut * BOKZMath::transToRadians);
    QS = sin(pointBeta * BOKZMath::transToRadians); QC = cos(pointBeta * BOKZMath::transToRadians);
    RS = sin(pointAlpha * BOKZMath::transToRadians); RC = cos(pointAlpha * BOKZMath::transToRadians);
    trMat[0][0] = -PC*RS-PS*RC*QS;
    trMat[0][1] = PC*RC-PS*RS*QS;
    trMat[0][2] = PS*QC;
    trMat[1][0] = PS*RS-PC*RC*QS;
    trMat[1][1] = -PS*RC-PC*RS*QS;
    trMat[1][2] = PC*QC;
    trMat[2][0] = QC*RC;
    trMat[2][1] = QC*RS;
    trMat[2][2] = QS;

    return trMat;

}

QVector <quint32> firstMinDistanceTable(RecognizedInfo** distMatrix, quint32 countOfMins, quint32 objectIndex, quint32 size)
{
    QVector <float> vec(countOfMins, std::numeric_limits <float>::max());
    QVector <quint32> minIndexes(countOfMins, std::numeric_limits <quint32>::max());

    for (quint32 i = 0; i < size; i++)
    {
        if (i == objectIndex) continue;

        float dist = 0;
        if (objectIndex > i) dist = distMatrix[objectIndex][i].distance;
        else dist = distMatrix[i][objectIndex].distance;
        //qDebug() << dist;

        for (int j = 0; j < countOfMins; j++)
        {
            if (dist < vec[j])
            {
                for (int k = countOfMins - 1; k >= j; k--)
                {
                    if (k == j)
                    {
                        vec[k] = dist;
                        minIndexes[k] = i;
                        break;
                    }
                    else
                    {
                        vec[k] = vec[k - 1];
                        minIndexes[k] = minIndexes[k - 1];
                    }
                }
                break;
            }
        }
    }
    return minIndexes;
}

}






