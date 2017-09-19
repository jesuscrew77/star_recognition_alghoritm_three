#include <QCoreApplication>
#include <starrecalgh.h>
int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    StarRecognizeThree alghoritm ("");
    alghoritm.makeTest();
    return a.exec();
}
