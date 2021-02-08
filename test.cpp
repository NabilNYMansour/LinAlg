#include "LinAlg.h"
int main(int argc, char const *argv[])
{
    int size = 10;
    int sizeY = size;
    Matrix<double> t(size, sizeY);
    int count = 1;
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < sizeY; ++j)
        {
            t.setValue(i, j, rand() % 10);
            ++count;
        }
    }

    // for (int i = 0; i < size; ++i)
    // {
    //     t.setValue(i, i, i + 1);
    // }
    // t.setValue(0, 0, 1);
    // t.setValue(0, 1, 6);
    // t.setValue(0, 2, 4);

    // t.setValue(1, 0, 2);
    // t.setValue(1, 1, 4);
    // t.setValue(1, 2, -1);

    // t.setValue(2, 0, -1);
    // t.setValue(2, 1, 2);
    // t.setValue(2, 2, 5);

    // t.rowOperation(0, 4, 's');

    t.print(',');
    cout << endl;

    // t.setValue(1, 2, -99);

    // cout << t.getValue(1, 2) << endl;

    // t = t.getTranspose();

    // t = t + t;

    // t.getRow(1).print(',');

    // t.getCol(3).print();

    // t = t.straightDivide(t);

    // t = t * t;

    // Vector<int> v(5);
    // v.setValue(0, 1);
    // v.setValue(1, 1);

    // v = t * v;
    // v = t.getDiagonal();
    t = t.getPower(4);
    t.print(',');
    cout << endl;
    // t.rowOperation(2, 0, 2, '+');
    // t.rowOperation(2, 2);

    // cout << t.getMinor(1, 1)<<endl;
    // t = t.getMinorMatrix(1, 1);
    // t.setValue(1, 1, 0);
    // t.setValue(2, 1, 0);
    // t.setValue(3, 1, 0);
    // t.setValue(2, 2, 22);
    // t.setValue(3, 3, 22);
    // t.orderMatrix();

    // t.getToEchelonOperations();
    // Matrix<float> t2 = t;
    // t2 = t.getEchelon();
    // t = t.getReducedEchelon();
    // t2 = t.getInverse(false);

    // t2.print(',');
    // cout << endl;

    // cout << t.getDeterminant(false) << endl;
    t = t.getInverse(true);

    // t = t.appendRows(t);
    // t = t.appendCols(t);

    t.print(',');

    // cout << t.getDeterminant(false) << endl;
    // cout << "max: " << *t.getArgMax() << ", " << *(t.getArgMax() + 1) << endl;
    // cout << "min: " << *t.getArgMin() << ", " << *(t.getArgMin() + 1) << endl;
    // cout << "find: " << *t.find(4) << ", " << *(t.find(4) + 1) << endl;
    // cout << "find: " << *t.find(17);

    // cout << t.getNullity() << endl;

    // cout << t.getEchelonType();
    // cout<<t.getLeading(4);
    // cout << endl;
    // v.print(',');

    // cout << t.getTriangular();
    //---------------------------------------------------------------------//

    Vector<float> test(3);
    test.setValue(0, 1);
    test.setValue(1, 2);
    test.setValue(2, 3);

    Vector<float> test2(3);
    test2.setValue(0, 4);
    test2.setValue(1, 5);
    test2.setValue(2, 6);
    // test = test2;

    Vector<float> test3;
    test3 = test.crossProduct(test2);

    test3.print();

    // cout << test.getMagnitude();

    // test = test * 5;

    // test = test.unitVector();

    // cout << test.getAngleSeperating(test2);

    // test.print();

    return 0;
}
