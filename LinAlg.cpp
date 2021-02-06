#ifndef LINALG_H
#define LINALG_H
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
using namespace std;

template <class T>
class Vector
{
public:
    T *values;
    int dimension;

    Vector(int dimension, T values[])
    {
        this->dimension = dimension;
        this->values = values;
    }

    Vector(int dimension)
    {
        this->dimension = dimension;
        values = new T[dimension];
        for (int i = 0; i < dimension; ++i)
        {
            values[i] = 0;
        }
    }

    Vector() {}

    void setDimension(int dimension)
    {
        this->dimension = dimension;
        values = new T[dimension];
        for (int i = 0; i < dimension; ++i)
        {
            values[i] = 0;
        }
    }

    ~Vector() {}

    void setValue(int index, T value)
    {
        if (index >= dimension)
        {
            throw "Index Error";
        }
        values[index] = value;
    }

    T getValue(int index)
    {
        return values[index];
    }

    void print(char seperator, int precision)
    {
        std::cout << std::fixed;
        std::cout << std::setprecision(precision);
        for (int i = 0; i < this->dimension; ++i)
        {
            if (i < this->dimension - 1)
            {
                std::cout << values[i] << seperator << ' ';
            }
            else
            {
                std::cout << values[i];
            }
        }
        std::cout << "\n";
    }

    void print(int precision)
    {
        std::cout << std::fixed;
        std::cout << std::setprecision(precision);
        for (int i = 0; i < this->dimension; ++i)
        {
            if (i < this->dimension - 1)
            {
                std::cout << values[i] << ' ';
            }
            else
            {
                std::cout << values[i];
            }
        }
        std::cout << "\n";
    }

    void print(char seperator)
    {
        for (int i = 0; i < this->dimension; ++i)
        {
            if (i < this->dimension - 1)
            {
                std::cout << values[i] << seperator << ' ';
            }
            else
            {
                std::cout << values[i];
            }
        }
        std::cout << "\n";
    }

    void print()
    {
        for (int i = 0; i < this->dimension; ++i)
        {
            if (i < this->dimension - 1)
            {
                std::cout << values[i] << ' ';
            }
            else
            {
                std::cout << values[i];
            }
        }
        std::cout << "\n";
    }

    template <class TT>
    Vector<T> operator+(const Vector<TT> &other)
    {
        int loopingLength;
        if (other.dimension >= this->dimension)
        {
            loopingLength = this->dimension;
        }
        else
        {
            loopingLength = other.dimension;
        }
        Vector<T> result(loopingLength);
        for (int i = 0; i < loopingLength; ++i)
        {
            result.setValue(i, this->values[i] + other.values[i]);
        }
        return result;
    }

    template <class TT>
    Vector<T> operator-(const Vector<TT> &other)
    {
        int loopingLength;
        if (other.dimension >= this->dimension)
        {
            loopingLength = this->dimension;
        }
        else
        {
            loopingLength = other.dimension;
        }
        Vector<T> result(loopingLength);
        T ac = 0;
        for (int i = 0; i < loopingLength; ++i)
        {
            ac = this->values[i] - other.values[i];
            if ((*typeid(this->getValue(0)).name() == 'f') && abs(ac) <= 0.00001)
            {
                ac = 0;
            }
            result.setValue(i, ac);
        }
        return result;
    }

    template <class TT>
    T operator*(const Vector<TT> &other) // DOT PRODUCT
    {
        if (other.dimension != this->dimension)
        {
            throw "Dimension error";
        }
        T result = 0;
        for (int i = 0; i < this->dimension; ++i)
        {
            result += this->values[i] * other.values[i];
        }
        return result;
    }

    template <class TT>
    Vector<T> operator*(TT scalar)
    {
        Vector<T> result(dimension);
        T ac = 0;
        for (int i = 0; i < this->dimension; ++i)
        {
            ac = this->values[i] * scalar;
            if (ac == -0)
            {
                ac = 0;
            }
            result.setValue(i, ac);
        }
        return result;
    }

    Vector<T> unitVector()
    {
        float mag = this->getMagnitude();
        Vector<T> result(dimension);
        for (int i = 0; i < this->dimension; ++i)
        {
            result.setValue(i, this->values[i] / mag);
        }
        return result;
    }

    template <class TT>
    Vector<T> crossProduct(const Vector<TT> &other);

    float getMagnitude()
    {
        float result = 0;
        for (int i = 0; i < this->dimension; ++i)
        {
            result += this->values[i] * this->values[i];
        }
        return sqrt(result);
    }

    template <class TT>
    float getAngleSeperating(Vector<TT> &other)
    {
        return acos(this->operator*(other) / (this->getMagnitude() * other.getMagnitude()));
    }
};

template <class T>
class Matrix
{
private:
    Vector<T> *rows;
    int row;
    int col;

public:
    /**
     * Construcor method for class.
     * 
     * @param row
     */
    Matrix(int row, int col, Vector<T> rows[])
    {
        this->row = row;
        this->col = col;
        this->rows = rows;
    }

    Matrix(int row, int col)
    {
        this->row = row;
        this->col = col;

        rows = new Vector<T>[row];
        for (int i = 0; i < row; ++i)
        {
            *(rows + i) = Vector<T>(col);
        }
    }

    ~Matrix() {}

    Vector<T> getRow(int row)
    {
        if (row >= this->row)
        {
            throw "Index error";
        }
        return *(rows + row);
    }

    Vector<T> getCol(int col)
    {
        if (col >= this->col)
        {
            throw "Index error";
        }
        Vector<T> v(row);
        for (int i = 0; i < row; ++i)
        {
            v.setValue(i, this->getValue(i, col));
        }
        return v;
    }

    void setValue(int row, int col, T value)
    {
        if (row >= this->row || col >= this->col)
        {
            throw "Dimension error";
        }
        Vector<T> v = *(rows + row);
        v.setValue(col, value);
    }

    T getValue(int row, int col)
    {
        Vector<T> v = *(rows + row);
        return v.getValue(col);
    }

    void print()
    {
        Vector<T> v;
        for (int i = 0; i < this->row; ++i)
        {
            v = *(rows + i);
            v.print();
        }
    }

    void print(char seperator)
    {
        Vector<T> v;
        for (int i = 0; i < this->row; ++i)
        {
            v = *(rows + i);
            v.print(seperator);
        }
    }

    void print(int precision)
    {
        Vector<T> v;
        for (int i = 0; i < this->row; ++i)
        {
            v = *(rows + i);
            v.print(precision);
        }
    }

    void print(char seperator, int precision)
    {
        Vector<T> v;
        for (int i = 0; i < this->row; ++i)
        {
            v = *(rows + i);
            v.print(seperator, precision);
        }
    }

    template <class TT>
    Matrix<T> operator+(Matrix<TT> &other)
    {
        if (this->row != other.row || this->col != other.col)
        {
            throw "Dimension error";
        }
        Matrix<T> result(this->row, this->col);
        for (int i = 0; i < row; ++i)
        {
            for (int j = 0; j < col; ++j)
            {
                result.setValue(i, j, this->getValue(i, j) + other.getValue(i, j));
            }
        }
        return result;
    }

    template <class TT>
    Matrix<T> operator*(Matrix<TT> &other)
    {
        if (this->row != other.col)
        {
            throw "Dimension error";
        }
        Matrix<T> result(this->row, other.col);
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = 0; j < other.col; ++j)
            {
                result.setValue(i, j, this->getRow(i) * other.getCol(j));
            }
        }
        return result;
    }

    template <class TT>
    Matrix<T> operator*(TT scalar)
    {
        Matrix<T> result(this->row, this->col);
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = 0; j < this->col; ++j)
            {
                result.setValue(i, j, this->getValue(i, j) * scalar);
            }
        }
        return result;
    }

    template <class TT>
    Vector<T> operator*(Vector<TT> &other)
    {
        if (this->row != other.dimension)
        {
            throw "Dimension error";
        }
        Vector<T> result(this->col);
        for (int i = 0; i < this->row; ++i)
        {
            result.setValue(i, this->getRow(i) * other);
        }
        return result;
    }

    template <class TT>
    Matrix<T> operator-(Matrix<TT> &other)
    {
        if (this->row != other.row || this->col != other.col)
        {
            throw "Dimension error";
        }
        Matrix<T> result(this->row, this->col);
        for (int i = 0; i < row; ++i)
        {
            for (int j = 0; j < col; ++j)
            {
                result.setValue(i, j, this->getValue(i, j) - other.getValue(i, j));
            }
        }
        return result;
    }

    template <class TT>
    Matrix<T> straightMultiply(Matrix<TT> &other)
    {
        if (this->row != other.row || this->col != other.col)
        {
            throw "Dimension error";
        }
        Matrix<T> result(this->row, this->col);
        for (int i = 0; i < row; ++i)
        {
            for (int j = 0; j < col; ++j)
            {
                result.setValue(i, j, this->getValue(i, j) * other.getValue(i, j));
            }
        }
        return result;
    }

    template <class TT>
    Matrix<T> straightDivide(Matrix<TT> &other)
    {
        if (this->row != other.row || this->col != other.col)
        {
            throw "Dimension error";
        }
        Matrix<T> result(this->row, this->col);
        for (int i = 0; i < row; ++i)
        {
            for (int j = 0; j < col; ++j)
            {
                result.setValue(i, j, (T)(this->getValue(i, j) / other.getValue(i, j)));
            }
        }
        return result;
    }

    Matrix<T> getTranspose()
    {
        Matrix<T> result(this->col, this->row);
        for (int i = 0; i < row; ++i)
        {
            for (int j = 0; j < col; ++j)
            {
                result.setValue(j, i, this->getValue(i, j));
            }
        }
        return result;
    }

    Vector<T> getDiagonal()
    {
        if (this->row != this->col)
        {
            throw "Non-square error";
        }
        Vector<T> v(this->row);
        for (int i = 0; i < row; ++i)
        {
            v.setValue(i, this->getValue(i, i));
        }
        return v;
    }

    int getTriangular()
    {
        if (this->row != this->col)
        {
            throw "Non-square error";
        }
        bool isUpper = true, isLower = true;
        Vector<T> diag = this->getDiagonal();
        for (int i = 0; i < this->row; ++i)
        {
            if (diag.getValue(i) == 0)
            {
                return Not;
            }
        }
        for (int i = 0; i < this->row - 1; ++i)
        {
            for (int j = 1 + i; j < this->col; ++j)
            {
                if (this->getValue(i, j) != 0)
                {
                    isUpper = false;
                    break;
                }
            }
        }
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = 0; j < i; ++j)
            {
                if (this->getValue(i, j) != 0)
                {
                    isLower = false;
                    break;
                }
            }
        }

        if (isLower && isUpper)
        {
            return Diagonal;
        }
        else if (isUpper)
        {
            return Upper;
        }
        else if (isLower)
        {
            return Lower;
        }
        else
        {
            return Not;
        }
    }

    enum Triangular
    {
        Upper,
        Lower,
        Diagonal,
        Not
    };

    Matrix<T> getPower(int power)
    {
        if (this->row != this->col)
        {
            throw "Non-square error";
        }
        if (power == 1)
        {
            return *this;
        }
        Matrix<T> m(this->row, this->col);
        if (this->getTriangular() == Diagonal)
        {
            for (int i = 0; i < 5; ++i)
            {
                m.setValue(i, i, pow(this->getValue(i, i), power));
            }
            return m;
        }
        else
        {
            m = this->getPower(power - 1);
            return this->operator*(m);
        }
    }

    T getLeading(int row)
    {
        for (int i = 0; i < this->col; ++i)
        {
            if (this->getValue(row, i) != 0)
            {
                return this->getValue(row, i);
            }
        }
        throw "Zero row";
    }

    int getLeadingIndex(int row)
    {
        for (int i = 0; i < this->col; ++i)
        {
            if (this->getValue(row, i) != 0)
            {
                return i;
            }
        }
        throw "Zero row";
    }

    bool isZeroRow(int row)
    {
        for (int i = 0; i < this->col; ++i)
        {
            if (this->getValue(row, i) != 0)
            {
                return false;
            }
        }
        return true;
    }

    template <class TT>
    void rowOperation(TT scalar, int selectRow)
    {
        if (selectRow >= this->row)
        {
            throw "Dimension error";
        }
        Vector<T> v = *(rows + selectRow);
        *(rows + selectRow) = v.operator*(scalar);
    }

    template <class TT>
    void rowOperation(TT scalar, int selectRow, int otherRow, char operation)
    {
        if (selectRow >= this->row)
        {
            throw "Dimension error";
        }
        if (operation == '+')
        {
            Vector<T> v = *(rows + selectRow);
            Vector<T> o = *(rows + otherRow);
            *(rows + selectRow) = o.operator*(scalar).operator+(v);
        }
        else if (operation == '-')
        {
            Vector<T> v = *(rows + selectRow);
            Vector<T> o = *(rows + otherRow);
            *(rows + selectRow) = o.operator*(scalar).operator-(v);
        }
        else
        {
            throw "Wrong operation error";
        }
    }

    void rowOperation(int selectRow, int otherRow, char operation)
    {
        if (selectRow >= this->row)
        {
            throw "Dimension error";
        }
        if (operation == '+')
        {
            Vector<T> v = *(rows + selectRow);
            Vector<T> o = *(rows + otherRow);
            *(rows + selectRow) = v.operator+(o);
        }
        else if (operation == '-')
        {
            Vector<T> v = *(rows + selectRow);
            Vector<T> o = *(rows + otherRow);
            *(rows + selectRow) = v.operator-(o);
        }
        else if (operation == 's')
        {
            Vector<T> v = *(rows + selectRow);
            Vector<T> o = *(rows + otherRow);
            *(rows + selectRow) = o;
            *(rows + otherRow) = v;
        }
        else
        {
            throw "Wrong operation error";
        }
    }

    int getRank()
    {
        int rank = 0;
        for (int i = 0; i < this->row; ++i)
        {
            try
            {
                this->getLeadingIndex(i);
                ++rank;
            }
            catch (const char *zr)
            {
                if (zr != "Zero row")
                {
                    throw zr;
                }
            }
        }
        return rank;
    }

    int getNullity()
    {
        return this->col - this->getRank();
    }

    Matrix<T> getMinorMatrix(int row, int col)
    {
        Matrix<T> minorMatrix(this->row - 1, this->col - 1);
        int rowCount = 0;
        int colCount;
        for (int i = 0; i < this->row; ++i)
        {
            colCount = 0;
            for (int j = 0; j < this->col; ++j)
            {
                if (i != row && j != col)
                {
                    minorMatrix.setValue(rowCount, colCount, this->getValue(i, j));
                    ++colCount;
                }
            }
            if (i != row)
            {
                ++rowCount;
            }
        }
        return minorMatrix;
    }

    T getMinor(int row, int col, bool useGaussian)
    {
        return this->getMinorMatrix(row, col).getDeterminant(useGaussian);
    }

    T getDeterminant(bool useGaussian)
    {
        if (this->row != this->col)
        {
            throw "Non-square error";
        }
        if (this->row == 2)
        {
            return this->getValue(0, 0) * this->getValue(1, 1) - this->getValue(0, 1) * this->getValue(1, 0);
        }
        else if (this->getTriangular() != Not)
        {
            Vector<T> diag = this->getDiagonal();
            T result = 1;
            for (int i = 0; i < this->col; ++i)
            {
                result *= diag.getValue(i);
            }
            return result;
        }
        else if (useGaussian)
        {
            if (*typeid(this->getValue(0, 0)).name() == 'i')
            {
                throw "Invalid (int) genertic type";
            }
            Matrix<T> copyMatrix(this->row, this->col);
            bool isZero = true;
            for (int i = 0; i < this->row; ++i)
            {
                for (int j = 0; j < this->col; ++j)
                {
                    copyMatrix.setValue(i, j, this->getValue(i, j));
                    if (isZero)
                    {
                        if (this->getValue(i, j) != 0)
                        {
                            isZero = false;
                        }
                    }
                }
            }
            if (isZero)
            {
                throw "Zero matrix error";
            }

            T nextValue;
            int nextValueIndex;
            int switchesNum;
            switchesNum = doEchelonRowOperations(copyMatrix, nextValue, nextValueIndex);
            Vector<T> diag = copyMatrix.getDiagonal();
            T result = 1;
            for (int i = 0; i < this->col; ++i)
            {
                result *= diag.getValue(i);
            }
            return result * pow(-1, switchesNum);
        }
        else
        {
            T minSum = 0;
            T sum = 0;
            int minRow = 0;
            int minCol = 0;
            bool checkRows = true;
            // Checking rows:
            for (int i = 0; i < this->row; ++i)
            {
                sum = 0;
                if (i == 0)
                {
                    for (int j = 0; j < this->col; ++j)
                    {
                        minSum += this->getValue(i, j);
                    }
                }
                else
                {
                    for (int j = 0; j < this->col; ++j)
                    {
                        sum += this->getValue(i, j);
                    }
                    if (sum < minSum)
                    {
                        minSum = sum;
                        minRow = i;
                    }
                }
            }
            // Checking cols:
            for (int j = 0; j < this->col; ++j)
            {
                sum = 0;
                for (int i = 0; i < this->row; ++i)
                {
                    sum += this->getValue(i, j);
                }
                if (sum < minSum)
                {
                    minSum = sum;
                    minCol = j;
                    checkRows = false;
                }
            }
            // Finding result:
            T result = 0;
            if (checkRows)
            {
                int signCheck = (minRow % 2);
                for (int j = 0; j < this->col; ++j)
                {
                    if (j % 2 == signCheck)
                    {
                        result += this->getValue(minRow, j) * this->getMinor(minRow, j, false);
                    }
                    else
                    {
                        result -= this->getValue(minRow, j) * this->getMinor(minRow, j, false);
                    }
                }
            }
            else
            {
                int signCheck = (minCol % 2);
                for (int i = 0; i < this->row; ++i)
                {
                    if (i % 2 == signCheck)
                    {
                        result += this->getValue(i, minCol) * this->getMinor(i, minCol, false);
                    }
                    else
                    {
                        result -= this->getValue(i, minCol) * this->getMinor(i, minCol, false);
                    }
                }
            }
            return result;
        }
    }

    Matrix<T> getAdjugate(bool useGaussian)
    {
        Matrix<T> adjugate(this->row, this->col);
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = 0; j < this->col; ++j)
            {
                if ((i + j) % 2 == 0)
                {
                    adjugate.setValue(i, j, this->getMinor(i, j, useGaussian));
                }
                else
                {
                    adjugate.setValue(i, j, -this->getMinor(i, j, useGaussian));
                }
            }
        }
        adjugate = adjugate.getTranspose();
        return adjugate;
    }

    Matrix<T> getInverse(bool useGaussian)
    {
        T det = this->getDeterminant(useGaussian);
        if (det != 0)
        {
            return this->getAdjugate(useGaussian).operator*((float)(1 / det));
        }
        else
        {
            throw "Non invertible error";
        }
    }

    T getMax()
    {
        T max = this->getValue(0, 0);
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = 0; j < this->col; ++j)
            {
                if (this->getValue(i, j) > max)
                {
                    max = this->getValue(i, j);
                }
            }
        }
        return max;
    }

    T getMin()
    {
        T min = this->getValue(0, 0);
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = 0; j < this->col; ++j)
            {
                if (this->getValue(i, j) < min)
                {
                    min = this->getValue(i, j);
                }
            }
        }
        return min;
    }

    int *getArgMax() // Will return a pointer.
    {                // Use "cout << *instance.getArgMax()<<", "<<*(instance.getArgMax() + 1);"
                     // to see the result.
        static int argMax[] = {0, 0};
        T max = this->getValue(0, 0);
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = 0; j < this->col; ++j)
            {
                if (this->getValue(i, j) > max)
                {
                    max = this->getValue(i, j);
                    argMax[0] = i;
                    argMax[1] = j;
                }
            }
        }
        return argMax;
    }

    int *getArgMin() // Will return a pointer.
    {                // Use "cout << *instance.getArgMin()<<", "<<*(instance.getArgMin() + 1);"
                     // to see the result.
        static int argMin[] = {0, 0};
        T min = this->getValue(0, 0);
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = 0; j < this->col; ++j)
            {
                if (this->getValue(i, j) < min)
                {
                    min = this->getValue(i, j);
                    argMin[0] = i;
                    argMin[1] = j;
                }
            }
        }
        return argMin;
    }

    int *find(T value) // Will throw an exception if value not in matrix.
    {                  // Will return a pointer.
                       // Use "cout << "find: " << *instance.find(value) << ", " << *(instance.find(value) + 1);"
                       // to see the result.
        static int arg[] = {0, 0};
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = 0; j < this->col; ++j)
            {
                if (this->getValue(i, j) == value)
                {
                    arg[0] = i;
                    arg[1] = j;
                    return arg;
                }
            }
        }
        throw "Cannot find value";
    }

    Matrix<T> getReducedEchelon() // Uses Gauss-Jordan Elimination
    {
        Matrix<T> copyMatrix = this->getEchelon();
        T nextValue;
        int nextValueIndex;
        doReducedEchelonRowOperations(copyMatrix, nextValue, nextValueIndex);
        return copyMatrix;
    }

    Matrix<T> getEchelon()
    {
        if (*typeid(this->getValue(0, 0)).name() == 'i')
        {
            throw "Invalid (int) genertic type";
        }
        Matrix<T> copyMatrix(this->row, this->col);
        bool isZero = true;
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = 0; j < this->col; ++j)
            {
                copyMatrix.setValue(i, j, this->getValue(i, j));
                if (isZero)
                {
                    if (this->getValue(i, j) != 0)
                    {
                        isZero = false;
                    }
                }
            }
        }
        if (isZero)
        {
            throw "Zero matrix error";
        }

        T nextValue;
        int nextValueIndex;
        doEchelonRowOperations(copyMatrix, nextValue, nextValueIndex);
        return copyMatrix;
    }

    int doEchelonRowOperations(Matrix<T> &copyMatrix, T &nextValue, int &nextValueIndex)
    {
        int switchesNum = 0;
        for (int i = 0; i < copyMatrix.row; ++i)
        {
            switchesNum += copyMatrix.orderMatrix();
            if (!copyMatrix.isZeroRow(i))
            {
                nextValue = copyMatrix.getLeading(i);
                nextValueIndex = copyMatrix.getLeadingIndex(i);
                for (int k = i + 1; k < copyMatrix.row; ++k)
                {
                    if (copyMatrix.getCol(nextValueIndex).getValue(k) != 0)
                    {
                        copyMatrix.rowOperation(copyMatrix.getValue(k, nextValueIndex) / copyMatrix.getValue(i, nextValueIndex), k, i, '-');
                    }
                }
            }
        }
        return switchesNum;
    }

    void doReducedEchelonRowOperations(Matrix<T> &copyMatrix, T &nextValue, int &nextValueIndex)
    {
        for (int i = copyMatrix.row - 1; i >= 0; --i)
        {
            copyMatrix.orderMatrix();
            if (!copyMatrix.isZeroRow(i))
            {
                nextValue = copyMatrix.getLeading(i);
                nextValueIndex = copyMatrix.getLeadingIndex(i);
                if (nextValue != 1)
                {
                    copyMatrix.rowOperation(1 / nextValue, i);
                }
                for (int k = i - 1; k >= 0; --k)
                {
                    if (copyMatrix.getCol(nextValueIndex).getValue(k) != 0)
                    {
                        copyMatrix.rowOperation(copyMatrix.getValue(k, nextValueIndex) / copyMatrix.getValue(i, nextValueIndex), k, i, '-');
                    }
                }
            }
        }
    }

    int orderMatrix() // Orders the matrix s.t. the matrix becomes as much echelon as possible.
    {
        int switchesNum = 0;
        bool brokeLoop;
        for (int i = 0; i < this->row; ++i)
        {
            brokeLoop = false;
            for (int j = 0; j < this->col; ++j)
            {
                if (this->getValue(i, j) == 0)
                {
                    for (int k = i; k < this->row; ++k)
                    {
                        if (this->getValue(k, j) != 0)
                        {
                            this->rowOperation(i, k, 's');
                            ++switchesNum;
                            brokeLoop = true;
                            break;
                        }
                    }
                }
                else
                {
                    break;
                }
                if (brokeLoop)
                {
                    break;
                }
            }
        }
        return switchesNum;
    }

    template <class TT>
    Matrix<T> appendRows(Matrix<TT> &other)
    {
        if (this->col != other.col)
        {
            throw "Mismatched matrices error";
        }
        Matrix<T> appendedMatrix(this->row + other.row, this->col);
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = 0; j < this->col; ++j)
            {
                appendedMatrix.setValue(i, j, this->getValue(i, j));
            }
        }
        for (int i = this->row; i < this->row + other.row; ++i)
        {
            for (int j = 0; j < other.col; ++j)
            {
                appendedMatrix.setValue(i, j, (T)other.getValue(i - this->row, j));
            }
        }
        return appendedMatrix;
    }

    template <class TT>
    Matrix<T> appendCols(Matrix<TT> &other)
    {
        if (this->row != other.row)
        {
            throw "Mismatched matrices error";
        }
        Matrix<T> appendedMatrix(this->row, this->col + other.col);
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = 0; j < this->col; ++j)
            {
                appendedMatrix.setValue(i, j, this->getValue(i, j));
            }
        }
        for (int i = 0; i < this->row; ++i)
        {
            for (int j = this->col; j < this->col + other.col; ++j)
            {
                appendedMatrix.setValue(i, j, (T)other.getValue(i, j - this->col));
            }
        }
        return appendedMatrix;
    }
};

int main(int argc, char const *argv[])
{
    int size = 100;
    int sizeY = size;
    Matrix<double> t(size, sizeY);
    // int count = 1;
    // for (int i = 0; i < size; ++i)
    // {
    //     for (int j = 0; j < sizeY; ++j)
    //     {
    //         t.setValue(i, j, rand() %10);
    //         ++count;
    //     }
    // }

    for (int i = 0; i < size; ++i)
    {
        t.setValue(i, i, i + 1);
    }
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

    // t = t.getPower(4);

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

    // Vector<float> test(3);
    // test.setValue(0, 55);
    // test.setValue(1, 0);
    // test.setValue(2, 0);

    // Vector<int> test2(3);
    // test2.setValue(0, 55);
    // test2.setValue(1, 0);
    // test2.setValue(2, 0);
    // test = test2;

    // cout << test.getMagnitude();

    // test = test * 5;

    // test = test.unitVector();

    // cout << test.getAngleSeperating(test2);

    // test.print();

    return 0;
}

#endif // LINALG_H