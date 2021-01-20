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
        for (int i = 0; i < loopingLength; ++i)
        {
            result.setValue(i, this->values[i] - other.values[i]);
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
        for (int i = 0; i < this->dimension; ++i)
        {
            result.setValue(i, this->values[i] * scalar);
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
};

int main(int argc, char const *argv[])
{
    Matrix<int> t(2, 4);
    int count = 0;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            t.setValue(i, j, count);
            ++count;
        }
    }
    t.setValue(1, 2, -99);

    // cout << t.getValue(1, 2) << endl;

    t = t.getTranspose();

    t = t + t;
    t.print(',');

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