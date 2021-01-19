#ifndef MATRIX_H
#define MATRIX_H
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

    Vector(T values[])
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

    ~Vector() {}

    void setValue(int index, T value)
    {
        if (index >= dimension)
        {
            throw "Index Error";
        }
        values[index] = value;
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

int main(int argc, char const *argv[])
{
    Vector<float> test(3);
    test.setValue(0, 55);
    test.setValue(1, 0);
    test.setValue(2, 0);

    Vector<int> test2(3);
    test2.setValue(0, -55);
    test2.setValue(1, 0);
    test2.setValue(2, 0);
    // test = test2;

    // cout << test.getMagnitude();

    // test = test * 5;

    // test = test.unitVector();

    cout<<test.getAngleSeperating(test2);
    // test.print();
    return 0;
}

/*
class matrix
{
private:
    int size;
    float *values;

public:
    matrix(int size);

    ~matrix();

    void print();

    void makeIdentity();

    void makeRandom(int from, int to);

    void setValue(int i, int j, float value);

    float getValue(int i, int j);

    matrix operator+(matrix other);
};

matrix::matrix(int size)
{
    this->size = size;
    this->values = new float[size * size];
    for (int i = 0; i < size * size; ++i)
    {
        values[i] = 0;
    }
}

matrix::~matrix()
{
}

void matrix::print()
{
    for (int i = 0; i < this->size; ++i)
    {
        cout << '{';
        for (int j = 0; j < this->size; ++j)
        {
            if (j < this->size - 1)
            {
                cout << values[(i * this->size) + j] << ", ";
            }
            else
            {
                cout << values[(i * this->size) + j] << '}' << endl;
            }
        }
    }
}

void matrix::makeIdentity()
{
    int count = 0;
    while (count < this->size * this->size)
    {
        for (int i = 0; i < this->size * this->size; ++i)
        {
            if (i == count)
            {
                this->values[i] = 1;
                count += this->size + 1;
            }
        }
    }
}

void matrix::makeRandom(int from, int to)
{
    for (int i = 0; i < this->size * this->size; ++i)
    {
        this->values[i] = rand() % to + from;
    }
}

void matrix::setValue(int i, int j, float value)
{
    this->values[(i * this->size) + j] = value;
}

float matrix::getValue(int i, int j)
{
    return this->values[(i * this->size) + j];
}

*/
#endif // MATRIX_H