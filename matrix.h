#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <iomanip>
#include <random>
using namespace std;

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

#endif // MATRIX_H