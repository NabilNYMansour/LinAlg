/**
    Vector: Vector class functions header.
    @author Nabil NY Mansour
*/

#ifndef VECTOR_H
#define VECTOR_H
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

template <class T>
class Vector
{
public:
    T *values;
    int dimension;

    /**
     * Construcor method for class Vector.
     * @param dimenstion: number of elements of the vector.
     */
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

    ~Vector() {}

    /**
     * Method for setting the value of an element in the vector.
     * @param index: the index of the element.
     * @param value: the value of the element.
     * @throws a string "Index Error" if index value is illegal.
     */
    void setValue(int index, T value)
    {
        if (index >= dimension)
        {
            throw "Index Error";
        }
        values[index] = value;
    }

    /**
     * Method for getting the value of an element in the vector.
     * @param index: the index of the element.
     * @returns the value of the desired element
     */
    T getValue(int index)
    {
        return values[index];
    }

    /**
     * Method for printing the vector.
     * @param seperator: the seperator char that is to be placed inbetween each entry.
     * @param precision: the precision level of float or double type values that is to be printed.
     */
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

    /**
     * Method for printing the vector.
     * @param precision: the precision level of float or double type values that is to be printed.
     */
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

    /**
     * Method for printing the vector.
     * @param seperator: the seperator char that is to be placed inbetween each entry.
     */
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

    /**
     * Method for printing the vector.
     */
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

    /**
     * Method for getting the unit vector that corrosponds to this vector.
     * @returns the unit vector in the form of a Vector object.
     */
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

    /**
     * Method for getting the cross product that corrosponds to this vector.
     * @throws a string "Dimension error" if the dimension of this vector is not 3.
     * @returns the resultnant cross produect vector in the form of a Vector object.
     */
    template <class TT>
    Vector<T> crossProduct(Vector<TT> &other)
    {
        if (dimension != 3)
        {
            throw "Dimension error";
        }
        Vector<T> returnVector(3);
        returnVector.setValue(0, this->getValue(1) * other.getValue(2) - this->getValue(2) * other.getValue(1));
        returnVector.setValue(1, this->getValue(2) * other.getValue(0) - this->getValue(0) * other.getValue(2));
        returnVector.setValue(2, this->getValue(0) * other.getValue(1) - this->getValue(1) * other.getValue(0));
        return returnVector;
    }

    /**
     * Method for getting the magnitude that corrosponds to this vector.
     * @returns the magnitude of this vector in the form of a float object.
     */
    float getMagnitude()
    {
        float result = 0;
        for (int i = 0; i < this->dimension; ++i)
        {
            result += this->values[i] * this->values[i];
        }
        return sqrt(result);
    }

    /**
     * Method for getting the angle seperating this vector with another vector.
     * @param other: the other vector.
     * @returns the angle in the form of a float object.
     */
    template <class TT>
    float getAngleSeperating(Vector<TT> &other)
    {
        return acos(this->operator*(other) / (this->getMagnitude() * other.getMagnitude()));
    }
};

#endif // VECTOR_H