/**
    LinAlg: Linear Algebra classes library
    @author Nabil NY Mansour
*/

#ifndef LINALG_H
#define LINALG_H
#include "Matrix.h"

/**
 * Method for connecting two matrices horizontally.
 * @param m1: the first matrix.
 * @param m2: the second matrix.
 * @throws a string "Mismatched matrices error" if the rows of the two matrices is not equal.
 * @returns the combined matrix in the form of a Matrix object.
*/
template <class T, class TT>
Matrix<T> connectMatricesHorizontally(Matrix<T> &m1, Matrix<TT> &m2)
{
    if (m1.row != m2.row)
    {
        throw "Mismatched matrices error";
    }
    Matrix<T> appendedMatrix(m1.row, m1.col + m2.col);
    for (int i = 0; i < m1.row; ++i)
    {
        for (int j = 0; j < m1.col; ++j)
        {
            appendedMatrix.setValue(i, j, m1.getValue(i, j));
        }
    }
    for (int i = 0; i < m1.row; ++i)
    {
        for (int j = m1.col; j < m1.col + m2.col; ++j)
        {
            appendedMatrix.setValue(i, j, (T)m2.getValue(i, j - m1.col));
        }
    }
    return appendedMatrix;
}

/**
 * Method for connecting two matrices vertically.
 * @param m1: the first matrix.
 * @param m2: the second matrix.
 * @throws a string "Mismatched matrices error" if the columns of the two matrices is not equal.
 * @returns the combined matrix in the form of a Matrix object.
*/
template <class T, class TT>
Matrix<T> connectMatricesVerticallly(Matrix<T> &m1, Matrix<TT> &m2)
{
    if (m1.col != m2.col)
    {
        throw "Mismatched matrices error";
    }
    Matrix<T> appendedMatrix(m1.row + m2.row, m1.col);
    for (int i = 0; i < m1.row; ++i)
    {
        for (int j = 0; j < m1.col; ++j)
        {
            appendedMatrix.setValue(i, j, m1.getValue(i, j));
        }
    }
    for (int i = m1.row; i < m1.row + m2.row; ++i)
    {
        for (int j = 0; j < m2.col; ++j)
        {
            appendedMatrix.setValue(i, j, (T)m2.getValue(i - m1.row, j));
        }
    }
    return appendedMatrix;
}

#endif // LINALG_H