/**
    Matrix: Matrix class functions header.
    @author Nabil NY Mansour
*/

#ifndef MATRIX_H
#define MATRIX_H
#include "Vector.h"

template <class T>
class Matrix
{
public:
    Vector<T> *rows;
    int row;
    int col;

    /**
     * Construcor method for class Matrix.
     * @param row: number of rows of the matrix.
     * @param col: number of columns of the matrix.
     */
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

    /**
     * Destrucor method for class Matrix.
     */
    ~Matrix() {}

    /**
     * Method for returning a row vector.
     * @param row: the row that is to be returned.
     * @throws a string "Index error" if the row input is illegal.
     * @returns an array of class "Vector" of that specific row.
     */
    Vector<T> getRow(int row)
    {
        if (row >= this->row)
        {
            throw "Index error";
        }
        return *(rows + row);
    }

    /**
     * Method for returning a column vector.
     * @param col: the column that is to be returned.
     * @throws a string "Index error" if the col input is illegal.
     * @returns an array of class "Vector" of that specific column.
     */
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

    /**
     * Method for setting or changing the value of an element in the matrix.
     * @param row: the row of the element that is to be changed.
     * @param col: the column of the element that is to be changed.
     * @throws a string "Dimension error" if the row or col inputs are illegal.
     */
    void setValue(int row, int col, T value)
    {
        if (row >= this->row || col >= this->col)
        {
            throw "Dimension error";
        }
        Vector<T> v = *(rows + row);
        v.setValue(col, value);
    }

    /**
     * Method for returning the value of an element.
     * @param row: the row of the element that is to be returned.
     * @param col: the column of the element that is to be returned.
     * @returns the value of the element corrosponding to the given row and column inputted.
     */
    T getValue(int row, int col)
    {
        Vector<T> v = *(rows + row);
        return v.getValue(col);
    }

    /**
     * Method for printing the matrix.
     * Will utilize the method for printing in class Vector.
     * @param seperator: the seperator char that is to be placed inbetween each entry.
     * @param precision: the precision level of float or double type values that is to be printed.
     */
    void print(int precision, char seperator)
    {
        Vector<T> v;
        for (int i = 0; i < this->row; ++i)
        {
            v = *(rows + i);
            v.print(seperator, precision);
        }
    }

    /**
     * Method for printing the matrix.
     * Will utilize the method for printing in class Vector.
     * @param seperator: the seperator char that is to be placed inbetween each entry.
     */
    void print(char seperator)
    {
        Vector<T> v;
        for (int i = 0; i < this->row; ++i)
        {
            v = *(rows + i);
            v.print(seperator);
        }
    }

    /**
     * Method for printing the matrix.
     * Will utilize the method for printing in class Vector.
     * @param precision: the precision level of float or double type values that is to be printed.
     */
    void print(int precision)
    {
        Vector<T> v;
        for (int i = 0; i < this->row; ++i)
        {
            v = *(rows + i);
            v.print(precision);
        }
    }

    /**
     * Method for printing the matrix.
     * Will utilize the method for printing in class Vector.
     */
    void print()
    {
        Vector<T> v;
        for (int i = 0; i < this->row; ++i)
        {
            v = *(rows + i);
            v.print();
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

    /**
     * Method for multiplying two matrixes where each element is multiplied with the element it
     * corrosponds to and not in the regular fashion way.
     * @param other: the other matrix.
     * @throws a string "Dimension error" if the dimensions of the other matrix is not equal to the
     * this matrix's dimentions.
     * @returns an object of class Matrix of the resultant matrix.
     */
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

    /**
     * Method for dividing two matrixes where each element is divided with the element it
     * corrosponds to.
     * @param other: the other matrix.
     * @throws a string "Dimension error" if the dimensions of the other matrix is not equal to the
     * this matrix's dimentions.
     * @returns an object of class Matrix of the resultant matrix.
     */
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

    /**
     * Method for getting the transpose of the matrix.
     * @returns the transpose matrix of this matrix.
     */
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

    /**
     * Method for getting the diagonal vector of the matrix.
     * @throws a string "Non-square error" if the matrix is not a square matrix.
     * @returns an object of class Vector that holds the values of the diagonal of this matrix.
     */
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

    /**
     * Method for getting type of triangular for this matrix.
     * @throws a string "Non-square error" if the matrix is not a square matrix.
     * @returns an int corrosponding to an enumerator where:
     * 0 -> matrix is Upper triangular
     * 1 -> matrix is Lower triangular
     * 2 -> matrix is both Lower and Upper triangular
     * 3 -> matrix is neither Lower and Upper triangular
     */
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

    /**
     * Method for getting a matrix to the power of a integer.
     * @param power: the desired power for the matrix.
     * @throws a string "Non-square error" if the matrix is not a square matrix.
     * @returns a matrix that corrosponds to this matrix to the power of the input integer.
     */
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

    /**
     * Method for getting the leading elemet given a row.
     * @param row: the row to which the leading value of is desired.
     * @throws a string "Zero row" if the row is a zero row.
     * @returns an element that is the leading element of the inputted row.
     */
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

    /**
     * Method for getting the leading elemet index given a row.
     * @param row: the row to which the leading value index of is desired.
     * @throws a string "Zero row" if the row is a zero row.
     * @returns the index of the element that is the leading element of the inputted row.
     */
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

    /**
     * Method for checking if a row is a zero row.
     * @param row: the row that is to be checked if it is a zero row.
     * @returns a boolean where:
     * true -> the row is a zero row.
     * false -> the row is a non-zero row.
     */
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

    /**
     * Method for performing row operations on a matrix given a row.
     * @param scalar: the amount to scale the row with.
     * @param selectRow: the row to have the operation applied on.
     * @throws a string "Dimension error" if the row input is illegal.
     */
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

    /**
     * Method for performing row operations on a matrix given a row.
     * @param scalar: the amount to scale the row with.
     * @param selectRow: the row to have the operation applied on.
     * @param otherRow: the other row that will be used in the row operation.
     * @param operation: the operation that is to be performed between the selectRow
     * and the otherRow.
     * @note: legal operation are "+" and "-".
     * @throws a string "Dimension error" if the matrix is not a square matrix.
     * @throws a string "Wrong operation error" if the operation char is illegal.
     * the selected row and the other row.
     */
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

    /**
     * Method for performing row operations on a matrix given a row.
     * @param selectRow: the row to have the operation applied on.
     * @param otherRow: the other row that will be used in the row operation.
     * @param operation: the operation that is to be performed between the selectRow
     * and the otherRow.
     * @throws a string "Dimension error" if the matrix is not a square matrix.
     * @throws a string "Wrong operation error" if the operation char is illegal.
     * the selected row and the other row.
     */
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

    /**
     * Method for getting the rank of the matrix. 
     * @throws a string "Zero row" if the whole matrix is a zero matrix.
     * @returns the rank of the matrix as an integer value.
     */
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

    /**
     * Method for getting the nullity of the matrix.
     * @returns the nullity of the matrix as an integer value.
     */
    int getNullity()
    {
        return this->col - this->getRank();
    }

    /**
     * Method for getting the minor matrix of the matrix given a row and a column.
     * @param row: the row of the desired minor matrix.
     * @param col: the column of hte desired minor matrix.
     * @returns the minor matrix of the specified row and column in the form of a Matrix object.
     */
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

    /**
     * Method for getting the determinant of the minor matrix of the matrix given a row and a column.
     * @param row: the row of the desired minor matrix.
     * @param col: the column of hte desired minor matrix.
     * @param useGaussian: the boolean value for whether the algorithm that is to be used for finding
     * the determinant is the Gaussian algorithm or otherwise.
     * @note: using Gaussian is better for efficiency.
     * @returns the determinant of the minor matrix of the specified row and column 
     * in the form of the template value used.
     */
    T getMinor(int row, int col, bool useGaussian)
    {
        return this->getMinorMatrix(row, col).getDeterminant(useGaussian);
    }

    /**
     * Method for getting the determinant the matrix.
     * @param useGaussian: the boolean value for whether the algorithm that is to be used for finding
     * the determinant is the Gaussian algorithm or otherwise.
     * @note: using Gaussian is better for efficiency.
     * @throws a string "Non-square error" if the matrix is not a square matrix.
     * @throws a string "Invalid (int) genertic type" if the matrix template is an int.
     * @throws a string "Zero matrix error" if the whole matrix is a zero matrix.
     * @returns the determinant of the matrix in the form of the template value used.
     */
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

    /**
     * Method for getting the adjugate of the matrix.
     * @param useGaussian: the boolean value for whether the algorithm that is to be used for finding
     * the adjugate is the Gaussian algorithm or otherwise.
     * @note: using Gaussian is better for efficiency.
     * @returns the adjugate of this matrix in the form of a Matrix object.
     */
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

    /**
     * Method for getting the inverse matrix the matrix.
     * @param useGaussian: the boolean value for whether the algorithm that is to be used for finding
     * the determinant is the Gaussian algorithm or otherwise.
     * @note: using Gaussian is better for efficiency.
     * @throws a string "Non invertible error" if the matrix is not invertible.
     * @returns the determinant of the matrix in the form of the template value used.
     */
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

    /**
     * Method for getting the maximum value of the matrix.
     * @returns the maximum value of the matrix in the form of the template value used.
     */
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

    /**
     * Method for getting the minimum value of the matrix.
     * @returns the minimum value of the matrix in the form of the template value used.
     */
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

    /**
     * Method for getting the location of the maximum value of the matrix.
     * @returns the x and y values of the maximum value of the matrix in the form
     * of a pointer to an array of two values.
     */
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

    /**
     * Method for getting the location of the minimum value of the matrix.
     * @returns the x and y values of the minimum value of the matrix in the form
     * of a pointer to an array of two values.
     */
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

    /**
     * Method for getting the location of a given value of the matrix.
     * @param value: the value to be found.
     * @throws a string "Cannot find value" if the value does not exist in the matrix.
     * @returns the x and y values of the minimum value of the matrix in the form
     * of a pointer to an array of two values.
     */
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

    /**
     * Method for getting the reduced row echelon form of this matrix.
     * @note: will use Gauss-Jordan Elimination algorithm.
     * @returns a reduced row echelon version of this matrix in the form of a Matrix object.
     */
    Matrix<T> getReducedEchelon() // Uses Gauss-Jordan Elimination
    {
        Matrix<T> copyMatrix = this->getEchelon();
        T nextValue;
        int nextValueIndex;
        doReducedEchelonRowOperations(copyMatrix, nextValue, nextValueIndex);
        return copyMatrix;
    }

    /**
     * Method for getting the row echelon form of this matrix.
     * @throws a string "Invalid (int) genertic type" if the matrix template is an int.
     * @throws a string "Zero matrix error" if the whole matrix is a zero matrix.
     * @returns a row echelon version of this matrix in the form of a Matrix object.
     */
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
};

#endif // MATRIX_H