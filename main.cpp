#include "Matrix.h"
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>

//A*x=y
Matrix<double> gauss(Matrix<double> A, Matrix<double> b) {
    int n{A.getRows()};
    Matrix<double> x(n, 1);
    int i, j, k, m, rowx;
    double xfac, temp, temp1, amax;

    rowx = 0;
    for (k = 0; k < n - 1; ++k) {
        amax = std::fabs(A(k, k));
        m = k;
        for (i = k + 1; i < n; i++) {
            xfac = std::fabs(A(i, k));
            if (xfac > amax) {
                amax = xfac;
                m = i;
            }
        }
        if (m != k) {
            rowx = rowx + 1;
            temp1 = b(k, 0);
            b(k, 0) = b(m, 0);
            b(m, 0) = temp1;
            for (j = k; j < n; j++) {
                temp = A(k, j);
                A(k, j) = A(m, j);
                A(m, j) = temp;
            }
        }
        for (i = k + 1; i < n; ++i) {
            xfac = A(i, k) / A(k, k);

            for (j = k + 1; j < n; ++j) {
                A(i, j) = A(i, j) - xfac * A(k, j);
            }
            b(i, 0) = b(i, 0) - xfac * b(k, 0);
        }

    }

    for (j = 1; j <= n; ++j) {
        k = n - j;
        x(k, 0) = b(k, 0);
        for (i = k + 1; i < n; ++i) {
            x(k, 0) = x(k, 0) - A(k, i) * x(i, 0);
        }
        x(k, 0) = x(k, 0) / A(k, k);
    }
    return x;
}

double norm(Matrix<double> v) {
    double cumulativeSum{0.0};
    for (int i = 0; i < v.getRows(); ++i)
        cumulativeSum += (v(i, 0) * v(i, 0));
    return sqrt(cumulativeSum);
}

double norm(const std::vector<double>& v) {
    double cumulativeSum{0.0};
    for (double i : v)
        cumulativeSum += (i * i);
    return sqrt(cumulativeSum);
}

void Normalize(Matrix<double> &v) {
    double vecNorm = norm(v);

    for (int i = 0; i < v.getRows(); ++i) {
        double temp = v(i, 0) * (1.0 / vecNorm);
        v(i, 0) = temp;
    }
}

Matrix<double> get_matrix(int s) {
    Matrix<double> A(s, s);
    for (int i = 0; i < s; ++i) {
        for (int k = 0; k < s; ++k) {
            A(i, k) = pow((static_cast<double>(std::abs(i - k) + 1)), -3.0)
                      + pow(static_cast<double>(i + k - s + 2), 2.0);
        }
    }
    return A;
}

double powermethod(const Matrix<double> &A) {
    std::vector<double> x_0(A.getCols(), 1.0);
    Matrix<double> x_init, x_temp, x_iter;
    x_init.pushBackColumn(x_0);
    x_temp = x_init;
    double eps{1e-16};
    double lambda{0.0};
    double old_lambda{1.0};

    while (std::fabs(1.0 - lambda / old_lambda) > eps) {

        old_lambda = lambda;
        x_iter = A * x_temp;
        lambda = x_iter(0, 0);
        x_iter /= lambda;
        x_temp = x_iter;
    }
    return lambda;
}


double inverse_powermethod(const Matrix<double> &A) {
    double eps{1e-16};
    std::vector<double> x_0(A.getCols(), 1.0);
    Matrix<double> x_temp, x_iter(A.getRows(), 1);
    x_temp.pushBackColumn(x_0);
    double lambda{0.0};
    double old_lambda{1.0};

    while (std::fabs(1.0 - lambda / old_lambda) > eps) {
        old_lambda = lambda;
        x_iter = gauss(A, x_temp);
        lambda = x_iter(0, 0);
        x_iter /= lambda;
        x_temp = x_iter;
    }

    return 1.0 / lambda;

}

double getValueFromVector(const Matrix<double> &A, Matrix<double> eigenVector) {
    double temp{eigenVector(0, 0)};
    eigenVector = A * eigenVector;
    return eigenVector(0, 0) / temp;
}

std::vector<double> generate_orthogonal(std::vector<std::vector<double>> vectors) {
    std::vector<double> approx(vectors[0].size(), 1.0);
    const std::vector<double> &orthogonal{approx};
    double s;
    vectors.push_back(approx);
    Matrix<double> q(vectors[0].size(), vectors.size()), r(vectors[0].size(), vectors.size());

    for (int i = 0; i < vectors.size(); ++i) {
        s = norm(vectors[i]);
        r(i, i) = s;
        std::vector<double> temp_q;
        for (int j = 0; j < vectors[0].size(); ++j) {
            q(j, i) = vectors[i][j] / s;
            temp_q.push_back(q(j, i));
        }
        for (int j = i + 1; j < vectors.size(); ++j) {
            s = std::inner_product(vectors[j].begin(), vectors[j].end(), temp_q.begin(), 0.0);
            r(j, i) = s;
            for (int k = 0; k < vectors[0].size(); ++k) {
                vectors[j][k] -= temp_q[k] * s;
            }
        }
    }
    return vectors.back();
}


std::vector<double> inverse_powermethod_modified(const Matrix<double> &A, double eigenValue) {
    const double eps{1e-9};
    std::vector<double> eigenVector(A.getCols(), 1.0);
    int count{0};
    double delta{1e6};
    Matrix<double> x_iter, x_temp, tempMatrix(A.getRows(), A.getCols());
    Matrix<double> I = Matrix<double>::createIdentity(A.getCols());
    x_iter.pushBackColumn(eigenVector);
    while (delta > eps && count < 300) {
        x_temp = x_iter;
        I *= eigenValue;
        tempMatrix = A - I;
        x_iter = gauss(tempMatrix, x_iter);
        Normalize(x_iter);
        x_temp *= x_iter - x_temp;
        delta = norm(x_iter - x_temp);

        count++;
    }

    for (int i = 0; i < A.getRows(); ++i) {
        eigenVector[i] = x_iter(i, 0);
    }

    return eigenVector;
}

Matrix<double> inverse_powermethod_modified(Matrix<double> A, double t_n, const std::vector<double> &eigenVectorInit) {
    double eps{1e-8};
    Matrix<double> x, x_n(A.getCols(), 1);
    x.pushBackColumn(eigenVectorInit);
    double old_t_n = t_n + 1.0;
    double lambda;
    int count = 0;

    Matrix<double> I = Matrix<double>::createIdentity(A.getCols());
    while (fabs(t_n / old_t_n - 1.0) > eps && count < 100) {
        for (int i = 0; i < A.getCols(); ++i) {
            A(i, i) -= t_n;
        }
        x_n = gauss(A, x);
        old_t_n = t_n;
        lambda = x_n(0, 0);
        t_n += 1 / lambda;
        x_n /= lambda;
        x = x_n;
        count += 1;
        Normalize(x_n);
    }

    return x;
}

std::vector<double> find_them_all(const Matrix<double> &A, const double &first) {
    std::vector<double> x_0(A.getCols(), 1.0);
    std::vector<double> eigenValues, eigenVector, orthogonal(A.getCols());
    Matrix<double> eigenVectorMatrix(A.getRows(), A.getCols()), zeros(A.getRows(), 1), initApproximation(A.getRows(),1);

    std::vector<std::vector<double>> eigenVectors;
    eigenValues.reserve(A.getCols());
    eigenValues.push_back(first);
    eigenVector = inverse_powermethod_modified(A, first);
    eigenVectors.push_back(eigenVector);
    for (int i = 1; i < A.getRows() - 1; ++i) {
        eigenVectors.push_back(generate_orthogonal(eigenVectors));
        Matrix<double> temp;
        temp.pushBackColumn(eigenVectors[i]);
        initApproximation = inverse_powermethod_modified(A, getValueFromVector(A, temp), eigenVectors[i]);
        for (int j = 0; j < A.getRows(); ++j) {
            eigenVectors[i][j] = initApproximation(j, 0);
        }
        eigenValues.push_back(getValueFromVector(A, initApproximation));

    }
    eigenValues.push_back(1e5);
    std::sort(eigenValues.begin(), eigenValues.end());
    return eigenValues;
}

int main() {
    Matrix<double> A = get_matrix(20);
    std::cout << "Matrix:\n";
    std::cout << std::fixed << A << '\n';
    std::cout << "||||||||||||||||||||||||||||||||\n\n";
    std::cout << "Eigenvalues:\n";
    double lambda_max = powermethod(A);
    double lambda_min = inverse_powermethod(A);
    std::cout << "lambda_min = " << std::setprecision(16) << std::fixed << lambda_min << '\t'
              << "lambda_max = " << std::setprecision(16) << std::fixed << lambda_max << '\n';
    std::cout << "||||||||||||||||||||||||||||||||\n\n";
    std::vector<double> values = find_them_all(A, lambda_min);
    values.back() = lambda_max;
    std::cout << "Values: \n";
    for (auto i : values) {
        std::cout << i << '\n';
    }
}