#include "Matrix.h"
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <random>
#include "Eigen/Eigenvalues"

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

//Matrix<double> tridiagonalMatrixMethod(Matrix<double> A, Matrix<double> b)
Matrix<double> solveTridiagonal(const Matrix<double> &A, const Matrix<double> &input_vector) {
    assert(A.getCols() == A.getRows());
    assert(A.getCols() == input_vector.getRows());
    assert(input_vector.getRows() >= 2);

    Matrix<double> modifiedMainDiagonal(input_vector.getRows(), 1);
    Matrix<double> result(input_vector.getRows(), 1);

    modifiedMainDiagonal(0, 0) = A(0, 0);
    result(0, 0) = input_vector(0, 0);

    for (int i = 1; i < input_vector.getRows(); ++i) {
        const double ratio = A(i, i - 1) / modifiedMainDiagonal(i - 1, 0);
        modifiedMainDiagonal(i, 0) = A(i, i) - ratio * A(i - 1, i);
        result(i, 0) = input_vector(i, 0) - ratio * result(i - 1, 0);
    }

    int n = input_vector.getRows() - 1;
    result(n, 0) /= modifiedMainDiagonal(n, 0);

    for (int i = n - 1; i >= 0; --i) {   // while not overflow
        result(i, 0) = (result(i, 0) - result(i + 1, 0) * A(i, i + 1)) / modifiedMainDiagonal(i, 0);
    }

    return result;
}

double norm(Matrix<double> v) {
    double cumulativeSum{0.0};
    for (int i = 0; i < v.getRows(); ++i)
        cumulativeSum += (v(i, 0) * v(i, 0));
    return sqrt(cumulativeSum);
}

double norm(std::vector<double> v) {
    double cumulativeSum{0.0};
    for (int i = 0; i < v.size(); ++i)
        cumulativeSum += (v[i] * v[i]);
    return sqrt(cumulativeSum);
}

void Normalize(Matrix<double> &v) {
    // Compute the vector norm.
    double vecNorm = norm(v);

    // Compute the elements of the normalized version of the vector.
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

//Matrix<double> get_matrix(int s, double alpha) {
//    Matrix<double> A(s, s);
//    for (int i = 0; i < s; ++i) {
//        for (int k = 0; k < s; ++k) {
//            if (k != i) {
//                A(i, k) = 0.2 / sqrt(static_cast<double>(i + k + 2))
//                          + 0.01 * sqrt(static_cast<double>(i + k + 2));
//            } else if (k == i) {
//                A(i, k) = 5.0 / (static_cast<double>(i + 1) + alpha);
//            }
//        }
//    }
//    return A;
//}

Matrix<double> givensturn(Matrix<double> A, int p, int q) {
    double cos_phi = A(p - 1, p) / sqrt(A(p - 1, p) * A(p - 1, p) + A(p - 1, q) * A(p - 1, q));
    double sin_phi = A(p - 1, q) / sqrt(A(p - 1, p) * A(p - 1, p) + A(p - 1, q) * A(p - 1, q));
    Matrix<double> O(A.getRows(), A.getCols());
    for (int i = 0; i < A.getCols(); ++i) {
        O(i, i) = 1.0;
    }
    O(p, p) = cos_phi;
    O(q, q) = cos_phi;
    O(q, p) = sin_phi;
    O(p, q) = -sin_phi;
    return O.transpose() * A * O;
}

Matrix<double> tridiag_with_givens(Matrix<double> A) {
    double eps{1.0e-10};
    for (int i = 1; i < A.getCols() - 1; ++i) {
        for (int j = i + 1; j < A.getCols(); ++j) {
            A = givensturn(A, i, j);
        }
    }
    for (int i = 0; i < A.getCols(); ++i) {
        for (int j = 0; j < A.getCols(); ++j) {
            if (A(i, j) < eps && A(i, j) > -1.0 * eps) {
                A(i, j) = 0.0;
            }
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

//double inverse_powermethod_with_varbias(const Matrix<double>& A, const Matrix<double> &x_0) {
//    double eps{1.0e-8};
//    bool rethist{false};
//    Matrix<double> x_temp, x_iter;
//    x_temp =  x_0;
//    int count = 0;
//
//    double lambda{0.0};
//    double old_lambda{1.0};
//
//    while (std::fabs(1.0 - lambda / old_lambda) > eps && count < 50) {
//        old_lambda = lambda;
//        x_iter = gauss(A, x_temp);
//        lambda = x_iter(0, 0);
//        x_iter /= lambda;
//        x_temp = x_iter;
//        count += 1;
//    }
//    return lambda;
//}

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
//    double eps{1.0e-8};
////    bool rethist{false};
//    Matrix<double> x_temp, x_iter;
//    x_temp = x_0;
//    int count = 0;
//
//    double lambda{0.0};
//    double old_lambda{1.0};
//
//    while (std::fabs(1.0 - lambda / old_lambda) > eps && count < 50) {
//        old_lambda = lambda;
//        x_iter = gauss(A, x_temp);
//        lambda = x_iter(0, 0);
//        x_iter /= lambda;
//        x_temp = x_iter;
//        count += 1;
//    }
    double temp{eigenVector(0, 0)};
    eigenVector = A * eigenVector;
//    std::cout << "A*v = \n" << eVector << "\n lambda = " << eVector(1,0) / v[1] << '\n';
    return eigenVector(0, 0) / temp;
}

//std::vector<double> initial_interval_net(const Matrix<double> &A) {
//    int s{A.getCols()};
//
//    double delta{0.0}; //opnorm(A, 2)
//    double lambda_max{powermethod(A)};
//    double lambda_min{inverse_powermethod(A)};
//    double a{lambda_min - delta};
//    double b{lambda_max + delta};
//
//    double step{(b - a) / (5 * s)};
//    std::vector<double> net;
//    net.reserve(5 * s + 2);
//    for (int i = 0; i < 5 * s + 2; ++i)
//        net.push_back(a + step * (i - 1));
//    return net;
//
//}

std::vector<double> generate_orthogonal(std::vector<std::vector<double>> vectors) {
    std::vector<double> approx(vectors[0].size(), 1.0);
    std::vector<double> orthogonal{approx};
    double s;
    vectors.push_back(approx);
    Matrix<double> q(vectors[0].size(), vectors.size()), r(vectors[0].size(), vectors.size());

    std::cout << vectors.size() << " - ";
    for (int i = 0; i < vectors.size(); ++i) {
//        std::cout << "You came to generate_orthogonal \n";
//        for (int j = 0; j < vectors[0].size(); ++j) {
//            std::cout << vectors[i][j] << '\n';
//        }
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
    int s{A.getCols()}, count{0};
    double lambda, delta{1e6};
    Matrix<double> x_iter, x_temp, tempMatrix(A.getRows(), A.getCols());
    Matrix<double> I = Matrix<double>::createIdentity(A.getCols());
    x_iter.pushBackColumn(eigenVector);
    while (delta > eps && count < 300) {
        x_temp = x_iter;
        I *= eigenValue;
        tempMatrix = A - I;
        x_iter = gauss(tempMatrix, x_iter);
        Normalize(x_iter);
//        std::cout << x_iter << '\n';
        x_temp *= x_iter - x_temp;
        delta = norm(x_iter - x_temp);

        count++;

//        x_iter = gauss(A, x_temp);
//        lambda = x_iter(0, 0);
//        x_iter /= lambda;
//        x_temp = x_iter;
    }

    for (int i = 0; i < A.getRows(); ++i) {
//        std::cout << "eigenVector[i] = " << x_iter(i,0) << '\n';
        eigenVector[i] = x_iter(i, 0);
    }

    return eigenVector;
}

Matrix<double> inverse_powermethod_modified(const Matrix<double> &A, const std::vector<double> &eigenVectorInit) {
    const double eps{1e-8};
    int count{0};
    double delta{1e6}, lambda, old_lambda;
    double diff{0.0}, new_diff{0.0};
    Matrix<double> eigenVector(A.getRows(), 1);
    Matrix<double> x_iter, x_temp;
    x_iter.pushBackColumn(eigenVectorInit);
    std::cout << '\n';
    while (delta > eps && count < 500) {
        x_temp = x_iter;
        old_lambda = lambda;
        x_iter = gauss(A, x_iter);
        Normalize(x_iter);
        delta = norm(x_iter - x_temp);
        lambda = getValueFromVector(A, x_iter);
        new_diff = std::fabs(lambda - old_lambda);
        std::cout << count << " & " << new_diff << " & " << getValueFromVector(A, x_iter) << "\\\\\n";
        std::cout << "\\hline \n";
        if (new_diff < eps || (new_diff > diff && count > 20)) {
            return x_iter;
        }
        diff = new_diff;
        count++;
    }

    return x_iter;
}

Matrix<double> inverse_powermethod_with_varbias(Matrix<double> A, double t_n, const std::vector<double> &eigenVectorInit) {
    double eps{1e-8}, delta{1e6};
    Matrix<double> x, x_n(A.getCols(), 1);
    x.pushBackColumn(eigenVectorInit);
    double diff{0.0}, new_diff{0.0};
    double old_t_n = t_n + 1.0;
    double lambda, old_lambda;
    int count = 0;
    Matrix<double> A_copy = A;

    Matrix<double> I = Matrix<double>::createIdentity(A.getCols());
    while (fabs(t_n / old_t_n - 1.0) > eps && count < 100) {
        for (int i = 0; i < A.getCols(); ++i) {
            A(i, i) -= t_n;
        }
        x_n = solveTridiagonal(A, x);
        old_t_n = t_n;
        lambda = x_n(0, 0);
//        lambda = getValueFromVector(A,x_n);
        t_n += 1 / lambda;
        x_n /= lambda;
        x = x_n;
        count += 1;
        Normalize(x_n);
        lambda = getValueFromVector(A_copy, x_n);
        new_diff = std::fabs(lambda - old_lambda);
        std::cout << count << " & " << new_diff << " & " << lambda << "\\\\\n";
        std::cout << "\\hline \n";
//        diff = new_diff;
        old_lambda = lambda;
//        delta = norm(x_n - x);
    }

    std::cout << "value = " << getValueFromVector(A_copy, x) << '\n';
    return x;
}

std::vector<double> find_them_all(const Matrix<double> &A, const double &first) {
    std::vector<double> x_0(A.getCols(), 1.0);
    std::vector<double> eigenValues, eigenVector, orthogonal(A.getCols());
    Matrix<double> eigenVectorMatrix(A.getRows(), A.getCols()), zeros(A.getRows(), 1), initApproximation(A.getRows(),
                                                                                                         1);
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(-10.0, 100.0);

    std::vector<std::vector<double>> eigenVectors;
    eigenValues.reserve(A.getCols());
    eigenValues.push_back(first);
    eigenVector = inverse_powermethod_modified(A, first);
    eigenVectors.push_back(eigenVector);
    for (int i = 1; i < A.getRows() - 1; ++i) {
//        for (int j = 0; j < A.getRows(); ++j) {
//            eigenVectorMatrix(i - 1, j) = eigenVector[j];
//        }
//        std::cout << "eigenVectorMatrix : \n" << eigenVectorMatrix << '\n';
//        eigenVector = generate_orthogonal(eigenVectors);                           ////TODO
        eigenVectors.push_back(generate_orthogonal(eigenVectors));
//        for (int j = 0; j < A.getRows(); ++j) {
//            initApproximation(j, 0) = eigenVector[j];
//        }
//        Normalize(initApproximation);
//        initApproximation = inverse_powermethod_modified(A, eigenVectors[i]);
        Matrix<double> temp;
        temp.pushBackColumn(eigenVectors[i]);
        initApproximation = inverse_powermethod_with_varbias(A, getValueFromVector(A, temp), eigenVectors[i]);
        std::cout << "Pushed value: " << getValueFromVector(A, initApproximation) << '\n';
//        std::cout << "Pushed eigenVector: \n";
        for (int j = 0; j < A.getRows(); ++j) {
            eigenVectors[i][j] = initApproximation(j, 0);
//            std::cout << eigenVectors[i][j] << '\n';
        }
//        std::cout << '\n';
//        std::cout << "initApproximation : \n" << initApproximation << '\n';
//        for (int j = 0; j < A.getRows(); ++j) {
//            initApproximation(j, 0) = -initApproximation(j, 0);             ////WHY?
//            eigenVector[j] = initApproximation(j, 0);
//        }
//        std::cout << "Pushed value: " << getValueFromVector(A, initApproximation) << '\n';
        eigenValues.push_back(getValueFromVector(A, initApproximation));

    }
    eigenValues.push_back(1e5);
    std::sort(eigenValues.begin(), eigenValues.end());
    return eigenValues;
}

using std::cout;
using std::endl;

int main() {
    Matrix<double> A = get_matrix(20);
    std::cout << "Matrix:\n";
    std::cout << std::fixed << A << '\n';
    std::cout << "||||||||||||||||||||||||||||||||\n\n";
    std::cout << "Tridiagonal form of matrix:\n";
    A = tridiag_with_givens(A);
    std::cout << std::fixed << A << '\n';
    std::cout << "||||||||||||||||||||||||||||||||\n\n";
    std::cout << "Eigenvalues:\n";
    std::cout << "lambda_max = " << std::setprecision(16) << std::fixed << powermethod(A) << '\t'
              << "lambda_min = " << std::setprecision(16) << std::fixed << inverse_powermethod(A) << '\n';
    double lambda_max = powermethod(A);
    std::cout << "||||||||||||||||||||||||||||||||\n\n";
    std::vector<double> values = find_them_all(A, inverse_powermethod(A));
    values.back() = lambda_max;
    std::cout << "Values: \n";
    for (auto i : values) {
        std::cout << i << '\n';
    }

//    Matrix<double> eVector;
//    std::vector<double> v{-0.00704204616,
//                          0.0207184355,
//                          -0.151905104,
//                          0.302183775,
//                          -0.417725541,
//                          0.479727084,
//                          -0.474703335,
//                          0.401082565,
//                          -0.27593617,
//                          0.13471083};
//    v = generate_orthogonal(v);
//    eVector = inverse_powermethod_modified(A, v);
//    eVector.pushBackColumn(v);
//    eVector = A * eVector;
//    std::cout << "vector = \n" << eVector << '\n';

}
//    /*
//     * TODO:
//     *       1) Найти следующее собственное число
//     *       2) Дописать отчёт*/

/* TEST */

//    Eigen::MatrixXd A(10,10);
//
//    A << 1.250000, 0.321084, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
//         0.321084, 1.402853, 0.134806, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
//         0.000000, 0.134806, 0.651579, 0.164810, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
//         0.000000, 0.000000, 0.164810, 0.591594, 0.147146, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
//         0.000000, 0.000000, 0.000000, 0.147146, 0.558219, 0.127510, 0.000000, 0.000000, 0.000000, 0.000000,
//         0.000000, 0.000000, 0.000000, 0.000000, 0.127510, 0.522435, 0.107190, 0.000000, 0.000000, 0.000000,
//         0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.107190, 0.486946, 0.086736, 0.000000, 0.000000,
//         0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.086736, 0.453409, 0.066012, 0.000000,
//         0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.066012, 0.422525, 0.043501,
//         0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.043501, 0.394443;
//    cout << "Here is a 10x10 matrix, A:" << endl << A << endl << endl;
//    Eigen::EigenSolver<Eigen::MatrixXd> es(A);
//    cout << "The eigenvalues of A are:" << endl << std::setprecision(9) << es.eigenvalues().real() << endl;
//    cout << "The matrix of eigenvectors, V, is:" << endl << std::setprecision(9) << es.eigenvectors().real() << endl << endl;
//
//    std::complex<double> lambda = es.eigenvalues()[0];
//    cout << "Consider the first eigenvalue, lambda = " << lambda.real() << endl;
//    Eigen::VectorXcd v = es.eigenvectors().col(0);
//    cout << "If v is the corresponding eigenvector, then lambda * v = " << endl << lambda * v << endl;
//    cout << "... and A * v = " << endl << A.cast<std::complex<double> >() * v << endl << endl;
//
//    Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
//    Eigen::MatrixXcd V = es.eigenvectors();
//    cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;
//}

/*
 *
1.6678972
1.01871523
0.821207147
0.681495297
0.578181502
0.305337326
0.498798237
0.341903097
0.384613524
0.435854438

1.6678972       1.01871523      0.821207147     0.681495297   0.578181502   0.305337326   0.498798237   0.341903097    0.384613524    0.435854438

0.605660162    -0.778333243     0.142430466    -0.0656488895  0.0385370767 -0.00704204616 0.0254576042 -0.00980685621 -0.0132028503  -0.017989998
0.788278725     0.560652735    -0.190209309     0.116236569  -0.0806328593  0.0207184355 -0.0595601062  0.027735968    0.0355843584   0.0456157174
0.107272056     0.256238311     0.481449728    -0.46562717    0.401480514  -0.151905104   0.338794774  -0.19492937    -0.237334944   -0.284364439
0.0167334603    0.112220229     0.65110601     -0.179596065  -0.112844319   0.302183775  -0.265349582   0.34358324     0.355337973    0.334901533
0.002247898     0.038743876     0.476769811     0.411795732  -0.439389989  -0.417725541  -0.21212571   -0.36469425    -0.234004667   -0.0359596235
0.000252385425  0.0104199765    0.231959566     0.605376013   0.0614325678  0.479727084   0.40506471    0.22219643    -0.0914600396  -0.351966413
2.30291142e-05  0.00215502136   0.0793930327    0.408463609   0.554634464  -0.474703335   0.163016425   0.0596008902   0.39596131     0.327070262
1.64944508e-06  0.00033500249   0.0193029456    0.168050787   0.507487048   0.401082565  -0.478310932  -0.374261223   -0.354133228    0.242307108
8.75346695e-08  3.72820537e-05  0.00323204191   0.0439553727  0.230468001  -0.27593617   -0.543076422   0.553881006   -0.151204873   -0.49418835
2.99017087e-09  2.59791567e-06  0.000329449079  0.00666116484 0.0545644403  0.13471083   -0.226384111  -0.458591967    0.669167206   -0.519124386

-2.4542955810107614
-2.4542956194392307  39.1383251756 11.6083794364 -2.4542955422 0.8811279103 0.8264630199
-2.4542955422
0.8264630156152148
0.8264630156152148
0.8264630199
0.8811279118009789
0.8811279118009772
0.8811279103
11.6083794731384806
11.6083794377122214
11.6083794364
39.1383251733118129
39.1383251733118129
39.1383251756
 * */