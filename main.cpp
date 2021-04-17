#include "Matrix.h"
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

//A*x=y
Matrix<double> gauss(Matrix<double> A, Matrix<double> b) {
    const double eps{1e-16};  // точность
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

//Matrix<double> get_matrix(int s) {
//    Matrix<double> A(s, s);
//    for (int i = 0; i < s; ++i) {
//        for (int k = 0; k < s; ++k) {
//            A(i, k) = pow((static_cast<double>(std::abs(i - k) + 1)), -3.0)
//                      + pow(static_cast<double>(i + k - s + 2), 2.0);
//        }
//    }
//    return A;
//}

Matrix<double> get_matrix(int s, double alpha) {
    Matrix<double> A(s, s);
    for (int i = 0; i < s; ++i) {
        for (int k = 0; k < s; ++k) {
            if (k != i) {
                A(i, k) = 0.2 / sqrt(static_cast<double>(i + k + 2))
                          + 0.01 * sqrt(static_cast<double>(i + k + 2));
            } else if (k == i) {
                A(i, k) = 5.0 / (static_cast<double>(i + 1) + alpha);
            }
        }
    }
    return A;
}

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

double inverse_powermethod_with_varbias(Matrix<double> A, double t_0, const Matrix<double> &x_0) {
    double eps{1e-8};
//    bool rethist{false};
    Matrix<double> x_init, x_temp, x_iter;
    x_init = x_0;
    x_temp = x_init;
    double old_t{t_0 + 1.0};
    double t = t_0;
    int count = 0;

    std::vector<double> t_history, c_history;
    t_history.push_back(t_0);
    c_history.push_back(0);

    while (std::fabs(static_cast<double>(t) / static_cast<double>(old_t) - 1.0) > eps) {
        std::cout << "\t " << count << '\t' << t << '\n';
        old_t = t;
        //создаём единичную диагональную матрицу размером A и покомпонентно умножаем на t
        Matrix<double> t_matrix(A.getRows(), A.getCols());
        for (int i = 0; i < A.getCols(); ++i) {
            t_matrix(i, i) = t;
        }
        A = A - t_matrix;
        x_iter = gauss(A, x_temp);
        double lambda = x_iter(0, 0);
        t += 1 / lambda;
        x_iter /= lambda;
        x_temp = x_iter;
        count += 1;

        t_history.push_back(t);
        c_history.push_back(count);
    }

//    if (rethist) {
//        puts("Implement printing history?");
//        return t;
//    } else
    return t;
}

std::vector<double> initial_interval_net(const Matrix<double> &A) {
    int s{A.getCols()};

    double delta{0.0}; //opnorm(A, 2)
    double lambda_max{powermethod(A)};
    double lambda_min{inverse_powermethod(A)};
    double a{lambda_min - delta};
    double b{lambda_max + delta};

    double step{(b - a) / (5 * s)};
    std::vector<double> net;
    net.reserve(5 * s + 2);
    for (int i = 0; i < 5 * s + 2; ++i)
        net.push_back(a + step * (i - 1));
    return net;

}

double inverse_powermethod_modified(Matrix<double> A, double t) {
    std::vector<double> x_init(A.getCols(), 1.0);
    int s = A.getCols();
    double lambda;
    Matrix<double> x_iter, x_temp;
    x_temp.pushBackColumn(x_init);
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < s; ++j) {
            A(j, j) -= t;
        }
        x_iter = gauss(A, x_temp);
        lambda = x_iter(0, 0);
        x_iter /= lambda;
        x_temp = x_iter;

    }
    double t_0 = t;
    lambda = inverse_powermethod_with_varbias(A, t_0, x_temp);

    return lambda;
}

std::vector<double> find_them_all(const Matrix<double> &A) {
    std::vector<double> x_0(A.getCols(), 1.0);
//    bool rethist{false};
    std::vector<double> net = initial_interval_net(A);
    std::vector<double> found;
    found.reserve(net.size());
    for (double i : net) {
        std::cout << i << " : finding\n";
        found.push_back(inverse_powermethod_modified(A, i));
    }
//    if (rethist)
//        return found;
//    else {
    std::sort(found.begin(), found.end());
    return found;
//    }
}

int main() {
    Matrix<double> A = get_matrix(10, 3);
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
    std::cout << "||||||||||||||||||||||||||||||||\n\n";
    std::vector<double> eigen = find_them_all(A);
    for (auto i : eigen) {
        std::cout << i << '\n';
    }
    /*
     * TODO:
     *       1) Понять в чём ошибка
     *       2) Начать писать моё задание
     *       3) Дописать отчёт*/

}