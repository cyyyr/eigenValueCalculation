#include "Matrix.h"
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>

Matrix<double> get_matrix(int s) {
    Matrix<double> A(s, s);
    for (int i = 0; i < s; ++i) {
        for (int k = 0; k < s; ++k) {
            A(i, k) = pow((static_cast<double>(std::abs(i - k) + 1)), -3.0)
                      + pow(static_cast<double>(i + k - s), 2.0);
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
    for (int i = 1; i < A.getCols() - 1; ++i) {
        for (int j = i + 1; j < A.getCols(); ++j) {
            A = givensturn(A, i, j);
        }
    }

    return A;
}

double powermethod(const Matrix<double> &A) {
    std::vector<double> x_0(A.getCols(), 1.0);
    double eps{1e-16};
    std::vector<double> x_n = x_0;
    double lambda{0.0};
    double old_lambda{1.0};

//    while (std::fabs(1.0 - lambda / old_lambda) > eps) {
//
//        old_lambda = lambda;
//        x_np1 = A * x_n;
//        lambda = x_np1[0];
//        x_n = x_np1 / lambda;
//    }
    return lambda;
}

double inverse_powermethod(const Matrix<double> &A) {
    std::vector<double> x_0(A.getCols(), 1.0);
    double eps{1e-16};
    std::vector<double> x_n = x_0;
    double lambda{0.0};
    double old_lambda{1.0};

//    while (std::fabs(1.0 - lambda / old_lambda) > eps) {
//
//        old_lambda = lambda;
//        x_np1 = A / x_n;
//        lambda = x_np1[0];
//        x_n = x_np1 / lambda;
//    }

    return 1.0 / lambda;

}

double inverse_powermethod_with_varbias(Matrix<double> A, double t_0, std::vector<double> x_0) {
    double eps{1e-16};
    bool rethist{false};

    std::vector<double> x_n = x_0;
    double old_t{t_0 + 1.0};
    double t = t_0;
    int count = 0;

    std::vector<double> t_history, c_history;
    t_history.push_back(t_0);
    c_history.push_back(0);

    while (std::fabs(static_cast<double>(t) / static_cast<double>(old_t) - 1.0) > eps) {

        old_t = t;
        //создаём диагональную матрицу размером A и покомпонентно умножаем на t
        Matrix<double> t_matrix(A.getRows(), A.getCols());
        for (int i = 0; i < A.getCols(); ++i) {
            t_matrix(i, i) = t;
        }
//        A -= t_matrix;
//        A.pushBackColumn(x_n);
//        Matrix<double> x_np1 = A.gaussianEliminate();
//        double lambda = x_np1[0];
//        t += 1 / lambda;
//        x_n = x_np1 / lambda;
        count += 1;

        t_history.push_back(t);
        c_history.push_back(count);
    }

    if (rethist) {
        puts("Implement printing history?");
        return t;
    } else
        return t;
}

std::vector<double> initial_interval_net(const Matrix<double> &A) {
    int s{A.getCols()};
    std::vector<double> x(A.getCols(), 1.0);
    double eps{1e-16};

    double delta{0.0}; //opnorm(A, 2)
    double lambda_max{powermethod(A)};
    double lambda_min{inverse_powermethod(A)};
    double a{lambda_min - delta};
    double b{lambda_max + delta};

    double step{(b - a) / (5 * s)};
    std::vector<double> net;
    for (int i = 0; i < 5 * s + 2; ++i)
        net.push_back(a + step * (i - 1));
    return net;

}

double inverse_powermethod_modified(Matrix<double> A, double t) {
    std::vector<double> x_0(A.getCols(), 1.0);
    double eps{1e-16};
    bool rethist{false};
    int s = A.getCols();
    std::vector<double> x_n{x_0};
    double lambda{0.0};
//    for (int i = 0; i < 10; ++i) {
//        x_np1 = (A .- t .* Diagonal(ones(Float64, s))) \ x_n;
//        lambda = x_np1[0];
//        x_n = x_np1 / lambda;
//
//    }
    double t_0 = t;
    lambda = inverse_powermethod_with_varbias(A, t_0, x_n);

    return lambda;
}

std::vector<double> find_them_all(const Matrix<double> &A) {
    std::vector<double> x_0(A.getCols(), 1.0);
    double eps{1e-16};
    bool rethist{false};
    int s = A.getCols();
    std::vector<double> net = initial_interval_net(A);
    std::vector<double> found;
    for (double i : net)
        found.push_back(inverse_powermethod_modified(A, i));
    if (rethist)
        return found;
    else {
        std::sort(found.begin(), found.end());
        return found;
    }
}

int main() {
    Matrix<double> A = get_matrix(3);
    std::vector<double> v{1.0, 2.0, 1.0};
    Matrix<double> b;
    b.pushBackColumn(v);
    std::cout << std::fixed << A << '\n';
    std::cout << std::fixed << b << '\n';
    Matrix C = Matrix<double>::solve(A, b);
    std::cout << std::fixed << C << '\n';
//    std::cout << std::fixed << tridiag_with_givens(A) << '\n';
}
/* TODO:
 *      1) Разобраться с решением СЛАУ
 *      2) Умножение на вектор или переделать все вектора в матрицы
 *      3) Получить похожие результаты
 *          3а) Взять его матрицу
 *      4) Сделать моё задание
 */