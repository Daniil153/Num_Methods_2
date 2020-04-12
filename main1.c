#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
std::vector <double> masx;
std::vector <double> masy;
std::vector <double> masx1;
std::vector <double> masy1;
std::vector <double> masansx;
std::vector <double> masansy;

double f0(double x){
    return 4 - x - 4 * exp(-x);
}

double f1(double x, double y){
    return 3 - y - x;
}

double f2(double x, double u, double v){
    return log(2 * x + sqrt(4 * x * x + v * v));
}

double f3(double x, double u, double v){
    return sqrt(4 * x * x + u * u);
}

double RK_2(double x, double y0, double low, double high, double n)
{
    if (low < x) {
        return 0;
    }
    double h = (high - low) / n;
    double k;
    double dy;
    for (int i = 1; i <= n; i++) {
        masx.push_back(x);
        masy.push_back(y0);
        printf("(%lf;%lf)\n", x, y0);
        k = f1(x, y0);
        dy = h * f1(x + h / 2, y0 + (h * k) / 2 );
        y0 = y0 + dy;
        x += h;
    }
    printf("(%lf;%lf)\n", x, y0);
    return 0;
}

double RK_4(double x, double y0, double low, double high, double n)
{
    if (low < x) {
        return 0;
    }
    double h = (high - low) / n;
    double k1, k2 ,k3, k4;
    double dy;
    for (int i = 1; i <= n; i++) {
        masx.push_back(x);
        masy.push_back(y0);
        printf("(%lf;%lf)\n", x, y0);
        k1 = f1(x, y0);
        k2 = f1(x + h / 2, y0 + (h * k1) / 2);
        k3 = f1(x + h / 2, y0 + (h * k2) / 2);
        k4 = f1(x + h, y0 + h * k3);
        dy = h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        y0 = y0 + dy;
        x += h;
    }
    printf("(%lf;%lf)\n", x, y0);
    return 0;
}

double RK_system_2(double x, double y0, double z0, double low, double high, double n)
{
    if (low < x) {
        return 0;
    }
    double h = (high - low) / n;
    double k, l;
    double dy, dz;
    for (int i = 1; i <= n; i++) {
        masx.push_back(x);
        masy.push_back(y0);
        masx1.push_back(x);
        masy1.push_back(z0);
        printf("(%lf;%lf)\n ", y0, z0);
        k = f2(x, y0, z0);
        l = f3(x, y0, z0);
        dy = h * f2(x + h / 2, y0 + (h * k) / 2, z0 + (h * l) / 2);
        dz = h * f3(x + h / 2, y0 + (h * k) / 2, z0 + (h * l) / 2);
        y0 = y0 + dy;
        z0 = z0 + dz;
        x += h;
    }
    printf("(%lf;%lf)\n ", y0, z0);
    return 0;
}

double RK_system_4(double x, double y0, double z0, double low, double high, double n)
{
    if (low < x) {
        return 0;
    }
    double h = (high - low) / n;
    double k1, k2 ,k3, k4;
    double l1, l2 ,l3, l4;
    double dy, dz;
    for (int i = 1; i <= n; i++) {
        masx.push_back(x);
        masy.push_back(y0);
        masx1.push_back(x);
        masy1.push_back(z0);
        printf("(%lf;%lf)\n", y0, z0);
        k1 = f2(x, y0, z0);
        l1 = f3(x, y0, z0);
        k2 = f2(x + h / 2, y0 + (h * k1) / 2, z0 + (h * l1) / 2);
        l2 = f3(x + h / 2, y0 + (h * k1) / 2, z0 + (h * l1) / 2);
        k3 = f2(x + h / 2, y0 + (h * k2) / 2, z0 + (h * l2) / 2);
        l3 = f3(x + h / 2, y0 + (h * k2) / 2, z0 + (h * l2) / 2);
        k4 = f2(x + h, y0 + h * k3, z0 + h * l3);
        l4 = f3(x + h, y0 + h * k3, z0 + h * l3);
        dy = h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        dz = h * (l1 + 2 * l2 + 2 * l3 + l4) / 6;
        y0 = y0 + dy;
        z0 = z0 + dz;
        x += h;
    }
    printf("(%lf;%lf)\n", y0, z0);
    return 0;
}
//Для смены теста поменяем функции.
double p1(double x) { return 1/x; }
double q1(double x) { return 2; }
double f1(double x) { return x; }

void sol(double A, double a, double a0, double a1, double B, double b, double b0, double b1, int iter) {
    double x_i = a;
    double h = (b - a) / iter;
    std::vector <double> alpha(iter + 5, 0.0), betta(iter + 5, 0.0), y(iter + 5, 0.0);
    alpha[1] = -a1 / (a0 * h - a1);
    betta[1] = (A * h) / (a0 * h - a1);
    double A_i, C_i, B_i, D_i;
    for (int i = 1; i < iter; ++i) {
        x_i += h;
        A_i = (1.0 / (h * h)) - (p1(x_i) / (2.0 * h));
        C_i = (1.0 / (h * h)) + (p1(x_i) / (2.0 * h));
        B_i = -2.0 / (h * h) + q1(x_i);
        D_i = f1(x_i);
        alpha[i + 1] = C_i / (-B_i - A_i * alpha[i]);
        betta[i + 1] = (A_i * betta[i] - D_i) / (-B_i - A_i * alpha[i]);
    }
    x_i = a;
    y[iter] = (b1 * betta[iter] + B * h) / (b1 * (1.0 - alpha[iter]) + b0 * h);
    for (int i = iter; i > 0; --i) {
        y[i - 1] = y[i] * alpha[i] + betta[i];
    }
    for (int i = 0; i <= iter; ++i) {
        masx.push_back(x_i);
        masy.push_back(y[i]);
        printf("x_i = %.5lf  y_i = %.5lf\n", x_i, y[i]);
        x_i += h;
    }
}

int main(int argc, char **argv) {
    int problem, low, high, n;
    printf("Введите интервал [a;b]\n");
    scanf("%d%d", &low, &high);
    double a1 = low;
    double shag = 0.01;
    while (a1 < high) {
        masansx.push_back(a1);
        masansy.push_back(f0(a1));
        a1 += shag;
    }
    printf("Введите разбиение отрезка n\n");
    scanf("%d", &n);
    printf("Если вы хотите решить задачу Коши - введите 1.\n Если вы хотите решать краевую задачу - введите 2\n");
    scanf("%d", &problem);
    int mode, type;
    if (problem == 1) {
        printf("Введите 0, если вы хотите решать одно уравнение.\n Введите 1, если вы хотите решать систему уравнений:\n");
        scanf("%d", &mode);
        printf("Введите 2, если вы хотите использовать метод порядка точности 2.\n Введите 4, если хотите использовать метод порядка точности 4:\n");
        scanf("%d", &type);
        if (!mode) {
            if (type == 2) {
                RK_4(0, 0, low, high, n);
                plt::plot(masx, masy, "-o");
                plt::plot(masansx, masansy, "y");
                masx.clear();
                masy.clear();
                RK_2(0, 0, low, high, n);
                plt::plot(masx, masy, "-o");
                plt::plot(masansx, masansy, "y");
                plt::show();
            }
        } else if (mode == 1) {
            if (type == 2) {
                RK_system_2(0, 0.5, 1, low, high, n);
                plt::plot(masx, masy, "-ob");
                plt::plot(masx1, masy1, "-*r");
                masx.clear();
                masy.clear();
                masx1.clear();
                masy1.clear();
                RK_system_4(0, 0.5, 1, low, high, n);
                plt::plot(masx, masy, "-oy");
                plt::plot(masx1, masy1, "-*g");
                plt::show();
            }
        }
    } else if (problem == 2) {
        //Для смены теста - меняем коэффициенты
        sol(0.5, 0.7, 1, 0, 1.2, 1, 2, 3, n);
        plt::plot(masx, masy);
        plt::show();
    }
    return 0;
}

