#include <iostream>
#include <cmath>
using namespace std;
const int MAX_NODES = 10;
//вычисления коэффициентов кубического сплайна
void Splinecb(double x[], double y[], int n, double a[], double b[], double c[], double d[]) {
    double h[MAX_NODES]; // Шаги между узлами
    double alpha[MAX_NODES]; //альфа-коэффициенты
    double l[MAX_NODES], mu[MAX_NODES], z[MAX_NODES]; //лемма для сплайнов
    for (int i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }
    //альфа-коэффициенты
    for (int i = 1; i < n - 1; i++) {
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
    }
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;
    for (int i = 1; i < n - 1; i++) {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    l[n - 1] = 1.0;
    z[n - 1] = 0.0;
    c[n - 1] = 0.0;
    for (int j = n - 2; j >= 0; j--) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / h[j];
    }
    for (int i = 0; i < n - 1; i++) {
        a[i] = y[i];
    }
}
//функция для интерполяции кубических сплайнов
double SplineInterpol(double x[], double y[], double a[], double b[], double c[], double d[], int n, double x_value) {
    int i = 0;
    while (i < n - 1 && x_value > x[i + 1]) {
        i++;
    }
    double dx = x_value - x[i];
    return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
}
int main() {
    setlocale(LC_ALL, "Russian");
    double x[MAX_NODES], y[MAX_NODES];
    double a[MAX_NODES], b[MAX_NODES], c[MAX_NODES], d[MAX_NODES]; // Коэффициенты сплайна
    int n;
    cout << "Введите количество узлов (максимум " << MAX_NODES << "): ";
    cin >> n;
    if (n < 2 || n > MAX_NODES) {
        cout << "Ошибка: количество узлов должно быть от 2 до " << MAX_NODES << "." << endl;
        return 1;
    }
    cout << "Введите узлы (x):" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = ";
        cin >> x[i];
    }
    cout << "Введите значения функции (y):" << endl;
    for (int i = 0; i < n; i++) {
        cout << "y[" << i << "] = ";
        cin >> y[i];
    }
    //коэффициент кубического сплайна
    Splinecb(x, y, n, a, b, c, d);
    double x_value;
    cout << "Введите значение x для интерполяции: ";
    cin >> x_value;
    double result = SplineInterpol(x, y, a, b, c, d, n, x_value);
    cout << "Интерполированное значение в x = " << x_value << " равно " << result << endl;
    return 0;
}