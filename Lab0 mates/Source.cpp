#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <Eigen/Dense>
#include <chrono>

using namespace std;
using namespace Eigen;


int main() {

    srand(static_cast<unsigned>(time(nullptr)));

    int sz;

    cout << "Choose your size: ";
    cin >> sz;

    int** Mat = new int* [sz];

    for (int i = 0; i < sz; i++) {
        Mat[i] = new int[sz + 1];
        for (int j = 0; j < sz + 1; j++) {
            Mat[i][j] = rand() % 100;
        }
    }

    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz + 1; j++) {
            std::cout << Mat[i][j] << " ";
        }
        std::cout << std::endl;
    }
    auto start = chrono::high_resolution_clock::now();
    int typeGauss = 0;
    cout << "Choose what type of loop you want (1,2,3): " << endl;
    cin >> typeGauss;

    if (typeGauss == 3) {
        for (int Mati = 0; Mati <= sz - 1; Mati++) {

            for (int Matj = Mati + 1; Matj <= sz - 1; Matj++) {
                double fact = Mat[Matj][Mati] / Mat[Matj][Matj];

                for (int Matk = Mati; Matk <= sz - 1; Matk++) {
                    Mat[Matj][Matk] -= fact * Mat[Mati][Matk];
                }
            }

        }
        double* solucion = new double[sz];
        for (int i = sz - 1; i >= 0; i--) {
            solucion[i] = Mat[i][sz];
            for (int j = i + 1; j < sz; j++) {
                solucion[i] -= Mat[i][j] * solucion[j];
            }
            solucion[i] /= Mat[i][i];
        }

        cout << "Solucion del sistema:" << endl;
        for (int i = 0; i < sz; i++) {
            cout << "x" << i + 1 << " = " << solucion[i] << endl;
        }

        delete[] solucion;

        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "duration: " << duration.count() << " milliseconds" << endl;

        for (int i = 0; i < sz; i++) {
            delete[] Mat[i];
        }
        delete[] Mat;

        return 0;
    }
    if (typeGauss == 2) {

        MatrixXd mat(sz, sz);
        for (int i = 0; i < sz; i++) {
            for (int j = 0; j < sz; j++) {
                mat(i, j) = i * sz + j + 1;
            }
        }
        double* solucion = new double[sz];
        for (int i = sz - 1; i >= 0; i--) {
            solucion[i] = Mat[i][sz];
            for (int j = i + 1; j < sz; j++) {
                solucion[i] -= Mat[i][j] * solucion[j];
            }
            solucion[i] /= Mat[i][i];
        }

        cout << "Solucion del sistema:" << endl;
        for (int i = 0; i < sz; i++) {
            cout << "x" << i + 1 << " = " << solucion[i] << endl;
        }

        delete[] solucion;

        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "duration: " << duration.count() << " milliseconds" << endl;
        for (int i = 0; i < sz; i++) {
            delete[] Mat[i];
        }
        delete[] Mat;

        return 0;

    }

    if (typeGauss == 1) {

        MatrixXd mat(sz, sz);
        Eigen::VectorXd b(sz);
        for (int col = 0; col < sz - 1; col++) {
            // Encontrar el elemento máximo en la columna actual
            int max_row = mat.block(col, col, sz - col, 1).cwiseAbs().maxCoeff() + col;
            int max_val = -1;
            
            // Intercambiar filas si es necesario
            if (max_row != col) {
                mat.row(max_row).swap(mat.row(col));
                swap(b(max_row), b(col));
            }

            // Realizar la eliminación de Gauss en la columna actual
            VectorXd ratios = mat.block(col + 1, col, sz - col - 1, 1) / mat(col, col);
            mat.block(col + 1, col, sz - col - 1, sz - col) -= ratios * mat.block(col, col, 1, sz - col);
            b.segment(col + 1, sz - col - 1) -= ratios * b(col);
        }
        double* solucion = new double[sz];
        for (int i = sz - 1; i >= 0; i--) {
            solucion[i] = Mat[i][sz];
            for (int j = i + 1; j < sz; j++) {
                solucion[i] -= Mat[i][j] * solucion[j];
            }
            solucion[i] /= Mat[i][i];
        }

        cout << "Solucion del sistema:" << endl;
        for (int i = 0; i < sz; i++) {
            cout << "x" << i + 1 << " = " << solucion[i] << endl;
        }

        delete[] solucion;

        for (int i = 0; i < sz; i++) {
            delete[] Mat[i];
        }
        delete[] Mat;

        return 0;
    }
}

