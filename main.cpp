#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> // ��� setw()
#include "Problems.h"

/*
� ������ ��������� ����������� ��������� ������� ������ ������� ��� ��������� ���������������� � ��������������.

������� �������:
��� x = xMin u(t, x, y) = mu1(t, y)
��� x = xMax u(t, x, y) = mu2(t, y)
��� y = yMin u(t, x, y) = mu3(t, x)
��� y = yMax u(t, x, y) = mu4(t, x)

��������� �������: u(0, x, y) = phi(x, y).

������� ����� �������-��������. ����� �������� ���������� ����������, ������� ����� �������� ����� ���� �� ����������� � �������

������� �� ���� t_n ������������ � ���� ���������� ������� w_n[i][j], ��� i �������� ����� ��� Y, � j �������� ����� ��� X
*/

/* ��������� �� ������� ����������� ��� ���������� ����������, ����� �� ������ �������� ���������
   ��� �����������, ������� ���������� ��� ������� */
double (*uExact)(double, double, double);
double (*f)(double, double, double);
double (*phi)(double, double);
double (*mu_xMin)(double, double, double); // �������� ������� �� ������� x = xMin
double (*mu_xMax)(double, double, double);
double (*mu_yMin)(double, double, double);
double (*mu_yMax)(double, double, double);

// ��������� �����. ������� g_ ���������� ���������� ����������
double g_xMin, g_xMax, g_yMin, g_yMax;        // ������� ��������������
int g_Nx, g_Ny;                               // ����� ����� ����� �� ��� X � ��� Y ��������������
int g_N;                                      // ����� ����� �����
double g_hx, g_hy;                            // ��� �� ��� X � ��� Y ��������������
double g_hx_sq, g_hy_sq;                      // ���� �� ��� X � ��� Y, ���������� � �������
double g_tau;                                 // ��� �� �������
double g_T;                                   // ������������ �������� ��������� t

void assignW_np12OnBorders(double** w_np12, double t);                            // ��������� ������� w_np12 �������� ��������� ������� � ������ ������� t
void assignWOnBorders(double** w, double t);                                      // ��������� ������� w_n �������� ��������� ������� � ������ ������� t
void initW(double** w);                                                           // ���������������� ������ ������� w_n � ������ ������� t = 0
void gridInfoToFile(std::string filename);                                        // �������� � ��������� ����� ���������� � �����
void writeStepResultsToStream(double** u, double** w_n, double** res,             // �������� � ����� ������ ostream ������� ������� ������� u, ������������� ������� w_n
	double maxRes, int N_rows, int N_cols, double t, std::ostream& streamOut);    // � ������� ����� ���� res � ������ ������� t
void getExactSolution(double** ar, double t);                                     // �������� � ������ ar �������� ������� �������
void getRes(double** res, double** a1, double** a2, int N_rows, int N_cols);      // �������� � ������ diffAr ������ �������� �������� a1 � a2
double getMax(double** ar, int N_rows, int N_cols);                               // ����� ����� ������� ����� � ������� ar
double getMaxInRow(double** ar, int N_cols, int rowId);

int main() {
	// ������� ������
	// ������� 1
	uExact = uExact1;
	f = f1;
	phi = phi1;
	mu_xMin = mu1;
	mu_xMax = mu1;
	mu_yMin = mu1;
	mu_yMax = mu1;
	g_xMin = 0;
	g_xMax = M_PI;
	g_yMin = 0;
	g_yMax = M_PI;

	// ������� 2
	/*uExact = uExact2;
	f = f2;
	phi = phi2;
	mu_xMin = mu_xMin2;
	mu_xMax = mu_xMax2;
	mu_yMin = mu_yMin2;
	mu_yMax = mu_yMax2;
	g_xMin = 3;
	g_xMax = 6;
	g_yMin = -2;
	g_yMax = 2;*/

	// ������������� ���������� �����
	g_Nx = 30;
	g_Ny = 40;
	g_N = g_Nx * g_Ny;
	g_hx = (g_xMax - g_xMin) / (g_Nx - 1.0);
	g_hy = (g_yMax - g_yMin) / (g_Ny - 1.0);
	g_hx_sq = g_hx * g_hx;
	g_hy_sq = g_hy * g_hy;
	g_tau = 0.05;
	g_T = 10;

	// ���������� � ������������� �������� w_n, w_np12, u � res � ������ ������� 0
	double** w_n = new double*[g_Ny];        // ������ ��� ������������� �������
	double** u = new double*[g_Ny];          // ������ ��� ������� �������
	double** res = new double*[g_Ny];        // ������ ��� ������� ������� � ������������� �������
	double** w_np12 = new double*[g_Ny];     // ������������ ������� �� ��������� ����. np12 ���������������� ��� n plus 1/2
	for (int i = 0; i < g_Ny; i++) {
		w_n[i] = new double[g_Nx];
		u[i] = new double[g_Nx];
		res[i] = new double[g_Nx];
		w_np12[i] = new double[g_Nx];
	}

	initW(w_n);
	getExactSolution(u, 0);
	getRes(res, w_n, u, g_Ny, g_Nx);

	// ����� � ���� ����� � ����������� ����������� � ������� ��� t = 0
	std::string filename = "Results.txt";
	gridInfoToFile(filename);
	std::ofstream outstream(filename, std::ios::app);
	outstream << std::fixed;
	writeStepResultsToStream(u, w_n, res, 0, g_Ny, g_Nx, 0, outstream);

	// ������� ������������� ��� ������ ��������
	double* alpha1 = new double[g_Nx - 2];
	double* beta1 = new double[g_Nx - 2];
	double* s1 = new double[g_Nx - 2];  // ������ ������ ����� ����, ������� ������� ��������� ����� w_np12 ��� ������� i
	double* alpha2 = new double[g_Ny - 2];
	double* beta2 = new double[g_Ny - 2];
	double* s2 = new double[g_Ny - 2];

	// ������� ������ �� ��������, ������� ������� ��������� ������������ ������� alpha �� ����� �� t
	double a1 = -g_tau / (2.0 * g_hx_sq);
	double b1 = a1;
	double c1 = g_tau / (g_hx_sq) + 1.0;
	alpha1[0] = -b1 / c1;
	for (int j = 1; j < g_Nx - 2; j++)
		alpha1[j] = -b1 / (c1 + alpha1[j-1]*a1);
	double a2 = -g_tau / (2 * g_hy_sq);
	double b2 = a2;
	double c2 = g_tau / (g_hy_sq) + 1.0;
	alpha2[0] = -b2 / c2;
	for (int i = 1; i < g_Ny - 2; i++)
		alpha2[i] = -b2 / (c2 + alpha2[i - 1] * a2);

	
	double tau2 = 1.0;  // ������ ���������� ������� � ����
	double t2 = 0.0;
	double t = 0.0;
	double x, y;
	double maxRes;
	printf("%10s%14s%14s\n", "t", "max(|Res|)", "Progress, %");
	while (t < g_T) {
		t += g_tau;
		//printf("t = %.4f of %.4f\n", t, g_T);
		//printf("%10.4f%14.2f\n", t, t / g_T * 100.0);

		// ���������� ������������� ������� �� ��������� ����
		assignW_np12OnBorders(w_np12, t);
		for (int i = 1; i <= g_Ny - 2; i++) {

			// ���������� ������� s1
			for (int j = 1; j <= g_Nx - 2; j++) {
				x = g_xMin + j * g_hx;
				y = g_yMin + i * g_hy;
				s1[j - 1] = g_tau / (2.0 * g_hy_sq) * (w_n[i-1][j] - 2*w_n[i][j] + w_n[i+1][j]);
				s1[j - 1] += g_tau / 2.0 * f(t + g_tau / 2.0, x, y);
				s1[j - 1] += w_n[i][j];
			}
			s1[0] -= a1 * w_np12[i][0];                   // ���� ��������� ������� ��� x = xMin
			s1[g_Nx - 3] -= a1 * w_np12[i][g_Nx - 1];     // ���� ��������� ������� ��� x = xMax

			// ���������� ������� beta1
			beta1[0] = s1[0] / c1;
			for (int i = 1; i < g_Nx - 2; i++)
				beta1[i] = (s1[i] - beta1[i - 1] * a1) / (c1 + alpha1[i - 1] * a1);

			// ��������� w_np12 ��� j �� 1 �� g_Nx - 2
			w_np12[i][g_Nx - 2] = beta1[g_Nx - 3];
			for (int j = g_Nx - 3; j >= 1; j--) {
				w_np12[i][j]= alpha1[j - 1] * w_np12[i][j+1] + beta1[j - 1];
			}
		}

		// ���������� ������������� ������� �� ����� ����
		assignWOnBorders(w_n, t);
		for (int j = 1; j <= g_Nx - 2; j++) {

			// ���������� ������� s2
			for (int i = 1; i <= g_Ny - 2; i++) {
				x = g_xMin + j * g_hx;
				y = g_yMin + i * g_hy;
				s2[i - 1] = g_tau / (2.0 * g_hx_sq) * ( w_np12[i][j-1] - 2*w_np12[i][j] + w_np12[i][j+1] );
				s2[i - 1] += g_tau / 2.0 * f(t + g_tau / 2.0, x, y);
				s2[i - 1] += w_np12[i][j];
			}
			s2[0] -= a2 * w_n[0][j];                                 // ���� ��������� ������� ��� y = yMin
			s2[g_Ny - 3] -= a2 * w_n[g_Ny - 1][j];                   // ���� ��������� ������� ��� y = yMax

			// ���������� ������� beta2
			beta2[0] = s2[0] / c2;
			for (int i = 1; i < g_Ny - 2; i++)
				beta2[i] = (s2[i] - beta2[i - 1] * a2) / (c2 + alpha2[i-1]*a2);

			// ��������� w_n �� ����� ���� �� �������
			w_n[g_Ny - 2][j] = beta2[g_Ny - 3];
			for (int i = g_Ny - 3; i >= 1; i--)
				w_n[i][j] = alpha2[i - 1] * w_n[i + 1][j] + beta2[i - 1];
		}

		if (t > t2) {
			maxRes = getMax(res, g_Ny, g_Nx);
			printf("%10.4f%14.4f%14.2f\n", t, maxRes, t / g_T * 100.0);
			//printf("%14.4f\n", maxRes);

			t2 += tau2;
			getExactSolution(u, t);
			getRes(res, w_n, u, g_Ny, g_Nx);
			writeStepResultsToStream(u, w_n, res, maxRes, g_Ny, g_Nx, t, outstream);
		}
	}

	// �������� ���������
	for (int i = 0; i < g_Ny; i++) {
		delete[] w_n[i];
		delete[] w_np12[i];
		delete[] u[i];
		delete[] res[i];
	}
	delete[] w_n;
	delete[] w_np12;
	delete[] u;
	delete[] res;
	
	delete[] alpha1;
	delete[] alpha2;
	delete[] beta1;
	delete[] beta2;
	delete[] s1;
	delete[] s2;

	outstream.close();
	return 0;
}

void assignW_np12OnBorders(double** w_np12, double t) {
	double y, yPast, yNew;
	double tNew = t + g_tau;
	double C = g_tau / (4.0 * g_hy_sq);
	for (int i = 1; i <= g_Ny - 2; i++) {
		y = g_yMin + i * g_hy;
		yNew = y + g_hy;
		yPast = y - g_hy;
		
		w_np12[i][0] = (mu_xMin(t, g_xMin, y) + mu_xMin(tNew, g_xMin, y)) / 2.0;
		w_np12[i][0] -= C * (mu_xMin(tNew, g_xMin, yPast) - 2*mu_xMin(tNew, g_xMin, y) + mu_xMin(tNew, g_xMin, yNew));
		w_np12[i][0] += C * ( mu_xMin(t, g_xMin, yPast) - 2*mu_xMin(t, g_xMin, y) + mu_xMin(t, g_xMin, yNew) );
		
		w_np12[i][g_Nx - 1] = (mu_xMax(t, g_xMax, y) + mu_xMax(tNew, g_xMax, y)) / 2.0;
		w_np12[i][g_Nx - 1] -= C * (mu_xMax(tNew, g_xMax, yPast) - 2*mu_xMax(tNew, g_xMax, y) + mu_xMax(tNew, g_xMax, yNew));
		w_np12[i][g_Nx - 1] += C * ( mu_xMax(t, g_xMax, yPast) - 2*mu_xMax(t, g_xMax, y) + mu_xMax(t, g_xMax, yNew) );
	}
}

void assignWOnBorders(double** w, double t) {
	double y;
	for (int i = 0; i < g_Ny; i++) {
		y = g_yMin + i * g_hy;
		w[i][0] = mu_xMin(t, g_xMin, y);                     // ������� ����� ������� x = xMin
		w[i][g_Nx - 1] = mu_xMax(t, g_xMax, y);              // ������� ����� ������� x = xMax
	}
	double x;
	for (int j = 0; j < g_Nx; j++) {
		x = g_xMin + j * g_hx;
		w[0][j] = mu_yMin(t, x, g_yMin);                     // ������� ����� ������� y = yMin
		w[g_Ny - 1][j] = mu_yMax(t, x, g_yMax);              // ������� ����� ������� y = yMax
	}
}

void initW(double** w) {
	assignWOnBorders(w, 0);

	double x, y;
	for (int i = 1; i < g_Ny - 1; i++) {
		y = g_yMin + g_hy * i;
		for (int j = 1; j < g_Nx - 1; j++) {
			x = g_xMin + g_hx * j;
			w[i][j] = phi(x, y);
		}
	}
}

void gridInfoToFile(std::string filename) {
	std::ofstream outstream(filename);

	outstream << "// xMin = " << g_xMin << "\n";
	outstream << "// xMax = " << g_xMax << "\n";
	outstream << "// yMin = " << g_yMin << "\n";
	outstream << "// yMax = " << g_yMax << "\n";
	outstream << "//\n";
	outstream << "// Nx = " << g_Nx << "\n";
	outstream << "// Ny = " << g_Ny << "\n";
	outstream << "// hx = " << g_hx << "\n";
	outstream << "// hy = " << g_hy << "\n";
	outstream << "//\n";
	outstream << "// tau = " << g_tau << "\n";
	outstream << "// T = " << g_T << "\n";
	outstream << "//\n";
	outstream << "// x = xMin + hx*j, y = yMin + hy*i, s = j + i*Nx\n";
	outstream << "//\n";
	outstream << "// For each gridnode program results are presented in the following style: \n";
	outstream << "//                  ...    ...    ...\n";
	outstream << "// max(res[i][:]) | ...   u[i][j] ...\n";
	outstream << "//                | ...   w[i][j] ...\n";
	outstream << "//                | ... res[i][j] ...\n";
	outstream << "//                  ...    ...    ...\n";
	outstream << "// where u - exact solution, w - approximate solution\n";
	outstream << "// res[i][j] = abs( u[i][j] - w[i][j] )\n";
	outstream << "// max(res[i][:]) - biggest res value among all elements in row i \n";

	outstream << "\n";
	outstream.close();
}

void writeStepResultsToStream(double** u, double** w_n, double** res, double maxRes, int N_rows, int N_cols, double t, std::ostream& streamOut)
{
	streamOut << "t = " << std::setprecision(4) << t << std::endl;
	streamOut << "max(|Res|) = " << maxRes << std::endl;
	for (int i = 0; i < N_rows; i++) {
		streamOut << std::setw(11) << std::setprecision(4) << getMaxInRow(res, N_cols, i) << " | ";
		for (int j = 0; j < N_cols; j++)
			streamOut << std::setw(14) << std::setprecision(4) << u[i][j];
		streamOut << std::endl;

		streamOut << std::setw(14) << "| ";
		for (int j = 0; j < N_cols; j++)
			streamOut << std::setw(14) << std::setprecision(4) << w_n[i][j];
		streamOut << std::endl;

		streamOut << std::setw(14) << "| ";
		for (int j = 0; j < N_cols; j++)
			streamOut << std::setw(14) << std::setprecision(4) << res[i][j];
		streamOut << "\n\n";

	}
	streamOut << "\n";
}

void getExactSolution(double** ar, double t) {
	double x, y;
	for (int i = 0; i < g_Ny; i++) {
		y = g_yMin + g_hy * i;
		for (int j = 0; j < g_Nx; j++) {
			x = g_xMin + g_hx * j;
			ar[i][j] = uExact(t, x, y);
		}
	}
}

void getRes(double** res, double** a1, double** a2, int N_rows, int N_cols) {
	for (int i = 0; i < N_rows; i++)
		for (int j = 0; j < N_cols; j++)
			res[i][j] = abs(a1[i][j] - a2[i][j]);
}

double getMax(double** ar, int N_rows, int N_cols) {
	double temp = -INFINITY;
	for (int i = 0; i < N_rows; i++)
		for (int j = 0; j < N_cols; j++)
			if (ar[i][j] > temp) temp = ar[i][j];
	return temp;
}

double getMaxInRow(double** ar, int N_cols, int rowId) {
	double temp = -INFINITY;
	for (int j = 0; j < N_cols; j++)
		if (ar[rowId][j] > temp) temp = ar[rowId][j];
	return temp;
}