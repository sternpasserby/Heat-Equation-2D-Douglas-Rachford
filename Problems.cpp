#include "Problems.h"

double uExact1(double t, double x, double y) {
	return (exp(-2 * t) + 2. * t - 1) / 4 * sin(x) * sin(y) + exp(-25 * t)*sin(4 * x) * sin(3 * y);
}

double f1(double t, double x, double y) {
	return t * sin(x) * sin(y);
}

double phi1(double x, double y) {
	return sin(4 * x) * sin(3 * y);
}

double mu1(double t, double x, double y) {
	return 0;
}

static double x01 = 3.5;
static double y01 = -1.5;
static double sigmaX1 = 1;
static double sigmaY1 = 3;
double uExact2(double t, double x, double y) {
	double u = 5.0 / (t + 1.0) * exp( -(x - x01)*(x - x01)/(2.0*sigmaX1*sigmaX1) - (y - y01)* (y - y01) / (2 * sigmaY1 *sigmaY1));
	return u;
}

double f2(double t, double x, double y) {
	double a = uExact2(t, x, y);
	double sx_sq = sigmaX1 * sigmaX1;
	double sy_sq = sigmaY1 * sigmaY1;
	//return a * (2.0 - 1.0 / (t + 1.0) - (x - x01) * (x - x01) / (sx_sq * sx_sq) - (y - y01) * (y - y01) / (sy_sq * sy_sq));
	return a * (1.0 / sx_sq + 1.0 / sy_sq - (x - x01) * (x - x01) / (sx_sq * sx_sq) - (y - y01) * (y - y01) / (sy_sq * sy_sq) - 1.0 / (t + 1.0));
	
}

double phi2(double x, double y) {
	return uExact2(0, x, y);
}

double mu_xMin2(double t, double x, double y) {
	return uExact2(t, x, y);
}

double mu_xMax2(double t, double x, double y) {
	return uExact2(t, x, y);
}

double mu_yMin2(double t, double x, double y) {
	return uExact2(t, x, y);
}

double mu_yMax2(double t, double x, double y) {
	return uExact2(t, x, y);
}
