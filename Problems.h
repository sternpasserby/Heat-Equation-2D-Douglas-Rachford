#pragma once
#define _USE_MATH_DEFINES //Для числа пи (константа M_PI)
#include <math.h>

double uExact1(double t, double x, double y);
double f1(double t, double x, double y);
double phi1(double x, double y);
double mu1(double t, double x, double y);

double uExact2(double t, double x, double y);
double f2(double t, double x, double y);
double phi2(double x, double y);
double mu_xMin2(double t, double x, double y);
double mu_xMax2(double t, double x, double y);
double mu_yMin2(double t, double x, double y);
double mu_yMax2(double t, double x, double y);