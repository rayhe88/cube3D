/**
 * @file   lagrange.h
 * @brief
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _LAGRANGE_H_
#define _LAGRANGE_H_

int getData(double, int, int, double *, double *, double *, double *);
int getLagParameters(double, double *, double *, double *, int);

double lagrange(double *xxmu, double delta, int i, int n);
double lagrange1(double *xxmu, double delta, int i, int n);
double lagrange2(double *xxmu, double delta, int i, int n);

double interpolador1D(double x0, double *xcon, double *f, int n);
double interpolador1Dder1(double x0, double *xcon, double *f, int n);
double interpolador1Dder2(double x0, double *xcon, double *f, int n);

int getDerivatives1D(double x0, double *xcon, double *f, int n, double *val);

double interpolador3D(double x0, double y0, double z0, double *xcon,
                      double *ycon, double *zcon, double *fun, int n);

int gradient3D(double x0, double y0, double z0, double *xcon, double *ycon,
               double *zcon, double *fun, int n, int trans, const double *,
               double *val);

int getDerivatives3D(double x0, double y0, double z0, double *xcon,
                     double *ycon, double *zcon, double *fun, int n, int,
                     const double *, double *val);

int gradient3DLog(double x0, double y0, double z0, double *xcon, double *ycon,
                  double *zcon, double *fun, int n, int trans, const double *,
                  double min, double *val);

int getDerivatives3DLog(double x0, double y0, double z0, double *xcon,
                        double *ycon, double *zcon, double *fun, int n, int,
                        const double *, double min, double *val);

#endif
