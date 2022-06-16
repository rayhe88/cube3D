/**
 * @file   fields.c
 * @brief
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */

#include "fields.h"
#include <omp.h>

/**
 * @brief
 * @param
 * @param
 * @return
 */
int myTernary(int i, int izq, int der, int n) {

    int outOfRange;
    n = n - 1;

    if (i + izq < 0 || i + der > n) {
        outOfRange = YES;
    } else {
        outOfRange = NOT;
    }

    return outOfRange;
}
///////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @param
 * @param
 * @return
 */
double getGrd(double *val) {
    double ret;
    ret = val[1] * val[1] + val[2] * val[2] + val[3] * val[3];

    ret = sqrt(ret);

    return ret;
}

/**
 * @brief
 * @param
 * @param
 * @return
 */
double getRed(double *val) {
    double ret;
    double den;
    ret = val[1] * val[1] + val[2] * val[2] + val[3] * val[3];
    ret = sqrt(ret);

    den = fabs(val[0]);

    if (den < 1.E-9)
        return 1.E9;

    den = pow(den, FOT);

    ret *= (CF * den);

    return ret;
}

/**
 * @brief
 * @param
 * @param
 * @return
 */
double getLap(double *val) { return val[4] + val[5] + val[6]; }

/**
 * @brief This function get the kinetic energy from density
 * and the first and second derivative. The expresion can be
 * found in the article:
 * Yu.A. Abramov, Acta Cryst. A53 1997 264–272.
 * On the Possibility of Kinetic Energy Density Evaluation
 * from the Experimental Electron-Density Distribution.
 * @param
 * @param
 * @return
 */
double getKin(double *val) {
    double term1, term2, term3;
    double grad2;

    grad2 = val[1] * val[1] + val[2] * val[2] + val[3] * val[3];

    term1 = 2.871234000188192 * pow(val[0], 1.66666666666);
    term2 = 0.01388888888888889 * grad2 / val[0];
    term3 = 0.1666666666666667 * (val[4] + val[5] + val[6]);

    return term1 + term2 + term3;
}

double getKEW(double *val) {
    double grad,rho,kew;
    grad = val[1] * val[1] + val[2] * val[2] + val[3] * val[3];
    grad = sqrt(grad);
    rho = sqrt(fabs(val[0]));

    kew = (double) 0.;
    // In the limit when rho^1/2 and grad go to zero, the value is zero too
    // this is a principal diference with reduced gradient the value goes to infinity.
    if( rho > 1.E-10 && grad > 1.E-10) 
        kew = grad / rho;

    return kew;
}

double getVir(double *val) {
    double lap = getLap(val);
    double kinG = getKin(val);

    return (0.25 * lap - 2. * kinG);
}
/**
 * @brief
 * @param
 * @param
 * @return
 */
int getFieldNoPer(dataCube cube, dataRun param, const double *matU,
                  double *field2) {
    int i, j, k;
    int nx, ny, nz;
    int n1, n2;
    int idxGlobal;
    int index[3];
    int vecn[3];
    double x, y, z, tmp;

    double (*gfunc)(double *);

    nx = cube.pts[0];
    ny = cube.pts[1];
    nz = cube.pts[2];

    vecn[0] = nx;
    vecn[1] = ny;
    vecn[2] = nz;

    n1 = ny * nz;
    n2 = nz;

    switch (param.task) {
    case RED:
        gfunc = &getRed;
        break;
    case GRA:
        gfunc = &getGrd;
        break;
    case LAP:
        gfunc = &getLap;
        break;
    case KIN:
        gfunc = &getKin;
        break;
    case KEW:
        gfunc = &getKEW;
        break;
    case VIR:
        gfunc = &getVir;
        break;
    }

    int npt = nx * ny * nz;

#pragma omp parallel private(idxGlobal, i, j, k, index, x, y, z, tmp)          \
    shared(npt, n1, n2, cube, vecn, param, matU, gfunc, field2)
    {

#pragma omp single
        printf(" Number of threads for evaluation : %4d\n",
               omp_get_num_threads());

#pragma omp barrier

#pragma omp for schedule(dynamic)
        for (idxGlobal = 0; idxGlobal < npt; idxGlobal++) {
            i = idxGlobal / n1;
            j = (idxGlobal - i * n1) / n2;
            k = idxGlobal % n2;
            index[0] = i;
            index[1] = j;
            index[2] = k;

            x = cube.min[0] + i * cube.hvec[0];
            y = cube.min[1] + j * cube.hvec[1];
            z = cube.min[2] + k * cube.hvec[2];

            tmp = campoNoPer(x, y, z, index, vecn, param, cube.hvec, cube.field,
                             matU, gfunc, cube.min);

            field2[idxGlobal] = tmp;
        }
    }

    return 0;
}

/**
 * @brief
 * @param
 * @param
 * @return
 */
int getFieldPer(dataCube cube, dataRun p, const double *matU, double *field2) {
    int i, j, k;
    int nx, ny, nz;
    int n1, n2;
    int idxGlobal;
    int index[3];
    int vecn[3];

    int (*fpol)(double, double, double, double *, double *);
    double (*gfunc)(double *);

    nx = cube.pts[0];
    ny = cube.pts[1];
    nz = cube.pts[2];

    vecn[0] = nx;
    vecn[1] = ny;
    vecn[2] = nz;

    n1 = ny * nz;
    n2 = nz;

    if (p.task == GRA || p.task == RED || p.task == KEW) {
        switch (p.pol) {
        case 1:
            fpol = &gradPol01;
            break;
        case 2:
            fpol = &gradPol02;
            break;
        case 3:
            fpol = &gradPol03;
            break;
        case 4:
            fpol = &gradPol04;
            break;
        case 5:
            fpol = &gradPol05;
            break;
        case 6:
            fpol = &gradPol06;
            break;
        default:
            fpol = &gradPol02;
        }
    }

    if (p.task == LAP || p.task == KIN || p.task == VIR) {
        switch (p.pol) {
        case 1:
            fpol = &hessPol01;
            break;
        case 2:
            fpol = &hessPol02;
            break;
        case 3:
            fpol = &hessPol03;
            break;
        case 4:
            fpol = &hessPol04;
            break;
        case 5:
            fpol = &hessPol05;
            break;
        case 6:
            fpol = &hessPol06;
            break;
        default:
            fpol = &hessPol02;
        }
    }

    switch (p.task) {
    case RED:
        gfunc = &getRed;
        break;
    case GRA:
        gfunc = &getGrd;
        break;
    case LAP:
        gfunc = &getLap;
        break;
    case KIN:
        gfunc = &getKin;
        break;
    case VIR:
        gfunc = &getVir;
        break;
    case KEW:
        gfunc = &getKEW;
        break;
    }

    int npt = nx * ny * nz;

#pragma omp parallel private(idxGlobal, i, j, k, index)                        \
    shared(field2, vecn, p, cube, matU)
    {
#pragma omp single
        printf(" Number of threads for evaluation : %4d\n",
               omp_get_num_threads());

#pragma omp barrier

#pragma omp for schedule(dynamic)
        for (idxGlobal = 0; idxGlobal < npt; idxGlobal++) {
            i = idxGlobal / n1;
            j = (idxGlobal - i * n1) / n2;
            k = idxGlobal % n2;
            index[0] = i;
            index[1] = j;
            index[2] = k;

            field2[idxGlobal] = campoPer(index, vecn, p, cube.hvec, cube.field,
                                         matU, fpol, gfunc);
        }

    } // fin de parallalel section

    return 0;
}

/**
 * @brief
 * @param
 * @param
 * @return
 */
double campoPer(int *index, int *vecn, dataRun param, double *h, double *field,
                const double *matU,
                int (*fpol)(double, double, double, double *, double *),
                double (*gfun)(double *)) {
    int i, j, k;
    int p, q, r, mu;
    int p2, q2, r2;
    int n1, n2;
    double f[param.size];
    double val[10];

    n1 = vecn[1] * vecn[2];
    n2 = vecn[2];

    i = index[0];
    j = index[1];
    k = index[2];

    mu = 0;
    for (p = i + param.izq; p <= i + param.der; p++) {
        p2 = getPeriodicIndex(p, vecn[0]);

        for (q = j + param.izq; q <= j + param.der; q++) {
            q2 = getPeriodicIndex(q, vecn[1]);

            for (r = k + param.izq; r <= k + param.der; r++) {
                r2 = getPeriodicIndex(r, vecn[2]);

                f[mu] = field[p2 * n1 + q2 * n2 + r2];
                mu++;
            }
        }
    }

    fpol(h[0], h[1], h[2], f, val);

    if (param.orth != YES)
        trans02(val, matU);

    return gfun(val);
}

/**
 * @brief
 * @param
 * @param
 * @return
 */
int getNCIPer(dataCube cube, dataRun p, const double *matU, double *red1,
              double *rho2) {
    int i, j, k;
    int nx, ny, nz;
    int n1, n2;
    int idxGlobal;
    int index[3];
    int vecn[3];
    double ret[2];

    nx = cube.pts[0];
    ny = cube.pts[1];
    nz = cube.pts[2];

    vecn[0] = nx;
    vecn[1] = ny;
    vecn[2] = nz;

    n1 = ny * nz;
    n2 = nz;

    int npt = nx * ny * nz;

#pragma omp parallel private(idxGlobal, i, j, k, index, ret)                   \
    shared(red1, rho2, vecn, p, cube, matU)
    {
#pragma omp single
        printf(" Number of threads for evaluation : %4d\n",
               omp_get_num_threads());

#pragma omp barrier

#pragma omp for schedule(dynamic)
        for (idxGlobal = 0; idxGlobal < npt; idxGlobal++) {
            i = idxGlobal / n1;
            j = (idxGlobal - i * n1) / n2;
            k = idxGlobal % n2;
            index[0] = i;
            index[1] = j;
            index[2] = k;

            NCIPer(index, vecn, p, cube.hvec, cube.field, matU, ret);

            red1[idxGlobal] = ret[0];
            rho2[idxGlobal] = ret[1];
        }
    }
    /*  for( i = 0; i < nx; i++ ){
        for( j = 0; j < ny; j++ ){
          for( k = 0; k < nz; k++ ){
            index[0] = i; index[1] = j; index[2] = k;

            idxGlobal = i*n1 + j*n2 + k;

            NCIPer (index,vecn,p,cube.hvec,cube.field,matU,ret);
            red1[idxGlobal] = ret[0];
            rho2[idxGlobal] = ret[1];


          }
        }
      }
     */

    return 0;
}

/**
 * @brief
 * @param
 * @param
 * @return
 */
int NCIPer(int *index, int *vecn, dataRun param, double *h, double *field,
           const double *matU, double *ret) {
    int i, j, k;
    int p, q, r, mu;
    int p2, q2, r2;
    int n1, n2;
    double f[param.size];
    double l2;
    double val[10];
    double rgd, rho, tmp;

    n1 = vecn[1] * vecn[2];
    n2 = vecn[2];

    i = index[0];
    j = index[1];
    k = index[2];

    rho = field[i * n1 + j * n2 + k];

    tmp = fabs(rho);

    if (tmp >= param.la2) {
        ret[0] = param.rgd;
        ret[1] = 100. * param.la2;
        return 0;
    }

    mu = 0;
    for (p = i + param.izq; p <= i + param.der; p++) {
        p2 = getPeriodicIndex(p, vecn[0]);

        for (q = j + param.izq; q <= j + param.der; q++) {
            q2 = getPeriodicIndex(q, vecn[1]);

            for (r = k + param.izq; r <= k + param.der; r++) {
                r2 = getPeriodicIndex(r, vecn[2]);

                f[mu] = field[p2 * n1 + q2 * n2 + r2];
                mu++;
            }
        }
    }

    switch (param.pol) {
    case 1:
        gradPol01(h[0], h[1], h[2], f, val);
        break;
    case 2:
        gradPol02(h[0], h[1], h[2], f, val);
        break;
    case 3:
        gradPol03(h[0], h[1], h[2], f, val);
        break;
    case 4:
        gradPol04(h[0], h[1], h[2], f, val);
        break;
    case 5:
        gradPol05(h[0], h[1], h[2], f, val);
        break;
    case 6:
        gradPol06(h[0], h[1], h[2], f, val);
        break;
    }

    if (param.orth != YES)
        trans01(val, matU);

    rho = val[0];
    rgd = getRed(val);

    tmp = fabs(rho);

    if (tmp < param.la2)
        if (rgd < param.rgd) {
            switch (param.pol) {
            case 1:
                hessPol01(h[0], h[1], h[2], f, val);
                break;
            case 2:
                hessPol02(h[0], h[1], h[2], f, val);
                break;
            case 3:
                hessPol03(h[0], h[1], h[2], f, val);
                break;
            case 4:
                hessPol04(h[0], h[1], h[2], f, val);
                break;
            case 5:
                hessPol05(h[0], h[1], h[2], f, val);
                break;
            case 6:
                hessPol06(h[0], h[1], h[2], f, val);
                break;
            }

            if (param.orth != YES)
                trans02(val, matU);

            l2 = eValNci(val);

            rgd = getRed(val);
            rho = l2 * 100. * val[0];

        } else {
            rgd = param.rgd;
            rho = 100. * param.la2;
        }

    ret[0] = rgd;
    ret[1] = rho;

    return 0;
}

/**
 * @brief
 * @param
 * @param
 * @return
 */
double campoNoPer(double x, double y, double z, int *index, int *vecn,
                  dataRun param, double *h, double *field, const double *matU,
                  double (*gfun)(double *), double *min) {
    int i, j, k;
    int p, q, r, mu;
    int n1, n2;
    double f[param.size];
    double xx[param.pol + 1];
    double yy[param.pol + 1];
    double zz[param.pol + 1];
    double val[10];

    n1 = vecn[1] * vecn[2];
    n2 = vecn[2];

    checkIndex(index, vecn, param);

    i = index[0];
    j = index[1];
    k = index[2];

    mu = 0;
    for (p = i + param.izq; p <= i + param.der; p++) {
        xx[mu] = min[0] + p * h[0];
        mu++;
    }
    mu = 0;
    for (q = j + param.izq; q <= j + param.der; q++) {
        yy[mu] = min[1] + q * h[1];
        mu++;
    }
    mu = 0;
    for (r = k + param.izq; r <= k + param.der; r++) {
        zz[mu] = min[2] + r * h[2];
        mu++;
    }

    mu = 0;
    for (p = i + param.izq; p <= i + param.der; p++) {

        for (q = j + param.izq; q <= j + param.der; q++) {

            for (r = k + param.izq; r <= k + param.der; r++) {
                f[mu] = field[p * n1 + q * n2 + r];

                mu++;
            }
        }
    }

    if (param.task == RED || param.task == GRA)
        gradient3D(x, y, z, xx, yy, zz, f, param.pol, param.orth, matU, val);
    else
        getDerivatives3D(x, y, z, xx, yy, zz, f, param.pol, param.orth, matU,
                         val);

    return gfun(val);
}

/**
 * @brief
 * @param
 * @param
 * @return
 */
int getNCINoPer(dataCube cube, dataRun param, const double *matU, double *red1,
                double *rho2) {
    int i, j, k;
    int nx, ny, nz;
    int n1, n2;
    int idxGlobal;
    int index[3];
    int vecn[3];
    int orx, ory, orz;
    double x, y, z;
    double ret[2];

    nx = cube.pts[0];
    ny = cube.pts[1];
    nz = cube.pts[2];

    vecn[0] = nx;
    vecn[1] = ny;
    vecn[2] = nz;

    n1 = ny * nz;
    n2 = nz;

    int npt = nx * ny * nz;

#pragma omp parallel private(idxGlobal, i, j, k, index, ret, orx, ory, orz)    \
    shared(red1, rho2, vecn, param, cube, nx, ny, nz, matU)
    {
#pragma omp single
        printf(" Number of threads for evaluation : %4d\n",
               omp_get_num_threads());

#pragma omp barrier

#pragma omp for schedule(dynamic)
        for (idxGlobal = 0; idxGlobal < npt; idxGlobal++) {
            i = idxGlobal / n1;
            j = (idxGlobal - i * n1) / n2;
            k = idxGlobal % n2;
            index[0] = i;
            index[1] = j;
            index[2] = k;

            orx = myTernary(i, param.izq, param.der, nx);
            ory = myTernary(j, param.izq, param.der, ny);
            orz = myTernary(k, param.izq, param.der, nz);

            if (orx == YES || ory == YES || orz == YES) {
                x = cube.min[0] + i * cube.hvec[0];
                y = cube.min[1] + j * cube.hvec[1];
                z = cube.min[2] + k * cube.hvec[2];

                NCINoPer(x, y, z, index, vecn, param, cube.hvec, cube.field,
                         matU, cube.min, ret);

            } else {
                NCIPer(index, vecn, param, cube.hvec, cube.field, matU, ret);
            }

            red1[idxGlobal] = ret[0];
            rho2[idxGlobal] = ret[1];
        }
    }

    /*for( i = 0; i < nx; i++ ){
      orx = myTernary( i, param.izq, param.der,nx);
      for( j = 0; j < ny; j++ ){
        ory = myTernary( j, param.izq, param.der,ny);
        for( k = 0; k < nz; k++ ){
          orz = myTernary( k, param.izq, param.der,nz);

          index[0] = i; index[1] = j; index[2] = k;

          idxGlobal = i*n1 + j*n2 + k;
          if( orx == YES || ory == YES || orz == YES ){
            x = cube.min[0] + i*cube.hvec[0];
            y = cube.min[1] + j*cube.hvec[1];
            z = cube.min[2] + k*cube.hvec[2];

            NCINoPer(x,y,z,index,vecn,param,cube.hvec,cube.field,matU,cube.min,ret);

          }else{
            NCIPer (index,vecn,param,cube.hvec,cube.field,matU,ret);
          }
          red1[idxGlobal] = ret[0];
          rho2[idxGlobal] = ret[1];


        }
      }
    }
    */

    return 0;
}
/**
 * @brief
 * @param
 * @param
 * @return
 */
int NCINoPer(double x, double y, double z, int *index, int *vecn, dataRun param,
             double *h, double *field, const double *matU, double *min,
             double *ret) {

    int i, j, k;
    int p, q, r, mu;
    int n1, n2;
    double f[param.size];
    double xx[param.pol + 1];
    double yy[param.pol + 1];
    double zz[param.pol + 1];
    double val[10];
    double l2, rgd, rho, tmp;

    n1 = vecn[1] * vecn[2];
    n2 = vecn[2];

    checkIndex(index, vecn, param);

    i = index[0];
    j = index[1];
    k = index[2];

    rho = field[i * n1 + j * n2 + k];
    tmp = fabs(rho);

    if (tmp >= param.la2) {
        ret[0] = param.rgd;
        ret[1] = 100. * param.la2;
        return 0;
    }

    mu = 0;
    for (p = i + param.izq; p <= i + param.der; p++) {
        xx[mu] = min[0] + p * h[0];
        mu++;
    }
    mu = 0;
    for (q = j + param.izq; q <= j + param.der; q++) {
        yy[mu] = min[1] + q * h[1];
        mu++;
    }
    mu = 0;
    for (r = k + param.izq; r <= k + param.der; r++) {
        zz[mu] = min[2] + r * h[2];
        mu++;
    }

    mu = 0;
    for (p = i + param.izq; p <= i + param.der; p++) {

        for (q = j + param.izq; q <= j + param.der; q++) {

            for (r = k + param.izq; r <= k + param.der; r++) {
                f[mu] = field[p * n1 + q * n2 + r];
                mu++;
            }
        }
    }

    gradient3D(x, y, z, xx, yy, zz, f, param.pol, param.orth, matU, val);

    rho = val[0];
    rgd = getRed(val);

    tmp = fabs(rho);

    if (tmp < param.la2)
        if (rgd < param.rgd) {
            getDerivatives3D(x, y, z, xx, yy, zz, f, param.pol, param.orth,
                             matU, val);

            l2 = eValNci(val);

            rgd = getRed(val);
            rho = 100. * l2 * val[0];
        } else {
            rgd = param.rgd;
            rho = 100. * param.la2;
        }

    ret[0] = rgd;
    ret[1] = rho;

    return 0;
}

/***********************************************************/

int gradientVec(int *index, int *vecn, dataRun param, double *h, double *field,
                const double *matU,
                int (*fpol)(double, double, double, double *, double *),
                double grad[]) {

    int i, j, k;
    int p, q, r, mu;
    int p2, q2, r2;
    int n1, n2;
    double f[param.size];
    double val[10];
    double gnorm;

    n1 = vecn[1] * vecn[2];
    n2 = vecn[2];

    i = index[0];
    j = index[1];
    k = index[2];

    mu = 0;
    for (p = i + param.izq; p <= i + param.der; p++) {
        p2 = getPeriodicIndex(p, vecn[0]);

        for (q = j + param.izq; q <= j + param.der; q++) {
            q2 = getPeriodicIndex(q, vecn[1]);

            for (r = k + param.izq; r <= k + param.der; r++) {
                r2 = getPeriodicIndex(r, vecn[2]);

                f[mu] = field[p2 * n1 + q2 * n2 + r2];
                mu++;
            }
        }
    }
    fpol(h[0], h[1], h[2], f, val);

    if (param.orth != YES)
        trans01(val, matU);

    grad[0] = val[1];
    grad[1] = val[2];
    grad[2] = val[3];

    gnorm = sqrt(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);

    grad[0] /= gnorm;
    grad[1] /= gnorm;
    grad[2] /= gnorm;

    printf(" VecGrad % 10.6lf % 10.6lf % 10.6lf ", grad[0], grad[1], grad[2]);
    printf(" Norm   : %10.6E\n",
           sqrt(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]));

    grad[0] *= 2 * h[0];
    grad[1] *= 2 * h[0];
    grad[2] *= 2 * h[0];
    printf(" VecGrad % 10.6lf % 10.6lf % 10.6lf", grad[0], grad[1], grad[2]);
    printf(" Norm   : %10.6E\n",
           sqrt(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]));

    return 0;
}
