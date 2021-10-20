/**
 * @file   mathTools.h
 * @brief
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */
#include "mathTools.h"

void scalarVector(double a, double r[3]){
    r[0] *= a;
    r[1] *= a;
    r[2] *= a;
}

double distance(double r1[3], double r2[3]) {

    double ret = (r2[0] - r1[0]) * (r2[0] - r1[0]);
    ret += (r2[1] - r1[1]) * (r2[1] - r1[1]);
    ret += (r2[2] - r1[2]) * (r2[2] - r1[2]);

    return sqrt(ret);
}

// Return  vector's norm for a 3D-vector
double getNorm(double v1, double v2, double v3) {

    double valor;

    valor = v1 * v1 + v2 * v2 + v3 * v3;

    valor = sqrt(valor);

    return valor;
}
// Return  vector's norm for a 3D-vector
double getNormVec(double vec[3]) {

    double valor;

    valor = vec[0] * vec[0];
    valor += vec[1] * vec[1];
    valor += vec[2] * vec[2];

    valor = sqrt(valor);

    return valor;
}

// Return  dot product between 2 3D-vectors
double dotProduct(double vecA[3], double vecB[3]) {

    double valor;

    valor = vecA[0] * vecB[0];
    valor += vecA[1] * vecB[1];
    valor += vecA[2] * vecB[2];

    return valor;
}

// This function calculates the cross product between 2 3D-vectors
void crossProduct(double vecA[3], double vecB[3], double vecOut[3]) {

    vecOut[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1];
    vecOut[1] = vecA[2] * vecB[0] - vecA[0] * vecB[2];
    vecOut[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
}

// This function calculates the product between vector and matrix 3x3.
void matVecProduct(double v[3], double m[9], double vecOut[3]) {

    vecOut[0] = v[0] * m[0] + v[1] * m[1] + v[2] * m[2];
    vecOut[1] = v[0] * m[3] + v[1] * m[4] + v[2] * m[5];
    vecOut[2] = v[0] * m[6] + v[1] * m[7] + v[2] * m[8];
}

double detMat(double m[9]) {

    double det;

    det = m[0] * (m[4] * m[8] - m[5] * m[7]);
    det += m[1] * (m[5] * m[6] - m[3] * m[8]);
    det += m[2] * (m[3] * m[7] - m[4] * m[7]);

    return det;
}
//
// Return eigen values and eigenvector for a matrix
//       | a  b  c |
// mat = | 0  d  e |
//       | 0  0  g |
//
//
int eigenVV(double mat[9], double val[3], double vec[9]) {

    if (fabs(mat[3]) > 1.E-7)
        return 3;
    if (fabs(mat[6]) > 1.E-7)
        return 6;
    if (fabs(mat[7]) > 1.E-7)
        return 7;

    double a, b, c, d, e, g;

    a = mat[0];
    b = mat[1];
    c = mat[2];
    d = mat[4];
    e = mat[5];
    g = mat[8];

    val[0] = a;
    val[1] = d;
    val[2] = g;

    vec[0] = (double)1.;
    vec[1] = (double)0.;
    vec[2] = (double)0.;

    vec[3] = b / (d - a);
    vec[4] = (double)1.;
    vec[5] = (double)0.;

    vec[6] = (c * d - b * e - c * g) / ((a - g) * (g - d));
    vec[7] = e / (g - d);
    vec[8] = (double)1.;

    if (fabs(b) < 1.E-7 && fabs(c) < 1.E-7 && fabs(e) < 1.E-7) {
        vec[0] = (double)1.;
        vec[1] = (double)0.;
        vec[2] = (double)0.;

        vec[3] = (double)0.;
        vec[4] = (double)1.;
        vec[5] = (double)0.;

        vec[6] = (double)0.;
        vec[7] = (double)0.;
        vec[8] = (double)1.;
    }

    return 0;
}

double determinant3(double matA[3][3]) {

    double det;

    det = matA[0][0] * (matA[1][1] * matA[2][2] - matA[1][2] * matA[2][1]);
    det += matA[0][1] * (matA[1][2] * matA[2][0] - matA[1][0] * matA[2][2]);
    det += matA[0][2] * (matA[1][0] * matA[2][2] - matA[1][1] * matA[2][0]);

    return det;
}

double mayor(double matAA[][N], int *p, int *q) {
    int i, j;

    *p = 1;
    *q = 0;
    double temp = fabs(matAA[*p][*q]);

    for (i = 2; i < N; i++)
        for (j = 0; j < i; j++)
            if (fabs(matAA[i][j]) > temp) {
                temp = fabs(matAA[i][j]);
                *p = i;
                *q = j;
            }

    return temp;
}

int jacobi(double matAA[N][N], double valores[N], double eigenvectors[N][N]) {

    int i, j, p, q, n = 0;
    double cot, sen, tan, cos, max, ip, iq, temp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            eigenvectors[i][j] = 0.0;
            eigenvectors[i][i] = 1.0;
        }
    }
    max = mayor(matAA, &p, &q);
    while (n < NMAX && max > TOL) {
        cot = (matAA[q][q] - matAA[p][p]) / (2. * matAA[p][q]);

        if (cot >= TOL)
            tan = -cot + sqrt(1. + cot * cot);
        else
            tan = -cot - sqrt(1. + cot * cot);

        cos = 1. / sqrt(1. + tan * tan);
        sen = tan * cos;

        for (j = 0; j < N; j++) {
            if (j != p && j != q) {
                temp = matAA[p][j];
                matAA[p][j] = matAA[p][j] * cos - matAA[q][j] * sen;
                matAA[j][p] = matAA[p][j];
                matAA[q][j] = temp * sen + matAA[q][j] * cos;
                matAA[j][q] = matAA[q][j];
            }
        } // For

        matAA[q][q] = matAA[q][q] + tan * matAA[p][q];
        matAA[p][p] = matAA[p][p] - tan * matAA[p][q];
        matAA[p][q] = 0.;
        matAA[q][p] = 0.;

        for (j = 0; j < N; j++) {
            ip = cos * eigenvectors[j][p] + sen * eigenvectors[j][q];
            iq = -sen * eigenvectors[j][p] + cos * eigenvectors[j][q];
            eigenvectors[j][p] = ip;
            eigenvectors[j][q] = iq;
        }
        n++;
        max = mayor(matAA, &p, &q);

    } // While

    if (n == NMAX)
        printf(" Maximum number of iterations\n");
    j = 0;
    for (i = 0; i < N; i++) {
        valores[i] = matAA[i][i];
        if (valores[i] > 0)
            j = i;
    }
    if (j != 0) {
        sen = eigenvectors[0][0]; // x
        cos = eigenvectors[0][1]; // y
        cot = eigenvectors[0][2]; // z
        eigenvectors[0][0] = eigenvectors[j][0];
        eigenvectors[0][1] = eigenvectors[j][1];
        eigenvectors[0][2] = eigenvectors[j][2];
        eigenvectors[j][0] = sen;
        eigenvectors[j][1] = cos;
        eigenvectors[j][2] = cot;

        sen = valores[0];
        valores[0] = valores[j];
        valores[j] = sen;
    }

    return 0;
}

int valoresPropios3x3(double *matA, double *eigenVal) {
    //         | a1  a2  a3 |
    //  matA = | a2  a4  a5 |
    //         | a3  a5  a6 |
    int i, j;
    double a1, a2, a3, a4, a5, a6;
    double ta, tb, tc;
    double q, r, sq, theta, aux;
    a1 = matA[0];
    a2 = matA[1];
    a3 = matA[2];
    a4 = matA[4];
    a5 = matA[5];
    a6 = matA[8];

    ta = -a1 - a4 - a6;
    tb = -a2 * a2 - a3 * a3 + a1 * a4 - a5 * a5 + a1 * a6 + a4 * a6;
    tc = a3 * a3 * a4 - 2. * a2 * a3 * a5 + a1 * a5 * a5 + a2 * a2 * a6 -
         a1 * a4 * a6;

    q = (ta * ta - 3. * tb) / 9.;
    r = (2. * ta * ta * ta - 9. * ta * tb + 27. * tc) / 54;
    theta = acos(r / sqrt(q * q * q));
    sq = sqrt(q);

    eigenVal[0] = -2. * sq * cos(theta / 3.) - ta / 3.;
    eigenVal[1] = -2. * sq * cos((theta + 2. * M_PI) / 3.) - ta / 3.;
    eigenVal[2] = -2. * sq * cos((theta - 2. * M_PI) / 3.) - ta / 3.;

    for (i = 0; i < 2; i++)
        for (j = i + 1; j < 3; j++) {
            if (eigenVal[i] > eigenVal[j]) {
                aux = eigenVal[i];
                eigenVal[i] = eigenVal[j];
                eigenVal[j] = aux;
            }
        }
}

int valoresPropios3x3_v0(double *matA, double *eigenVal) {

    //        | a  b  c |
    //  matA =| b  d  e |
    //        | c  e  g |

    double a, b, c, d, e, g;
    double aux;
    double det;
    double p, eig1, eig2, eig3, q, detB, r, phi;
    int i, j;

    a = matA[0];
    b = matA[1];
    c = matA[2];
    d = matA[3];
    e = matA[4];
    g = matA[5];

    det = a * d * g + 2. * b * c * e - (c * c * d + a * e * e + b * b * g);
    if (fabs(det) < 1.E-18) {
        eigenVal[0] = 0.;
        eigenVal[1] = 0.;
        eigenVal[2] = 0.;
        return 0;
    }

    p = b * b + c * c + e * e;
    if (p == 0) {
        eig1 = a;
        eig2 = d;
        eig3 = g;
    } else {
        q = (a + d + g) / 3.;
        p = (a - q) * (a - q) + (d - q) * (d - q) + (g - q) * (g - q) + 2. * p;
        p = sqrt(p / 6.);
        detB = (a - q) * ((d - q) * (g - q) - e * e);
        detB = detB - b * (b * (g - q) - e * c);
        detB = detB + c * (b * e - c * (d - q));
        detB = detB / (p * p * p);
        r = detB / 2.;
        if (r <= -1.)
            phi = M_PI / 3.;
        else {
            if (r >= 1)
                phi = 0.;
            else
                phi = acos(r) / 3.;
        }
        eig1 = q + 2. * p * cos(phi);
        eig3 = q + 2. * p * cos(phi + M_PI * (2. / 3.));
        eig2 = 3. * q - eig1 - eig3;
    }
    eigenVal[0] = eig1;
    eigenVal[1] = eig2;
    eigenVal[2] = eig3;

    for (i = 0; i < 2; i++)
        for (j = i + 1; j < 3; j++) {
            if (eigenVal[i] > eigenVal[j]) {
                aux = eigenVal[i];
                eigenVal[i] = eigenVal[j];
                eigenVal[j] = aux;
            }
        }

    return 0;
}

int eValNci(double *val) {
    int ret;
    double matH[6];
    double eval[3];

    matH[0] = val[4];
    matH[1] = val[7];
    matH[2] = val[8];
    matH[3] = val[5];
    matH[4] = val[9];
    matH[5] = val[6];
    valoresPropios3x3_v0(matH, eval);
    if (eval[1] >= 0)
        ret = 1;
    else
        ret = -1;
    return ret;
}
