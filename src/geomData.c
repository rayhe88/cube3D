#include "geomData.h"
#include "cubeIndex.h"
#include "fields.h"
#include "file.h"
#include "lagrange2.h"
#include "numBondPath.h"
#include "tableP.h"
#include "transU.h"
#include "utils.h"

int getFieldInLine(double min0, dataCube cube, dataRun param, dataRC config,
                   const double *matU, char *name) {

    printBanner("Field-line", stdout);

    int at1, at2;
    int i, n;
    int idx[3];
    double dist;
    double h;
    double r1[3], r2[3];
    double q1[3], q2[3];
    double ux, uy, uz, norm;

    double f[param.size];
    double xx[param.pol + 1];
    double yy[param.pol + 1];
    double zz[param.pol + 1];
    double matT[9];
    double r[3], q[3], rin[3];
    double val[10];
    double f0;

    FILE *out;
    char nameOut[128];

    getMatInv(matU, matT); // The inverse transformation matrix is calculated
    /**********************************************************************************/
    for (i = 0; i < cube.natm; i++) {
        rin[0] = cube.coor[3 * i];
        rin[1] = cube.coor[3 * i + 1];
        rin[2] = cube.coor[3 * i + 2];

        trans00(rin, matU);

        cube.coor[3 * i] = rin[0];
        cube.coor[3 * i + 1] = rin[1];
        cube.coor[3 * i + 2] = rin[2];
    }
    rin[0] = cube.min[0];
    rin[1] = cube.min[1];
    rin[2] = cube.min[2];
    trans00(rin, matU);
    cube.min[0] = rin[0];
    cube.min[1] = rin[1];
    cube.min[2] = rin[2];

    rin[0] = cube.max[0];
    rin[1] = cube.max[1];
    rin[2] = cube.max[2];
    trans00(rin, matU);
    cube.max[0] = rin[0];
    cube.max[1] = rin[1];
    cube.max[2] = rin[2];
    /**********************************************************************************/

    at1 = param.geom[0];
    at2 = param.geom[1];
    if (at1 > cube.natm || at2 > cube.natm)
        return 1;

    sprintf(nameOut, "%s_line%d-%d.dat", name, at1, at2);
    openFile(&out, nameOut, "w+");

    at1--;
    at2--;

    // These coordinates are in the orthogonal system  q
    q1[0] = cube.coor[3 * at1];
    q1[1] = cube.coor[3 * at1 + 1];
    q1[2] = cube.coor[3 * at1 + 2];

    q2[0] = cube.coor[3 * at2];
    q2[1] = cube.coor[3 * at2 + 1];
    q2[2] = cube.coor[3 * at2 + 2];

    // Here  the coordinates qi are transformed to ri
    getRiU(q1, matU, r1);
    getRiU(q2, matU, r2);
    //

    fprintf(out, "# Coordinates atom %4d:  % 10.6lf % 10.6lf % 10.6lf\n",
            at1 + 1, r1[0], r1[1], r1[2]);
    fprintf(out, "# Coordinates atom %4d:  % 10.6lf % 10.6lf % 10.6lf\n",
            at2 + 1, r2[0], r2[1], r2[2]);
    dist = distance(r1, r2);

    ux = r2[0] - r1[0];
    uy = r2[1] - r1[1];
    uz = r2[2] - r1[2];

    norm = getNorm(ux, uy, uz);

    ux /= norm;
    uy /= norm;
    uz /= norm;

    n = ceil(dist * (config.geom_npua0));
    h = dist / (double)(n - 1);

    double (*gfun)(double *);
    int (*gNum)(double, double, double, double *, double *, double *, double *,
                int, int, const double *, double, double *);

    switch (param.task) {
    case RED:
        gfun = &getRed;
        gNum = &gradient3DLog;
        break;
    case GRA:
        gfun = &getGrd;
        gNum = &gradient3DLog;
        break;
    case KEW:
        gfun = &getKEW;
        gNum = &gradient3DLog;
        break;
    case LAP:
        gfun = &getLap;
        gNum = &getDerivatives3DLog;
        break;
    case KIN:
        gfun = &getKin;
        gNum = &getDerivatives3DLog;
        break;
    case VIR:
        gfun = &getVir;
        gNum = &getDerivatives3DLog;
        break;
    default:
        gfun = &getGrd;
        gNum = &gradient3DLog;
        break;
    }

    printf("  Line between the  atoms %3d and %3d\n", at1 + 1, at2 + 1);
    printf("   Number of points      :  % 10d\n", n);
    printf("   Distance in Angstrom  :  % 10.6lf\n", dist * B2A);
    printf("   Step in Angstrom      :  % 10.6lf\n", h * B2A);

    fprintf(out, "# Distance between atoms %3d and %3d is : %10.6lf Angstrom\n",
            at1 + 1, at2 + 2, dist * B2A);
    fprintf(out,
            "# Unitary vector between atoms : % 10.6lf i + % 10.6lf j + % "
            "10.6lf k\n",
            ux, uy, uz);
    fprintf(out, "# Number of points  % 4d  step : % 10.6lf\n", n, h);
    ux *= h;
    uy *= h;
    uz *= h;
    fprintf(
        out,
        "#    h    vector between atoms : % 10.6lf i + % 10.6lf j + % 10.6lf "
        "k\n#\n",
        ux, uy, uz);
    fprintf(out, "#        r         field        field0\n");

    for (i = 0; i < n; i++) {
        r[0] = r1[0] + i * ux;
        r[1] = r1[1] + i * uy;
        r[2] = r1[2] + i * uz;

        getRiU(r, matT, q);
        getIndex3D(cube.pts, q, cube.min, cube.hvec, idx);
        loadLocalField(idx, cube, param, xx, yy, zz, f);
        gNum(q[0], q[1], q[2], xx, yy, zz, f, param.pol, param.orth, matU, min0,
             val);
        f0 = gfun(val);

        fprintf(out, " % 12.6lf % 12.6lf % 12.6lf\n", distance(r1, r) * B2A, f0,
                val[0]);
    }

    fclose(out);

    printBar(stdout);
    printf("  File %s was generated\n", nameOut);

    return 0;
}

int getFieldInPlane(double min0, dataCube cube, dataRun param, dataRC config,
                    const double *matU, char *name) {

    printBanner("Field-plane", stdout);

    int n;
    int i, j;
    int idx[3];
    int at1, at2, at3;
    double val[10];
    double f0;
    double ejeL, ejeM;
    double cgL, cgM;
    double dist, h;
    double r[3], r0[3], r1[3], r2[3], r3[3]; // "real" space.
    double q[3], q1[3], q2[3], q3[3];        // "virtual" orthogonal space.
    double vL[3], vM[3], vN[3], rcg[3];
    double matT[9];
    double f[param.size];
    double xx[param.pol + 1];
    double yy[param.pol + 1];
    double zz[param.pol + 1];
    char sym1[6];
    char sym2[6];
    char sym3[6];

    double rin[3];

    FILE *out;
    char nameOut[128];

    getMatInv(matU, matT); // The inverse transformation matrix is calculated

    /**********************************************************************************/
    for (i = 0; i < cube.natm; i++) {
        rin[0] = cube.coor[3 * i];
        rin[1] = cube.coor[3 * i + 1];
        rin[2] = cube.coor[3 * i + 2];

        trans00(rin, matU);

        cube.coor[3 * i] = rin[0];
        cube.coor[3 * i + 1] = rin[1];
        cube.coor[3 * i + 2] = rin[2];
    }
    rin[0] = cube.min[0];
    rin[1] = cube.min[1];
    rin[2] = cube.min[2];
    trans00(rin, matU);
    cube.min[0] = rin[0];
    cube.min[1] = rin[1];
    cube.min[2] = rin[2];

    rin[0] = cube.max[0];
    rin[1] = cube.max[1];
    rin[2] = cube.max[2];
    trans00(rin, matU);
    cube.max[0] = rin[0];
    cube.max[1] = rin[1];
    cube.max[2] = rin[2];
    /**********************************************************************************/

    at1 = param.geom[0];
    at2 = param.geom[1];
    at3 = param.geom[2];

    if (at1 > cube.natm || at2 > cube.natm || at3 > cube.natm)
        return 1;

    sprintf(nameOut, "%s_plane%d-%d-%d.dat", name, at1, at2, at3);
    openFile(&out, nameOut, "w+");

    printf("  Plane conformed by the atoms %3d, %3d and %3d\n", at1, at2, at3);
    reOrderAtoms(&at1, &at2, &at3, cube, matU);

    at1--;
    at2--;
    at3--;

    // These coordinates are in the orthogonal system  q
    q1[0] = cube.coor[3 * at1];
    q1[1] = cube.coor[3 * at1 + 1];
    q1[2] = cube.coor[3 * at1 + 2];

    q2[0] = cube.coor[3 * at2];
    q2[1] = cube.coor[3 * at2 + 1];
    q2[2] = cube.coor[3 * at2 + 2];

    q3[0] = cube.coor[3 * at3];
    q3[1] = cube.coor[3 * at3 + 1];
    q3[2] = cube.coor[3 * at3 + 2];

    // Here  the coordinates qi are transformed to ri
    getRiU(q1, matU, r1);
    getRiU(q2, matU, r2);
    getRiU(q3, matU, r3);

    printf("   Centre %3d:  % 10.6lf % 10.6lf % 10.6lf [A]\n", at1 + 1,
           r1[0] * B2A, r1[1] * B2A, r1[2] * B2A);
    printf("   Centre %3d:  % 10.6lf % 10.6lf % 10.6lf [A]\n", at2 + 1,
           r2[0] * B2A, r2[1] * B2A, r2[2] * B2A);
    printf("   Centre %3d:  % 10.6lf % 10.6lf % 10.6lf [A]\n", at3 + 1,
           r3[0] * B2A, r3[1] * B2A, r3[2] * B2A);

    dist = getVecLMN(r1, r2, r3, r0, rcg, vL, vM, vN);
    n = ceil(dist * (config.geom_npua1));
    h = dist / (double)(n - 1);

    printf("   Baricentre:  % 10.6lf % 10.6lf % 10.6lf [A]\n", rcg[0] * B2A,
           rcg[1] * B2A, rcg[2] * B2A);

    printf("   Number of points      :  % 10d x %-10d\n", n, n);
    printf("   Distance in Angstrom  :        % 10.6lf\n", n * h * B2A);
    printf("   Step in Angstrom      :        % 10.6lf\n", h * B2A);

    double (*gfun)(double *);
    int (*gNum)(double, double, double, double *, double *, double *, double *,
                int, int, const double *, double, double *);

    switch (param.task) {
    case RED:
        gfun = &getRed;
        gNum = &gradient3DLog;
        break;
    case GRA:
        gfun = &getGrd;
        gNum = &gradient3DLog;
        break;
    case KEW:
        gfun = &getKEW;
        gNum = &gradient3DLog;
        break;
    case LAP:
        gfun = &getLap;
        gNum = &getDerivatives3DLog;
        break;
    case KIN:
        gfun = &getLap;
        gNum = &getDerivatives3DLog;
        break;
    case VIR:
        gfun = &getLap;
        gNum = &getDerivatives3DLog;
        break;
    default:
        gfun = &getGrd;
        gNum = &gradient3DLog;
        break;
    }

    cgL = dotProduct(rcg, vL);
    cgM = dotProduct(rcg, vM);

    if (at1 < cube.natm)
        getAtomicSymbol(cube.zatm[at1], 6, sym1);
    else
        strcpy(sym1, "NNA");

    if (at2 < cube.natm)
        getAtomicSymbol(cube.zatm[at2], 6, sym2);
    else
        strcpy(sym2, "NNA");

    if (at3 < cube.natm)
        getAtomicSymbol(cube.zatm[at3], 6, sym3);
    else
        strcpy(sym3, "NNA");

    fprintf(out, "#  File of information in 2D for the plane %d-%d-%d\n",
            at1 + 1, at2 + 1, at3 + 1);
    fprintf(out,
            "#                         3D Info                   2D Info\n");

    trans3Dto2D(r1, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(out,
            "#  %5d %6s (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
            at1 + 1, sym1, r1[0] * B2A, r1[1] * B2A, r1[2] * B2A, ejeM * B2A,
            ejeL * B2A);

    trans3Dto2D(r2, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(out,
            "#  %5d %6s (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
            at2 + 1, sym2, r2[0] * B2A, r2[1] * B2A, r2[2] * B2A, ejeM * B2A,
            ejeL * B2A);

    trans3Dto2D(r3, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(out,
            "#  %5d %6s (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
            at3 + 1, sym3, r3[0] * B2A, r3[1] * B2A, r3[2] * B2A, ejeM * B2A,
            ejeL * B2A);

    trans3Dto2D(rcg, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(
        out,
        "#  Baricentre   (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
        rcg[0] * B2A, rcg[1] * B2A, rcg[2] * B2A, ejeM * B2A, ejeL * B2A);

    trans3Dto2D(r0, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(
        out,
        "#  Origin Plane (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
        r0[0] * B2A, r0[1] * B2A, r0[2] * B2A, ejeM * B2A, ejeL * B2A);

    fprintf(out, "#     L axis       M axis      field eval        field 0\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            r[0] = r0[0] + i * h * vL[0] + j * h * vM[0];
            r[1] = r0[1] + i * h * vL[1] + j * h * vM[1];
            r[2] = r0[2] + i * h * vL[2] + j * h * vM[2];

            trans3Dto2D(r, vL, vM, cgL, cgM, &ejeL, &ejeM);

            getRiU(r, matT, q);
            getIndex3D(cube.pts, q, cube.min, cube.hvec, idx);
            loadLocalField(idx, cube, param, xx, yy, zz, f);
            gNum(q[0], q[1], q[2], xx, yy, zz, f, param.pol, param.orth, matU,
                 min0, val);
            f0 = gfun(val);

            fprintf(out, " % 12.6lf % 12.6lf  % 14.8lf % 14.8lf\n", ejeM * B2A,
                    ejeL * B2A, f0, val[0]);
        }
        fprintf(out, "\n");
    }

    fclose(out);

    printBar(stdout);
    printf("  File %s was generated\n", nameOut);

    return 0;
}

double getVecLMN(double r1[3], double r2[3], double r3[3], double r0[3],
                 double cg[3], double vecL[3], double vecM[3], double vecN[3]) {

    double tmp1[3], tmp2[3];
    double vecACm[3];
    double dist[3], dst, norm;

    cg[0] = (r1[0] + r2[0] + r3[0]) / (double)3.;
    cg[1] = (r1[1] + r2[1] + r3[1]) / (double)3.;
    cg[2] = (r1[2] + r2[2] + r3[2]) / (double)3.;

    // getCircumcentre(r1, r2, r3, cg);
    getIncentre(r1, r2, r2, cg);
    dist[0] = distance(r1, cg);
    dist[1] = distance(r2, cg);
    dist[2] = distance(r3, cg);
    if (dist[0] >= dist[1] && dist[0] >= dist[2])
        dst = dist[0];
    else if (dist[1] > dist[2])
        dst = dist[1];
    else
        dst = dist[2];

    dst *= 1.75;

    tmp1[0] = (r1[0] - r2[0]);
    tmp1[1] = (r1[1] - r2[1]);
    tmp1[2] = (r1[2] - r2[2]);

    tmp2[0] = (r3[0] - r2[0]);
    tmp2[1] = (r3[1] - r2[1]);
    tmp2[2] = (r3[2] - r2[2]);

    vecACm[0] = r2[0] - 0.5 * (r1[0] + r3[0]);
    vecACm[1] = r2[1] - 0.5 * (r1[1] + r3[1]);
    vecACm[2] = r2[2] - 0.5 * (r1[2] + r3[2]);

    crossProduct(tmp1, tmp2, vecN);
    norm = getNormVec(vecN);
    vecN[0] /= norm;
    vecN[1] /= norm;
    vecN[2] /= norm;

    crossProduct(vecN, vecACm, vecM);
    norm = getNormVec(vecM);
    vecM[0] /= (-1. * norm);
    vecM[1] /= (-1. * norm);
    vecM[2] /= (-1. * norm);

    /*crossProduct(vecN,vecM,vecL);
    norm = getNormVec(vecL);
    vecL[0] /= norm;
    vecL[1] /= norm;
    vecL[2] /= norm;*/
    norm = getNormVec(vecACm);
    vecL[0] = vecACm[0] / norm;
    vecL[1] = vecACm[1] / norm;
    vecL[2] = vecACm[2] / norm;

    r0[0] = cg[0] - dst * vecL[0] - dst * vecM[0];
    r0[1] = cg[1] - dst * vecL[1] - dst * vecM[1];
    r0[2] = cg[2] - dst * vecL[2] - dst * vecM[2];

    double p1[3];

    p1[0] = cg[0] - dst * vecL[0] + dst * vecM[0];
    p1[1] = cg[1] - dst * vecL[1] + dst * vecM[1];
    p1[2] = cg[2] - dst * vecL[2] + dst * vecM[2];

    double m[3], l[3], minmax[2];

    l[0] = dotProduct(r1, vecL);
    l[1] = dotProduct(r2, vecL);
    l[2] = dotProduct(r3, vecL);

    m[0] = dotProduct(r1, vecM);
    m[1] = dotProduct(r2, vecM);
    m[2] = dotProduct(r3, vecM);

    sortCoor(l, m, minmax);

    return distance(r0, p1);
}

void sortCoor(double l[3], double m[3], double minmax[2]) {

    int i, j, tmp;
    double coorm[3], coorl[3];

    double extra = 1.5;

    for (i = 0; i < 3; i++) {
        coorm[i] = m[i];
        coorl[i] = l[i];
    }

    for (i = 0; i < 3; i++)
        for (j = i + 1; j < 3; j++) {
            if (coorm[i] > coorm[j]) {
                tmp = coorm[i];
                coorm[i] = coorm[j];
                coorm[j] = tmp;
            }
        }

    for (i = 0; i < 3; i++)
        for (j = i + 1; j < 3; j++) {
            if (coorl[i] > coorl[j]) {
                tmp = coorl[i];
                coorl[i] = coorl[j];
                coorl[j] = tmp;
            }
        }

    coorm[0] -= extra;
    coorl[0] -= extra;

    coorm[2] += extra;
    coorl[2] += extra;

    if (coorm[0] < coorl[0])
        minmax[0] = coorm[0];
    else
        minmax[0] = coorl[0];

    if (coorm[2] > coorl[2])
        minmax[1] = coorm[2];
    else
        minmax[1] = coorl[2];
}

void origin(double vecL[3], double vecM[3], double vecN[3], double minmax[2],
            double r0[3]) {
    double m0, n0, l0;
    double rx, ry, rz;
    double mx, my, mz;
    double nx, ny, nz;
    double lx, ly, lz;
    double det;

    l0 = minmax[0];
    m0 = minmax[0];
    n0 = 0.;

    lx = vecL[0];
    ly = vecL[1];
    lz = vecL[2];
    mx = vecM[0];
    my = vecM[1];
    mz = vecM[2];
    nx = vecN[0];
    ny = vecN[1];
    nz = vecN[2];

    det = mx * (ny * lz - nz * ly);
    det += my * (nz * lx - nx * lz);
    det += mz * (nx * ly - ny * lx);

    rx = -lz * my * n0 + ly * mz * n0 + lz * m0 * ny - l0 * mz * ny -
         ly * m0 * nz + l0 * my * nz;
    ry = lz * mx * n0 - lx * mz * n0 - lz * m0 * nx + l0 * mz * nx +
         lx * m0 * nz - l0 * mx * nz;
    rz = -ly * mx * n0 + lx * my * n0 + ly * m0 * nx - l0 * my * nx -
         lx * m0 * ny + l0 * mx * ny;

    r0[0] = rx / det;
    r0[1] = ry / det;
    r0[2] = rz / det;
}

void reOrderAtoms(int *at1, int *at2, int *at3, dataCube cube,
                  const double *matU) {
    int i;
    int atA, atB, atC;

    double coorA[3], coorB[3], coorC[3];
    double vecAB[3], vecAC[3], vecBC[3];
    double nAB, nAC, nBC;

    atA = (*at1);
    atB = (*at2);
    atC = (*at3);

    atA--;
    atB--;
    atC--;

    for (i = 0; i < 3; i++) {
        coorA[i] = cube.coor[3 * atA + i];
        coorB[i] = cube.coor[3 * atB + i];
        coorC[i] = cube.coor[3 * atC + i];
    }

    // itrans00(coorA,matU);
    // itrans00(coorB,matU);
    // itrans00(coorC,matU);

    for (i = 0; i < 3; i++) {
        vecAB[i] = coorA[i] - coorB[i];
        vecAC[i] = coorA[i] - coorC[i];
        vecBC[i] = coorB[i] - coorC[i];
    }

    nAB = getNormVec(vecAB);
    nAC = getNormVec(vecAC);
    nBC = getNormVec(vecBC);

    if (nAB > nAC && nAB > nBC) {
        (*at1) = atA + 1;
        (*at2) = atC + 1;
        (*at3) = atB + 1;
    } else {
        if (nAC > nAB && nAC > nBC) {
            (*at1) = atA + 1;
            (*at2) = atB + 1;
            (*at3) = atC + 1;
        } else {
            if (nBC > nAB && nBC > nAC) {
                (*at1) = atB + 1;
                (*at2) = atA + 1;
                (*at3) = atC + 1;
            } else {
                (*at1) = atA + 1;
                (*at2) = atB + 1;
                (*at3) = atC + 1;
            }
        }
    }
}
void getCircumcentre(double rA[3], double rB[3], double rC[3], double c0[3]) {
    // https://en.wikipedia.org/wiki/Circumscribed_circle
    double cc[3];
    double tmpa[3], tmpb[3];
    double tmp1[3], tmp2[3];
    double normA, normB;
    double normcrossAB;
    double crossAB[3];

    tmpa[0] = rA[0] - rC[0];
    tmpa[1] = rA[1] - rC[1];
    tmpa[2] = rA[2] - rC[2];
    normA = getNormVec(tmpa);
    normA = normA * normA;

    tmpb[0] = rB[0] - rC[0];
    tmpb[1] = rB[1] - rC[1];
    tmpb[2] = rB[2] - rC[2];
    normB = getNormVec(tmpb);
    normB = normB * normB;

    crossProduct(tmpa, tmpb, crossAB);
    normcrossAB = getNormVec(crossAB);
    normcrossAB = normcrossAB * normcrossAB;

    tmp1[0] = normA * tmpb[0] - normB * tmpa[0];
    tmp1[1] = normA * tmpb[1] - normB * tmpa[1];
    tmp1[2] = normA * tmpb[2] - normB * tmpa[2];

    crossProduct(tmp1, crossAB, tmp2);

    cc[0] = tmp2[0] / (2. * normcrossAB) + rC[0];
    cc[1] = tmp2[1] / (2. * normcrossAB) + rC[1];
    cc[2] = tmp2[2] / (2. * normcrossAB) + rC[2];

    if (distance(cc, c0) < 2) {
        c0[0] = cc[0];
        c0[1] = cc[1];
        c0[2] = cc[2];
    }
}
void getIncentre(double rA[3], double rB[3], double rC[3], double c0[3]) {
    double a, b, c, den;

    a = distance(rB, rC);
    b = distance(rA, rC);
    c = distance(rA, rB);

    den = a + b + c;

    c0[0] = (a * rA[0] + b * rB[0] + c * rC[0]) / den;
    c0[1] = (a * rA[1] + b * rB[1] + c * rC[1]) / den;
    c0[2] = (a * rA[2] + b * rB[2] + c * rC[2]) / den;
}

// New feature1
int getGradLinesInPlane1(double min0, dataCube cube, dataRun param,
                         dataRC config, const double *matU, char *name) {

    printBanner("GradLines-plane", stdout);

    int n;
    int i, j;
    int idx[3];
    int at1, at2, at3;
    double val[10];
    double ejeL, ejeM;
    double cgL, cgM;
    double dist, h;
    double r0[3], r1[3], r2[3], r3[3]; // "real" space.
    double q1[3], q2[3], q3[3];        // "virtual" orthogonal space.
    double vL[3], vM[3], vN[3], rcg[3];
    double ri[4][3];
    double matT[9];
    double f[param.size];
    double xx[param.pol + 1];
    double yy[param.pol + 1];
    double zz[param.pol + 1];
    char sym1[6];
    char sym2[6];
    char sym3[6];

    double rin[3];

    FILE *out;
    char nameOut[128];

    getMatInv(matU, matT); // The inverse transformation matrix is calculated

    /**********************************************************************************/
    for (i = 0; i < cube.natm; i++) {
        rin[0] = cube.coor[3 * i];
        rin[1] = cube.coor[3 * i + 1];
        rin[2] = cube.coor[3 * i + 2];

        trans00(rin, matU);

        cube.coor[3 * i] = rin[0];
        cube.coor[3 * i + 1] = rin[1];
        cube.coor[3 * i + 2] = rin[2];
    }
    rin[0] = cube.min[0];
    rin[1] = cube.min[1];
    rin[2] = cube.min[2];
    trans00(rin, matU);
    cube.min[0] = rin[0];
    cube.min[1] = rin[1];
    cube.min[2] = rin[2];

    rin[0] = cube.max[0];
    rin[1] = cube.max[1];
    rin[2] = cube.max[2];
    trans00(rin, matU);
    cube.max[0] = rin[0];
    cube.max[1] = rin[1];
    cube.max[2] = rin[2];
    /**********************************************************************************/

    at1 = param.geom[0];
    at2 = param.geom[1];
    at3 = param.geom[2];

    if (at1 > cube.natm || at2 > cube.natm || at3 > cube.natm)
        return 1;

    sprintf(nameOut, "%s_plane%d-%d-%d-vec.dat", name, at1, at2, at3);
    openFile(&out, nameOut, "w+");

    printf("  Plane conformed by the atoms %3d, %3d and %3d\n", at1, at2, at3);
    reOrderAtoms(&at1, &at2, &at3, cube, matU);

    at1--;
    at2--;
    at3--;

    // These coordinates are in the orthogonal system  q
    q1[0] = cube.coor[3 * at1];
    q1[1] = cube.coor[3 * at1 + 1];
    q1[2] = cube.coor[3 * at1 + 2];

    q2[0] = cube.coor[3 * at2];
    q2[1] = cube.coor[3 * at2 + 1];
    q2[2] = cube.coor[3 * at2 + 2];

    q3[0] = cube.coor[3 * at3];
    q3[1] = cube.coor[3 * at3 + 1];
    q3[2] = cube.coor[3 * at3 + 2];

    // Here  the coordinates qi are transformed to ri
    getRiU(q1, matU, r1);
    getRiU(q2, matU, r2);
    getRiU(q3, matU, r3);

    ri[0][0] = r1[0];
    ri[0][1] = r1[1];
    ri[0][2] = r1[2];
    ri[1][0] = r2[0];
    ri[1][1] = r2[1];
    ri[1][2] = r2[2];
    ri[2][0] = r3[0];
    ri[2][1] = r3[1];
    ri[2][2] = r3[2];

    printf("   Centre %3d:  % 10.6lf % 10.6lf % 10.6lf [A]\n", at1 + 1,
           r1[0] * B2A, r1[1] * B2A, r1[2] * B2A);
    printf("   Centre %3d:  % 10.6lf % 10.6lf % 10.6lf [A]\n", at2 + 1,
           r2[0] * B2A, r2[1] * B2A, r2[2] * B2A);
    printf("   Centre %3d:  % 10.6lf % 10.6lf % 10.6lf [A]\n", at3 + 1,
           r3[0] * B2A, r3[1] * B2A, r3[2] * B2A);

    dist = getVecLMN(r1, r2, r3, r0, rcg, vL, vM, vN);
    n = ceil(dist * (config.geom_npua2));
    h = dist / (double)(n - 1);
    ri[3][0] = rcg[0];
    ri[3][1] = rcg[1];
    ri[3][2] = rcg[2];

    printf("   Baricentre:  % 10.6lf % 10.6lf % 10.6lf [A]\n", rcg[0] * B2A,
           rcg[1] * B2A, rcg[2] * B2A);

    printf("   Number of points      :  % 10d x %-10d\n", n, n);
    printf("   Distance in Angstrom  :        % 10.6lf\n", n * h * B2A);
    printf("   Step in Angstrom      :        % 10.6lf\n", h * B2A);

    int (*gNum)(double, double, double, double *, double *, double *, double *,
                int, int, const double *, double, double *);

    gNum = &gradient3DLog;

    cgL = dotProduct(rcg, vL);
    cgM = dotProduct(rcg, vM);

    if (at1 < cube.natm)
        getAtomicSymbol(cube.zatm[at1], 6, sym1);
    else
        strcpy(sym1, "NNA");

    if (at2 < cube.natm)
        getAtomicSymbol(cube.zatm[at2], 6, sym2);
    else
        strcpy(sym2, "NNA");

    if (at3 < cube.natm)
        getAtomicSymbol(cube.zatm[at3], 6, sym3);
    else
        strcpy(sym3, "NNA");

    fprintf(out, "#  File of information in 2D for the plane %d-%d-%d\n",
            at1 + 1, at2 + 1, at3 + 1);
    fprintf(out,
            "#                         3D Info                   2D Info\n");

    trans3Dto2D(r1, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(out,
            "#  %5d %6s (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
            at1 + 1, sym1, r1[0] * B2A, r1[1] * B2A, r1[2] * B2A, ejeM * B2A,
            ejeL * B2A);

    trans3Dto2D(r2, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(out,
            "#  %5d %6s (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
            at2 + 1, sym2, r2[0] * B2A, r2[1] * B2A, r2[2] * B2A, ejeM * B2A,
            ejeL * B2A);

    trans3Dto2D(r3, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(out,
            "#  %5d %6s (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
            at3 + 1, sym3, r3[0] * B2A, r3[1] * B2A, r3[2] * B2A, ejeM * B2A,
            ejeL * B2A);

    trans3Dto2D(rcg, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(
        out,
        "#  Baricentre   (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
        rcg[0] * B2A, rcg[1] * B2A, rcg[2] * B2A, ejeM * B2A, ejeL * B2A);

    trans3Dto2D(r0, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(
        out,
        "#  Origin Plane (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
        r0[0] * B2A, r0[1] * B2A, r0[2] * B2A, ejeM * B2A, ejeL * B2A);
    double d2min[2];
    double d2max[2];

    d2min[0] = ejeL;
    d2min[1] = ejeM;
    d2max[0] = -ejeL;
    d2max[1] = -ejeM;

    fprintf(out, "#     L axis       M axis      field eval        field 0\n");

    double gnorm;

    double eps = 0.2 * B2A;
    double theta = 2 * M_PI / config.geom_nang;
    double rj[3], qi[3];
    int flag;
    int iter;
    double cM, cL;

    for (i = 0; i < 3; i++) {
        rj[0] = ri[i][0];
        rj[1] = ri[i][1];
        rj[2] = ri[i][2];
        trans3Dto2D(rj, vL, vM, cgL, cgM, &ejeL, &ejeM);
        cM = ejeM * B2A;
        cL = ejeL * B2A;
        for (j = 0; j < config.geom_nang; j++) {
            rj[0] = eps * vL[0] * sin(j * theta) +
                    eps * vM[0] * cos(j * theta) + ri[i][0];
            rj[1] = eps * vL[1] * sin(j * theta) +
                    eps * vM[1] * cos(j * theta) + ri[i][1];
            rj[2] = eps * vL[2] * sin(j * theta) +
                    eps * vM[2] * cos(j * theta) + ri[i][2];

            trans3Dto2D(rj, vL, vM, cgL, cgM, &ejeL, &ejeM);

            fprintf(out, " % 12.6lf % 12.6lf  %2d\n", cM, cL, i);
            fprintf(out, " % 12.6lf % 12.6lf  %2d\n", ejeM * B2A, ejeL * B2A,
                    i);
            flag = 0;
            iter = 0;
            while (flag == 0 && iter < 10000) {

                getRiU(rj, matT, qi);
                getIndex3D(cube.pts, qi, cube.min, cube.hvec, idx);
                loadLocalField(idx, cube, param, xx, yy, zz, f);
                gNum(qi[0], qi[1], qi[2], xx, yy, zz, f, param.pol, param.orth,
                     matU, min0, val);

                gnorm = getNorm(val[1], val[2], val[3]);

                if (gnorm < 1.E-10) {
                    iter = 20000;
                    flag = 1;
                    break;
                } else {
                    if (iter == 0) {
                        rj[0] -= 0.1 * val[1] / gnorm;
                        rj[1] -= 0.1 * val[2] / gnorm;
                        rj[2] -= 0.1 * val[3] / gnorm;
                    } else {
                        rj[0] -= 0.05 * val[1] / gnorm;
                        rj[1] -= 0.05 * val[2] / gnorm;
                        rj[2] -= 0.05 * val[3] / gnorm;
                    }

                    trans3Dto2D(rj, vL, vM, cgL, cgM, &ejeL, &ejeM);
                    fprintf(out, " % 12.6lf % 12.6lf  %2d\n", ejeM * B2A,
                            ejeL * B2A, i);

                    if (ejeL < d2min[0] || ejeL > d2max[0]) {
                        flag = 1;
                    } else {
                        if (ejeM < d2min[1] || ejeM > d2max[1]) {
                            flag = 1;
                        }
                    }
                    iter++;
                }
            }
            fprintf(out, "\n");
        }
        fprintf(out, "\n");
    }

    fclose(out);

    printBar(stdout);
    printf("  File %s was generated\n", nameOut);

    return 0;
}
// New feature2
int getGradVectorsInPlane(double min0, dataCube cube, dataRun param,
                          dataRC config, const double *matU, char *name) {

    printBanner("Gradient Vectorial - Plane", stdout);

    int n;
    int i, j;
    int idx[3];
    int at1, at2, at3;
    double val[10];
    double f0;
    double ejeL, ejeM;
    double cgL, cgM;
    double dist, h;
    double r[3], r0[3], r1[3], r2[3], r3[3]; // "real" space.
    double q[3], q1[3], q2[3], q3[3];        // "virtual" orthogonal space.
    double vL[3], vM[3], vN[3], rcg[3];
    double matT[9];
    double f[param.size];
    double xx[param.pol + 1];
    double yy[param.pol + 1];
    double zz[param.pol + 1];
    char sym1[6];
    char sym2[6];
    char sym3[6];

    double rin[3];

    FILE *out;
    char nameOut[128];

    getMatInv(matU, matT); // The inverse transformation matrix is calculated

    /**********************************************************************************/
    for (i = 0; i < cube.natm; i++) {
        rin[0] = cube.coor[3 * i];
        rin[1] = cube.coor[3 * i + 1];
        rin[2] = cube.coor[3 * i + 2];

        trans00(rin, matU);

        cube.coor[3 * i] = rin[0];
        cube.coor[3 * i + 1] = rin[1];
        cube.coor[3 * i + 2] = rin[2];
    }
    rin[0] = cube.min[0];
    rin[1] = cube.min[1];
    rin[2] = cube.min[2];
    trans00(rin, matU);
    cube.min[0] = rin[0];
    cube.min[1] = rin[1];
    cube.min[2] = rin[2];

    rin[0] = cube.max[0];
    rin[1] = cube.max[1];
    rin[2] = cube.max[2];
    trans00(rin, matU);
    cube.max[0] = rin[0];
    cube.max[1] = rin[1];
    cube.max[2] = rin[2];
    /**********************************************************************************/

    at1 = param.geom[0];
    at2 = param.geom[1];
    at3 = param.geom[2];

    if (at1 > cube.natm || at2 > cube.natm || at3 > cube.natm)
        return 1;

    sprintf(nameOut, "%s_plane%d-%d-%d-vField.dat", name, at1, at2, at3);
    openFile(&out, nameOut, "w+");

    printf("  Plane conformed by the atoms %3d, %3d and %3d\n", at1, at2, at3);
    reOrderAtoms(&at1, &at2, &at3, cube, matU);

    at1--;
    at2--;
    at3--;

    // These coordinates are in the orthogonal system  q
    q1[0] = cube.coor[3 * at1];
    q1[1] = cube.coor[3 * at1 + 1];
    q1[2] = cube.coor[3 * at1 + 2];

    q2[0] = cube.coor[3 * at2];
    q2[1] = cube.coor[3 * at2 + 1];
    q2[2] = cube.coor[3 * at2 + 2];

    q3[0] = cube.coor[3 * at3];
    q3[1] = cube.coor[3 * at3 + 1];
    q3[2] = cube.coor[3 * at3 + 2];

    // Here  the coordinates qi are transformed to ri
    getRiU(q1, matU, r1);
    getRiU(q2, matU, r2);
    getRiU(q3, matU, r3);
    // We get the baricentre and the vectors that made the plane
    dist = getVecLMN(r1, r2, r3, r0, rcg, vL, vM, vN);
    n = ceil(dist * (config.geom_npua2));
    h = dist / (double)(n - 1);

    printf("   Centre %3d:  % 10.6lf % 10.6lf % 10.6lf [A]\n", at1 + 1,
           r1[0] * B2A, r1[1] * B2A, r1[2] * B2A);
    printf("   Centre %3d:  % 10.6lf % 10.6lf % 10.6lf [A]\n", at2 + 1,
           r2[0] * B2A, r2[1] * B2A, r2[2] * B2A);
    printf("   Centre %3d:  % 10.6lf % 10.6lf % 10.6lf [A]\n", at3 + 1,
           r3[0] * B2A, r3[1] * B2A, r3[2] * B2A);
    printf("   Baricentre:  % 10.6lf % 10.6lf % 10.6lf [A]\n", rcg[0] * B2A,
           rcg[1] * B2A, rcg[2] * B2A);
    printf("   Number of points      :  % 10d x %-10d\n", n, n);
    printf("   Distance in Angstrom  :        % 10.6lf\n", n * h * B2A);
    printf("   Step in Angstrom      :        % 10.6lf\n", h * B2A);

    cgL = dotProduct(rcg, vL);
    cgM = dotProduct(rcg, vM);

    if (at1 < cube.natm)
        getAtomicSymbol(cube.zatm[at1], 6, sym1);
    else
        strcpy(sym1, "NNA");

    if (at2 < cube.natm)
        getAtomicSymbol(cube.zatm[at2], 6, sym2);
    else
        strcpy(sym2, "NNA");

    if (at3 < cube.natm)
        getAtomicSymbol(cube.zatm[at3], 6, sym3);
    else
        strcpy(sym3, "NNA");

    fprintf(out, "#  File of information in 2D for the plane %d-%d-%d\n",
            at1 + 1, at2 + 1, at3 + 1);
    fprintf(out,
            "#                         3D Info                   2D Info\n");

    trans3Dto2D(r1, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(out,
            "#  %5d %6s (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
            at1 + 1, sym1, r1[0] * B2A, r1[1] * B2A, r1[2] * B2A, ejeM * B2A,
            ejeL * B2A);

    trans3Dto2D(r2, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(out,
            "#  %5d %6s (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
            at2 + 1, sym2, r2[0] * B2A, r2[1] * B2A, r2[2] * B2A, ejeM * B2A,
            ejeL * B2A);

    trans3Dto2D(r3, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(out,
            "#  %5d %6s (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
            at3 + 1, sym3, r3[0] * B2A, r3[1] * B2A, r3[2] * B2A, ejeM * B2A,
            ejeL * B2A);

    trans3Dto2D(rcg, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(
        out,
        "#  Baricentre   (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
        rcg[0] * B2A, rcg[1] * B2A, rcg[2] * B2A, ejeM * B2A, ejeL * B2A);

    trans3Dto2D(r0, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(
        out,
        "#  Origin Plane (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
        r0[0] * B2A, r0[1] * B2A, r0[2] * B2A, ejeM * B2A, ejeL * B2A);

    fprintf(out, "#     L axis       M axis      field eval        field 0\n");

    double gnorm;
    double gvec[3];
    double gL, gM;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            r[0] = r0[0] + i * h * vL[0] + j * h * vM[0];
            r[1] = r0[1] + i * h * vL[1] + j * h * vM[1];
            r[2] = r0[2] + i * h * vL[2] + j * h * vM[2];
            trans3Dto2D(r, vL, vM, cgL, cgM, &ejeL, &ejeM);

            getRiU(r, matT, q);
            getIndex3D(cube.pts, q, cube.min, cube.hvec, idx);
            loadLocalField(idx, cube, param, xx, yy, zz, f);
            gradient3DLog(q[0], q[1], q[2], xx, yy, zz, f, param.pol,
                          param.orth, matU, min0, val);
            f0 = getGrd(val);

            gnorm = getNorm(val[1], val[2], val[3]);

            gvec[0] = val[1] / gnorm;
            gvec[1] = val[2] / gnorm;
            gvec[2] = val[3] / gnorm;

            gL = dotProduct(gvec, vL);
            gM = dotProduct(gvec, vM);

            gnorm = sqrt(gL * gL + gM * gM);

            gL *= h / gnorm;
            gM *= h / gnorm;

            gL *= B2A;
            gM *= B2A;
            ejeL *= B2A;
            ejeM *= B2A;

            fprintf(out,
                    " % 12.6lf % 12.6lf % 12.6lf % 12.6lf  % 14.8lf % 14.8lf\n",
                    ejeM, ejeL, gM, gL, f0, val[0]);
        }
        fprintf(out, "\n");
    }

    fclose(out);

    printBar(stdout);
    printf("  File %s was generated\n", nameOut);

    return 0;
}

void transformCube1(dataCube *cube, const double *matU) {
    int i;
    double rin[3];
    for (i = 0; i < cube->natm; i++) {
        rin[0] = cube->coor[3 * i];
        rin[1] = cube->coor[3 * i + 1];
        rin[2] = cube->coor[3 * i + 2];

        trans00(rin, matU);

        cube->coor[3 * i] = rin[0];
        cube->coor[3 * i + 1] = rin[1];
        cube->coor[3 * i + 2] = rin[2];
    }

    rin[0] = cube->min[0];
    rin[1] = cube->min[1];
    rin[2] = cube->min[2];
    trans00(rin, matU);
    cube->min[0] = rin[0];
    cube->min[1] = rin[1];
    cube->min[2] = rin[2];

    rin[0] = cube->max[0];
    rin[1] = cube->max[1];
    rin[2] = cube->max[2];
    trans00(rin, matU);
    cube->max[0] = rin[0];
    cube->max[1] = rin[1];
    cube->max[2] = rin[2];
}

void printInfoCube1(dataCube cube) {
    int i;
    printf(" Min   :  % 12.6lf  % 12.6lf  % 12.6lf\n", cube.min[0], cube.min[1],
           cube.min[2]);
    printf(" Max   :  % 12.6lf  % 12.6lf  % 12.6lf\n", cube.max[0], cube.max[1],
           cube.max[2]);
    for (i = 0; i < cube.natm; i++)
        printf(" Atm%2d :  % 12.6lf  % 12.6lf  % 12.6lf\n", i, cube.coor[3 * i],
               cube.coor[3 * i + 1], cube.coor[3 * i + 2]);
}

int getStreamLinesInPlane(double min0, dataCube cube, dataRun param,
                          dataRC config, const double *matU, char *name) {

    printBanner("StreamLines plane", stdout);

    int n;
    int i, j;
    int idx[3];
    int at1, at2, at3;
    double val[10];
    double ejeL, ejeM;
    double cgL, cgM, cgN;
    double dist, h;
    double r[3], r0[3], r1[3], r2[3], r3[3]; // "real" space.
    double q[3], q1[3], q2[3], q3[3];        // "virtual" orthogonal space.
    double vL[3], vM[3], vN[3], rcg[3];
    double matT[9];
    double f[param.size];
    double xx[param.pol + 1];
    double yy[param.pol + 1];
    double zz[param.pol + 1];
    char sym1[6];
    char sym2[6];
    char sym3[6];

    double rin[3];

    FILE *out;
    FILE *vec;
    char nameOut[128];
    char nameVec[128];

    getMatInv(matU, matT); // The inverse transformation matrix is calculated

    /**********************************************************************************/
    for (i = 0; i < cube.natm; i++) {
        rin[0] = cube.coor[3 * i];
        rin[1] = cube.coor[3 * i + 1];
        rin[2] = cube.coor[3 * i + 2];

        trans00(rin, matU);

        cube.coor[3 * i] = rin[0];
        cube.coor[3 * i + 1] = rin[1];
        cube.coor[3 * i + 2] = rin[2];
    }
    rin[0] = cube.min[0];
    rin[1] = cube.min[1];
    rin[2] = cube.min[2];
    trans00(rin, matU);
    cube.min[0] = rin[0];
    cube.min[1] = rin[1];
    cube.min[2] = rin[2];

    rin[0] = cube.max[0];
    rin[1] = cube.max[1];
    rin[2] = cube.max[2];
    trans00(rin, matU);
    cube.max[0] = rin[0];
    cube.max[1] = rin[1];
    cube.max[2] = rin[2];
    /**********************************************************************************/

    at1 = param.geom[0];
    at2 = param.geom[1];
    at3 = param.geom[2];

    if (at1 > cube.natm || at2 > cube.natm || at3 > cube.natm)
        return 1;

    sprintf(nameOut, "%s_plane%d-%d-%d-gLine.dat", name, at1, at2, at3);
    openFile(&out, nameOut, "w+");

    sprintf(nameVec, "%s_plane%d-%d-%d-gLineVec.dat", name, at1, at2, at3);
    openFile(&vec, nameVec, "w+");

    printf("  Plane conformed by the atoms %3d, %3d and %3d\n", at1, at2, at3);
    reOrderAtoms(&at1, &at2, &at3, cube, matU);

    at1--;
    at2--;
    at3--;

    // These coordinates are in the orthogonal system  q
    q1[0] = cube.coor[3 * at1];
    q1[1] = cube.coor[3 * at1 + 1];
    q1[2] = cube.coor[3 * at1 + 2];

    q2[0] = cube.coor[3 * at2];
    q2[1] = cube.coor[3 * at2 + 1];
    q2[2] = cube.coor[3 * at2 + 2];

    q3[0] = cube.coor[3 * at3];
    q3[1] = cube.coor[3 * at3 + 1];
    q3[2] = cube.coor[3 * at3 + 2];

    // Here  the coordinates qi are transformed to ri
    getRiU(q1, matU, r1);
    getRiU(q2, matU, r2);
    getRiU(q3, matU, r3);

    // We get the baricentre and the vectors that conform the plane
    dist = getVecLMN(r1, r2, r3, r0, rcg, vL, vM, vN);
    n = ceil(dist * (config.geom_npua2));
    h = dist / (double)(n - 1);

    printf("   Centre %3d:  % 10.6lf % 10.6lf % 10.6lf [A]\n", at1 + 1,
           r1[0] * B2A, r1[1] * B2A, r1[2] * B2A);
    printf("   Centre %3d:  % 10.6lf % 10.6lf % 10.6lf [A]\n", at2 + 1,
           r2[0] * B2A, r2[1] * B2A, r2[2] * B2A);
    printf("   Centre %3d:  % 10.6lf % 10.6lf % 10.6lf [A]\n", at3 + 1,
           r3[0] * B2A, r3[1] * B2A, r3[2] * B2A);
    printf("   Baricentre:  % 10.6lf % 10.6lf % 10.6lf [A]\n", rcg[0] * B2A,
           rcg[1] * B2A, rcg[2] * B2A);

    printf("   Number of points      :  % 10d x %-10d\n", n, n);
    printf("   Distance in Angstrom  :        % 10.6lf\n", n * h * B2A);
    printf("   Step in Angstrom      :        % 10.6lf\n", h * B2A);

    cgL = dotProduct(rcg, vL);
    cgM = dotProduct(rcg, vM);
    cgN = dotProduct(rcg, vN);

    if (at1 < cube.natm)
        getAtomicSymbol(cube.zatm[at1], 6, sym1);
    else
        strcpy(sym1, "NNA");

    if (at2 < cube.natm)
        getAtomicSymbol(cube.zatm[at2], 6, sym2);
    else
        strcpy(sym2, "NNA");

    if (at3 < cube.natm)
        getAtomicSymbol(cube.zatm[at3], 6, sym3);
    else
        strcpy(sym3, "NNA");

    fprintf(out, "#  File of information in 2D for the plane %d-%d-%d\n",
            at1 + 1, at2 + 1, at3 + 1);
    fprintf(out,
            "#                         3D Info                   2D Info\n");

    trans3Dto2D(r1, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(out,
            "#  %5d %6s (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
            at1 + 1, sym1, r1[0] * B2A, r1[1] * B2A, r1[2] * B2A, ejeM * B2A,
            ejeL * B2A);

    trans3Dto2D(r2, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(out,
            "#  %5d %6s (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
            at2 + 1, sym2, r2[0] * B2A, r2[1] * B2A, r2[2] * B2A, ejeM * B2A,
            ejeL * B2A);

    trans3Dto2D(r3, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(out,
            "#  %5d %6s (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
            at3 + 1, sym3, r3[0] * B2A, r3[1] * B2A, r3[2] * B2A, ejeM * B2A,
            ejeL * B2A);

    trans3Dto2D(rcg, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(
        out,
        "#  Baricentre   (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
        rcg[0] * B2A, rcg[1] * B2A, rcg[2] * B2A, ejeM * B2A, ejeL * B2A);

    trans3Dto2D(r0, vL, vM, cgL, cgM, &ejeL, &ejeM);
    fprintf(
        out,
        "#  Origin Plane (% 7.3lf % 7.3lf % 7.3lf) -> (% 7.3lf, % 7.3lf) [A]\n",
        r0[0] * B2A, r0[1] * B2A, r0[2] * B2A, ejeM * B2A, ejeL * B2A);

    fprintf(vec, "#  File aux of information in 2D for the plane %d-%d-%d\n",
            at1 + 1, at2 + 1, at3 + 1);
    fprintf(vec,
            "#      coor l       coor m       vec l        vec m        fun ");

    double min2[2], max2[2];
    min2[0] = ejeL;
    min2[1] = ejeM;
    max2[0] = -ejeL;
    max2[1] = -ejeM;

    fprintf(out, "#     L axis       M axis      field eval\n");

    int id, nline;

    dataSLine *matrix;

    matrix = (dataSLine *)malloc(n * n * sizeof(dataSLine));
    if (matrix == NULL) {
        printf(" Failed to allocate memory: [%s]\n", "matrix");
        exit(EXIT_FAILURE);
    }

    double basisLMN[4][3];

    basisLMN[0][0] = vL[0];
    basisLMN[0][1] = vL[1];
    basisLMN[0][2] = vL[2];
    basisLMN[1][0] = vM[0];
    basisLMN[1][1] = vM[1];
    basisLMN[1][2] = vM[2];
    basisLMN[2][0] = vN[0];
    basisLMN[2][1] = vN[1];
    basisLMN[2][2] = vN[2];
    basisLMN[3][0] = cgL;
    basisLMN[3][1] = cgM;
    basisLMN[3][2] = cgN;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            id = i * n + j;

            r[0] = r0[0] + i * h * vL[0] + j * h * vM[0];
            r[1] = r0[1] + i * h * vL[1] + j * h * vM[1];
            r[2] = r0[2] + i * h * vL[2] + j * h * vM[2];
            trans3Dto2D(r, vL, vM, cgL, cgM, &ejeL, &ejeM);

            getRiU(r, matT, q);
            getIndex3D(cube.pts, q, cube.min, cube.hvec, idx);
            loadLocalField(idx, cube, param, xx, yy, zz, f);
            gradient3DLog(q[0], q[1], q[2], xx, yy, zz, f, param.pol,
                          param.orth, matU, min0, val);

            matrix[id].i = i;
            matrix[id].j = j;
            matrix[id].f = getGrd(val);
            matrix[id].sline = 0;
            matrix[id].r[0] = r[0];
            matrix[id].r[1] = r[1];
            matrix[id].r[2] = r[2];
            matrix[id].s[0] = ejeL;
            matrix[id].s[1] = ejeM;

            if (matrix[id].f < 1.E-7 && fabs(val[0]) < 1.E-7)
                matrix[id].sline = -1;
        }
    }

    nline = 1;
    id = firstSeed(n, matrix);
    getStreamLine(id, n, nline, min0, cube, param, config, matU, matT, basisLMN,
                  min2, max2, matrix, out, vec);

    int full = checkFullMatrix(n, matrix);
    int jter = 0;
    while (full != 0 && jter < n * n) {
        id = findNeighbours(nline, n, matrix);
        getStreamLine(id, n, nline, min0, cube, param, config, matU, matT,
                      basisLMN, min2, max2, matrix, out, vec);
        full = checkFullMatrix(n, matrix);
        jter++;
        // printf(" Full Matrix : % 10d\n",full);
    }
    // printMatrixGeom( n, matrix);

    free(matrix);
    fclose(out);
    fclose(vec);

    printBar(stdout);
    printf("  File %s was generated\n", nameOut);
    printBar(stdout);
    printf("  File %s was generated\n", nameVec);

    return 0;
}

int firstSeed(int n, dataSLine *matrix) {

    int i, j;
    int id;
    srand(time(NULL));
    do {
        i = rand() % n;
        j = rand() % n;

        id = i * n + j;
    } while (matrix[id].f > 0.5);

    return id;
}

double EulerMethod(double rin[3], double rout[3], dataCube cube, dataRun param,
                   const double *matU, const double *matT, double min0,
                   double eps, int sign) {

    int i;
    double q[3];
    double val[10], vec[3];

    if (abs(sign) > 1)
        sign /= abs(sign);

    getRiU(rin, matT, q);
    numCritical01Vec(q, cube, param, matU, min0, val);
    getKnRungeKuta(vec, val);

    for (i = 0; i < 3; i++) {
        rout[i] = rin[i] + eps * sign * vec[i];
    }

    getRiU(rout, matT, q);
    numCritical01Vec(q, cube, param, matU, min0, val);

    return val[0];
}

double RungeKutta4_fast(double rin[3], double rout[3], dataCube cube,
                        dataRun param, const double *matU, const double *matT,
                        double min0, double eps, int sign) {

    int i;
    double r0[3], r[3], q[3];
    double k1[3], k2[3], k3[3], k4[3];
    double val[10], vec[3];

    if (abs(sign) > 1)
        sign /= abs(sign);

    cpyVec3(rin, r0);

    getRiU(r0, matT, q);
    numCritical01Vec(q, cube, param, matU, min0, val);
    getKnRungeKuta(k1, val);

    for (i = 0; i < 3; i++)
        r[i] = r0[i] + sign * eps / 2.;
    getRiU(r, matT, q);
    numCritical01Vec(q, cube, param, matU, min0, val);
    getKnRungeKuta(k2, val);

    for (i = 0; i < 3; i++)
        r[i] = r0[i] + sign * eps / 2.;
    getRiU(r, matT, q);
    numCritical01Vec(q, cube, param, matU, min0, val);
    getKnRungeKuta(k3, val);

    for (i = 0; i < 3; i++)
        r[i] = r0[i] + sign * eps;
    getRiU(r, matT, q);
    numCritical01Vec(q, cube, param, matU, min0, val);
    getKnRungeKuta(k4, val);

    for (i = 0; i < 3; i++) {
        vec[i] = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.;
        rout[i] = rin[i] + eps * sign * vec[i];
    }

    getRiU(rout, matT, q);
    numCritical01Vec(q, cube, param, matU, min0, val);

    return val[0];
}

double RungeKutta4(double rin[3], double rout[3], dataCube cube, dataRun param,
                   const double *matU, const double *matT, double min0,
                   double eps, int sign) {

    int i;
    double r0[3], r[3], q[3];
    double k1[3], k3[3], k4[3], k6[3];
    double a3 = 0.300 * eps;
    double a4 = 0.600 * eps;
    double a6 = 0.875 * eps;
    double c1 = 0.097883598;
    double c3 = 0.40257649;
    double c4 = 0.21043771;
    double c6 = 0.289102202;
    double val[10], vec[3];

    if (abs(sign) > 1)
        sign /= abs(sign);

    cpyVec3(rin, r0);

    getRiU(r0, matT, q);
    numCritical01Vec(q, cube, param, matU, min0, val);
    getKnRungeKuta(k1, val);
    for (i = 0; i < 3; i++)
        r[i] = r0[i] + sign * a3;

    getRiU(r, matT, q);
    numCritical01Vec(q, cube, param, matU, min0, val);
    getKnRungeKuta(k3, val);
    for (i = 0; i < 3; i++)
        r[i] = r0[i] + sign * a4;

    getRiU(r, matT, q);
    numCritical01Vec(q, cube, param, matU, min0, val);
    getKnRungeKuta(k4, val);
    for (i = 0; i < 3; i++)
        r[i] = r0[i] + sign * a6;

    getRiU(r, matT, q);
    numCritical01Vec(q, cube, param, matU, min0, val);
    getKnRungeKuta(k6, val);

    for (i = 0; i < 3; i++) {
        vec[i] = c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c6 * k6[i];
        rout[i] = rin[i] + eps * sign * vec[i];
    }

    getRiU(rout, matT, q);
    numCritical01Vec(q, cube, param, matU, min0, val);

    return val[0];
}

int getStreamLine(int id0, int n, int nline, double min0, dataCube cube,
                  dataRun param, dataRC config, const double *matU,
                  const double *matT, double basisLMN[4][3], double min[2],
                  double max[2], dataSLine *matrix, FILE *out, FILE *vec) {

    int goon;
    int iter;
    int sign;
    double val[10];
    double ri[3], rout[3];
    double coorL, coorM;
    double fi, fj;

    double q[3], gvec[3], gnorm, gL, gM;
    double h = 1. / (config.geom_npua2);

    fi = matrix[id0].f;
    ri[0] = matrix[id0].r[0];
    ri[1] = matrix[id0].r[1];
    ri[2] = matrix[id0].r[2];

    trans3Dto2D(ri, basisLMN[0], basisLMN[1], basisLMN[3][0], basisLMN[3][1],
                &coorL, &coorM);
    fprintf(out, " % 12.6lf % 12.6lf % 14.8lf\n", coorM * B2A, coorL * B2A, fi);

    goon = 0;
    iter = 0;
    while (goon == 0 && iter < config.geom_maxiter) {

        fj = RungeKutta4(ri, rout, cube, param, matU, matT, min0,
                         config.geom_eps, 1);

        trans3Dto2D(rout, basisLMN[0], basisLMN[1], basisLMN[3][0],
                    basisLMN[3][1], &coorL, &coorM);

        if (iter == 2)
            sign = (fj - fi) / fabs(fj - fi);

        if (iter > 2 && sign * (fj - fi) < 0) {
            iter = 2 * config.geom_maxiter;
            break;
        }
        fi = fj;
        fprintf(out, " % 12.6lf % 12.6lf % 14.8lf\n", coorM * B2A, coorL * B2A,
                fi);

        trans2Dto3D(coorL, coorM, basisLMN, ri);

        goon = check(coorL, coorM, min, max, val);
        goon += checkMatrix(n, matrix, coorL, coorM, nline);
        iter++;

        if (iter % 50 == 0) {
            getRiU(rout, matT, q);

            numCritical01Vec(q, cube, param, matU, min0, val);
            gnorm = getNorm(val[1], val[2], val[3]);

            gvec[0] = val[1] / gnorm;
            gvec[1] = val[2] / gnorm;
            gvec[2] = val[3] / gnorm;

            gL = dotProduct(gvec, basisLMN[0]);
            gM = dotProduct(gvec, basisLMN[1]);

            gnorm = sqrt(gL * gL + gM * gM);

            gL *= h / gnorm;
            gM *= h / gnorm;

            fprintf(vec, " % 12.6lf % 12.6lf % 12.6lf % 12.6lf % 14.8lf\n",
                    coorM * B2A, coorL * B2A, gM * B2A, gL * B2A, val[0]);
        }
    }
    fprintf(out, "\n");

    fi = matrix[id0].f;
    ri[0] = matrix[id0].r[0];
    ri[1] = matrix[id0].r[1];
    ri[2] = matrix[id0].r[2];

    trans3Dto2D(ri, basisLMN[0], basisLMN[1], basisLMN[3][0], basisLMN[3][1],
                &coorL, &coorM);
    fprintf(out, " % 12.6lf % 12.6lf % 14.8lf\n", coorM * B2A, coorL * B2A, fi);

    goon = 0;
    iter = 0;
    while (goon == 0 && iter < config.geom_maxiter) {

        fj = RungeKutta4(ri, rout, cube, param, matU, matT, min0,
                         config.geom_eps, -1);

        trans3Dto2D(rout, basisLMN[0], basisLMN[1], basisLMN[3][0],
                    basisLMN[3][1], &coorL, &coorM);

        if (iter == 2)
            sign = (fj - fi) / fabs(fj - fi);

        if (iter > 2 && sign * (fj - fi) < 0) {
            iter = 2 * config.geom_maxiter;
            break;
        }
        fi = fj;
        fprintf(out, " % 12.6lf % 12.6lf % 14.8lf\n", coorM * B2A, coorL * B2A,
                fi);

        trans2Dto3D(coorL, coorM, basisLMN, ri);

        goon = check(coorL, coorM, min, max, val);
        goon += checkMatrix(n, matrix, coorL, coorM, nline);
        iter++;

        if (iter % 50 == 0) {
            getRiU(rout, matT, q);

            numCritical01Vec(q, cube, param, matU, min0, val);
            gnorm = getNorm(val[1], val[2], val[3]);

            gvec[0] = val[1] / gnorm;
            gvec[1] = val[2] / gnorm;
            gvec[2] = val[3] / gnorm;

            gL = dotProduct(gvec, basisLMN[0]);
            gM = dotProduct(gvec, basisLMN[1]);

            gnorm = sqrt(gL * gL + gM * gM);

            gL *= h / gnorm;
            gM *= h / gnorm;

            fprintf(vec, " % 12.6lf % 12.6lf % 12.6lf % 12.6lf % 14.8lf\n",
                    coorM * B2A, coorL * B2A, gM * B2A, gL * B2A, val[0]);
        }
    }
    fprintf(out, "\n");

    return 0;
}

// Check the point given by cooL and coorM is inside the region
// [min:max][min:max] and the value of the gradient norm will be greater
// than 1.E-7
int check(double coorL, double coorM, double min[2], double max[2],
          double val[10]) {
    int ret = 0;

    if (coorL < min[0] || coorL > max[0])
        ret += 1;
    if (coorM < min[1] || coorM > max[1])
        ret += 1;

    if (getNorm(val[1], val[2], val[3]) < 1.E-7)
        ret += 1;

    return ret;
}

int checkMatrix(int n, dataSLine *matrix, double coorL, double coorM, int id0) {
    int id;
    int val;
    double dist, distmin;
    int idmin;
    idmin = 0;
    distmin = 100;

    for (id = 0; id < n * n; id++) {
        dist = sqrt(pow(coorL - matrix[id].s[0], 2) +
                    pow(coorM - matrix[id].s[1], 2));
        if (dist < distmin) {
            idmin = id;
            distmin = dist;
        }
    }

    val = 0;

    if (matrix[idmin].sline == 0 || matrix[idmin].sline == id0)
        matrix[idmin].sline = id0;
    else
        val = 1;

    return val;
}
void printMatrixGeom(int n, dataSLine *matrix) {
    int i, j, id;

    for (i = n - 1; i >= 0; i--) {
        for (j = 0; j < n; j++) {
            id = i * n + j;
            if (matrix[id].sline == 0)
                printf("[  ]");
            else
                printf("[%2d]", matrix[id].sline % 10);
        }
        printf("\n");
    }
}

int checkFullMatrix(int n, dataSLine *matrix) {
    int id;
    int occ = 0;

    for (id = 0; id < n * n; id++) {
        if (matrix[id].sline != 0) {
            occ++;
        }
    }

    return (n * n - occ);
}

int findNeighbours(int nline, int n, dataSLine *matrix) {

    int i, j, id;
    int nl, nwId;
    int mu, nu;
    int val;

    nl = 1;
    int flag = 0;
    while (nl <= nline && flag == 0) {
        i = 1;
        while (i < n - 1 && flag == 0) {
            j = 1;
            while (j < n - 1 && flag == 0) {
                id = i * n + j;

                if (matrix[id].sline == nl) {
                    mu = matrix[id].i;
                    nu = matrix[id].j;
                    nwId = getIdfromMatrix(mu, nu, n, matrix);
                    if (nwId >= 0) {
                        val = nwId;
                        flag = 1;
                    }
                }
                j++;
            }
            i++;
        }
        nline++;
    }
    return val;
}

int getIdfromMatrix(int mu, int nu, int n, dataSLine *matrix) {
    int i;
    int id[8];

    id[0] = (mu - 1) * n + (nu - 1);
    id[1] = (mu - 1) * n + (nu);
    id[2] = (mu - 1) * n + (nu + 1);
    id[3] = (mu)*n + (nu - 1);
    id[4] = (mu)*n + (nu + 1);
    id[5] = (mu + 1) * n + (nu - 1);
    id[6] = (mu + 1) * n + (nu);
    id[7] = (mu + 1) * n + (nu + 1);

    for (i = 0; i < 8; i++) {
        if (matrix[id[i]].sline == 0)
            return id[i];
    }
    return -1;
}

void trans3Dto2D(double ri[3], double vL[3], double vM[3], double cgL,
                 double cgM, double *coorL, double *coorM) {

    double cL, cM;

    cL = dotProduct(ri, vL);
    cL -= cgL;
    cM = dotProduct(ri, vM);
    cM -= cgM;

    (*coorL) = cL;
    (*coorM) = cM;
}

void trans2Dto3D(double coorL, double coorM, double basisLMN[4][3],
                 double rout[3]) {

    double vL[3], vM[3], vN[3], cg[3];

    vL[0] = basisLMN[0][0];
    vL[1] = basisLMN[0][1];
    vL[2] = basisLMN[0][2];
    vM[0] = basisLMN[1][0];
    vM[1] = basisLMN[1][1];
    vM[2] = basisLMN[1][2];
    vN[0] = basisLMN[2][0];
    vN[1] = basisLMN[2][1];
    vN[2] = basisLMN[2][2];
    cg[0] = basisLMN[3][0];
    cg[1] = basisLMN[3][1];
    cg[2] = basisLMN[3][2];

    double alpha = coorL + cg[0];
    double beta = coorM + cg[1];
    double gamma = 0. + cg[2];

    double den = vN[0] * (vL[2] * vM[1] - vL[1] * vM[2]) -
                 vN[1] * (vL[2] * vM[0] - vL[0] * vM[2]) +
                 vN[2] * (vL[1] * vM[0] - vL[0] * vM[1]);

    double a = alpha * (vM[1] * vN[2] - vM[2] * vN[1]) +
               beta * (vL[2] * vN[1] - vL[1] * vN[2]) +
               gamma * (vL[1] * vM[2] - vL[2] * vM[1]);

    double b = alpha * (vM[2] * vN[0] - vM[0] * vN[2]) +
               beta * (vL[0] * vN[2] - vL[2] * vN[0]) +
               gamma * (vL[2] * vM[0] - vL[0] * vM[2]);

    double c = alpha * (vM[0] * vN[1] - vM[1] * vN[0]) +
               beta * (vL[1] * vN[0] - vL[0] * vN[1]) +
               gamma * (vL[0] * vM[1] - vL[1] * vM[0]);

    rout[0] = -a / den;
    rout[1] = -b / den;
    rout[2] = -c / den;
}
