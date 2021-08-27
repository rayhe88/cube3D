#include "numGradLines.h"

/**
 * @file   numGradLines.h
 * @brief
 * @author Raymundo Hernández-Esparza.
 * @date   June 2021.
 */

#include "array.h"
#include "fields.h"
#include "file.h"
#include "findCrit.h"
#include "jacobi.h"
#include "lagrange2.h"
#include "mathTools.h"
#include "numBondPath.h"
#include "struct.h"
#include "tableP.h"
#include "transU.h"
#include "utils.h"
#include "version.h"

#include "lebedev.h"
#include "numGradLines.h"

#include <omp.h>

void transVec(dataVec in, double out[3]) {
    out[0] = in.x;
    out[1] = in.y;
    out[2] = in.z;
}

int gradientLines(int ncp, dataCritP *nnucCrit, dataCube cube, dataRun param,
                  double min0, const double *matU, char *name) {
    int np;
    int i, j, k, iter;
    int attractors;
    double val[10];
    char nameOut[128], tmpname[128];

    double rad;
    double ri[3], qi[3], rn[3], qn[3];
    double qc[3];

    double *coorAttr;
    double rho;
    FILE *tmp;
    FILE *out;

    dataVec *rleb;

    /**********************************************************
      TODO  this is only a form to repair the gradient lines
            (bpaths). In the future may do a structure for
            matrix of transformation and rotation.
    **********************************************************/
    double matT[9];
    getMatInv(matU, matT); // The inverse of U is the original matrix T

    sprintf(nameOut, "%sGLine.xyz", name);

    tmpFile(&tmp, ".c3dgline_", tmpname, "w+");

    attractors = ncp + cube.natm;

    createArrayDou(3 * attractors, &coorAttr, "GL01");

    for (i = 0; i < cube.natm; i++) {
        coorAttr[3 * i] = cube.coor[3 * i];
        coorAttr[3 * i + 1] = cube.coor[3 * i + 1];
        coorAttr[3 * i + 2] = cube.coor[3 * i + 2];
    }
    for (j = 0; j < ncp; j++) {
        coorAttr[3 * i] = nnucCrit[j].x;
        coorAttr[3 * i + 1] = nnucCrit[j].y;
        coorAttr[3 * i + 2] = nnucCrit[j].z;
        i++;
    }

    printBar(stdout);

    printf("  Number of threads in bond path  : %6d\n", omp_get_num_threads());

    for (i = 0; i < attractors; i++) {
        np = 86;
        rad = getChemProperty(1, "RADIUS");

        if (i < cube.natm) {
            rad = getChemProperty(cube.zatm[i], "RADIUS");
            if (cube.zatm[i] >= 2) {
                np = 230;
            }
        }
        createArrayVec(np, &rleb, "Lebedev points");
        if (np == 230)
            Lebedev0230(rleb);
        if (np == 86)
            Lebedev0086(rleb);

        qc[0] = coorAttr[3 * i];
        qc[1] = coorAttr[3 * i + 1];
        qc[2] = coorAttr[3 * i + 2];

        numCritical01Vec(qc, cube, param, matU, min0, val);
        rho = val[0];
        if (rho > 0.2)
            rad *= 0.5;

        printf(" Attractor [%3d] (% 10.6lf % 10.6lf % 10.6lf) = % 10.6lf Rad % "
               "10.6lf\n",
               i, qc[0], qc[1], qc[2], rho, rad);

        perfectCube(param.pbc, qc, cube.min, cube.max);
        translateAndChange(np, rleb, qc, rad);

        for (k = 0; k < np; k++) {

            transVec(rleb[k], ri);
            getRiU(ri, matT, qi);
            iter = 0;
            numCritical01Vec(qi, cube, param, matU, min0, val);
            rho = val[0];
            while (iter < MAXPTS && rho > 1.E-7) {
                perfectCube(param.pbc, qi, cube.min, cube.max);

                numCritical01Vec(qi, cube, param, matU, min0, val);
                rho = val[0];

                rn[0] = ri[0] - 0.02 * val[1];
                rn[1] = ri[1] - 0.02 * val[2];
                rn[2] = ri[2] - 0.02 * val[3];

                getRiU(rn, matT, qn);
                perfectCube(param.pbc, qn, cube.min, cube.max);
                cpyVec3(qn, qi);
                getRiU(qi, matU, ri);

                if (myIsNanInf_V3(qn) == 0) {
                    getRiU(qn, matU, rn);
                    fprintf(tmp, " GL%04d  % 10.6lf   % 10.6lf  % 10.6lf\n",
                            i + 1, rn[0] * B2A, rn[1] * B2A, rn[2] * B2A);
                } else {
                    iter = MAXPTS;
                }
                iter++;
            }
        }

        free(rleb);
    }

    rewind(tmp);

    int npoints = 0;
    char c;

    while ((c = fgetc(tmp)) != EOF) {
        if (c == '\n')
            npoints++;
    }
    rewind(tmp);

    openFile(&out, nameOut, "w+");

    fprintf(out, " %10d\n", npoints);
    fprintf(out, " Gradient Lines (basins) file in Aangstrom\n");

    while ((c = fgetc(tmp)) != EOF) {
        fprintf(out, "%c", c);
    }

    fclose(tmp);
    fclose(out);

    remove(tmpname);
    printBar(stdout);
    printf("  File %s was generated\n", nameOut);
    return 0;
}

int gradientLines2(int bcp, dataCritP *bondCrit, dataCube cube, dataRun param,
                   double min0, const double *matU, char *name) {
    int np;
    int i, k, iter;
    double val[10];
    char nameOut[128], tmpname[128];

    double rad;
    double ri[3], qi[3], rn[3], qn[3];
    double qc[3];

    double rho;
    FILE *tmp;
    FILE *out;

    dataVec *rleb;

    /**********************************************************
      TODO  this is only a form to repair the gradient lines
            (bpaths). In the future may do a structure for
            matrix of transformation and rotation.
    **********************************************************/
    double matT[9];
    getMatInv(matU, matT); // The inverse of U is the original matrix T

    sprintf(nameOut, "%sGLine2.xyz", name);

    tmpFile(&tmp, ".c3dgline_", tmpname, "w+");

    printBar(stdout);

    printf("  Number of threads in bond path  : %6d\n", omp_get_num_threads());

    for (i = 0; i < bcp; i++) {
        np = 86;
        rad = getChemProperty(1, "RADIUS");

        createArrayVec(np, &rleb, "Lebedev points");
        if (np == 86)
            Lebedev0086(rleb);

        qc[0] = bondCrit[i].x;
        qc[1] = bondCrit[i].y;
        qc[2] = bondCrit[i].z;

        rad *= 0.05;

        printf(" Attractor [%3d] (% 10.6lf % 10.6lf % 10.6lf) = % 10.6lf Rad % "
               "10.6lf\n",
               i, qc[0], qc[1], qc[2], rho, rad);

        perfectCube(param.pbc, qc, cube.min, cube.max);
        translateAndChange(np, rleb, qc, rad);

        for (k = 0; k < np; k++) {

            transVec(rleb[k], ri);
            getRiU(ri, matT, qi);
            iter = 0;
            numCritical01Vec(qi, cube, param, matU, min0, val);
            rho = val[0];
            while (iter < MAXPTS && rho > 1.E-7) {
                perfectCube(param.pbc, qi, cube.min, cube.max);

                numCritical01Vec(qi, cube, param, matU, min0, val);
                rho = val[0];

                rn[0] = ri[0] - 0.02 * val[1];
                rn[1] = ri[1] - 0.02 * val[2];
                rn[2] = ri[2] - 0.02 * val[3];

                getRiU(rn, matT, qn);
                perfectCube(param.pbc, qn, cube.min, cube.max);
                cpyVec3(qn, qi);
                getRiU(qi, matU, ri);

                if (myIsNanInf_V3(qn) == 0) {
                    getRiU(qn, matU, rn);
                    fprintf(tmp, " GL%04d  % 10.6lf   % 10.6lf  % 10.6lf\n",
                            i + 1, rn[0] * B2A, rn[1] * B2A, rn[2] * B2A);
                } else {
                    iter = MAXPTS;
                }
                iter++;
            }
        }

        free(rleb);
    }

    rewind(tmp);

    int npoints = 0;
    char c;

    while ((c = fgetc(tmp)) != EOF) {
        if (c == '\n')
            npoints++;
    }
    rewind(tmp);

    openFile(&out, nameOut, "w+");

    fprintf(out, " %10d\n", npoints);
    fprintf(out, " Gradient Lines (basins) file in Aangstrom\n");

    while ((c = fgetc(tmp)) != EOF) {
        fprintf(out, "%c", c);
    }

    fclose(tmp);
    fclose(out);

    remove(tmpname);
    printBar(stdout);
    printf("  File %s was generated\n", nameOut);
    return 0;
}
int gradientLines3(int bcp, int rcp, int ccp, dataCritP *bondCrit,
                   dataCritP *ringCrit, dataCritP *cageCrit, dataCube cube,
                   dataRun param, double min0, const double *matU, char *name) {
    int np;
    int i, k, iter;
    double val[10];
    char nameOut[128], tmpname[128];

    double ri[3], qi[3], rn[3], qn[3];
    double qc[3];

    FILE *tmp;
    FILE *out;
    double matH[9], rc[3];
    double eval[3], evec[9], vec1[3], vec2[3], nv1, nv2;

    dataVec *rleb;

    /**********************************************************
      TODO  this is only a form to repair the gradient lines
            (bpaths). In the future may do a structure for
            matrix of transformation and rotation.
    **********************************************************/
    double matT[9];
    getMatInv(matU, matT); // The inverse of U is the original matrix T

    sprintf(nameOut, "%sGLine3.xyz", name);

    tmpFile(&tmp, ".c3dgline_", tmpname, "w+");

    printBar(stdout);

    printf("  Number of threads in bond path  : %6d\n", omp_get_num_threads());
    double theta, dtheta, dr;

    np = 72;
    dtheta = (2. * M_PI) / (double)np;
    dr = 0.1;

    for (i = 0; i < rcp; i++) {
        qc[0] = ringCrit[i].x;
        qc[1] = ringCrit[i].y;
        qc[2] = ringCrit[i].z;

        getRiU(qc, matU, rc);

        numCritical02Vec(qc, cube, param, matU, min0, val);

        matH[0] = val[4];
        matH[1] = val[7];
        matH[2] = val[8];
        matH[3] = val[7];
        matH[4] = val[5];
        matH[2] = val[9];
        matH[6] = val[8];
        matH[7] = val[9];
        matH[2] = val[6];

        JacobiNxN(matH, eval, evec);
        // The 1st eivec is taken
        vec1[0] = evec[3];
        vec1[1] = evec[4];
        vec1[2] = evec[5];
        vec2[0] = evec[6];
        vec2[1] = evec[7];
        vec2[2] = evec[8];

        nv1 = getNormVec(vec1);
        nv2 = getNormVec(vec2);

        vec1[0] /= nv1;
        vec1[1] /= nv1;
        vec1[2] /= nv1;
        vec2[0] /= nv2;
        vec2[1] /= nv2;
        vec2[2] /= nv2;

        for (k = 0; k < np; k++) {
            theta = k * dtheta;
            ri[0] =
                rc[0] + dr * cos(theta) * vec1[0] + dr * sin(theta) * vec2[0];
            ri[1] =
                rc[1] + dr * cos(theta) * vec1[1] + dr * sin(theta) * vec2[1];
            ri[2] =
                rc[2] + dr * cos(theta) * vec1[2] + dr * sin(theta) * vec2[2];

            iter = 0;
            getRiU(ri, matT, qi);
            while (iter < MAXPTS) {
                perfectCube(param.pbc, qi, cube.min, cube.max);

                numCritical01Vec(qi, cube, param, matU, min0, val);

                rn[0] = ri[0] + 0.1 * val[1];
                rn[1] = ri[1] + 0.1 * val[2];
                rn[2] = ri[2] + 0.1 * val[3];

                getRiU(rn, matT, qn);
                perfectCube(param.pbc, qn, cube.min, cube.max);
                cpyVec3(qn, qi);
                getRiU(qi, matU, ri);

                if (myIsNanInf_V3(qn) == 0) {
                    getRiU(qn, matU, rn);
                    fprintf(tmp, " GR%04d  % 10.6lf   % 10.6lf  % 10.6lf\n",
                            i + 1, rn[0] * B2A, rn[1] * B2A, rn[2] * B2A);
                } else {
                    iter = MAXPTS;
                }
                iter++;
            }
        }

        free(rleb);
    }

    rewind(tmp);

    int npoints = 0;
    char c;

    while ((c = fgetc(tmp)) != EOF) {
        if (c == '\n')
            npoints++;
    }
    rewind(tmp);

    openFile(&out, nameOut, "w+");

    fprintf(out, " %10d\n", npoints);
    fprintf(out, " Gradient Lines (basins) file in Aangstrom\n");

    while ((c = fgetc(tmp)) != EOF) {
        fprintf(out, "%c", c);
    }

    fclose(tmp);
    fclose(out);

    remove(tmpname);
    printBar(stdout);
    printf("  File %s was generated\n", nameOut);
    return 0;
}

int gradientLines4(int bcp, int rcp, int ccp, dataCritP *bondCrit,
                   dataCritP *ringCrit, dataCritP *cageCrit, dataCube cube,
                   dataRun param, double min0, const double *matU, char *name) {
    int i, j, iter;
    double val[10];
    char nameOut[128], tmpname[128];

    double ri[3], qi[3], rn[3], qn[3];
    double rj[3], qj[3];
    double qc[3];

    FILE *tmp;
    FILE *out;
    double matH[9], rc[3];
    double eval[3], evec[9], vec0[3], vec1[3], vec2[3], nv0, nv1, nv2;
    double vij[3], nv;

    /**********************************************************
      TODO  this is only a form to repair the gradient lines
            (bpaths). In the future may do a structure for
            matrix of transformation and rotation.
    **********************************************************/
    double matT[9];
    double step;
    getMatInv(matU, matT); // The inverse of U is the original matrix T

    sprintf(nameOut, "%sGLine4.xyz", name);

    tmpFile(&tmp, ".c3dgline_", tmpname, "w+");

    printBar(stdout);

    printf("  Number of threads in bond path  : %6d\n", omp_get_num_threads());
    double cop;
    int k;

    for (i = 0; i < rcp; i++) {
        qc[0] = ringCrit[i].x;
        qc[1] = ringCrit[i].y;
        qc[2] = ringCrit[i].z;

        getRiU(qc, matU, rc);
        numCritical02Vec(qc, cube, param, matU, min0, val);

        matH[0] = val[4];
        matH[1] = val[7];
        matH[2] = val[8];
        matH[3] = val[7];
        matH[4] = val[5];
        matH[5] = val[9];
        matH[6] = val[8];
        matH[7] = val[9];
        matH[8] = val[6];

        JacobiNxN(matH, eval, evec);

        vec0[0] = evec[0];
        vec0[1] = evec[1];
        vec0[2] = evec[2];
        vec1[0] = evec[3];
        vec1[1] = evec[4];
        vec1[2] = evec[5];
        vec2[0] = evec[6];
        vec2[1] = evec[7];
        vec2[2] = evec[8];

        nv0 = getNormVec(vec0);
        nv1 = getNormVec(vec1);
        nv2 = getNormVec(vec2);

        vec0[0] /= nv0;
        vec0[1] /= nv0;
        vec0[2] /= nv0;
        vec1[0] /= nv1;
        vec1[1] /= nv1;
        vec1[2] /= nv1;
        vec2[0] /= nv2;
        vec2[1] /= nv2;
        vec2[2] /= nv2;
        /*
            for(k=0;k<2;k++){
              rn[0] = rc[0] + k * 0.75 * vec0[0];
              rn[1] = rc[1] + k * 0.75 * vec0[1];
              rn[2] = rc[2] + k * 0.75 * vec0[2];
              fprintf(tmp," V1%04d  % 10.6lf   % 10.6lf
           % 10.6lf\n",i+1,rn[0]*B2A,rn[1]*B2A,rn[2]*B2A);
            }
            for(k=0;k<2;k++){
              rn[0] = rc[0] - k * 0.75 * vec0[0];
              rn[1] = rc[1] - k * 0.75 * vec0[1];
              rn[2] = rc[2] - k * 0.75 * vec0[2];
              fprintf(tmp," V1%04d  % 10.6lf   % 10.6lf
           % 10.6lf\n",i+1,rn[0]*B2A,rn[1]*B2A,rn[2]*B2A);
            }
            for(k=0;k<2;k++){
              rn[0] = rc[0] + k * 0.75 * vec1[0];
              rn[1] = rc[1] + k * 0.75 * vec1[1];
              rn[2] = rc[2] + k * 0.75 * vec1[2];
              fprintf(tmp," V2%04d  % 10.6lf   % 10.6lf
           % 10.6lf\n",i+1,rn[0]*B2A,rn[1]*B2A,rn[2]*B2A);
            }
            for(k=0;k<2;k++){
              rn[0] = rc[0] - k * 0.75 * vec1[0];
              rn[1] = rc[1] - k * 0.75 * vec1[1];
              rn[2] = rc[2] - k * 0.75 * vec1[2];
              fprintf(tmp," V2%04d  % 10.6lf   % 10.6lf
           % 10.6lf\n",i+1,rn[0]*B2A,rn[1]*B2A,rn[2]*B2A);
            }
            for(k=0;k<2;k++){
              rn[0] = rc[0] + k * 0.75 * vec2[0];
              rn[1] = rc[1] + k * 0.75 * vec2[1];
              rn[2] = rc[2] + k * 0.75 * vec2[2];
              fprintf(tmp," V3%04d  % 10.6lf   % 10.6lf
           % 10.6lf\n",i+1,rn[0]*B2A,rn[1]*B2A,rn[2]*B2A);
            }
            for(k=0;k<2;k++){
              rn[0] = rc[0] - k * 0.75 * vec2[0];
              rn[1] = rc[1] - k * 0.75 * vec2[1];
              rn[2] = rc[2] - k * 0.75 * vec2[2];
              fprintf(tmp," V3%04d  % 10.6lf   % 10.6lf
           % 10.6lf\n",i+1,rn[0]*B2A,rn[1]*B2A,rn[2]*B2A);
            }
        */
        for (j = 0; j < bcp; j++) {
            qj[0] = bondCrit[j].x;
            qj[1] = bondCrit[j].y;
            qj[2] = bondCrit[j].z;
            getRiU(qj, matU, rj);

            if (distance(rc, rj) < 5.66918) { // 3 Angstrom = 5.66918 Bohr

                vij[0] = rj[0] - rc[0];
                vij[1] = rj[1] - rc[1];
                vij[2] = rj[2] - rc[2];

                nv = getNormVec(vij);
                vij[0] /= nv;
                vij[1] /= nv;
                vij[2] /= nv;

                cop = dotProduct(vec0, vij);
                // printf(" Cop  %3d %3d distance = % 10.6lf = % 10.6lf =
                // % 10.6lf\n",i,j, distance(rc,rj),cop,acos(cop)*180./M_PI);

                if (cop > -0.5 && cop < 0.5) {

                    ri[0] = rc[0] + 0.1 * vij[0];
                    ri[1] = rc[1] + 0.1 * vij[1];
                    ri[2] = rc[2] + 0.1 * vij[2];

                    getRiU(ri, matT, qi);

                    iter = 0;
                    step = 5;
                    while (iter < MAXPTS) {
                        perfectCube(param.pbc, qi, cube.min, cube.max);
                        numCritical01Vec(qi, cube, param, matU, min0, val);

                        rn[0] = ri[0] + 0.2 * val[1];
                        rn[1] = ri[1] + 0.2 * val[2];
                        rn[2] = ri[2] + 0.2 * val[3];

                        getRiU(rn, matT, qn);
                        perfectCube(param.pbc, qn, cube.min, cube.max);
                        cpyVec3(qn, qi);
                        getRiU(qi, matU, ri);
                        if (myIsNanInf_V3(qn) == 0) {
                            getRiU(qn, matU, rn);
                            if (step == 5) {
                                fprintf(
                                    tmp,
                                    " GR%04d  % 10.6lf   % 10.6lf  % 10.6lf\n",
                                    i + 1, rn[0] * B2A, rn[1] * B2A,
                                    rn[2] * B2A);
                                step = 0;
                            }

                            if (distance(rj, rn) < 0.05)
                                iter = 2 * MAXPTS;
                        } else {
                            iter = 2 * MAXPTS;
                        }
                        iter++;
                        step++;
                    }
                }
            }
        }
    }

    rewind(tmp);

    int npoints = 0;
    char c;

    while ((c = fgetc(tmp)) != EOF) {
        if (c == '\n')
            npoints++;
    }
    rewind(tmp);

    openFile(&out, nameOut, "w+");

    fprintf(out, " %10d\n", npoints);
    fprintf(out, " Gradient Lines (basins) file in Aangstrom\n");

    while ((c = fgetc(tmp)) != EOF) {
        fprintf(out, "%c", c);
    }

    fclose(tmp);
    fclose(out);

    remove(tmpname);
    printBar(stdout);
    printf("  File %s was generated\n", nameOut);
    return 0;
}
void cpyCritP(int np, dataCritP *in, dataCritP *out) {
    int i;
    for (i = 0; i < np; i++)
        out[i] = in[i];
}
int sortBondCritP(int np, double qi[3], dataCritP *crit) {

    int i, j;
    double qj[3];
    dataCritP tmp;

    for (j = 0; j < np; j++) {
        qj[0] = crit[j].x;
        qj[1] = crit[j].y;
        qj[2] = crit[j].z;
        crit[j].fun = distance(qi, qj);
    }

    for (i = 0; i < np - 1; i++)
        for (j = i + 1; j < np; j++)
            if (crit[j].fun > crit[i].fun) {
                tmp = crit[i];
                crit[i] = crit[j];
                crit[j] = tmp;
            }
}
