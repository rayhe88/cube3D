#include "replicate.h"
#include "array.h"
#include "file.h"
#include "lectura.h"
#include "struct.h"
#include "utils.h"

int replicate(dataCube cubeOld, int rep[3], const char *name) {

    char text[64];
    int i, npt2, check;
    int tamanio;
    int newxp = ceil(cubeOld.pts[0] * rep[0]);
    int newyp = ceil(cubeOld.pts[1] * rep[1]);
    int newzp = ceil(cubeOld.pts[2] * rep[2]);
    int *zatm2;
    double *coor2;
    double *fieldNew;

    FILE *out;
    dataCube cubeNew;

    sprintf(text, " Cube replicate %3d x %3d x %3d times", rep[0], rep[1],
            rep[2]);

    openFile(&out, name, "w+");

    printf("%s", TRBI);
    printf(" Original points\n");
    printf("  xn = % 6d\n", cubeOld.pts[0]);
    printf("  yn = % 6d\n", cubeOld.pts[1]);
    printf("  zn = % 6d\n", cubeOld.pts[2]);
    printf("%s", TRST);

    checkBoundaryCond(cubeOld, &check);

    printf(" Check value for boundary conditions: %d\n", check);

    newxp -= (rep[0] - 1) * check;
    newyp -= (rep[1] - 1) * check;
    newzp -= (rep[2] - 1) * check;
    printf("%s", TRBI);
    printf(" Replicate points\n");
    printf("  xn = % 6d\n", newxp);
    printf("  yn = % 6d\n", newyp);
    printf("  zn = % 6d\n", newzp);
    printf("%s", TRST);

    npt2 = newxp * newyp * newzp;

    createArrayDou(npt2, &fieldNew, "Field2");

    cubeNew.pts[0] = newxp;
    cubeNew.pts[1] = newyp;
    cubeNew.pts[2] = newzp;
    cubeNew.npt = npt2;

    cubeNew.natm = cubeOld.natm;
    for (i = 0; i < 3; i++) {
        cubeNew.min[i] = cubeOld.min[i];
        cubeNew.hvec[i] = cubeOld.hvec[i];
    }
    for (i = 0; i < 9; i++)
        cubeNew.mvec[i] = cubeOld.mvec[i];

    for (i = 0; i < npt2; i++) {
        fieldNew[i] =
            cubeOld.field[getIndexOld(check, i, cubeOld.pts, cubeNew.pts)];
    }

    tamanio = cubeOld.natm;
    tamanio *= ceil(rep[0]);
    tamanio *= ceil(rep[1]);
    tamanio *= ceil(rep[2]);

    createArrayInt(tamanio, &zatm2, "New zatm");
    createArrayDou(3 * tamanio, &coor2, "New coor");

    replicateCoor(check, cubeOld, zatm2, coor2, rep);

    cubeNew.natm = tamanio;

    cubeNew.zatm = zatm2;
    cubeNew.coor = coor2;
    cubeNew.field = fieldNew;

    printCube(text, cubeNew, out);

    unloadData(&cubeNew, &zatm2, &coor2, &fieldNew);

    // checkData(
    // cubeNew.natm,cubeNew.pts,zatm2,cubeNew.min,cubeNew.mvec,coor2,fieldNew,"prueba
    // replicada",out);

    fclose(out);

    return 0;
}

int getIndexOld(int check, int idx, int oldpts[3], int newpts[3]) {

    int oldx, oldy, oldz;
    int newx, newy, newz;
    int oldnpx, oldnpy, oldnpz;
    int newnpy, newnpz;

    int oldIdx;

    oldnpx = oldpts[0];
    oldnpy = oldpts[1];
    oldnpz = oldpts[2];

    newnpy = newpts[1];
    newnpz = newpts[2];

    newx = (idx / newnpz) / newnpy;
    newy = (idx / newnpz) % newnpy;
    newz = (idx % newnpz);

    // NOTA USAR -1 para cubos de Crystal
    //            0 para cubos de GPAW
    // esto viene dado por check

    oldx = newx % (oldnpx - check);
    oldy = newy % (oldnpy - check);
    oldz = newz % (oldnpz - check);

    oldIdx = oldx * oldnpz * oldnpy + oldy * oldnpz + oldz;

    return oldIdx;
}

int replicateCoor(int check, dataCube cubeOld, int *zatm2, double *coor2,
                  int *rep) {

    int i, j, k, l, n;
    int nrx, nry, nrz;

    int natm = cubeOld.natm;
    int npx = cubeOld.pts[0] - check;
    int npy = cubeOld.pts[1] - check;
    int npz = cubeOld.pts[2] - check;

    double hrx[3], hry[3], hrz[3];

    hrx[0] = npx * (cubeOld.mvec[0]);
    hrx[1] = npx * (cubeOld.mvec[1]);
    hrx[2] = npx * (cubeOld.mvec[2]);

    hry[0] = npy * (cubeOld.mvec[3]);
    hry[1] = npy * (cubeOld.mvec[4]);
    hry[2] = npy * (cubeOld.mvec[5]);

    hrz[0] = npz * (cubeOld.mvec[6]);
    hrz[1] = npz * (cubeOld.mvec[7]);
    hrz[2] = npz * (cubeOld.mvec[8]);

    nrx = (int)rep[0];
    nry = (int)rep[1];
    nrz = (int)rep[2];

    for (i = 0; i < nrx; i++) {
        for (j = 0; j < nry; j++) {
            for (k = 0; k < nrz; k++) {
                for (n = 0; n < natm; n++) {
                    l = i * nry * nrz * natm;
                    l += j * nrz * natm;
                    l += k * natm;
                    l += n;
                    zatm2[l] = cubeOld.zatm[n];
                    coor2[3 * l] = cubeOld.coor[3 * n] + i * hrx[0] +
                                   j * hry[0] + k * hrz[0];
                    coor2[3 * l + 1] = cubeOld.coor[3 * n + 1] + i * hrx[1] +
                                       j * hry[1] + k * hrz[1];
                    coor2[3 * l + 2] = cubeOld.coor[3 * n + 2] + i * hrx[2] +
                                       j * hry[2] + k * hrz[2];
                }
            }
        }
    }

    return 0;
}
