#include "basins.h"
#include "fields.h"
#include "utils.h"

#include <omp.h>

int Sign ( double dx ){
    if( dx > 0 )
        return +1;
    else if( dx < 0 )
            return -1;
         else
            return 0;
}

void createArrayCells(int n, dataCells **ptr, const char *mess){
    int i;
    if( n < 1){
        printf(" There is an error in the size for [%s]\n",mess);
        exit(EXIT_FAILURE);
    }

    *ptr = (dataCells*) malloc (n * sizeof(dataCells));
    if( *ptr == NULL){
        printf(" Failed to allocate memory: [%s]\n",mess);
        exit(EXIT_FAILURE);
    }


    for(i=0; i<n; i++)
        (*ptr)[i].attr = SET_NULL;

}

void evalBasins     (dataCube cube, dataRun param, double *matU, char *name){

    printf("%sEntramos aqui\n",TGBI);

    if ( param.pbc == YES )
        getCellsPer   (cube,param,matU);
    else
        getCellsNoPer (cube,param,matU);

    printf("Salimos de evalBasins%s \n",TRST);

}


void getCellsPer   (dataCube cube, dataRun param, double *matU){

    printf("Parte periodica\n");
    

}

int screening      (int npt, dataCells *cellsIn, dataCells *cellsOut){
    int i,j;
    j=0;
    for(i=0;i<npt; i++){
        if( cellsIn[i].attr == SET_VOID){
            cellsOut[j] = cellsIn[i];
            j++;
        }
    }

    return j-1;
}

void getCellsNoPer (dataCube cube, dataRun param, double *matU){
    printf("PArte no periodica\n");
    
    int i,j,k;
    int nx = cube.pts[0];
    int ny = cube.pts[1];
    int nz = cube.pts[2];
    int n1 = ny*nz;
    int n2 = nz;

    int idxGlobal;

    dataCells *cbasins;
    dataCells *cbasins2;

    int npt = nx * ny * nz;
    int npt2 = 0;
    double gnorm,val[10];

    createArrayCells(npt,&cbasins,"Basins cells");

#pragma omp parallel private(idxGlobal, i,j,k,gnorm,val)\
                     shared(npt,npt2,n1,n2,cube,param,matU,nx,ny,nz)
{
#pragma omp single
    printf(" Number of threads for evaluation of fields : %4d\n",omp_get_num_threads());
#pragma omp barrier

#pragma omp for schedule (dynamic)
    for(idxGlobal = 0; idxGlobal < npt; idxGlobal++){
        i = idxGlobal/n1;
        j = (idxGlobal - i*n1)/n2;
        k = idxGlobal%n2;


        if(    ( i > 0 && i < nx - 1) 
            && ( j > 0 && j < ny - 1) 
            && ( k > 0 && k < nz - 1) ){ 
            
            loadFieldBasins(i,j,k,n1,n2,cube.hvec,cube.field,val);

            if( param.orth != YES)
                trans02(val,matU);

            gnorm = getGrd(val);

            cbasins[idxGlobal].idx = idxGlobal;
            cbasins[idxGlobal].fun0 = val[0];
            cbasins[idxGlobal].fun1 = getKin(val);
            cbasins[idxGlobal].fun2 = getLap(val);
            cbasins[idxGlobal].dx   = cube.hvec[0]*(val[1] / gnorm);
            cbasins[idxGlobal].dy   = cube.hvec[1]*(val[2] / gnorm);
            cbasins[idxGlobal].dz   = cube.hvec[2]*(val[3] / gnorm);
            cbasins[idxGlobal].norm = gnorm;
            cbasins[idxGlobal].attr = SET_VOID;

            cbasins[idxGlobal].r[0] = cube.min[0] + i*cube.hvec[0];
            cbasins[idxGlobal].r[1] = cube.min[1] + j*cube.hvec[1];
            cbasins[idxGlobal].r[2] = cube.min[2] + k*cube.hvec[2];

            if( cbasins[idxGlobal].fun0 < TOLERANCE ){
                cbasins[idxGlobal].attr = SET_NULL;
            }
        }
#pragma omp critical
{
        if( cbasins[idxGlobal].attr != SET_NULL )
            npt2++;
}
            
        
    }
}//end omp
    
    npt2--;

    printf("npt = % 10d \n",npt);
    printf("npt = % 10d \n",npt2);
    createArrayCells(npt2,&cbasins2,"Basins cells");

    int npOut = screening(npt,cbasins,cbasins2);

    for(i=0;i<10;i++)
        printf(" %4d  %4d % 10.6lf \n",cbasins2[i].attr,cbasins2[i].idx,cbasins2[i].fun0);

    if( npOut == npt2 )
        printf("Vamos bien\n");

    mergeSort(cbasins2,0,npt2);

    for(i=0;i<10;i++)
        printf(" %4d  %4d % 10.6lf \n",cbasins2[i].attr,cbasins2[i].idx,cbasins2[i].fun0);


    assignAttr(npt2,cbasins2,cube);

    free(cbasins);
    free(cbasins2);
}


void loadFieldBasins(int i, int j, int k, int n1, int n2,
                     double *hvec, double *field, double *val){
    int p,q,r,mu;
    double f[27];

    mu = 0;
    for( p= i - 1; p <= i + 1; p++)
        for( q= j - 1; q <= j + 1; q++)
            for( r= k - 1; r <= k + 1; r++){
                f[mu] = field[p*n1 + q*n2 + r];
                mu++;
            }

    hessPol02(hvec[0],hvec[1],hvec[2],f,val);

}


void mergeSort(dataCells *array, int l, int r ){
    if( l < r){
        int m = l + (r-l)/2;

        mergeSort(array,l,m);
        mergeSort(array,m+1,r);

        merge(array,l,m,r);
    }
}

void merge(dataCells *array,int l, int m, int r){
    int i,j,k;
    int n1 = m - l + 1;
    int n2 = r - m;

    //temporal arrays
    dataCells *arrayLeft;
    dataCells *arrayRight;
    createArrayCells(n1,&arrayLeft ,"Basins cells");
    createArrayCells(n2,&arrayRight,"Basins cells");


    for( i = 0; i < n1; i++)
        arrayLeft[i] = array[l+i];
    for( j = 0; j < n2; j++)
        arrayRight[j] = array[m + 1 + j];

    i = 0;
    j = 0;
    k = l;
    while( i < n1 && j < n2){
        if( arrayLeft[i].fun0 >= arrayRight[j].fun0){
            array[k] = arrayLeft[i];
            i++;
        }else{
            array[k] = arrayRight[j];
            j++;
        }
        k++;
    }

    while( j < n2){
        array[k] = arrayRight[j];
        j++;
        k++;
    }

    while( i < n1){
        array[k] = arrayLeft[i];
        i++;
        k++;
    }

    free(arrayRight);
    free(arrayLeft);

}
void assignAttr(int npt, dataCells *cbasins,dataCube cube){
    double r[3];
    double hx,hy,hz;
    double x0,y0,z0;
    int  i,j,k;
    int at;
    int idx;
    int n1 = cube.pts[1]*cube.pts[2];
    int n2 = cube.pts[2];

    x0 = cube.min[0];
    y0 = cube.min[1];
    z0 = cube.min[2];

    hx = cube.hvec[0];
    hy = cube.hvec[1];
    hz = cube.hvec[2];

    double h = (hx + hy + hz)/3;
    double dij;

    idx = cbasins[0].idx;

    i = idx/n1;
    j = (idx - i*n1)/n2;
    k = idx%n2;

    int nwi,nwj,nwk,nwidx;

    at = getAttr( cube.natm,0,h,cbasins[0].r,cube.coor);
    printf("\n");
    cbasins[0].attr = at;


    for(i=1; i<npt; i++){
        idx = cbasins[i].idx;

        nwi = idx/n1;
        nwj = (idx - nwi*n1)/n2;
        nwk = idx%n2;

        //r[0] = x0 + nwi*hx + cbasins[i].dx;
        //r[1] = y0 + nwj*hx + cbasins[i].dy;
        //r[2] = z0 + nwk*hx + cbasins[i].dz;

        for(j=0; j < i; j++){

            dij = distance(cbasins[i].r,cbasins[j].r);

            if( dij<= 1.2*h ){
                cbasins[i].attr = cbasins[j].attr;
                break;
            }
            
        }
        if( cbasins[i].attr < 0) 
            cbasins[i].attr = getAttr(cube.natm,0,h,cbasins[i].r,cube.coor);

    }


    for(i=0;i<npt;i++){
        printf(" BAS%d % 10.6lf % 10.6lf % 10.6lf\n",cbasins[i].attr,cbasins[i].r[0]*B2A,
        cbasins[i].r[1]*B2A, cbasins[i].r[2]*B2A);
    }

}


int getAttr( int natm,int nna,double h, double r[3], double *rcoor){
    int attr;
    double rc[3];
    double min_dist = 1.E7;
    double dmu;
    int mu;

    for(mu=0;mu<natm; mu++){
        rc[0] = rcoor[3*mu];
        rc[1] = rcoor[3*mu+1];
        rc[2] = rcoor[3*mu+2];
        dmu = distance(r,rc);
        if( dmu < min_dist){
            min_dist = dmu;
        }
        if( dmu < 1.25*h){
            attr = mu;
            break;
        }
    }

    return attr;
}


