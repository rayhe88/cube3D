#include "basins.h"
#include "fields.h"
#include "utils.h"
#include "findCrit.h"

#include <omp.h>

#define FALSE 0
#define TRUE !FALSE


void decomposeIdx(int idx, int n1, int n2,int *i,int *j,int *k){
    int ti,tj,tk;
    ti = idx/n1;
    tj = (idx - ti*n1)/n2;
    tk = idx%n2;
    (*i) = ti;
    (*j) = tj;
    (*k) = tk;
}

int composeIdx(int i, int j, int k, int n1, int n2){
    return i*n1 + j*n2 + k;
}

void getHash(dataCells *cells, int *hash, int npt){
    int i;
//#pragma omp parallel private(i) shared(npt,hash)
//{
//#pragma omp for schedule (dynamic)
    for(i=0;i<npt;i++)
        hash[i] = 0;
//}

//#pragma omp parallel private(i) shared(npt,hash,cells)
//{
//#pragma omp for schedule (dynamic)
    for(i=0;i<npt;i++){
        hash[cells[i].idx] = i;
    }
//}
}

void checkHash(int npt, dataCells *cells, int *hash){
    int i,sum = 0;
    for(i=0;i<npt;i++){
        if ( hash[cells[i].idx] != i )
            sum++;
    }
    if (sum != 0 )
        printf("[1] Error en la tabla hash en %7d entradas\n",sum);

    sum = 0;
    for(i=0;i<npt;i++){
        if ( cells[hash[i]].idx != i )
            sum++;
    }
    if (sum != 0 )
        printf("[2] Error en la tabla hash en %7d entradas\n",sum);

}

int getIdxSingle(double r, double r0, double h){
    
    int i_up,i_down;
    double id,d_up,d_down;

    id     = (r - r0)/h;
    i_up   = (int) ceil(id);
    i_down = (int) floor(id);

    d_down = fabs(r - (r0 + h * i_down));
    d_up   = fabs(r - (r0 + h * i_up  ));



    if( d_up < d_down )
        return i_up;
    return  i_down;
}

int getIdx(dataCube cube, double r[3]){
    int i,j,k;

    i = getIdxSingle(r[0],cube.min[0],cube.hvec[0]);
    j = getIdxSingle(r[1],cube.min[1],cube.hvec[1]);
    k = getIdxSingle(r[2],cube.min[2],cube.hvec[2]);

    return composeIdx(i,j,k,cube.pts[1]*cube.pts[2],cube.pts[2]);
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

    for(i=0; i<n; i++){
        (*ptr)[i].attr = SET_NULL;
        (*ptr)[i].max  = FALSE;
        (*ptr)[i].idx  = 0;
        (*ptr)[i].fun0 = (double)0.;
        (*ptr)[i].fun1 = (double)0.;
        (*ptr)[i].fun2 = (double)0.;
        (*ptr)[i].dx   = (double)0.;
        (*ptr)[i].dy   = (double)0.;
        (*ptr)[i].dz   = (double)0.;
        (*ptr)[i].norm = (double)0.;
    }

}

void evalBasins     (dataCube cube, dataRun param, double *matU, char *name){

    printf("%s\n We enter in \'evalBains\' \n",TGBI);

    if ( param.pbc == YES )
        getCellsPer   (cube,param,matU);
    else
        getCellsNoPer (cube,param,matU);

    printf(" We out of \'evalBasins\'\n%s\n",TRST);

}


void getCellsPer   (dataCube cube, dataRun param, double *matU){

    printf(" Periodic subroutine\n");
    

}

void getCellsNoPer (dataCube cube, dataRun param, double *matU){
    printf(" Non-periodic subroutine\n");
    
    int i,j,k;
    int npt = cube.npt;
    int n1 = cube.pts[1]*cube.pts[2];
    int n2 = cube.pts[2];
    int natt, max, index;
    int *hash;
    double gnorm,val[10];

    dataCells *cells;

    createArrayCells (npt, &cells, "Basins cells");
    createArrayInt   (npt, &hash ," Basins cells");
    natt = 0;

//#pragma omp parallel private()\\
                     shared()
//{
//#pragma omp single
    printf(" Number of threads for evaluation of fields : %4d\n",omp_get_num_threads());
//#pragma omp barrier

//#pragma omp for schedule (dynamic)
    double qcharge=(double)0.;
    double qcharge0 = (double) 0.;
    for(index = 0; index < npt; index++){
        
        qcharge0 += cube.field[index];

        decomposeIdx(index,n1,n2,&i,&j,&k);

        cells[index].idx = index;
        cells[index].r[0] = cube.min[0] + i*cube.hvec[0];
        cells[index].r[1] = cube.min[1] + j*cube.hvec[1];
        cells[index].r[2] = cube.min[2] + k*cube.hvec[2];
        cells[index].attr = SET_NULL;

        if( ( i > 0 && i < cube.pts[0]-1) && 
            ( j > 0 && j < cube.pts[1]-1) && 
            ( k > 0 && k < cube.pts[2]-1) ){ 
            
            max = loadFieldBasins(i,j,k,n1,n2,cube.hvec,cube.field,val);

            if( param.orth != YES) trans02(val,matU);

            gnorm = getGrd(val);

            cells[index].fun0 = val[0];
            cells[index].fun1 = getKin(val);
            cells[index].fun2 = getLap(val);
            cells[index].dx   = 1.5 * cube.hvec[0]*(val[1] / gnorm);
            cells[index].dy   = 1.5 * cube.hvec[1]*(val[2] / gnorm);
            cells[index].dz   = 1.5 * cube.hvec[2]*(val[3] / gnorm);
            cells[index].norm = gnorm;
            cells[index].attr = SET_VOID;
            cells[index].max  = max;

            qcharge += val[0];

            if( max == 26){
                natt++;
            }

        }
    }
//}//end omp
    
    printf("Qcharge % 10.6lf\n",qcharge);

    mergeSort(cells,0,npt);

    getHash(cells,hash,npt);
    checkHash(npt,cells,hash);

    assignAttr(npt,cells,cube,hash,param,matU);

    double v0;
    double vol[10],q[10],lap[10];
    double vecA[3],vecB[3],vecC[3];
    double tmp[3];
    vecA[0] = cube.mvec[0]; vecA[1] = cube.mvec[1]; vecA[2] = cube.mvec[2];
    vecB[0] = cube.mvec[3]; vecB[1] = cube.mvec[4]; vecB[2] = cube.mvec[5];
    vecC[0] = cube.mvec[6]; vecC[1] = cube.mvec[7]; vecC[2] = cube.mvec[8];
    
    crossProduct(vecA,vecB,tmp);
    v0 = fabs( dotProduct(vecC,tmp));

    for(i=0;i<10;i++){
        vol[i] = (double) 0.;
        lap[i] = (double) 0.;
        q[i]   = (double) 0.;
    }
    double totalvol = (double) 0.;
    double totalq   = (double) 0.;
    double totallap = (double) 0.;
    printf("v0 = %lf\n",v0);

    for(index = 0; index < npt; index++){
        if( cells[index].attr >= 0 ){
            q[cells[index].attr] += v0 * cells[index].fun0;
          vol[cells[index].attr] += v0;
          lap[cells[index].attr] += v0 * cells[index].fun2;
        }
    }

    for(i=0; i< 3; i++){
        printf(" Attr %2d  Vol = % 10.6E e = % 10.6lf q = % 10.6lf Lap = % 10.6lf\n",i+1,
        vol[i],q[i],cube.zatm[i] - q[i] ,lap[i]);
        totalvol += vol[i];
        totallap += lap[i];
        totalq   += q[i];
    }
    printf("-------------------------------------------------------------\n");
    printf("   TOTAL         % 10.6lf     % 10.6lf        % 10.6lf\n",totalvol,totalq,totallap);
    printf("  Q  TOTAL       % 10.6lf     % 10.6lf        % 10.6lf\n",qcharge,v0,qcharge*v0);
    printf("  Q0 TOTAL       % 10.6lf     % 10.6lf        % 10.6lf\n",qcharge0,v0,qcharge0*v0);


    for(i=0;i<10;i++){
        vol[i] = (double) 0.;
        lap[i] = (double) 0.;
        q[i]   = (double) 0.;
    }
    totalvol = (double) 0.;
    totalq   = (double) 0.;
    totallap = (double) 0.;
    printf("v0 = %lf\n",v0);

#pragma omp parallel for shared (cells) private(index)reduction(+: q,vol,lap) 
    for(index = 0; index < npt; index++){
        if( cells[index].attr >= 0 ){
          q[cells[index].attr]   += v0 * cells[index].fun0;
          vol[cells[index].attr] += v0;
          lap[cells[index].attr] += v0 * cells[index].fun2;
        }
    }

    for(i=0; i< 10; i++){
        printf(" Attr %2d  Vol = % 10.6lf Q = % 10.6lf  Lap = % 10.6lf\n",i+1,
        vol[i],q[i],lap[i]);
        totalvol += vol[i];
        totallap += lap[i];
        totalq   += q[i];
    }
    printf("-------------------------------------------------------------\n");
    printf("   TOTAL         % 10.6lf     % 10.6lf        % 10.6lf\n",totalvol,totalq,totallap);

    
    free(hash);
    free(cells);
}


int loadFieldBasins(int i, int j, int k, int n1, int n2,
                     double *hvec, double *field, double *val){
    int p,q,r,mu;
    int sum;
    double f[27];

    mu = 0;
    for( p= i - 1; p <= i + 1; p++)
        for( q= j - 1; q <= j + 1; q++)
            for( r= k - 1; r <= k + 1; r++){
                f[mu] = field[p*n1 + q*n2 + r];
                mu++;
            }

    hessPol02(hvec[0],hvec[1],hvec[2],f,val);
    sum = 0;
    for( mu = 0; mu < 27; mu++){
        if(f[mu] < f[13] )
            sum  += 1;
    }

    return sum;
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
void assignAttr(int npt, dataCells *cells,dataCube cube, int* hash, dataRun param, double *matU){
    double r[3];
    int  i,j,k;
    int idx;
    int n1 = cube.pts[1]*cube.pts[2];
    int n2 = cube.pts[2];

    double coor[10][3];

    idx = cells[0].idx;
    printf(" cells[0].idx = %3d\n",cells[0].idx);
    printf(" hash[cells[0].idx] = %3d\n",hash[cells[0].idx]);


    int new_idx;

    int nnuc=0;
    int ii,jj,kk;
    for(idx=0;idx<npt;idx++){
        if( cells[idx].max == 26 ){
            coor[nnuc][0] = cells[idx].r[0];
            coor[nnuc][1] = cells[idx].r[1];
            coor[nnuc][2] = cells[idx].r[2];
            cells[idx].attr = nnuc;
            nnuc++;
        }
    }

    printf(" Nucleos totales %d\n",nnuc);
    printf(" Coordenadas \n");
    for(i=0;i<nnuc;i++)
        printf("%4d % 10.6lf % 10.6lf % 10.6lf\n",i,coor[i][0],coor[i][1],coor[i][2]);


//#pragma omp parallel private(idxGlobal,idx,nwi,nwj,nwk,r,dij)\
                     shared(npt,cells,n1,n2,h,nnuc,coor,r0,cube,have,hash)
//{
//#pragma omp single
    printf(" Number of threads to assigned attractors   : %4d\n",omp_get_num_threads());
//#pragma omp barrier

//#pragma omp for schedule (dynamic)
    for(idx = 0; idx < npt; idx++){
        if( cells[idx].attr == -5  && cells[idx].fun0 > TOLERANCE){
            r[0] = cells[idx].r[0] + cells[idx].dx;
            r[1] = cells[idx].r[1] + cells[idx].dy;
            r[2] = cells[idx].r[2] + cells[idx].dz;

            new_idx = getIdx(cube,r);

//            decomposeIdx(cells[idx].idx,n1,n2,&i,&j,&k);
//            decomposeIdx(new_idx,n1,n2,&ii,&jj,&kk);

            /*printf("(%3d,%3d,%3d)rho(%10.6E) ---> (%3d,%3d,%3d)rho(%10.6E)  ",
                    i,j,k,cells[idx].fun0,ii,jj,kk,cells[hash[new_idx]].fun0);

            printf("(%4d) <- (%4d) ",cells[idx].attr,cells[hash[new_idx]].attr);

            if( cells[idx].fun0 > cells[hash[new_idx]].fun0 )
                printf(" We need check it! % 10.6lf % 10.6lf % 10.6lf  % 10.6lf %4d\n",cells[idx].r[0]*B2A,cells[idx].r[1]*B2A,cells[idx].r[2]*B2A, cells[idx].norm,cells[idx].max);
            else
                printf("\n");
            */
            if( cells[idx].fun0 > cells[hash[new_idx]].fun0){
                new_idx = ascendingGLine( idx, hash, r, cells, cube, param, matU);
            }else{
                if( cells[hash[new_idx]].attr == -5)
                    new_idx = ascendingGLine( idx, hash, r, cells, cube, param, matU);

                cells[idx].attr = cells[hash[new_idx]].attr;
            }



          }

    }
//}
   
   int nbas=0;
    for(i=0;i<npt;i++){
        if(cells[i].attr >= 0)
        nbas++;
    }
    
    FILE *file;
    openFile(&file,"basins.xyz","w+");
    fprintf(file," %5d\n Basins \n",nbas);
    for(i=0;i<npt;i++){
        if(cells[i].attr >= 0)
            fprintf(file," BAS%d % 10.6lf % 10.6lf % 10.6lf\n",(cells[i].attr)+1,cells[i].r[0]*B2A,
            cells[i].r[1]*B2A, cells[i].r[2]*B2A);
    }
    fclose(file);

}

int ascendingGLine( int idx, int *hash, double r[3], dataCells *cells, 
                    dataCube cube, dataRun param, double *matU){
    int iter,flag=0;
    int new_idx;
    double val[10];
    double gnorm;

    numCritical01VecBasins(r,cube,param,matU,val);
    gnorm = getGrd(val);

    
    printf(" Entramos a ascendingGLine  ");
    iter = 0;
    do {
        r[0] += 0.05*cube.hvec[0] * val[1];
        r[1] += 0.05*cube.hvec[1] * val[2];
        r[2] += 0.05*cube.hvec[2] * val[3];
        
        new_idx = getIdx(cube,r);
        if( cells[hash[new_idx]].attr >= 0 ){
            flag = 1;
        }
        iter++;

    } while(iter < 2000 && flag == 0);
    printf(" [%5d] \n",iter);

    return new_idx;

}

int getAttr(int *nnuc, double coor[10][3],double r0[3],double h[3], int i, int j , int k,  int n1,int n2,dataCube  cube){

    int p,q,r;
    int mu;
    double f[27];
    mu = 0;
    int sum;
    for(p=i-1;p<= i+1; p++){
        for(q=j-1;q<= j+1; q++){
            for(r=k-1;r<= k+1; r++){
                f[mu] = cube.field[p*n1+q*n2+r];
                mu++;
            }
        }
    }


    for( mu = 0; mu < 27; mu++){
        if(f[mu] < f[13] )
            sum  += 1;
    }

    return sum;
}
