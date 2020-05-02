#include "basins2.h"
#include "fields.h"
#include "utils.h"
#include "findCrit.h"
#include "mergeSort.h"
#include "kernels.h"
#include "file.h"
#include <omp.h>

#define FALSE 0
#define TRUE !FALSE

double getMax3( double val[3]){
    
    if(val[1] > val[2]){
        if( val[1] > val[3])
            return val[1];
        else
            return val[3];
    }else{
        if(val[2] > val[3])
            return val[2];
        else
            return val[3];
    }
}

void getHash(dataCells2 *cells, int *hash, int npt){
    int i;
    for(i=0;i<npt;i++)
        hash[i] = 0;

    for(i=0;i<npt;i++){
        hash[cells[i].idx] = i;
    }

}

double getVolumenCell( double mvec[9] ){
    double vecA[3],vecB[3],vecC[3];
    double tmp[3];

    vecA[0] = mvec[0]; vecA[1] = mvec[1]; vecA[2] = mvec[2];
    vecB[0] = mvec[3]; vecB[1] = mvec[4]; vecB[2] = mvec[5];
    vecC[0] = mvec[6]; vecC[1] = mvec[7]; vecC[2] = mvec[8];
    
    crossProduct(vecA,vecB,tmp);
    return fabs( dotProduct(vecC,tmp));
}

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

void createArrayCells(int n, dataCells2 **ptr, const char *mess){
    int i;
    if( n < 1){
        printf(" There is an error in the size for [%s]\n",mess);
        exit(EXIT_FAILURE);
    }

    *ptr = (dataCells2*) malloc (n * sizeof(dataCells2));
    if( *ptr == NULL){
        printf(" Failed to allocate memory: [%s]\n",mess);
        exit(EXIT_FAILURE);
    }

    for(i=0; i<n; i++){
        (*ptr)[i].attr = SET_NULL;
        (*ptr)[i].max  = FALSE;
        (*ptr)[i].idx  = -1;
        (*ptr)[i].fun0 = (double)0.;
        (*ptr)[i].fun1 = (double)0.;
        (*ptr)[i].fun2 = (double)0.;
    }

}

void evalBasins    (dataCube cube, dataRun param, double *matU, double min0, char *name){

    printf("%s\n We enter in \'evalBains\' \n",TGBI);

    if ( param.pbc == YES )
        getCellsPer   (cube,param,min0,matU,name);
    else
        getCellsNoPer (cube,param,min0,matU,name);

    printf(" We out of \'evalBasins\'\n%s\n",TRST);

}


void getCellsPer   (dataCube cube, dataRun param, double min0, double *matU, char *name){

    printf(" Periodic subroutine\n");
    

}

void getCellsNoPer (dataCube cube, dataRun param, double min0, double *matU, char *name){
    
    int i,j,k;
    int npt = cube.npt;
    int n1 = cube.pts[1]*cube.pts[2];
    int n2 = cube.pts[2];
    int natt, max, index;
    int maxatt;
    int *hash;
    double val[10];
    double qcharge;
    double qcharge0;
    double *field2;
    double totalvol = (double) 0.;
    double totalq   = (double) 0.;
    double totallap = (double) 0.;
    double *vol,*lap,*ne,*q;
    double **coorAt;
    double v0 = getVolumenCell(cube.mvec);
    double hmax = getMax3(cube.hvec);
    dataCells2 *cells;

    FILE *file;
    char nameOut[128];
    sprintf(nameOut,"%sBasins.log",name);
    openFile(&file,nameOut,"w+");

    createArrayCells (npt, &cells , "Basins cells");
    createArrayDou   (npt, &field2, "Basins cells");
    createArrayInt   (npt, &hash  , "Basins cells");

    qcharge0 = (double) 0.;
    for(index = 0; index < npt; index++){
        field2[index] = exp(cube.field[index]) + min0 - DELTA;
        qcharge0 += field2[index];
    }


    natt = 0;
    qcharge = (double) 0.;
#pragma omp parallel private(index,i,j,k,val,max)                                   \
                     shared(qcharge,qcharge0,field2,npt,n1,n2,cells,cube,matU,natt)
{
#pragma omp single
    printf(" Number of threads for evaluation of fields : %4d\n",omp_get_num_threads());
#pragma omp barrier

#pragma omp for schedule (dynamic)
    for(index = 0; index < npt; index++){
        

        decomposeIdx(index,n1,n2,&i,&j,&k);

        cells[index].idx = index;
        cells[index].r[0] = cube.min[0] + i*cube.hvec[0];
        cells[index].r[1] = cube.min[1] + j*cube.hvec[1];
        cells[index].r[2] = cube.min[2] + k*cube.hvec[2];
        cells[index].attr = SET_NULL;

        if( ( i > 0 && i < cube.pts[0]-1) && 
            ( j > 0 && j < cube.pts[1]-1) && 
            ( k > 0 && k < cube.pts[2]-1) ){ 
            
            max = loadFieldBasins(i,j,k,n1,n2,cube.hvec,field2,val);

            if( param.orth != YES) trans02(val,matU);

            cells[index].fun0 = val[0];
            cells[index].fun1 = getKin(val);
            cells[index].fun2 = getLap(val);
            cells[index].attr = SET_VOID;
            cells[index].max  = max;
#pragma omp critical
{
            qcharge += cells[index].fun0;
            if( max == 26){
                natt++;
            }
}

        }
    }
}//end omp

    /**********************************************************************/
    maxatt = natt;
    if( natt > cube.natm){
        fprintf(file," More local maximums were found than nuclei in the system,");
        fprintf(file,"\n possible there are/is [%5d] NNACPs.\n",natt-cube.natm);
        maxatt = natt;
    }else{
        if( natt < cube.natm){
            fprintf(file," Fewer local maximums were found than nuclei in the system,");
            fprintf(file," possibly having a poor cube.");
            maxatt = cube.natm;
        }
    }
    /////////////   memory block for a matrix /////////////////////////////
    coorAt = (double**) malloc(maxatt * sizeof(double*));
    for(index=0; index<maxatt; index++)
        coorAt[index] = (double*) malloc(3 * sizeof(double));
    ///////////////////////////////////////////////////////////////////////


    //  this works when natt = cube.natm  or natt < cube.natm
    double ratm[3];
    double dij;
    double mind;
    for(index=0; index<npt; index++){
        if( cells[index].max  == 26 ){
            mind = 1.E9;
            i = 0;
            while(i < cube.natm){
                ratm[0] = cube.coor[3*i];
                ratm[1] = cube.coor[3*i+1];
                ratm[2] = cube.coor[3*i+2];
                dij = distance(ratm,cells[index].r);
                if( dij < mind ){
                    mind = dij;
                    if( dij < 3*hmax){
                        cells[index].attr = i;
                        i = cube.natm;
                    }
                }
                i++;
            }
        }
    }
    for(i=0;i<cube.natm;i++){
        coorAt[i][0] = cube.coor[3*i];
        coorAt[i][1] = cube.coor[3*i+1];
        coorAt[i][2] = cube.coor[3*i+2];
    }
    //  this works when natt > cube.natm
    j = cube.natm;
    if( natt > cube.natm ){
        for(index=0; index<npt; index++){
            if( cells[index].max  == 26 && cells[index].attr < 0){
                coorAt[j][0] = cells[index].r[0];
                coorAt[j][1] = cells[index].r[1];
                coorAt[j][2] = cells[index].r[2];
                cells[index].attr = j;
                j++;
            }
        }
    }
#ifdef DEBUG
    for(index=0; index<npt; index++){
        if( cells[index].attr >= 0){
            printf(" Cells[%9d]  with Atractor %4d  and  % 8.3lf % 8.3lf % 8.3lf\n",
            index,cells[index].attr,cells[index].r[0],cells[index].r[1],cells[index].r[2]);
        }

    }

#endif
    /**********************************************************************/
    fprintf(file," Qcharge     : % 10.6lf\n",qcharge*v0);
    fprintf(file," Qcharge0    : % 10.6lf\n",qcharge0*v0);
    fprintf(file," Attractors  :    % 10d\n",natt);
    fprintf(file," volumen cell: % 10.6lf\n",v0);
    /**********************************************************************/

    for(index=0;index<natt;index++)
        fprintf(file," %5d  % 10.6lf % 10.6lf % 10.6lf\n",index,coorAt[index][0]*B2A,coorAt[index][1]*B2A,coorAt[index][2]*B2A);

    /**********************************************************************/
    mergeSort(cells,0,npt);

    getHash  (cells,hash,npt);

    natt = assignAttr(npt,natt,cells,cube,param,min0,matU,coorAt);
    /**********************************************************************/


    createArrayDou   (natt, &q  , "Basins cells");
    createArrayDou   (natt, &ne , "Basins cells");
    createArrayDou   (natt, &lap, "Basins cells");
    createArrayDou   (natt, &vol, "Basins cells");

    for(i=0;i<natt;i++){
        q[i]   = (double) 0.;
        ne[i]  = (double) 0.;
        lap[i] = (double) 0.;
        vol[i] = (double) 0.;
    }
    totalvol = (double) 0.;
    totalq   = (double) 0.;
    totallap = (double) 0.;
    double totalcha = (double) 0.;

    for(index = 0; index < npt; index++){
        if( cells[index].attr >= 0 ){
          q[cells[index].attr]   += cells[index].fun0;
          vol[cells[index].attr] += 1.;
          lap[cells[index].attr] += cells[index].fun2;
        }
    }
    for(i=0; i< natt; i++){
        vol[i] *= v0;
        lap[i] *= v0;
        q[i]   *= v0;
        if( i < cube.natm )
            ne[i] = cube.zatm[i]  - q[i];
       
    }

    fprintf(file,"-------------------------------------------------------------\n");
    fprintf(file,"  Attractor      Volumen     nElectron  Charge  Laplacian    \n");
    fprintf(file,"-------------------------------------------------------------\n");
    for(i=0; i< natt; i++){
    fprintf(file,"  Attr %3d  % 10.6E % 10.6lf % 10.6lf % 10.6E\n",i+1,
        vol[i],q[i],ne[i],lap[i]);
        totalvol += vol[i];
        totallap += lap[i];
        totalq   += q[i];
        totalcha += ne[i];
    }
    fprintf(file,"-------------------------------------------------------------\n");
    fprintf(file,"   TOTAL     % 10.6E  % 10.6lf %10.6lf % 10.6lf \n",totalvol,totalq,totalcha,totallap);
    fprintf(file,"-------------------------------------------------------------\n");
    fprintf(file,"  Q  TOTAL   % 10.6E  % 10.6lf  % 10.6lf \n",qcharge,v0,qcharge*v0);
    fprintf(file,"  Q0 TOTAL   % 10.6E  % 10.6lf  % 10.6lf \n",qcharge0,v0,qcharge0*v0);
    fprintf(file,"-------------------------------------------------------------\n");
    printBar(stdout);
    printf("  File %s was generated\n",nameOut);


    printXYZBasins( npt,cells,cube,hash,name);

    for(index=0;index<maxatt;index++)
        free(coorAt[index]);
    free(coorAt);

   
    free(q);
    free(ne);
    free(lap);
    free(vol);
    free(hash);
    free(field2);
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

int assignAttr(int npt, int natt, dataCells2 *cells,dataCube cube, dataRun param, double min0, double *matU, double **coorAt){
    int index;
    double hav;

    hav  = cube.hvec[0] * cube.hvec[0];
    hav += cube.hvec[1] * cube.hvec[1];
    hav += cube.hvec[2] * cube.hvec[2];
    hav = sqrt(hav);


#pragma omp parallel private(index)                                   \
                     shared(cube, cells, natt, hav, coorAt, param, min0, matU)
{
#pragma omp single
    printf(" Number of threads to assign attractors to cells : %4d\n",omp_get_num_threads());
#pragma omp barrier

#pragma omp for schedule (dynamic)
    for(index = 0; index < cube.npt; index++){
        if (cells[index].attr == SET_VOID && cells[index].fun0 > TOLERANCE){
            cells[index]. attr = ascendingGLine( natt, cells[index].r,
                                                 hav, coorAt, cube, param, 
                                                 min0, matU, 1000);
            if( cells[index].attr == -1 )
                cells[index]. attr = ascendingGLine( natt, cells[index].r,
                                                     hav, coorAt, cube, param, 
                                                     min0, matU, 3000);
        }

    }
}


    


    
    return natt;
}

int ascendingGLine( int nattr, double r[3], double haverage, double **coorAt,
                    dataCube cube, dataRun param, double min0,double *matU,int nmax){
   
    int atractor = -1;
    int i,iter,flag;
    double q[3];
    double val[10];
    double gnorm;
    double dist;

    q[0] = r[0];  q[1] = r[1]; q[2] = r[2];

    numCritical01Vec(q,cube,param,matU,min0,val);
    gnorm = getGrd(val);
    
    iter = 0; flag = 0;

    while(iter < nmax && flag ==0){
        q[0] += 0.5 * cube.hvec[0] * val[1]/gnorm;
        q[1] += 0.5 * cube.hvec[1] * val[2]/gnorm;
        q[2] += 0.5 * cube.hvec[2] * val[3]/gnorm;
        
        iter++;

        for(i=0;i<nattr;i++){
            dist = distance(coorAt[i],q);
            if ( dist < 3.* haverage){
                atractor = i;
                flag = 1;
                break;
            }
        }
        numCritical01Vec(q,cube,param,matU,min0,val);
        gnorm = getGrd(val);
    }


    return atractor;
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
void printXYZBasins( int n, dataCells2 *cells, dataCube cube, int *hash, char *name){
    char nameOut[128];
    int idx,new;
    int i,j,k;
    FILE *file;
    int nx = cube.pts[0];
    int ny = cube.pts[1];
    int nz = cube.pts[2];
    int n1 = ny * nz;
    int n2 = nz;
    int at0,at1,at2,at3,at4,at5,at6;
    int id01,id02,id03,id04,id05,id06;

    sprintf(nameOut,"%sBasins.xyz",name);
    openFile(&file,nameOut,"w+");

    for(idx=0;idx<n;idx++){
        cells[idx].max = 0;
        decomposeIdx(cells[idx].idx,n1,n2,&i,&j,&k);
        if( i == 0 || i == nx - 1 ) cells[idx].max = 1;
        if( j == 0 || j == ny - 1 ) cells[idx].max = 1;
        if( k == 0 || k == nz - 1 ) cells[idx].max = 1;

        if( cells[idx].attr >= -1 ){
           id01 = getCubeIdx(i-1, j, k, nx, ny, nz);
           id02 = getCubeIdx(i+1, j, k, nx, ny, nz);
           id03 = getCubeIdx(i, j-1, k, nx, ny, nz);
           id04 = getCubeIdx(i, j+1, k, nx, ny, nz);
           id05 = getCubeIdx(i, j, k-1, nx, ny, nz);
           id06 = getCubeIdx(i, j, k+1, nx, ny, nz);

           if( id01 >= 0 && id02 >= 0 && id03 >= 0 &&
               id04 >= 0 && id05 >= 0 && id06 >= 0){

               at0 = cells[idx       ].attr;
               at1 = cells[hash[id01]].attr;
               at2 = cells[hash[id02]].attr;
               at3 = cells[hash[id03]].attr;
               at4 = cells[hash[id04]].attr;
               at5 = cells[hash[id05]].attr;
               at6 = cells[hash[id06]].attr;

               if( at0 == at1 && at0 == at2 &&
                   at0 == at3 && at0 == at4 &&
                   at0 == at5 && at0 == at6){
                   cells[idx].max = 0;
               }else{
                   cells[idx].max = 1;
               }

           }
        }
    }
    new=0;
    for(idx=0;idx<n;idx++){
        if(cells[idx].attr >= -1 && cells[idx].max == 1)
            new++;
    }

    mergeSortbyAtr(cells,0,n);

    fprintf(file," %6d\n",new-1);
    fprintf(file," prueba\n");
    for(idx=0;idx<n;idx++){
        if( cells[idx].attr >= -1  && cells[idx].max == 1){
            fprintf(file,"BAS%04d % 10.6lf % 10.6lf % 10.6lf\n",
                            cells[idx].attr + 1,
                            cells[idx].r[0] * B2A,
                            cells[idx].r[1] * B2A,
                            cells[idx].r[2] * B2A);
        }
    }

    fclose(file);

    printBar(stdout);
    printf("  File %s was generated\n",nameOut);
}

