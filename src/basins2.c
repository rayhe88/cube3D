#include "basins2.h"
#include "fields.h"
#include "utils.h"
#include "findCrit.h"
#include "mergeSort.h"
#include "kernels.h"
#include "file.h"
#include "transU.h"
#include "numBondPath.h"
#include "jacobi.h"
#include "tableP.h"
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
        (*ptr)[i].max  = -1;
        (*ptr)[i].idx  = -1;
        (*ptr)[i].fun0 = (double)0.;
        (*ptr)[i].fun1 = (double)0.;
        (*ptr)[i].fun2 = (double)0.;
    }

}

void evalBasins    (dataCube cube, dataRun param, double *matU, double min0, char *name){

    printf("%s\n We enter in \'evalBains\' \n",TGBI);

    double matT[9];
    getMatInv(matU,matT);
    double r0[3],q0[3];
    int i, zsum=0;
    double max[3],min[3];

    for(i=0;i<cube.natm;i++){
        zsum += cube.zatm[i];
    }
    printf(" The system has %5d protons\n",zsum);
    
    if( param.orth != YES){
        r0[0] = cube.min[0];
        r0[1] = cube.min[1]; 
        r0[2] = cube.min[2]; 

        getRiU(r0,matT,q0);

        cube.min[0] = q0[0];
        cube.min[1] = q0[1];
        cube.min[2] = q0[2];
        
        for(i=0;i<cube.natm;i++){
            r0[0] = cube.coor[3*i];
            r0[1] = cube.coor[3*i+1];
            r0[2] = cube.coor[3*i+2];
            
            getRiU(r0,matT,q0);

            cube.coor[3*i]   = q0[0];
            cube.coor[3*i+1] = q0[1];
            cube.coor[3*i+2] = q0[2];
        }
    }

    if( param.pbc == YES ){

        min[0] = cube.min[0];
        min[1] = cube.min[1];
        min[2] = cube.min[2];

        max[0] = min[0] + (cube.pts[0] - 1)*cube.hvec[0];
        max[1] = min[1] + (cube.pts[1] - 1)*cube.hvec[1];
        max[2] = min[2] + (cube.pts[2] - 1)*cube.hvec[2];
        
        printf(" X RANGE : % 10.6lf  % 10.6lf\n",min[0] * B2A, max[0] * B2A);
        printf(" Y RANGE : % 10.6lf  % 10.6lf\n",min[1] * B2A, max[1] * B2A);
        printf(" Z RANGE : % 10.6lf  % 10.6lf\n",min[2] * B2A, max[2] * B2A);
        
        printf(" MAX FROM CUBE \n");
        printf(" X RANGE : % 10.6lf\n", cube.max[0] * B2A );
        printf(" Y RANGE : % 10.6lf\n", cube.max[1] * B2A );
        printf(" Z RANGE : % 10.6lf\n", cube.max[2] * B2A );
        cube.max[0] = max[0];
        cube.max[1] = max[1];
        cube.max[2] = max[2];

        for(i=0;i<cube.natm;i++){
            r0[0] = cube.coor[3*i];
            r0[1] = cube.coor[3*i+1];
            r0[2] = cube.coor[3*i+2];
            
            perfectCube(param.pbc, r0, min,max);

            cube.coor[3*i]   = r0[0];
            cube.coor[3*i+1] = r0[1];
            cube.coor[3*i+2] = r0[2];
        }

    }

    
    getCells (cube,param,min0,matU,name);

    printf(" We out of \'evalBasins\'\n%s\n",TRST);

}

void getCells (dataCube cube, dataRun param, double min0, double *matU, char *name){
    
    int i,j;
    int npt = cube.npt;
    int natt, index;
    int maxatt;
    int *hash;
    double qcharge;
    double qcharge0;
    double *field2;
    double totalvol = (double) 0.;
    double totalq   = (double) 0.;
    double totallap = (double) 0.;
    double *vol,*lap,*ne,*q;
    double *coorAt;
    double v0 = getVolumenCell(cube.mvec);
    double hmax = getMax3(cube.hvec);
    dataCells2 *cells;

    FILE *file;
    char nameOut[128];
    sprintf(nameOut,"%sBasins.log",name);
#ifdef DEBUG
    file = stdout;
#else
    openFile(&file,nameOut,"w+");
#endif

    createArrayCells (npt, &cells , "Basins cells");
    createArrayDou   (npt, &field2, "Basins cells");
    createArrayInt   (npt, &hash  , "Basins cells");

    qcharge0 = (double) 0.;
    for(index = 0; index < npt; index++){
        field2[index] = exp(cube.field[index]) + min0 - DELTA;
        qcharge0 += field2[index];
    }

    if( param.pbc == YES){
        natt = kernelBasinsLoadPer   (cells,cube,param,field2,matU,&qcharge);
    }else{
        natt = kernelBasinsLoadNoPer (cells,cube,param,field2,matU,&qcharge);
    }

    /**********************************************************************/
    maxatt = natt;
    printf(" %5d local maximums were found\n",natt);
    if( natt > cube.natm){
        fprintf(file," More local maximums were found than nuclei in the system,");
        fprintf(file,"\n possible there are/is [%5d] NNACPs.\n",natt-cube.natm);
        maxatt = natt;
    }else{
        if( natt < cube.natm){
            fprintf(file," Fewer local maximums were found than nuclei in the system,");
            fprintf(file," possibly having a poor cube.\n");
            maxatt = cube.natm;
        }
    }
    createArrayDou ( 3*maxatt,&coorAt, "Basins cells");

    //  this works when natt = cube.natm  or natt < cube.natm
    double ratm[3];
    double dij;
    double mind;
    double matT[9];
    getMatInv(matU,matT);

    for(i=0;i<cube.natm;i++){
        ratm[0] = cube.coor[3*i];
        ratm[1] = cube.coor[3*i+1];
        ratm[2] = cube.coor[3*i+2];

        ratm[0] *= B2A;
        ratm[1] *= B2A;
        ratm[2] *= B2A;
        printf(" %d  % 10.6lf % 10.6lf % 10.6lf\n",cube.zatm[i],ratm[0],ratm[1],ratm[2]);
    }

    printf(" Cual es el valor de hmax % 10.5lf\n",hmax);
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

                        printf(" cells[%9d] has %5d as attractor\n",index,cells[index].attr);
                    }
                }
                i++;
            }
        }
    }
    for(index=0; index<npt; index++){
        if( cells[index].max  == 26 && cells[index].attr >= 0){
            j = cells[index].attr;
            coorAt[3*j  ] = cells[index].r[0];
            coorAt[3*j+1] = cells[index].r[1];
            coorAt[3*j+2] = cells[index].r[2];
        }
    }

    j = cube.natm;
    if( natt > cube.natm ){
        for(index=0; index<npt; index++){
            if( cells[index].max  == 26 && cells[index].attr < 0){
                coorAt[3*j  ] = cells[index].r[0];
                coorAt[3*j+1] = cells[index].r[1];
                coorAt[3*j+2] = cells[index].r[2];
                cells[index].attr = j;
                j++;
            }
        }
    }
//#ifndef DEBUG
    for(index=0; index<npt; index++){
        if( cells[index].attr >= 0){
            printf(" Cells[%9d]  with Atractor %4d  and  % 8.3lf % 8.3lf % 8.3lf\n",
            index,cells[index].attr,cells[index].r[0],cells[index].r[1],cells[index].r[2]);
        }

    }

//#endif
    /**********************************************************************/
    fprintf(file," Qcharge     : % 12.6lf\n",qcharge*v0);
    fprintf(file," Qcharge0    : % 12.6lf\n",qcharge0*v0);
    fprintf(file," Attractors  : % 12d\n",natt);
    fprintf(file," volumen cell: % 12.6lf\n",v0);
    /**********************************************************************/

    for(index=0;index<natt;index++)
        fprintf(file," %5d  % 10.6lf % 10.6lf % 10.6lf\n",index,
                coorAt[3*index]*B2A,coorAt[3*index+1]*B2A,coorAt[3*index+2]*B2A);

    /**********************************************************************/
    printf(" 0\n");
    mergeSort(cells,0,npt);

    printf(" 1\n");

    getHash  (cells,hash,npt);

    printf(" 2\n");
    natt = assignAttr(npt,natt,cells,cube,param,min0,matU,coorAt);
    /**********************************************************************/


    printf(" 3\n");
    createArrayDou   (natt, &q  , "Basins cells");
    createArrayDou   (natt, &ne , "Basins cells");
    createArrayDou   (natt, &lap, "Basins cells");
    createArrayDou   (natt, &vol, "Basins cells");
    printf(" 4\n");

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
        fprintf(file,"  Attr %3d  %9.3E % 10.5lf  % 10.5lf  % 9.4E\n",i+1,
        vol[i],q[i],ne[i],lap[i]);
        totalvol += vol[i];
        totallap += lap[i];
        totalq   += q[i];
        totalcha += ne[i];
    }
    fprintf(file,"-------------------------------------------------------------\n");
    fprintf(file,"   TOTAL    %9.3E % 10.5lf  % 10.5lf  % 9.4E\n",totalvol,totalq,totalcha,totallap);
    fprintf(file,"-------------------------------------------------------------\n");
    fprintf(file,"  Q  TOTAL   % 10.6lf \n",qcharge*v0);
    fprintf(file,"-------------------------------------------------------------\n");
    printBar(stdout);
    printf("  File %s was generated\n",nameOut);


    printXYZBasins( npt,cells,cube,hash,matU,name);

    free(q);
    free(ne);
    free(lap);
    free(vol);
    free(hash);
    free(cells);
    free(coorAt);
    free(field2);
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

int getPerIndex(int i, int nx){
    int ret;
    ret = i;
    if( i < 0 )  ret += (nx-1);
    if( i >= nx) ret -= (nx-1);

    return ret;
}
int loadFieldBasinsPer( int i, int j, int k, int n1, int n2,
                        double *hvec, double *field, double *val,int *np){
    int p,q,r,mu;
    int sum;
    int nwp,nwq,nwr;
    double f[27];

    mu = 0;
    for( p= i - 1; p <= i + 1; p++)
        for( q= j - 1; q <= j + 1; q++)
            for( r= k - 1; r <= k + 1; r++){
                nwp = getPerIndex(p,np[0]);
                nwq = getPerIndex(q,np[1]);
                nwr = getPerIndex(r,np[2]);
                
                f[mu] = field[nwp*n1 + nwq*n2 + nwr];
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

int assignAttr(int npt, int natt, dataCells2 *cells,dataCube cube, dataRun param, double min0, double *matU, double *coorAt){
    int index;
    double hav;
    double rout[3];
    double matT[9];

    getMatInv(matU,matT);

    hav  = cube.hvec[0] * cube.hvec[0];
    hav += cube.hvec[1] * cube.hvec[1];
    hav += cube.hvec[2] * cube.hvec[2];
    hav = sqrt(hav);

#pragma omp parallel private(index,rout)                                   \
                     shared(cube, cells, natt, hav, coorAt, param, min0, matU,matT)
{
#pragma omp single
    printf(" Number of threads to assign attractors to cells : %4d\n",omp_get_num_threads());
#pragma omp barrier

#pragma omp for schedule (dynamic)
    for(index = 0; index < cube.npt; index++){
        if (cells[index].attr == SET_VOID && cells[index].fun0 > TOLERANCE){
            cells[index]. attr = ascendingGLine( natt, cells[index].r,
                                                 hav, coorAt, cube, param, 
                                                 min0, matU, matT, 3000,rout);

            if( cells[index].attr == -1)
                cells[index]. attr = ascendingGLine2( natt, cells[index].r,
                                                     hav, coorAt, cube, param, 
                                                     min0, matU, matT, 3000,rout);

        }

    }
}


    


    
    return natt;
}

int ascendingGLine( int nattr, double qc[3], double haverage, double *coorAt,
                    dataCube cube, dataRun param, double min0,double *matU,
                    double *matT, int nmax, double qout[3]){
   
    int atractor = -1;
    int i,iter,flag;
    double qi[3],qn[3];
    double ri[3],rn[3];
    double val[10],vec[3];
    double gnorm;
    double dist;

    qi[0] = qc[0];  qi[1] = qc[1]; qi[2] = qc[2];
    
    iter = 0; flag = 0;

    while(iter < nmax && flag == 0){
        getRiU(qi,matU,ri);
        numCritical01Vec(qi,cube,param,matU,min0,val);
        gnorm = getGrd(val);
        vec[0] = val[1]/gnorm; vec[1] = val[2]/gnorm; vec[2] = val[3]/gnorm;
        if( myIsNanInf_V3(vec)!= 0){
            vec[0] = vec[1] = vec[2] = 0.577350269;
        }
        rn[0] = ri[0] + 0.5 * cube.hvec[0] * vec[0];
        rn[1] = ri[1] + 0.5 * cube.hvec[1] * vec[1];
        rn[2] = ri[2] + 0.5 * cube.hvec[2] * vec[2];

        getRiU(rn,matT,qn);
        perfectCube(param.pbc,qn,cube.min,cube.max);
        cpyVec3(qn,qi);
        

        for(i=0;i<nattr;i++){
            dist = distance2(coorAt[3*i],coorAt[3*i+1],coorAt[3*i+2],qi);
            if ( dist < 1.* haverage){
                atractor = i;
                flag = 1;
                break;
            }
        }
        iter++;
    }
    qout[0] = qn[0]; qout[1] = qn[1]; qout[2] = qn[2];


    return atractor;
}
int ascendingGLine2( int nattr, double qc[3], double haverage, double *coorAt,
                    dataCube cube, dataRun param, double min0,double *matU,
                    double *matT, int nmax, double qout[3]){
   
    int atractor = -1;
    int i,iter,flag;
    double qi[3],qn[3];
    double ri[3],rn[3];
    double val[10];
    double matH[9];
    double eval[3],evec[9];
    double vec2[3];
    double gnorm;
    double dist;

    qi[0] = qc[0];  qi[1] = qc[1]; qi[2] = qc[2];
    
    iter = 0; flag = 0;
    
    numCritical02Vec(qi,cube,param,matU,min0,val);
    matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
    matH[3] = val[7]; matH[4] = val[5]; matH[5] = val[9];
    matH[6] = val[8]; matH[7] = val[9]; matH[8] = val[6];
    JacobiNxN(matH,eval,evec);
    vec2[0] =evec[6];  vec2[1] =evec[7]; vec2[2] =evec[8]; 
    gnorm = getNormVec(vec2);

    qi[0] = qc[0] + 0.02 * vec2[0];
    qi[1] = qc[1] + 0.02 * vec2[1];
    qi[2] = qc[2] + 0.02 * vec2[2];

    while(iter < nmax && flag == 0){
        getRiU(qi,matU,ri);
        numCritical01Vec(qi,cube,param,matU,min0,val);
        gnorm = getGrd(val);

        rn[0] = ri[0] + 0.5 * cube.hvec[0] * val[1]/gnorm;
        rn[1] = ri[1] + 0.5 * cube.hvec[1] * val[2]/gnorm;
        rn[2] = ri[2] + 0.5 * cube.hvec[2] * val[3]/gnorm;

        getRiU(rn,matT,qn);
        perfectCube(param.pbc,qn,cube.min,cube.max);
        cpyVec3(qn,qi);
        

        for(i=0;i<nattr;i++){
            dist = distance2(coorAt[3*i],coorAt[3*i+1],coorAt[3*i+2],qi);
            if ( dist < 1.* haverage){
                atractor = i;
                flag = 1;
                break;
            }
        }
        iter++;
    }
    qout[0] = qn[0]; qout[1] = qn[1]; qout[2] = qn[2];


    return atractor;
}

void printXYZBasins( int n, dataCells2 *cells, dataCube cube, int *hash, double *matU, char *name){
    char nameOut[128];
    char symbol[6];
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
    double r[3];

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
        if( cells[idx].attr >= -1 && cells[idx].max == 1 ){
            new++;
        }
    }

    mergeSortbyAtr(cells,0,n);

    fprintf(file," %6d\n",new);
    fprintf(file," prueba\n");
    for(idx=0;idx<n;idx++){
        if( cells[idx].attr >= -1  && cells[idx].max == 1){
            getRiU(cells[idx].r,matU,r);
            r[0] *= B2A;
            r[1] *= B2A;
            r[2] *= B2A;
            
            if( cells[idx].attr >= 0 && cells[idx].attr < cube.natm)
                getAtomicSymbol(cube.zatm[cells[idx].attr],4,symbol);
            else
                strcpy(symbol,"NNA");

            fprintf(file,"%s-%04d % 10.6lf % 10.6lf % 10.6lf\n",symbol,
                            cells[idx].attr + 1,r[0],r[1],r[2]);
        }
    }

    fclose(file);

    printBar(stdout);
    printf("  File %s was generated\n",nameOut);
}

int kernelBasinsLoadNoPer(dataCells2 *cells,dataCube cube, dataRun param, double *field2,double *matU, double *charge){
    int index,i,j,k,max;
    double val[10];

    int natt = 0;
    double qcharge  = (double) 0.;
    int n1 = cube.pts[1] * cube.pts[2];
    int n2 = cube.pts[2];
#pragma omp parallel private(index,i,j,k,val,max)                                   \
                     shared(qcharge,field2,n1,n2,cells,cube,matU,natt)
{
#pragma omp single
    printf(" Number of threads for evaluation of fields : %4d\n",omp_get_num_threads());
#pragma omp barrier

#pragma omp for schedule (dynamic)
    for(index = 0; index < cube.npt; index++){

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

    (*charge) = qcharge;
    
    return natt;

}
int kernelBasinsLoadPer(dataCells2 *cells,dataCube cube, dataRun param, double *field2,double *matU, double *charge){
    int index,i,j,k,max;
    double val[10];

    int natt = 0;
    double qcharge  = (double) 0.;
    int n1 = cube.pts[1] * cube.pts[2];
    int n2 = cube.pts[2];
#pragma omp parallel private(index,i,j,k,val,max)                                   \
                     shared(qcharge,field2,n1,n2,cells,cube,matU,natt)
{
#pragma omp single
    printf(" Number of threads for evaluation of fields : %4d\n",omp_get_num_threads());
#pragma omp barrier

#pragma omp for schedule (dynamic)
    for(index = 0; index < cube.npt; index++){
        
        decomposeIdx(index,n1,n2,&i,&j,&k);

        max = loadFieldBasinsPer(i,j,k,n1,n2,cube.hvec,field2,val,cube.pts);
        
        if( param.orth != YES) trans02(val,matU);

        cells[index].idx = index;
        cells[index].r[0] = cube.min[0] + i*cube.hvec[0];
        cells[index].r[1] = cube.min[1] + j*cube.hvec[1];
        cells[index].r[2] = cube.min[2] + k*cube.hvec[2];
        cells[index].max  = max;
        cells[index].attr = SET_VOID;

        if( i < cube.pts[0] - 1 &&
            j < cube.pts[1] - 1 &&
            k < cube.pts[2] - 1){
            cells[index].fun0 = val[0];
            cells[index].fun1 = getKin(val);
            cells[index].fun2 = getLap(val);

        }

#pragma omp critical
{
            qcharge += cells[index].fun0;
            if( cells[index].max == 26){
                natt++;
            }
}

    }
}//end omp

    (*charge) = qcharge;
    
    return natt;

}

int ascendingGLinePrint( int nattr, double qc[3], double haverage, double *coorAt,
                    dataCube cube, dataRun param, double min0,double *matU,
                    double *matT, int nmax, double qout[3]){

    static int flag1 = 0;
    if( flag1 > 3 ) {
        return ascendingGLine(nattr,qc, haverage, coorAt, cube, param, min0, matU,
                       matT, nmax, qout);
        
    }
   
    int atractor = -1;
    int i,iter,flag;
    double qi[3],qn[3];
    double ri[3],rn[3];
    double val[10];
    double gnorm;
    double dist;

    qi[0] = qc[0];  qi[1] = qc[1]; qi[2] = qc[2];
    
    iter = 0; flag = 0;

    while(iter < nmax && flag == 0){
        getRiU(qi,matU,ri);
        numCritical01Vec(qi,cube,param,matU,min0,val);
        gnorm = getGrd(val);

        rn[0] = ri[0] + 0.5 * cube.hvec[0] * val[1]/gnorm;
        rn[1] = ri[1] + 0.5 * cube.hvec[1] * val[2]/gnorm;
        rn[2] = ri[2] + 0.5 * cube.hvec[2] * val[3]/gnorm;

        getRiU(rn,matT,qn);
        perfectCube(param.pbc,qn,cube.min,cube.max);
        cpyVec3(qn,qi);
        

        for(i=0;i<nattr;i++){
            dist = distance2(coorAt[3*i],coorAt[3*i+1],coorAt[3*i+2],qi);
            printf("                 dij = % 10.6E  , atm %4d\n",dist,i+1);
            if ( dist < 3.* haverage){
                atractor = i;
                flag = 1;
                break;
            }
        }
        iter++;
    }
    qout[0] = qn[0]; qout[1] = qn[1]; qout[2] = qn[2];


    flag1++;
    return atractor;
}
