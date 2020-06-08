/** 
 * @file   numBondPath.c
 * @brief
 * @author Raymundo Hernández-Esparza.
 * @date   August 2018.
 */

#include "file.h"
#include "array.h"
#include "fields.h"
#include "jacobi.h"
#include "version.h"
#include "findCrit.h"
#include "mathTools.h"
#include "numBondPath.h"
#include "lagrange2.h"
#include "utils.h"
#include "tableP.h"
#include "transU.h"

#include <omp.h>

void cpyVec3(double in[3], double out[3]){
    out[0] = in[0];
    out[1] = in[1];
    out[2] = in[2];
}

void getRiU ( double q[3], const double *matU, double r[3]){
    cpyVec3(q,r);
    itrans00(r,matU);
}

void getKnRungeKuta( double k[3], double val[10]){
  double norm;

  norm = getNorm(val[1],val[2],val[3]);

  k[0] = val[1]/norm;
  k[1] = val[2]/norm;
  k[2] = val[3]/norm;
}

int bondPath( int bcp, dataCritP *bondCrit, int ncp, dataCritP *nnucCrit,int *bonding,
              int *cellInfo, dataCube cube, dataRun param, double min0,
              const double *matU, char *name){


  int i,j,k,step;
  int iterp,itern,amico;
  int attractors,flagNeg,flagPos;
  int nucleo1,nucleo2,itmp;
  int cell1,cell2;
  double ratm[3],rij;
  double dist,norm,difmin,dij;
  double val[10];
  double matH[9],eval[3],evec[9];
  char nameOut[128],tmpname[128];

  double vec[3],vec2[3];
  double ri[3],qi[3],rn[3],qn[3];
  double qc[3],q[3],rc[3];

  double k1[3],k3[3],k4[3],k6[3];
  double a3 = 0.300*BPATH_EPS;
  double a4 = 0.600*BPATH_EPS;
  double a6 = 0.875*BPATH_EPS;
  double c1 = 0.097883598;
  double c3 = 0.40257649;
  double c4 = 0.21043771;
  double c6 = 0.289102202;

  double *coorAttr;
  FILE *tmp;
  FILE *out;
  /**********************************************************
    TODO  this is only a form to repair the gradient lines 
          (bpaths). In the future may do a structure for 
          matrix of transformation and rotation.
  **********************************************************/
  double matT[9];
  getMatInv(matU,matT);//The inverse of U is the original matrix T

  sprintf(nameOut,"%sBPath.xyz",name);

  tmpFile(&tmp,".c3dBpath_",tmpname,"w+");

  attractors = ncp + cube.natm;

  createArrayDou(3*attractors,&coorAttr,"BP01");

  for(i=0;i<cube.natm;i++){
    coorAttr[3*i]   = cube.coor[3*i];
    coorAttr[3*i+1] = cube.coor[3*i+1];
    coorAttr[3*i+2] = cube.coor[3*i+2];
  }
  for(j=0;j<ncp;j++){
    coorAttr[3*i]   = nnucCrit[j].x;
    coorAttr[3*i+1] = nnucCrit[j].y;
    coorAttr[3*i+2] = nnucCrit[j].z;
    i++;
  }

  printBar(stdout);

#pragma omp parallel private(i,nucleo1,nucleo2,val,matH,eval,evec,vec,vec2,    \
                             norm,ri,rn,qi,qn,iterp,itern,dist,qc,q,k,k1,k3,k4,\
                             k6,amico,difmin,step,ratm,rij,flagPos,flagNeg,    \
                             cell1, cell2,rc)                                  \
                     shared(bcp,bondCrit,cube,param,matU,matT,min0,a3,a4,      \
                            a6,c1,c3,c4,c6,attractors,coorAttr,tmp,cellInfo)
{

#pragma omp single 
{
  printf("  Number of threads in bond path  : %6d\n",omp_get_num_threads());
}
#pragma omp barrier 


#pragma omp for schedule (dynamic) 
  for(i=0; i < bcp; i++){
    nucleo1 = nucleo2 = 0;

    qc[0] = bondCrit[i].x;
    qc[1] = bondCrit[i].y;
    qc[2] = bondCrit[i].z;
    
    getRiU(qc,matU,rc);
    dij = 0.;
    fprintf(tmp,"%d % 10.6lf % 10.6lf % 10.6lf % 10.6lf\n",
                  i+1,rc[0]*B2A,rc[1]*B2A,rc[2]*B2A,dij);

    numCritical02Vec(qc,cube,param,matU,min0,val); 

    matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
    matH[3] = val[7]; matH[4] = val[5]; matH[5] = val[9];
    matH[6] = val[8]; matH[7] = val[9]; matH[8] = val[6];

    JacobiNxN(matH,eval,evec);
/*  The 3rd eigenvector is taken */
    vec2[0] = evec[6];  vec2[1] = evec[7];  vec2[2] = evec[8]; 
    
    norm = getNormVec(vec2);
    vec2[0] /= norm; vec2[1] /= norm; vec2[2] /= norm;

    ri[0] = rc[0] + 2.0 * BPATH_EPS * vec2[0];
    ri[1] = rc[1] + 2.0 * BPATH_EPS * vec2[1];
    ri[2] = rc[2] + 2.0 * BPATH_EPS * vec2[2];
    getRiU(ri,matT,qi);

    cellInfo[2*i  ] = 0;
    cellInfo[2*i+1] = 0;
   
    step = iterp= 0; dist = 6.;
    flagPos = 0;
    flagPos += perfectCube ( param.pbc,qi,cube.min,cube.max);
    while( iterp < MAXPTS && dist > 0.0 ){

      getRiU(qi,matU,ri);
      numCritical01Vec(qi,cube,param,matU,min0,val);
      getKnRungeKuta(k1,val);   

      q[0] = qi[0] + a3;
      q[1] = qi[1] + a3;
      q[2] = qi[2] + a3;

      flagPos += perfectCube ( param.pbc,q,cube.min,cube.max);
      numCritical01Vec(q,cube,param,matU,min0,val);

      getKnRungeKuta(k3,val);

      q[0] = qi[0] + a4;
      q[1] = qi[1] + a4;
      q[2] = qi[2] + a4;

      flagPos += perfectCube ( param.pbc,q,cube.min,cube.max);
      numCritical01Vec(q,cube,param,matU,min0,val);

      getKnRungeKuta(k4,val);

      q[0] = qi[0] + a6;
      q[1] = qi[1] + a6;
      q[2] = qi[2] + a6;

      flagPos += perfectCube ( param.pbc,q,cube.min,cube.max);
      numCritical01Vec(q,cube,param,matU,min0,val);

      getKnRungeKuta(k6,val);
      
      vec[0] = c1*k1[0] + c3*k3[0] + c4*k4[0] + c6*k6[0];
      vec[1] = c1*k1[1] + c3*k3[1] + c4*k4[1] + c6*k6[1];
      vec[2] = c1*k1[2] + c3*k3[2] + c4*k4[2] + c6*k6[2];

      //if( param.orth != YES) trans01(vec,matU);
      
      
      rn[0] = ri[0] + BPATH_EPS * vec[0];
      rn[1] = ri[1] + BPATH_EPS * vec[1];
      rn[2] = ri[2] + BPATH_EPS * vec[2];

      getRiU(rn,matT,qn);

      flagPos += perfectCube ( param.pbc,qn,cube.min,cube.max);
       
      step++;

      amico = 0;
      difmin = 1.E4;

      cpyVec3(qn,qi);

      for(k=0; k < attractors; k++){
        ratm[0] = coorAttr[3*k];
        ratm[1] = coorAttr[3*k+1];
        ratm[2] = coorAttr[3*k+2];

        rij = distance(ratm,qn);

        if( rij < difmin){
          difmin = rij;
          nucleo1 = k;
        }
        if( rij <= TOL_DIST_ATM){
          amico = 1;
        }

      }// end for k

      if( amico == 0 ){
        iterp++;
        dist = 6.;
        if( step == NSTEP ){
          if( myIsNanInf_V3(qn) == 0 ){
            getRiU(qn,matU,rn);
            dij = distance(rn,rc);
            fprintf(tmp,"%d % 10.6lf % 10.6lf % 10.6lf % 10.6lf\n",
                            i+1,rn[0]*B2A,rn[1]*B2A,rn[2]*B2A,dij);
          }
          step = 0;
        }
      }else{
        dist = -1.;
      }
    }//end while


 // At this moment we change the direction
    ri[0] = rc[0] - 2.0 * BPATH_EPS * vec2[0]; 
    ri[1] = rc[1] - 2.0 * BPATH_EPS * vec2[1];
    ri[2] = rc[2] - 2.0 * BPATH_EPS * vec2[2];
    getRiU(ri,matT,qi);

    step = itern = 0; dist = 6.;
    flagNeg = 0;
    flagNeg += perfectCube ( param.pbc,qi,cube.min,cube.max);
    while( itern < MAXPTS && dist > 0. ){

      getRiU(qi,matU,ri);
      numCritical01Vec(qi,cube,param,matU,min0,val);
      getKnRungeKuta(k1,val);

      q[0] = qi[0] + a3;
      q[1] = qi[1] + a3;
      q[2] = qi[2] + a3;

      flagNeg += perfectCube ( param.pbc,q,cube.min,cube.max);
      numCritical01Vec(q,cube,param,matU,min0,val);
      getKnRungeKuta(k3,val);

      q[0] = qi[0] + a4;
      q[1] = qi[1] + a4;
      q[2] = qi[2] + a4;

      flagNeg += perfectCube ( param.pbc,q,cube.min,cube.max);
      numCritical01Vec(q,cube,param,matU,min0,val);
      getKnRungeKuta(k4,val);

      q[0] = qi[0] + a6;
      q[1] = qi[1] + a6;
      q[2] = qi[2] + a6;

      flagNeg += perfectCube ( param.pbc,q,cube.min,cube.max);
      numCritical01Vec(q,cube,param,matU,min0,val);
      getKnRungeKuta(k6,val);

      vec[0] = c1*k1[0] + c3*k3[0] + c4*k4[0] + c6*k6[0];
      vec[1] = c1*k1[1] + c3*k3[1] + c4*k4[1] + c6*k6[1];
      vec[2] = c1*k1[2] + c3*k3[2] + c4*k4[2] + c6*k6[2];
      
      rn[0] = ri[0] + BPATH_EPS * vec[0];
      rn[1] = ri[1] + BPATH_EPS * vec[1];
      rn[2] = ri[2] + BPATH_EPS * vec[2];
 
      getRiU(rn,matT,qn);

      flagNeg += perfectCube ( param.pbc,qn,cube.min,cube.max);

      step++;

      amico = 0;
      difmin = 1.E4;

      cpyVec3(qn,qi);

      for(k=0; k < attractors; k++){
        if ( k ==  nucleo1 )
            continue;
        ratm[0] = coorAttr[3*k];
        ratm[1] = coorAttr[3*k+1];
        ratm[2] = coorAttr[3*k+2];

        rij = distance(ratm,qn);

        if( rij < difmin){
          difmin = rij;
          nucleo2 = k;
        }
        if( rij <= TOL_DIST_ATM){
          amico = 1;
        }

      }// end for k

      if( amico == 0 ){
        itern++;
        dist = 6.;

        if( step == NSTEP ){
          if( myIsNanInf_V3(qn) == 0 ){
            getRiU(qn,matU,rn);
            dij = -distance(rn,rc);
            fprintf(tmp,"%d % 10.6lf % 10.6lf % 10.6lf % 10.6lf\n",
                          i+1,rn[0]*B2A,rn[1]*B2A,rn[2]*B2A,dij);

          }
         step = 0;
        }
        //}
      }else{
        dist = -1.;
      }
    }//end while

    cell1 = flagPos;
    cell2 = flagNeg;

    if( nucleo1 > nucleo2){
      itmp = nucleo1;
      nucleo1 = nucleo2;
      nucleo2 = itmp;

      itmp = cell1;
      cell1 = cell2;
      cell2 = itmp;
    }

    bonding [2*i  ] = nucleo1;
    bonding [2*i+1] = nucleo2;
    cellInfo[2*i  ] = cell1;
    cellInfo[2*i+1] = cell2;
  }//fin for
}//fin OMP

  printBar(stdout);
  printBanner("Connectivity of Attractors by BCP",stdout);
  char symb1[6],symb2[6];
  char cc1,cc2;
  for( i = 0; i < bcp ; i++ ){
    
    cc1=' ';
    cc2=' ';
    if( cellInfo[2*i] != 0 )
        cc1='*';
    if( cellInfo[2*i+1] != 0 )
        cc2='*';

    nucleo1 = bonding[2*i];
    nucleo2 = bonding[2*i+1];

    if( nucleo1 < cube.natm )
      getAtomicSymbol(cube.zatm[nucleo1],4,symb1);
    else
      strcpy(symb1,"NNA");
    if( nucleo2 < cube.natm )
      getAtomicSymbol(cube.zatm[nucleo2],4,symb2);
    else
      strcpy(symb2,"NNA");
  
    printf("     Critical point %4d between  : %3d %c%-4s %3d %c%-4s\n",i+1,nucleo1+1,cc1,symb1,nucleo2+1,cc2,symb2);
  }

  rewind(tmp);

  int npoints=0;
  char c;

  while( (c = fgetc(tmp)) != EOF){
    if( c == '\n' )
      npoints++;
  }
  rewind(tmp);


  openFile(&out,nameOut,"w+");

  printBpaths(npoints,bcp,tmp,out);

  fclose(tmp);
  fclose(out);

  remove(tmpname);
  printBar(stdout);
  printf("  File %s was generated\n",nameOut);
  return 0;
}


int perfectCube ( int bcp, double *r, double *min, double *max){
  
  int flag = 0;
  if( bcp == YES ){
    double dx = fabs( max[0] - min[0]);
    double dy = fabs( max[1] - min[1]);
    double dz = fabs( max[2] - min[2]);

    if ( r[0] <= min[0] ){ r[0] += dx; flag += 1;}
    if ( r[0]  > max[0] ){ r[0] -= dx; flag += 1;}
        
    if ( r[1] <= min[1] ){ r[1] += dy; flag += 1;}
    if ( r[1]  > max[1] ){ r[1] -= dy; flag += 1;}

    if ( r[2] <= min[2] ){ r[2] += dz; flag += 1;}
    if ( r[2]  > max[2] ){ r[2] -= dz; flag += 1;}
  }

  return flag;
}

int logFile( int bcp, int rcp, int ccp, int ncp, dataCritP *bondCrit,
             dataCritP *ringCrit, dataCritP *cageCrit, dataCritP *nnucCrit,
             dataCube cube, dataRun param, double min0, int *bonding, int *cells,
             const double *matU, char *name){
  int i;
  int at1,at2;
  char nameOut[128];
  char symb1[6],symb2[6];
  double x,y,z;
  double fun, lap, val[10];
  double l1,l2,l3, matH[9];
  double eval[3];
  double r[3];
  double kinG,kinK,virial,eneH;
  double ngrad;
  char cc1,cc2;
  FILE *out;

  sprintf(nameOut,"%sCritP.log",name);
  openFile(&out,nameOut,"w+");

  printLog1(bcp,rcp,ccp,ncp,cube.natm,out);


  /////////////// IMPRIMIMOS INFORMACION DE LOS NUCLEOS //////////
  for(i=0; i < cube.natm; i++){
    x = cube.coor[3*i];
    y = cube.coor[3*i+1];
    z = cube.coor[3*i+2];
    numCritical02(x,y,z,cube,param,matU,min0,val);
    fun = val[0];
    lap = getLap(val);
    getAtomicSymbol(cube.zatm[i],4,symb1);
    fprintf(out,"%5d  %5s",i+1,symb1);
    r[0] = x*B2A; r[1] = y*B2A; r[2] = z*B2A;
    itrans00(r,matU);
    fprintf(out," % 10.6lf % 10.6lf % 10.6lf % 14.8lf % 16.7lf\n",r[0],r[1],r[2],fun,lap);
  }

  //// IMPRIMIMOS INFORMACION DE LOS ATTRACTORES NO NUCLEARES //////
  if( ncp > 0 ) {
    printBar82(out);
    centraMess("Non nuclear Attractor Critical Points",out);
    printBar82(out);
  }

  for(i=0;i<ncp;i++){
    x = nnucCrit[i].x;
    y = nnucCrit[i].y;
    z = nnucCrit[i].z;
    numCritical02(x,y,z,cube,param,matU,min0,val);
    fun = val[0];
    lap = getLap(val);
    fprintf(out,"%5d  %5s",i+1,"NNACP");
    r[0] = x*B2A; r[1] = y*B2A; r[2] = z*B2A;
    itrans00(r,matU);
    fprintf(out," % 10.6lf % 10.6lf % 10.6lf % 14.8lf % 16.7lf\n",r[0],r[1],r[2],fun,lap);
  }
  //// IMPRIMIMOS INFORMACION DE LOS DEMAS PUNTOS CRITICOS//////
  printBar82(out);
  centraMess("Bond, Ring and Cage critical points",out);
  printBar82(out);
  for(i=0;i<bcp;i++){
    x = bondCrit[i].x;
    y = bondCrit[i].y;
    z = bondCrit[i].z;
    numCritical02(x,y,z,cube,param,matU,min0,val);
    fun = val[0];
    lap = getLap(val);
    fprintf(out,"%5d  %5s",i+1,"BCP");
    r[0] = x*B2A; r[1] = y*B2A; r[2] = z*B2A;
    itrans00(r,matU);
    fprintf(out," % 10.6lf % 10.6lf % 10.6lf % 14.8lf % 16.7lf\n",r[0],r[1],r[2],fun,lap);
  }

  if( rcp > 0 ) 
    printBar82(out);
  for(i=0;i<rcp;i++){
    x = ringCrit[i].x;
    y = ringCrit[i].y;
    z = ringCrit[i].z;
    numCritical02(x,y,z,cube,param,matU,min0,val);
    fun = val[0];
    lap = getLap(val);
    fprintf(out,"%5d  %5s",i+1,"RCP");
    r[0] = x*B2A; r[1] = y*B2A; r[2] = z*B2A;
    itrans00(r,matU);
    fprintf(out," % 10.6lf % 10.6lf % 10.6lf % 14.8lf % 16.7lf\n",r[0],r[1],r[2],fun,lap);
  }
  if( ccp > 0 ) 
    printBar82(out);
  for(i=0;i<ccp;i++){
    x = cageCrit[i].x;
    y = cageCrit[i].y;
    z = cageCrit[i].z;
    numCritical02(x,y,z,cube,param,matU,min0,val);
    fun = val[0];
    lap = getLap(val);
    fprintf(out,"%5d  %5s",i+1,"CCP");
    r[0] = x*B2A; r[1] = y*B2A; r[2] = z*B2A;
    itrans00(r,matU);
    fprintf(out," % 10.6lf % 10.6lf % 10.6lf % 14.8lf % 16.7lf\n",r[0],r[1],r[2],fun,lap);
  }

  /////////////////////////////////////////INFORMACION DETALLADA 
  if( ncp > 0 ){
    printBar82(out);
    centraMess("No-nuclear Attractor Critical Points",out);
    for(i=0;i<ncp;i++){
      printBar82(out);
      x = nnucCrit[i].x;
      y = nnucCrit[i].y;
      z = nnucCrit[i].z;
      numCritical02(x,y,z,cube,param,matU,min0,val);

      fun  = val[0];
	   ngrad = getGrd(val);
	   lap   = getLap(val);

	   getEnergies(fun,lap,&kinG,&kinK,&virial);
	   eneH = kinG + virial;
	
      matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
      matH[3] = val[7]; matH[4] = val[5]; matH[5] = val[9];
      matH[6] = val[8]; matH[7] = val[9]; matH[8] = val[6];

      //JacobiNxN(matH,eval,evec);
      valoresPropios3x3(matH,eval);

	   l1 = eval[0]; l2 = eval[1]; l3 = eval[2];

      r[0] = x*B2A; r[1] = y*B2A; r[2] = z*B2A;
      itrans00(r,matU);

      fprintf(out,"\n   Critical Point  %7d of type (3,-3)\n\n",i+1);
      fprintf(out,"   Coordinates [A]    % 12.8E  % 12.8E  % 12.8E  \n\n",
                   r[0],r[1],r[2]);

	   fprintf(out,"   Density            % 12.8E\n",fun);
	   fprintf(out,"   NormGrad           % 12.8E\n",ngrad);
	   fprintf(out,"   Laplacian          % 12.8E\n\n",lap);

	   fprintf(out,"   Kinetic Energy G   % 12.8E\n",kinG);
	   fprintf(out,"   Kinetic Energy K   % 12.8E\n",kinK);
       fprintf(out,"   Virial field V     % 12.8E\n",virial);
	   fprintf(out,"   Total energy H     % 12.8E\n\n",eneH);

       fprintf(out,"                     | % 12.8E  % 12.8E  % 12.8E |\n"  ,matH[0],matH[1],matH[2]);
       fprintf(out,"   Hessian Matrix  = | % 12.8E  % 12.8E  % 12.8E |\n"  ,matH[3],matH[4],matH[5]);
       fprintf(out,"                     | % 12.8E  % 12.8E  % 12.8E |\n\n",matH[6],matH[7],matH[8]);


	   fprintf(out,"   Eigenvalues        % 12.8E  % 12.8E  % 12.8E\n\n",l1,l2,l3);
    }
  }


  if( bcp > 0 ){
    printBar82(out);
    centraMess("Bond Critical Points",out);
    for(i=0;i<bcp;i++){
      
      printBar82(out);
      x = bondCrit[i].x;
      y = bondCrit[i].y;
      z = bondCrit[i].z;
      numCritical02(x,y,z,cube,param,matU,min0,val);

      fun  = val[0];
	  ngrad = getGrd(val);
	  lap   = getLap(val);

	  getEnergies(fun,lap,&kinG,&kinK,&virial);
	  eneH = kinG + virial;
	
      matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
      matH[3] = val[7]; matH[4] = val[5]; matH[5] = val[9];
      matH[6] = val[8]; matH[7] = val[9]; matH[8] = val[6];

      //JacobiNxN(matH,eval,evec);
      valoresPropios3x3(matH,eval);
	  l1 = eval[0]; l2 = eval[1]; l3 = eval[2];
      
      at1 = bonding[2*i];
      at2 = bonding[2*i+1];

      if( at1 < cube.natm )
        getAtomicSymbol(cube.zatm[at1],4,symb1);
      else
        strcpy(symb1,"NNACP");

      if( at2 < cube.natm )
        getAtomicSymbol(cube.zatm[at2],4,symb2);
      else
        strcpy(symb2,"NNACP");

      at1++;
      at2++;

      cc1 = ' ';
      if( cells[2*i]   != 0 ) cc1 = '*';
      cc2 = ' '; 
      if( cells[2*i+1] != 0 ) cc2 = '*';
      

      r[0] = x*B2A; r[1] = y*B2A; r[2] = z*B2A;
      itrans00(r,matU);

      fprintf(out,"\n   Critical Point  %7d of type (3,-1)\n\n",i+1);
      fprintf(out,"   Between the nucleous : %2d %c%s and %2d %c%s\n",
                  at1,cc1,symb1,at2,cc2,symb2);
      fprintf(out,"   Coordinates [A]    % 12.8E  % 12.8E  % 12.8E  \n\n",
                   r[0],r[1],r[2]);

	   fprintf(out,"   Density            % 12.8E\n",fun);
	   fprintf(out,"   NormGrad           % 12.8E\n",ngrad);
	   fprintf(out,"   Laplacian          % 12.8E\n\n",lap);

	   fprintf(out,"   Kinetic Energy G   % 12.8E\n",kinG);
	   fprintf(out,"   Kinetic Energy K   % 12.8E\n",kinK);
      fprintf(out,"   Virial field V     % 12.8E\n",virial);
	   fprintf(out,"   Total energy H     % 12.8E\n\n",eneH);

      fprintf(out,"                     | % 12.8E  % 12.8E  % 12.8E |\n"  ,matH[0],matH[1],matH[2]);
      fprintf(out,"   Hessian Matrix  = | % 12.8E  % 12.8E  % 12.8E |\n"  ,matH[3],matH[4],matH[5]);
      fprintf(out,"                     | % 12.8E  % 12.8E  % 12.8E |\n\n",matH[6],matH[7],matH[8]);


	   fprintf(out,"   Eigenvalues        % 12.8E  % 12.8E  % 12.8E\n\n",l1,l2,l3);
	   fprintf(out,"   Bond ellipticity   %12.8lf\n",(l1/l2)-1.);
	   fprintf(out,"   Eta index          % 12.8E\n\n",fabs(l1)/l3);
    }
  }


  if( rcp > 0 ){
    printBar82(out);
    centraMess("Ring Critical Points",out);
    for(i=0;i<rcp;i++){
      printBar82(out);
      x = ringCrit[i].x;
      y = ringCrit[i].y;
      z = ringCrit[i].z;
      numCritical02(x,y,z,cube,param,matU,min0,val);

      fun  = val[0];
	   ngrad = getGrd(val);
	   lap   = getLap(val);

	   getEnergies(fun,lap,&kinG,&kinK,&virial);
	   eneH = kinG + virial;
	
      matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
      matH[3] = val[7]; matH[4] = val[5]; matH[5] = val[9];
      matH[6] = val[8]; matH[7] = val[9]; matH[8] = val[6];

      //JacobiNxN(matH,eval,evec);
      valoresPropios3x3(matH,eval);
	   l1 = eval[0]; l2 = eval[1]; l3 = eval[2];

      r[0] = x*B2A; r[1] = y*B2A; r[2] = z*B2A;
      itrans00(r,matU);

      fprintf(out,"\n   Critical Point  %7d of type (3, 1)\n\n",i+1);
      fprintf(out,"   Coordinates [A]    % 12.8E  % 12.8E  % 12.8E  \n\n",
                   r[0],r[1],r[2]);

	   fprintf(out,"   Density            % 12.8E\n",fun);
	   fprintf(out,"   NormGrad           % 12.8E\n",ngrad);
	   fprintf(out,"   Laplacian          % 12.8E\n\n",lap);

	   fprintf(out,"   Kinetic Energy G   % 12.8E\n",kinG);
	   fprintf(out,"   Kinetic Energy K   % 12.8E\n",kinK);
      fprintf(out,"   Virial field V     % 12.8E\n",virial);
	   fprintf(out,"   Total energy H     % 12.8E\n\n",eneH);

      fprintf(out,"                     | % 12.8E  % 12.8E  % 12.8E |\n"  ,matH[0],matH[1],matH[2]);
      fprintf(out,"   Hessian Matrix  = | % 12.8E  % 12.8E  % 12.8E |\n"  ,matH[3],matH[4],matH[5]);
      fprintf(out,"                     | % 12.8E  % 12.8E  % 12.8E |\n\n",matH[6],matH[7],matH[8]);


	   fprintf(out,"   Eigenvalues        % 12.8E  % 12.8E  % 12.8E\n\n",l1,l2,l3);
    }
  }

  if( ccp > 0 ){
    printBar82(out);
    centraMess("Cage Critical Points",out);
    for(i=0;i<ccp;i++){
      printBar82(out);
      x = cageCrit[i].x;
      y = cageCrit[i].y;
      z = cageCrit[i].z;
      numCritical02(x,y,z,cube,param,matU,min0,val);

      fun  = val[0];
	   ngrad = getGrd(val);
	   lap   = getLap(val);

	   getEnergies(fun,lap,&kinG,&kinK,&virial);
	   eneH = kinG + virial;
	
      matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
      matH[3] = val[7]; matH[4] = val[5]; matH[5] = val[9];
      matH[6] = val[8]; matH[7] = val[9]; matH[8] = val[6];

      //JacobiNxN(matH,eval,evec);
      valoresPropios3x3(matH,eval);
	   l1 = eval[0]; l2 = eval[1]; l3 = eval[2];

      r[0] = x*B2A; r[1] = y*B2A; r[2] = z*B2A;
      itrans00(r,matU);

      fprintf(out,"\n   Critical Point  %7d of type (3, 3)\n\n",i+1);
      fprintf(out,"   Coordinates [A]    % 12.8E  % 12.8E  % 12.8E  \n\n",
                   r[0],r[1],r[2]);

	   fprintf(out,"   Density            % 12.8E\n",fun);
	   fprintf(out,"   NormGrad           % 12.8E\n",ngrad);
	   fprintf(out,"   Laplacian          % 12.8E\n\n",lap);

	   fprintf(out,"   Kinetic Energy G   % 12.8E\n",kinG);
	   fprintf(out,"   Kinetic Energy K   % 12.8E\n",kinK);
      fprintf(out,"   Virial field V     % 12.8E\n",virial);
	   fprintf(out,"   Total energy H     % 12.8E\n\n",eneH);

      fprintf(out,"                     | % 12.8E  % 12.8E  % 12.8E |\n"  ,matH[0],matH[1],matH[2]);
      fprintf(out,"   Hessian Matrix  = | % 12.8E  % 12.8E  % 12.8E |\n"  ,matH[3],matH[4],matH[5]);
      fprintf(out,"                     | % 12.8E  % 12.8E  % 12.8E |\n\n",matH[6],matH[7],matH[8]);


	   fprintf(out,"   Eigenvalues        % 12.8E  % 12.8E  % 12.8E\n\n",l1,l2,l3);
    }
  }

  printBar82(out);
  centraMess("End of File",out);
  printBar82(out);

  fclose(out);

  printBar(stdout);
  printf("  File %s was generated\n",nameOut);

  return 0;
}

  
void getEnergies(double rho, double lap, double *kinG,double *kinK,double *vir){

  double kG,kK;

  kG = (3./10.)*pow(3.*M_PI*M_PI,2./3.)*pow(rho,5./3.);
  kG += (1./6.)*lap;

  kK = kG - 0.25*lap;


  (*vir) = 0.25*lap - 2.*kG;

  (*kinG) = kG;
  (*kinK) = kK;


}

int myIsNanInf(double val){

  int ret=0;

  ret += isnan(val);
  ret += isinf(val);

  return ret; 
}

int myIsNanInf_V3(double vec[3]){

  int ret = 0; 
  
  ret += myIsNanInf(vec[0]);
  ret += myIsNanInf(vec[1]);
  ret += myIsNanInf(vec[2]);

  return ret; 

}

void printLog1(int nbc, int nrc, int ncc, int nna, int nnu, FILE *out){
  
  char tchar[128];

  time_t t= time(NULL);
  struct tm *tlocal = localtime(&t);

  strftime(tchar,128,"%H:%M:%S  %d/%m/%y",tlocal);


  printBar82(out);

  centraMess("File log for cube3d",out);
  printBar82(out);
  fprintf(out,"%82s\n",tchar);
  printBar82(out);
  
  fprintf(out," This file contains information about critical points and properties, these\n");
  fprintf(out," information are determineted by code Cube3D-%s\n",VERSION);
  fprintf(out," The units in this  file for  distance are  Angstrom and  the units for fields are\n");
  fprintf(out," atomic units.\n\n");
  fprintf(out,"  NNACP               Non-nuclear attractor critical point\n");
  fprintf(out,"  BCP                 Bond critical point\n");
  fprintf(out,"  RCP                 Ring critical point\n");
  fprintf(out,"  CCP                 Cage critical point\n");
  fprintf(out,"  Kinetic Energy G    Kinetic Energy Density in the form Lagrangian\n");
  fprintf(out,"                      (Abramov's approximation)\n");
  fprintf(out,"  Kinetic Energy K    Kinetic Energy Density in the form Hamiltonian\n");
  fprintf(out,"  Virial Field        Or Potential Energy Density  V = -K - G\n");
  //fprintf(out,"  Lagrangian Density  -(1/4) Laplacian or K - G\n");
 // fprintf(out,"  KenergyG/Density    Kinetic Energy Density per electron\n");
  
  fprintf(out,"  Eigenvalues         Eigenvalues for Hessian matrix: lambda_1<lambda_2<lambda_3\n");
  fprintf(out,"  Bond Ellipiticy     lambda_1/lambda_2 - 1\n");
  fprintf(out,"  Eta index           |lambda_1|/lambda_3\n\n");
  printBar82(out);
  centraMess("Information about the system",out);
  printBar82(out);
  fprintf(out,"    Nuclei  (3,-3)         : %15d\n",nnu);
  fprintf(out,"    NNACP   (3,-3)         : %15d\n",nna);
  fprintf(out,"    BCP     (3,-1)         : %15d\n",nbc);
  fprintf(out,"    RCP     (3,+1)         : %15d\n",nrc);
  fprintf(out,"    CCP     (3,+3)         : %15d\n",ncc);
  fprintf(out,"    Total Critical Points  : %15d\n",nna+nbc+nrc+ncc);
  printBar82(out);
  fprintf(out,"    Poincare-Hopf rule\n");
  fprintf(out,"    nuclei + NNACP + RCP - BCP - CCP = %5d\n",nnu+nna+nrc-nbc-ncc);
  printBar82(out);

  centraMess("Nuclear and Non-nuclear Attractors",out);
  printBar82(out);
  fprintf(out,"  Atractor                Coordinates [A]            Density         Laplacian\n");
  printBar82(out);



}

void centraMess(char *mess,FILE *out){
  int nmax = 82;
  int i, n = strlen(mess);
  int spacio;
  if( n%2 == 1)
    n++;
  if(nmax%2 == 1)
    nmax++;

  spacio = nmax/2 - n/2;

  if( spacio > 0 )
    for(i=0;i<spacio;i++)
      fprintf(out," ");
  fprintf(out,"%s\n",mess);

}

int axesCrit( int bcp, int rcp, int ccp, int ncp, dataCritP *bondCrit,
              dataCritP *ringCrit, dataCritP *cageCrit, dataCritP *nnucCrit,
             dataCube cube, dataRun param, double min0,  const double *matU, char *name){


  int i,j;
  char nameOut[128];

  FILE *out;
  sprintf(nameOut,"%sAxes.xyz",name);

  openFile(&out,nameOut,"w+");

  double x0,y0,z0;
  double v1,v2,v3,norm;
  double val[10];
  double matH[9];
  double eval[3];
  double evec[9];
  double r[3];

  for(i=0; i < bcp; i++){
    x0 = bondCrit[i].x;
    y0 = bondCrit[i].y;
    z0 = bondCrit[i].z;
    numCritical02(x0,y0,z0,cube,param,matU,min0,val);
  
    matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
    matH[3] = val[7]; matH[4] = val[5]; matH[5] = val[9];
    matH[6] = val[8]; matH[7] = val[9]; matH[8] = val[6];
  
    JacobiNxN(matH,eval,evec);

    r[0] = x0;
    r[1] = y0;
    r[2] = z0;
    itrans00(r,matU);

    fprintf(out," BCP  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);

    v1 = evec[6]; v2 = evec[7]; v3 = evec[8];
    norm = getNorm(v1,v2,v3);
    v1 = v1*0.05/norm; v2 = v2*0.05/norm; v3 = v3*0.05/norm;

    for(j=1;j<=5;j++){
      r[0] = x0 + j*v1;
      r[1] = y0 + j*v2; 
      r[2] = z0 + j*v3;
      itrans00(r,matU);
      fprintf(out," LA3B  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    for(j=1;j<=5;j++){
      r[0] = x0 - j*v1;
      r[1] = y0 - j*v2; 
      r[2] = z0 - j*v3;
      itrans00(r,matU);
      fprintf(out," LA3B  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }

    v1 = evec[3]; v2 = evec[4]; v3 = evec[5];
    norm = getNorm(v1,v2,v3);
    v1 = v1*0.05/norm; v2 = v2*0.05/norm; v3 = v3*0.05/norm;
    for(j=1;j<=5;j++){
      r[0] = x0 + j*v1;
      r[1] = y0 + j*v2; 
      r[2] = z0 + j*v3;
      itrans00(r,matU);
      fprintf(out," LA2B  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    for(j=1;j<=5;j++){
      r[0] = x0 - j*v1;
      r[1] = y0 - j*v2; 
      r[2] = z0 - j*v3;
      itrans00(r,matU);
      fprintf(out," LA2B  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    v1 = evec[0]; v2 = evec[1]; v3 = evec[2];
    norm = getNorm(v1,v2,v3);
    v1 = v1*0.05/norm; v2 = v2*0.05/norm; v3 = v3*0.05/norm;
    for(j=1;j<=5;j++){
      r[0] = x0 + j*v1;
      r[1] = y0 + j*v2; 
      r[2] = z0 + j*v3;
      itrans00(r,matU);
      fprintf(out," LA1B  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    for(j=1;j<=5;j++){
      r[0] = x0 - j*v1;
      r[1] = y0 - j*v2; 
      r[2] = z0 - j*v3;
      itrans00(r,matU);
      fprintf(out," LA1B  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
  
  } // fin BCP

  for(i=0; i < rcp; i++){
    x0 = ringCrit[i].x;
    y0 = ringCrit[i].y;
    z0 = ringCrit[i].z;
    numCritical02(x0,y0,z0,cube,param,matU,min0,val);
  
    matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
    matH[3] = val[7]; matH[4] = val[5]; matH[5] = val[9];
    matH[6] = val[8]; matH[7] = val[9]; matH[8] = val[6];
  
    JacobiNxN(matH,eval,evec);

    r[0] = x0;
    r[1] = y0;
    r[2] = z0;
    itrans00(r,matU);

    fprintf(out," RCP  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);

    v1 = evec[6]; v2 = evec[7]; v3 = evec[8];
    norm = getNorm(v1,v2,v3);
    v1 = v1*0.05/norm; v2 = v2*0.05/norm; v3 = v3*0.05/norm;

    for(j=1;j<=5;j++){
      r[0] = x0 + j*v1;
      r[1] = y0 + j*v2; 
      r[2] = z0 + j*v3;
      itrans00(r,matU);
      fprintf(out," LA3R  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    for(j=1;j<=5;j++){
      r[0] = x0 - j*v1;
      r[1] = y0 - j*v2; 
      r[2] = z0 - j*v3;
      itrans00(r,matU);
      fprintf(out," LA3R  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }

    v1 = evec[3]; v2 = evec[4]; v3 = evec[5];
    norm = getNorm(v1,v2,v3);
    v1 = v1*0.05/norm; v2 = v2*0.05/norm; v3 = v3*0.05/norm;
    for(j=1;j<=5;j++){
      r[0] = x0 + j*v1;
      r[1] = y0 + j*v2; 
      r[2] = z0 + j*v3;
      itrans00(r,matU);
      fprintf(out," LA2R  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    for(j=1;j<=5;j++){
      r[0] = x0 - j*v1;
      r[1] = y0 - j*v2; 
      r[2] = z0 - j*v3;
      itrans00(r,matU);
      fprintf(out," LA2R  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    v1 = evec[0]; v2 = evec[1]; v3 = evec[2];
    norm = getNorm(v1,v2,v3);
    v1 = v1*0.05/norm; v2 = v2*0.05/norm; v3 = v3*0.05/norm;
    for(j=1;j<=5;j++){
      r[0] = x0 + j*v1;
      r[1] = y0 + j*v2; 
      r[2] = z0 + j*v3;
      itrans00(r,matU);
      fprintf(out," LA1R  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    for(j=1;j<=5;j++){
      r[0] = x0 - j*v1;
      r[1] = y0 - j*v2; 
      r[2] = z0 - j*v3;
      itrans00(r,matU);
      fprintf(out," LA1R  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
  
  } // fin RCP

  for(i=0; i < ccp; i++){
    x0 = cageCrit[i].x;
    y0 = cageCrit[i].y;
    z0 = cageCrit[i].z;
    numCritical02(x0,y0,z0,cube,param,matU,min0,val);
  
    matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
    matH[3] = val[7]; matH[4] = val[5]; matH[5] = val[9];
    matH[6] = val[8]; matH[7] = val[9]; matH[8] = val[6];
  
    JacobiNxN(matH,eval,evec);

    r[0] = x0;
    r[1] = y0;
    r[2] = z0;
    itrans00(r,matU);

    fprintf(out," CCP  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);

    v1 = evec[6]; v2 = evec[7]; v3 = evec[8];
    norm = getNorm(v1,v2,v3);
    v1 = v1*0.05/norm; v2 = v2*0.05/norm; v3 = v3*0.05/norm;

    for(j=1;j<=5;j++){
      r[0] = x0 + j*v1;
      r[1] = y0 + j*v2; 
      r[2] = z0 + j*v3;
      itrans00(r,matU);
      fprintf(out," LA3C  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    for(j=1;j<=5;j++){
      r[0] = x0 - j*v1;
      r[1] = y0 - j*v2; 
      r[2] = z0 - j*v3;
      itrans00(r,matU);
      fprintf(out," LA3C  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }

    v1 = evec[3]; v2 = evec[4]; v3 = evec[5];
    norm = getNorm(v1,v2,v3);
    v1 = v1*0.05/norm; v2 = v2*0.05/norm; v3 = v3*0.05/norm;
    for(j=1;j<=5;j++){
      r[0] = x0 + j*v1;
      r[1] = y0 + j*v2; 
      r[2] = z0 + j*v3;
      itrans00(r,matU);
      fprintf(out," LA2C  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    for(j=1;j<=5;j++){
      r[0] = x0 - j*v1;
      r[1] = y0 - j*v2; 
      r[2] = z0 - j*v3;
      itrans00(r,matU);
      fprintf(out," LA2C  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    v1 = evec[0]; v2 = evec[1]; v3 = evec[2];
    norm = getNorm(v1,v2,v3);
    v1 = v1*0.05/norm; v2 = v2*0.05/norm; v3 = v3*0.05/norm;
    for(j=1;j<=5;j++){
      r[0] = x0 + j*v1;
      r[1] = y0 + j*v2; 
      r[2] = z0 + j*v3;
      itrans00(r,matU);
      fprintf(out," LA1C  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    for(j=1;j<=5;j++){
      r[0] = x0 - j*v1;
      r[1] = y0 - j*v2; 
      r[2] = z0 - j*v3;
      itrans00(r,matU);
      fprintf(out," LA1C  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
  
  } // fin CCP

  for(i=0; i < ncp; i++){
    x0 = nnucCrit[i].x;
    y0 = nnucCrit[i].y;
    z0 = nnucCrit[i].z;
    numCritical02(x0,y0,z0,cube,param,matU,min0,val);
  
    matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
    matH[3] = val[7]; matH[4] = val[5]; matH[5] = val[9];
    matH[6] = val[8]; matH[7] = val[9]; matH[8] = val[6];
  
    JacobiNxN(matH,eval,evec);

    r[0] = x0;
    r[1] = y0;
    r[2] = z0;
    itrans00(r,matU);

    fprintf(out," NCP  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);

    v1 = evec[6]; v2 = evec[7]; v3 = evec[8];
    norm = getNorm(v1,v2,v3);
    v1 = v1*0.05/norm; v2 = v2*0.05/norm; v3 = v3*0.05/norm;

    for(j=1;j<=5;j++){
      r[0] = x0 + j*v1;
      r[1] = y0 + j*v2; 
      r[2] = z0 + j*v3;
      itrans00(r,matU);
      fprintf(out," LA3N  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    for(j=1;j<=5;j++){
      r[0] = x0 - j*v1;
      r[1] = y0 - j*v2; 
      r[2] = z0 - j*v3;
      itrans00(r,matU);
      fprintf(out," LA3N  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }

    v1 = evec[3]; v2 = evec[4]; v3 = evec[5];
    norm = getNorm(v1,v2,v3);
    v1 = v1*0.05/norm; v2 = v2*0.05/norm; v3 = v3*0.05/norm;
    for(j=1;j<=5;j++){
      r[0] = x0 + j*v1;
      r[1] = y0 + j*v2; 
      r[2] = z0 + j*v3;
      itrans00(r,matU);
      fprintf(out," LA2N  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    for(j=1;j<=5;j++){
      r[0] = x0 - j*v1;
      r[1] = y0 - j*v2; 
      r[2] = z0 - j*v3;
      itrans00(r,matU);
      fprintf(out," LA2N  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    v1 = evec[0]; v2 = evec[1]; v3 = evec[2];
    norm = getNorm(v1,v2,v3);
    v1 = v1*0.05/norm; v2 = v2*0.05/norm; v3 = v3*0.05/norm;
    for(j=1;j<=5;j++){
      r[0] = x0 + j*v1;
      r[1] = y0 + j*v2; 
      r[2] = z0 + j*v3;
      itrans00(r,matU);
      fprintf(out," LA1N  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
    for(j=1;j<=5;j++){
      r[0] = x0 - j*v1;
      r[1] = y0 - j*v2; 
      r[2] = z0 - j*v3;
      itrans00(r,matU);
      fprintf(out," LA1N  % 10.6lf % 10.6lf % 10.6lf\n",r[0]*B2A,r[1]*B2A,r[2]*B2A);
    }
  
  } // fin NCP

}


int logFileCSV( int ncp, dataCritP *crit, dataCube cube, 
                dataRun param, double min0, 
                const double *matU, char *name, char *string){

  int i;
  char nameOut[128];
  double fun,lap,val[10];
  double ngrad;
  double l1,l2,l3,matH[9];
  double eval[3],evec[9];
  double r[3],q[3],ellep;
  double kinG,kinK,virial,eneH,eneHb;
  FILE *out;

  sprintf(nameOut,"%sCritP_%s.csv",name,string);
  openFile(&out,nameOut,"w+");
  
  fprintf(out,"Index, Type, Coor X, Coor Y, Coor Z, Density, Ngrad, ");
  fprintf(out,"Laplacian, Kinetic G, Kinetic K, Virial, energy H, ");
  fprintf(out,"lambda_1, lambda_2, lambda_3, ellepticiy, energyInt(kcal/mol)\n");

  for(i=0;i<ncp;i++){
      q[0] = crit[i].x;
      q[1] = crit[i].y;
      q[2] = crit[i].z;
      numCritical02Vec_exp(q,cube,param,matU,min0,val);

      fun  = val[0];
	  lap   = getLap(val);
      ngrad = sqrt( val[1]*val[1] + val[2]*val[2] + val[3]*val[3]);
	  getEnergies(fun,lap,&kinG,&kinK,&virial);
	  eneH = kinG + virial;
	
      matH[0] = val[4]; matH[1] = val[7]; matH[2] = val[8];
      matH[3] = val[7]; matH[4] = val[5]; matH[5] = val[9];
      matH[6] = val[8]; matH[7] = val[9]; matH[8] = val[6];

      JacobiNxN(matH,eval,evec);
      //valoresPropios3x3(matH,eval);
	   l1 = eval[0]; l2 = eval[1]; l3 = eval[2];
      ellep = (l1/l2)-1.;
      eneHb = 0.5 * virial * 627.509;
      
      getRiU(q,matU,r);
      r[0] *= B2A;
      r[1] *= B2A;
      r[2] *= B2A;
      fprintf(out,"%3d,%s,%10.6lf,%10.6lf,%10.6lf,",i+1,string,r[0],r[1],r[2]);
      fprintf(out,"%10.6lf,%10.6lf,%10.6lf,%10.6lf,%10.6lf,",fun,ngrad,lap,kinG,kinK);
      fprintf(out,"%10.6lf,%10.6lf,",virial,eneH);
      fprintf(out,"%10.6lf,%10.6lf,%10.6lf,",l1,l2,l3);
      fprintf(out,"%10.6lf,%10.6lf\n",ellep,eneHb);

  }

  fclose(out);

  printBar(stdout);
  printf("  File %s was generated\n",nameOut);

  return 0;
}

void printBpaths(int npoints, int nbcp, FILE *tmp, FILE *out){
    int *type=NULL;
    dataBpath *bpaths=NULL;
    double dij,rx,ry,rz;
    int i,bond;

    rewind(tmp);

    createArrayInt(nbcp,&type,"print Bond Paths");
    createArrayDataBpath(npoints,&bpaths,"bond paths");

    for(i=0;i<nbcp;i++)
        type[i] = 0;

    for(i=0;i<npoints;i++){
        fscanf(tmp,"%d %lf %lf %lf %lf",&bond,&rx,&ry,&rz,&dij);
        bpaths[i].bcp  = bond;
        bpaths[i].r[0] = rx;
        bpaths[i].r[1] = ry;
        bpaths[i].r[2] = rz;
        bpaths[i].dij  = dij;
        type[bond-1]++;
    }

    sortBpaths(npoints,nbcp,type,bpaths);

    fprintf(out," %6d\n",npoints);
    fprintf(out," Bond Path file Coordinates in Angstroms \n");
    for(i=0;i<npoints;i++)
        fprintf(out," BP%04d  % 10.6lf % 10.6lf % 10.6lf \n",
                    bpaths[i].bcp,
                    bpaths[i].r[0], bpaths[i].r[1], bpaths[i].r[2]);
                 //   bpaths[i].r[0], bpaths[i].r[1], bpaths[i].r[2],bpaths[i].dij);


    free(type);
    free(bpaths);
}

void sortBpaths(int np, int nbcp, int* type, dataBpath* bpaths){
    int k;
    int min,max;

    mergeBPathBCP(bpaths,0,np-1);
    
    min = 0;
    for(k=0;k<nbcp;k++){
        max = min + type[k];
        mergeBPathDist(bpaths,min,max-1);
        min = max;
    }
}

void createArrayDataBpath(int n, dataBpath **ptr, char *mess){
    if( n <= 0 ){
        printf(" There is an error in the size (%d) for [%s] \n",n, mess);
        exit(EXIT_FAILURE);
    }

    (*ptr) = (dataBpath*) malloc((n)*sizeof(dataBpath));
    if( (*ptr) == NULL){
        printf(" Failed to allocate memory: [%s]\n",mess);
        exit(EXIT_FAILURE);
    }
}

void mergeBPathDist(dataBpath *array, int l, int r){
    if(l < r ){
        int m  = l + (r-l)/2;

        mergeBPathDist(array,l,m);
        mergeBPathDist(array,m+1,r);

        sortBPathDist(array,l,m,r);

    }

}

void sortBPathDist( dataBpath *array, int l, int m, int r){
    int i, j, k; 
    int n1 = m - l + 1; 
    int n2 =  r - m; 
    dataBpath *leftArray, *rightArray;
  
    createArrayDataBpath(n1+1,&leftArray,"Temporal array L");
    createArrayDataBpath(n2+1,&rightArray,"Temporal array R");
  
    for (i = 0; i < n1; i++) 
        leftArray[i] = array[l + i]; 
    for (j = 0; j < n2; j++) 
        rightArray[j] = array[m + 1+ j]; 
  
    i = 0; 
    j = 0; 
    k = l;
    while (i < n1 && j < n2) 
    { 
        if (leftArray[i].dij <= rightArray[j].dij) 
        { 
            array[k] = leftArray[i]; 
            i++; 
        } 
        else
        { 
            array[k] = rightArray[j]; 
            j++; 
        } 
        k++; 
    } 
  
    while (i < n1) 
    { 
        array[k] = leftArray[i]; 
        i++; 
        k++; 
    } 
  
    while (j < n2) 
    { 
        array[k] = rightArray[j]; 
        j++; 
        k++; 
    } 

    free(rightArray);
    free(leftArray);
}

void mergeBPathBCP(dataBpath *array, int l, int r){

    if ( l < r ){
        int m  = l + (r - l) / 2;
        mergeBPathBCP(array,l,m);
        mergeBPathBCP(array,m+1,r);

        sortBPathBCP(array,l,m,r);

    }

}

void sortBPathBCP( dataBpath *array, int l, int m, int r){
    int i, j, k; 
    int n1 =  m - l + 1; 
    int n2 =  r - m; 
    dataBpath *leftArray;
    dataBpath *rightArray;
  
    createArrayDataBpath(n1,&leftArray,"Temporal array L");

    createArrayDataBpath(n2,&rightArray,"Temporal array R");
  
    for (i = 0; i < n1; i++) 
        leftArray[i] = array[l + i]; 

    for (j = 0; j < n2; j++) 
        rightArray[j] = array[m + 1 + j]; 
  
    i = 0; 
    j = 0; 
    k = l;
    while (i < n1 && j < n2){ 
        if (leftArray[i].bcp <= rightArray[j].bcp){
            array[k] = leftArray[i]; 
            i++; 
        }else{ 
            array[k] = rightArray[j]; 
            j++; 
        } 
        k++; 
    } 
  
    while (i < n1){ 
        array[k] = leftArray[i]; 
        i++; 
        k++; 
    } 
  
    while (j < n2){ 
        array[k] = rightArray[j]; 
        j++; 
        k++; 
    } 


    free(leftArray);
    free(rightArray);
}
