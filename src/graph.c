#include "graph.h"
#include "file.h"

 void  InitPlot( FILE *gp, char *name, char *background,
               int type, int npx, int npy, int key){

  fprintf(gp,"set term qt size %d,%d enhanced font \'Helvetica,11\' title \"%s\"\n",npx,npy,name);
  fprintf(gp,"set object rectangle from screen 0,0 to screen 1,1 behind fillcolor rgb \"%s\" fillstyle solid noborder\n",background);
  fprintf(gp,"set encoding utf8\n");
  if( key == 0 ) 
    fprintf(gp,"unset key\n");
 
  fprintf(gp,"set style line 1 lc rgb '#0000DE' pt 7 ps 1.5 lt 1 lw 2\n"); 
  fprintf(gp,"set style line 2 lc rgb '#006400' pt 7 ps 1.5 lt 1 lw 2\n"); 
  fprintf(gp,"set style line 3 lc rgb '#B0BEC5' lt 1 lw 0.1 \n");
  fprintf(gp,"set style line 4 lc rgb '#0000DE' lt 1 lw 2\n"); 

  if ( key == 4){
    fprintf(gp,"set style line 101 lc rgb '#FFFFFF' lt 1 lw 1.25\n"); 
    fprintf(gp,"set key textcolor rgb '#FFFFFF'\n");
  }else
    fprintf(gp,"set style line 101 lc rgb \'#263238\' lt 1 lw 1.25\n"); 

  fprintf(gp,"set style line 201 lt 1 pt 5 ps 0.75 lw 1\n");
  fprintf(gp,"set style line 301 lt 1 pt 7 ps 0.25 lw 1\n");
  fprintf(gp,"set style line 401 lt 1 pt 7 ps 0.25 lw 2 lc rgb \'#FFF59D\'\n");

 
  fprintf(gp,"set tics nomirror out scale 0.75\n");
 
  if(type == _3D)
    fprintf(gp,"set border 21 front ls 101\n");
  else
    if(type == _2D)
      fprintf(gp,"set border 3 front ls 101\n");
    else
      fprintf(gp,"set border  front ls 101\n");

    
  fflush(gp);
 
}

void InitPlotWx( FILE *gp, char *name, char *background,
                 int type, int npx, int npy, int key){

  fprintf(gp,"set term wxt size %d,%d enhanced font \'Helvetica,11\' title \"%s\"\n",npx,npy,name);
  fprintf(gp,"set object rectangle from screen 0,0 to screen 1,1 behind fillcolor rgb \"%s\" fillstyle solid noborder\n",background);
  fprintf(gp,"set encoding utf8\n");
  if( key == 0 ) 
    fprintf(gp,"unset key\n");
 
  fprintf(gp,"set style line 1 lc rgb '#0000DE' pt 7 ps 1.5 lt 1 lw 2\n"); 
  fprintf(gp,"set style line 2 lc rgb '#006400' pt 7 ps 1.5 lt 1 lw 2\n"); 
  fprintf(gp,"set style line 3 lc rgb '#B0BEC5' lt 1 lw 0.1 \n");
  fprintf(gp,"set style line 4 lc rgb '#0000DE' lt 1 lw 2\n"); 

  if ( key == 4){
    fprintf(gp,"set style line 101 lc rgb '#FFFFFF' lt 1 lw 1.25\n"); 
    fprintf(gp,"set key textcolor rgb '#FFFFFF'\n");
  }else
    fprintf(gp,"set style line 101 lc rgb \'#263238\' lt 1 lw 1.25\n"); 

  fprintf(gp,"set style line 201 lt 1 pt 5 ps 0.75 lw 1\n");
  fprintf(gp,"set style line 301 lt 1 pt 7 ps 0.25 lw 1\n");
  fprintf(gp,"set style line 401 lt 1 pt 7 ps 0.25 lw 2 lc rgb \'#FFF59D\'\n");

 
  fprintf(gp,"set tics nomirror out scale 0.75\n");
 
  if(type == _3D)
    fprintf(gp,"set border 21 front ls 101\n");
  else
    if(type == _2D)
      fprintf(gp,"set border 3 front ls 101\n");
    else
      fprintf(gp,"set border  front ls 101\n");

    
  fflush(gp);
 
}



void plotNCIdat(char *name,dataRun param ){
  
  FILE *gp;
  double ymax = param.rgd;
  double xmax = param.la2;

  gp = popen(GNUPLOT_PATH,"w");

  InitPlot(gp,"Plot  NCI","#0F0F0F",_2D,1.5*NPX,1.5*NPY,4);
  fprintf(gp," set title  'Non-covalentet Interactions index' tc rgb 'white'\n");
  fprintf(gp," set xlabel 'sign({/Symbol l}_2){/Symbol r}(r)' tc rgb 'white'\n");
  fprintf(gp," set ylabel 's(r)' tc rgb 'white' \n");
  fprintf(gp," unset key \n");
  fprintf(gp," set palette \n");
  fprintf(gp," set palette defined ( 1 '#0000FF', 2 '#00FF00', 3 '#FF0000' ) \n");
  fprintf(gp," unset colorbox \n");
  fprintf(gp," set xrange[%f:%f]\n",-xmax,xmax);
  fprintf(gp," set yrange[%f:%f]\n",0.,ymax);
  fprintf(gp," set mxtics 10 \n");
  fprintf(gp," set mytics 4 \n");
  fprintf(gp," set format x \"%c4.3f\"\n",'%');
  fprintf(gp," set format y \"%c4.2f\"\n",'%');

  fprintf(gp," plot '%s.dat' u ($1/100.):2:1 w p pt 7 ps 0.45 lc palette\n",name);

  fflush(gp);
  fclose(gp);

}

void getNCIgnu(char *name, dataRun param){

  char tmpname[120];
  FILE *tmp;
  double ymax = fabs(param.rgd);
  double xmax = fabs(param.la2);
  
  

  sprintf(tmpname,"plot_NCI_%s.gp",name);
  openFile(&tmp,tmpname,"w+");

  fprintf(tmp,"#!/usr/bin/gnuplot\n");
  fprintf(tmp," reset\n");

  fprintf(tmp," set term post enhanced color eps size 15cm,10cm font 'Arial,30' \n");
  fprintf(tmp," set output '%s.eps'\n",name);
  fprintf(tmp," set encoding utf8\n");
  fprintf(tmp," unset key\n");
  fprintf(tmp," set style line 500 lc rgb '#303030' lt 1 lw 1.25\n"); 

  fprintf(tmp," set tics font 'Arial-Italic,16'\n");
  fprintf(tmp," set xlabel font 'Arial,16'\n");
  fprintf(tmp," set ylabel font 'Arial,16'\n");

  fprintf(tmp," set title  'Non-covalentet Interactions index' \n");
  fprintf(tmp," set xlabel 'sign({/Symbol l}_2){/Symbol r}(r)' \n");
  fprintf(tmp," set ylabel 's(r)' \n");
  fprintf(tmp," set palette \n");
  fprintf(tmp," set palette defined ( 1 '#0000FF', 2 '#00FF00', 3 '#FF0000' ) \n");
  fprintf(tmp," unset colorbox \n");
  fprintf(tmp," set xrange[%f:%f]\n",-xmax,xmax);
  fprintf(tmp," set yrange[%f:%f]\n",0.,ymax);
  fprintf(tmp," set mxtics 4\n");
  fprintf(tmp," set mytics 4\n");
  
  fprintf(tmp," set tics nomirror out scale 0.75\n");
  fprintf(tmp," set border 3 front ls 500\n");
  fprintf(tmp," set format x \"%c4.3f\"\n",'%');
  fprintf(tmp," set format y \"%c4.2f\"\n",'%');

  fprintf(tmp," plot '%s.dat' u ($1/100.):2:1 w p pt 7 ps 0.5 lc palette\n",name);
 
  fclose(tmp);


}

