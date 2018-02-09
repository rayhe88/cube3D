#include "tableP.h"


char chemSymbol[TPNMAX][3] = {
"H",                                                                                 "HE",
"LI","BE",                                                  "B" ,"C" ,"N" ,"O" ,"F" ,"NE",
"NA","MG",                                                  "AL","SI","P" ,"S" ,"CL","AR",
"K" ,"CA","SC","TI","V" ,"CR","MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR","KR",
"RB","SR","Y" ,"ZR","NB","MO","TC","RU","RH","PD","AG","CD","ID","SN","SB","TE","I" ,"XE",
"CS","BA","LA","HF","TA","W" ,"RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","AT","RN",
"FR","RA","AC","RF","DB","SG","BH","HS","MT","DS","RG","CN","--","--","--","--","--","--",
"CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YT","LU",
"TH","PA","U" ,"NP","PU","AM","CM","BK","CF","ES","FM","MD","NO","LR",
"BQ"
 };
char chemSymbol2[TPNMAX][3] = {
"H",                                                                                 "He",
"Li","Be",                                                  "B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
"Na","Mg",                                                  "Al","Si","P" ,"S" ,"Cl","Ar",
"K" ,"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
"Rb","Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","Id","Sn","Sb","Te","I" ,"Xe",
"Cs","Ba","La","Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
"Fr","Ra","Ac","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","--","--","--","--","--","--",
"Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yt","Lu",
"Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
"Bq"
 };

int chemAtomicNumber[TPNMAX] = { 
 1,                                                                  2,
 3,  4,                                          5,  6,  7,  8,  9, 10,
11, 12,                                         13, 14, 15, 16, 17, 18,
19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,
37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
55, 56, 57, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86,
87, 88, 89,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,
58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101,102,103,
0
};

float chemAtomicRadius[TPNMAX] = {
 37,                                                                 37,
157,112,                                         80, 77, 74, 74, 71, 71,
191,160,                                        143,118,110,103, 99, 99,
235,197,164,147,135,129,137,126,125,125,128,137,153,139,120,116,114,114,
250,215,182,160,147,140,135,134,134,137,144,152,167,158,161,143,133,133,
272,224,188,159,147,141,137,135,136,139,144,155,171,175,182,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
183,182,181,181,180,199,179,176,175,174,173,173,194,172,
180,161,139,140,151,140,  0,  0,  0,  0,  0,  0,  0,  0,
37
};
float chemAtomicRadVdW[TPNMAX] = {
120,                                                                140,
182,153,                                        192,170,155,152,147,154,
227,173,                                        184,210,180,180,175,188,
275,231,211,211,211,211,211,211,211,163,140,139,187,211,185,190,185,202,
303,249,249,249,249,249,249,249,249,163,172,158,193,217,206,206,192,216,
272,224,188,159,147,141,137,135,136,139,144,155,171,175,182,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
183,182,181,181,180,199,179,176,175,174,173,173,194,172,
180,161,139,140,151,140,  0,  0,  0,  0,  0,  0,  0,  0,
120
};
//http://www.educaplus.org/elementos-quimicos/propiedades/radio-atomico.html
//http://www.educaplus.org/elementos-quimicos/propiedades/radio-atomico.html
//datos en pm
// 1 u.a. = 52.918 pm
//
float chemAtomicMass[TPNMAX] = {
  1.0079,                                                                                                                                                  4.0026,
  6.9410,  9.0122,                                                                                           10.8110, 12.0107, 14.0067, 15.9994, 18.9984, 20.1797,
 22.9898, 24.3050,                                                                                           26.9815, 28.0855, 30.9737, 32.0650, 35.4530, 39.9480,
 39.0983, 40.0780, 44.9559, 47.8670, 50.9415, 51.9961, 54.9380, 55.8450, 58.9332, 58.6934, 63.5460, 65.3800, 69.7230, 72.6400, 74.9216, 78.9600, 79.9040, 83.7980,
 85.4678, 87.6200, 88.9059, 91.2240, 92.9064, 95.9600, 98.0000,101.0700,102.9055,106.4200,107.8682,112.4110,114.8180,118.7100,121.7600,127.6000,126.9045,131.2930,
132.91  ,137.33  ,138.91  ,178.49  ,180.95  ,183.84  ,186.21  ,190.23  ,192.22  ,195.08  ,196.97  ,200.59  ,204.38  ,207.2   ,208.98  ,209.    ,210.    ,222.    ,
223.    ,226.    ,227.    ,267.    ,268.    ,271.    ,272.    ,277.    ,276.    ,281.    ,280.    ,285.    ,  0.    ,287.    ,  0.    ,291.    ,  0.    ,   0.   ,
140.12  ,140.91  ,144.24  ,145.    ,150.36  ,151.96  ,157.25  ,158.93  ,162.50  ,164.93  ,167.26  ,168.93  ,173.05  ,174.97,
232.04  ,231.04  ,238.03  ,237.    ,244.    ,243.    ,247.    ,247.    ,251.    ,252.    ,257.    ,258.    ,259.    ,262.,
  0.00
};

// en unidades de masa atomica (Dalton) 1 Da = 1.660538921*10^-27 Kg
//
float getChemProperty(int n, const char *prop){

  int idx,flag=0;
  idx = 0;
  float data;

  data = 0;
  while( idx < TPNMAX && flag == 0 ){
    if( n == chemAtomicNumber[idx] ){
      flag = 1;
    }
    idx++;
  }
  idx--;

  if( !strncmp(prop,"MASS",4) )
    data = chemAtomicMass[idx];
  if( !strncmp(prop,"RADIUS",6) )
    data = chemAtomicRadius[idx]/52.918;
  if( !strncmp(prop,"VANDERWALLS",6))
    data = chemAtomicRadVdW[idx]/52.918;

  return data;
}

int printSymbol( int n ,FILE *file){
  int idx,flag=0;
  idx = 0;
  while( idx < TPNMAX && flag == 0){
    if( n == chemAtomicNumber[idx] ){
      flag = 1;
    }
    idx++;
  }
  idx--;

  fprintf(file,"%3s",chemSymbol2[idx]);
  return 0;
}

int getAtomicSymbol(int n,int m,char *symb){

  int i,idx,flag=0;
  idx = 0;
  while( idx < TPNMAX && flag == 0){
    if( n == chemAtomicNumber[idx] ){
      flag = 1;
    }
    idx++;
  }
  idx--;

  for(i=0;i<m;i++)
    symb[i]='\0';

  strcpy(symb,chemSymbol2[idx]);

  return 0;
}

int getAtomicNumber(const char *symb){
  int i,j,n,z;
  char tmp[LSTRING];
  n = strlen(symb);

  for(i=0,j=0;i<n;i++)
    if( isalpha(symb[i])){
      tmp[j] = toupper(symb[i]);
      j++;
    }
  j++;

  i=0;
  while( i < TPNMAX ){
    if( !strncmp(tmp,chemSymbol[i],j) ){
      z = chemAtomicNumber[i];
      i= 200;
    }
   i++;
  }

  clearString(tmp);
  return z;
}

int printProperties(int n,FILE *out){

  int idx,flag=0;
  idx = 0;

  while( idx < TPNMAX && flag == 0 ){
    if( n == chemAtomicNumber[idx] ){
      flag = 1;
    }
    idx++;
  }
  idx--;

  fprintf(out,"        Symbol :  %s\n",chemSymbol[idx]);
  fprintf(out," Atomic number : %4d\n",n);
  fprintf(out," Number of e-  : %4d\n",n);
  fprintf(out," Number of p-  : %4d\n",n);
  fprintf(out," Number of n0  : %4.0f\n",chemAtomicMass[idx]-n);
  fprintf(out," Atomic mass   : %6.3f Da\n",chemAtomicMass[idx]);
  fprintf(out," Atomic radius : %6.3f pm\n",chemAtomicRadius[idx]);
  return 0;

}

int clearString(char *string){
  int i,n;

  n = strlen(string);
  for(i=0;i<n;i++)
    string[i] = '\0';
}

