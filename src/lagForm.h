/** @file   lagForm.h
 *  @brief  
 *  @author Raymundo Hernández-Esparza.
 *  @date   August 2018.
 */
#include <stdio.h>
#include <stdlib.h>
// val[0] campo en r0
//
// val[1] derivada  x evaluada en r0
// val[2] derivada  y evaluada en r0
// val[3] derivada  z evaluada en r0
//
// val[4] derivada xx evaluada en r0
// val[5] derivada yy evaluada en r0
// val[6] derivada zz evaluada en r0
//
// val[7] derivada xy evaluada en r0
// val[8] derivada xz evaluada en r0
// val[9] derivada yz evaluada en r0

#ifndef _LAGFORM_H_
 #define _LAGFORM_H_
 void setVal(double *val);
 int Ind(int idx,int idy, int idz,int n);
/***************************************************************/
 int gradPol01(double,double,double,double *,double *);
 int lapPol01 (double,double,double,double *,double *);
 int hessPol01(double,double,double,double *,double *);
/***************************************************************/
 int gradPol02(double,double,double,double *,double *);
 int lapPol02 (double,double,double,double *,double *);
 int hessPol02(double,double,double,double *,double *);
/***************************************************************/
 int gradPol03(double,double,double,double *,double *);
 int lapPol03 (double,double,double,double *,double *);
 int hessPol03(double,double,double,double *,double *);
/***************************************************************/
 int gradPol04(double,double,double,double *,double *);
 int lapPol04 (double,double,double,double *,double *);
 int hessPol04(double,double,double,double *,double *);
/***************************************************************/
 int gradPol05(double,double,double,double *,double *);
 int lapPol05 (double,double,double,double *,double *);
 int hessPol05(double,double,double,double *,double *);
/***************************************************************/
 int gradPol06(double,double,double,double *,double *);
 int lapPol06 (double,double,double,double *,double *);
 int hessPol06(double,double,double,double *,double *);
#endif
