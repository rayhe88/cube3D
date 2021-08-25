/**
 * @file    struct.h
 * @brief   structure for storing data cube.
 * @authors Raymundo Hernández-Esparza.
 * @date    August 2018.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _STRUCT_H_
#define _STRUCT_H_

#define NOT 0 /**< Macro for NOT.                  */
#define YES 1 /**< Macro for YES.                  */

#define RED 0 /**< Macro for Reduced gradient.     */
#define GRA 1 /**< Macro for gradient.             */
#define LAP 2 /**< Macro for laplacian.            */
#define KIN 3 /**< Macro for Kinetic (Abramov) E.  */
#define VIR 4 /**< Macro for Potencial Energy.     */
#define KEW 5 /**< Macro for Kin Ene Den Weissaker.*/
#define NCI 6 /**< Macro for NCI index.            */
#define CRI 7 /**< Macro for critical points.      */
#define VOI 8 /**< Macro for voids in density.     */
#define REP 9 /**< Macro for replication.          */

#define LIN 100 /**< Macro for Geometry Line.        */
#define PLA 200 /**< Macro for Geometry Plane.       */

#define PLA_F 210 /**< Macro for Geom Plane and field. */
#define PLA_V 220 /**< Macro for Geom Plane and vector.*/
#define PLA_S 230 /**< Macro for Geom Plane and Stream.*/

#define DELTA 0.001 /**< Macro for delta in the log-func.*/

#define B2A 0.529177249 /**< Macro to convert Bohr to Angstr.*/

typedef struct {

    unsigned int natm;   /**< Number of atoms in the cube.    */
    unsigned int npt;    /**< Number of total points.         */
    unsigned int pts[3]; /**< Array for points in x, y and z. */
    double min[3];       /**< Array for \f$ \vec{r}_0 \f$.    */
    double max[3];       /**< Array for \f$ \vec{r}_{n-1} \f$.*/
    double hvec[3];      /**< Array for steps in hx,hy,hz.    */
    double mvec[9];      /**< Matrix of transformation.       */
    int *zatm;           /**< Pointer to atomic numbers.      */
    double *field;       /**< Pointer to field data           */
    double *coor;        /**< Pointer to coordinates          */

} dataCube; /**< Structure for cube information. */

typedef struct {

    int pol;  /**< Polynomial's degree.            */
    int pbc;  /**< Periodic boundary conditions.   */
    int orth; /**< Orthogonality of the system.    */
    int task; /**< Task of the running.            */
    int izq;  /**< Points to the left.             */
    int der;  /**< Points to the right.            */
    int size; /**< Size of the array for the values
                   of the field                    */
    int geoTask;
    int geoProp;
    int rep[3];  /**< Replica in the 3 direcctions.   */
    int geom[3]; /**< Atoms for geometry lin or plan. */
    double la2;  /**< Cutoff for the lambda2 value.   */
    double rgd;  /**< Cutoff for the reduced gradient.*/
    double vac;  /**< Cutoff for density in voids.    */

} dataRun; /**< Structure for run information.  */

typedef struct {
    double tolfun0;
    double tolfun;
    double tolgrd;
    double tolnrm;
    int maxiter1;
    int maxiter2;
    double percent;
    double toldist1;
    double toldist2;
    double toldist3;
    int npua0;
    int npua1;
    int npua2;
    int nang;
    int maxiter;
    double eps;
    int kindCenter;
    char gnuplot_path;
    char gnuplot_path_wx;
    double jacobi_eps;
    int jacobi_maxiter;
    double math_tol;
    int math_nmax;
    double bpath_eps;
    int bpath_nstep;
    double bpath_tol_dist_atm;
    int maxpts;
    double gradlines_step;
} RunCommands;

#endif
