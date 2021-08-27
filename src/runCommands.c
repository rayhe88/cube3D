/**
 *  @file   loadRunCommands.c
 *  @brief  File for loading run commands.
 *
 *  @author Raymundo Hernández-Esparza.
 *  @date   August 2021.
 */
#include "runCommands.h"

int checkFile(dataRC *config){
    char path[120];
    char filename[132];
    uid_t uid = getuid();
    struct passwd *pw = getpwuid(uid);

    if( pw == NULL ){
        printf(" Failed to acces home directory!\n");
        exit(EXIT_FAILURE);
    }

    strcpy(path, pw -> pw_dir);

//    printf("path : %s\n", path);
    sprintf(filename,"%s/.cube3drc",path);

//    printf("file : %s\n",filename);

    defaultRunCommands(config);

    if(!access(filename, F_OK) ){
        //printf("The file %s was found\n",filename);
        readFileRC(path, config);

    }else{
        //printf("The file %s was not found\n",filename);
        createFileRC(filename, *config);
    }





    return 0;

}

int createFileRC(char* namerc, dataRC config){

    FILE *filerc;

    openFile(&filerc, namerc, "w+");

    fprintf(filerc, "# File of configuration for Cube3D projec.\n");
    fprintf(filerc, "# This file contains the value of many parameters\n");
    fprintf(filerc, "# for the correct function of the program.\n#\n");
    fprintf(filerc, "#  FIND CRITICAL POINTS SECTION\n");
    fprintf(filerc, "CRIT_TOLERANCE_FUN0 %E\n", config.crit_tolfun0);
    fprintf(filerc, "CRIT_TOLERANCE_FUN  %E\n", config.crit_tolfun);
    fprintf(filerc, "CRIT_TOLERANCE_GRD  %E\n", config.crit_tolgrd);
    fprintf(filerc, "CRIT_TOLERANCE_NRM  %E\n", config.crit_tolnrm);
    fprintf(filerc, "CRIT_MAXITER_STEP1  %d\n", config.crit_maxiter1);
    fprintf(filerc, "CRIT_MAXITER_STEP2  %d\n", config.crit_maxiter2);
    fprintf(filerc, "CRIT_PERCENT_OVERL  %E\n", config.crit_percent);
    fprintf(filerc, "CRIT_TOL_DISTANCE1  %E\n", config.crit_toldist1);
    fprintf(filerc, "CRIT_TOL_DISTANCE2  %E\n", config.crit_toldist2);
    fprintf(filerc, "CRIT_TOL_DISTANCE3  %E\n", config.crit_toldist3);
    fprintf(filerc, "#  GEOMETRICAL DATA SECTION\n");
    fprintf(filerc, "GEOM_MAXITER        %d\n", config.geom_maxiter);
    fprintf(filerc, "GEOM_EPSILON        %E\n", config.geom_eps);
    fprintf(filerc, "GEOM_NUM_POINTS_1D  %d\n", config.geom_npua0);
    fprintf(filerc, "GEOM_NUM_POINTS_2D  %d\n", config.geom_npua1);
    fprintf(filerc, "GEOM_NUM_POINTS_2Dv %d\n", config.geom_npua2);
    fprintf(filerc, "GEOM_NUM_ANGPOINTS  %d\n", config.geom_nang);
    fprintf(filerc, "#  JACOBI DIAGONALIZATION SECTION\n");
    fprintf(filerc, "JACOBI_MAXITER      %d\n", config.jacobi_maxiter);
    fprintf(filerc, "JACOBI_EPSILON      %E\n", config.jacobi_eps);
    fprintf(filerc, "#  NUMERICAL BOND PATH SECTION\n");
    fprintf(filerc, "BOND_PATH_EPSILON   %E\n", config.bpath_eps);
    fprintf(filerc, "BOND_PATH_NSTEP     %d\n", config.bpath_nstep);
    fprintf(filerc, "BOND_PATH_TOL_DIST  %E\n", config.bpath_tol_dist_atm);
    fprintf(filerc, "BOND_PATH_MAXPTS    %d\n", config.bpath_maxpts);
    fprintf(filerc, "#  NUMERICAL GRADIENT LINES SECTION\n");
    fprintf(filerc, "GRAD_LINES_EPSILON  %E\n", config.gradlines_step);
    fprintf(filerc, "# END OF FILE\n");

    fclose(filerc);
    return 0;
}

int readFileRC(char *path, dataRC *config){

    int i, val, idata;
    char namerc[120];
    char nametmp[120];
    char buffer[120];
    double fdata;
    FILE *inp;
    FILE *tmp;


    sprintf(namerc, "%s/.cube3drc", path);

    openFile(&inp, namerc, "r");

    tmpFile(&tmp, ".c3drc", nametmp, "w+");

    while( !feof(inp)){
        fgets(buffer, 100, inp);

        if(buffer[0] != '#'){
            i = 0;
            while( i < 100 && !feof(inp) && buffer[i] != '\n'){
                fprintf(tmp, "%c", buffer[i]);
                i++;
            }
            fprintf(tmp, "\n");
        }
    }

    rewind(inp);
    rewind(tmp);

    while( !feof(tmp)){
        fgets(buffer, 100, tmp);

        if( !strncmp(buffer, "CRIT_TOLERANCE_FUN0", strlen("CRIT_TOLERANCE_FUN0") - 1 )){
            val = sscanf(buffer, "CRIT_TOLERANCE_FUN0 %lf", &fdata);
            if ( val == 1 )
                config -> crit_tolfun0 = fdata;
         }

        if( !strncmp(buffer, "CRIT_TOLERANCE_FUN", strlen("CRIT_TOLERANCE_FUN") - 1 )){
            val = sscanf(buffer, "CRIT_TOLERANCE_FUN %lf", &fdata);
            if ( val == 1 )
                config -> crit_tolfun = fdata;
         }

        if( !strncmp(buffer, "CRIT_TOLERANCE_GRD", strlen("CRIT_TOLERANCE_GRD") - 1 )){
            val = sscanf(buffer, "CRIT_TOLERANCE_GRD %lf", &fdata);
            if ( val == 1 )
                config -> crit_tolgrd = fdata;
         }

        if( !strncmp(buffer, "CRIT_TOLERANCE_NRM", strlen("CRIT_TOLERANCE_NRM") - 1 )){
            val = sscanf(buffer, "CRIT_TOLERANCE_NRM %lf", &fdata);
            if ( val == 1 )
                config -> crit_tolnrm = fdata;
         }

        if( !strncmp(buffer, "CRIT_MAXITER_STEP1", strlen("CRIT_MAXITER_STEP1") - 1 )){
            val = sscanf(buffer, "CRIT_MAXITER_STEP1 %d", &idata);
            if ( val == 1 )
                config -> crit_maxiter1 = idata;
         }

        if( !strncmp(buffer, "CRIT_MAXITER_STEP2", strlen("CRIT_MAXITER_STEP2") - 1 )){
            val = sscanf(buffer, "CRIT_MAXITER_STEP2 %d", &idata);
            if ( val == 1 )
                config -> crit_maxiter2 = idata;
         }

        if( !strncmp(buffer, "CRIT_PERCENT_OVERL", strlen("CRIT_PERCENT_OVERL") - 1 )){
            val = sscanf(buffer, "CRIT_PERCENT_OVERL %lf", &fdata);
            if ( val == 1 )
                config -> crit_percent = fdata;
         }

        if( !strncmp(buffer, "CRIT_TOL_DISTANCE1", strlen("CRIT_TOL_DISTANCE1") - 1 )){
            val = sscanf(buffer, "CRIT_TOL_DISTANCE1 %lf", &fdata);
            if ( val == 1 )
                config -> crit_toldist1 = fdata;
         }

        if( !strncmp(buffer, "CRIT_TOL_DISTANCE2", strlen("CRIT_TOL_DISTANCE2") - 1 )){
            val = sscanf(buffer, "CRIT_TOL_DISTANCE2 %lf", &fdata);
            if ( val == 1 )
                config -> crit_toldist2 = fdata;
         }

        if( !strncmp(buffer, "CRIT_TOL_DISTANCE3", strlen("CRIT_TOL_DISTANCE3") - 1 )){
            val = sscanf(buffer, "CRIT_TOL_DISTANCE3 %lf", &fdata);
            if ( val == 1 )
                config -> crit_toldist3 = fdata;
         }

        if( !strncmp(buffer, "GEOM_MAXITER", strlen("GEOM_MAXITER") - 1 )){
            val = sscanf(buffer, "GEOM_MAXITER %d", &idata);
            if ( val == 1 )
                config -> geom_maxiter = idata;
         }

        if( !strncmp(buffer, "GEOM_EPSILON", strlen("GEOM_EPSILON") - 1 )){
            val = sscanf(buffer, "GEOM_EPSILON %lf", &fdata);
            if ( val == 1 )
                config -> geom_eps = fdata;
         }

        if( !strncmp(buffer, "GEOM_NUM_POINTS_1D", strlen("GEOM_NUM_POINTS_1D") - 1 )){
            val = sscanf(buffer, "GEOM_NUM_POINTS_1D %d", &idata);
            if ( val == 1 )
                config -> geom_npua0 = idata;
         }

        if( !strncmp(buffer, "GEOM_NUM_POINTS_2D", strlen("GEOM_NUM_POINTS_2D") - 1 )){
            val = sscanf(buffer, "GEOM_NUM_POINTS_2D %d", &idata);
            if ( val == 1 )
                config -> geom_npua1 = idata;
         }

        if( !strncmp(buffer, "GEOM_NUM_POINTS_2Dv", strlen("GEOM_NUM_POINTS_2Dv") - 1 )){
            val = sscanf(buffer, "GEOM_NUM_POINTS_2Dv %d", &idata);
            if ( val == 1 )
                config -> geom_npua2 = idata;
         }

        if( !strncmp(buffer, "GEOM_NUM_ANGPOINTS", strlen("GEOM_NUM_ANGPOINTS") - 1 )){
            val = sscanf(buffer, "GEOM_NUM_ANGPOINTS %d", &idata);
            if ( val == 1 )
                config -> geom_nang = idata;
         }

        if( !strncmp(buffer, "JACOBI_MAXITER", strlen("JACOBI_MAXITER") - 1 )){
            val = sscanf(buffer, "JACOBI_MAXITER %d", &idata);
            if ( val == 1 )
                config -> jacobi_maxiter = idata;
         }

        if( !strncmp(buffer, "JACOBI_EPSILON", strlen("JACOBI_EPSILON") - 1 )){
            val = sscanf(buffer, "JACOBI_EPSILON %lf", &fdata);
            if ( val == 1 )
                config -> jacobi_eps = fdata;
         }

        if( !strncmp(buffer, "BOND_PATH_EPSILON", strlen("BOND_PATH_EPSILON") - 1 )){
            val = sscanf(buffer, "BOND_PATH_EPSILON %lf", &fdata);
            if ( val == 1 )
                config -> bpath_eps = fdata;
         }

        if( !strncmp(buffer, "BOND_PATH_NSTEP", strlen("BOND_PATH_NSTEP") - 1 )){
            val = sscanf(buffer, "BOND_PATH_NSTEP %d", &idata);
            if ( val == 1 )
                config -> bpath_nstep = idata;
         }

        if( !strncmp(buffer, "BOND_PATH_TOL_DIST", strlen("BOND_PATH_TOL_DIST") - 1 )){
            val = sscanf(buffer, "BOND_PATH_TOL_DIST %lf", &fdata);
            if ( val == 1 )
                config -> bpath_tol_dist_atm = fdata;
         }

        if( !strncmp(buffer, "BOND_PATH_MAXPTS", strlen("BOND_PATH_MAXPTS") - 1 )){
            val = sscanf(buffer, "BOND_PATH_MAXPTS %d", &idata);
            if ( val == 1 )
                config -> bpath_maxpts = idata;
         }

        if( !strncmp(buffer, "GRAD_LINES_EPSILON", strlen("GRAD_LINES_EPSILON") - 1 )){
            val = sscanf(buffer, "GRAD_LINES_EPSILON %lf", &fdata);
            if ( val == 1 )
                config -> gradlines_step = fdata;
         }


    }



    fclose(inp);
    fclose(tmp);
    remove(nametmp);

    return 0;
}

int defaultRunCommands(dataRC *config){
    // CRIT_TOLERANCE_FUN0
    config -> crit_tolfun0 = 1.E-7;
    // CRIT_TOLERANCE_FUN
    config -> crit_tolfun = 1.E-7;
    // CRIT_TOLERANCE_GRD
    config -> crit_tolgrd = 1.E-5;
    // CRIT_TOLERANCE_NRM
    config -> crit_tolnrm = 100.;
    // CRIT_MAXITER_STEP1
    config -> crit_maxiter1 = 30;
    // CRIT_MAXITER_STEP2
    config -> crit_maxiter2 = 1000;
    // CRIT_PERCENT_OVERL
    config -> crit_percent = 1.5;
    // CRIT_TOL_DISTANCE1
    config -> crit_toldist1 = 1.72E-1;
    // CRIT_TOL_DISTANCE2
    config -> crit_toldist2 = 0.2;
    // CRIT_TOL_DISTANCE3
    config -> crit_toldist3 = 0.7;
    // GEOM_MAXITER
    config -> geom_maxiter = 2000;
    // GEOM_EPSILON
    config -> geom_eps = 0.0125;
    // GEOM_NUM_POINTS_1D
    config -> geom_npua0 = 75;
    // GEOM_NUM_POINTS_2D
    config -> geom_npua1 = 50;
    // GEOM_NUM_POINTS_2Dv
    config -> geom_npua2 = 8;
    // GEOM_NUM_ANGPOINTS
    config -> geom_nang = 36;
    // JACOBI_MAXITER
    config -> jacobi_maxiter = 500;
    // JACOBI_EPSILON
    config -> jacobi_eps = 1.E-8;
    // BOND_PATH_EPSILON
    config -> bpath_eps = 0.005;
    // BOND_PATH_NSTEP
    config -> bpath_nstep = 5;
    // BOND_PATH_TOL_DIST
    config -> bpath_tol_dist_atm= 0.1;
    // BOND_PATH_MAXPTS
    config -> bpath_maxpts = 3000;
    // GRAD_LINES_EPSILON
    config -> gradlines_step = 0.005;

    return 0;
}
