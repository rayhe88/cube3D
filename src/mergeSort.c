#include "basins2.h"
#include "mergeSort.h"

void mergeSort(dataCells2 *array, int l, int r ){
    if( l < r){
        int m = l + (r-l)/2;

        mergeSort(array,l,m);
        mergeSort(array,m+1,r);

        merge(array,l,m,r);
    }
}

void merge(dataCells2 *array,int l, int m, int r){
    int i,j,k;
    int n1 = m - l + 1;
    int n2 = r - m;

    //temporal arrays
    dataCells2 *arrayLeft;
    dataCells2 *arrayRight;
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
void mergeSortbyAtr(dataCells2 *array, int l, int r ){
    if( l < r){
        int m = l + (r-l)/2;

        mergeSortbyAtr(array,l,m);
        mergeSortbyAtr(array,m+1,r);

        mergebyAtr(array,l,m,r);
    }
}

void mergebyAtr(dataCells2 *array,int l, int m, int r){
    int i,j,k;
    int n1 = m - l + 1;
    int n2 = r - m;

    //temporal arrays
    dataCells2 *arrayLeft;
    dataCells2 *arrayRight;
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
        if( arrayLeft[i].attr <= arrayRight[j].attr){
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
