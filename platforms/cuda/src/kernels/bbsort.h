/*
 * Authored by: Chen, Shifu
 * 
 * Email: chen@gmtk.org
 *
 * Website: http://www.gmtk.org/gsort
 *
 * The code is distributed under BSD license, you are allowed to use, modify or sell this code, but a statement is required if you used this code any where.
 * 
 */
#ifndef _BBSORT_H_
#define _BBSORT_H_
#include "windowsExportCuda.h"

#define BLOCK_SIZE 512

#define DISORDERLY 0
#define NEARLY_SORTED 1
#define AUTO_EVALUATE 2

template <typename T>
void OPENMMCUDA_EXPORT bbSort(T* dData,int number,int listOrder=AUTO_EVALUATE);

#endif // _BBSORT_H_
