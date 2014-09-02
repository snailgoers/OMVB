//
//  geom.h
//  GCL
//
//  Created by Jia Yonglei on 14-8-29.
//  Copyright (c) 2014年 Jia Yonglei. All rights reserved.
//

#ifndef __GCL__geom__
#define __GCL__geom__

#define __cminpack_double__



#include "DateStruct.h"
#include "cminpack134/cminpack.h"

// geom base
bool GetInverseMatrix(double *A, double *B, const int len);
void GetRotateMatrix(gclVector vectorSrc, gclVector vectorDst, double *matrix);


// geom distance
double geom_dist_point_point(gclPoint pt1, gclPoint pt2);
double geom_dist_point_lseg(gclPoint pt, gclLseg lseg);
double geom_dist_point_elipse(gclElipse elipse, gclPoint point);
double geom_dist_point_plane(gclPlane plane, gclPoint point);


//geom fitting
bool geom_fit_lseg_3D(gclLseg &lseg, gclPoint *points, int len);
bool geom_fit_elipse(gclElipse &elipse, gclPoint *points, int len);
bool geom_fit_plane(gclPlane plane, gclPoint *points, int len);


#endif /* defined(__GCL__geom__) */
