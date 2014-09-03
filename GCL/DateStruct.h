//
//  DateStruct.h
//  GCL
//
//  Created by Jia Yonglei on 14-8-29.
//  Copyright (c) 2014年 Jia Yonglei. All rights reserved.
//

#ifndef GCL_DateStruct_h
#define GCL_DateStruct_h

#include <iostream>
#include <math.h>
#include <float.h>

using namespace std;

#define PI 3.1415926
#define FTOL    1.0E-16
#define XTOL    1.0E-16
#define GTOL    0.0
#define MAXFEV  100
#define EPSFCN  1.0E-16
#define MODE    1
#define FACTOR  100.0
#define EPS     1.0E-16

struct gclPoint
{
    double x, y, z;
    gclPoint()
    {
        x = 0; y = 0; z = 0;
    }
    gclPoint(double x1, double y1, double z1)
    {
        x = x1; y = y1; z = z1;
    }
    void Rotate(double *matrix)
    {
        x = matrix[0] * x + matrix[1] * y + matrix[2] * z;
        y = matrix[3] * x + matrix[4] * y + matrix[5] * z;
        z = matrix[6] * x + matrix[7] * y + matrix[8] * z;
    }
    void operator = (gclPoint pt)
    {
        x = pt.x; y = pt.y; z = pt.z;
    }
};
struct gclLseg{
    gclPoint start, end;
    gclLseg()
    {
        start.x = 0; start.y = 0; start.z = 0;
        end.x = 0; end.y = 0; end.z = 0;
    }
    gclLseg(gclPoint pt1, gclPoint pt2)
    {
        start = pt1;
        end = pt2;
    }
};
struct gclVector{
    double a , b , c;
    gclVector()
    {
        a = 0; b = 0; c = 0;
    }
    gclVector(double a1, double b1, double c1)
    {
        a = a1; b = b1; c = c1;
    }
    gclVector(gclPoint pt1, gclPoint pt2)
    {
        a = pt2.x - pt1.x;
        b = pt2.y - pt1.y;
        c = pt2.z - pt1.z;
    }
    double vLength()
    {
        return sqrt(a * a + b * b + c * c);
    }
    double vPointMult(gclVector v1)
    {
        return v1.a * a + v1. b * c + v1.c * c;
    }
    gclVector vCrossMult(gclVector v1)
    {
        gclVector temp;
        temp.a = b * v1.c - c * v1.b;
        temp.b = c * v1.a - a * v1.c;
        temp.c = a * v1.b - b * v1.a;
        return temp;
    }
    double vAngle(gclVector v1)
    {
        return acos(this->vPointMult(v1) / v1.vLength() / this->vLength());
    }
    void vNorm()
    {
        double len = this->vLength();
        if (len == 0) {
            return;
        }
        a /= len;
        b /= len;
        c /= len;
    }
    void Rotate(double *matrix)
    {
        a = matrix[0] * a + matrix[1] * b + matrix[2] * c;
        b = matrix[3] * a + matrix[4] * b + matrix[5] * c;
        c = matrix[6] * a + matrix[7] * b + matrix[8] * c;
    }
};
struct gclLine{
    double a, b, c;
    gclLine()
    {
        a = 0; b = 0; c = 0;
    }
};

struct gclElipse{
    double a, b;
    gclPoint center;
    gclVector norVector;
    gclVector mainVector;
};

struct gclPlane{
    gclPoint point;
    gclVector norVector;
};

struct gclSphere{
    gclPoint center;
    double r;
    gclSphere()
    {
        center.x = 0; center.y = 0; center.z = 0; r = 0;
    }
};
struct minpackData{
    const gclPoint *points;
    int len;
};

#endif
