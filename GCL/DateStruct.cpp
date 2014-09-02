//
//  DateStruct.cpp
//  GCL
//
//  Created by Jia Yonglei on 14-8-29.
//  Copyright (c) 2014年 Jia Yonglei. All rights reserved.
//

#include "DateStruct.h"

gclPoint::gclPoint()
{
    x = 0; y = 0; z = 0;
}
gclPoint::gclPoint(double x1, double y1, double z1)
{
    x = x1; y = y1; z = z1;
}
double gclPoint::getX()
{
    return x;
}
double gclPoint::getY()
{
    return y;
}
double gclPoint::getZ()
{
    return z;
}
void gclPoint::setX(double x1)
{
    x = x1;
}
void gclPoint::setY(double y1)
{
    y = y1;
}
void gclPoint::setZ(double z1)
{
    z = z1;
}