//
//  testFile.cpp
//  GCL
//
//  Created by Jia Yonglei on 14-8-29.
//  Copyright (c) 2014年 Jia Yonglei. All rights reserved.
//

#include "testFile.h"

int main(void)
{
    // read the data
    int len = 11;
    gclPoint *points = new gclPoint[len];
    ifstream file("../../../GCL/testdata/line3Dfitting2.txt");
    for (int i = 0; i < len; i++) {
        file >> points[i].x >> points[i].y  >> points[i].z;
    }
    file.close();
    
    for (int i = 0; i < len; i++) {
        printf("%f %f %f\n", points[i].x, points[i].y, points[i].z);
    }
    
    // test for fitting lseg
    gclLseg lseg;
    if (geom_fit_lseg_3D(lseg, points, len)) {
        printf("%f %f %f\n%f %f %f\n", lseg.start.x, lseg.start.y, lseg.start.z, lseg.end.x, lseg.end.y, lseg.end.z);
    }
    
    
    
    delete []points;
    return 0;
}
