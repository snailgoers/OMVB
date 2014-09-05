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
    // test for opencv
    Mat src = imread("../../../GCL/testdata/1.tif", 0);
    imshow("src", src);
    waitKey();
    
    
    
    
    // read the data
    int len = 1296;
    gclPoint *points = new gclPoint[len];
    ifstream file("../../../GCL/testdata/sphereFit.txt");
    for (int i = 0; i < len; i++) {
        file >> points[i].x >> points[i].y  >> points[i].z;
    }
    file.close();
    
    for (int i = 0; i < len; i++) {
        printf("%f %f %f\n", points[i].x, points[i].y, points[i].z);
    }
    
    // test for fitting lseg
//    gclLseg lseg;
//    if (geom_fit_lseg_3D(lseg, points, len)) {
//        printf("%f %f %f\n%f %f %f\n", lseg.start.x, lseg.start.y, lseg.start.z, lseg.end.x, lseg.end.y, lseg.end.z);
//    }
    
    // test for elipse fitting
//    gclElipse elipse;
//    if (geom_fit_elipse(elipse, points, len)) {
//        printf("a= %f b= %f\np.x = %f p.y = %f p.z = %f\nnv.a = %f nv.b = %f nv.c = %f\nmv.a = %f mv.b = %f mv.c = %f\n", elipse.a, elipse.b, elipse.center.x, elipse.center.y, elipse.center.z, elipse.norVector.a, elipse.norVector.b, elipse.norVector.c, elipse.mainVector.a, elipse.mainVector.b, elipse.mainVector.c);
//    }
    
    // test for sphere fitting
    // true value:center =[-10, 5, 20],R = 30;
    gclSphere sphere;
    if (geom_fit_sphere(sphere, points, len)) {
        printf("%f %f %f\n%f\n", sphere.center.x, sphere.center.y, sphere.center.z, sphere.r);
    }
    
    
    delete []points;
    return 0;
}
