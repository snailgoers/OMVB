//
//  geom.cpp
//  GCL
//
//  Created by Jia Yonglei on 14-8-29.
//  Copyright (c) 2014年 Jia Yonglei. All rights reserved.
//

#include "geom.h"

void GetRotateMatrix(gclVector vectorSrc, gclVector vectorDst, double *matrix)
{
    double angle = vectorSrc.vAngle(vectorDst);
    gclVector roteAxis = vectorSrc.vCrossMult(vectorDst);
    roteAxis.vNorm();
    double cost = cos(angle);
    double sint = sin(angle);
    double x = roteAxis.a;
    double y = roteAxis.b;
    double z = roteAxis.c;
    matrix[0] = (1 - cost) * x * x + cost;
    matrix[1] = (1 - cost) * x * y - sint * z;
    matrix[2] = (1 - cost) * x * z + sint * y;
    matrix[3] = (1 - cost) * x * y + sint * z;
    matrix[4] = (1 - cost) * y * y + cost;
    matrix[5] = (1 - cost) * y * z - sint * x;
    matrix[6] = (1 - cost) * x * z - sint * y;
    matrix[7] = (1 - cost) * y * z + sint * x;
    matrix[8] = (1 - cost) * z * z + cost;
}
bool GetInverseMatrix(double *A, double *B, const int len)
{
    double **t = new double *[len];
    for (int i = 0; i < len; i++) {
        t[i] = new double[len];
    }
    for (int i = 0; i < len; i++) {
        for (int j =0; j < len; j++) {
            t[i][j] = A[i * len + j];
            // init B as the unit matrix
            B[i * len + j] = (i ==j) ? (double)1 : 0;
        }
    }
    
    for (int i = 0; i < len; i++) {
        //寻找主元
        double maxValue = t[i][i];
        int k = i;
        for (int j = i + 1; j < len; j++) {
            if (fabs(t[i][j]) > fabs(maxValue)) {
                maxValue = t[i][j];
                k = j;
            }
        }
        // 如果主元所在行不是第i行，进行换行
        if (k != i) {
            for (int j = 0; j < len; j++) {
                double temp = t[i][j];
                t[i][j] = t[k][j];
                t[k][j] = temp;
                // B伴随交换
                temp = B[i];
                B[i] = B[k];
                B[k] = temp;
            }
        }
        // 判断主元是否为零，若是，则矩阵不满足满秩，不存在逆矩阵
        if (t[i][i] == 0) {
            return false;
        }
        // 消去A的第i列除去i行以外的各行原属

        for (int j = 0; j < len; j++) {
            t[i][j] /= t[i][i];
            B[i * len + j] /= t[i][i];
        }
        // 消元
        for (int j = 0; j < len; j++) {
            if (j != i) {
                double temp = t[j][i];
                for (int k = 0; k < len; k++) {
                    t[j][k] -= t[i][k] * temp;
                    B[j * len + k] -= B[i * len + k] * temp;
                }
            }
        }
    }
    
    
    for (int i = 0; i < len; i++) {
        delete []t[i];
    }
    delete []t;
    return true;
}
double geom_dist_point_point(gclPoint pt1, gclPoint pt2)
{
    double dist = sqrt((pt1.x - pt2.x) * (pt1.x - pt2.x) + (pt1.y - pt2.y) * (pt1.y - pt2.y) + (pt1.z - pt2.z) * (pt1.z - pt2.z));
    return dist;
}

double geom_dist_point_lseg(gclPoint pt, gclLseg lseg)
{
    gclVector v(pt, lseg.start);
    gclVector vLine(lseg.start, lseg.end);
    gclVector temp = v.vCrossMult(vLine);
    
    gclVector v1(lseg.start, pt);
    gclVector v2(lseg.start, lseg.end);
    gclVector v3(lseg.end, pt);
    gclVector v4(lseg.end, lseg.start);
    
    if (v1.vAngle(v2) > PI / 4) {
        return v1.vLength();
    }
    else if (v3.vAngle(v4) > PI / 4)
    {
        return v3.vLength();
    }
    else
    {
        return temp.vLength() / vLine.vLength();
    }
}
double geom_dist_point_plane(gclPlane plane, gclPoint point)
{
    gclVector v = plane.norVector;
    gclPoint pt = plane.point;
    return fabs(v.a * (point.x - pt.x) + v.b * (point.y - pt.y) + v.c * (point.z - pt.z)) / sqrt(v.a * v.a + v.b * v.b + v.c * v.c);
}
// 第一象限中点到标准椭圆距离
double DistancePointEllipseSpecial(const double e[2], const double y[2], double x[2])
{
    double distance;
    if (y[1] > 0)
    {
        if (y[0] > 0)
        {
            // bisect to compute the root of F(t) for t >= -e1*d1
            double esqr[2] = {e[0] * e[0], e[1] * e[1]};
            double ey[2] = {e[0] * y[0], e[1] * y[1]};
            double t0 = -esqr[1] + ey[1];
            double t1 = -esqr[1] + sqrtf(ey[0]*ey[0] + ey[1]*ey[1]);
            double t = t0;
            for (int i = 0; ; ++i)
            {
                t = 0.5*(t0 + t1);
                if (t == t0 || t == t1)
                {
                    break;
                }
                double r[2] = { ey[0]/(t + esqr[0]), ey[1]/(t + esqr[1]) };
                double f = r[0]*r[0] + r[1]*r[1] - 1;
                if (f > 0)
                {
                    t0 = t;
                }
                else if (f < 0)
                {
                    t1 = t;
                }
                else
                {
                    break;
                }
            }
            x[0] = esqr[0]*y[0]/(t + esqr[0]);
            x[1] = esqr[1]*y[1]/(t + esqr[1]);
            double d[2] = { x[0] - y[0], x[1] - y[1] };
            distance = sqrtf(d[0]*d[0] + d[1]*d[1]);
        }
        else // y0 == 0
        {
            x[0] = 0;
            x[1] = e[1];
            distance = fabs(y[1] - e[1]);
        }
    }
    else // y1 == 0
    {
        double denom0 = e[0]*e[0] - e[1]*e[1];
        double e0y0 = e[0]*y[0];
        if (e0y0 < denom0)
        {
            // y0 is inside the subinterval.
            double x0de0 = e0y0/denom0;
            double x0de0sqr = x0de0*x0de0;
            x[0] = e[0]*x0de0;
            x[1] = e[1]*sqrtf(fabs(1 - x0de0sqr));
            double d0 = x[0] - y[0];
            distance = sqrtf(d0*d0 + x[1]*x[1]);
        }
        else
        {
            // y0 is outside the subinterval. The closest ellipse point has
            // x1 == 0 and is on the domain-boundary interval (x0/e0)^2 = 1.
            x[0] = e[0];
            x[1] = 0;
            distance = fabs(y[0] - e[0]);
        }
    }
    return distance;
}

//任意象限点到椭圆距离
double DistancePointEllipse (const double e[2], const double y[2], double x[2])
{
    // Determine reflections for y to the first quadrant.
    bool reflect[2];
    int i, j;
    for (i = 0; i < 2; ++i)
    {
        reflect[i] = y[i] < 0;
    }
    // Determine the axis order for decreasing extents.
    int permute[2];
    if (e[0] < e[1])
    {
        permute[0] = 1; permute[1] = 0;
    }
    else
    {
        permute[0] = 0; permute[1] = 1;
    }
    int invpermute[2];
    for (i = 0; i < 2; ++i)
    {
        invpermute[permute[i]] = i;
    }
    double locE[2], locY[2];
    for (i = 0; i < 2; ++i)
    {
        j = permute[i];
        locE[i] = e[j];
        locY[i] = y[j];
        if (reflect[j])
        {
            locY[i] = -locY[i];
        }
    }
    double locX[2];
    double distance = DistancePointEllipseSpecial(locE, locY, locX);
    // Restore the axis order and reflections.
    for (i = 0; i < 2; ++i)
    {
        j = invpermute[i];
        if (reflect[i])
        {
            locX[j] = -locX[j];
        }
        x[i] = locX[j];
    }
    return distance;
}
double geom_dist_point_elipseHoriranttol(gclElipse elipse, gclPoint point)
{
    // 计算点到法向量为［0，0，1］得椭圆得距离
    double yin[2], xout[2];
    double ein2[2] = {elipse.a, elipse.b};
    
    gclVector v1(1, 0, 0);
    double theta = acos(v1.vPointMult(elipse.mainVector) / elipse.mainVector.vLength());
    double cost = cos(theta);
    double sint = sin(theta);
    gclPoint pt;
    pt.x -= elipse.center.x;
    pt.y -= elipse.center.y;
    pt.z -= elipse.center.z;
    yin[0] = cost * pt.x + sint * pt.y;
    yin[1] = -sint * pt.x + cost * pt.y;
    double distance = DistancePointEllipse(ein2, yin, xout);
    
    pt.x = cost * xout[0] - sint * xout[1] + elipse.center.x;
    pt.y = sint * xout[0] + cost * xout[1] + elipse.center.y;
    pt.z = 0;
    distance = sqrt((pt.x - point.x) * (pt.x - point.x) + (pt.y - point.y) * (pt.y - point.y));
    
    return distance;
}

double geom_dist_point_elipse(gclElipse elipse, gclPoint point)
{
    // dist = sqrt(disV ^ 2 + disH ^2)
    gclPlane plane;
    plane.point = elipse.center;
    plane.norVector = elipse.norVector;
    double distHorirantol = geom_dist_point_plane(plane, point);
    
    // 二维 点到椭圆距离
    double *matrix = new double[9];
    gclVector vsrc = elipse.norVector;
    vsrc. vNorm();
    gclVector vdst(0, 0, 1);
    GetRotateMatrix(vsrc, vdst, matrix);
    
    point.Rotate(matrix);
    gclVector mvector = elipse.mainVector;
    mvector.Rotate(matrix);
    
    gclPoint pcenter = elipse.center;
    pcenter.Rotate(matrix);
    
    gclElipse elipseV;
    elipseV.a = elipse.a;
    elipseV.b = elipse.b;
    elipseV.norVector = vdst;
    elipseV.mainVector = mvector;
    
    double distVertical = geom_dist_point_elipseHoriranttol(elipseV, point);
    
    delete []matrix;
    return sqrt(distHorirantol * distHorirantol + distVertical * distVertical);
}

/// geom fitting

// plane fitting
bool geom_fit_plane_nonVertical(gclPlane &plane, gclPoint *points, int len)
{
    if (points == NULL || len < 3) {
        return false;
    }
    double sxx = 0, sxy = 0, sx = 0, syy = 0, sy = 0, sxz = 0, syz = 0, sz = 0;
    for (int i = 0; i < len; i++) {
        sxx += points[i].x * points[i].x;
        sxy += points[i].x * points[i].y;
        sx += points[i].x;
        syy += points[i].y * points[i].y;
        sy += points[i].y;
        sxz += points[i].x * points[i].z;
        syz += points[i].y * points[i].z;
        sz += points[i].z;
    }
    double AA[] = {sxx, sxy, sx, sxy, syy, sy, sx, sy, double(len)};
    double invAA[9];
    if (GetInverseMatrix(AA, invAA, 3)) {
        return false;
    }
    double x[] = {sxz, syz, sz};
    double a[3];
    for (int i = 0; i < 3; i++) {
        double sumtemp = 0.0;
        for (int j = 0; j < 3; j++) {
            sumtemp += invAA[i * 3 + j] * x[j];
        }
        a[i] = sumtemp;
    }
    plane.norVector.a =  -a[0];
    plane.norVector.b = -a[1];
    plane.norVector.c = 1;
    plane.point.x = sx / len;
    plane.point.y = sy / len;
    plane.point.z = a[0] * plane.point.x + a[1] * plane.point.y + a[2];
    
    return true;
}

int Minpack_CostFun_Plane(void *pspo, int m, int n, const __cminpack_real__ *x, __cminpack_real__ *fvec, int iflag)
{
    const gclPoint *points = ((minpackData *)pspo)->points;
    int len = ((minpackData*)pspo)->len;
    gclPlane plane;
    plane.norVector.a = x[0];
    plane.norVector.b = x[1];
    plane.norVector.c = x[2];
    plane.point.z = x[3];
    for (int i = 0; i < len; i++) {
        fvec[i] = geom_dist_point_plane(plane, points[i]);
    }
    return 0;
}
bool geom_fit_plane(gclPlane plane, gclPoint *points, int len)
{
    if (points == NULL || len < 3) {
        return  false;
    }
    // init the plane
    double *dist = new double[len * len];
    memset(dist, 0, sizeof(double) * len * len);
    double maxDis = -DBL_MAX;
    int index0, index1, index2;
    for (int i = 0; i < len - 1; i++) {
        for (int j = i + 1; j < len; j++) {
            dist[i * len + j] =fabs(points[i].x - points[j].x) + fabs(points[i].y - points[i].y) + fabs(points[i].z - points[j].z);
            if (dist[i * len + j] > maxDis) {
                maxDis = dist[i * len +j];
                index0 = i;
                index1 = j;
            }
            dist[j * len + i] = dist[i * len + j];
        }
    }
    maxDis = -DBL_MAX;
    gclPoint pt;
    for (int i = 0; i < len; i++) {
        double distemp  = dist[index0 * len + i] + dist[index1 * len + i];
        if (distemp > maxDis) {
            maxDis = distemp;
            index2 = i;
        }
        pt.x += points[i].x;
        pt.y += points[i].y;
        pt.z += points[i].z;
    }
    pt.x /= len;
    pt.y /= len;
    pt.z /= len;
    
    delete []dist;
    
    plane.point = pt;
    gclVector v1(points[index0], points[index1]);
    gclVector v2(points[index0], points[index2]);
    plane.norVector = v1.vCrossMult(v2);
    plane.norVector.vNorm();
    
    if (plane.norVector.c > EPS) {
        geom_fit_plane_nonVertical(plane, points, len);
        return true;
    }

    // fitting the plane
    int m = len;
    int n = 4;
    double *fvec = new double[m];
    double *x = new double[n];
    double *fjac = new double[m * n];
    double *wa4 = new double[m];
    double *diag = new double[n];
    double *qtf = new double[n];
    double *wa1 = new double[n];
    double *wa2 = new double[n];
    double *wa3 = new double[n];
    int *ipvt = new int[n];
    double ftol = FTOL;
    double xtol = XTOL;
    double gtol = GTOL;
    int maxfev = MAXFEV * n;
    double epsfcn = EPSFCN;
    int mode = MODE;
    double factor = FACTOR;
    int nprint = 0;
    int info = 0;
    int nfev;
    int ldfjac = m;
    
    // init unknown
    minpackData data;
    data.len = len;
    data.points = points;
    
    x[0] = plane.norVector.a;
    x[1] = plane.norVector.b;
    x[2] = plane.norVector.c;
    x[4] = plane.point.z;
    
    info = __cminpack_func__(lmdif)(Minpack_CostFun_Plane, &data, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor, nprint, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4);
    
    plane.norVector.a = x[0];
    plane.norVector.b = x[1];
    plane.norVector.c = x[2];
    plane.point.z = x[3];

    delete[] fvec;
    delete []x;
    delete []fjac;
    delete []diag;
    delete []ipvt;
    delete []qtf;
    delete []wa1;
    delete []wa2;
    delete []wa3;
    delete []wa4;
    
    return true;
}

// fit the lseg
int Minpack_CostFun_Lseg(void *pspo, int m, int n, const __cminpack_real__ *x, __cminpack_real__ *fvec, int iflag)
{
    const gclPoint *points = ((minpackData *)pspo)->points;
    int len = ((minpackData *)pspo)->len;
    
    gclLseg lseg;
    lseg.start.x = x[0];
    lseg.start.y = x[1];
    lseg.start.z = x[2];
    lseg.end.x = x[3];
    lseg.end.y = x[4];
    lseg.end.z = x[5];
    for (int i = 0; i < len; i++) {
        fvec[i] = geom_dist_point_lseg(points[i], lseg);
    }
    return 0;
}

bool geom_fit_lseg_3D(gclLseg &lseg, gclPoint *points, int len)
{
    if (points == NULL || len < 2) {
        return false;
    }
    // init the lseg
    int index1 = 0;
    int index2 = len - 1;
    double disMax = -DBL_MAX;
    for (int i = 0; i < len - 1; i++) {
        for (int j = i + 1; j < len; j++) {
            gclPoint pt;
            pt.x = points[i].x - points[j].x;
            pt.y = points[i].y - points[j].y;
            pt.z = points[i].z - points[j].z;
            double dis = sqrt(pt.x * pt.x + pt.y * pt.y + pt.z * pt.z);
            if (dis > disMax) {
                disMax = dis;
                index1 = i;
                index2 = j;
            }
        }
    }
    lseg.start = points[index1];
    lseg.end   = points[index2];
    
    int m = len;
    int n = 6; // number of unknown
    double *fvec = new double[m];
    double *x = new double[n];
    double *fjac = new double[m * n];
    double *wa4 = new double[m];
    double *diag = new double[n];
    double *qtf = new double[n];
    double *wa1 = new double[n];
    double *wa2 = new double[n];
    double *wa3 = new double[n];
    int *ipvt = new int[n];
    double ftol = FTOL;
    double xtol = XTOL;
    double gtol = GTOL;
    int maxfev = MAXFEV * n;
    double epsfcn = EPSFCN;
    int mode = MODE;
    double factor = FACTOR;
    int nprint = 0;
    int info = 0;
    int nfev;
    int ldfjac = m;
    
    // init unknown
    
    minpackData data;
    data.len = len;
    data.points = points;
    
    x[0] = lseg.start.x;
    x[1] = lseg.start.y;
    x[2] = lseg.start.z;
    x[3] = lseg.end.x - lseg.start.x;
    x[4] = lseg.end.y - lseg.start.y;
    x[5] = lseg.end.z - lseg.start.z;
    
    
    info = __cminpack_func__(lmdif)(Minpack_CostFun_Lseg, &data, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor, nprint, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4);
    
    lseg.start.x = x[0];
    lseg.start.y = x[1];
    lseg.start.z = x[2];
    lseg.end.x = x[0] + x[3];
    lseg.end.y = x[1] + x[4];
    lseg.end.z = x[2] + x[5];
    
    delete[] fvec;
    delete []x;
    delete []fjac;
    delete []diag;
    delete []ipvt;
    delete []qtf;
    delete []wa1;
    delete []wa2;
    delete []wa3;
    delete []wa4;
    
    return true;
}

// elipse fitting
int Minpack_CostFun_Elipse(void *pspo, int m, int n, const __cminpack_real__ *x, __cminpack_real__ *fvec, int iflag)
{
    const gclPoint *points = ((minpackData*)pspo)->points;
    int len =((minpackData*)pspo)->len;
    gclElipse elipse;
    elipse.a = x[0];
    elipse.b = x[1];
    elipse.center.x = x[2];
    elipse.center.y = x[3];
    elipse.center.z = x[4];
    elipse.norVector.a = x[5];
    elipse.norVector.b = x[6];
    elipse.norVector.c = x[7];
    elipse.mainVector.a = x[8];
    elipse.mainVector.b = x[9];
    elipse.mainVector.c = x[10];
    
    for (int i = 0; i < len; i++) {
        fvec[i] = geom_dist_point_elipse(elipse, points[i]);
    }
    return 0;
}

bool geom_fit_elipse(gclElipse &elipse, gclPoint *points, int len)
{
    if (points == NULL || len < 11) {
        return false;
    }
    // init the elipse
    gclPlane plane;
    geom_fit_plane(plane, points, len);
    gclPoint pt, pta;
    for (int i = 0; i < len; i++) {
        pt.x += points[i].x;
        pt.y += points[i].y;
        pt.z += points[i].z;
    }
    pt.x /= len;
    pt.y /= len;
    pt.z /= len;
    
    double maxDis = -DBL_MAX;
    double minDis = DBL_MAX;
    for (int i = 0; i < len; i++) {
        double dis = fabs(pt.x - points[i].x) + fabs(pt.y - points[i].y) + fabs(pt.z - points[i].z);
        if (dis > maxDis) {
            maxDis = dis;
            pta.x = points[i].x;
            pta.y = points[i].y;
            pta.z = points[i].z;
        }
        if (dis < minDis) {
            minDis = dis;
        }
    }
    gclVector mainVector(pt, pta);
    elipse.a = maxDis;
    elipse.b = minDis;
    elipse.center = pt;
    elipse.norVector = plane.norVector;
    elipse.mainVector = mainVector;
    
    
    int m = len;
    int n = 11; // number of unknown
    double *fvec = new double[m];
    double *x = new double[n];
    double *fjac = new double[m * n];
    double *wa4 = new double[m];
    double *diag = new double[n];
    double *qtf = new double[n];
    double *wa1 = new double[n];
    double *wa2 = new double[n];
    double *wa3 = new double[n];
    int *ipvt = new int[n];
    double ftol = FTOL;
    double xtol = XTOL;
    double gtol = GTOL;
    int maxfev = MAXFEV * n;
    double epsfcn = EPSFCN;
    int mode = MODE;
    double factor = FACTOR;
    int nprint = 0;
    int info = 0;
    int nfev;
    int ldfjac = m;
    
    // init unknown
    
    minpackData data;
    data.len = len;
    data.points = points;
    
    x[0] = elipse.a;
    x[1] = elipse.b;
    x[2] = elipse.center.x;
    x[3] = elipse.center.y;
    x[4] = elipse.center.z;
    x[5] = elipse.norVector.a;
    x[6] = elipse.norVector.b;
    x[7] = elipse.norVector.c;
    x[8] = elipse.mainVector.a;
    x[9] = elipse.mainVector.b;
    x[10] = elipse.mainVector.c;
    
    
    info = __cminpack_func__(lmdif)(Minpack_CostFun_Elipse, &data, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor, nprint, &nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4);
    
    elipse.a = x[0];
    elipse.b = x[1];
    elipse.center.x = x[2];
    elipse.center.y = x[3];
    elipse.center.z = x[4];
    elipse.norVector.a = x[5];
    elipse.norVector.b = x[6];
    elipse.norVector.c = x[7];
    elipse.mainVector.a = x[8];
    elipse.mainVector.b = x[9];
    elipse.mainVector.c = x[10];
    
    delete[] fvec;
    delete []x;
    delete []fjac;
    delete []diag;
    delete []ipvt;
    delete []qtf;
    delete []wa1;
    delete []wa2;
    delete []wa3;
    delete []wa4;
    
    return true;
}




