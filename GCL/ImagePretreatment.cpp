//
//  ImagePretreatment.cpp
//  GCL
//
//  Created by Jia Yonglei on 14-9-9.
//  Copyright (c) 2014年 Jia Yonglei. All rights reserved.
//

#include "ImagePretreatment.h"

/// get colormap

void GetColorMap(int level, double *map)
{
    int index = 0;
    int nScale[] = { int(level * 0.125), int(level * 0.375), int(level * 0.625), int(level * 0.875), level};
    double step = 1 / (level / 4.0);
    double iniValue = (nScale[0] + 1) * step;
    for (int i = 1; i < nScale[0]; i++)
    {
        map[index ++] = iniValue + (i - 1) * step;
        map[index ++] = 0;
        map[index ++] = 0;
    }
    for (int i = nScale[0]; i < nScale[1]; i++)
    {
        map[index ++] = 1;
        map[index ++] = (i - nScale[0]) * step;
        map[index ++] = 0;
    }
    for (int i = nScale[1]; i < nScale[2]; i++)
    {
        map[index ++] = 1 - (i - nScale[1]) * step;
        map[index ++] = 1;
        map[index ++] = (i - nScale[1]) * step;
    }
    for (int i = nScale[2]; i < nScale[3]; i++)
    {
        map[index ++] = 0;
        map[index ++] = 1 - (i - nScale[2]) * step;
        map[index ++] = 1;
    }
    for (int i = nScale[3]; i <= nScale[4]; i++)
    {
        map[index ++] = 0;
        map[index ++] = 0;
        map[index ++] = 1 - (i - nScale[3]) * step;
    }
}

/// 归一化[0,1]
void Normalization(double *dataIn, double *dateOut, int size)
{
    double maxV = -DBL_MAX;
    double minV = DBL_MAX;
    for (int i = 0; i < size; i++)
    {
        if (dataIn[i] > maxV)
        {
            maxV = dataIn[i];
        }
        if (dataIn[i] < minV)
        {
            minV = dataIn[i];
        }
    }
    for (int i = 0; i < size; i++)
    {
        dateOut[i] = (dataIn[i] - minV) / (maxV - minV);
    }
}

/// 卷积操作
bool Conv2(uchar *dataIn, double *dataOut, int width, int height, double *kernel, int len)
{
    if (dataIn == NULL)
    {
        return false;
    }
    for (int i = 0; i < width * height; i++)
    {
        dataOut[i] = dataIn[i] / 255.0;
    }
    int slen = sqrt(double(len));
    int halfLen = slen / 2;
    for (int i = halfLen; i < height - halfLen; i++)
    {
        for (int j = halfLen; j < width - halfLen; j++)
        {
            double temp = 0;
            for (int m = -halfLen; m <= halfLen; m++)
            {
                for (int n = -halfLen; n <= halfLen; n++)
                {
                    temp += dataIn[(i + m) * width + j + n] / 255.0 * kernel[(m + halfLen) * slen + n + halfLen];
                }
            }
            dataOut[i * width + j] = temp;
        }
    }
    return true;
}
/// 彩色空间转换

bool RGB2HSI(uchar *imgIn, int nchannels, double *Hvector, double *Svector, double *Ivector, int width, int height)
{
    if (imgIn == NULL || nchannels != 3)
    {
        return false;
    }
    int index = 0;
    for (int i = 0; i < width * height * nchannels; i = i + 3)
    {
        double b = imgIn[i + 0] / 255.0;
        double g = imgIn[i + 1] / 255.0;
        double r = imgIn[i + 2] / 255.0;
        double num = 0.5 * (r - g + r - b);
        double den = sqrt((r - g) * (r - g) + (r - b) * (g - b));
        double theta = acos(num / (den + EPSMIN));
        double hh = theta;
        if (b > g)
        {
            hh = 2 * PI - hh;
        }
        hh /= 2 * PI;
        num = gclMin(gclMin(r, g), b);
        den = r + g + b;
        double ss = 1 - 3 * num / (den + EPSMIN);
        if (ss == 0)
        {
            hh = 0;
        }
        double ii = (r + g + b) / 3;
        Hvector[index] = hh;
        Svector[index] = ss;
        Ivector[index++] = ii;
    }
    return true;
}

bool HSI2RGB(double *Hvector, double *Svector, double *Ivector, uchar *imgOut, int width, int height)
{
    if (Hvector == NULL || Svector == NULL || Ivector == NULL)
    {
        return false;
    }
    int index = 0;
    for (int i = 0; i < height * width; i++)
    {
        double hh = Hvector[i] * 2 * PI;
        double ss = Svector[i];
        double ii = Ivector[i];
        double B = 0, G = 0, R = 0;
        if (hh >= 0 && hh < 2 * PI / 3)
        {
            B = ii * (1 - ss);
            R = ii * (1 + ss * cos(hh) / cos(PI / 3 - hh));
            G = 3 * ii - R - B;
        }
        else if (hh >= 2 * PI / 3 && hh < 4 * PI / 3)
        {
            R = ii * (1 - ss);
            G = ii * (1 + ss * cos(hh - 2 * PI / 3) / cos(PI - hh));
            B = 3 * ii - R - G;
        }
        else if (hh >= 4 * PI / 3 && hh <= 2 * PI)
        {
            G = ii * (1 - ss);
            B = ii * (1 + ss * cos(hh - 4 * PI / 3) / cos(5 * PI / 3 - hh));
            R = 3 * ii - G - B;
        }
        R = gclMax(0, gclMin(1, R));
        G = gclMax(0, gclMin(1, G));
        B = gclMax(0, gclMin(1, B));
        
        imgOut[index++] = uchar(B * 255);
        imgOut[index++] = uchar(G * 255);
        imgOut[index++] = uchar(R * 255);
    }
    return true;
}
/// 平滑图像

bool ImAverageFilter(uchar *imgIn, uchar *imgOut, int width, int height, int nchannels, int size)
{
    if (imgIn == NULL || size % 2 == 0 || size < 3 || !(nchannels == 1 || nchannels == 3))
    {
        return false;
    }
    int halfsize = size / 2;
    memcpy(imgOut, imgIn, sizeof(uchar) * width * height * nchannels);
    if (nchannels == 1)
    {
        for (int i = halfsize; i < height - halfsize; i++)
        {
            for (int j = halfsize; j < width - halfsize; j++)
            {
                double temp = 0;
                for (int m = -halfsize; m <= halfsize; m++)
                {
                    for (int n = -halfsize; n <= halfsize; n++)
                    {
                        temp += imgIn[(i + m) * width + j + n];
                    }
                }
                temp /= size * size;
                imgOut[i * width + j] = uchar(temp);
            }
        }
    }
    else if (nchannels == 3)
    {
        for (int i = halfsize; i < height - halfsize; i++)
        {
            for (int j = halfsize; j < width - halfsize; j++)
            {
                double tempb = 0, tempg = 0, tempr = 0;
                for (int m = -halfsize; m <= halfsize; m++)
                {
                    for (int n = -halfsize; n <= halfsize; n++)
                    {
                        tempb += imgIn[3 * ((i + m) * width + j + n) + 0];
                        tempg += imgIn[3 * ((i + m) * width + j + n) + 1];
                        tempr += imgIn[3 * ((i + m) * width + j + n) + 2];
                    }
                }
                tempb /= size * size;
                tempg /= size * size;
                tempr /= size * size;
                imgOut[3 * (i * width + j) + 0] = uchar(tempb);
                imgOut[3 * (i * width + j) + 1] = uchar(tempg);
                imgOut[3 * (i * width + j) + 2] = uchar(tempr);
            }
        }
    }
    return true;
}
/// 对比度增强
bool ImAdjust(uchar *imgIn, uchar *imgOut, int width, int height, int nchannels,
              
              double low_in, double high_in, double low_out, double high_out, double gamma)
{
    if (imgIn == NULL || !(nchannels == 1 || nchannels == 3) || low_in < 0 || low_in > 1 ||
        high_in < low_in || high_in > 1 || low_out < 0 || low_out > 1 || high_out < low_out || high_out > 1)
    {
        return false;
    }
    if ( nchannels == 1)
    {
        for (int i = 0; i < height * width; i++)
        {
            if (imgIn[i] < low_in * 255)
            {
                imgOut[i] = uchar(low_out * 255);
            }
            else if (imgIn[i] > high_in * 255)
            {
                imgOut[i] = uchar(high_out * 255);
            }
            else
            {
                double temp = pow((double(imgIn[i]) - low_in * 255) / (255 *(high_in - low_in)), gamma) *
                (high_out - low_out) * 255 + low_out * 255;
                imgOut[i] = uchar(temp);
            }
        }
    }
    else if (nchannels == 3)
    {
        double *Hvector = new double[width * height];
        double *Svector = new double[width * height];
        double *Ivector = new double[width * height];
        if (!RGB2HSI(imgIn, nchannels, Hvector, Svector, Ivector, width, height))
        {
            delete[] Hvector;
            delete[] Svector;
            delete[] Ivector;
            return false;
        }
        // enhancement the I vector
        for (int i = 0; i < width * height; i++)
        {
            if (Ivector[i] < low_in)
            {
                Ivector[i] = low_out;
            }
            else if (Ivector[i] > high_in)
            {
                Ivector[i] = high_out;
            }
            else
            {
                double temp = pow((Ivector[i] - low_in) / (high_in - low_in), gamma) *
                (high_out - low_out) + low_out;
                if (temp < low_out)
                {
                    temp = low_out;
                }
                if (temp > high_out)
                {
                    temp = high_out;
                }
                Ivector[i] = temp;
            }
        }
        HSI2RGB(Hvector, Svector, Ivector, imgOut, width, height);
        delete[] Hvector;
        delete[] Svector;
        delete[] Ivector;
    }
    return true;
}
/// 图像锐化滤波
bool ImSharpening(uchar *imgIn, uchar *imgOut, int widht, int height, int nchannels)
{
    if (imgIn == NULL || !(nchannels == 1 || nchannels == 3))
    {
        return false;
    }
    double kernel[] = {-1, -1, -1, -1, 8, -1, -1, -1, -1};
    double *temp = new double[widht * height];
    double *norm = new double[widht * height];
    if (nchannels == 1)
    {
        Conv2(imgIn, temp, widht, height, kernel, 9);
        Normalization(temp, norm, widht * height);
        for (int i = 0; i < widht * height; i++)
        {
            norm[i] += imgIn[i] / 255.0;
        }
        Normalization(norm, temp, widht * height);
        for (int i = 0; i < widht * height; i++)
        {
            imgOut[i] = uchar(temp[i] * 255);
        }
    }
    else if(nchannels == 3)
    {
        uchar *R = new uchar[widht * height];
        uchar *G = new uchar[widht * height];
        uchar *B = new uchar[widht * height];

    
        int index = 0;
        for (int i = 0; i < widht * height * nchannels; i = i + 3)
        {
            B[index] = imgIn[i + 0];
            G[index] = imgIn[i + 1];
            R[index ++] = imgIn[i + 2];
        }

        // B
        Conv2(B, temp, widht, height, kernel, 9);
        Normalization(temp, norm, widht * height);
        for (int i = 0; i < widht * height; i++)
        {
            norm[i] += B[i] / 255.0;
        }
        Normalization(norm, temp, widht * height);
        for (int i = 0; i < widht * height; i++)
        {
            B[i] = uchar(255 * temp[i]);
        }
        // G
        Conv2(G, temp, widht, height, kernel, 9);
        Normalization(temp, norm, widht * height);

        for (int i = 0; i < widht * height; i++)
        {
            norm[i] += G[i] / 255.0;
        }
        Normalization(norm, temp, widht * height);
        for (int i = 0; i < widht * height; i++)
        {
            G[i] = uchar(255 * temp[i]);
        }
        // R
        Conv2(R, temp, widht, height, kernel, 9);
        Normalization(temp, norm, widht * height);
        for (int i = 0; i < widht * height; i++)
        {
            norm[i] += R[i] / 255.0;
        }
        Normalization(norm, temp, widht * height);
        for (int i = 0; i < widht * height; i++)
        {
            R[i] = uchar(255 * temp[i]);
        }
        
        index = 0;

        for (int i = 0; i < widht * height; i++)
        {
            imgOut[index ++] = B[i];
            imgOut[index ++] = G[i];
            imgOut[index ++] = R[i];
        }
        delete[] R;
        delete[] G;
        delete[] B;
    }
    delete[] temp;
    delete[] norm;
    return true;
}
/// 伪彩色滤

bool ImColorMap(uchar *imgIn, uchar *imgOut, int widht, int height, int nchannels, double low, double high, int level)
{
    if (level < 8 || level % 8 != 0 || imgIn == NULL || nchannels != 1 || low < 0 || high > 1 || low >= high)
    {
        return false;
    }
    double *colorMap = new double[level * 3];

    GetColorMap(level, colorMap);
    
    uchar lowP = uchar(255 * low);
    uchar highP = uchar(255 * high);
    if (lowP >= highP)
    {
        return false;
    }
    
    int index = 0;
    for (int i = 0; i < widht * height; i++)
    {
        if (imgIn[i] < lowP)
        {
            imgOut[index ++] = uchar(255 * colorMap[0]);
            imgOut[index ++] = uchar(255 * colorMap[1]);
            imgOut[index ++] = uchar(255 * colorMap[2]);
        }
        else if (imgIn[i] > highP)
        {
            imgOut[index ++] = uchar(255 * colorMap[level * 3 - 3]);
            imgOut[index ++] = uchar(255 * colorMap[level * 3 - 2]);
            imgOut[index ++] = uchar(255 *colorMap[level * 3 - 1]);
        }
        else
        {
            int colorDex = int((imgIn[i] - lowP) / double(highP - lowP) * (level - 1));
            imgOut[index ++] = uchar( 255 *colorMap[colorDex * 3 + 0]);
            imgOut[index ++] = uchar(255 * colorMap[colorDex * 3 + 1]);
            imgOut[index ++] = uchar(255 * colorMap[colorDex * 3 + 2]);
        }
    }
    delete[] colorMap;
    return true;
}

/// 灰度图像直方图均衡化
void Histeq(uchar *imgIn, uchar *imgOut, int width, int height)
{
    const int level = 256;
    double *bin = new double[level];
    double *cunsumBin = new double[level];
    memset(bin, 0, sizeof(double) * level);
    memset(cunsumBin, 0, sizeof(double) * level);
    for (int i = 0; i < width * height; i++)
    {
        bin[imgIn[i]]++;
    }
    for (int i = 0; i < level; i++)
    {
        bin[i] /= width * height;
    }
    for (int i = 0; i < level; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            cunsumBin[i] += bin[j];
        }
    }
    for (int i = 0; i < width * height; i++)
    {
        imgOut[i] = uchar(255 * cunsumBin[imgIn[i]] + 0.5);
    }
    
    delete[] cunsumBin;
    delete[] bin;
    
}

/// 直方图均衡化，解决过曝及欠曝
bool ImHisteq(uchar *imgIn, uchar *imgOut, int width, int height, int nchannels)
{
    if (imgIn == NULL || !(nchannels == 1 || nchannels == 3))
    {
        return false;
    }
    if (nchannels == 1)
    {
        Histeq(imgIn, imgOut, width, height);
    }
    else if (nchannels == 3)
    {
        uchar *R = new uchar[width * height];
        uchar *G = new uchar[width * height];
        uchar *B = new uchar[width * height];
        uchar *r = new uchar[width * height];
        uchar *g = new uchar[width * height];
        uchar *b = new uchar[width * height];
        int index = 0;
        for (int i = 0; i < width * height * nchannels; i = i + 3)
        {
            B[index] = imgIn[i + 0];
            G[index] = imgIn[i + 1];
            R[index++] = imgIn[i + 2];
        }
        Histeq(R, r, width, height);
        Histeq(G, g, width, height);
        Histeq(B, b, width, height);
        
        index = 0;
        for (int i = 0; i < width * height; i++)
        {
            imgOut[index ++] = b[i];
            imgOut[index ++] = g[i];
            imgOut[index ++] = r[i];
        }
        delete[] R;
        delete[] G;
        delete[] B;
        delete[] r;
        delete[] g;
        delete[] b;
    }
    return true;
}
