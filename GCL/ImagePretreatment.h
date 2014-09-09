//
//  ImagePretreatment.h
//  GCL
//
//  Created by Jia Yonglei on 14-9-9.
//  Copyright (c) 2014年 Jia Yonglei. All rights reserved.
//

#ifndef __GCL__ImagePretreatment__
#define __GCL__ImagePretreatment__

#include <stdio.h>
#include "DateStruct.h"

typedef unsigned char uchar;

bool ImAverageFilter(uchar *imgIn, uchar *imgOut, int width, int height, int nchannels, int size);
bool ImAdjust(uchar *imgIn, uchar *imgOut, int width, int height, int nchannels,
              double low_in, double high_in, double low_out, double high_out, double gamma);
bool ImSharpening(uchar *imgIn, uchar *imgOut, int widht, int height, int nchannels);
bool ImColorMap(uchar *imgIn, uchar *imgOut, int widht, int height, int nchannels, double low, double high, int level);

#endif /* defined(__GCL__ImagePretreatment__) */
