#include <cmath>
#include "image.h"

using namespace std;

// HW1 #1
// float x,y: inexact coordinates
// int c: channel
// returns the nearest neibor to pixel (x,y,c)
float Image::pixel_nearest(float x, float y, int c) const
  {
  // Since you are inside class Image you can
  // use the member function pixel(a,b,c)
  
  // TODO: Your code here
  
	int near_x = roundf(x);
	int near_y = roundf(y);
  
  
	return clamped_pixel(near_x, near_y, c);
  }

// HW1 #1
// float x,y: inexact coordinates
// int c: channel
// returns the bilinearly interpolated pixel (x,y,c)
float Image::pixel_bilinear(float x, float y, int c) const
  {
  // Since you are inside class Image you can
  // use the member function pixel(a,b,c)
  
  
  // TODO: Your code here
  
	int x1 = floor(x);
	int x2 = floor(x)+1;
	int y1 = floor(y);
	int y2 = floor(y)+1;
	float v1 = clamped_pixel(x1, y1, c);
	float v2 = clamped_pixel(x2, y1, c);
	float v3 = clamped_pixel(x1, y2, c);
	float v4 = clamped_pixel(x2, y2, c);
  
  
  return v1*(x2-x)*(y2-y)+v2*(x-x1)*(y2-y)+v3*(x2-x)*(y-y1)+v4*(x-x1)*(y-y1);
  }

// HW1 #1
// int w,h: size of new image
// const Image& im: input image
// return new Image of size (w,h,im.c)
Image nearest_resize(const Image& im, int w, int h)
  {
  Image ret(w,h,im.c);
  
  // TODO: Your code here
  float scale_w = (im.w * 1.0) / w;
  float scale_h = (im.h * 1.0) / h;
  for (int c = 0; c < im.c; c++) {
	  for (int i = 0; i < w; i++) {
		  for (int j = 0; j < h; j++) {
			  ret(i, j, c) = im.pixel_nearest(scale_w*(i+0.5)-0.5, scale_h*(j+0.5)-0.5, c);
		  }
	  }
  }
  
  
  
  return ret;
  }


// HW1 #1
// int w,h: size of new image
// const Image& im: input image
// return new Image of size (w,h,im.c)
Image bilinear_resize(const Image& im, int w, int h)
  {
	Image ret(w, h, im.c);
  // TODO: Your code here
  
	float scale_w = (im.w * 1.0) / w;
	float scale_h = (im.h * 1.0) / h;
	for (int c = 0; c < im.c; c++) {
		for (int i = 0; i < w; i++) {
			for (int j = 0; j < h; j++) {
				ret(i, j, c) = im.pixel_bilinear(scale_w*(i + 0.5) - 0.5, scale_h*(j + 0.5) - 0.5, c);
			}
		}
	}
  
  
  return ret;
  }


