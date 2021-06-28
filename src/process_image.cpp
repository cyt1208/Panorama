#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include "image.h"

using namespace std;

// HW0 #3
// const Image& im: input image
// return the corresponding grayscale image
Image rgb_to_grayscale(const Image& im)
  {
  assert(im.c == 3); // only accept RGB images
  Image gray(im.w,im.h,1); // create a new grayscale image (note: 1 channel)
  
  for (int i = 0; i < im.w; i++) {
	  for (int j = 0; j < im.h; j++) {
		  gray(i, j) = 0.299*im(i, j, 0) + 0.587*im(i, j, 1)+0.114*im(i, j, 2);
	  }
  }
  
  return gray;
  }



// Example function that changes the color of a grayscale image
Image grayscale_to_rgb(const Image& im, float r, float g, float b)
  {
  assert(im.c == 1);
  Image rgb(im.w,im.h,3);
  
  for(int q2=0;q2<im.h;q2++)for(int q1=0;q1<im.w;q1++)
    {
    rgb(q1,q2,0)=r*im(q1,q2);
    rgb(q1,q2,1)=g*im(q1,q2);
    rgb(q1,q2,2)=b*im(q1,q2);
    }
  
  return rgb;
  }




// HW0 #4
// Image& im: input image to be modified in-place
// int c: which channel to shift
// float v: how much to shift
void shift_image(Image& im, int c, float v)
  {
  assert(c>=0 && c<im.c); // needs to be a valid channel
  
  for (int i = 0; i < im.w; i++) {
	  for (int j = 0; j < im.h; j++) {
		  im(i, j, c) += v;
	  }
  }
  
  }

// HW0 #8
// Image& im: input image to be modified in-place
// int c: which channel to scale
// float v: how much to scale
void scale_image(Image& im, int c, float v)
  {
  assert(c>=0 && c<im.c); // needs to be a valid channel
  
  for (int i = 0; i < im.w; i++) {
	  for (int j = 0; j < im.h; j++) {
		  im(i, j, c) *= v;
	  }
  }
  
  }


// HW0 #5
// Image& im: input image to be modified in-place
void clamp_image(Image& im)
  {
	for (int i = 0; i < im.w; i++) {
		for (int j = 0; j < im.h; j++) {
			for (int k = 0; k < im.c; k++) {
				im(i, j, k) = im(i, j, k) < 0 ? 0 : im(i, j, k);
				im(i, j, k) = im(i, j, k) > 1 ? 1 : im(i, j, k);
			}
		}
  }
  
  }

// These might be handy
float max(float a, float b, float c)
  {
  return max({a,b,c});
  }

float min(float a, float b, float c)
  {
  return min({a,b,c});
  }

float computeH(float R, float G, float B, float V, float C) {
	if (C == 0)
		return 0;
	float h = 0;
	if (V == R)
		h = (G - B) / C;
	if (V == G)
		h = (B - R) / C + 2;
	if (V == B)
		h = (R - G) / C + 4;
	return h < 0 ? h / 6 + 1 : h / 6;
}

// HW0 #6
// Image& im: input image to be modified in-place
void rgb_to_hsv(Image& im)
  {
  assert(im.c==3 && "only works for 3-channels images");
  
  for (int i = 0; i < im.w; i++) {
	  for (int j = 0; j < im.h; j++) {
		  float R = im(i, j, 0);
		  float G = im(i, j, 1);
		  float B = im(i, j, 2);
		  float V = max(R, G, B);
		  float C = V - min(R, G, B);
		  float S = V == 0 ? 0 : C / V;
		  float H = computeH(R, G, B, V, C);
		  im(i, j, 0) = H;
		  im(i, j, 1) = S;
		  im(i, j, 2) = V;
	  }
  }
  
  }

// HW0 #7
// Image& im: input image to be modified in-place
void hsv_to_rgb(Image& im)
  {
  assert(im.c==3 && "only works for 3-channels images");
  
  for (int i = 0; i < im.w; i++) {
	  for (int j = 0; j < im.h; j++) {
		  float H = im(i, j, 0);
		  float S = im(i, j, 1);
		  float V = im(i, j, 2);
		  float C = V * S;
		  float X = C * (1 - abs(fmod(6 * H, 2.0) - 1));
		  float m = V - C;
		  float temp = 1.0 / 6.0;
		  if (H >= 0 && H < temp) {
			  im(i, j, 0) = C + m;
			  im(i, j, 1) = X + m;
			  im(i, j, 2) = m;
		  }
		  if (H >= temp && H < 2 * temp) {
			  im(i, j, 0) = X + m;
			  im(i, j, 1) = C + m;
			  im(i, j, 2) = m;
		  }
		  if (H >= 2 * temp && H < 3 * temp) {
			  im(i, j, 0) = m;
			  im(i, j, 1) = C + m;
			  im(i, j, 2) = X + m;
		  }
		  if (H >= 3 * temp && H < 4 * temp) {
			  im(i, j, 0) = m;
			  im(i, j, 1) = X + m;
			  im(i, j, 2) = C + m;
		  }
		  if (H >= 4 * temp && H < 5 * temp) {
			  im(i, j, 0) = X + m;
			  im(i, j, 1) = m;
			  im(i, j, 2) = C + m;
		  }
		  if (H >= 5 * temp && H < 1) {
			  im(i, j, 0) = C + m;
			  im(i, j, 1) = m;
			  im(i, j, 2) = X + m;
		  }
	  }
  }
  
  }

// HW0 #9
// Image& im: input image to be modified in-place
void rgb_to_lch(Image& im)
{
	assert(im.c == 3 && "only works for 3-channels images");

	for (int i = 0; i < im.w; i++) {
		for (int j = 0; j < im.h; j++) {
			//gamma decompression while gamma = 2.2.
			float r = pow(im(i, j, 0), 2.2);
			float g = pow(im(i, j, 1), 2.2);
			float b = pow(im(i, j, 2), 2.2);

			//I use the RGB Working Space Matrix 
			//where the RGB working space is Adobe RGB(1998) 
			//and the reference white is D65. The reference link is:
			//http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
			float X = 0.5767309*r + 0.1855540*g + 0.1881852*b;
			float Y = 0.2973769*r + 0.6273491*g + 0.0752741*b;
			float Z = 0.0270343*r + 0.0706872*g + 0.9911085*b;

			//For reference white, I use D65 with x, y, z value: 0.9504, 1.0000, 1.0888
			//Link: https://www.mathworks.com/help/images/ref/whitepoint.html
			//I set k = 903.3, e = 0.008856.
			float x = 0.9504;
			float y = 1.0000;
			float z = 1.0888;
			float yr = Y / y;
			float uu = 4 * X / (X + 15 * Y + 3 * Z);
			float vv = 9 * Y / (X + 15 * Y + 3 * Z);
			float ur = 4 * x / (x + 15 * y + 3 * z);
			float vr = 9 * y / (x + 15 * y + 3 * z);
			float L = yr > 0.008856 ? 116 * pow(yr, 1.0 / 3.0) - 16 : 903.3*yr;
			float u = 13 * L*(uu - ur);
			float v = 13 * L*(vv - vr);

			float C = sqrt(pow(u, 2) + pow(v, 2));
			float H = atan2(v, u);
			im(i, j, 0) = L;
			im(i, j, 1) = C;
			im(i, j, 2) = H;
		}
	}
}


// HW0 #9
// Image& im: input image to be modified in-place
void lch_to_rgb(Image& im)
  {
  assert(im.c==3 && "only works for 3-channels images");
  
  for (int i = 0; i < im.w; i++) {
	  for (int j = 0; j < im.h; j++) {
		  float L = im(i, j, 0);
		  float u = im(i, j, 1) * cos(im(i, j, 2));
		  float v = im(i, j, 1) * sin(im(i, j, 2));

		  //all the parameters are the same as them in the rgb_to_lch function.
		  float x = 0.9504;
		  float y = 1.0000;
		  float z = 1.0888;
		  float Y = L > 0.008856*903.3 ? pow((L + 16) / 116, 3) : L / 903.3;
		  float u0 = 4 * x / (x + 15 * y + 3 * z);
		  float v0 = 9 * y / (x + 15 * y + 3 * z);
		  float a = (52 * L / (u + 13 * L * u0) - 1) / 3;
		  float b = -5 * Y;
		  float c = -1.0 / 3.0;
		  float d = Y * (39 * L / (v + 13 * L*v0) - 5);
		  float X = (d - b) / (a - c);
		  float Z = X * a + b;

		  //I also use the inverse matrix of the matrix I used above.
		  float R = 2.0413690*X - 0.5649464*Y - 0.3446944*Z;
		  float G = -0.9692660*X + 1.8760108*Y + 0.0415560*Z;
		  float B = 0.0134474*X - 0.1183897*Y + 1.0154096*Z;
		  im(i, j, 0) = pow(R, 1 / 2.2);
		  im(i, j, 1) = pow(G, 1 / 2.2);
		  im(i, j, 2) = pow(B, 1 / 2.2);
		  im(i, j, 0) = im(i, j, 0) != im(i, j, 0) ? 0 : im(i, j, 0);
		  im(i, j, 1) = im(i, j, 1) != im(i, j, 1) ? 0 : im(i, j, 1);
		  im(i, j, 2) = im(i, j, 2) != im(i, j, 2) ? 0 : im(i, j, 2);
	  }
  }
  
  }



// Implementation of member functions
void Image::clamp(void) { clamp_image(*this); }
void Image::shift(int c, float v) { shift_image(*this,c,v); }
void Image::scale(int c, float v) { scale_image(*this,c,v); }

void Image::HSVtoRGB(void) { hsv_to_rgb(*this); }
void Image::RGBtoHSV(void) { rgb_to_hsv(*this); }
void Image::LCHtoRGB(void) { lch_to_rgb(*this); }
void Image::RGBtoLCH(void) { rgb_to_lch(*this); }
