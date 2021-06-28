#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"

#define M_PI 3.14159265358979323846

// HW1 #2.1
// Image& im: image to L1-normalize
void l1_normalize(Image& im)
  {
  
  // TODO: Normalize each channel
	for (int c = 0; c < im.c; c++) {
		float sum = 0;
		for (int i = 0; i < im.w; i++) {
			for (int j = 0; j < im.h; j++) {
				sum += im(i, j, c);
			}
		}
		for (int i = 0; i < im.w; i++) {
			for (int j = 0; j < im.h; j++) {
				im(i, j, c) /= sum;
			}
		}
	}
  
  
  }

// HW1 #2.1
// int w: size of filter
// returns the filter Image of size WxW
Image make_box_filter(int w)
  {
  assert(w%2); // w needs to be odd
  
  // TODO: Implement the filter
  Image im = Image(w, w, 1);
  for (int i = 0; i < w; i++) {
	  for (int j = 0; j < w; j++) {
		  im(i, j) = 1.0 / (w*w);
	  }
  }
  
  return im;
  }

// HW1 #2.2
// const Image&im: input image
// const Image& filter: filter to convolve with
// bool preserve: whether to preserve number of channels
// returns the convolved image
Image convolve_image(const Image& im, const Image& filter, bool preserve)
  {
  assert(filter.c==1);
  Image ret;
  // This is the case when we need to use the function clamped_pixel(x,y,c).
  // Otherwise you'll have to manually check whether the filter goes out of bounds
  int x = filter.w / 2;
  int y = filter.h / 2;
  if (preserve) {
	  ret = Image(im.w, im.h, im.c);
	  for (int c = 0; c < im.c; c++) {
		  for (int i = 0; i < im.w; i++) {
			  for (int j = 0; j < im.h; j++) {
				  float q = 0;
				  for (int m = 0; m < filter.w; m++) {
					  for (int n = 0; n < filter.h; n++) {
						  q += filter(m, n)*im.clamped_pixel(i+m-x, j+n-y, c);
					  }
				  }
				  ret(i, j, c) = q;
			  }
		  }
	  }
  }
  else {
	  ret = Image(im.w, im.h, 1);
	  for (int i = 0; i < im.w; i++) {
		  for (int j = 0; j < im.h; j++) {
			  float q = 0;
			  for (int c = 0; c < im.c; c++) {
				  for (int m = 0; m < filter.w; m++) {
					  for (int n = 0; n < filter.h; n++) {
						  q += filter(m, n)*im.clamped_pixel(i+m-x, j+n-y, c);
					  }
				  }
			  }
			  ret(i, j) = q;
		  }
	  }
  }
  // TODO: Make sure you set the sizes of ret properly. Use ret=Image(w,h,c) to reset ret
  // TODO: Do the convolution operator

  
  // Make sure to return ret and not im. This is just a placeholder
  return ret;
  }

// HW1 #2.3
// returns basic 3x3 high-pass filter
Image make_highpass_filter()
  {
  // TODO: Implement the filter
	Image filter = make_box_filter(3);
	filter(0, 0) = filter(2, 0) = filter(0, 2) = filter(2, 2) = 0;
	filter(1, 0) = filter(0, 1) = filter(2, 1) = filter(1, 2) = -1;
	filter(1, 1) = 4;
  
  return filter;
  
  }

// HW1 #2.3
// returns basic 3x3 sharpen filter
Image make_sharpen_filter()
  {
  // TODO: Implement the filter
	Image filter = make_box_filter(3);
	filter(0, 0) = filter(2, 0) = filter(0, 2) = filter(2, 2) = 0;
	filter(1, 0) = filter(0, 1) = filter(2, 1) = filter(1, 2) = -1;
	filter(1, 1) = 5;

	return filter;
  
  }

// HW1 #2.3
// returns basic 3x3 emboss filter
Image make_emboss_filter()
  {
  // TODO: Implement the filter
	Image filter = make_box_filter(3);
	filter(0, 0) = -2;
	filter(2, 0) = filter(0, 2) = 0;
	filter(1, 0) = filter(0, 1) = -1;
	filter(1, 1) = filter(2, 1) = filter(1, 2) = 1;
	filter(2, 2) = 2;

	return filter;
  
  }

float compute_gaussian(float sigma, int x, int y) {
	return expf(-(pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2))) / (2 * M_PI * pow(sigma, 2));
}

float compute_1d_gaussian(float sigma, float dist) {
	return expf(-(powf(dist, 2)) / (2 * pow(sigma, 2))) / (2 * M_PI * pow(sigma, 2));
}

// HW1 #2.4
// float sigma: sigma for the gaussian filter
// returns basic gaussian filter
Image make_gaussian_filter(float sigma)
  {
  // TODO: Implement the filter
	int w = ceil(6 * sigma);
	w = w % 2 == 0 ? w + 1 : w;
	Image filter = make_box_filter(w);
	int offset = w / 2;
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < w; j++) {
			filter(i, j) = compute_gaussian(sigma, (i-offset), (j-offset));
		}
	}
	l1_normalize(filter);
  
  return filter;
  
  }


// HW1 #3
// const Image& a: input image
// const Image& b: input image
// returns their sum
Image add_image(const Image& a, const Image& b)
  {
  assert(a.w==b.w && a.h==b.h && a.c==b.c); // assure images are the same size
  int w = a.w;
  int h = a.h;
  int c = a.c;
  Image im = Image(w, h, c);
  // TODO: Implement addition
  for (int k = 0; k < c; k++) {
	  for (int i = 0; i < w; i++) {
		  for (int j = 0; j < h; j++) {
			  im(i, j, k) = a(i, j, k) + b(i, j, k);
		  }
	  }
  }
  
  return im;
  
  }

// HW1 #3
// const Image& a: input image
// const Image& b: input image
// returns their difference res=a-b
Image sub_image(const Image& a, const Image& b)
  {
  assert(a.w==b.w && a.h==b.h && a.c==b.c); // assure images are the same size
  int w = a.w;
  int h = a.h;
  int c = a.c;
  Image im = Image(w, h, c);
  // TODO: Implement subtraction
  for (int k = 0; k < c; k++) {
	  for (int i = 0; i < w; i++) {
		  for (int j = 0; j < h; j++) {
			  im(i, j, k) = a(i, j, k) - b(i, j, k);
		  }
	  }
  }
  
  return im;
  
  }

// HW1 #4.1
// returns basic GX filter
Image make_gx_filter()
  {
  // TODO: Implement the filter
	Image f = make_box_filter(3);
	f(0, 0) = f(0, 2) = -1;
	f(1, 0) = f(1, 1) = f(1, 2) = 0;
	f(2, 0) = f(2, 2) = 1;
	f(0, 1) = -2;
	f(2, 1) = 2;
  
  return f;
  }

// HW1 #4.1
// returns basic GY filter
Image make_gy_filter()
  {
  // TODO: Implement the filter
	Image f = make_box_filter(3);
	f(0, 0) = f(2, 0) = -1;
	f(0, 1) = f(1, 1) = f(2, 1) = 0;
	f(0, 2) = f(2, 2) = 1;
	f(1, 0) = -2;
	f(1, 2) = 2;

	return f;
  }

// HW1 #4.2
// Image& im: input image
void feature_normalize(Image& im)
  {
  assert(im.w*im.h); // assure we have non-empty image
  
  // TODO: Normalize the features for each channel
  float l = 1000000000000000;
  float r = -1000000000000000;
  for (int c = 0; c < im.c; c++) {
	  for (int i = 0; i < im.w; i++) {
		  for (int j = 0; j < im.h; j++) {
			  l = min(l, im(i, j, c));
			  r = max(r, im(i, j, c));
		  }
	  }
  }
  float range = r - l;
  for (int c = 0; c < im.c; c++) {
	  for (int i = 0; i < im.w; i++) {
		  for (int j = 0; j < im.h; j++) {
			  im(i, j, c) = range == 0 ? 0 : (im(i, j, c) - l) / range;
		  }
	  }
  }
  
  }


// Normalizes features across all channels
void feature_normalize_total(Image& im)
  {
  assert(im.w*im.h*im.c); // assure we have non-empty image
  
  int nc=im.c;
  im.c=1;im.w*=nc;
  
  feature_normalize(im);
  
  im.w/=nc;im.c=nc;
  
  }


// HW1 #4.3
// Image& im: input image
// return a pair of images of the same size
pair<Image,Image> sobel_image(const Image& im)
  {
  // TODO: Your code here
	Image m_x = convolve_image(im, make_gx_filter(), false);
	Image m_y = convolve_image(im, make_gy_filter(), false);
	Image mag = Image(im.w, im.h);
	Image grad = Image(im.w, im.h);
	for (int i = 0; i < im.w; i++) {
		for (int j = 0; j < im.h; j++) {
			mag(i, j) = sqrtf(powf(m_x(i, j), 2) + powf(m_y(i, j), 2));
			grad(i, j) = atan2(m_y(i, j), m_x(i, j));
		}
	}
  return {mag,grad};
  }


// HW1 #4.4
// const Image& im: input image
// returns the colorized Sobel image of the same size
Image colorize_sobel(const Image& im)
  {
  
  // TODO: Your code here
	pair<Image, Image> sobel_res = sobel_image(im);
	Image mag = sobel_res.first;
	Image grad = sobel_res.second;
	Image hsv_img = Image(mag.w, mag.h, 3);

	//Nomalize the magnitude:
	feature_normalize(mag);
	//Nomalize the angle and get the hsv image
	for (int i = 0; i < grad.w; i++) {
		for (int j = 0; j < grad.h; j++) {
			grad(i, j) = grad(i, j) / (2 * M_PI) + 0.5;
			hsv_img(i, j, 0) = grad(i, j);
			hsv_img(i, j, 1) = hsv_img(i, j, 2) = mag(i, j);
		}
	}

	//Gaussian smooth the image
	hsv_to_rgb(hsv_img);
	Image filter = make_gaussian_filter(4);
  return convolve_image(hsv_img, filter, true);
  }


// HW1 #4.5
// const Image& im: input image
// float sigma1,sigma2: the two sigmas for bilateral filter
// returns the result of applying bilateral filtering to im
Image bilateral_filter(const Image& im, float sigma1, float sigma2)
  {
	int w = ceil(6 * sigma1);
	w = w % 2 == 0 ? w + 1 : w;
	int k = w / 2;
	Image bf = make_box_filter(w);
	Image res = Image(im.w, im.h, im.c);
	//
	Image spatial = make_gaussian_filter(sigma1);

	for (int c = 0; c < im.c; c++) {
		for (int x = 0; x < im.w; x++) {
			for (int y = 0; y < im.h; y++) {
				//Compute the weights
				float N = 0;
				for (int i = 0; i < w; i++) {
					for (int j = 0; j < w; j++) {
						bf(i, j) = spatial(i, j) * compute_1d_gaussian(sigma2, (im.clamped_pixel(x, y, c) - im.clamped_pixel(x + i - k, y + j - k, c)));
						N += bf(i, j);
					}
				}
				l1_normalize(bf);
				//Convolution
				float q = 0;
				for (int m = 0; m < bf.w; m++) {
					for (int n = 0; n < bf.w; n++) {
						q += bf(m, n)*im.clamped_pixel(x + m - k, y + n - k, c);
					}
				}
				res(x, y, c) = q;
			 }
		}
	}
  
  // TODO: Your bilateral code
  
  
  return res;
  }



// HELPER MEMBER FXNS

void Image::feature_normalize(void) { ::feature_normalize(*this); }
void Image::feature_normalize_total(void) { ::feature_normalize_total(*this); }
void Image::l1_normalize(void) { ::l1_normalize(*this); }

Image operator-(const Image& a, const Image& b) { return sub_image(a,b); }
Image operator+(const Image& a, const Image& b) { return add_image(a,b); }
