// Source file for image class


#define _USE_MATH_DEFINES

// Include files 

#include <iostream>
#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include <math.h>
#include <vector>
#include <array>
#include <algorithm>
#include <valarray>
#include <iterator>

using namespace std;
using std::vector;
using std::sort;
using std::slice;


////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////

struct Feature
{
	int centerX;
	int centerY;
	R2Pixel HarrisValue;

	Feature()
	{
		centerX = -1;
		centerY = -1;
	}

	Feature(int x, int y, R2Pixel val)
	{
		centerX = x;
		centerY = y;
		HarrisValue = val;
	}

	bool operator < (const Feature& feature) const {
		double valueIntensity = HarrisValue[0] + HarrisValue[1] + HarrisValue[2];
		double featureIntensity = feature.HarrisValue[0] + feature.HarrisValue[1] + feature.HarrisValue[2];

		return valueIntensity < featureIntensity;
	}

	bool operator == (const Feature& feature) const {
		double valueIntensity = HarrisValue[0] + HarrisValue[1] + HarrisValue[2];
		double featureIntensity = feature.HarrisValue[0] + feature.HarrisValue[1] + feature.HarrisValue[2];

		return valueIntensity == featureIntensity;
	}

};

R2Image::
R2Image(void)
	: pixels(NULL),
	npixels(0),
	width(0),
	height(0)
{
}



R2Image::
R2Image(const char *filename)
	: pixels(NULL),
	npixels(0),
	width(0),
	height(0)
{
	// Read image
	Read(filename);
}



R2Image::
R2Image(int width, int height)
	: pixels(NULL),
	npixels(width * height),
	width(width),
	height(height)
{
	// Allocate pixels
	pixels = new R2Pixel[npixels];
	assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
	: pixels(NULL),
	npixels(width * height),
	width(width),
	height(height)
{
	// Allocate pixels
	pixels = new R2Pixel[npixels];
	assert(pixels);

	// Copy pixels 
	for (int i = 0; i < npixels; i++)
		pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
	: pixels(NULL),
	npixels(image.npixels),
	width(image.width),
	height(image.height)

{
	// Allocate pixels
	pixels = new R2Pixel[npixels];
	assert(pixels);

	// Copy pixels 
	for (int i = 0; i < npixels; i++)
		pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
	// Free image pixels
	if (pixels) delete[] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
	// Delete previous pixels
	if (pixels) { delete[] pixels; pixels = NULL; }

	// Reset width and height
	npixels = image.npixels;
	width = image.width;
	height = image.height;

	// Allocate new pixels
	pixels = new R2Pixel[npixels];
	assert(pixels);

	// Copy pixels 
	for (int i = 0; i < npixels; i++)
		pixels[i] = image.pixels[i];

	// Return image
	return *this;
}

double** R2Image::
svdTest(double* xcoords, double *ycoords, int size)
{
	double **matrixH = dmatrix(1, 3, 1, 3);

	int n = size;

	 /*int *xcoords = (int*)malloc(8 * sizeof(int));
	 int *ycoords = (int*)malloc(8 * sizeof(int));*/

	// Matrix A is a 2n x 9 matrix 
	double** matrixA = dmatrix(1, 2 * n, 1, 9);
	int j = 1;
	for (int i = 0; i < 2 * n; i+=2) 
	{

		matrixA[j][1] = 0;
		matrixA[j][2] = 0;
		matrixA[j][3] = 0;
		matrixA[j][4] = -xcoords[i];
		matrixA[j][5] = -ycoords[i];
		matrixA[j][6] = -1;
		matrixA[j][7] = ycoords[i + 1] * xcoords[i];
		matrixA[j][8] = ycoords[i + 1] * ycoords[i];
		matrixA[j][9] = ycoords[i + 1];

		matrixA[j + 1][1] = xcoords[i];
		matrixA[j + 1][2] = ycoords[i];
		matrixA[j + 1][3] = 1;
		matrixA[j + 1][4] = 0;
		matrixA[j + 1][5] = 0;
		matrixA[j + 1][6] = 0;
		matrixA[j + 1][7] = -xcoords[i + 1] * xcoords[i];
		matrixA[j + 1][8] = -xcoords[i + 1] * ycoords[i];
		matrixA[j + 1][9] = -xcoords[i + 1];

		j += 2;
	}

	// compute the SVD
	double singularValues[10]; // 1..9
	double** nullspaceMatrix = dmatrix(1, 9, 1, 9);
	svdcmp(matrixA, 8, 9, singularValues, nullspaceMatrix);

	// get the result
	// printf("\n Singular values: %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", singularValues[1], singularValues[2], singularValues[3], singularValues[4], singularValues[5], singularValues[6], singularValues[7], singularValues[8], singularValues[9], singularValues[10]);

	// find the smallest singular value:
	int smallestIndex = 1;
	for (int i = 2; i < 10; i++) if (singularValues[i] < singularValues[smallestIndex]) smallestIndex = i;

	// solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)

	/*printf("Matrix H: \n[%f, %f, %f, \n", nullspaceMatrix[1][smallestIndex], nullspaceMatrix[2][smallestIndex], nullspaceMatrix[3][smallestIndex]);
	printf("%f, %f, %f, \n", nullspaceMatrix[4][smallestIndex], nullspaceMatrix[5][smallestIndex], nullspaceMatrix[6][smallestIndex]);
	printf("%f, %f, %f] \n", nullspaceMatrix[7][smallestIndex], nullspaceMatrix[8][smallestIndex], nullspaceMatrix[9][smallestIndex]);*/

	// Save the new matrix H

	matrixH[1][1] = nullspaceMatrix[1][smallestIndex];
	matrixH[1][2] = nullspaceMatrix[2][smallestIndex];
	matrixH[1][3] = nullspaceMatrix[3][smallestIndex];

	matrixH[2][1] = nullspaceMatrix[4][smallestIndex];
	matrixH[2][2] = nullspaceMatrix[5][smallestIndex];
	matrixH[2][3] = nullspaceMatrix[6][smallestIndex];

	matrixH[3][1] = nullspaceMatrix[7][smallestIndex];
	matrixH[3][2] = nullspaceMatrix[8][smallestIndex];
	matrixH[3][3] = nullspaceMatrix[9][smallestIndex];

	free(xcoords);
	free(ycoords);
	return matrixH;
}


////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
	// Brighten the image by multiplying each pixel component by the factor.
	// This is implemented for you as an example of how to access and set pixels
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			Pixel(i, j) *= factor;
			Pixel(i, j).Clamp();
		}
	}
}

void R2Image::
SobelX(void)
{
	// Convert RGB to Grayscale 
	
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			const R2Pixel& pixel = Pixel(i, j);

			double r = (double)(pixel.Red());
			double g = (double)(pixel.Green());
			double b = (double)(pixel.Blue());

			double gray = (r + g + b) / 3;

			R2Pixel newPixel(gray, gray, gray, 1.0);
			// newPixel.Clamp();
			(*this).SetPixel(i, j, newPixel);
		}
	}
	
	// Apply SobelX 
	R2Image temp(*this);
	int sobel[] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };

	for (int y = 1; y < height - 1; y++) {
		for (int x = 1; x < width - 1; x++) {
			R2Pixel val(0.0, 0.0, 0.0, 1.0);
			int count = 0;
			for (int ly = -1; ly <= 1; ly++) {
				val += Pixel(x, y + ly) * sobel[count];
				count++;
			}
			// val.Clamp();
			temp.SetPixel(x, y, val);
		}
	}
	(*this) = temp;
}

void R2Image::
SobelY(void)
{
	// Apply the Sobel oprator to the image in Y direction

	// Convert RGB to Grayscale
	
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			const R2Pixel& pixel = Pixel(i, j);

			double r = (double)(pixel.Red());
			double g = (double)(pixel.Green());
			double b = (double)(pixel.Blue());

			double gray = (r + g + b) / 3;

			R2Pixel newPixel(gray, gray, gray, 1.0);
			// newPixel.Clamp();
			(*this).SetPixel(i, j, newPixel);
		}
	}
	
	// Apply SobelY
	R2Image temp(*this);
	int sobel[] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };

	for (int y = 1; y < height - 1; y++) {
		for (int x = 1; x < width - 1; x++) {
			R2Pixel val(0.0, 0.0, 0.0, 1.0);
			int count = 0;
			for (int lx = -1; lx <= 1; lx++) {
				val += Pixel(x + lx, y) * sobel[count];
				count++;
			}
			// val.Clamp();
			temp.SetPixel(x, y, val);
		}
	}
	(*this) = temp;
}

void R2Image::
LoG(void)
{
	// Apply the LoG oprator to the image

	// FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
	fprintf(stderr, "LoG() not implemented\n");
}


void R2Image::
ChangeSaturation(double factor)
{
	// Changes the saturation of an image

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			const R2Pixel& pixel = Pixel(i, j);

			double r = (double)(pixel.Red());
			double g = (double)(pixel.Green());
			double b = (double)(pixel.Blue());

			double p = sqrt(r * r * .299 + g * g * .587 + b * b * .114);

			r = p + (r - p) * factor;
			g = p + (g - p) * factor;
			b = p + (b - p) * factor;

			R2Pixel newPixel(r, g, b, 1.0);
			newPixel.Clamp();
			(*this).SetPixel(i, j, newPixel);

		}
	}
}

// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
	// Gaussian blur of the image. 

	// Gaussian kernel
	int size = 6 * sigma + 1;
	double *kernel_h = (double*)malloc(size * sizeof(double));

	double sum = 0.0;
	int w = 3 * sigma;

	for (int i = 0; i < size; i++) {
		kernel_h[i] = exp(-0.5*(pow((i - size / 2) / sigma, 2.0))) / sqrt(2.0*M_PI*pow(sigma, 2.0));
		sum += kernel_h[i];
	}

	// Normalize the kernel so that all the weights add up to one
	for (int i = 0; i < size; i++) {
		kernel_h[i] /= sum;
	}

	R2Image temp(*this);

	// Apply the Gaussian kernel in the X direction
	for (int y = w; y < height - w; y++) {
		for (int x = w; x < width - w; x++) {
			R2Pixel val(0.0, 0.0, 0.0, 1.0);
			int count = 0;
			for (int i = -w; i <= w; i++) {
				val += kernel_h[count] * Pixel(x + i, y);
				count++;
			}
			// val.Clamp();
			temp.SetPixel(x, y, val);
		}
	}
	// Apply the Gaussian kernel in the Y direction
	for (int y = w; y < height - w; y++) {
		for (int x = w; x < width - w; x++) {
			R2Pixel val(0.0, 0.0, 0.0, 1.0);
			int count = 0;
			for (int i = -w; i <= w; i++) {
				val += kernel_h[count] * temp.Pixel(x, y + i);
				count++;
			}
			// val.Clamp();
			(*this).SetPixel(x, y, val);
		}
	}
}


std::vector<Feature> features;

void R2Image::
Harris(double sigma)
{
	// Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
	// Output should be 50% grey at flat regions, white at corners and black/dark near edges

	R2Image img1(*this);
	R2Image img2(*this);
	R2Image img3(*this);

	img1.SobelX();
	img2.SobelY();

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			img3.Pixel(x, y) = img1.Pixel(x, y) * img2.Pixel(x, y);
		}
	}

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			img1.Pixel(x, y) *= img1.Pixel(x, y);
			img2.Pixel(x, y) *= img2.Pixel(x, y);
		}
	}

	// Using Gaussian blur
	img1.Blur(sigma);
	img2.Blur(sigma);
	img3.Blur(sigma);

	R2Image final(*this);
	R2Pixel val(0.5, 0.5, 0.5, 1.0);

	// Compute Harris scores
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			final.Pixel(x, y) = img1.Pixel(x, y) * img2.Pixel(x, y) - img3.Pixel(x, y) * img3.Pixel(x, y) -
				0.04 * ((img1.Pixel(x, y) + img2.Pixel(x, y)) * (img1.Pixel(x, y) + img2.Pixel(x, y)));
			final.Pixel(x, y) *= 150.0;
			final.Pixel(x, y) += val;
			final.Pixel(x, y).Clamp();
		}
	}
	*this = final;

	// Sort Harris scores 
	std::vector<Feature> featureVec;

	for (int y = 0; y < height; y+= 10) {
		for (int x = 0; x < width; x+=10) {
			featureVec.push_back(Feature(x, y, Pixel(x,y)));
		}
	}
	std::sort(featureVec.begin(), featureVec.end());

	int w = 10;

	// Pick 150 with the highest scores and draw a box around each one 
	for (int i = featureVec.size() - 150; i < featureVec.size(); i++) {
		features.push_back(featureVec[i]);
		for (int j = 0; j <= w; j++) {
			Pixel(featureVec[i].centerX - w/2 + j, featureVec[i].centerY - w/2) = R2red_pixel;
			Pixel(featureVec[i].centerX - w/2, featureVec[i].centerY - w/2 + j) = R2red_pixel;
			Pixel(featureVec[i].centerX + w/2, featureVec[i].centerY - w/2 + j) = R2red_pixel;
			Pixel(featureVec[i].centerX - w/2 + j, featureVec[i].centerY + w/2) = R2red_pixel;
		}
	}
}


void R2Image::
Sharpen()
{
	// Sharpen an image using a linear filter. Use a kernel of your choosing.
	int sharpen[] = { 0, -1, 0 , -1, 5, -1 , 0, -1, 0 };
	R2Image temp(*this);

	for (int x = 1; x < width - 1; x++) {
		for (int y = 1; y < height - 1; y++) {

			int count = 0;
			double rval = 0.0;
			double gval = 0.0;
			double bval = 0.0;
			for (int lx = -1; lx <= 1; lx++) {
				for (int ly = -1; ly <= 1; ly++) {
					rval += Pixel(x + lx, y + ly).Red() * sharpen[count];
					gval += Pixel(x + lx, y + ly).Green() * sharpen[count];
					bval += Pixel(x + lx, y + ly).Blue() * sharpen[count];
					count++;
				}
			}

			if (rval > 255) {
				rval = 255;
			}
			if (gval > 255) {
				gval = 255;
			}
			if (bval > 255) {
				bval = 255;
			}
			if (rval < 0) {
				rval = 0;
			}
			if (gval < 0) {
				gval = 0;
			}
			if (bval < 0) {
				bval = 0;
			}

			R2Pixel newPixel(rval, gval, bval, 1.0);
			temp.SetPixel(x, y, newPixel);
		}
	}
	(*this) = temp;
}

std::vector <Feature> maybeinliers;
std::vector <Feature> basepoints;

void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching translation (pixel precision is OK), and blend the translated "otherImage" 
	// into this image with a 50% opacity.
	// fprintf(stderr, "fit other image using translation and blend imageB over imageA\n");


	R2Image A(*this);
	R2Image B(*otherImage);
	double sigma = 2.0;
	double Sx = width / 5;
	double Sy = height / 5;


	// detect 150 features on image A
	A.Harris(sigma);
	basepoints = features;

	std::vector <Feature> good_matches;

	// Runs a local search & compute SSD at each location
	for (int i = 0; i < features.size(); i++) {
		R2Pixel ssd;
		double featureX = features[i].centerX;
		double featureY = features[i].centerY;
		double Fx = featureX;
		double Fy = featureY;
		double w = width * 0.2;
		double h = height * 0.2;
		double score_min = DBL_MAX;


		for (double x = featureX - w; x <= featureX + w; x++) {
			for (double y = featureY - h; y <= featureY + h; y++) {

				if (x > 0 && y > 0 && x < width && y < height) {
					double score = 0;

					for (int offX = -3; offX <= 3; offX++) {
						for (int offY = -3; offY <= 3; offY++) {

							if ((x + offX) > 0 && (y + offY) > 0 && (featureX + offX) > 0 && (featureY + offY) > 0
								&& (x + offX) < width && (y + offY) < height && (featureX + offX) < width & (featureY + offY) < height) {

								ssd = (*this).Pixel(featureX + offX, featureY + offY) - (*otherImage).Pixel(x + offX, y + offY);
								ssd *= ssd;
								score += ssd.Red() + ssd.Green() + ssd.Blue();
							}
							else {
								// edge cases 
								score = score_min;
							}
						}
					}

					if (score < score_min) {
						// update current min
						score_min = score;
						Fx = x;
						Fy = y;
					}
				}

			}
		}

		good_matches.push_back(Feature(Fx, Fy, Pixel(Fx, Fy)));
	}

	// save the good tracks
	maybeinliers = good_matches;

	// calling the draw line function
	/*for (int i = 0; i < features.size(); i++) {
		Feature f = features[i];
		line(f.centerX, good_matches[i].centerX, f.centerY, good_matches[i].centerY, 0.0, 0.0, 1.0);
	}*/

	return;
}

void R2Image::
line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
	if (x0>x1)
	{
		int x = y1;
		y1 = y0;
		y0 = x;

		x = x1;
		x1 = x0;
		x0 = x;
	}
	int deltax = x1 - x0;
	int deltay = y1 - y0;
	float error = 0;
	float deltaerr = 0.0;
	if (deltax != 0) deltaerr = fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
																		// note that this division needs to be done in a way that preserves the fractional part
	int y = y0;
	for (int x = x0; x <= x1; x++)
	{
		Pixel(x, y).Reset(r, g, b, 1.0);
		error = error + deltaerr;
		if (error >= 0.5)
		{
			if (deltay>0) y = y + 1;
			else y = y - 1;

			error = error - 1.0;
		}
	}
	if (x0>3 && x0<width - 3 && y0>3 && y0<height - 3)
	{
		for (int x = x0 - 3; x <= x0 + 3; x++)
		{
			for (int y = y0 - 3; y <= y0 + 3; y++)
			{
				Pixel(x, y).Reset(r, g, b, 1.0);
			}
		}
	}
}

void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.

	(*this).blendOtherImageTranslated(otherImage);

	int k = 1000;		
	double thresh = 2;

	std::vector <Feature> bestinliers;
	std::vector <Feature> bestoutliers;
	std::vector <Feature> withinliers;
	std::vector <Feature> withoutliers;

	double** bestH = NULL;
	
	double maxInliers = 0;
	int count = 0;

	while (count <= k) {
		double** matrixH = dmatrix(1, 3, 1, 3);

		double *xcoords = (double*)malloc(8 * sizeof(double));
		double *ycoords = (double*)malloc(8 * sizeof(double));

		// pick 4 random points
		int randomIndex1 = rand() % 150;
		int randomIndex2 = rand() % 150;
		int randomIndex3 = rand() % 150;
		int randomIndex4 = rand() % 150;

		xcoords[0] = basepoints[randomIndex1].centerX;
		xcoords[1] = maybeinliers[randomIndex1].centerX;
		xcoords[2] = basepoints[randomIndex2].centerX;
		xcoords[3] = maybeinliers[randomIndex2].centerX;
		xcoords[4] = basepoints[randomIndex3].centerX;
		xcoords[5] = maybeinliers[randomIndex3].centerX;
		xcoords[6] = basepoints[randomIndex4].centerX;
		xcoords[7] = maybeinliers[randomIndex4].centerX;

		ycoords[0] = basepoints[randomIndex1].centerY;
		ycoords[1] = maybeinliers[randomIndex1].centerY;
		ycoords[2] = basepoints[randomIndex2].centerY;
		ycoords[3] = maybeinliers[randomIndex2].centerY;
		ycoords[4] = basepoints[randomIndex3].centerY;
		ycoords[5] = maybeinliers[randomIndex3].centerY;
		ycoords[6] = basepoints[randomIndex4].centerY;
		ycoords[7] = maybeinliers[randomIndex4].centerY;


		// use the SVD to compute the homography matrix
		matrixH = svdTest(xcoords, ycoords, 4);

		double currentInliers = 0;

		std::vector <Feature> inliers;
		std::vector <Feature> outliers;
		std::vector <Feature> baseinliers;
		std::vector <Feature> baseoutliers;

		for (int i = 0; i < 150; i++) {
			// compute HA
			double newX = matrixH[1][1] * basepoints[i].centerX + matrixH[1][2] * basepoints[i].centerY + matrixH[1][3];
			double newY = matrixH[2][1] * basepoints[i].centerX + matrixH[2][2] * basepoints[i].centerY + matrixH[2][3];
			double newZ = matrixH[3][1] * basepoints[i].centerX + matrixH[3][2] * basepoints[i].centerY + matrixH[3][3];

			// normalize the new point
			double normX = newX / newZ;
			double normY = newY / newZ;

			// compute the vector from A to the new point
			double thisX = basepoints[i].centerX - normX;
			double thisY = basepoints[i].centerY - normY;

			// check if HA is close enough to B
			double otherX = basepoints[i].centerX - maybeinliers[i].centerX;
			double otherY = basepoints[i].centerY - maybeinliers[i].centerY;

			double diffX = thisX - otherX;
			double diffY = thisY - otherY;
			double error = sqrt(diffX * diffX + diffY * diffY);

			// fprintf(stderr, "error: %g \n", error);

			if (error < thresh) {
				currentInliers++;
				inliers.push_back(maybeinliers[i]);	
				baseinliers.push_back(basepoints[i]);
			}
			else {
				outliers.push_back(maybeinliers[i]);
				baseoutliers.push_back(basepoints[i]);
			}
		}

		if (currentInliers > maxInliers) {
			bestH = matrixH;
			maxInliers = currentInliers;

			bestinliers = inliers;
			bestoutliers = outliers;

			withinliers = baseinliers;
			withoutliers = baseoutliers;
		}
		count++;
	}

	// Recalculate H -- this currently gives incorrect results 

	/*int n = bestinliers.size();
	fprintf(stderr, "size of inlier set: %d \n", n);
	f
	double *xcoords = (double*)malloc(2 * n * sizeof(double));
	double *ycoords = (double*)malloc(2 * n * sizeof(double));

	for (int i = 0; i < 2 * n; i+=2) {
		xcoords[i + 1] = withinliers[i].centerX;
		xcoords[i] = bestinliers[i].centerX;
	
		ycoords[i + 1] = withinliers[i].centerY;
		ycoords[i] = bestinliers[i].centerY;
	}

	bestH = svdTest(xcoords, ycoords, n);*/


	R2Image warpedImage(*this);

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {

			double newX = bestH[1][1] * x + bestH[1][2] * y + bestH[1][3];
			double newY = bestH[2][1] * x + bestH[2][2] * y + bestH[2][3];
			double newZ = bestH[3][1] * x + bestH[3][2] * y + bestH[3][3]; 

			double normX = newX / newZ;
			double normY = newY / newZ;

			if (normX < width && normY < height && normX > 0 && normY > 0) {
				R2Pixel value = ((otherImage ->Pixel(normX, normY)) + warpedImage.Pixel(x, y)) / 2;
				SetPixel(x, y, value);
			}
			else {
				SetPixel(x, y, (R2Pixel(0, 0, 0, 0) + warpedImage.Pixel(x, y)) / 2);
			}
		}
	}

	return;
}

R2Point R2Image::
detectFeature(R2Image *featureImage) {

	// get size of the feature image
	double w = featureImage->width;
	double h = featureImage->height;

	int minX = 0;
	int minY = 0;

	for (int y = h; y < height - h; y++) {
		for (int x = w; x <= width - w; x++) {

			// Runs a local search & compute SSD at each location
			R2Pixel ssd;
			double score_min = DBL_MAX;

			if (x > 0 && y > 0 && x < width && y < height) {
				fprintf(stderr, "entering loop... \n");
				double score = 0;

				for (int offX = -3; offX <= 3; offX++) {
					for (int offY = -3; offY <= 3; offY++) {

						// check upper bounds (of the current image)
						if ((x + offX) > 0 && (y + offY) > 0 && (x + offX) < width && (y + offY) < height)
						{
							// check lower bounds (of the corner image)
							if ((w + offX) > 0 && (h + offY) > 0 && (w + offX) < w && (h + offY) < h) 
							{
								ssd = (*this).Pixel(x + offX, y + offY) - (*featureImage).Pixel(w + offX, h + offY);
								ssd *= ssd;
								score += ssd.Red() + ssd.Green() + ssd.Blue();
								fprintf(stderr, "printing ssd score: %g \n", score);
							}
							else
							{
								fprintf(stderr, "out of (lower) bounds \n");
							}
						}
						else 
						{
							fprintf(stderr, "out of (upper) bounds \n");
							score = score_min;	// make sure we don't update score_min
						}
					}
				}

				if (score < score_min) {
					// update current min
					fprintf(stderr, "updating current min \n");
					score_min = score;
					minX = x;
					minY = y;
				}
			}
		}
	}
	

	// draw a box around the feature -- testing 

	double s = 2000;
	for (int j = 0; j <= w; j++) {
		(*this).Pixel(minX - s / 2 + j, minY - s / 2) = R2red_pixel;
		(*this).Pixel(minX - s / 2, minY - s / 2 + j) = R2red_pixel;
		(*this).Pixel(minX + s / 2, minY - s / 2 + j) = R2red_pixel;
		(*this).Pixel(minX - s / 2 + j, minY + s / 2) = R2red_pixel;
	}

	return R2Point(minX, minY);
}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
	// Initialize everything
	if (pixels) { delete[] pixels; pixels = NULL; }
	npixels = width = height = 0;

	// Parse input filename extension
	char *input_extension;
	if (!(input_extension = (char*)strrchr(filename, '.'))) {
		fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
		return 0;
	}

	// Read file of appropriate type
	if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
	else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
	else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
	else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);

	// Should never get here
	fprintf(stderr, "Unrecognized image file extension");
	return 0;
}



int R2Image::
Write(const char *filename) const
{
	// Parse input filename extension
	char *input_extension;
	if (!(input_extension = (char*)strrchr(filename, '.'))) {
		fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
		return 0;
	}

	// Write file of appropriate type
	if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
	else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
	else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
	else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

	// Should never get here
	fprintf(stderr, "Unrecognized image file extension");
	return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
	unsigned short int bfType;
	unsigned int bfSize;
	unsigned short int bfReserved1;
	unsigned short int bfReserved2;
	unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
	unsigned int biSize;
	int biWidth;
	int biHeight;
	unsigned short int biPlanes;
	unsigned short int biBitCount;
	unsigned int biCompression;
	unsigned int biSizeImage;
	int biXPelsPerMeter;
	int biYPelsPerMeter;
	unsigned int biClrUsed;
	unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
	unsigned char rgbtBlue;
	unsigned char rgbtGreen;
	unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
	unsigned char rgbBlue;
	unsigned char rgbGreen;
	unsigned char rgbRed;
	unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
	// Read a unsigned short int from a file in little endian format 
	unsigned short int lsb, msb;
	lsb = getc(fp);
	msb = getc(fp);
	return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
	// Write a unsigned short int to a file in little endian format
	unsigned char lsb = (unsigned char)(x & 0x00FF); putc(lsb, fp);
	unsigned char msb = (unsigned char)(x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
	// Read a unsigned int word from a file in little endian format 
	unsigned int b1 = getc(fp);
	unsigned int b2 = getc(fp);
	unsigned int b3 = getc(fp);
	unsigned int b4 = getc(fp);
	return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
	// Write a unsigned int to a file in little endian format 
	unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
	unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
	unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
	unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
	// Read a int word from a file in little endian format 
	int b1 = getc(fp);
	int b2 = getc(fp);
	int b3 = getc(fp);
	int b4 = getc(fp);
	return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
	// Write a int to a file in little endian format 
	char b1 = (x & 0x000000FF); putc(b1, fp);
	char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
	char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
	char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
	// Open file
	FILE *fp = fopen(filename, "rb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s\n", filename);
		return 0;
	}

	/* Read file header */
	BITMAPFILEHEADER bmfh;
	bmfh.bfType = WordReadLE(fp);
	bmfh.bfSize = DWordReadLE(fp);
	bmfh.bfReserved1 = WordReadLE(fp);
	bmfh.bfReserved2 = WordReadLE(fp);
	bmfh.bfOffBits = DWordReadLE(fp);

	/* Check file header */
	assert(bmfh.bfType == BMP_BF_TYPE);
	/* ignore bmfh.bfSize */
	/* ignore bmfh.bfReserved1 */
	/* ignore bmfh.bfReserved2 */
	assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);

	/* Read info header */
	BITMAPINFOHEADER bmih;
	bmih.biSize = DWordReadLE(fp);
	bmih.biWidth = LongReadLE(fp);
	bmih.biHeight = LongReadLE(fp);
	bmih.biPlanes = WordReadLE(fp);
	bmih.biBitCount = WordReadLE(fp);
	bmih.biCompression = DWordReadLE(fp);
	bmih.biSizeImage = DWordReadLE(fp);
	bmih.biXPelsPerMeter = LongReadLE(fp);
	bmih.biYPelsPerMeter = LongReadLE(fp);
	bmih.biClrUsed = DWordReadLE(fp);
	bmih.biClrImportant = DWordReadLE(fp);

	// Check info header 
	assert(bmih.biSize == BMP_BI_SIZE);
	assert(bmih.biWidth > 0);
	assert(bmih.biHeight > 0);
	assert(bmih.biPlanes == 1);
	assert(bmih.biBitCount == 24);  /* RGB */
	assert(bmih.biCompression == BI_RGB);   /* RGB */
	int lineLength = bmih.biWidth * 3;  /* RGB */
	if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
	assert(bmih.biSizeImage == (unsigned int)lineLength * (unsigned int)bmih.biHeight);

	// Assign width, height, and number of pixels
	width = bmih.biWidth;
	height = bmih.biHeight;
	npixels = width * height;

	// Allocate unsigned char buffer for reading pixels
	int rowsize = 3 * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
	int nbytes = bmih.biSizeImage;
	unsigned char *buffer = new unsigned char[nbytes];
	if (!buffer) {
		fprintf(stderr, "Unable to allocate temporary memory for BMP file");
		fclose(fp);
		return 0;
	}

	// Read buffer 
	fseek(fp, (long)bmfh.bfOffBits, SEEK_SET);
	if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
		fprintf(stderr, "Error while reading BMP file %s", filename);
		return 0;
	}

	// Close file
	fclose(fp);

	// Allocate pixels for image
	pixels = new R2Pixel[width * height];
	if (!pixels) {
		fprintf(stderr, "Unable to allocate memory for BMP file");
		fclose(fp);
		return 0;
	}

	// Assign pixels
	for (int j = 0; j < height; j++) {
		unsigned char *p = &buffer[j * rowsize];
		for (int i = 0; i < width; i++) {
			double b = (double) *(p++) / 255;
			double g = (double) *(p++) / 255;
			double r = (double) *(p++) / 255;
			R2Pixel pixel(r, g, b, 1);
			SetPixel(i, j, pixel);
		}
	}

	// Free unsigned char buffer for reading pixels
	delete[] buffer;

	// Return success
	return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
	// Open file
	FILE *fp = fopen(filename, "wb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s\n", filename);
		return 0;
	}

	// Compute number of bytes in row
	int rowsize = 3 * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

	// Write file header 
	BITMAPFILEHEADER bmfh;
	bmfh.bfType = BMP_BF_TYPE;
	bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
	bmfh.bfReserved1 = 0;
	bmfh.bfReserved2 = 0;
	bmfh.bfOffBits = BMP_BF_OFF_BITS;
	WordWriteLE(bmfh.bfType, fp);
	DWordWriteLE(bmfh.bfSize, fp);
	WordWriteLE(bmfh.bfReserved1, fp);
	WordWriteLE(bmfh.bfReserved2, fp);
	DWordWriteLE(bmfh.bfOffBits, fp);

	// Write info header 
	BITMAPINFOHEADER bmih;
	bmih.biSize = BMP_BI_SIZE;
	bmih.biWidth = width;
	bmih.biHeight = height;
	bmih.biPlanes = 1;
	bmih.biBitCount = 24;       /* RGB */
	bmih.biCompression = BI_RGB;    /* RGB */
	bmih.biSizeImage = rowsize * (unsigned int)bmih.biHeight;  /* RGB */
	bmih.biXPelsPerMeter = 2925;
	bmih.biYPelsPerMeter = 2925;
	bmih.biClrUsed = 0;
	bmih.biClrImportant = 0;
	DWordWriteLE(bmih.biSize, fp);
	LongWriteLE(bmih.biWidth, fp);
	LongWriteLE(bmih.biHeight, fp);
	WordWriteLE(bmih.biPlanes, fp);
	WordWriteLE(bmih.biBitCount, fp);
	DWordWriteLE(bmih.biCompression, fp);
	DWordWriteLE(bmih.biSizeImage, fp);
	LongWriteLE(bmih.biXPelsPerMeter, fp);
	LongWriteLE(bmih.biYPelsPerMeter, fp);
	DWordWriteLE(bmih.biClrUsed, fp);
	DWordWriteLE(bmih.biClrImportant, fp);

	// Write image, swapping blue and red in each pixel
	int pad = rowsize - width * 3;
	for (int j = 0; j < height; j++) {
		for (int i = 0; i < width; i++) {
			const R2Pixel& pixel = (*this)[i][j];
			double r = 255.0 * pixel.Red();
			double g = 255.0 * pixel.Green();
			double b = 255.0 * pixel.Blue();
			if (r >= 255) r = 255;
			if (g >= 255) g = 255;
			if (b >= 255) b = 255;
			fputc((unsigned char)b, fp);
			fputc((unsigned char)g, fp);
			fputc((unsigned char)r, fp);
		}

		// Pad row
		for (int i = 0; i < pad; i++) fputc(0, fp);
	}

	// Close file
	fclose(fp);

	// Return success
	return 1;
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
	// Open file
	FILE *fp = fopen(filename, "rb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s\n", filename);
		return 0;
	}

	// Read PPM file magic identifier
	char buffer[128];
	if (!fgets(buffer, 128, fp)) {
		fprintf(stderr, "Unable to read magic id in PPM file");
		fclose(fp);
		return 0;
	}

	// skip comments
	int c = getc(fp);
	while (c == '#') {
		while (c != '\n') c = getc(fp);
		c = getc(fp);
	}
	ungetc(c, fp);

	// Read width and height
	if (fscanf(fp, "%d%d", &width, &height) != 2) {
		fprintf(stderr, "Unable to read width and height in PPM file");
		fclose(fp);
		return 0;
	}

	// Read max value
	double max_value;
	if (fscanf(fp, "%lf", &max_value) != 1) {
		fprintf(stderr, "Unable to read max_value in PPM file");
		fclose(fp);
		return 0;
	}

	// Allocate image pixels
	pixels = new R2Pixel[width * height];
	if (!pixels) {
		fprintf(stderr, "Unable to allocate memory for PPM file");
		fclose(fp);
		return 0;
	}

	// Check if raw or ascii file
	if (!strcmp(buffer, "P6\n")) {
		// Read up to one character of whitespace (\n) after max_value
		int c = getc(fp);
		if (!isspace(c)) putc(c, fp);

		// Read raw image data 
		// First ppm pixel is top-left, so read in opposite scan-line order
		for (int j = height - 1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				double r = (double)getc(fp) / max_value;
				double g = (double)getc(fp) / max_value;
				double b = (double)getc(fp) / max_value;
				R2Pixel pixel(r, g, b, 1);
				SetPixel(i, j, pixel);
			}
		}
	}
	else {
		// Read asci image data 
		// First ppm pixel is top-left, so read in opposite scan-line order
		for (int j = height - 1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				// Read pixel values
				int red, green, blue;
				if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
					fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
					fclose(fp);
					return 0;
				}

				// Assign pixel values
				double r = (double)red / max_value;
				double g = (double)green / max_value;
				double b = (double)blue / max_value;
				R2Pixel pixel(r, g, b, 1);
				SetPixel(i, j, pixel);
			}
		}
	}

	// Close file
	fclose(fp);

	// Return success
	return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
	// Check type
	if (ascii) {
		// Open file
		FILE *fp = fopen(filename, "w");
		if (!fp) {
			fprintf(stderr, "Unable to open image file: %s\n", filename);
			return 0;
		}

		// Print PPM image file 
		// First ppm pixel is top-left, so write in opposite scan-line order
		fprintf(fp, "P3\n");
		fprintf(fp, "%d %d\n", width, height);
		fprintf(fp, "255\n");
		for (int j = height - 1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				const R2Pixel& p = (*this)[i][j];
				int r = (int)(255 * p.Red());
				int g = (int)(255 * p.Green());
				int b = (int)(255 * p.Blue());
				fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
				if (((i + 1) % 4) == 0) fprintf(fp, "\n");
			}
			if ((width % 4) != 0) fprintf(fp, "\n");
		}
		fprintf(fp, "\n");

		// Close file
		fclose(fp);
	}
	else {
		// Open file
		FILE *fp = fopen(filename, "wb");
		if (!fp) {
			fprintf(stderr, "Unable to open image file: %s\n", filename);
			return 0;
		}

		// Print PPM image file 
		// First ppm pixel is top-left, so write in opposite scan-line order
		fprintf(fp, "P6\n");
		fprintf(fp, "%d %d\n", width, height);
		fprintf(fp, "255\n");
		for (int j = height - 1; j >= 0; j--) {
			for (int i = 0; i < width; i++) {
				const R2Pixel& p = (*this)[i][j];
				int r = (int)(255 * p.Red());
				int g = (int)(255 * p.Green());
				int b = (int)(255 * p.Blue());
				fprintf(fp, "%c%c%c", r, g, b);
			}
		}

		// Close file
		fclose(fp);
	}

	// Return success
	return 1;
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
extern "C" {
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
};
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
	// Open file
	FILE *fp = fopen(filename, "rb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s\n", filename);
		return 0;
	}

	// Initialize decompression info
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, fp);
	jpeg_read_header(&cinfo, TRUE);
	jpeg_start_decompress(&cinfo);

	// Remember image attributes
	width = cinfo.output_width;
	height = cinfo.output_height;
	npixels = width * height;
	int ncomponents = cinfo.output_components;

	// Allocate pixels for image
	pixels = new R2Pixel[npixels];
	if (!pixels) {
		fprintf(stderr, "Unable to allocate memory for BMP file");
		fclose(fp);
		return 0;
	}

	// Allocate unsigned char buffer for reading image
	int rowsize = ncomponents * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
	int nbytes = rowsize * height;
	unsigned char *buffer = new unsigned char[nbytes];
	if (!buffer) {
		fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
		fclose(fp);
		return 0;
	}

	// Read scan lines 
	// First jpeg pixel is top-left, so read pixels in opposite scan-line order
	while (cinfo.output_scanline < cinfo.output_height) {
		int scanline = cinfo.output_height - cinfo.output_scanline - 1;
		unsigned char *row_pointer = &buffer[scanline * rowsize];
		jpeg_read_scanlines(&cinfo, &row_pointer, 1);
	}

	// Free everything
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);

	// Close file
	fclose(fp);

	// Assign pixels
	for (int j = 0; j < height; j++) {
		unsigned char *p = &buffer[j * rowsize];
		for (int i = 0; i < width; i++) {
			double r, g, b, a;
			if (ncomponents == 1) {
				r = g = b = (double) *(p++) / 255;
				a = 1;
			}
			else if (ncomponents == 1) {
				r = g = b = (double) *(p++) / 255;
				a = 1;
				p++;
			}
			else if (ncomponents == 3) {
				r = (double) *(p++) / 255;
				g = (double) *(p++) / 255;
				b = (double) *(p++) / 255;
				a = 1;
			}
			else if (ncomponents == 4) {
				r = (double) *(p++) / 255;
				g = (double) *(p++) / 255;
				b = (double) *(p++) / 255;
				a = (double) *(p++) / 255;
			}
			else {
				fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
				return 0;
			}
			R2Pixel pixel(r, g, b, a);
			SetPixel(i, j, pixel);
		}
	}

	// Free unsigned char buffer for reading pixels
	delete[] buffer;

	// Return success
	return 1;
#else
	fprintf(stderr, "JPEG not supported");
	return 0;
#endif
}




int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
	// Open file
	FILE *fp = fopen(filename, "wb");
	if (!fp) {
		fprintf(stderr, "Unable to open image file: %s\n", filename);
		return 0;
	}

	// Initialize compression info
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_stdio_dest(&cinfo, fp);
	cinfo.image_width = width; 	/* image width and height, in pixels */
	cinfo.image_height = height;
	cinfo.input_components = 3;		/* # of color components per pixel */
	cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
	cinfo.dct_method = JDCT_ISLOW;
	jpeg_set_defaults(&cinfo);
	cinfo.optimize_coding = TRUE;
	jpeg_set_quality(&cinfo, 100, TRUE);
	jpeg_start_compress(&cinfo, TRUE);

	// Allocate unsigned char buffer for reading image
	int rowsize = 3 * width;
	if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
	int nbytes = rowsize * height;
	unsigned char *buffer = new unsigned char[nbytes];
	if (!buffer) {
		fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
		fclose(fp);
		return 0;
	}

	// Fill buffer with pixels
	for (int j = 0; j < height; j++) {
		unsigned char *p = &buffer[j * rowsize];
		for (int i = 0; i < width; i++) {
			const R2Pixel& pixel = (*this)[i][j];
			int r = (int)(255 * pixel.Red());
			int g = (int)(255 * pixel.Green());
			int b = (int)(255 * pixel.Blue());
			if (r > 255) r = 255;
			if (g > 255) g = 255;
			if (b > 255) b = 255;
			*(p++) = r;
			*(p++) = g;
			*(p++) = b;
		}
	}



	// Output scan lines
	// First jpeg pixel is top-left, so write in opposite scan-line order
	while (cinfo.next_scanline < cinfo.image_height) {
		int scanline = cinfo.image_height - cinfo.next_scanline - 1;
		unsigned char *row_pointer = &buffer[scanline * rowsize];
		jpeg_write_scanlines(&cinfo, &row_pointer, 1);
	}

	// Free everything
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);

	// Close file
	fclose(fp);

	// Free unsigned char buffer for reading pixels
	delete[] buffer;

	// Return number of bytes written
	return 1;
#else
	fprintf(stderr, "JPEG not supported");
	return 0;
#endif
}






