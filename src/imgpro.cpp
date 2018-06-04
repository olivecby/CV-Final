

// Computer Vision for Digital Post-Production
// Lecturer: Gergely Vass - vassg@vassg.hu
//
// Skeleton Code for programming assigments
// 
// Code originally from Thomas Funkhouser
// main.c
// original by Wagner Correa, 1999
// modified by Robert Osada, 2000
// modified by Renato Werneck, 2003
// modified by Jason Lawrence, 2004
// modified by Jason Lawrence, 2005
// modified by Forrester Cole, 2006
// modified by Tom Funkhouser, 2007
// modified by Chris DeCoro, 2007
//



// Include files

#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"



// Program arguments

static char options[] =
"  -help\n"
"  -svdTest\n"
"  -sobelX\n"
"  -sobelY\n"
"  -log\n"
"  -harris <real:sigma>\n"
"  -saturation <real:factor>\n"
"  -brightness <real:factor>\n"
"  -blur <real:sigma>\n"
"  -sharpen \n"
"  -matchTranslation <file:other_image>\n"
"  -matchHomography <file:other_image>\n"
"  -video\n";


static void 
ShowUsage(void)
{
  // Print usage message and exit
  fprintf(stderr, "Usage: imgpro input_image output_image [  -option [arg ...] ...]\n");
  fprintf(stderr, options);
  exit(EXIT_FAILURE);
}



static void 
CheckOption(char *option, int argc, int minargc)
{
  // Check if there are enough remaining arguments for option
  if (argc < minargc)  {
    fprintf(stderr, "Too few arguments for %s\n", option);
    ShowUsage();
    exit(-1);
  }
}



static int 
ReadCorrespondences(char *filename, R2Segment *&source_segments, R2Segment *&target_segments, int& nsegments)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open correspondences file %s\n", filename);
    exit(-1);
  }

  // Read number of segments
  if (fscanf(fp, "%d", &nsegments) != 1) {
    fprintf(stderr, "Unable to read correspondences file %s\n", filename);
    exit(-1);
  }

  // Allocate arrays for segments
  source_segments = new R2Segment [ nsegments ];
  target_segments = new R2Segment [ nsegments ];
  if (!source_segments || !target_segments) {
    fprintf(stderr, "Unable to allocate correspondence segments for %s\n", filename);
    exit(-1);
  }

  // Read segments
  for (int i = 0; i <  nsegments; i++) {

    // Read source segment
    double sx1, sy1, sx2, sy2;
    if (fscanf(fp, "%lf%lf%lf%lf", &sx1, &sy1, &sx2, &sy2) != 4) { 
      fprintf(stderr, "Error reading correspondence %d out of %d\n", i, nsegments);
      exit(-1);
    }

    // Read target segment
    double tx1, ty1, tx2, ty2;
    if (fscanf(fp, "%lf%lf%lf%lf", &tx1, &ty1, &tx2, &ty2) != 4) { 
      fprintf(stderr, "Error reading correspondence %d out of %d\n", i, nsegments);
      exit(-1);
    }

    // Add segments to list
    source_segments[i] = R2Segment(sx1, sy1, sx2, sy2);
    target_segments[i] = R2Segment(tx1, ty1, tx2, ty2);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int 
main(int argc, char **argv)
{
  // Look for help
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "-help")) {
      ShowUsage();
    }
  /*if (!strcmp(argv[i], "-svdTest")) {
      R2Image *image = new R2Image();
    image->svdTest();
    return 0;
    }*/
	else if (!strcmp(argv[i], "-video")) {
		printf("Video processing started\n");

		char inputName[100] = "../videoinput/v2p%02d.jpg";
		char outputName[100] = "../videooutput/output%02d.jpg";

		R2Image *mainImage = new R2Image();
		char currentFilename[100];
		char currentOutputFilename[100];

		if (!mainImage) {
			fprintf(stderr, "Unable to allocate image\n");
			exit(-1);
		}


		// read very first frame

		sprintf(currentFilename, inputName, 1);
		printf("Filename: %s\n", currentFilename);
		if (!mainImage->Read(currentFilename)) {
			fprintf(stderr, "Unable to read first image\n");
			exit(-1);
		}

		


		//--------------------------read markers---------------------------

		char *upperLeftName = "../markersv5/upperleft.jpg";
		char *upperRightName = "../markersv5/upperright.jpg";
		char *lowerLeftName = "../markersv5/lowerleft.jpg";
		char *lowerRightName = "../markersv5/lowerright.jpg";

		R2Image *upperLeftMarkerImage = new R2Image();
		R2Image *upperRightMarkerImage = new R2Image();
		R2Image *lowerLeftMarkerImage = new R2Image();
		R2Image *lowerRightMarkerImage = new R2Image();


		if (!upperLeftMarkerImage) {
			fprintf(stderr, "Unable to allocate upper left marker image\n");
			exit(-1);
		}
		if (!upperLeftMarkerImage->Read(upperLeftName)) {
			fprintf(stderr, "Unable to read upperleft marker image file\n");
			exit(-1);
		}

		if (!upperRightMarkerImage) {
			fprintf(stderr, "Unable to allocate upper right marker image\n");
			exit(-1);
		}
		if (!upperRightMarkerImage->Read(upperRightName)) {
			fprintf(stderr, "Unable to read upperright marker image file\n");
			exit(-1);
		}

		if (!lowerRightMarkerImage) {
			fprintf(stderr, "Unable to allocate lower right marker image\n");
			exit(-1);
		}
		if (!lowerRightMarkerImage->Read(lowerRightName)) {
			fprintf(stderr, "Unable to read upperright marker image file\n");
			exit(-1);
		}

		if (!lowerLeftMarkerImage) {
			fprintf(stderr, "Unable to allocate lower left marker image\n");
			exit(-1);
		}
		if (!lowerLeftMarkerImage->Read(lowerLeftName)) {
			fprintf(stderr, "Unable to read upperleft marker image file\n");
			exit(-1);
		}


    char replacementFileName[100] = "../replacement/v1p%02d.jpg";


		// =============== VIDEO PROCESSING ===============

		//mainImage->Blur(3.0f);
		// here you could call mainImage->FirstFrameProcessing( ); 

    R2Point prevUpperLeftMarker(0,0);
    R2Point prevUpperRightMarker(0,0);
    R2Point prevLowerRightMarker(0,0);
    R2Point prevLowerLeftMarker(0,0);

		int end = 110;
		for (int i = 1; i < end; i++)

		{


      // read replacement image
    char replacementName[100];
    sprintf(replacementName,replacementFileName, i);
    R2Image *replacementImage = new R2Image();

    if (!replacementImage) {
      fprintf(stderr, "Unable to allocate replacement image\n");
      exit(-1);
    }
    if (!replacementImage->Read(replacementName)) {
      fprintf(stderr, "Unable to read replacement image file\n");
      exit(-1);
    }



			R2Image *currentImage = new R2Image();
			if (!currentImage) {
				fprintf(stderr, "Unable to allocate image %d\n", i);
				exit(-1);
			}

			sprintf(currentFilename, inputName, i+9);
			sprintf(currentOutputFilename, outputName, i);

			printf("Processing file %s\n", currentFilename);
			if (!currentImage->Read(currentFilename)) {
				fprintf(stderr, "Unable to read image %d\n", i+9);
				exit(-1);
			}

			// detect 4 corner markers on current + next frame

			R2Point upperLeftMarker = currentImage->detectFeature(upperLeftMarkerImage, prevUpperLeftMarker);
			R2Point upperRightMarker = currentImage->detectFeature(upperRightMarkerImage, prevUpperRightMarker);
			R2Point lowerLeftMarker = currentImage->detectFeature(lowerLeftMarkerImage, prevLowerLeftMarker);
			R2Point lowerRightMarker = currentImage->detectFeature(lowerRightMarkerImage, prevLowerRightMarker);

      prevUpperLeftMarker = upperLeftMarker;
      prevUpperRightMarker = upperRightMarker;
      prevLowerRightMarker = lowerRightMarker;
      prevLowerLeftMarker = lowerLeftMarker;


			//fprintf(stderr, "upper left marker: (%g, %g)\n", upperLeftMarker.X(), upperLeftMarker.Y());
			// fprintf(stderr, "upper right marker: (%g, %g)\n", upperRightMarker.X(), upperRightMarker.Y());
			// fprintf(stderr, "lower left marker: (%g, %g)\n", lowerLeftMarker.X(), lowerLeftMarker.Y());
			// fprintf(stderr, "lower right marker: (%g, %g)\n", lowerRightMarker.X(), lowerRightMarker.Y());


			// get size of replacement image
			double w = replacementImage->Width();
			double h = replacementImage->Height();

			// save feature points
			double *xcoords = (double*)malloc(8 * sizeof(double));
			double *ycoords = (double*)malloc(8 * sizeof(double));


			/*

			(0,h) ---- (w,h))
			  |			 |
			  |			 |
			  |			 |
			(0,0) ---- (w,0)

			*/


			xcoords[0] = lowerLeftMarker.X();
			xcoords[1] = 0;
			xcoords[2] = lowerRightMarker.X();
			xcoords[3] = w;
			xcoords[4] = upperRightMarker.X();
			xcoords[5] = w;
			xcoords[6] = upperLeftMarker.X();
			xcoords[7] = 0;

			ycoords[0] = lowerLeftMarker.Y();
			ycoords[1] = 0;
			ycoords[2] = lowerRightMarker.Y();
			ycoords[3] = 0;
			ycoords[4] = upperRightMarker.Y();
			ycoords[5] = h;
			ycoords[6] = upperLeftMarker.Y();
			ycoords[7] = h;


			// replace frame 
			currentImage->replaceFrameContent(xcoords, ycoords, replacementImage);

      printf("**************************\n");

			// write result to file
			if (!currentImage->Write(currentOutputFilename)) {
				fprintf(stderr, "Unable to write %s\n", currentOutputFilename);
				exit(-1);
			}
			delete currentImage;
      delete replacementImage;
		}
		delete mainImage;
		
		delete upperLeftMarkerImage;
		delete upperRightMarkerImage;
		delete lowerLeftMarkerImage;
		delete lowerRightMarkerImage;

		// Return success
		return EXIT_SUCCESS;
	}
  }

  // Read input and output image filenames
  if (argc < 3)  ShowUsage();
  argv++, argc--; // First argument is program name
  char *input_image_name = *argv; argv++, argc--; 
  char *output_image_name = *argv; argv++, argc--; 

  // Allocate image
  R2Image *image = new R2Image();
  if (!image) {
    fprintf(stderr, "Unable to allocate image\n");
    exit(-1);
  }


  // Read input image
  if (!image->Read(input_image_name)) {
    fprintf(stderr, "Unable to read image from %s\n", input_image_name);
    exit(-1);
  }

  // Initialize sampling method
  int sampling_method = R2_IMAGE_POINT_SAMPLING;

  // Parse arguments and perform operations 
  while (argc > 0) {
    if (!strcmp(*argv, "-brightness")) {
      CheckOption(*argv, argc, 2);
      double factor = atof(argv[1]);
      argv += 2, argc -=2;
      image->Brighten(factor);
    }
  else if (!strcmp(*argv, "-sobelX")) {
      argv++, argc--;
      image->SobelX();
    }
  else if (!strcmp(*argv, "-sobelY")) {
      argv++, argc--;
      image->SobelY();
    }
  else if (!strcmp(*argv, "-log")) {
      argv++, argc--;
      image->LoG();
    }
    else if (!strcmp(*argv, "-saturation")) {
      CheckOption(*argv, argc, 2);
      double factor = atof(argv[1]);
      argv += 2, argc -= 2;
      image->ChangeSaturation(factor);
    }
  else if (!strcmp(*argv, "-harris")) {
      CheckOption(*argv, argc, 2);
      double sigma = atof(argv[1]);
      argv += 2, argc -= 2;
      image->Harris(sigma);
    }
    else if (!strcmp(*argv, "-blur")) {
      CheckOption(*argv, argc, 2);
      double sigma = atof(argv[1]);
      argv += 2, argc -= 2;
      image->Blur(sigma);
    }
    else if (!strcmp(*argv, "-sharpen")) {
      argv++, argc--;
      image->Sharpen();
    }
    else if (!strcmp(*argv, "-matchTranslation")) {
      CheckOption(*argv, argc, 2);
      R2Image *other_image = new R2Image(argv[1]);
      argv += 2, argc -= 2;
      image->blendOtherImageTranslated(other_image);
      delete other_image;
    }
    else if (!strcmp(*argv, "-matchHomography")) {
      CheckOption(*argv, argc, 2);
      R2Image *other_image = new R2Image(argv[1]);
      argv += 2, argc -= 2;
      image->blendOtherImageHomography(other_image);
      delete other_image;
    }
    
    else {
      // Unrecognized program argument
      fprintf(stderr, "image: invalid option: %s\n", *argv);
      ShowUsage();
    }
  }

  // Write output image
  if (!image->Write(output_image_name)) {
    fprintf(stderr, "Unable to read image from %s\n", output_image_name);
    exit(-1);
  }

  // Delete image
  delete image;

  // Return success
  return EXIT_SUCCESS;}