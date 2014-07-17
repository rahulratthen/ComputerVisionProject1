#include <cv.h>
#include <highgui.h>
#include<math.h>
#include <iostream>
using namespace cv;
using namespace std;



double* xyYToXYZ(double x,double y,double Y)
{
	double XYZ[3];

	//Convert xyY to XYZ
	double cX=0,cY=0,cZ=0;
	if(y==0)
		XYZ[0] = 0;
	else
		XYZ[0] = x*Y/y;
	XYZ[1] = Y;
	
	if(y==0)
		XYZ[2] = 0;
	else
		XYZ[2] = (1-x-y)/y*Y;

	if(y==0)
		XYZ[0] = XYZ[1] = XYZ[2] = 0;
	/*
	if(XYZ[0] < 0)
		XYZ[0] = 0;
	else if(XYZ[0] > 1.0)
		XYZ[0] = 1.0;

	if(XYZ[1] < 0)
		XYZ[1] = 0;
	else if(XYZ[1] > 1.0)
		XYZ[1] = 1.0;

	if(XYZ[2] < 0)
		XYZ[2] = 0;
	else if(XYZ[2] > 1.0)
		XYZ[2] = 1.0;
		*/
	return XYZ;
}

double* LuvToXYZ(double L, double u, double v)
{
	double XYZ[3], X=0,Y=0,Z=0;
	//Section 6.2 in Lecture notes
	double uw=0,vw=0,Xw=0.95,Yw=1.0,Zw=1.09, uPrime=0, vPrime=0;
	uw = (4*Xw) / (Xw + 15 * Yw + 3 *Zw);
	vw = (9*Yw) / (Xw + 15 * Yw + 3 *Zw);

	uPrime = (u + 13 * uw * L) / 13 * L;
	vPrime = (v + 13 * vw * L) / 13 * L;

	//Compute Y
	if(L > 7.9996)
		Y = pow((L+16)/116,3)*Yw;
	else
		Y = L/903.3*Yw;

	//Compute X,Z
	if(vPrime == 0)
	{
		X = 0;
		Z = 0;
	}
	else
	{
		X = Y * 2.25 * uPrime/vPrime;
		Z = Y * (3 - 0.75*uPrime - 5*vPrime)/vPrime;
	}

	XYZ[0] = X;
	XYZ[1] = Y;
	XYZ[2] = Z;
	/*
	if(XYZ[0] < 0)
		XYZ[0] = 0;
	else if(XYZ[0] > 1.0)
		XYZ[0] = 1.0;

	if(XYZ[1] < 0)
		XYZ[1] = 0;
	else if(XYZ[1] > 1.0)
		XYZ[1] = 1.0;

	if(XYZ[2] < 0)
		XYZ[2] = 0;
	else if(XYZ[2] > 1.0)
		XYZ[2] = 1.0;
 */
	return XYZ;
}

double GammaCorrection(double D)
{
	double result;
	if(D<0.00304)
		result = 12.92*D;
	else
		result = 1.055*pow(D,1.0/2.4) - 0.055;

	return result;
}


double* XYZTosRGB(double cX, double cY, double cZ)
{
	double sRGB[3];
	//Convert the XYZ to Linear RGB using Formula 5a in Handouts
	double lr=0,lg=0,lb=0;
	lr = (3.240479 * cX) + (-1.53715 * cY) + (-0.498535 * cZ);
	lg = (-0.969259 * cX) + (1.875991 * cY) + (0.041556 * cZ);
	lb = (0.055648 * cX) + (-0.204043 * cY) + (1.057311 * cZ);

	//Now do gamma correction to get nonLinear sRGB
	sRGB[0] = GammaCorrection(lr)*255;
	sRGB[1] = GammaCorrection(lg)*255;
	sRGB[2] = GammaCorrection(lb)*255;

	if(sRGB[0] < 0)
		sRGB[0] = 0;
	else if(sRGB[0] > 255)
		sRGB[0] = 255;
	
	if(sRGB[1] < 0)
		sRGB[1] = 0;
	else if(sRGB[1] > 255)
		sRGB[1] = 255;
	
	if(sRGB[2] < 0)
		sRGB[2] = 0;
	else if(sRGB[2] > 255)
		sRGB[2] = 255;
	
	return sRGB;
}


int main(int argc, char** argv) {
	if(argc != 3) 
	{
		cout << argv[0] << ": "<< "got " << argc-1 << " arguments. Expecting two: width height." << endl ;
		return(-1);
	}

	int width = atoi(argv[1]);
	int height = atoi(argv[2]);
	int** RED1 = new int*[height];
	int** GREEN1 = new int*[height];
	int** BLUE1 = new int*[height];
	int** RED2 = new int*[height];
	int** GREEN2 = new int*[height];
	int** BLUE2 = new int*[height];

	for(int i = 0 ; i < height ; i++) 
	{
		RED1[i] = new int[width];
		GREEN1[i] = new int[width];
		BLUE1[i] = new int[width];
		RED2[i] = new int[width];
		GREEN2[i] = new int[width];
		BLUE2[i] = new int[width];
	}
  
	for(int i = 0 ; i < height ; i++)
	{
		for(int j = 0 ; j < width ; j++)
		{
			int r1, g1, b1;
			int r2, g2, b2;

			double x = (double)j/(double)width;
			double y = (double)i/(double)height;
			double Y = 1.0;

			double L = 90;
			double u = x * 512 - 255;
			double v = y * 512 - 255;


			/* Your code should be placed here
			   It should translate xyY to byte sRGB
			   and Luv to byte sRGB
			
			r1 = (int) (x * 255);
			g1 = (int) (y * 255);
			b1 = (int) (1.0 * 255);

			r2 = (int) (1.0 * 255);
			g2 = (int) (x * 255);
			b2 = (int) (y * 255);
			*/

			//xyY -> XYZ -> nonLinear sRGB
			double *tempXYZ = xyYToXYZ(x,y,Y);
			double *converted = XYZTosRGB(tempXYZ[0],tempXYZ[1],tempXYZ[2]);

			r1 = (int) converted[0];
			g1 = (int) converted[1];
			b1 = (int) converted[2];

			//Luv -> XYZ -> nonLinear sRGB
			double *temp1XYZ = LuvToXYZ(L,u,v);
			double *converted1 = XYZTosRGB(temp1XYZ[0],temp1XYZ[1],temp1XYZ[2]);
			r2 = (int) converted1[0];
			g2 = (int) converted1[1];
			b2 = (int) converted1[2];
			// this is the end of your code

			RED1[i][j] = r1;
			GREEN1[i][j] = g1;
			BLUE1[i][j] = b1;
			RED2[i][j] = r2;
			GREEN2[i][j] = g2;
			BLUE2[i][j] = b2;
			  
		}
	}


	Mat R1(height, width, CV_8UC1);
	Mat G1(height, width, CV_8UC1);
	Mat B1(height, width, CV_8UC1);

	Mat R2(height, width, CV_8UC1);
	Mat G2(height, width, CV_8UC1);
	Mat B2(height, width, CV_8UC1);

	for(int i = 0 ; i < height ; i++)
	{
		for(int j = 0 ; j < width ; j++) 
		{

			R1.at<uchar>(i,j) = RED1[i][j];
			G1.at<uchar>(i,j) = GREEN1[i][j];
			B1.at<uchar>(i,j) = BLUE1[i][j];

			R2.at<uchar>(i,j) = RED2[i][j];
			G2.at<uchar>(i,j) = GREEN2[i][j];
			B2.at<uchar>(i,j) = BLUE2[i][j];
	    }
	}

	Mat xyY;
	Mat xyY_planes[] = {B1, G1, R1};
	merge(xyY_planes, 3, xyY);
	namedWindow("xyY",CV_WINDOW_AUTOSIZE);
	imshow("xyY", xyY);

	Mat Luv;
	Mat Luv_planes[] = {B2, G2, R2};
	merge(Luv_planes, 3, Luv);
	namedWindow("Luv",CV_WINDOW_AUTOSIZE);
	imshow("Luv", Luv);
	waitKey(0); // Wait for a keystroke
	return(0);
}