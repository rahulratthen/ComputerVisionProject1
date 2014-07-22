#include <cv.h>
#include <highgui.h>
#include<math.h>
#include <iostream>
using namespace cv;
using namespace std;

double InverseGammaCorrection(double v)
{
	if(v<0.03928)
		return v/12.92;
	else
		return pow((v+0.055)/1.055,2.4);
}

double* sRGBToXYZ(double r, double g, double b)
{
	//Section 5c in Lecture notes
	double XYZ[3];
	double R=0,G=0,B=0,RPrime=0,GPrime=0,BPrime=0;
	
	RPrime = (double)r/255.0;
	GPrime = (double)g/255.0;
	BPrime = (double)b/255.0;

	R = InverseGammaCorrection(RPrime);
	G = InverseGammaCorrection(GPrime);
	B = InverseGammaCorrection(BPrime);

	XYZ[0] = 0.412453*R + 0.35758*G + 0.180423*B;
	XYZ[1] = 0.212671*R + 0.71516*G + 0.072169*B;
	XYZ[2] = 0.019334*R + 0.119193*G + 0.950227*B;

	return XYZ;
}

double* XYZToLuv(double X, double Y, double Z)
{
	double Luv[3];
	double t,uw=0,vw=0,Xw=0.95,Yw=1.0,Zw=1.09,L=0,u=0,v=0,d=0,uPrime=0,vPrime=0;

	uw = (4*Xw) / (Xw + 15 * Yw + 3 *Zw);
	vw = (9*Yw) / (Xw + 15 * Yw + 3 *Zw);

	t=Y/Yw;

	if(t>0.008856)
		L = 116 * pow(t,(double)1/3.0) - 16;
	else
		L = 903.3*t;

	if(L<0)
		L=0;
	else if(L>100)
		L=100;

	d = X + 15*Y + 3*Z;
	uPrime = 4*X/d;
	vPrime = 9*Y/d;

	Luv[0] = L;
	Luv[1] = 13*L*(uPrime - uw);
	Luv[2] = 13*L*(vPrime - vw);

	return Luv;

}

double* LuvToXYZ(double L, double u, double v)
{
	//http://www.brucelindbloom.com/index.html?Eqn_Luv_to_XYZ.html
	double XYZ[3];
	double Y,a,b,c,d,Xw=0.95,Yw=1.0,Zw=1.09,uw,vw;

	uw = (4*Xw) / (Xw + 15 * Yw + 3 *Zw);
	vw = (9*Yw) / (Xw + 15 * Yw + 3 *Zw);

	if(L > 7.9996)
		Y = pow((double)(L+16)/116,3)*Yw;
	else
		Y = L/903.3*Yw;

	a = 0.333 * (((double)52*L/(u+(13*L*uw)))-1);
	b = -5*Y;
	c = -0.333;
	d = Y * ( ((double)39*L/(v+13*L*vw))-5);

	XYZ[0] = (double)(d-b)/(a-c);
	XYZ[1] = Y;
	XYZ[2] = XYZ[0]*a + b;

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


void runOnWindow(int W1,int H1, int W2,int H2, Mat inputImage, char *outName) 
{
	int rows = inputImage.rows;
	int cols = inputImage.cols;

	vector<Mat> i_planes;
	split(inputImage, i_planes);
	Mat iB = i_planes[0];
	Mat iG = i_planes[1];
	Mat iR = i_planes[2];

	// dynamically allocate RGB arrays of size rows x cols
	int** R = new int*[rows];
	int** G = new int*[rows];
	int** B = new int*[rows];
	
	for(int i = 0 ; i < rows ; i++) 
	{
		R[i] = new int[cols];
		G[i] = new int[cols];
		B[i] = new int[cols];
	}

	for(int i = 0 ; i < rows ; i++)
	{
		for(int j = 0 ; j < cols ; j++) 
		{
			R[i][j] = iR.at<uchar>(i,j);
			G[i][j] = iG.at<uchar>(i,j);
			B[i][j] = iB.at<uchar>(i,j);
		}
	}
	
	/*
	//TEsting module
			double *t5XYZ = sRGBToXYZ(100,100,100);
			double *t5Luv = XYZToLuv(t5XYZ[0],t5XYZ[1],t5XYZ[2]);

			double *t6XYZ = LuvToXYZ2(t5Luv[0],t5Luv[1],t5Luv[2]);
			double *NewsRGB1 = XYZTosRGB(t6XYZ[0],t6XYZ[1],t6XYZ[2]);
			NewsRGB1[0];
			NewsRGB1[1];
			NewsRGB1[2];

*/

	//	   The transformation should be based on the
	//	   historgram of the pixels in the W1,W2,H1,H2 range.
	//	   The following code goes over these pixels


	//Go over the window, convert pixel to Luv, find the minimum L value
	double min = std::numeric_limits<double>::max();
	double max = std::numeric_limits<double>::min();

	for(int i = H1 ; i <= H2 ; i++) 
	{
		for(int j = W1 ; j <= W2 ; j++) 
		{
			double r = R[i][j];
			double g = G[i][j];
			double b = B[i][j];
			
			//Converting XYZ to Luv
			double *tXYZ = sRGBToXYZ(r,g,b);
			double *Luv = XYZToLuv(tXYZ[0],tXYZ[1],tXYZ[2]);

			if(Luv[0]<min)
				min = Luv[0];

			if(Luv[0]>max)
				max = Luv[0];
			
			//int gray = (int) (0.3*r + 0.6*g + 0.1*b + 0.5);
		
			//R[i][j] = G[i][j] = B[i][j] = gray;
		}
	}

	for(int i = 0 ; i < rows ; i++)
	{
		for(int j = 0 ; j < cols ; j++) 
		{
			double r = R[i][j];
			double g = G[i][j];
			double b = B[i][j];
			
			//Converting sRGB to Luv
			double *tXYZ = sRGBToXYZ(r,g,b);
			double *Luv = XYZToLuv(tXYZ[0],tXYZ[1],tXYZ[2]);
			
			//Linear scaling in L domain
			
			if(Luv[0] <= min)
				Luv[0] = 0;
			else if(Luv[0] >= max)
				Luv[0] = 100;
			else
			{
				//Perform Linear Scaling
				Luv[0] = ((Luv[0]-min)*(100)/(max-min));
			}
				

			//Converting back to sRGB
			double *t1XYZ = LuvToXYZ(Luv[0],Luv[1],Luv[2]);
			double *NewsRGB = XYZTosRGB(t1XYZ[0],t1XYZ[1],t1XYZ[2]);

			R[i][j] = NewsRGB[0];
			G[i][j] = NewsRGB[1];
			B[i][j] = NewsRGB[2];


		}
	}

	Mat oR(rows, cols, CV_8UC1);
	Mat oG(rows, cols, CV_8UC1);
	Mat oB(rows, cols, CV_8UC1);
	
	for(int i = 0 ; i < rows ; i++)
	{
		for(int j = 0 ; j < cols ; j++) 
		{
			oR.at<uchar>(i,j) = R[i][j];
			oG.at<uchar>(i,j) = G[i][j];
			oB.at<uchar>(i,j) = B[i][j];
		}
	}

	Mat o_planes[] = {oB, oG, oR};
	Mat outImage;
	merge(o_planes, 3, outImage);
  
	namedWindow("output", CV_WINDOW_AUTOSIZE);
	imshow("output", outImage);
	imwrite(outName, outImage);
}

int main(int argc, char** argv) 
{
	if(argc != 7) 
	{
		cerr << argv[0] << ": "<< "got " << argc-1 << " arguments. Expecting six: w1 h1 w2 h2 ImageIn ImageOut." << endl ;
		cerr << "Example: proj1b 0.2 0.1 0.8 0.5 fruits.jpg out.bmp" << endl;
		return(-1);
	}
	double w1 = atof(argv[1]);
	double h1 = atof(argv[2]);
	double w2 = atof(argv[3]);
	double h2 = atof(argv[4]);
	char *inputName = argv[5];
	char *outputName = argv[6];

	if(w1<0 || h1<0 || w2<=w1 || h2<=h1 || w2>1 || h2>1) 
	{
		cerr << " arguments must satisfy 0 <= w1 < w2 <= 1"<< " ,  0 <= h1 < h2 <= 1" << endl;
		return(-1);
	}

	Mat inputImage = imread(inputName, CV_LOAD_IMAGE_UNCHANGED);
	if(inputImage.empty()) 
	{
		cout <<  "Could not open or find the image " << inputName << endl;
		return(-1);
	}
 
	string windowInput("input: ");
	windowInput += inputName;
  
	namedWindow(windowInput, CV_WINDOW_AUTOSIZE);
	imshow(windowInput, inputImage);

	if(inputImage.type() != CV_8UC3) 
	{
		cout <<  inputName << " is not a standard color image  " << endl;
		return(-1);
	}

	int rows = inputImage.rows;
	int cols = inputImage.cols;
	int W1 = (int) (w1*(cols-1));
	int H1 = (int) (h1*(rows-1));
	int W2 = (int) (w2*(cols-1));
	int H2 = (int) (h2*(rows-1));

	runOnWindow(W1, H1, W2, H2, inputImage, outputName);

	waitKey(0); // Wait for a keystroke
	return(0);
}