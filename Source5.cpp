#include <cv.h>
#include <highgui.h>
#include<math.h>
#include <iostream>
using namespace cv;
using namespace std;

void runOnWindow(int W1,int H1, int W2,int H2, Mat inputImage, char *outName) 
{
	int rows = inputImage.rows;
	int cols = inputImage.cols;

	Mat LuvImage(inputImage);
	cvtColor(LuvImage,LuvImage, CV_BGR2Luv);


	vector<Mat> i_planes;
	split(LuvImage, i_planes);
	Mat iL = i_planes[0];
	Mat iu = i_planes[1];
	Mat iv = i_planes[2];

	// dynamically allocate Luv arrays of size rows x cols
	int** L = new int*[rows];
	int** u = new int*[rows];
	int** v = new int*[rows];
	
	for(int i = 0 ; i < rows ; i++) 
	{
		L[i] = new int[cols];
		u[i] = new int[cols];
		v[i] = new int[cols];
	}

	for(int i = 0 ; i < rows ; i++)
	{
		for(int j = 0 ; j < cols ; j++) 
		{
			L[i][j] = iL.at<uchar>(i,j);
			u[i][j] = iu.at<uchar>(i,j);
			v[i][j] = iv.at<uchar>(i,j);
		}
	}
	

	//Go over the window, convert pixel to Luv, find the minimum L value
	double min = std::numeric_limits<double>::max();
	double max = std::numeric_limits<double>::min();

	for(int i = H1 ; i <= H2 ; i++) 
	{
		for(int j = W1 ; j <= W2 ; j++) 
		{
			double tL = L[i][j];
			double tu = u[i][j];
			double tv = v[i][j];
			
			if(tL<min)
				min = tL;

			if(tL>max)
				max = tL;
			
			//int gray = (int) (0.3*r + 0.6*g + 0.1*b + 0.5);
		
			//R[i][j] = G[i][j] = B[i][j] = gray;
		}
	}

	for(int i = 0 ; i < rows ; i++)
	{
		for(int j = 0 ; j < cols ; j++) 
		{
			if(L[i][j] <= min)
				L[i][j] = 0;
			else if(L[i][j] >= max)
				L[i][j] = 100;
			else
			{
				//Perform Linear Scaling
				L[i][j] = ((L[i][j]-min)*(100)/(max-min));
			}
				

		}
	}

	Mat oL(rows, cols, CV_8UC1);
	Mat ou(rows, cols, CV_8UC1);
	Mat ov(rows, cols, CV_8UC1);
	
	for(int i = 0 ; i < rows ; i++)
	{
		for(int j = 0 ; j < cols ; j++) 
		{
			oL.at<uchar>(i,j) = L[i][j];
			ou.at<uchar>(i,j) = u[i][j];
			ov.at<uchar>(i,j) = v[i][j];
		}
	}

	Mat o_planes[] = {oL, ou, ov};
	Mat outLuvImage;
	merge(o_planes, 3, outLuvImage);

	//Luv to sRGB
	Mat outRGBImage(outLuvImage);
	cvtColor(outRGBImage, outRGBImage, CV_Luv2BGR);
  
	namedWindow("output", CV_WINDOW_AUTOSIZE);
	imshow("output", outRGBImage);
	imwrite(outName, outRGBImage);
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