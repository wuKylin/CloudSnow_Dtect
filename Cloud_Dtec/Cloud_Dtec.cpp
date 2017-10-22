#include"Cloud_Dtec.h"
#include <windows.h>
#include <winbase.h>
#include<thread>
#include<mutex>

int GetTiffFile(const string& directory, vector<string>* entries)
{
	int count = 0;//统计tiff文件的数量
	char dir[MAX_PATH];
	string  str_dir = directory;

	strncpy(dir, str_dir.c_str(), MAX_PATH);

	DIR  *p_dir;
	struct dirent *p_dirent;

	p_dir = opendir(dir);

	while (p_dirent = readdir(p_dir))
	{
		string  str_fn = p_dirent->d_name;

		//if (str_fn != "." && str_fn != "..")  entries->push_back(str_fn);  读取所有文件
		size_t len = str_fn.length();
		if (len > 5 && str_fn.substr(len - 5, len) == ".tiff")
		{
			entries->push_back(str_fn);                //读取tiff文件
			count++;
		}
	}
	return count;

} 

void ReadImgData(byte *img_src)
{
	hasBoundaries_ = false;
	if (hasBoundaries_)
	{
		hasBoundaries_ = false;
		boundaries_->CleanData();
	}
	if (img_src != NULL == NULL)
		return;

	int width_, height_;
	width_ = width;
	height_ = height;

	byte *tpImg;

	tpImg = new byte[width*height*3];
	int i;
	for (i = 0; i < height*width*3; i++)
	{
		tpImg[i] = img_src[i];
	}

	/*                                 */
	cbgImage_ = new BgImage();
	whiteImage_ = new BgImage();
	filtImage_ = new BgImage();
	segmImage_ = new BgImage();
	boundaries_ = new BgPointSet();
	regionPts_ = new BgPointSet();

	if (iProc)
	{
		delete iProc;
		iProc = NULL;
	}
	iProc = new msImageProcessor;
	//	sigmaS		= 7;
	sigmaS = 7;
	//	sigmaR		= float(5.5);
	sigmaR = float(5.5);
	aij = float(0.8);
	epsilon = float(0.6);
	//	minRegion	= 1000;
	//	minRegion	= 400;
	minRegion = 300;
	kernelSize = 2;

	cbgImage_->SetImageFromRGB(tpImg, width_, height_, true);

	//set cbgImage
	if (tpImg)
	{
		delete[]tpImg;
	}
	tpImg = NULL;

	// 图像的显示，这四个变量有利于进行图像的显示。
	/************************************************************************/
	/*
	cbgImage_ hasImage_	;  // original image
	filtImage_ hasFilter_	;  // 平滑后的图像 filtered image
	segmImage_ hasSegment_	;  // 分割后的数据 segmented image
	whiteImage_ hasBoundaries_ ;  // 是否显示分割的轮廓 contour image
	*/
	/************************************************************************/
	hasImage_ = 1;
	hasFilter_ = 0;
	hasSegment_ = 0;
	hasBoundaries_ = 0;

	//set white image
	unsigned char *img = new unsigned char[(width_)*(height_)* 3];
	memset(img, 255, (width_)*(height_)* 3);
	whiteImage_->SetImageFromRGB(img, width_, height_, true);
	delete[] img;

	//delete edge maps if they exists...
	if (customMap_)	delete[] customMap_;
	if (gradMap_)	delete[] gradMap_;
	if (confMap_)	delete[] confMap_;
	if (weightMap_)	delete[] weightMap_;

	//set edge maps to NULL...
	customMap_ = (float *)NULL;
	gradMap_ = (float *)NULL;
	confMap_ = (float *)NULL;
	weightMap_ = (float *)NULL;

	hasImage_ = 1;
}

void GetBoundariesPtIdx()
{
	// 进行分割后
	if (hasSegment_)
	{
		//clean boundary data if any exists
		boundaries_->CleanData();
		int	width, height;
		width = cbgImage_->x_; // 原始影像的宽度
		height = cbgImage_->y_; // 原始影像的高度
		int			*boundaryIndeces = regionList->GetBoundaryIndeces(0); // 获取多边形的指数
		int			numRegions = regionList->GetNumRegions();  // 多边形的个数
		int			boundaryPointCount = 0; // 边界的点数

		//calculate the number of boundary points stored by region list
		//class
		int i;

		// 计算边界的点数，通过多边形表来存储
		for (i = 0; i < numRegions; i++)
		{
			boundaryPointCount += regionList->GetBourdaryCount(i); // 获取每个多边形的边界数，将边界的点数进行累加
		}

		//create a point set using calculated boundary point count...
		boundaries_->x_ = new int[boundaryPointCount]; // 边界的x序号
		boundaries_->y_ = new int[boundaryPointCount]; // 边界的y序号
		boundaries_->n_ = boundaryPointCount;  // 边界的点数



		for (i = 0; i < boundaryPointCount; i++)
		{
			boundaries_->x_[i] = boundaryIndeces[i] % width; // 列数
			boundaries_->y_[i] = boundaryIndeces[i] / width; // 行数
		}




		//set point type to point (= 1)
		boundaries_->type_ = 1;
		//set has boundaries to true
		hasBoundaries_ = true;
	}
	else
	{
		hasBoundaries_ = false;
		cerr << "没有打开影像，或没有进行分割！" << endl;
	}
}

void GetRegionsPtIdx()			// 获取所有区域所包含的像素
{
	if (hasFilter_)
	{
		//clean boundary data if any exists
		regionPts_->CleanData();  // 清除数据
		int	width, height;
		width = cbgImage_->x_; // 原始影像的宽度
		height = cbgImage_->y_; // 原始影像的高度

		// 问题的关键出在regionIndeces上
		int			*regionIndeces = regionList->GetRegionPtIndeces(0); // 获取多边形的指数
		int			numRegions = regionList->GetNumRegions();  // 多边形的个数
		int			regionPointCount = 0; // 所有区域内包含的总点数，其实就是整幅影像的总点数

		//calculate the number of region points stored by region list
		//class
		int i;

		regionPointCount = width*height;

		//create a point set using calculated boundary point count...
		regionPts_->x_ = new int[regionPointCount]; // 点的x序号
		regionPts_->y_ = new int[regionPointCount]; // 点的y序号

		regionPts_->n_ = regionPointCount;  // 点数

		for (i = 0; i < regionPointCount; i++)
		{
			regionPts_->x_[i] = regionIndeces[i] % width; // 列数
			regionPts_->y_[i] = regionIndeces[i] / width; // 行数
		}

		regionPts_->type_ = 1;
	}
}

void imageseg()
{
	if (hasImage_ == 0)	{ return; }   //hasImage_为0表示没有正确的完成ReadImgData()
	//obtain image dimensions
	int	width, height;
	width = cbgImage_->x_; // 原始影像的宽度
	height = cbgImage_->y_; // 原始影像的高度

	//obtain image type (color or grayscale)
	imageType 	gtype;
	if (cbgImage_->colorIm_)
		gtype = COLOR;
	else
		gtype = GRAYSCALE;

	buseWeightMap_ = true;
	if (buseWeightMap_)
	{
		//if the weight map has already been defined
		//then find out if it needs to be recomputed
		edgeParamsHaveChanged_ = false;
		// 如果参数改变且权重图已经计算，则要重置权重图，根据置信图与梯度图，进行计算权重图
		if ((weightMap_) && (edgeParamsHaveChanged_))
		{
			delete[] confMap_; // 置信度图
			delete[] gradMap_; // 梯度图
			delete[] weightMap_; // 权重图
			weightMap_ = (float *)NULL;
			//indicate that the change has been recognized...
			edgeParamsHaveChanged_ = false;
		}
		//if the weight map has not been computed or discarded
		//then recompute it...
		if (!weightMap_) // 当权重图重置或权重图为空时，对权重图进行计算
		{
			//allocate memory for gradient and confidence maps
			confMap_ = new float[width*height];
			gradMap_ = new float[width*height];

			//compute gradient and confidence maps
			BgEdgeDetect	edgeDetector(kernelSize);
			/*************************************************************/
			cerr << " 根据原始影像计算置信图和梯度图..." << endl;
			//pFrame->UpdateInfo(strInfo);
			/*************************************************************/
			edgeDetector.ComputeEdgeInfo(cbgImage_, confMap_, gradMap_);  // 根据原始影像计算置信图，梯度图（边缘检测来进行计算）

			//compute weight map...

			//allocate memory for weight map
			weightMap_ = new float[width*height];
			/*************************************************************/
			cerr << " 根据置信图和梯度图计算权重图..." << endl;
			//pFrame->UpdateInfo(strInfo);
			/*************************************************************/
			//compute weight map using gradient and confidence maps
			int i;
			for (i = 0; i<width*height; i++)
			{
				if (gradMap_[i] > 0.02)
					weightMap_[i] = aij*gradMap_[i] + (1 - aij)*confMap_[i]; //
				else
					weightMap_[i] = 0; // 如果梯度过小，则权重赋予0
			}
		}
	}
	//determine operation (filtering or segmentation)
	/******************************************************************************************************************/
	/* 此处默认操作为3号操作，但是我们希望主要是进行以下几步:
	operation	= 1 //filter
	operation	= 2 //fuse
	operation	= 3 //segment
	*/
	/************************************************************************************************/
	/*	int	operation	= 3;*/

	//create instance of image processor class
	//msImageProcessor *iProc = new msImageProcessor();

	//define an input image using the image under consideration
	//(if filtering or segmentation has taken place, then use this
	// result upon performing fusing...)
	hasFilter_ = false;
	if ((operation == 2) && (hasFilter_))
		iProc->DefineImage(filtImage_->im_, gtype, height, width); // 如果已经滤波以及融合，则对设置对平滑的图像进行处理
	else
		iProc->DefineImage(cbgImage_->im_, gtype, height, width); // 如果是没有进行平滑处理，就对原始图像进行处理

	bool useCustomMap = false;
	//set the weight map (if one was specified and a custom map is not being utilized)
	if ((buseWeightMap_) && (weightMap_) && (!useCustomMap))	iProc->SetWeightMap(weightMap_, epsilon);

	//set the custom map (if one was provided)
	if ((buseWeightMap_) && (customMap_) && (useCustomMap))	iProc->SetWeightMap(customMap_, epsilon);

	// 进行图像的分割和图像的平滑
	//perform image segmentation or filtering....
	iProc->SetSpeedThreshold(speedUpThreshold_);
	speedUpLevel_ = MED_SPEEDUP;

	//filter the image...
	/*************************************************************/
	cerr << " 正在平滑影像...\r\n" << endl;
	//pFrame->UpdateInfo(strInfo);
	/*************************************************************/
	iProc->Filter(sigmaS, sigmaR, speedUpLevel_);

	//filter the image....
	int dim;
	if (cbgImage_->colorIm_)
		dim = 3;
	else
		dim = 1;
	unsigned char *tempImage = new unsigned char[dim*height*width];
	iProc->GetResults(tempImage);

	//fuse regions...
	/*************************************************************/
	cerr << " 正在融合影像...\r\n" << endl;
	//pFrame->UpdateInfo(strInfo);
	/*************************************************************/
	// 区域生长，合并，修剪，FuseRegions就是比segment少了一个filter，分割可以看作先filter后fuseRegions
	iProc->FuseRegions(sigmaR, minRegion);


	//obtain the segmented and filtered image...
	filtImage_->Resize(width, height, cbgImage_->colorIm_);  // 
	memcpy(filtImage_->im_, tempImage, dim*height*width*sizeof(unsigned char));
	delete[] tempImage;
	segmImage_->Resize(width, height, cbgImage_->colorIm_);
	iProc->GetResults(segmImage_->im_); // 存储分割数据

	//indicate that the segmented image has been computed...
	hasFilter_ = 1;
	hasSegment_ = 1;

	//**************************************这段代码需要判断是否为NULL*********************************************
	/************************************************************************/
	/*  defineboudaries modepointcounts                                                                    */
	/************************************************************************/
	regionList = iProc->GetBoundariesAndRegions(); // 获取边界 ,
	regionLables = iProc->GetLabels(); // 获取多边形的标识
	regionPC = iProc->GetMPC(); // 获取多边形的所包含的像素个数
	regionData = iProc->GetModeData(); // 获取每个多边形的数据
	regionCount = regionList->GetNumRegions();
	//**************************************这段代码需要判断是否为NULL*********************************************

	// 获取边界，boundaries_存储边界点坐标。
	GetBoundariesPtIdx();
	// 获取每个区域的像素，regionPts_存储区域点坐标。
	GetRegionsPtIdx();

	/**************************分割已经完成!*******************************/

	cout << " 影像被分为" << regionCount << "块"<<endl;
}

void TF_TOAto255(float *TOA,byte *DN,int height,int width,int nband) 
{
	for (int j = 0; j < height*width; j++)
	{
		DN[height*width * 0 + j] = (int)((TOA[height*width * 0 + j] - 0.013455) / (0.368287 - 0.013455)*255.0);
		DN[height*width * 1 + j] = (int)((TOA[height*width * 1 + j] - 0.044422) / (0.315785 - 0.044422)*255.0);
		DN[height*width * 2 + j] = (int)((TOA[height*width * 2 + j] - 0.017165) / (0.401000 - 0.017165)*255.0);
		DN[height*width * 3 + j] = (int)((TOA[height*width * 3 + j] - 0.000000) / (0.568269 - 0.000000)*255.0);
		if (DN[height*width * 0 + j] < 0)DN[height*width * 0 + j] = 0;
		if (DN[height*width * 0 + j] > 255)DN[height*width * 0 + j] = 255;
		if (DN[height*width * 1 + j] < 0)DN[height*width * 0 + j] = 0;
		if (DN[height*width * 1 + j] > 255)DN[height*width * 0 + j] = 255;
		if (DN[height*width * 2 + j] < 0)DN[height*width * 0 + j] = 0;
		if (DN[height*width * 2 + j] > 255)DN[height*width * 0 + j] = 255;
		if (DN[height*width * 3 + j] < 0)DN[height*width * 0 + j] = 0;
		if (DN[height*width * 3 + j] > 255)DN[height*width * 0 + j] = 255;
	}
}

int Calcu_Ostu(double **m_gray_Arr, int m_height, int m_width)
{
	int i, j, nThresh = 0, m_gray;

	int nHistogram[256];
	//规一化直方图   
	double fStdHistogram[256];
	double  fGrayAccu[256];
	double  fGrayAve[256];
	double  fAverage = 0;
	double  fTemp, fMax = 0;
	//初始化   
	for (i = 0; i <= 255; i++)
	{
		nHistogram[i] = 0;
		fStdHistogram[i] = 0;
		fGrayAccu[i] = 0;
		fGrayAve[i] = 0;
	}
	//统计直方图   
	// 每行   
	for (i = 0; i < m_height; i++)
	{
		// 每列   
		for (j = 0; j < m_width; j++)
		{
			m_gray = m_gray_Arr[i][j];
			nHistogram[m_gray]++;
		}
	}
	int pnum = m_height*m_width;
	for (i = 0; i <= 255; i++)
	{
		fStdHistogram[i] = nHistogram[i] / (double)pnum;
	}
	for (i = 0; i <= 255; i++)
	{
		for (j = 0; j <= i; j++)
		{
			fGrayAccu[i] += fStdHistogram[j];
			fGrayAve[i] += j*fStdHistogram[j];
		}
		fAverage += i*fStdHistogram[i];
	}
	//计算OSTU   
	for (i = 0; i <= 255; i++)
	{
		fTemp = (fAverage*fGrayAccu[i] - fGrayAve[i])*(fAverage*fGrayAccu[i] - fGrayAve[i]) / (fGrayAccu[i] * (1 - fGrayAccu[i]));
		if (fTemp>fMax)
		{
			fMax = fTemp;
			nThresh = i;
		}
	}
	return nThresh;
}

void savetiffsample(BYTE*SrcImgData1, string sample_dir, int regionNum, Mat boudryMat, BOOL **IsInterestArea, byte *img_src)
{
	int j;
	int count = 0;
	for (int i = 0; i < regionNum; i++)
	{
		isCS[i]=false;
		int num_CS = 0;
		double p_CS;
		int x_max = -1, x_min = width + 1, y_max = -1, y_min = height + 1;
		int x_dis, y_dis;
		string filename,bmpname,tem,file_dir,bmp_dir;
		stringstream tt;
	
		mkdir(sample_dir.c_str());

		int *regionindex = regionList->GetRegionPtIndeces(i);
		int pos, x_t, y_t;
		for (j = 0; j < regionPC[i]; j++)//regionPC[i]每个多边形包含点个数
		{
			pos = *(regionindex + j);
			x_t = pos % width;//列数
			y_t = pos / width;//行数
			if (x_t > x_max)  x_max = x_t;
			if (x_t < x_min)  x_min = x_t;
			if (y_t > y_max)  y_max = y_t;
			if (y_t < y_min)  y_min = y_t;
			if (IsInterestArea[y_t][x_t])
			{
				num_CS++;
			}
		}
		p_CS = (double)num_CS / regionPC[i];
		if (p_CS > thres_p_CS)
		{
			isCS[i] = true;
		}
		if (!isCS[i])
		{
			for (j = 0; j < regionPC[i]; j++)//regionPC[i]每个多边形包含点个数
			{
				pos = *(regionindex + j);
				x_t = pos % width;//列数
				y_t = pos / width;//行数
				boudryMat.ptr(y_t)[x_t * 3] = boudryMat.ptr(y_t)[x_t * 3 + 1] = boudryMat.ptr(y_t)[x_t * 3 + 2] = 0;
			} 
			continue;
		}

		count++;
		tt << count;
		tem = tt.str();
		filename = sample_dir.substr(sample_dir.length() - 1, sample_dir.length()) + "_" + tem + ".tiff";
		bmpname = sample_dir.substr(sample_dir.length() - 1, sample_dir.length()) + "_" + tem + ".bmp";
		file_dir = sample_dir + "\\" + filename;
		bmp_dir = sample_dir + "\\" + bmpname;

		x_dis = x_max - x_min + 1;
		y_dis = y_max - y_min + 1;
		
		BYTE *temp;//保存样本四波段DN值
		Mat tempbmp(y_dis, x_dis, CV_8UC3, Scalar(0));
		temp = new BYTE[x_dis*y_dis * 4];
		for (j = 0; j < x_dis*y_dis*4; j++)
		{
			temp[j] = 0.0;
		}
		for (j = 0; j < regionPC[i]; j++)
		{
			pos = *(regionindex + j);
			x_t = pos % width;//列数
			y_t = pos / width;//行数
			temp[(y_t-y_min)*x_dis + x_t-x_min] = SrcImgData1[y_t*width + x_t];
			temp[x_dis*y_dis + (y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[height*width + y_t*width + x_t];
			temp[x_dis*y_dis * 2 + (y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[height*width * 2 + y_t*width + x_t];
			temp[x_dis*y_dis * 3 + (y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[height*width * 3 + y_t*width + x_t];

			tempbmp.ptr(y_t-y_min)[(x_t-x_min)* 3] = img_src[3 * (y_t*width + x_t)];
			tempbmp.ptr(y_t - y_min)[(x_t - x_min) * 3 + 1] = img_src[3 * (y_t*width + x_t) + 1];
			tempbmp.ptr(y_t - y_min)[(x_t - x_min) * 3 + 2] = img_src[3 * (y_t*width + x_t) + 2];
		}
		imwrite(bmp_dir, tempbmp);
		GDALAllRegister();
		GDALDriver*poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
		if (poDriver == NULL)
		{
			cerr<<"保存样本错误"<<endl;
		} 
		char **papszOptions = NULL;
		GDALDataset *WriteDataSet = poDriver->Create(file_dir.c_str(),
			x_dis, y_dis, 4, GDT_Byte, NULL);
		WriteDataSet->RasterIO(GF_Write, 0, 0, x_dis, y_dis, temp, x_dis, y_dis, GDT_Byte, 4, NULL, 0, 0, 0);
		GDALClose(WriteDataSet);
	}
}

void ct_svm_train()
{
	ifstream file;
	file.open("CloudT.txt", ifstream::in);
	if (!file.is_open())
	{
		return;
	}

	Vector<float> Feature_vec;
	float temp;
	while (file >> temp)
	{
		Feature_vec.push_back(temp);
	}
	file.close();

	size_t size1 = Feature_vec.size() / FeatureNumCT;

	Mat labelMat(size1, 1, CV_32S);
	Mat trainMat(size1, FeatureNumCT, CV_32F);

	auto trans1 = Feature_vec.begin();
	for (int i = 0; i < size1; i++)
	{
		if (i < 110)
		{
			labelMat.at<int>(i, 0) = 1;
		}
		else{
			labelMat.at<int>(i, 0) = 2;
		}
		for (int j = 0; j < FeatureNumCT; j++)
		{
			trainMat.at<float>(i, j) = *(trans1++);
		}
	}

	//设置支持向量机的参数  
	ctparams.svm_type = CvSVM::C_SVC;//SVM类型：使用C支持向量机  
	CvTermCriteria criteria;
	ctparams.kernel_type = CvSVM::RBF;
	ctparams.gamma = 1.0 / FeatureNumCT;
	//params.C=100;
	ctparams.term_crit = cvTermCriteria(CV_TERMCRIT_ITER, 1e+3, 1e-6);

	//训练SVM  
	//建立一个SVM类的实例  
	//训练模型，参数为：输入数据、响应、XX、XX、参数（前面设置过）  
	ctSVM.train(trainMat, labelMat, Mat(), Mat(), ctparams);
}
void cl_svm_train()
{
	ifstream file;
	file.open("CloudL.txt", ifstream::in);
	if (!file.is_open())
	{
		return;
	}

	Vector<float> Feature_vec;
	float temp;
	while (file >> temp)
	{
		Feature_vec.push_back(temp);
	}
	file.close();

	size_t size1 = Feature_vec.size() / FeatureNumCL;

	Mat labelMat(size1, 1, CV_32S);
	Mat trainMat(size1, FeatureNumCL, CV_32F);

	auto trans1 = Feature_vec.begin();
	for (int i = 0; i < size1; i++)
	{
		if (i < 110)
		{
			labelMat.at<int>(i, 0) = 1;
		}
		else{
			labelMat.at<int>(i, 0) = 3;
		}
		for (int j = 0; j < FeatureNumCL; j++)
		{
			trainMat.at<float>(i, j) = *(trans1++);
		}
	}

	//设置支持向量机的参数  
	clparams.svm_type = CvSVM::C_SVC;//SVM类型：使用C支持向量机  
	CvTermCriteria criteria;
	clparams.kernel_type = CvSVM::RBF;
	clparams.gamma = 1.0 / FeatureNumCL;
	//params.C=100;
	clparams.term_crit = cvTermCriteria(CV_TERMCRIT_ITER, 1e+3, 1e-6);

	//训练SVM  
	//建立一个SVM类的实例  
	//训练模型，参数为：输入数据、响应、XX、XX、参数（前面设置过）  
	clSVM.train(trainMat, labelMat, Mat(), Mat(), clparams);
}
void cs_svm_train()
{
	ifstream file;
	file.open("CloudS.txt", ifstream::in);
	if (!file.is_open())
	{
		return;
	}

	Vector<float> Feature_vec;
	float temp;
	while (file >> temp)
	{
		Feature_vec.push_back(temp);
	}
	file.close();

	size_t size1 = Feature_vec.size() / FeatureNumCS;

	Mat labelMat(size1, 1, CV_32S);
	Mat trainMat(size1,FeatureNumCS,CV_32F);

	auto trans1 = Feature_vec.begin();
	for (int i = 0; i < size1;i++)
	{
		if (i < 110)
		{
			labelMat.at<int>(i,0) = 1;
		}
		else{
			labelMat.at<int>(i,0) = 4;
		}
		for (int j = 0; j < FeatureNumCS; j++)
		{
			trainMat.at<float>(i,j) = *(trans1++);
		}
	}

	//设置支持向量机的参数  
	csparams.svm_type = CvSVM::C_SVC;//SVM类型：使用C支持向量机  
	CvTermCriteria criteria;
	csparams.kernel_type = CvSVM::RBF;
	csparams.gamma = 1.0/FeatureNumCS;
	//params.C=100;
	csparams.term_crit = cvTermCriteria(CV_TERMCRIT_ITER, 1e+3, 1e-6);

	//训练SVM  
	//建立一个SVM类的实例  
	//训练模型，参数为：输入数据、响应、XX、XX、参数（前面设置过）  
	csSVM.train(trainMat, labelMat, Mat(), Mat(), csparams);

}
void tl_svm_train()
{
	ifstream file;
	file.open("ThinL.txt", ifstream::in);
	if (!file.is_open())
	{
		return;
	}

	Vector<float> Feature_vec;
	float temp;
	while (file >> temp)
	{
		Feature_vec.push_back(temp);
	}
	file.close();

	size_t size1 = Feature_vec.size() / FeatureNumTL;

	Mat labelMat(size1, 1, CV_32S);
	Mat trainMat(size1, FeatureNumTL, CV_32F);

	auto trans1 = Feature_vec.begin();
	for (int i = 0; i < size1; i++)
	{
		if (i < 100)
		{
			labelMat.at<int>(i, 0) = 2;
		}
		else{
			labelMat.at<int>(i, 0) = 3;
		}
		for (int j = 0; j < FeatureNumTL; j++)
		{
			trainMat.at<float>(i, j) = *(trans1++);
		}
	}

	//设置支持向量机的参数  
	tlparams.svm_type = CvSVM::C_SVC;//SVM类型：使用C支持向量机  
	CvTermCriteria criteria;
	tlparams.kernel_type = CvSVM::RBF;
	tlparams.gamma = 1.0 / FeatureNumTL;
	//params.C=100;
	tlparams.term_crit = cvTermCriteria(CV_TERMCRIT_ITER, 1e+3, 1e-6);

	//训练SVM  
	//建立一个SVM类的实例  
	//训练模型，参数为：输入数据、响应、XX、XX、参数（前面设置过）  
	tlSVM.train(trainMat, labelMat, Mat(), Mat(), tlparams);
}
void ts_svm_train()
{
	ifstream file;
	file.open("ThinS.txt", ifstream::in);
	if (!file.is_open())
	{
		return;
	}

	Vector<float> Feature_vec;
	float temp;
	while (file >> temp)
	{
		Feature_vec.push_back(temp);
	}
	file.close();

	size_t size1 = Feature_vec.size() / FeatureNumTS;

	Mat labelMat(size1, 1, CV_32S);
	Mat trainMat(size1, FeatureNumTS, CV_32F);

	auto trans1 = Feature_vec.begin();
	for (int i = 0; i < size1; i++)
	{
		if (i < 100)
		{
			labelMat.at<int>(i, 0) = 2;
		}
		else{
			labelMat.at<int>(i, 0) = 4;
		}
		for (int j = 0; j < FeatureNumTS; j++)
		{
			trainMat.at<float>(i, j) = *(trans1++);
		}
	}

	//设置支持向量机的参数  
	tsparams.svm_type = CvSVM::C_SVC;//SVM类型：使用C支持向量机  
	CvTermCriteria criteria;
	tsparams.kernel_type = CvSVM::RBF;
	tsparams.gamma = 1.0 / FeatureNumTS;
	//params.C=100;
	tsparams.term_crit = cvTermCriteria(CV_TERMCRIT_ITER, 1e+3, 1e-6);

	//训练SVM  
	//建立一个SVM类的实例  
	//训练模型，参数为：输入数据、响应、XX、XX、参数（前面设置过）  
	tsSVM.train(trainMat, labelMat, Mat(), Mat(), tsparams);
}
void ls_svm_train()
{
	ifstream file;
	file.open("LandS.txt", ifstream::in);
	if (!file.is_open())
	{
		return;
	}

	Vector<float> Feature_vec;
	float temp;
	while (file >> temp)
	{
		Feature_vec.push_back(temp);
	}
	file.close();

	size_t size1 = Feature_vec.size() / FeatureNumLS;

	Mat labelMat(size1, 1, CV_32S);
	Mat trainMat(size1, FeatureNumLS, CV_32F);

	auto trans1 = Feature_vec.begin();
	for (int i = 0; i < size1; i++)
	{
		if (i < 100)
		{
			labelMat.at<int>(i, 0) = 3;
		}
		else{
			labelMat.at<int>(i, 0) = 4;
		}
		for (int j = 0; j < FeatureNumLS; j++)
		{
			trainMat.at<float>(i, j) = *(trans1++);
		}
	}

	//设置支持向量机的参数  
	lsparams.svm_type = CvSVM::C_SVC;//SVM类型：使用C支持向量机  
	CvTermCriteria criteria;
	lsparams.kernel_type = CvSVM::RBF;
	lsparams.gamma = 1.0 / FeatureNumLS;
	//params.C=100;
	lsparams.term_crit = cvTermCriteria(CV_TERMCRIT_ITER, 1e+3, 1e-6);

	//训练SVM  
	//建立一个SVM类的实例  
	//训练模型，参数为：输入数据、响应、XX、XX、参数（前面设置过）  
	lsSVM.train(trainMat, labelMat, Mat(), Mat(), lsparams);
}

void calFeature(BYTE *SrcImgData1, vector<Mat> &VecCT, vector<Mat> &VecCL, vector<Mat> &VecCS, vector<Mat> &VecTL, vector<Mat> &VecTS, vector<Mat> &VecLS, BOOL**IsInterestArea)
{
	/*以下是对每块样本分别提取FeatureNum维向量 */
	isCS = new BOOL[regionList->GetNumRegions()];
	for (int n = 0; n < regionList->GetNumRegions(); n++)
	{
		Mat TempMat(1, 50, CV_32FC1);
		isCS[n] = false;
		int num_CS = 0;
		float p_CS;
		int x_max = -1, x_min = width + 1, y_max = -1, y_min = height + 1;
		int x_dis, y_dis;

		int *regionindex = regionList->GetRegionPtIndeces(n);
		int pos, x_t, y_t;
		for (int j = 0; j < regionPC[n]; j++)//regionPC[i]每个多边形包含点个数
		{
			pos = *(regionindex + j);
			x_t = pos % width;//列数
			y_t = pos / width;//行数
			if (x_t > x_max)  x_max = x_t;
			if (x_t < x_min)  x_min = x_t;
			if (y_t > y_max)  y_max = y_t;
			if (y_t < y_min)  y_min = y_t;
			if (IsInterestArea[y_t][x_t])
			{
				num_CS++;
			}
		}
		p_CS = (double)num_CS / regionPC[n];
		if (p_CS > thres_p_CS)
		{
			isCS[n] = true;
		}
		if (!isCS[n])
		{
			continue;
		}
		x_dis = x_max - x_min + 1;
		y_dis = y_max - y_min + 1;

		float *pSmpdata = new float[y_dis*x_dis * 4];//保存样本四波段DN值
		SmpFeature SmpFeature;
		for (int j = 0; j < x_dis*y_dis * 4; j++)
		{
			pSmpdata[j] = 0.0;
		}
		for (int j = 0; j < regionPC[n]; j++)//regionPC[i]每个多边形包含点个数
		{
			pos = *(regionindex + j);
			x_t = pos % width;//列数
			y_t = pos / width;//行数
			pSmpdata[(y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[y_t*width + x_t];
			pSmpdata[x_dis*y_dis + (y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[height*width + y_t*width + x_t];
			pSmpdata[x_dis*y_dis * 2 + (y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[height*width * 2 + y_t*width + x_t];
			pSmpdata[x_dis*y_dis * 3 + (y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[height*width * 3 + y_t*width + x_t];
		}

		for (int i = 0; i < x_dis*y_dis; i++)
		{
			if (pSmpdata[i] != 0 || pSmpdata[x_dis*y_dis + i] != 0 || pSmpdata[x_dis*y_dis * 2 + i] != 0 || pSmpdata[x_dis*y_dis * 3 + i] != 0)
			{
				SmpFeature.area_pnum++;
				SmpFeature.NIR += pSmpdata[x_dis*y_dis * 3 + i];
				SmpFeature.R += pSmpdata[x_dis*y_dis * 2 + i];
				SmpFeature.G += pSmpdata[x_dis*y_dis * 1 + i];
				SmpFeature.B += pSmpdata[x_dis*y_dis * 0 + i];
			}
		}
		SmpFeature.NIR = SmpFeature.NIR / SmpFeature.area_pnum;
		SmpFeature.R = SmpFeature.R / SmpFeature.area_pnum;
		SmpFeature.G = SmpFeature.G / SmpFeature.area_pnum;
		SmpFeature.B = SmpFeature.B / SmpFeature.area_pnum;

		float v2r, v2g, v2b, v2nir;
		v2r = v2g = v2b = v2nir = 0;
		for (int i = 0; i < x_dis*y_dis; i++)
		{
			if (pSmpdata[i] != 0 || pSmpdata[x_dis*y_dis + i] != 0 || pSmpdata[x_dis*y_dis * 2 + i] != 0 || pSmpdata[x_dis*y_dis * 3 + i] != 0)
			{
				v2nir += (pSmpdata[x_dis*y_dis * 3 + i] - SmpFeature.NIR)* (pSmpdata[x_dis*y_dis * 3 + i] - SmpFeature.NIR);
				v2r += (pSmpdata[x_dis*y_dis * 2 + i] - SmpFeature.R)*(pSmpdata[x_dis*y_dis * 2 + i] - SmpFeature.R);
				v2g += (pSmpdata[x_dis*y_dis * 1 + i] - SmpFeature.G)*(pSmpdata[x_dis*y_dis * 1 + i] - SmpFeature.G);
				v2b += (pSmpdata[x_dis*y_dis * 0 + i] - SmpFeature.B)*(pSmpdata[x_dis*y_dis * 0 + i] - SmpFeature.B);
			}
		}
		TempMat.at<float>(0, 0) = SmpFeature.B;
		TempMat.at<float>(0, 1) = SmpFeature.G;
		TempMat.at<float>(0, 2) = SmpFeature.R;
		TempMat.at<float>(0, 3) = SmpFeature.NIR;
		
		TempMat.at<float>(0, 4) = sqrt(v2r / SmpFeature.area_pnum);
		TempMat.at<float>(0, 5) = sqrt(v2g / SmpFeature.area_pnum);
		TempMat.at<float>(0, 6) = sqrt(v2b / SmpFeature.area_pnum);
		TempMat.at<float>(0, 7) = SmpFeature.v2nir = sqrt(v2nir / SmpFeature.area_pnum);


		SmpFeature.mgray_b = new int*[y_dis];
		SmpFeature.mgray_g = new int*[y_dis];
		SmpFeature.mgray_r = new int*[y_dis];
		SmpFeature.mgray_nr = new int*[y_dis];
		SmpFeature.mgray_h = new double*[y_dis];
		SmpFeature.mgray_i = new double*[y_dis];
		SmpFeature.mgray_s = new double*[y_dis];
		double maxI, minI, maxS, minS;
		int tempCount = 0;
		for (int i = 0; i < y_dis; i++)
		{
			SmpFeature.mgray_b[i] = new int[x_dis];
			SmpFeature.mgray_g[i] = new int[x_dis];
			SmpFeature.mgray_r[i] = new int[x_dis];
			SmpFeature.mgray_nr[i] = new int[x_dis];
			SmpFeature.mgray_h[i] = new double[x_dis];
			SmpFeature.mgray_i[i] = new double[x_dis];
			SmpFeature.mgray_s[i] = new double[x_dis];
			for (int j = 0; j < x_dis; j++)
			{
				SmpFeature.mgray_b[i][j] = pSmpdata[i*x_dis + j];
				SmpFeature.mgray_g[i][j] = pSmpdata[x_dis*y_dis + i*x_dis + j];
				SmpFeature.mgray_r[i][j] = pSmpdata[x_dis*y_dis * 2 + i*x_dis + j];
				SmpFeature.mgray_nr[i][j] = pSmpdata[x_dis*y_dis * 3 + i * x_dis + j];

				transRGBtoHIS(SmpFeature.mgray_r[i][j], SmpFeature.mgray_g[i][j], SmpFeature.mgray_b[i][j],
					SmpFeature.mgray_h[i][j], SmpFeature.mgray_i[i][j], SmpFeature.mgray_s[i][j]);
				if (SmpFeature.mgray_b[i][j] != 0 || SmpFeature.mgray_g[i][j] != 0 || SmpFeature.mgray_r[i][j] != 0)
				{
					if (tempCount == 0)
					{
						maxI = minI = SmpFeature.mgray_i[i][j];
						maxS = minS = SmpFeature.mgray_s[i][j];
						tempCount++;
					}
					else{
						if (maxI < SmpFeature.mgray_i[i][j])maxI = SmpFeature.mgray_i[i][j];
						if (minI > SmpFeature.mgray_i[i][j])minI = SmpFeature.mgray_i[i][j];
						if (maxS < SmpFeature.mgray_s[i][j])maxS = SmpFeature.mgray_s[i][j];
						if (minS > SmpFeature.mgray_s[i][j])minS = SmpFeature.mgray_s[i][j];
					}
				}
			}
		}
		/*cal SF image*/
		tempCount = 0;
		double minSF, maxSF;
		SmpFeature.SF = new double*[y_dis];
		for (int i = 0; i < y_dis; i++)
		{
			SmpFeature.SF[i] = new double[x_dis];
			for (int j = 0; j < x_dis; j++)
			{
				SmpFeature.SF[i][j] = 0;
				if (SmpFeature.mgray_b[i][j] != 0 || SmpFeature.mgray_g[i][j] != 0 || SmpFeature.mgray_r[i][j] != 0)
				{
					SmpFeature.mgray_i[i][j] = (SmpFeature.mgray_i[i][j] - minI) / (maxI - minI);
					SmpFeature.mgray_s[i][j] = (SmpFeature.mgray_s[i][j] - minS) / (maxS - minS);
					SmpFeature.SF[i][j] = (SmpFeature.mgray_i[i][j] + 1.0) / (SmpFeature.mgray_s[i][j] + 1.0);
					if (tempCount == 0)
					{
						maxSF = minSF = SmpFeature.SF[i][j];
						tempCount++;
					}
					else{
						if (maxSF < SmpFeature.SF[i][j])maxSF = SmpFeature.SF[i][j];
						if (minSF > SmpFeature.SF[i][j])minSF = SmpFeature.SF[i][j];
					}
				}
			}
		}
		tempCount = 0;
		double valueST = 0;
		for (int i = 0; i < y_dis; i++)
		{
			for (int j = 0; j < x_dis; j++)
			{
				if (SmpFeature.mgray_b[i][j] != 0 || SmpFeature.mgray_g[i][j] != 0 || SmpFeature.mgray_r[i][j] != 0)
				{
					SmpFeature.SF[i][j] = (SmpFeature.SF[i][j] - minSF) / (maxSF - minSF)*255.0;
					valueST += SmpFeature.mgray_nr[i][j] - (1.73*SmpFeature.SF[i][j] - 219.2);
					tempCount++;
				}
			}
		}
		valueST /= tempCount;
		TempMat.at<float>(0, 8) = valueST;

		Cal_IWCS_LTP(SmpFeature.mgray_r, y_dis, x_dis, SmpFeature.LTP_features);
		for (int i = 0; i < 18; i++){ TempMat.at<float>(0, 9 + i) = SmpFeature.LTP_features[i]; }

		calColorGLCM(pSmpdata, y_dis, x_dis, SmpFeature.CGLCM_features);
		for (int i = 0; i < 9; i++){ TempMat.at<float>(0, 27 + i) = SmpFeature.CGLCM_features[i]; }

		float CurHistogram[10];
		calEdgeCurHis(SmpFeature.mgray_b, y_dis, x_dis, CurHistogram);
		TempMat.at<float>(0, 36) = CurHistogram[0] + CurHistogram[1] + CurHistogram[2];

		double aValueMCLFD[32];
		calMCLFD(SmpFeature.mgray_b, y_dis, x_dis, aValueMCLFD);
		for (int i = 0; i < 8; i++){
			TempMat.at<float>(0, 37 + i) = aValueMCLFD[4 * i];
		}

		Mat MatCT(1, FeatureNumCT, CV_32FC1);
		Mat MatCL(1, FeatureNumCL, CV_32FC1);
		Mat MatCS(1, FeatureNumCS, CV_32FC1);
		Mat MatTL(1, FeatureNumTL, CV_32FC1);
		Mat MatTS(1, FeatureNumTS, CV_32FC1);
		Mat MatLS(1, FeatureNumLS, CV_32FC1);
	    
		for (int i = 0; i < 4; i++)
		{
			MatCT.at<float>(0, i) = TempMat.at<float>(0, i);
			MatCL.at<float>(0, i) = TempMat.at<float>(0, i);
			MatTS.at<float>(0, i) = TempMat.at<float>(0, i);
			MatLS.at<float>(0, i) = TempMat.at<float>(0, i);
		}
		for (int i = 0; i < 4; i++)
		{
			MatCL.at<float>(0, 4 + i) = TempMat.at<float>(0, 4 + i);
			MatTL.at<float>(0, i) = TempMat.at<float>(0, 4 + i);
		}
		MatLS.at<float>(0, 4) = TempMat.at<float>(0, 5);

		MatCS.at<float>(0, 0) = TempMat.at<float>(0, 7);
		MatTS.at<float>(0, 4) = TempMat.at<float>(0, 7);

		MatCT.at<float>(0, 4) = TempMat.at<float>(0, 8);
		MatTL.at<float>(0, 4) = TempMat.at<float>(0, 8);


		for (int i = 0; i < 18; i++)
		{
			MatTS.at<float>(0, 5 + i) = TempMat.at<float>(0, 9 + i);
		}

		for (int i = 0; i < 9; i++)
		{
			MatCT.at<float>(0, 5 + i) = TempMat.at<float>(0, 27 + i);
			MatCS.at<float>(0, 1 + i) = TempMat.at<float>(0, 27 + i);
			MatTL.at<float>(0, 5 + i) = TempMat.at<float>(0, 27 + i);
			MatTS.at<float>(0, 23 + i) = TempMat.at<float>(0, 27 + i);
			MatLS.at<float>(0, 5 + i) = TempMat.at<float>(0, 27 + i);
		}

		MatCL.at<float>(0, 8) = TempMat.at<float>(0, 36);
		MatCS.at<float>(0, 10) = TempMat.at<float>(0, 36);
		MatTL.at<float>(0, 14) = TempMat.at<float>(0, 36);
	//	MatTS.at<float>(0, 32) = TempMat.at<float>(0, 36);

		for (int i = 0; i < 8; i++)
		{
			MatCL.at<float>(0, 9 + i) = TempMat.at<float>(0, 37 + i);
			MatCS.at<float>(0, 11 + i) = TempMat.at<float>(0, 37 + i);
			MatTL.at<float>(0, 15 + i) = TempMat.at<float>(0, 37 + i);
	//		MatTS.at<float>(0, 33 + i) = TempMat.at<float>(0, 37 + i);
		}

		VecCT.push_back(MatCT);
		VecCL.push_back(MatCL);
		VecCS.push_back(MatCS);
		VecTL.push_back(MatTL);
		VecTS.push_back(MatTS);
		VecLS.push_back(MatLS);

		delete[] pSmpdata;
		for (int i = 0; i < y_dis; i++)
		{
			delete[]SmpFeature.mgray_b[i];delete[]SmpFeature.mgray_g[i];
			delete[]SmpFeature.mgray_r[i];delete[]SmpFeature.mgray_nr[i];
			delete[]SmpFeature.mgray_h[i];delete[]SmpFeature.mgray_i[i];
			delete[]SmpFeature.mgray_s[i]; delete[]SmpFeature.SF[i];
		}
		delete[]SmpFeature.mgray_b; delete[]SmpFeature.mgray_g;
		delete[]SmpFeature.mgray_r; delete[]SmpFeature.mgray_nr;
		delete[]SmpFeature.mgray_h; delete[]SmpFeature.mgray_i;
		delete[]SmpFeature.mgray_s; delete[]SmpFeature.SF;
	}

	int length = VecCT.size();

	/*特征标准化*/
	ifstream file;
	float ParamsCT[FeatureNumCT][2], ParamsCL[FeatureNumCL][2], ParamsCS[FeatureNumCS][2], ParamsTL[FeatureNumTL][2], ParamsTS[FeatureNumTS][2], ParamsLS[FeatureNumLS][2];
	file.open("ParamsCT.txt");
	for (int i = 0; i < FeatureNumCT; i++)
	{
		file >> ParamsCT[i][0];
		file >> ParamsCT[i][1];
	}
	file.close();

	file.open("ParamsCL.txt");
	for (int i = 0; i < FeatureNumCL; i++)
	{
		file >> ParamsCL[i][0];
		file >> ParamsCL[i][1];
	}
	file.close();

	file.open("ParamsCS.txt");
	for (int i = 0; i < FeatureNumCS; i++)
	{
		file >> ParamsCS[i][0];
		file >> ParamsCS[i][1];
	}
	file.close();

	file.open("ParamsTL.txt");
	for (int i = 0; i < FeatureNumTL; i++)
	{
		file >> ParamsTL[i][0];
		file >> ParamsTL[i][1];
	}
	file.close();

	file.open("ParamsTS.txt");
	for (int i = 0; i < FeatureNumTS; i++)
	{
		file >> ParamsTS[i][0];
		file >> ParamsTS[i][1];
	}
	file.close();

	file.open("ParamsLS.txt");
	for (int i = 0; i < FeatureNumLS; i++)
	{
		file >> ParamsLS[i][0];
		file >> ParamsLS[i][1];
	}
	file.close();

	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < FeatureNumCT; j++)
		{
			VecCT[i].at<float>(0, j) = (VecCT[i].at<float>(0, j) - ParamsCT[j][0]) / ParamsCT[j][1];
		}
		for (int j = 0; j < FeatureNumCL; j++)
		{
			VecCL[i].at<float>(0, j) = (VecCL[i].at<float>(0, j) - ParamsCL[j][0]) / ParamsCL[j][1];
		}
		for (int j = 0; j < FeatureNumCS; j++)
		{
			VecCS[i].at<float>(0, j) = (VecCS[i].at<float>(0, j) - ParamsCS[j][0]) / ParamsCS[j][1];
		}
		for (int j = 0; j < FeatureNumTL; j++)
		{
			VecTL[i].at<float>(0, j) = (VecTL[i].at<float>(0, j) - ParamsTL[j][0]) / ParamsTL[j][1];
		}
		for (int j = 0; j < FeatureNumTS; j++)
		{
			VecTS[i].at<float>(0, j) = (VecTS[i].at<float>(0, j) - ParamsTS[j][0]) / ParamsTS[j][1];
		}
		for (int j = 0; j < FeatureNumLS; j++)
		{
			VecLS[i].at<float>(0, j) = (VecLS[i].at<float>(0, j) - ParamsLS[j][0]) / ParamsLS[j][1];
		}
	}

	
}

//void calFeatureCS(BYTE *SrcImgData1,Vector<Mat> &FeatureVec)
//{
//	/*以下是对每块样本分别提取FeatureNum维向量 */
//	isCS = new BOOL[regionList->GetNumRegions()];
//	for (int n = 0; n < regionList->GetNumRegions(); n++)
//	{
//		Mat TempMat(1, 100, CV_32FC1);
//		isCS[n] = false;
//		int num_CS = 0;
//		float p_CS;
//		int x_max = -1, x_min = width + 1, y_max = -1, y_min = height + 1;
//		int x_dis, y_dis;
//	
//		int *regionindex = regionList->GetRegionPtIndeces(n);
//		int pos, x_t, y_t;
//		for (int j = 0; j < regionPC[n]; j++)//regionPC[i]每个多边形包含点个数
//		{
//			pos = *(regionindex + j);
//			x_t = pos % width;//列数
//			y_t = pos / width;//行数
//			if (x_t > x_max)  x_max = x_t;
//			if (x_t < x_min)  x_min = x_t;
//			if (y_t > y_max)  y_max = y_t;
//			if (y_t < y_min)  y_min = y_t;
//			if (IsInterestArea[y_t][x_t])
//			{
//				num_CS++;
//			}
//		}
//		p_CS = (double)num_CS / regionPC[n];
//		if (p_CS > 0.3)
//		{
//			isCS[n] = true;
//		}
//		if (!isCS[n])
//		{
//			continue;
//		}
//		x_dis = x_max - x_min + 1;
//		y_dis = y_max - y_min + 1;
//
//		float *pSmpdata = new float[y_dis*x_dis * 4];//保存样本四波段DN值
//		SmpFeature SmpFeature;
//		for (int j = 0; j < x_dis*y_dis * 4; j++)
//		{
//			pSmpdata[j] = 0.0;
//		}
//		for (int j = 0; j < regionPC[n]; j++)//regionPC[i]每个多边形包含点个数
//		{
//			pos = *(regionindex + j);
//			x_t = pos % width;//列数
//			y_t = pos / width;//行数
//			pSmpdata[(y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[y_t*width + x_t];
//			pSmpdata[x_dis*y_dis + (y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[height*width + y_t*width + x_t];
//			pSmpdata[x_dis*y_dis * 2 + (y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[height*width * 2 + y_t*width + x_t];
//			pSmpdata[x_dis*y_dis * 3 + (y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[height*width * 3 + y_t*width + x_t];
//		}
//
//		SmpFeature.mgray_b = new int*[y_dis];
//		SmpFeature.mgray_g = new int*[y_dis];
//		SmpFeature.mgray_r = new int*[y_dis];
//		SmpFeature.mgray_nr = new int*[y_dis];
//		for (int i = 0; i < y_dis; i++)
//		{
//			SmpFeature.mgray_b[i] = new int[x_dis];
//			SmpFeature.mgray_g[i] = new int[x_dis];
//			SmpFeature.mgray_r[i] = new int[x_dis];
//			SmpFeature.mgray_nr[i] = new int[x_dis];
//			for (int j = 0; j < x_dis; j++)
//			{
//				SmpFeature.mgray_b[i][j] = pSmpdata[i*x_dis + j];
//				SmpFeature.mgray_g[i][j] = pSmpdata[x_dis*y_dis + i*x_dis + j];
//				SmpFeature.mgray_r[i][j] = pSmpdata[x_dis*y_dis * 2 + i*x_dis + j];
//				SmpFeature.mgray_nr[i][j] = pSmpdata[x_dis*y_dis * 3 + i * x_dis + j];
//			}
//		}
//		for (int i = 0; i < x_dis*y_dis; i++)
//		{
//			if (pSmpdata[i] != 0 || pSmpdata[x_dis*y_dis + i] != 0 || pSmpdata[x_dis*y_dis * 2 + i] != 0 || pSmpdata[x_dis*y_dis * 3 + i] != 0)
//			{
//				SmpFeature.area_pnum++;
//				SmpFeature.NIR += pSmpdata[x_dis*y_dis * 3 + i];
//				SmpFeature.R += pSmpdata[x_dis*y_dis * 2 + i];
//				SmpFeature.G += pSmpdata[x_dis*y_dis * 1 + i];
//				SmpFeature.B += pSmpdata[x_dis*y_dis * 0 + i];
//			}
//		}
//		SmpFeature.NIR = SmpFeature.NIR / SmpFeature.area_pnum;
//		SmpFeature.R = SmpFeature.R / SmpFeature.area_pnum;
//		SmpFeature.G = SmpFeature.G / SmpFeature.area_pnum;
//		SmpFeature.B = SmpFeature.B / SmpFeature.area_pnum;
//
//		float v2r, v2g, v2b, v2nir;
//		v2r = v2g = v2b = v2nir = 0;
//
//		for (int i = 0; i < x_dis*y_dis; i++)
//		{
//			if (pSmpdata[i] != 0 || pSmpdata[x_dis*y_dis + i] != 0 || pSmpdata[x_dis*y_dis * 2 + i] != 0 || pSmpdata[x_dis*y_dis * 3 + i] != 0)
//			{
//				v2nir += (pSmpdata[x_dis*y_dis * 3 + i] - SmpFeature.NIR)* (pSmpdata[x_dis*y_dis * 3 + i] - SmpFeature.NIR);
//				v2r += (pSmpdata[x_dis*y_dis * 2 + i] - SmpFeature.R)*(pSmpdata[x_dis*y_dis * 2 + i] - SmpFeature.R);
//				v2g += (pSmpdata[x_dis*y_dis * 1 + i] - SmpFeature.G)*(pSmpdata[x_dis*y_dis * 1 + i] - SmpFeature.G);
//				v2b += (pSmpdata[x_dis*y_dis * 0 + i] - SmpFeature.B)*(pSmpdata[x_dis*y_dis * 0 + i] - SmpFeature.B);
//			}
//		}
//		TempMat.at<float>(0, 0) = SmpFeature.B;
//		TempMat.at<float>(0, 1) = SmpFeature.G;
//		TempMat.at<float>(0, 2) = SmpFeature.R;
//		TempMat.at<float>(0, 3) = SmpFeature.NIR;
//		TempMat.at<float>(0, 3)=SmpFeature.v2nir = sqrt(v2nir / SmpFeature.area_pnum);
//	/*	TempMat.at<float>(0, 2) = sqrt(v2r / SmpFeature.area_pnum);
//		TempMat.at<float>(0, 3) = sqrt(v2g / SmpFeature.area_pnum);
//		TempMat.at<float>(0, 4) = sqrt(v2b / SmpFeature.area_pnum);
//	*/
//	    double b_FractalD, g_FractalD, r_FractalD, nr_FractalD;
//	    calFractalDimension(SmpFeature.mgray_b, y_dis, x_dis, b_FractalD);
//	    calFractalDimension(SmpFeature.mgray_g, y_dis, x_dis, g_FractalD);
//	    calFractalDimension(SmpFeature.mgray_r, y_dis, x_dis, r_FractalD);
//	    calFractalDimension(SmpFeature.mgray_nr, y_dis, x_dis, nr_FractalD);
//		TempMat.at<float>(0, 0) = SmpFeature.fractalDimension = (b_FractalD + g_FractalD + r_FractalD + nr_FractalD) / 4;
//
//	
//	    calColorGLCM(pSmpdata, y_dis, x_dis, SmpFeature.CGLCM_features);
//
//		for (int i = 0; i < 9; i++){ TempMat.at<float>(0,4 + i) = SmpFeature.CGLCM_features[i]; }
//		
//		float CurHistogram[10];
//		calEdgeCurHis(SmpFeature.mgray_b, y_dis, x_dis,CurHistogram);
//		TempMat.at<float>(0, 13) = CurHistogram[0] + CurHistogram[1] + CurHistogram[2];
//		for (int i = 0; i < 10; i++){
//			TempMat.at<float>(0, 14 + i) = CurHistogram[i];
//		}
//
//		double aValueMCLFD[32];
//		calMCLFD(SmpFeature.mgray_b, y_dis, x_dis, aValueMCLFD);
//		for (int i = 0; i < 8; i++){
//			TempMat.at<float>(0, 14 + i) = aValueMCLFD[4 * i];
//		}
//
//		Cal_IWCS_LTP(SmpFeature.mgray_r, y_dis, x_dis, SmpFeature.LTP_features);
//		for (int i = 0; i < 18; i++){ TempMat.at<float>(0, 5 + i) = SmpFeature.LTP_features[i]; }
//
//
//		Mat TrainMat(1,FeatureNumCS , CV_32FC1);
//		for (int i = 0; i < FeatureNumCS; i++)
//		{
//			TrainMat.at<float>(0, i) = TempMat.at<float>(0, i);
//		}
//
//		FeatureVec.push_back(TrainMat);
//     }
//
//	 int length = FeatureVec.size();
//	
//	 /*特征标准化*/
//	 ifstream file;
//	 file.open("NormalParamsCS.txt");
//	 for (int i = 0; i < FeatureNumCS; i++)
//	 {
//		 file >> Normal[i][0];
//		 file >> Normal[i][1];
//	 }
//	 file.close();
//
//	 for (int i = 0; i < FeatureNumCS; i++)
//	 {
//		 for (int j = 0; j < length; j++)
//		 {
//			 FeatureVec[j].at<float>(0, i) = (FeatureVec[j].at<float>(0, i) - Normal[i][0]) / Normal[i][1];
//		 }
//	 }
//}
//
//void calFeatureCL(BYTE *SrcImgData1, Vector<Mat> &FeatureVec)
//{
//	/*以下是对每块样本分别提取FeatureNum维向量 */
//	isCS = new BOOL[regionList->GetNumRegions()];
//	for (int n = 0; n < regionList->GetNumRegions(); n++)
//	{
//		Mat TempMat(1, 100, CV_32FC1);
//		isCS[n] = false;
//		int num_CS = 0;
//		float p_CS;
//		int x_max = -1, x_min = width + 1, y_max = -1, y_min = height + 1;
//		int x_dis, y_dis;
//
//		int *regionindex = regionList->GetRegionPtIndeces(n);
//		int pos, x_t, y_t;
//		for (int j = 0; j < regionPC[n]; j++)//regionPC[i]每个多边形包含点个数
//		{
//			pos = *(regionindex + j);
//			x_t = pos % width;//列数
//			y_t = pos / width;//行数
//			if (x_t > x_max)  x_max = x_t;
//			if (x_t < x_min)  x_min = x_t;
//			if (y_t > y_max)  y_max = y_t;
//			if (y_t < y_min)  y_min = y_t;
//			if (IsInterestArea[y_t][x_t])
//			{
//				num_CS++;
//			}
//		}
//		p_CS = (double)num_CS / regionPC[n];
//		if (p_CS > 0.3)
//		{
//			isCS[n] = true;
//		}
//		if (!isCS[n])
//		{
//			continue;
//		}
//		x_dis = x_max - x_min + 1;
//		y_dis = y_max - y_min + 1;
//
//		float *pSmpdata = new float[y_dis*x_dis * 4];//保存样本四波段DN值
//		SmpFeature SmpFeature;
//		for (int j = 0; j < x_dis*y_dis * 4; j++)
//		{
//			pSmpdata[j] = 0.0;
//		}
//		for (int j = 0; j < regionPC[n]; j++)//regionPC[i]每个多边形包含点个数
//		{
//			pos = *(regionindex + j);
//			x_t = pos % width;//列数
//			y_t = pos / width;//行数
//			pSmpdata[(y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[y_t*width + x_t];
//			pSmpdata[x_dis*y_dis + (y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[height*width + y_t*width + x_t];
//			pSmpdata[x_dis*y_dis * 2 + (y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[height*width * 2 + y_t*width + x_t];
//			pSmpdata[x_dis*y_dis * 3 + (y_t - y_min)*x_dis + x_t - x_min] = SrcImgData1[height*width * 3 + y_t*width + x_t];
//		}
//
//		SmpFeature.mgray_b = new int*[y_dis];
//		SmpFeature.mgray_g = new int*[y_dis];
//		SmpFeature.mgray_r = new int*[y_dis];
//		SmpFeature.mgray_nr = new int*[y_dis];
//		for (int i = 0; i < y_dis; i++)
//		{
//			SmpFeature.mgray_b[i] = new int[x_dis];
//			SmpFeature.mgray_g[i] = new int[x_dis];
//			SmpFeature.mgray_r[i] = new int[x_dis];
//			SmpFeature.mgray_nr[i] = new int[x_dis];
//			for (int j = 0; j < x_dis; j++)
//			{
//				SmpFeature.mgray_b[i][j] = pSmpdata[i*x_dis + j];
//				SmpFeature.mgray_g[i][j] = pSmpdata[x_dis*y_dis + i*x_dis + j];
//				SmpFeature.mgray_r[i][j] = pSmpdata[x_dis*y_dis * 2 + i*x_dis + j];
//				SmpFeature.mgray_nr[i][j] = pSmpdata[x_dis*y_dis * 3 + i * x_dis + j];
//			}
//		}
//		for (int i = 0; i < x_dis*y_dis; i++)
//		{
//			if (pSmpdata[i] != 0 || pSmpdata[x_dis*y_dis + i] != 0 || pSmpdata[x_dis*y_dis * 2 + i] != 0 || pSmpdata[x_dis*y_dis * 3 + i] != 0)
//			{
//				SmpFeature.area_pnum++;
//				SmpFeature.NIR += pSmpdata[x_dis*y_dis * 3 + i];
//				SmpFeature.R += pSmpdata[x_dis*y_dis * 2 + i];
//				SmpFeature.G += pSmpdata[x_dis*y_dis * 1 + i];
//				SmpFeature.B += pSmpdata[x_dis*y_dis * 0 + i];
//			}
//		}
//		SmpFeature.NIR = SmpFeature.NIR / SmpFeature.area_pnum;
//		SmpFeature.R = SmpFeature.R / SmpFeature.area_pnum;
//		SmpFeature.G = SmpFeature.G / SmpFeature.area_pnum;
//		SmpFeature.B = SmpFeature.B / SmpFeature.area_pnum;
//
//	/*	float v2r, v2g, v2b, v2nir;
//		v2r = v2g = v2b = v2nir = 0;
//
//		for (int i = 0; i < x_dis*y_dis; i++)
//		{
//			if (pSmpdata[i] != 0 || pSmpdata[x_dis*y_dis + i] != 0 || pSmpdata[x_dis*y_dis * 2 + i] != 0 || pSmpdata[x_dis*y_dis * 3 + i] != 0)
//			{
//				v2nir += (pSmpdata[x_dis*y_dis * 3 + i] - SmpFeature.NIR)* (pSmpdata[x_dis*y_dis * 3 + i] - SmpFeature.NIR);
//				v2r += (pSmpdata[x_dis*y_dis * 2 + i] - SmpFeature.R)*(pSmpdata[x_dis*y_dis * 2 + i] - SmpFeature.R);
//				v2g += (pSmpdata[x_dis*y_dis * 1 + i] - SmpFeature.G)*(pSmpdata[x_dis*y_dis * 1 + i] - SmpFeature.G);
//				v2b += (pSmpdata[x_dis*y_dis * 0 + i] - SmpFeature.B)*(pSmpdata[x_dis*y_dis * 0 + i] - SmpFeature.B);
//			}
//		}*/
//		TempMat.at<float>(0, 0) = SmpFeature.B;
//		TempMat.at<float>(0, 1) = SmpFeature.G;
//		TempMat.at<float>(0, 2) = SmpFeature.R;
//		TempMat.at<float>(0, 3) = SmpFeature.NIR;
//	    TempMat.at<float>(0, 4) = SmpFeature.v2nir = sqrt(v2nir / SmpFeature.area_pnum);
//		/*	TempMat.at<float>(0, 2) = sqrt(v2r / SmpFeature.area_pnum);
//		TempMat.at<float>(0, 3) = sqrt(v2g / SmpFeature.area_pnum);
//		TempMat.at<float>(0, 4) = sqrt(v2b / SmpFeature.area_pnum);
//		*/
//		   double b_FractalD, g_FractalD, r_FractalD, nr_FractalD;
//		   calFractalDimension(SmpFeature.mgray_b, y_dis, x_dis, b_FractalD);
//		   calFractalDimension(SmpFeature.mgray_g, y_dis, x_dis, g_FractalD);
//		   calFractalDimension(SmpFeature.mgray_r, y_dis, x_dis, r_FractalD);
//		   calFractalDimension(SmpFeature.mgray_nr, y_dis, x_dis, nr_FractalD);
//		TempMat.at<float>(0, 0) = SmpFeature.fractalDimension = (b_FractalD + g_FractalD + r_FractalD + nr_FractalD) / 4;
//
//
//		calColorGLCM(pSmpdata, y_dis, x_dis, SmpFeature.CGLCM_features);
//
//		for (int i = 0; i < 9; i++){ TempMat.at<float>(0, 5 + i) = SmpFeature.CGLCM_features[i]; }
//
//		float CurHistogram[10];
//		calEdgeCurHis(SmpFeature.mgray_b, y_dis, x_dis, CurHistogram);
//		TempMat.at<float>(0, 0) = CurHistogram[0] + CurHistogram[1] + CurHistogram[2];
//		for (int i = 0; i < 10; i++){
//			TempMat.at<float>(0, 14 + i) = CurHistogram[i];
//		}
//
//		double aValueMCLFD[32];
//		calMCLFD(SmpFeature.mgray_b, y_dis, x_dis, aValueMCLFD);
//		for (int i = 0; i < 8; i++){
//			TempMat.at<float>(0, 1 + i) = aValueMCLFD[4 * i];
//		}
//
//		Cal_IWCS_LTP(SmpFeature.mgray_r, y_dis, x_dis, SmpFeature.LTP_features);
//		for (int i = 0; i < 18; i++){ TempMat.at<float>(0, 5 + i) = SmpFeature.LTP_features[i]; }
//
//
//		Mat TrainMat(1, FeatureNumCL, CV_32FC1);
//		for (int i = 0; i < FeatureNumCL; i++)
//		{
//			TrainMat.at<float>(0, i) = TempMat.at<float>(0, i);
//		}
//
//		FeatureVec.push_back(TrainMat);
//	}
//
//	int length = FeatureVec.size();
//
//	/*特征标准化*/
//	ifstream file;
//	file.open("NormalParamsCL.txt");
//	for (int i = 0; i < FeatureNumCL; i++)
//	{
//		file >> Normal[i][0];
//		file >> Normal[i][1];
//	}
//	file.close();
//
//	for (int i = 0; i < FeatureNumCL; i++)
//	{
//		for (int j = 0; j < length; j++)
//		{
//			FeatureVec[j].at<float>(0, i) = (FeatureVec[j].at<float>(0, i) - Normal[i][0]) / Normal[i][1];
//		}
//	}
//}


const string Src_dir("F:\\Project\\TotalT");
string PrecisionTxtDir = Src_dir + "\\Precision.txt";
int countPrecision = 0, countPrecisionCloud = 0, countPrecisionSnow = 0;
double overall_cloud_accuracy = 0.0, overall_cloud_precision = 0.0, overall_cloud_recall = 0.0;
double overall_snow_accuracy = 0.0, overall_snow_precision = 0.0, overall_snow_recall = 0.0;
mutex mu;

void CloudSnow_Dtec(vector<string> &All_imgfile,int threadlabel)
{
	vector<string>::iterator beg, endd;
	if (threadlabel == 0)
	{
		beg = All_imgfile.begin();
		endd = All_imgfile.begin() + (All_imgfile.end() - All_imgfile.begin())/2;
	}
	else if (threadlabel == 1)
	{
		beg = All_imgfile.begin() + (All_imgfile.end() - All_imgfile.begin()) / 2;
		endd = All_imgfile.end();
	}
	for (vector<string>::iterator itr = beg; itr != endd; itr++)
	{
	

		string Src_img;
		Src_img = Src_dir + "\\" + *itr;

		

		if (threadlabel == 0)
			cout << "thread01" << endl;
		else if (threadlabel == 1)
			cout << "thread02" << endl;

		/*load data*/
		int nband;
		GDALAllRegister();
		GDALDataset *poDataset = NULL;
		poDataset = (GDALDataset*)GDALOpen(Src_img.c_str(), GA_ReadOnly);
		if (poDataset == NULL)
		{
			GDALClose(poDataset);
			return;
		}

		mu.lock();
		width = poDataset->GetRasterXSize();
		height = poDataset->GetRasterYSize();
		nband = poDataset->GetRasterCount();
		mu.unlock();

		cout << *itr << " processing:" << endl;
		cout << "Drivers:" << poDataset->GetDriver()->GetDescription() << poDataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME) << endl;
		cout << "Size is " << poDataset->GetRasterXSize() << "," << poDataset->GetRasterYSize() << "," << poDataset->GetRasterCount() << endl;

		float *pImgdata;
		pImgdata = NULL;
		pImgdata = new float[width*height*nband];
		poDataset->RasterIO(
			GF_Read, 0, 0, width, height, pImgdata, width, height, GDT_Float32, nband, NULL, 0, 0, 0);

		/*              image segmenting             */
		byte *pImgdata1 = new byte[width*height*nband];
		TF_TOAto255(pImgdata, pImgdata1, height, width, nband);

		/*SF图像*/
		double **SF, **HOT;
		SF = new double*[height];
		HOT = new double*[height];
		for (int i = 0; i < height; i++)
		{
			SF[i] = new double[width];
			HOT[i] = new double[width];
			for (int j = 0; j < width; j++)
			{
				SF[i][j] = 0;
				HOT[i][j] = pImgdata[i*width + j] - 0.5*pImgdata[height*width * 2 + i*width + j];
			}
		}
		calSF(pImgdata1, height, width, SF);
		int ostuSF = Calcu_Ostu(SF, height, width);


		/*生成HOT,SF图像*/
		//string HOT_des, seg_des;
		//size_t len = Src_img.length();
		//HOT_des = Src_img.substr(0, len - 5) + "_HOT.bmp";
		//seg_des = Src_img.substr(0, len - 5) + "_seg.bmp";
		//Mat SFimg(height, width, CV_8UC1, Scalar(0));
		//Mat HOTimg(height, width, CV_8UC1, Scalar(0));
		//Mat AREAimg(height, width, CV_8UC1, Scalar(0));
		//for (int i = 0; i < height; i++)
		//{
		//	for (int j = 0; j < width; j++)
		//	{
		//		if (HOT[i][j]>0.08)
		//		{
		//			HOTimg.at<uchar>(i, j) = 255;
		//		}
		//	}
		//}
		//for (int i = 0; i < height; i++)
		//{
		//	for (int j = 0; j < width; j++)
		//	{
		//		if (SF[i][j]>ostuSF)
		//		{
		//			SFimg.at<uchar>(i, j) = 255;
		//		}
		//	}
		//}

		BOOL **IsInterestArea;
		IsInterestArea = new BOOL*[height];
		for (int i = 0; i < height; i++)
		{
			IsInterestArea[i] = new BOOL[width];
			for (int j = 0; j < width; j++)
			{
				if (HOT[i][j]>0.08 && SF[i][j]>ostuSF)
				{
					IsInterestArea[i][j] = true;
				}
				else{
					IsInterestArea[i][j] = false;
				}
			}
		}
		//for (int i = 0; i < height; i++)
		//{
		//	for (int j = 0; j < width; j++)
		//	{
		//		if (IsInterestArea[i][j])
		//		{
		//			AREAimg.at<uchar>(i, j) = 255;
		//		}
		//	}
		//}
		//string SF_des,IntrestArea;
		//SF_des = Src_img.substr(0, len - 5) + "_SF.bmp";
		//IntrestArea = Src_img.substr(0, len - 5) + "_InArea.bmp";
		//imwrite(HOT_des, HOTimg);
		//imwrite(SF_des, SFimg);
		//imwrite(IntrestArea, AREAimg);
		byte *img_src;
		img_src = new byte[width*height * 4];
		Mat img11(height, width, CV_8UC3, Scalar(0));

		for (int i = 0; i < height*width; i++)
		{
			/*432*/
			img_src[3 * i] = pImgdata1[height*width * 1 + i];
			img_src[3 * i + 1] = pImgdata1[height*width * 2 + i];
			img_src[3 * i + 2] = pImgdata1[height*width * 3 + i];
		}
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				img11.ptr(i)[j * 3] = img_src[3 * (i*width + j)];
				img11.ptr(i)[j * 3 + 1] = img_src[3 * (i*width + j) + 1];
				img11.ptr(i)[j * 3 + 2] = img_src[3 * (i*width + j) + 2];
			}
		}

		if (iProc)
		{
			delete iProc;
			iProc = NULL;
		}
		iProc = new msImageProcessor;

		if (img_src != NULL)
		{
			mu.lock();
			ReadImgData(img_src);
			operation = 3;
			imageseg();
		}
		else
		{
			cerr << "请先打开影像！" << endl;
		}

		/***                     样本保存                     ***/
		//isCS = new BOOL[regionList->GetNumRegions()];
		//savetiffsample(pImgdata1,Src_img.substr(0, len - 5), regionList->GetNumRegions(),img11,IsInterestArea,img_src);
		//cout << "样本保存完成。" <<'\n'<<'\n'<< endl;
		//for (size_t i = 0; i < boundaries_->n_; i++)
		//{
		//	int boundaryLable;
		//	boundaryLable = regionLables[boundaries_->y_[i] * width + boundaries_->x_[i]];
		//	if (isCS[boundaryLable])
		//	{
		//		img11.ptr(boundaries_->y_[i])[boundaries_->x_[i] * 3] = 0;
		//		img11.ptr(boundaries_->y_[i])[boundaries_->x_[i] * 3 + 1] = 0;
		//		img11.ptr(boundaries_->y_[i])[boundaries_->x_[i] * 3 + 2] = 255;
		//	}
		//}
		//imwrite(seg_des, img11);
		//for (size_t i = 0; i < boundaries_->n_; i++)
		//{
		//		img11.ptr(boundaries_->y_[i])[boundaries_->x_[i] * 3] = 0;
		//		img11.ptr(boundaries_->y_[i])[boundaries_->x_[i] * 3 + 1] = 0;
		//		img11.ptr(boundaries_->y_[i])[boundaries_->x_[i] * 3 + 2] = 255;
		//}
		//imwrite(seg_des, img11);
		//continue;

		/***                      Feature calculating                     ***/
		vector<Mat> VecCT, VecCL, VecCS, VecTL, VecTS, VecLS;
		calFeature(pImgdata1, VecCT, VecCL, VecCS, VecTL, VecTS, VecLS, IsInterestArea);

		/*	ofstream file;
		file.open("regionFeature.txt");
		for (auto iter = FeatureVec.begin(); iter != FeatureVec.end(); iter++)
		{
		for (int i = 0; i < 18; i++)
		{
		file << (*iter).at<float>(0,i) << " ";
		}
		file << endl;
		}
		file.close();*/

		/***                      SVM predict                     ***/
		string regionP_dir(Src_img.substr(0, Src_img.length() - 5));
		regionP_dir += "_regionP.txt";
		ofstream file;
		file.open(regionP_dir);
		int length;
		length = VecCT.size();
		float *VClassCT = new float[length];		float *VScoreCT = new float[length];
		float *VClassCL = new float[length];		float *VScoreCL = new float[length];
		float *VClassCS = new float[length];		float *VScoreCS = new float[length];
		float *VClassTL = new float[length];		float *VScoreTL = new float[length];
		float *VClassTS = new float[length];		float *VScoreTS = new float[length];
		float *VClassLS = new float[length];		float *VScoreLS = new float[length];
		Mat ClassScore(length, 5, CV_32F, Scalar(0));
		for (int i = 0; i<length; i++)
		{
			VClassCT[i] = ctSVM.predict(VecCT[i], false); VScoreCT[i] = ctSVM.predict(VecCT[i], true);
			VClassCL[i] = clSVM.predict(VecCL[i], false); VScoreCL[i] = clSVM.predict(VecCL[i], true);
			VClassCS[i] = csSVM.predict(VecCS[i], false); VScoreCS[i] = csSVM.predict(VecCS[i], true);
			VClassTL[i] = tlSVM.predict(VecTL[i], false); VScoreTL[i] = tlSVM.predict(VecTL[i], true);
			VClassTS[i] = tsSVM.predict(VecTS[i], false); VScoreTS[i] = tsSVM.predict(VecTS[i], true);
			VClassLS[i] = lsSVM.predict(VecLS[i], false); VScoreLS[i] = lsSVM.predict(VecLS[i], true);

			ClassScore.at<float>(i, (int)VClassCT[i]) += fabs(VScoreCT[i]); ClassScore.at<float>(i, (int)VClassCL[i]) += fabs(VScoreCL[i]);
			ClassScore.at<float>(i, (int)VClassCS[i]) += fabs(VScoreCS[i]); ClassScore.at<float>(i, (int)VClassTL[i]) += fabs(VScoreTL[i]);
			ClassScore.at<float>(i, (int)VClassTS[i]) += fabs(VScoreTS[i]); ClassScore.at<float>(i, (int)VClassLS[i]) += 0.5*fabs(VScoreLS[i]);
			file << VClassCT[i] << " " << VScoreCT[i] << " " << VClassCL[i] << " " << VScoreCL[i] << " ";
			file << VClassCS[i] << " " << VScoreCS[i] << " " << VClassTL[i] << " " << VScoreTL[i] << " ";
			file << VClassTS[i] << " " << VScoreTS[i] << " " << VClassLS[i] << " " << VScoreLS[i] << " ";
			file << ClassScore.at<float>(i, 1) << " " << ClassScore.at<float>(i, 2) << " "
				<< ClassScore.at<float>(i, 3) << " " << ClassScore.at<float>(i, 4) << endl;
		}
		vector<Mat>(VecCT).swap(VecCT); vector<Mat>(VecCT).swap(VecCL);
		vector<Mat>(VecCT).swap(VecCS); vector<Mat>(VecCT).swap(VecTL);
		vector<Mat>(VecCT).swap(VecTS); vector<Mat>(VecCT).swap(VecLS);
		file.close();

		Mat cloud_mask(height, width, CV_8UC1, Scalar(GC_BGD));
		Mat snow_mask(height, width, CV_8UC1, Scalar(GC_BGD));
		Mat ResultImg(height, width, CV_8UC3, Scalar(0));
		int count = 0;

		for (int i = 0; i < regionList->GetNumRegions(); i++)
		{
			if (isCS[i] == true)
			{
				int *regionindex = regionList->GetRegionPtIndeces(i);
				int pos, x_t, y_t, ClassMax = 1;

				for (int j = 2; j < 5; j++)
				{
					if (ClassScore.at<float>(count, ClassMax) < ClassScore.at<float>(count, j))
					{
						ClassMax = j;
					}
				}

				if (ClassMax == 1)
				{
					for (int j = 0; j < regionPC[i]; j++)//regionPC[i]每个多边形包含点个数
					{
						pos = *(regionindex + j);
						x_t = pos % width;//列数
						y_t = pos / width;//行数
						ResultImg.at<uchar>(y_t, 3 * x_t) = 255;
						ResultImg.at<uchar>(y_t, 3 * x_t + 1) = 255;
						ResultImg.at<uchar>(y_t, 3 * x_t + 2) = 255;
						if (IsInterestArea[y_t][x_t] == true)
						{
							cloud_mask.at<uchar>(y_t, x_t) = GC_FGD;
						}
						else{
							cloud_mask.at<uchar>(y_t, x_t) = GC_PR_FGD;
						}

					}
					count++;
					continue;
				}
				else if (ClassMax == 2)
				{
					for (int j = 0; j < regionPC[i]; j++)//regionPC[i]每个多边形包含点个数
					{
						pos = *(regionindex + j);
						x_t = pos % width;//列数
						y_t = pos / width;//行数
						ResultImg.at<uchar>(y_t, 3 * x_t) = 255;
						ResultImg.at<uchar>(y_t, 3 * x_t + 1) = 255;
						ResultImg.at<uchar>(y_t, 3 * x_t + 2) = 0;
						if (IsInterestArea[y_t][x_t] == true)
						{
							cloud_mask.at<uchar>(y_t, x_t) = GC_FGD;
						}
						else{
							cloud_mask.at<uchar>(y_t, x_t) = GC_PR_FGD;
						}
					}
					count++;
					continue;
				}
				else if (ClassMax == 3)
				{
					for (int j = 0; j < regionPC[i]; j++)//regionPC[i]每个多边形包含点个数
					{
						pos = *(regionindex + j);
						x_t = pos % width;//列数
						y_t = pos / width;//行数
						ResultImg.at<uchar>(y_t, 3 * x_t) = 0;
						ResultImg.at<uchar>(y_t, 3 * x_t + 1) = 255;
						ResultImg.at<uchar>(y_t, 3 * x_t + 2) = 255;
						if (IsInterestArea[y_t][x_t] == true)
						{
							snow_mask.at<uchar>(y_t, x_t) = GC_FGD;
						}
						else{
							snow_mask.at<uchar>(y_t, x_t) = GC_PR_FGD;
						}
					}
					count++;
					continue;
				}
				else if (ClassMax == 4)
				{
					for (int j = 0; j < regionPC[i]; j++)//regionPC[i]每个多边形包含点个数
					{
						pos = *(regionindex + j);
						x_t = pos % width;//列数
						y_t = pos / width;//行数
						ResultImg.at<uchar>(y_t, 3 * x_t) = 127;
						ResultImg.at<uchar>(y_t, 3 * x_t + 1) = 127;
						ResultImg.at<uchar>(y_t, 3 * x_t + 2) = 127;
						if (IsInterestArea[y_t][x_t] == true)
						{
							snow_mask.at<uchar>(y_t, x_t) = GC_FGD;
						}
						else{
							snow_mask.at<uchar>(y_t, x_t) = GC_PR_FGD;
						}
					}
					count++;
					continue;
				}
			}
			else{
				int *regionindex = regionList->GetRegionPtIndeces(i);
				int pos, x_t, y_t;
				for (int j = 0; j < regionPC[i]; j++)//regionPC[i]每个多边形包含点个数
				{
					pos = *(regionindex + j);
					x_t = pos % width;//列数
					y_t = pos / width;//行数
					if (IsInterestArea[y_t][x_t] == true)
					{
						cloud_mask.at<uchar>(y_t, x_t) = GC_PR_FGD;
						snow_mask.at<uchar>(y_t, x_t) = GC_PR_FGD;
					}
				}
			}
		}
		count = 0;

		for (int i = 0; i < regionList->GetNumRegions(); i++)
		{
			if (isCS[i] == true)
			{
				int text_x = 0, text_y = 0;
				int *regionindex = regionList->GetRegionPtIndeces(i);
				int pos, x_t, y_t;
				for (int j = 0; j < regionPC[i]; j++)//regionPC[i]每个多边形包含点个数
				{
					pos = *(regionindex + j);
					x_t = pos % width;//列数
					y_t = pos / width;//行数
					text_x += x_t;
					text_y += y_t;
				}
				text_x /= regionPC[i];
				text_y /= regionPC[i];
				Scalar red = Scalar(0, 0, 255);
				stringstream word;
				string temp;
				word << count + 1;
				temp = word.str();
				putText(ResultImg, temp, cvPoint(text_x, text_y), FONT_HERSHEY_COMPLEX, 0.3, red);//在图片中输出字符 
				count++;
			}
		}
		delete  cbgImage_, filtImage_, whiteImage_, segmImage_, boundaries_, regionPts_,regionList,regionLables;
		delete[] isCS;
		mu.unlock();

		string ResultImgDir, ResultTxtDir, SrcImgDir;
		ResultImgDir = Src_img.substr(0, Src_img.length() - 5) + "_ClassU.bmp";
		SrcImgDir = Src_img.substr(0, Src_img.length() - 5) + "_Src.bmp";
		imwrite(ResultImgDir, ResultImg);
		imwrite(SrcImgDir, img11);
		Rect rect(0, 0, width, height);
		Mat bgdModel, fgdModel;
		grabCut(img11, cloud_mask, rect, bgdModel, fgdModel, 3, GC_INIT_WITH_MASK);
		grabCut(img11, snow_mask, rect, bgdModel, fgdModel, 3, GC_INIT_WITH_MASK);
		Mat cloud_mask_t = cloud_mask & 1;
		Mat snow_mask_t = snow_mask & 1;

		Mat cloud_result, snow_result, final_result(height, width, CV_8UC3, Scalar(0));
		threshold(cloud_mask_t, cloud_result, 0, 255, 0);
		threshold(snow_mask_t, snow_result, 0, 255, 0);
		ResultImgDir = Src_img.substr(0, Src_img.length() - 5) + "_Grabcut_cloud.bmp";
		//	imwrite(ResultImgDir, cloud_result);
		ResultImgDir = Src_img.substr(0, Src_img.length() - 5) + "_Grabcut_snow.bmp";
		//	imwrite(ResultImgDir, snow_result);

		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				if (cloud_result.at<uchar>(i, j) == 255 && snow_result.at<uchar>(i, j) != 255)
				{
					final_result.ptr(i)[3 * j] = final_result.ptr(i)[3 * j + 1] = final_result.ptr(i)[3 * j + 2] = 255;
				}
				else if (snow_result.at<uchar>(i, j) == 255 && cloud_result.at<uchar>(i, j) != 255)
				{
					final_result.ptr(i)[3 * j] = final_result.ptr(i)[3 * j + 1] = final_result.ptr(i)[3 * j + 2] = 127;
				}
				else{
					final_result.ptr(i)[3 * j] = final_result.ptr(i)[3 * j + 1] = final_result.ptr(i)[3 * j + 2] = 0;
				}
			}
		}
		ResultImgDir = Src_img.substr(0, Src_img.length() - 5) + "_result.bmp";
		imwrite(ResultImgDir, final_result);

		string GroudTruthDir;
		GroudTruthDir = Src_img.substr(0, Src_img.length() - 5) + "_gt.bmp";
		Mat gtImg;
		gtImg = imread(GroudTruthDir, 1);

		/*精度评价*/
		Mat CloudPreMat(height, width, CV_8UC3, Scalar(0));
		Mat SnowPreMat(height, width, CV_8UC3, Scalar(0));
		if (gtImg.dims != 0)
		{
			unsigned c_TP, c_FP, c_TN, c_FN, s_TP, s_FP, s_TN, s_FN;
			c_TP = c_FP = c_TN = c_FN = s_TP = s_FP = s_TN = s_FN = 0;
			for (int i = 0; i < height; i++)
			{
				for (int j = 0; j < width; j++)
				{
					if (gtImg.ptr(i)[3 * j] == 255 && final_result.ptr(i)[3 * j] == 255)
					{
						c_TP++;
						CloudPreMat.ptr(i)[3 * j] = CloudPreMat.ptr(i)[3 * j + 1] = CloudPreMat.ptr(i)[3 * j + 2] = 255;
					}
					if (gtImg.ptr(i)[3 * j] != 255 && final_result.ptr(i)[3 * j] == 255)
					{
						c_FP++;
						CloudPreMat.ptr(i)[3 * j] = 0;
						CloudPreMat.ptr(i)[3 * j + 1] = 0;
						CloudPreMat.ptr(i)[3 * j + 2] = 255;
					}
					if (gtImg.ptr(i)[3 * j] == 255 && final_result.ptr(i)[3 * j] != 255)
					{
						c_FN++;
						CloudPreMat.ptr(i)[3 * j] = 0;
						CloudPreMat.ptr(i)[3 * j + 1] = 255;
						CloudPreMat.ptr(i)[3 * j + 2] = 255;
					}
					if (gtImg.ptr(i)[3 * j] != 255 && final_result.ptr(i)[3 * j] != 255)
					{
						c_TN++;
					}

					if (gtImg.ptr(i)[3 * j] == 127 && final_result.ptr(i)[3 * j] == 127)
					{
						s_TP++;
						SnowPreMat.ptr(i)[3 * j] = SnowPreMat.ptr(i)[3 * j + 1] = SnowPreMat.ptr(i)[3 * j + 2] = 127;
					}
					if (gtImg.ptr(i)[3 * j] != 127 && final_result.ptr(i)[3 * j] == 127)
					{
						s_FP++;
						SnowPreMat.ptr(i)[3 * j] = 0;
						SnowPreMat.ptr(i)[3 * j + 1] = 0;
						SnowPreMat.ptr(i)[3 * j + 2] = 255;
					}
					if (gtImg.ptr(i)[3 * j] == 127 && final_result.ptr(i)[3 * j] != 127)
					{
						s_FN++;
						SnowPreMat.ptr(i)[3 * j] = 0;
						SnowPreMat.ptr(i)[3 * j + 1] = 255;
						SnowPreMat.ptr(i)[3 * j + 2] = 255;
					}
					if (gtImg.ptr(i)[3 * j] != 127 && final_result.ptr(i)[3 * j] != 127)
					{
						s_TN++;
					}
				}
			}
			ResultImgDir = Src_img.substr(0, Src_img.length() - 5) + "_cloudPre.bmp";
			imwrite(ResultImgDir, CloudPreMat);
			ResultImgDir = Src_img.substr(0, Src_img.length() - 5) + "_snowPre.bmp";
			imwrite(ResultImgDir, SnowPreMat);

			double c_accuracy, c_precision, c_recall, s_accuracy, s_precision, s_recall;
			c_accuracy = c_precision = c_recall = s_accuracy = s_precision = s_recall = 0.0;

			mu.lock();
			if (countPrecision == 0){
				PrecisionTxt.open(PrecisionTxtDir);
				countPrecision++;
			}
			else{
				PrecisionTxt.open(PrecisionTxtDir, ofstream::app);
			}
			PrecisionTxt << *itr << " :" << endl;
			PrecisionTxt << "  " << "Cloud:";
			if ((c_TP + c_FN) != 0){
				c_accuracy = (double)(c_TP + c_TN) / (c_TP + c_FP + c_TN + c_FN);
				c_precision = (double)c_TP / (c_TP + c_FP);
				c_recall = (double)c_TP / (c_TP + c_FN);
				overall_cloud_accuracy += c_accuracy;
				overall_cloud_precision += c_precision;
				overall_cloud_recall += c_recall;
				countPrecisionCloud++;
				PrecisionTxt << c_accuracy << " " << c_precision << " " << c_recall << endl;
			}
			else{
				PrecisionTxt << "NULL" << endl;
			}
			PrecisionTxt << "  " << "Snow :";
			if ((s_TP + s_FN) != 0){
				s_accuracy = (double)(s_TP + s_TN) / (s_TP + s_FP + s_TN + s_FN);
				s_precision = (double)s_TP / (s_TP + s_FP);
				s_recall = (double)s_TP / (s_TP + s_FN);
				overall_snow_accuracy += s_accuracy;
				overall_snow_precision += s_precision;
				overall_snow_recall += s_recall;
				countPrecisionSnow++;
				PrecisionTxt << s_accuracy << " " << s_precision << " " << s_recall << endl;
			}
			else{
				PrecisionTxt << "NULL" << endl;
			}
			PrecisionTxt << endl;
			PrecisionTxt.close();
			mu.unlock();
		}
		else{
			cout << *itr << " lacks Groud truth image." << endl;
		}

		/***                     end                     ***/
		cout << endl << endl;
		delete[] pImgdata; delete[] pImgdata1; delete[] img_src; 
		delete[] VClassCT; delete[] VClassCL; delete[] VClassCS;
		delete[] VClassTL; delete[] VClassTS; delete[] VClassLS;
		for (int i = 0; i < height; i++)
		{
			delete[]IsInterestArea[i];
			delete[]HOT[i];
			delete[]SF[i];
		}
		delete[] IsInterestArea;
		delete[] HOT;
		delete[] SF;
		
		GDALClose(poDataset);	
	}
}

int main()
{
	vector<string> All_imgfile;	
	GetTiffFile(Src_dir, &All_imgfile);
	/***                      SVM training                     ***/
	ct_svm_train();
	cl_svm_train();
	cs_svm_train();
	tl_svm_train();
	ts_svm_train();
	ls_svm_train();

	double t1 = GetTickCount();
	thread t01(CloudSnow_Dtec, All_imgfile, 0);
	thread t02(CloudSnow_Dtec, All_imgfile, 1);
	t01.join();
	t02.join();

	PrecisionTxt.open(PrecisionTxtDir, ofstream::app);
	PrecisionTxt << "Overall precision:" << endl;
	PrecisionTxt << "  Cloud:" << overall_cloud_accuracy / countPrecisionCloud << " " << overall_cloud_precision / countPrecisionCloud << " " << overall_cloud_recall / countPrecisionCloud << " " << countPrecisionCloud << endl;
	PrecisionTxt << "  Snow :" << overall_snow_accuracy / countPrecisionSnow << " " << overall_snow_precision / countPrecisionSnow << " " << overall_snow_recall / countPrecisionSnow << " " << countPrecisionSnow << endl;
	PrecisionTxt << "time cost：" <<(double)(GetTickCount() - t1) / 1000 << " seconds" << endl;
	PrecisionTxt.close();

	
	return 0;
}

/*计算样本特征值*/
//int main() 
//{
//	int count;
//	vector<string> All_SampleFile;
//	const string Sample_dir("C:\\Users\\Kylin\\OneDrive\\CS_sample\\snow");
//	string Smp_img;
//	count=GetTiffFile(Sample_dir, &All_SampleFile);
//	
//	ofstream file;
//	if (Sample_dir.substr(Sample_dir.length() - 4, 4) != "loud")
//	{
//		file.open("svmfeature.txt", ofstream::app);
//		
//	}
//	else{
//		file.open("svmfeature.txt");
//	}
//
//	/*以下是对每块样本分别提取FeatureNum维向量 */
//	for (vector<string>::iterator itr = All_SampleFile.begin(); itr != All_SampleFile.end(); itr++)
//	{
//		SmpFeature SmpFeature;
//		Smp_img = Sample_dir + "\\" + *itr;
//		int nband;
//		GDALAllRegister();
//		GDALDataset *poDataset = NULL;
//		poDataset = (GDALDataset*)GDALOpen(Smp_img.c_str(), GA_ReadOnly);
//		if (poDataset == NULL)
//		{
//			GDALClose(poDataset);
//			return 0;
//		}
//		width = poDataset->GetRasterXSize();
//		height = poDataset->GetRasterYSize();
//		nband = poDataset->GetRasterCount();
//		float *pSmpdata = new float[width*height*nband];
//		poDataset->RasterIO(
//			GF_Read, 0, 0, width, height, pSmpdata, width, height, GDT_Float32, nband, NULL, 0, 0, 0);
//	
//		for (int i = 0; i < width*height; i++)
//		{
//			if (pSmpdata[i] != 0 || pSmpdata[width*height + i] != 0 || pSmpdata[width*height * 2 + i] != 0 || pSmpdata[width*height * 3 + i] != 0)
//			{
//				SmpFeature.area_pnum++;
//				SmpFeature.NIR += pSmpdata[width*height * 3 + i];
//				SmpFeature.R += pSmpdata[width*height * 2 + i];
//				SmpFeature.G += pSmpdata[width*height * 1 + i];
//				SmpFeature.B += pSmpdata[width*height * 0 + i];
//			}
//		}
//		SmpFeature.NIR = SmpFeature.NIR / SmpFeature.area_pnum;
//		SmpFeature.R = SmpFeature.R / SmpFeature.area_pnum;
//		SmpFeature.G = SmpFeature.G / SmpFeature.area_pnum;
//		SmpFeature.B = SmpFeature.B / SmpFeature.area_pnum;
//
//		file << SmpFeature.B << " " << SmpFeature.G << " " << SmpFeature.R << " " << SmpFeature.NIR << " ";
//
//		float v2r, v2g, v2b, v2nir;
//		v2r = v2g = v2b = v2nir = 0;
//
//		for (int i = 0; i < width*height; i++)
//		{
//			if (pSmpdata[i] != 0 || pSmpdata[width*height + i] != 0 || pSmpdata[width*height * 2 + i] != 0 || pSmpdata[width*height * 3 + i] != 0)
//			{
//				v2nir += (pSmpdata[width*height * 3 + i] - SmpFeature.NIR)* (pSmpdata[width*height * 3 + i] - SmpFeature.NIR);
//				v2r += (pSmpdata[width*height * 2 + i] - SmpFeature.R)*(pSmpdata[width*height * 2 + i] - SmpFeature.R);
//				v2g += (pSmpdata[width*height * 1 + i] - SmpFeature.G)*(pSmpdata[width*height * 1 + i] - SmpFeature.G);
//				v2b += (pSmpdata[width*height * 0 + i] - SmpFeature.B)*(pSmpdata[width*height * 0 + i] - SmpFeature.B);
//			}
//		}
//		SmpFeature.v2nir =sqrt( v2nir / SmpFeature.area_pnum);
//		SmpFeature.v2r = sqrt(v2r / SmpFeature.area_pnum);
//		SmpFeature.v2g = sqrt(v2g / SmpFeature.area_pnum);
//		SmpFeature.v2b = sqrt(v2b / SmpFeature.area_pnum);
//		GDALClose(poDataset);
//
//		SmpFeature.mgray_b = new int*[height];
//		SmpFeature.mgray_g = new int*[height];
//		SmpFeature.mgray_r = new int*[height];
//		SmpFeature.mgray_nr = new int*[height];
//		SmpFeature.mgray_h = new double*[height];
//		SmpFeature.mgray_i = new double*[height];
//		SmpFeature.mgray_s = new double*[height];
//		double maxI, minI, maxS, minS;
//		int tempCount = 0;
//		for (int i = 0; i < height;i++)
//		{
//			SmpFeature.mgray_b[i] = new int[width];
//			SmpFeature.mgray_g[i] = new int[width];
//			SmpFeature.mgray_r[i] = new int[width];
//			SmpFeature.mgray_nr[i] = new int[width];
//			SmpFeature.mgray_h[i] = new double[width];
//			SmpFeature.mgray_i[i] = new double[width];
//			SmpFeature.mgray_s[i] = new double[width];
//			for (int j = 0; j < width;j++)
//			{
//				SmpFeature.mgray_b[i][j] = pSmpdata[i*width + j];
//				SmpFeature.mgray_g[i][j] = pSmpdata[width*height+i*width + j];
//				SmpFeature.mgray_r[i][j] = pSmpdata[width*height*2+i*width + j];
//				SmpFeature.mgray_nr[i][j] = pSmpdata[width*height * 3+i * width + j];
//				transRGBtoHIS(SmpFeature.mgray_r[i][j], SmpFeature.mgray_g[i][j], SmpFeature.mgray_b[i][j],
//					SmpFeature.mgray_h[i][j], SmpFeature.mgray_i[i][j], SmpFeature.mgray_s[i][j]);
//				if (SmpFeature.mgray_b[i][j] != 0 || SmpFeature.mgray_g[i][j] != 0 || SmpFeature.mgray_r[i][j] != 0)
//				{
//					if (tempCount == 0)
//					{
//						maxI = minI = SmpFeature.mgray_i[i][j];
//						maxS = minS = SmpFeature.mgray_s[i][j];
//						tempCount++;
//					}
//					else{
//						if (maxI < SmpFeature.mgray_i[i][j])maxI = SmpFeature.mgray_i[i][j];
//						if (minI > SmpFeature.mgray_i[i][j])minI = SmpFeature.mgray_i[i][j];
//						if (maxS < SmpFeature.mgray_s[i][j])maxS = SmpFeature.mgray_s[i][j];
//						if (minS > SmpFeature.mgray_s[i][j])minS = SmpFeature.mgray_s[i][j];
//					}
//				}
//			}
//		}
//
//		/*cal SF image*/
//		tempCount = 0;
//		double minSF, maxSF;
//		SmpFeature.SF = new double*[height];
//		for (int i = 0; i < height; i++)
//		{
//			SmpFeature.SF[i] = new double[width];
//			for (int j = 0; j < width; j++)
//			{
//				SmpFeature.SF[i][j] = 0;
//				if (SmpFeature.mgray_b[i][j] != 0 || SmpFeature.mgray_g[i][j] != 0 || SmpFeature.mgray_r[i][j] != 0)
//				{
//					SmpFeature.mgray_i[i][j] = (SmpFeature.mgray_i[i][j] - minI) / (maxI - minI);
//					SmpFeature.mgray_s[i][j] = (SmpFeature.mgray_s[i][j] - minS) / (maxS - minS);
//					SmpFeature.SF[i][j] = (SmpFeature.mgray_i[i][j] + 1.0) / (SmpFeature.mgray_s[i][j] + 1.0);
//					if (tempCount == 0)
//					{
//						maxSF = minSF = SmpFeature.SF[i][j];
//						tempCount++;
//					}
//					else{
//						if (maxSF < SmpFeature.SF[i][j])maxSF = SmpFeature.SF[i][j];
//						if (minSF > SmpFeature.SF[i][j])minSF = SmpFeature.SF[i][j];
//					}
//				}
//			}
//		}
//		tempCount = 0;
//		double valueST=0;
//		for (int i = 0; i < height; i++)
//		{
//			for (int j = 0; j < width; j++)
//			{
//				if (SmpFeature.mgray_b[i][j] != 0 || SmpFeature.mgray_g[i][j] != 0 || SmpFeature.mgray_r[i][j] != 0)
//				{
//					SmpFeature.SF[i][j] = (SmpFeature.SF[i][j] - minSF) / (maxSF - minSF)*255.0;
//					valueST += SmpFeature.mgray_nr[i][j] - (1.73*SmpFeature.SF[i][j] - 219.2);
//					tempCount++;
//				}
//			}
//		}
//		valueST /= tempCount;
//
//
//		//double b_FractalD, g_FractalD, r_FractalD, nr_FractalD;
//		//calFractalDimension(SmpFeature.mgray_b, height, width, b_FractalD);
//		//calFractalDimension(SmpFeature.mgray_g, height, width, g_FractalD);
//		//calFractalDimension(SmpFeature.mgray_r, height, width, r_FractalD);
//		//calFractalDimension(SmpFeature.mgray_nr, height, width, nr_FractalD);
//		//SmpFeature.fractalDimension = (b_FractalD + g_FractalD + r_FractalD + nr_FractalD) / 4;
//		//file << SmpFeature.fractalDimension << " ";
//		file << SmpFeature.v2b << " " << SmpFeature.v2g << " " << SmpFeature.v2r << " " << SmpFeature.v2nir << " ";
//		file << valueST << " ";
//
//		/*cal LTP*/
//		Cal_IWCS_LTP(SmpFeature.mgray_r, height, width, SmpFeature.LTP_features);
//		for (int i = 0; i < 18; i++){ file << SmpFeature.LTP_features[i] << " "; }
//
//		calColorGLCM(pSmpdata, height, width, SmpFeature.CGLCM_features);
//		for (int i = 0; i < 9; i++)
//		{
//			file << SmpFeature.CGLCM_features[i] << " ";
//		}
//
//		float CurHistogram[10];
//		calEdgeCurHis(SmpFeature.mgray_b, height, width, CurHistogram);
//	    /*	for (int i = 0; i < 10; i++)
//		{
//			file << CurHistogram[i] << " ";
//		}*/
//		file << CurHistogram[0] + CurHistogram[1] + CurHistogram[2] << " ";
//
//		double aValueMCLFD[8*4];
//		calMCLFD(SmpFeature.mgray_b, height, width, aValueMCLFD);
//		for (int i = 0; i < 8; i++)
//		{
//			file << aValueMCLFD[i*4] << " ";
//		}
//
//		delete[]pSmpdata;
//		file << endl; 
//	}
//	file.close();
//
//
//	if (Sample_dir.substr(Sample_dir.length() - 4, 4) == "snow")
//	{
//		ifstream file1;
//		file1.open("svmfeature.txt");
//		float **OFeatureD;
//		int nNum=410;
//		OFeatureD = new float*[nNum];
//		for (int i = 0; i < nNum; i++)
//		{
//			OFeatureD[i] = new float[FeatureNum];
//			for (int j = 0; j < FeatureNum; j++)
//			{
//				file1 >> OFeatureD[i][j];
//			}
//		}
//		file1.close();
//		ofstream file2;
//		file2.open("NormalParams.txt");
//		for (int i = 0; i < FeatureNum; i++)
//		{
//			float average = 0;
//			float Sd=0;
//			for (int j = 0; j < nNum; j++)
//			{
//				average += OFeatureD[j][i];
//			}
//			average /= nNum;
//			for (int j = 0; j < nNum; j++)
//			{
//				Sd += (OFeatureD[j][i] - average)*(OFeatureD[j][i] - average);
//			}
//			Sd = sqrt(Sd / nNum);
//			file2 << average << " " << Sd << endl;
//			for (int j = 0; j < nNum; j++)
//			{
//				OFeatureD[j][i] = (OFeatureD[j][i] - average) / Sd;
//			}
//		}
//		file2.close();
//		file2.open("svmfeatureU.txt");
//		for (int i = 0; i < nNum; i++)
//		{
//			for (int j = 0; j < FeatureNum; j++)
//			{
//				file2 << OFeatureD[i][j] << " ";
//			}
//			file2 << endl;
//		}
//		file2.close();
//	}
//	return 0;
//}
