#include "gdal_priv.h"
#include "gdal.h"
#include <iostream>
#include<opencv2/opencv.hpp>  
using namespace std;
using namespace cv;

int main()
{
	int width, height,nband;
	GDALAllRegister();
	GDALDataset *poDataset = NULL;
	poDataset = (GDALDataset*)GDALOpen("E:\\CloudDtec\\CSDiscrimination\\data\\sample\\1.tiff", GA_ReadOnly);
	if (poDataset == NULL)
	{
		GDALClose(poDataset);
		return false;
	}
	width = poDataset->GetRasterXSize();
	height = poDataset->GetRasterYSize();
	nband = poDataset->GetRasterCount();
	cout << "Drivers:" << poDataset->GetDriver()->GetDescription() << poDataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME) << endl;
	cout << "Size is " << poDataset->GetRasterXSize()<<"," << poDataset->GetRasterYSize()<< "," << poDataset->GetRasterCount() << endl;

	double        adfGeoTransform[6];
	if (poDataset->GetGeoTransform(adfGeoTransform) == CE_None)
		    {
		         printf("Origin = (%.6f,%.6f)\n",
			                adfGeoTransform[0], adfGeoTransform[3]);
			        printf("Pixel Size = (%.6f,%.6f)\n",
			                 adfGeoTransform[1], adfGeoTransform[5]);
		    }

	float *pImgdata;
	int *pImgdataband1;
	pImgdata = NULL;
	pImgdata = new float[width*height*nband];
	pImgdataband1 = new int[width*height];


	poDataset->RasterIO(
		GF_Read, 0, 0, width, height, pImgdata, width, height, GDT_Float32, nband, NULL, 0, 0, 0);

 	for (int i = 0; i < width*height; i++)
	{
		if (pImgdata[i] - 0.5*pImgdata[width*height * 2 + i] - 0.08>0)
		{
		/*	pImgdataband1[i] = 255;
			pImgdataband1[i*3 +1] = 255;
			pImgdataband1[i*3 + 2] = 255;*/
		
		}
		else
		{
			pImgdataband1[i] = 0;
			pImgdataband1[i * 3 + 1] = 0;
			pImgdataband1[i * 3 + 2] = 0;
		}
	}

	GDALDriver*poDriver = GetGDALDriverManager()->GetDriverByName("BMP");
	if (poDriver == NULL)
	{
		return false;
	}

	char **papszOptions = NULL;
	GDALDataset *WriteDataSet = poDriver->Create("E:\\CloudDtec\\CSDiscrimination\\data\\sample\\band2.bmp",
		width, height, 1, GDT_Byte, NULL);
	WriteDataSet->RasterIO(GF_Write, 0, 0, width, height, pImgdataband1, width, height, GDT_Byte, 1, NULL, 0, 0, 0);
	
	return 0;

}

