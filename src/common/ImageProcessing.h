/*
 * ImageProcessing.h
 *
 *  Created on: Sep 5, 2018
 *      Author: oscar
 */

#ifndef SRC_COMMON_IMAGEPROCESSING_H_
#define SRC_COMMON_IMAGEPROCESSING_H_


#include "itkOtsuThresholdImageFilter.h"



namespace ip
{

template<typename imageType>
typename imageType::Pointer otsuThreshold(typename imageType::Pointer inputImage, unsigned insideValue = 255, unsigned outsideValue = 0, unsigned bins=250)
{

	typedef itk::OtsuThresholdImageFilter<imageType, imageType> FilterType;
	typename FilterType::Pointer otsuFilter = FilterType::New();
	otsuFilter->SetInput(inputImage);

	otsuFilter->SetInsideValue(insideValue); //0
	otsuFilter->SetOutsideValue(outsideValue); //1

	otsuFilter->SetNumberOfHistogramBins(10);

	otsuFilter->Update();

	//IO::print("Otsu Threshold", 1);

	return otsuFilter->GetOutput();

}


}



#endif /* SRC_COMMON_IMAGEPROCESSING_H_ */
