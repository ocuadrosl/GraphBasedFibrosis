/*
 * ImageProcessing.h
 *
 *  Created on: Sep 5, 2018
 *      Author: oscar
 */

#ifndef SRC_COMMON_IMAGEPROCESSING_H_
#define SRC_COMMON_IMAGEPROCESSING_H_

#include "itkOtsuThresholdImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkKittlerIllingworthThresholdImageFilter.h"
#include "itkHistogramThresholdImageFilter.h"
#include "itkHuangThresholdImageFilter.h"
#include "itkLiThresholdImageFilter.h"
#include "itkIntermodesThresholdImageFilter.h"
#include "itkIsoDataThresholdImageFilter.h"
#include "itkMaximumEntropyThresholdImageFilter.h"
#include "itkMomentsThresholdImageFilter.h"
#include "itkRenyiEntropyThresholdImageFilter.h"
#include "itkShanbhagThresholdImageFilter.h"
#include "itkTriangleThresholdImageFilter.h"
#include "itkYenThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

namespace ip
{

template<typename imageType>
typename imageType::Pointer simpleThreshold(typename imageType::Pointer inputImage, unsigned lowerThreshold = 245, unsigned upperThreshold = 255, unsigned insideValue = 255, unsigned outsideValue =
		0)
{

	using BinaryThresholdImageFilterType = itk::BinaryThresholdImageFilter<imageType, imageType>;

	typename BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
	thresholdFilter->SetInput(inputImage);
	thresholdFilter->SetLowerThreshold(lowerThreshold);
	thresholdFilter->SetUpperThreshold(upperThreshold);
	thresholdFilter->SetInsideValue(insideValue);
	thresholdFilter->SetOutsideValue(outsideValue);
	thresholdFilter->Update();
	return thresholdFilter->GetOutput();

}

template<typename imageType>
typename imageType::Pointer histogramThreshold(typename imageType::Pointer inputImage, const std::string& algorithmName = "otsu", unsigned bins = 100, unsigned insideValue = 255,
		unsigned outsideValue = 0)
{

	typename itk::HistogramThresholdImageFilter<imageType, imageType>::Pointer filter;

	if (algorithmName == "otsu")
	{
		filter = itk::OtsuThresholdImageFilter<imageType, imageType>::New();

	}
	if (algorithmName == "huang")
	{
		filter = itk::HuangThresholdImageFilter<imageType, imageType>::New();

	}
	if (algorithmName == "intermodes")
	{
		filter = itk::IntermodesThresholdImageFilter<imageType, imageType>::New();

	}

	if (algorithmName == "isoData")
	{
		filter = itk::IsoDataThresholdImageFilter<imageType, imageType>::New();

	}
	if (algorithmName == "li")
	{
		filter = itk::LiThresholdImageFilter<imageType, imageType>::New();

	}
	if (algorithmName == "maximumEntropy")
	{
		filter = itk::MaximumEntropyThresholdImageFilter<imageType, imageType>::New();

	}
	if (algorithmName == "moments")
	{
		filter = itk::MomentsThresholdImageFilter<imageType, imageType>::New();

	}
	if (algorithmName == "renyi")
	{
		filter = itk::RenyiEntropyThresholdImageFilter<imageType, imageType>::New();

	}
	if (algorithmName == "shanbhag")
	{
		filter = itk::ShanbhagThresholdImageFilter<imageType, imageType>::New();

	}
	if (algorithmName == "triangle")
	{
		filter = itk::TriangleThresholdImageFilter<imageType, imageType>::New();

	}
	if (algorithmName == "yen")
	{
		filter = itk::YenThresholdImageFilter<imageType, imageType>::New();

	}

	filter->SetInput(inputImage);

	filter->SetInsideValue(insideValue); //255
	filter->SetOutsideValue(outsideValue); //0

	filter->SetNumberOfHistogramBins(bins);

	filter->Update();

	return filter->GetOutput();

}

template<typename imageType>
typename imageType::Pointer kittlerIllingworthThreshold(typename imageType::Pointer inputImage, unsigned bins = 256, unsigned insideValue = 255, unsigned outsideValue = 0)
{

	using KittlerIllingworthFilterType = itk::KittlerIllingworthThresholdImageFilter<imageType, imageType >;

	typename KittlerIllingworthFilterType::Pointer kittlerFilter = KittlerIllingworthFilterType::New();
	kittlerFilter->SetInput(inputImage);

	kittlerFilter->SetInsideValue(insideValue); //255
	kittlerFilter->SetOutsideValue(outsideValue); //0

	kittlerFilter->SetNumberOfHistogramBins(bins);

	kittlerFilter->Update();

	//IO::print("Otsu Threshold", 1);

	return kittlerFilter->GetOutput();

}

template<typename imageType>
typename imageType::Pointer otsuMultipleThreshold(typename imageType::Pointer inputImage, unsigned insideValue = 255, unsigned outsideValue = 0, unsigned bins = 256)
{

	typedef itk::OtsuMultipleThresholdsImageFilter<imageType, imageType> FilterType;
	typename FilterType::Pointer otsuFilter = FilterType::New();
	otsuFilter->SetInput(inputImage);

	//otsuFilter->SetInsideValue(insideValue); //255
	//otsuFilter->SetOutsideValue(outsideValue); //0

	otsuFilter->SetNumberOfHistogramBins(bins);
	otsuFilter->SetNumberOfThresholds(2);
	otsuFilter->SetLabelOffset(2);

	using RescaleType = itk::RescaleIntensityImageFilter< imageType, imageType >;
	typename RescaleType::Pointer rescaler = RescaleType::New();
	rescaler->SetInput(otsuFilter->GetOutput());
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);

	rescaler->Update();

	return rescaler->GetOutput();

}

template<typename imageType>
typename imageType::Pointer otsuThreshold(typename imageType::Pointer inputImage, unsigned bins = 256, unsigned insideValue = 255, unsigned outsideValue = 0)
{

	typedef itk::OtsuThresholdImageFilter<imageType, imageType> FilterType;
	typename FilterType::Pointer otsuFilter = FilterType::New();
	otsuFilter->SetInput(inputImage);

	otsuFilter->SetInsideValue(insideValue); //255
	otsuFilter->SetOutsideValue(outsideValue); //0

	otsuFilter->SetNumberOfHistogramBins(bins);

	otsuFilter->Update();

	//IO::print("Otsu Threshold", 1);

	return otsuFilter->GetOutput();

}

template<typename imageType>
typename imageType::Pointer imageGradient(typename imageType::Pointer image)
{

	using FilterType = itk::GradientMagnitudeImageFilter<imageType, imageType>;
	typename FilterType::Pointer filter = FilterType::New();
	filter->SetInput(image);
	filter->Update();
	return filter->GetOutput();

}

}

#endif /* SRC_COMMON_IMAGEPROCESSING_H_ */
