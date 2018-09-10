/*
 * Graph.h
 *
 *  Created on: Sep 3, 2018
 *      Author: oscar
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include "itkImage.h"
#include "common/Type.h"
#include "itkImageRegionIterator.h"
#include "common/Math.h"
#include "itkRGBToLuminanceImageFilter.h"
#include "common/ImageProcessing.h"
#include "common/InputOutput.h"
#include "itkCastImageFilter.h"

extern "C"
{
#include <igraph.h>
}

class Graph
{

private:
	//typedefs

	using componentType = unsigned char;
	using rgbPixelType = itk::RGBPixel< componentType >;
	using rgbImageType = itk::Image< rgbPixelType, 2 >;

	using grayType = double;
	using grayImageType = itk::Image<double,2>;
	using grayIteratorType = itk::ImageRegionIterator<grayImageType>;
	using grayImagePointer = typename grayImageType::Pointer;
	using grayIndexType = grayImageType::IndexType;

public:
	Graph();
	~Graph()
	{
	}
	;

	void setImage(rgbImageType::Pointer rgbImage);
	void buildOld();
	void build();

private:

	grayImagePointer image;

	//Auxiliary matrices to create the igraph, because it is faster than itkImage
	std::vector<std::vector<double> > imageIntensity;

	double laplaceWeigth;
	unsigned radius;
	igraph_t graph;

	//methods
	inline double laplacianWeight(grayIndexType currentIndex, grayIndexType neighborIndex, grayIteratorType currentIt, grayIteratorType neighborIt, double b);
	inline double stimateParameterB(grayIteratorType imageIt, grayImageType::IndexType currentIndex, grayIteratorType neighborhoodIt);

	void segmentBackground();

};

Graph::Graph()
{
	laplaceWeigth = 2.0;
	radius = 1;

}

/*
 * Separate background and foreground, background is set to 0.
 *
 * */

void Graph::segmentBackground()
{

	grayImagePointer segmented = ip::otsuThreshold<grayImageType>(image);

	grayIteratorType imIt(image, image->GetLargestPossibleRegion()); //image iterator
	grayIteratorType segIt(segmented, image->GetLargestPossibleRegion()); //segmented image iterator

	for (imIt.GoToBegin(), segIt.GoToBegin(); !imIt.IsAtEnd(); ++imIt, ++segIt)
	{

		if (segIt.Get() == 0)
		{
			imIt.Set(0);

		}

	}

	/*typedef itk::CastImageFilter<grayImageType, rgbImageType > CastFilterType;
	 CastFilterType::Pointer castFilter = CastFilterType::New();
	 castFilter->SetInput(image);
	 castFilter->Update();
	 io::writeImage<rgbImageType>(castFilter->GetOutput(), "/home/oscar/MEGA/post-doc/src/input/laudos/otsu.png");
	 */

}

void Graph::setImage(rgbImageType::Pointer rgbImage)
{

	using luminanceImageFilterType = itk::RGBToLuminanceImageFilter< rgbImageType, grayImageType >;
	luminanceImageFilterType::Pointer luminanceImageFilter = luminanceImageFilterType::New();
	luminanceImageFilter->SetInput(rgbImage);

	luminanceImageFilter->Update();
	this->image = luminanceImageFilter->GetOutput();

	segmentBackground();

}

inline double Graph::stimateParameterB(grayIteratorType imageIt, grayImageType::IndexType currentIndex, grayIteratorType neighborhoodIt)
{

	type::grayImage::SizeType size = image->GetLargestPossibleRegion().GetSize(); //aux

	double sumIntensity = 0.0;
	double sumIndex = 0.0;
	unsigned count = 0;

	neighborhoodIt.GoToBegin();

	grayImageType::IndexType loweIndex = neighborhoodIt.GetRegion().GetIndex();
	grayImageType::IndexType upperIndex = neighborhoodIt.GetRegion().GetUpperIndex();

	for (unsigned i = loweIndex[0]; i <= upperIndex[0]; ++i)
	{
		for (unsigned j = loweIndex[1]; j <= upperIndex[1]; ++j, ++neighborhoodIt)
		{
			//if()

			sumIntensity += std::abs(neighborhoodIt.Get() - imageIt.Get());
			sumIndex += std::sqrt(std::pow(static_cast<double>(i - currentIndex[0]), 2) + std::pow(static_cast<double>(j - currentIndex[1]), 2));
			++count;

		}

	}

	double b = (sumIndex + sumIntensity) / static_cast<double>(count);

	return (b == 0) ? 0.0000000000001 : b; //to avoid division by zero*/
}

inline double Graph::laplacianWeight(grayIndexType currentIndex, grayIndexType neighborIndex, grayIteratorType currentIt, grayIteratorType neighborIt, double b)
{

	double indexTerm;
	double intensityTerm;

	indexTerm = math::euclideanDistance<grayIndexType>(neighborIndex, currentIndex, 2);
	intensityTerm = std::abs(neighborIt.Get() - currentIt.Get());

	double numerator = std::abs(intensityTerm + indexTerm);

	//std::cin>>b;
	return (1.0 / (2.0 * b)) * std::exp(-laplaceWeigth * numerator / b); //2.0

}

void Graph::build()
{

	grayIteratorType imageIt(image, image->GetLargestPossibleRegion());

	grayImageType::RegionType neighborhood;

	grayImageType::IndexType lowerIndex;
	grayImageType::IndexType upperIndex;
	grayImageType::IndexType currentIndex;
	grayImageType::IndexType neighborIndex;

	grayImageType::SizeType imageSize = image->GetLargestPossibleRegion().GetSize();
	std::cout << imageSize << std::endl;

	double b;
	double weight;
	unsigned vectorIndex;

	igraph_vector_t weights;
	igraph_vector_t edges;
	igraph_vector_init(&edges, 0);
	igraph_vector_init(&weights, 0);

	//imageIt.GetIndex() is a very expensive operation,
	//size the matrix position is used, the traditional
	//two for loops strategy is faster than computing imageIt
	imageIt.GoToBegin();
	for (unsigned width = 0; width < imageSize[0]; ++width)
	{
		currentIndex[0] = width;

		lowerIndex[0] = (static_cast<signed>(width - radius) < 0) ? 0 : width - radius;
		upperIndex[0] = (width + radius >= imageSize[0]) ? imageSize[0] - 1 : width + radius;

		for (unsigned height = 0; height < imageSize[1]; ++height, ++imageIt)
		{
			currentIndex[1] = height;
			lowerIndex[1] = (static_cast<signed>(height - radius) < 0) ? 0 : height - radius;
			upperIndex[1] = (height + radius >= imageSize[1]) ? imageSize[1] - 1 : height + radius;

			neighborhood.SetIndex(lowerIndex);
			neighborhood.SetUpperIndex(upperIndex);

			grayIteratorType neighborIt(image, neighborhood);

			//Estimating parameter b
			b = stimateParameterB(imageIt, currentIndex, neighborIt);
			neighborIt.GoToBegin();
			vectorIndex = height * imageSize[0] + width;
			for (unsigned i = lowerIndex[0]; i <= upperIndex[0]; ++i)
			{
				neighborIndex[0] = i;
				for (unsigned j = lowerIndex[1]; j <= upperIndex[1]; ++j, ++neighborIt)
				{
					neighborIndex[1] = j;
					weight = laplacianWeight(currentIndex, neighborIndex, imageIt, neighborIt, b);

					igraph_vector_resize(&edges, igraph_vector_size(&edges) + 2);
					igraph_vector_set(&edges, igraph_vector_size(&edges) - 1, vectorIndex);
					igraph_vector_set(&edges, igraph_vector_size(&edges) - 2, j * imageSize[0] + i);

					igraph_vector_resize(&weights, igraph_vector_size(&weights) + 1);
					igraph_vector_set(&weights, igraph_vector_size(&weights) - 1, weight);

				}
			}

		}

	}
	igraph_create(&graph, &edges, 0, 0);
}

#endif /* GRAPH_H_ */
