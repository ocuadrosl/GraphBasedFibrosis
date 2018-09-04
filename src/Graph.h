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

extern "C"
{
#include <igraph.h>
}
class Graph
{

private:
	//typedefs
	typedef itk::ImageRegionIterator<type::grayImage> iteratorType;
public:
	Graph();
	~Graph()
	{
	}
	;

	void setImage(type::grayImagePointer image);
	void buildOld();
	void build();

private:

	type::grayImagePointer image;

	//Auxiliary matrices,
	std::vector<std::vector<double> > imageIntensity;
	//std::vector<std::vector<std::vector<type::grayImage::IndexType> > > imageIndex;

	double laplaceWeigth;
	unsigned radius;
	igraph_t graph;

	//methods
	inline double laplacianWeigh(const iteratorType& it1, const iteratorType& it2);
	inline double laplacianWeigh(const std::vector<unsigned>& index1, const std::vector<unsigned>& index2);
	double stimateParameterB(iteratorType inputIt);

};

Graph::Graph()
{
	laplaceWeigth = 2.0;
	radius = 1;

}

void Graph::setImage(type::grayImagePointer image)
{
	this->image = image;

	iteratorType it(image, image->GetLargestPossibleRegion());

	type::grayImage::SizeType imageSize = image->GetLargestPossibleRegion().GetSize();

	std::vector<double> rows(imageSize[1], 0.0);
	std::vector<std::vector<double>> imageAux(imageSize[0], rows);

	unsigned i = 0;
	for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++i)
	{
		type::grayImage::IndexType index = it.GetIndex();
		imageAux[index[0]][index[1]] = it.Get();
	}

	imageIntensity = imageAux;

}

double Graph::stimateParameterB(iteratorType inputIt)
{

	type::grayImage::RegionType neighborhood;

	type::grayImage::IndexType index = inputIt.GetIndex(); //lower index

	index[0] = (index[0] - radius < 0) ? 0 : index[0] - radius;
	index[1] = (index[1] - radius < 0) ? 0 : index[1] - radius;

	type::grayImage::IndexType upper = inputIt.GetIndex();

	type::grayImage::SizeType size = inputIt.GetRegion().GetSize(); //aux

	upper[0] = (upper[0] + radius >= size[0]) ? size[0] - 1 : upper[0] + radius;
	upper[1] = (upper[1] + radius >= size[1]) ? size[1] - 1 : upper[1] + radius;

	neighborhood.SetIndex(index);
	neighborhood.SetUpperIndex(upper);

	iteratorType it(image, neighborhood);

	double sumIntensity = 0.0;
	double sumIndex = 0.0;
	double intensity = inputIt.Get();
	index = inputIt.GetIndex();

	unsigned count = 0;

	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		sumIntensity += std::abs(it.Get() - intensity);
		sumIndex += math::euclideanDistance<type::grayImage::IndexType>(it.GetIndex(), index);
		++count;
	}

	double b = (sumIndex + sumIntensity) / static_cast<double>(count);

	return (b == 0) ? 0.0000000000001 : b; //to avoid division by zero

}

inline double Graph::laplacianWeigh(const iteratorType& it1, const iteratorType& it2)
{

	double indexTerm;
	double intensityTerm;

	indexTerm = math::euclideanDistance<type::rgbImage::IndexType>(it1.GetIndex(), it2.GetIndex());
	intensityTerm = math::euclideanDistance<type::rgbImage::PixelType>(it1.Get(), it2.Get());

	double numerator = std::abs(intensityTerm + indexTerm);

	double b = 1;

	return (1.0 / (2.0 * b)) * std::exp(-laplaceWeigth * numerator / b); //2.0

}

inline double Graph::laplacianWeigh(const std::vector<unsigned>& index1, const std::vector<unsigned>& index2)
{

	double indexTerm;
	double intensityTerm;

	indexTerm = math::euclideanDistance<std::vector<unsigned>>(index1, index2, 2);
	intensityTerm = std::abs(imageIntensity[index1[0]][index1[1]] - imageIntensity[index2[0]][index2[1]]);

	double numerator = std::abs(intensityTerm + indexTerm);

	double b = 1;

	return (1.0 / (2.0 * b)) * std::exp(-laplaceWeigth * numerator / b); //2.0

}

void Graph::build()
{

	type::grayImage::SizeType size = image->GetLargestPossibleRegion().GetSize();

	for (unsigned i = 0; i < size[0]; ++i)
	{

		for (unsigned j = 0; j < size[1]; ++j)
		{
			//get neighbors



		}
	}

}

void Graph::buildOld()
{

	iteratorType it1(image, image->GetLargestPossibleRegion());
	iteratorType it2(image, image->GetLargestPossibleRegion());

	for (it1.GoToBegin(); !it1.IsAtEnd(); ++it1)
	{
		double b = stimateParameterB(it1);

		for (it2.GoToBegin(); !it2.IsAtEnd(); ++it2)
		{
			laplacianWeigh(it1, it2);
		}

	}

}

#endif /* GRAPH_H_ */
