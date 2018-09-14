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
#include "itkImageRegionIteratorWithIndex.h"

extern "C"
{
#include <igraph.h>
}

class Graph
{

private:
	//standard class typedefs
	using componentType = unsigned char;
	using rgbPixelType = itk::RGBPixel< componentType >;
	using rgbImageType = itk::Image< rgbPixelType, 2 >;

	using grayType = double;
	using grayImageType = itk::Image<double,2>;
	using grayIteratorType = itk::ImageRegionIterator<grayImageType>;
	using grayIteratorIndexType = itk::ImageRegionIteratorWithIndex<grayImageType>;
	using grayImagePointer = typename grayImageType::Pointer;
	using grayIndexType = grayImageType::IndexType;
	const int BACKGROUND = -1; //background value

public:
	Graph();
	~Graph()
	{
	}
	;

	void setImage(rgbImageType::Pointer rgbImage);
	void buildOld();
	void build();

	igraph_t getGraph() const;
	igraph_vector_t getWeights() const;
	std::vector<int> getVertexIdMap() const;

	void setRadius(unsigned radius);

private:

	grayImagePointer image;

	double laplaceWeigth;
	unsigned radius;
	igraph_t graph;
	igraph_vector_t weights;
	std::vector<int> vertexIdMap;

	//methods
	inline double laplacianWeight(grayIteratorIndexType imageIt, grayIteratorIndexType neighborIt, double b);
	inline double stimateParameterB(grayIteratorIndexType imageIt, grayIteratorIndexType neighborhoodIt);
	inline std::vector<double> computeMeans(grayImageType::RegionType region);
	void segmentBackground();
	void imageGradient();
	void createVertexIdMap();

};

Graph::Graph()
{
	laplaceWeigth = 1.0;
	radius = 1;

}

void Graph::setRadius(unsigned radius)
{
	this->radius = radius;

}
std::vector<int> Graph::getVertexIdMap() const
{
	return this->vertexIdMap;
}
/*
 * compute the median for index and intensity
 *
 * */
inline std::vector<double> Graph::computeMeans(grayImageType::RegionType region)
{
	//[3] intensity mean
	std::vector<double> means(3, 0.0);

	grayIteratorType it(image, region);
	grayIndexType index;

	double count = 0;
	for (it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		index = it.GetIndex();
		means[0] += index[0];
		means[1] += index[1];
		means[3] += it.Get();
		++count;

	}

	means[0] = means[0] / count;
	means[1] = means[1] / count;
	means[2] = means[2] / count;

	return means;

}

/* returns the igraph_t object
 *
 * */
igraph_t Graph::getGraph() const
{

	return this->graph;

}
igraph_vector_t Graph::getWeights() const
{

	return this->weights;

}

void Graph::imageGradient()
{

	grayImagePointer gradient = ip::imageGradient<grayImageType>(image);

	typedef itk::CastImageFilter<grayImageType, rgbImageType> CastFilterType;
	CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(gradient);
	castFilter->Update();
	io::writeImage<rgbImageType>(castFilter->GetOutput(), "/home/oscar/MEGA/post-doc/src/output/gradient.png");

	grayIteratorType imIt(image, image->GetLargestPossibleRegion()); //image iterator
	grayIteratorType segIt(gradient, image->GetLargestPossibleRegion()); //segmented image iterator

	for (imIt.GoToBegin(), segIt.GoToBegin(); !imIt.IsAtEnd(); ++imIt, ++segIt)
	{

		if (segIt.Get() <= 0)
		{
			imIt.Set(BACKGROUND);
		}

	}

}

/*
 * Separate background and foreground, background is set to 0.
 *
 * */

void Graph::segmentBackground()
{

	grayImagePointer segmented = ip::histogramThreshold<grayImageType>(image, "triangle", 100);
	//grayImagePointer segmented = ip::simpleThreshold<grayImageType>(image, 0, 240);

	typedef itk::CastImageFilter<grayImageType, rgbImageType> CastFilterType;
	CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(segmented);
	castFilter->Update();
	io::writeImage<rgbImageType>(castFilter->GetOutput(), "/home/oscar/MEGA/post-doc/src/output/threshold.png");

	grayIteratorType imIt(image, image->GetLargestPossibleRegion()); //image iterator
	grayIteratorType segIt(segmented, image->GetLargestPossibleRegion()); //segmented image iterator

	for (imIt.GoToBegin(), segIt.GoToBegin(); !imIt.IsAtEnd(); ++imIt, ++segIt)
	{

		if (segIt.Get() == 0) //outside default value = 0
		{
			imIt.Set(BACKGROUND);
		}

	}

	castFilter->SetInput(image);
	castFilter->Update();
	io::writeImage<rgbImageType>(castFilter->GetOutput(), "/home/oscar/MEGA/post-doc/src/output/segmented.png");

}

void Graph::setImage(rgbImageType::Pointer rgbImage)
{

	using luminanceImageFilterType = itk::RGBToLuminanceImageFilter< rgbImageType, grayImageType >;
	luminanceImageFilterType::Pointer luminanceImageFilter = luminanceImageFilterType::New();
	luminanceImageFilter->SetInput(rgbImage);

	luminanceImageFilter->Update();
	this->image = luminanceImageFilter->GetOutput();

	segmentBackground();
	//imageGradient();
	createVertexIdMap();

}

inline double Graph::stimateParameterB(grayIteratorIndexType imageIt, grayIteratorIndexType neighborIt)
{

	double count = 0;

	//compute mean
	grayIndexType meanIndex;
	meanIndex.Fill(0);
	double meanIntensity = 0;

	for (neighborIt.GoToBegin(); !neighborIt.IsAtEnd(); ++neighborIt)
	{
		meanIntensity += neighborIt.Get();
		meanIndex[0] += neighborIt.GetIndex()[0];
		meanIndex[1] += neighborIt.GetIndex()[1];
		++count;
	}

	meanIntensity /= count;
	meanIndex[0] /= count;
	meanIndex[1] /= count;

	double sumIntensity = 0.0;
	double sumIndex = 0.0;

	for (neighborIt.GoToBegin(); !neighborIt.IsAtEnd(); ++neighborIt)
	{
		sumIntensity += std::abs(meanIntensity - imageIt.Get());
		sumIndex += math::euclideanDistance<grayIndexType>(meanIndex, neighborIt.GetIndex());

	}

	double b = (sumIndex + sumIntensity) / count;

	return (b == 0) ? 0.0000000000001 : b; //to avoid division by zero

}

inline double Graph::laplacianWeight(grayIteratorIndexType imageIt, grayIteratorIndexType neighborIt, double b)
{

	double indexTerm = 0;
	double intensityTerm = 0;

	indexTerm = math::euclideanDistance<grayIndexType>(neighborIt.GetIndex(), imageIt.GetIndex());
	intensityTerm = std::abs(neighborIt.Get() - imageIt.Get());

	double numerator = std::abs(intensityTerm + indexTerm);

	return (1.0 / (2.0 * b)) * std::exp(-(laplaceWeigth * numerator) / b); //2.0

}

/*
 * Since not all pixels are mapped as graph vertices, a vector with its vector
 * index is created to preserve their original index after building the graph
 * */
void Graph::createVertexIdMap()
{

	grayImageType::SizeType imageSize = image->GetLargestPossibleRegion().GetSize();

	vertexIdMap = std::vector<int>(imageSize[1] * imageSize[0], 0);

	grayIteratorType it(image, image->GetLargestPossibleRegion());

	unsigned vertexCount = 0;
	unsigned vertexId = 0;
	int x;
	for (it.GoToBegin(); !it.IsAtEnd(); ++it, ++vertexCount)
	{
		if (it.Get() != BACKGROUND)
		{
			this->vertexIdMap[vertexCount] = vertexId;
			++vertexId;

		}
		else
		{
			this->vertexIdMap[vertexCount] = BACKGROUND;

		}

	}

}

void Graph::build()
{

//creating iterators
	grayIteratorIndexType imageIt(image, image->GetLargestPossibleRegion());
	grayIteratorIndexType neighborIt(image, image->GetLargestPossibleRegion());

//auxiliary declarations
	grayImageType::IndexType lowerIndex;
	grayImageType::IndexType upperIndex;
	grayImageType::IndexType currentIndex;
	grayImageType::IndexType neighborIndex;
	grayImageType::SizeType imageSize = image->GetLargestPossibleRegion().GetSize();
	grayImageType::RegionType neighborhood;
	unsigned imageVecIndex = 0;
	unsigned neigborVecIndex = 0;

//Laplace variables
	double b;
	double weight;

//igraph variables
	igraph_vector_t edges;
	igraph_vector_init(&edges, 0);
	igraph_vector_init(&weights, 0);

//reserving some memory to avoid expensive resize operations O(n)
	igraph_vector_reserve(&edges, imageSize[1] * imageSize[0] * radius);
	igraph_vector_reserve(&weights, imageSize[1] * imageSize[0] * radius);

//creating an auxiliary matrix to verify if a new edge already exits
	std::vector<bool> heigthAux(imageSize[1] * imageSize[0], false);
	std::vector<std::vector<bool> > edgeAux(imageSize[1] * imageSize[0], heigthAux);

	for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt, ++imageVecIndex)
	{

		if (vertexIdMap[imageVecIndex] == BACKGROUND) //background
		{
			continue;
		}

		currentIndex = imageIt.GetIndex();

		//width
		lowerIndex[0] = (static_cast<signed>(currentIndex[0] - radius) < 0) ? 0 : currentIndex[0] - radius;
		upperIndex[0] = (currentIndex[0] + radius >= imageSize[0]) ? imageSize[0] - 1 : currentIndex[0] + radius;

		//height
		lowerIndex[1] = (static_cast<signed>(currentIndex[1] - radius) < 0) ? 0 : currentIndex[1] - radius;
		upperIndex[1] = (currentIndex[1] + radius >= imageSize[1]) ? imageSize[1] - 1 : currentIndex[1] + radius;

		neighborhood.SetIndex(lowerIndex);
		neighborhood.SetUpperIndex(upperIndex);

		grayIteratorIndexType neighborIt(image, neighborhood);
		b = stimateParameterB(imageIt, neighborIt);

		for (neighborIt.GoToBegin(); !neighborIt.IsAtEnd(); ++neighborIt)
		{
			neighborIndex = neighborIt.GetIndex();
			neigborVecIndex = (neighborIndex[1] * imageSize[0]) + neighborIndex[0];

			if (vertexIdMap[neigborVecIndex] == BACKGROUND) //background
			{
				continue;
			}

			if (edgeAux[imageVecIndex][neigborVecIndex] == false && edgeAux[neigborVecIndex][imageVecIndex] == false && imageVecIndex != neigborVecIndex)
			{

				weight = laplacianWeight(imageIt, neighborIt, b);

				igraph_vector_resize(&edges, igraph_vector_size(&edges) + 2);
				igraph_vector_set(&edges, igraph_vector_size(&edges) - 1, vertexIdMap[imageVecIndex]);
				igraph_vector_set(&edges, igraph_vector_size(&edges) - 2, vertexIdMap[neigborVecIndex]);

				igraph_vector_resize(&weights, igraph_vector_size(&weights) + 1);
				igraph_vector_set(&weights, igraph_vector_size(&weights) - 1, weight);
				edgeAux[imageVecIndex][neigborVecIndex] = true;
				edgeAux[neigborVecIndex][imageVecIndex] = true;
			}

		}

	}

	igraph_create(&graph, &edges, 0, 0);

//igraph_simplify(&graph, true, true, 0);

	igraph_vector_destroy(&edges);

	std::cout << igraph_vector_size(&weights) << std::endl;

	io::print("Building graph", true);

}

#endif /* GRAPH_H_ */
