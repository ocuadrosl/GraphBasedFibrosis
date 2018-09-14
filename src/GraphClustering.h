/*
 * GraphClustering.h
 *
 *  Created on: Sep 12, 2018
 *      Author: oscar
 */

#ifndef SRC_GRAPHCLUSTERING_H_
#define SRC_GRAPHCLUSTERING_H_

#include "common/InputOutput.h"
#include "itkImage.h"

extern "C"
{
#include <igraph.h>
}

class GraphClustering
{
private:

	using componentType = unsigned char;
	using rgbPixelType = itk::RGBPixel< componentType >;
	using rgbImageType = itk::Image< rgbPixelType, 2 >;

public:
	GraphClustering();
	~GraphClustering()
	{
	}
	;

	void fastGreedy(const igraph_t& graph, const igraph_vector_t& weights);
	void labelPropagation(const igraph_t& graph, const igraph_vector_t& weights);
	void multilevel(const igraph_t& graph, const igraph_vector_t& weights);
	void eigenvector(const igraph_t& graph, const igraph_vector_t& weights);

	rgbImageType::Pointer membershipToImage(unsigned imageWigth, unsigned imageHeight, const std::vector<int>& vertexIdMap);

private:

	igraph_vector_t membership;

	std::vector<rgbPixelType> clusterColorGenerator(unsigned clusterNumber);

};

GraphClustering::GraphClustering()
{

	igraph_vector_init(&this->membership, 0);

}

void GraphClustering::multilevel(const igraph_t& graph, const igraph_vector_t& weights)
{

	igraph_community_multilevel(&graph, &weights, &this->membership, 0, 0);
	io::print("Multilevel", true);

}

std::vector<GraphClustering::rgbPixelType> GraphClustering::clusterColorGenerator(unsigned clusterNumber)
{

	std::vector<rgbPixelType> clusterColors(clusterNumber);
	std::srand(time(nullptr));
	for (unsigned i = 0; i < clusterNumber; ++i)
	{
		clusterColors[i][0] = std::rand() % 255;
		clusterColors[i][1] = std::rand() % 255;
		clusterColors[i][2] = std::rand() % 255;
	}

	return clusterColors;

}

GraphClustering::rgbImageType::Pointer GraphClustering::membershipToImage(unsigned imageWigth, unsigned imageHeight, const std::vector<int>& vertexIdMap)
{

	unsigned clusterNumber = igraph_vector_max(&this->membership);

	std::vector<rgbPixelType> clusterColors = clusterColorGenerator(clusterNumber);
	unsigned clusterIndex;
	std::vector<unsigned> matrixIndexTmp;
	rgbImageType::IndexType index;

	//creating a new image
	rgbImageType::Pointer newImage = rgbImageType::New();
	rgbImageType::RegionType region;
	rgbImageType::SizeType size;
	size[0] = imageWigth;
	size[1] = imageHeight;
	region.SetSize(size);
	newImage->SetRegions(region);
	newImage->Allocate();
	rgbPixelType black;
	black.Fill(0);
	newImage->FillBuffer(black);

	std::vector<int>::const_iterator vertexId;
	for (unsigned i = 0; i < igraph_vector_size(&this->membership); ++i)
	{
		clusterIndex = igraph_vector_e(&this->membership, i);

		vertexId = std::find(vertexIdMap.begin(), vertexIdMap.end(), static_cast<int>(i));
		matrixIndexTmp = math::vectorIndexToMatrixIndex(imageWigth, distance(vertexIdMap.begin(), vertexId));
		index[0] = matrixIndexTmp[0];
		index[1] = matrixIndexTmp[1];
		newImage->SetPixel(index, clusterColors[clusterIndex]);

	}

	return newImage;

}

void GraphClustering::fastGreedy(const igraph_t& graph, const igraph_vector_t& weights)
{

	igraph_vector_t modularity;
	igraph_matrix_t merges;
	igraph_vector_init(&modularity, 0);
	igraph_matrix_init(&merges, 0, 0);

	igraph_community_fastgreedy(&graph, &weights, &merges, &modularity, &this->membership);

	igraph_vector_destroy(&modularity);
	igraph_matrix_destroy(&merges);

	io::print("Fast greedy", true);

}

void GraphClustering::labelPropagation(const igraph_t& graph, const igraph_vector_t& weights)
{

	igraph_community_label_propagation(&graph, &this->membership, &weights, 0, 0, 0);

	io::print("Label propagation", true);

}

void GraphClustering::eigenvector(const igraph_t& graph, const igraph_vector_t& weights)
{

	igraph_vector_t modularity;
	igraph_matrix_t merges;
	igraph_vector_init(&modularity, 0);
	igraph_matrix_init(&merges, 0, 0);

	igraph_arpack_options_t options;
	igraph_arpack_options_init(&options);

	igraph_community_leading_eigenvector(&graph, &weights, &merges, &this->membership, igraph_vcount(&graph), &options, 0, false, 0, 0, 0, 0, 0);

}

#endif /* SRC_GRAPHCLUSTERING_H_ */
