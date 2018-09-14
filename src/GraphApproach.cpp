//============================================================================
// Name        : GraphApproach.cpp
// Author      : Oscar Cuadros Linares
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "itkImage.h"
#include "Graph.h"
#include "common/InputOutput.h"
#include "common/ImageProcessing.h"
#include "GraphClustering.h"

int main()
{

	std::cout << "START" << std::endl;

	using componentType = unsigned char;
	using rgbPixelType = itk::RGBPixel< componentType >;
	using rgbImageType = itk::Image< rgbPixelType, 2 >;

	type::rgbImagePointer image = io::readImage<rgbImageType>("/home/oscar/MEGA/post-doc/src/input/rp/patient_2/small_3.jpg");
	//type::rgbImagePointer image = io::readImage<rgbImageType>("/home/oscar/MEGA/post-doc/src/input/laudos/d3.png");
	Graph graph;
	graph.setImage(image);
	graph.setRadius(10);
	graph.build();

	GraphClustering graphClustering;

	graphClustering.fastGreedy(graph.getGraph(), graph.getWeights());
	//graphClustering.labelPropagation(graph.getGraph(), graph.getWeights());
	//graphClustering.multilevel(graph.getGraph(), graph.getWeights());
	//graphClustering.eigenvector(graph.getGraph(), graph.getWeights());

	unsigned width = image->GetLargestPossibleRegion().GetSize()[0];
	unsigned height = image->GetLargestPossibleRegion().GetSize()[1];

	io::writeImage<type::rgbImage>(graphClustering.membershipToImage(width, height, graph.getVertexIdMap()), "/home/oscar/MEGA/post-doc/src/output/clusters.png");

	std::cout << "DONE" << std::endl; // prints !!!Hello World!!!

	return 0;
}
