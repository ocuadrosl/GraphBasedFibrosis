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

	//type::rgbImagePointer image = io::readImage<rgbImageType>("/home/oscar/MEGA/post-doc/src/input/rp/patient_1/00529 (3).jpg");
	type::rgbImagePointer image = io::readImage<rgbImageType>("/home/oscar/MEGA/post-doc/src/input/laudos/e3.png");

	//io::writeImage<type::grayImage>(ip::otsuThreshold<type::grayImage>(image), "/home/oscar/MEGA/post-doc/src/input/laudos/outsu.png");

	Graph graph;
	graph.setImage(image);
	graph.setRadius(3);
	graph.build();

	GraphClustering graphClustering;

	graphClustering.fastGreedy(graph.get(), graph.getWeights());

	std::vector<unsigned> imageSize(2);
	imageSize[0] = image->GetLargestPossibleRegion().GetSize()[0];
	imageSize[1] = image->GetLargestPossibleRegion().GetSize()[1];

	io::writeImage<type::rgbImage>(graphClustering.getMembershipAsImage(imageSize), "/home/oscar/MEGA/post-doc/src/input/laudos/clusters.png");

	std::cout << "DONE" << std::endl; // prints !!!Hello World!!!

	return 0;
}
