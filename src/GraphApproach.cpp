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

int main() {

	type::grayImagePointer image = io::readImage<type::grayImage>(
			"/home/oscar/MEGA/post-doc/src/input/laudos/a3.png");

	Graph graph;
	graph.setImage(image);
	graph.build();

	std::cout << "DONE" << std::endl; // prints !!!Hello World!!!
	return 0;
}
