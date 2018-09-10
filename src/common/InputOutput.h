/*
 * InputOutput.h
 *
 *  Created on: Sep 3, 2018
 *      Author: oscar
 */

#ifndef COMMON_INPUTOUTPUT_H_
#define COMMON_INPUTOUTPUT_H_

#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include <iomanip>
#include "Type.h"

//#include "QuickView.h"

namespace io
{

void print(const std::string& message, const bool & ok, short size = 50)
{

	short diff = size - message.length();

	if (ok)
	{
		std::cout << message << std::setfill(' ') << std::setw(diff) << "[OK]" << std::endl;
	}
	else
	{
		//std::cout << message << std::setfill(' ') << std::setw(diff) << "[FAIL]" << std::endl;
	}

}

template<typename imageType>
typename imageType::Pointer readImage(const std::string& file_name)
{

	typedef itk::ImageFileReader<imageType> readerType;

	typename readerType::Pointer reader = readerType::New();

	reader->SetFileName(file_name);

	try
	{
		reader->Update();

	} catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception reading: " << file_name << " null pointer is returned" << std::endl;
		typename imageType::Pointer nullImage;
		return nullImage;
	}



	return reader->GetOutput();

}

template<typename imageType>
void writeImage(const typename imageType::Pointer &image, const std::string &fileName, const std::string& message = "")
{

	typedef itk::ImageFileWriter<imageType> writer_t;
	typename writer_t::Pointer writer = writer_t::New();
	writer->SetFileName(fileName);
	writer->SetInput(image);
	writer->Update();

	print("Writing image " + message, true);
}


template<typename imageType>
void quickView(typename imageType::Pointer image)
{
	//QuickView viewer;

	//viewer.AddRGBImage(image);
	//viewer.AddImage(rescaleFilter->GetOutput());
	//viewer.Visualize();
}

}

#endif /* COMMON_INPUTOUTPUT_H_ */
