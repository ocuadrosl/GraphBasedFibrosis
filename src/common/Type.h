/*
 * Type.h
 *
 *  Created on: Sep 3, 2018
 *      Author: oscar
 */

#ifndef COMMON_TYPE_H_
#define COMMON_TYPE_H_

#include "itkImage.h"
#include <itkRGBPixel.h>
namespace type
{

typedef double gray;
typedef itk::Image<double, 2> grayImage;
typedef grayImage::Pointer grayImagePointer;

typedef double componentType;
typedef itk::RGBPixel<double> rgbPixel;
typedef itk::Image<rgbPixel, 2> rgbImage;
typedef rgbImage::Pointer rgbImagePointer;

}

#endif /* COMMON_TYPE_H_ */
