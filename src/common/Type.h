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

using gray = unsigned char;
using grayImage = itk::Image<gray,2>;
typedef grayImage::Pointer grayImagePointer;

typedef unsigned char componentType;
typedef itk::RGBPixel<componentType> rgbPixel;
typedef itk::Image<rgbPixel, 2> rgbImage;
typedef rgbImage::Pointer rgbImagePointer;

}

#endif /* COMMON_TYPE_H_ */
