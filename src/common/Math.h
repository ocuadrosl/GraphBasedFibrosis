/*
 * Math.h
 *
 *  Created on: Sep 3, 2018
 *      Author: oscar
 */

#ifndef COMMON_MATH_H_
#define COMMON_MATH_H_

namespace math
{
template<typename type> //vector type
inline double squaredEuclideanDistance(type const& p1, type const& p2, unsigned size = 2)
{
	double squaredSum = 0.0;
	for (unsigned i = 0; i < size; ++i)
	{
		squaredSum += std::pow(static_cast<double>(p2[i]) - static_cast<double>(p1[i]), 2.0);
	}

	return squaredSum;
}

template<typename type>
inline double euclideanDistance(type const& p1, type const& p2, unsigned size = 2)
{
	return std::sqrt(squaredEuclideanDistance<type>(p1, p2));

}


}

#endif /* COMMON_MATH_H_ */
