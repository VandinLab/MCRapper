/** @file NumericTools.h
 * Imported file from NetworKit, for the logically_equal function.
 *
 * This file is a combination and modification of two files taken from the
 * NetworKit project on 20170417. Their location in the NetworKit tree was
 * networkit/cpp/auxiliary/NumericTools.{h,cpp}
 *
 * @date 2017 06 17
 *
 * @copyright
 *
 * Copyright (c) 2013 Christian Staudt, Henning Meyerhenke
 *
 * MIT License (http://opensource.org/licenses/MIT)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial
 * portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 * NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
 * OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */
/*
 * NumericTools.h
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef MATTEOUTILS_NUMERICTOOLS_H_
#define MATTEOUTILS_NUMERICTOOLS_H_

#include <algorithm>
#include <cmath>
#include <limits>

namespace MatteoUtils {

/** Tests numeric values for equality within the smallest possible error.
*
* @param a,b           the values to compare.
* @param error_factor  the error factor, used in the approximate equality test.
*
* @returns  true if @a a == @a b or if |a-b| < |min(a,b)| * error_factor * eps,
*           where eps is the minimum positive representable double.  false
*           otherwise.
*/
inline bool logically_equal(const double a, const double b,
        const double error_factor=1.0) {
    return a == b ||
        std::abs(a - b) < std::abs(std::min(a, b)) *
        std::numeric_limits<double>::epsilon() * error_factor;
}

} // namespace MatteoUtils

#endif // MATTEOUTILS_NUMERICTOOLS_H_
