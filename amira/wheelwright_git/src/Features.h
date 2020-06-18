/** @file Features.h
 * Features for Projection.
 *
 * @author Matteo Riondato
 * @date 2017 10 20
 *
 * @copyright
 *
 * Copyright 2017 Matteo Riondato <riondato@acm.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef WHEELWRIGHT_FEATURES_H_
#define WHEELWRIGHT_FEATURES_H_

namespace Wheelwright {

/** The features of a Projection (and Vector).
 *
 * See documentation for Projection for an explaination of the different
 * features.
 */
enum Features : int {
    /** See documentation for Projection. */
    DISCRETE = 0b001,
    /** See documentation for Projection. */
    ONEDRAW = 0b0010,
    /** See documentation for Projection. */
    VANILLA = 0b0100,
    /** See documentation for Projection. */
    ZERO = 0b000
};

} // namespace Wheelwright

#endif // WHEELWRIGHT_FEATURES_H_
