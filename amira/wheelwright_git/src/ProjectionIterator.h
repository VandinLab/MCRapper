/** @file ProjectionIterator.h
 * Iterator for the vectors in a Projection object.
 *
 * @author Matteo Riondato
 * @date 2017 06 19
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

#ifndef WHEELWRIGHT_PROJECTIONITERATOR_H
#define WHEELWRIGHT_PROJECTIONITERATOR_H

#include <iterator>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>

namespace Wheelwright::impl {

/** Iterator through the vectors of a Projection.
 *
 * @tparam Iterator  the type of the Iterator through the vectors. It is
 *                   expected to behave like a std::map (or std::unordered_map)
 *                   iterator, i.e., to return pairs when dereferenced.
 */
template<typename Iterator> class ProjectionIterator :
        public boost::iterator_facade<ProjectionIterator<Iterator>,
            typename
            std::iterator_traits<Iterator>::value_type::first_type::element_type,
            boost::bidirectional_traversal_tag>
    {
    // TODO: ProjectionIterator: Enable this class only if Iterator behaves
    // like a std::map or std::unordered_map iterator where the keys are
    // shared_ptr.
    private:
        struct enabler {};
        friend class boost::iterator_core_access;
        template<class OtherIterator> bool equal(
                ProjectionIterator<OtherIterator> const& o) const {
            return it == o.it;
        }
        void decrement() { --it; }
        void increment() { ++it; }
        typename
            std::iterator_traits<Iterator>::value_type::first_type::element_type & dereference()
            const { return *(it->first); }
        /** The iterator through the vectors of the Projection. */
        Iterator it;
    public:
        // TODO: Wheelwright::ProjectionIterator::ProjectionIterator: In an
        // ideal world, this constructor would be private, and we would have
        // Projection as a friend class. But Projection has many template
        // arguments that we do not care about here, so we make this function a
        // public function. The downside is that potentially classes different
        // than Projection may create objects of this type (and maybe that's not
        // an issue, because all this class really does is returning Value's
        // from an iterator over std::pair<std::shared_ptr<Value>, foo>).
        /** Constructor. */
        ProjectionIterator<Iterator>(Iterator it) : it {it} {}

        /** Copy constructor. */
        template <typename OtherIterator> ProjectionIterator(
                ProjectionIterator<OtherIterator> const& o,
                typename boost::enable_if<
                    boost::is_convertible<typename
                    std::iterator_traits<OtherIterator>::value_type::first_type::element_type *,
                    typename
            std::iterator_traits<Iterator>::value_type::first_type::element_type *>,
                    enabler>::type = enabler()) : it {o.it} {}

}; // class ProjectionIterator

} // namespace Wheelwright::impl

#endif // WHEELWRIGHT_PROJECTIONITERATOR_H
