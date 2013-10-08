#ifndef REORDER_HPP__
#define REORDER_HPP__
#include <iterator>
#include <algorithm>
namespace Utility
{
  //! Reorders a vector according to a container of indices
  /*!
    It implements an algorithm that reorders a container in place.

    @param order_begin. Iterator to the begin of container of indexes.
    @param order_end. Iterator to the begin of container of indexes.
    @param v Iterator to the begin of container of values to be reordered.
  */
  template< typename order_iterator, typename value_iterator >
  void reorder ( order_iterator order_begin,
                 order_iterator order_end,
                 value_iterator v )
  {
    typedef typename std::iterator_traits< value_iterator >::value_type value_t;
    typedef typename std::iterator_traits< order_iterator >::value_type index_t;
    typedef typename std::iterator_traits< order_iterator >::difference_type diff_t;

    diff_t remaining = order_end - 1 - order_begin;
    for ( index_t s = index_t(), d; remaining > 0; ++ s )
    {
      for ( d = order_begin[s]; d > s; d = order_begin[d] ) ;
      if ( d == s )
      {
        -- remaining;
        value_t temp = v[s];
        while ( d = order_begin[d], d != s )
        {
          std::swap ( temp, v[d] );
          -- remaining;
        }
        v[s] = temp;
      }
    }
  }// end reorder
}// end namespace Utility
#endif
