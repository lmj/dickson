/////////////////////////////////////////////////////////////////////////
// 
// Copyright (C) 2006 James M. Lawrence.  All rights reserved.
// 
// software: dickson-1.0.3
// file: test/common.hxx
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 
// 2. The origin of this software must not be misrepresented; you must
//    not claim that you wrote the original software.  If you use this
//    software in a product, an acknowledgment in the product 
//    documentation would be appreciated but is not required.
// 
// 3. Altered source versions must be plainly marked as such, and must
//    not be misrepresented as being the original software.
// 
// 4. The name of the author may not be used to endorse or promote
//    products derived from this software without specific prior written
//    permission.
// 
// THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY,
// AND WITH NO CLAIM AS TO ITS SUITABILITY FOR ANY PURPOSE.  IN NO EVENT
// SHALL THE AUTHOR BE LIABLE FOR ANY CLAIM WHATSOEVER MADE IN CONNECTION
// TO THIS SOFTWARE.
// 
/////////////////////////////////////////////////////////////////////////

#ifndef test_common_hxx
#define test_common_hxx

#include "dickson.hxx"
#include <cstdlib>
#include <limits>
#include <cmath>

/////////////////////////////////////////////////////////
// within
/////////////////////////////////////////////////////////

template< typename scalar >
bool within( scalar x, scalar y, scalar epsilon )
{
    return std::abs(x - y) <= epsilon ;
}

template< typename subalg >
bool within( const dickson::algebra<subalg> & x,
             const dickson::algebra<subalg> & y,
             typename dickson::traits<subalg>::scalar epsilon )
{
    for( unsigned int i = 0 ;
         i != dickson::algebra<subalg>::dimension ;
         ++i )
    {
        if( !within(x[i], y[i], epsilon) )
        {
            return false ;
        }
    }
    
    return true ;
}

template< typename subalg >
bool within( typename dickson::traits<subalg>::scalar x,
             const dickson::algebra<subalg> & y,
             typename dickson::traits<subalg>::scalar epsilon )
{
    return within(dickson::algebra<subalg>(x), y, epsilon) ;
}

template< typename subalg >
bool within( const dickson::algebra<subalg> & x,
             typename dickson::traits<subalg>::scalar y,
             typename dickson::traits<subalg>::scalar epsilon )
{
    return within(y, x, epsilon) ;
}

/////////////////////////////////////////////////////////
// rand
/////////////////////////////////////////////////////////

template< typename scalar, bool is_integer >
struct select_rand ;

template< typename scalar >
struct select_rand<scalar, false>
{
    static scalar rand()
    {
        return scalar(2)*(scalar(std::rand())/RAND_MAX) - scalar(1) ;
    }
} ;

template< typename scalar >
struct select_rand<scalar, true>
{
    static scalar rand()
    {
        // small range to prevent overflow.
        return scalar(double(100)
                      *
                      (double(2)
                       *
                       (double(std::rand())/RAND_MAX) - double(1))) ;
    }
} ;

template< typename scalar >
struct generate_rand
{
    static scalar rand()
    {
        return
            select_rand<scalar, std::numeric_limits<scalar>::is_integer>::
            rand() ;
    }
} ;

template< typename subalg >
struct generate_rand< dickson::algebra<subalg> >
{
    static dickson::algebra<subalg> rand()
    {
        typedef typename dickson::traits<subalg>::scalar scalar ;
        dickson::algebra<subalg> res ;
        
        for( unsigned int i = 0 ;
             i != dickson::algebra<subalg>::dimension ;
             ++i )
        {
            res[i] =
                select_rand<scalar, std::numeric_limits<scalar>::is_integer>::
                rand() ;
        }
        
        return res ;
    }
} ;

template< typename algebra >
algebra rand()
{
    return generate_rand<algebra>::rand() ;
}

#endif
