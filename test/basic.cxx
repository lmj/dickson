/////////////////////////////////////////////////////////////////////////
// 
// Copyright (C) 2006 James M. Lawrence.  All rights reserved.
// 
// software: dickson-1.0.3
// file: test/basic.cxx
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
//
// A complete or near-complete test of dickson.hxx
//

#include "common.hxx"
#include <iostream>
#include <sstream>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <ctime>

typedef long double real ;

const real EPSILON = 1e-8 ;

using namespace dickson::scalar_func ;

namespace {

int g_assertion_count = 0 ;

} // anon namespace

/////////////////////////////////////////////////////////
// assertions
////////////////////////////////////////////////////////

#define show(a) { std::cout << (a) << " <== " << #a << std::endl ; }

#define assert_close(a, b) \
    if( !within((a), (b), EPSILON) ) \
    { \
        std::cout \
            << "Assertion failed. Not within " \
            << EPSILON \
            << ":" \
            << std::endl ; \
        show(a) ; \
        show(b) ; \
        std::cout << __FILE__ << ":" << __LINE__ << std::endl ; \
        exit(1) ; \
    } \
    else \
    { \
        g_assertion_count += 1 ; \
    } \

#define assert_equal(a, b) \
    if( !((a) == (b)) ) \
    { \
        std::cout \
            << "Assertion failed. Not equal:" \
            << std::endl ; \
        show(a) ; \
        show(b) ; \
        std::cout << __FILE__ << ":" << __LINE__ << std::endl ; \
        exit(1) ; \
    } \
    else \
    { \
        g_assertion_count += 1 ; \
    } \

/////////////////////////////////////////////////////////
// factorial (approx)
/////////////////////////////////////////////////////////

template< typename scalar >
scalar factorial( unsigned int x )
{
    scalar result = scalar(1) ;
    
    for( unsigned int t = 2 ; t <= x ; ++t )
    {
        result *= scalar(t) ;
    }
    
    return result ;
}

/////////////////////////////////////////////////////////
// core
/////////////////////////////////////////////////////////

template< typename algebra >
void core()
{
    typedef typename algebra::scalar scalar ;
    typedef typename algebra::subalg subalg ;
    const unsigned int dim = algebra::dimension ;

    algebra a = rand<algebra>() ;
    algebra b = rand<algebra>() ;
    scalar s = rand<scalar>() ;

    //// Construct (0, 0, 0, ..., 0).
    //algebra() ;
    {
        algebra a ;
        for( unsigned int i = 0 ; i < dim ; ++i )
        {
            assert_equal(a[i], scalar(0)) ;
        }
    }

    //// Construct (s, 0, 0, ..., 0).
    //explicit algebra( scalar s ) ;
    {
        algebra a = algebra(s) ;

        assert_equal(a[0], s) ;
        
        for( unsigned int i = 1 ; i != dim ; ++i )
        {
            assert_equal(a[i], scalar(0)) ;
        }
    }

    //// Construct (a + b i) where 'i' is the appropriate unit element
    //// orthogonal to a.
    //algebra( const subalg & a, const subalg & b ) ;

    //// The first subalgebra; corresponds to 'a' in algebra(a, b).
    //const subalg & first() const ;
    //      subalg & first() ;

    //// The second subalgebra; corresponds to 'b' in algebra(a, b).
    //const subalg & second() const ;
    //      subalg & second() ;

    {
        subalg a = rand<subalg>() ;
        subalg b = rand<subalg>() ;
        algebra c = algebra(a, b) ;

        assert_equal(a, static_cast<const algebra*>(&c)->first()) ;
        assert_equal(b, static_cast<const algebra*>(&c)->second()) ;
        assert_equal(a, static_cast<algebra*>(&c)->first()) ;
        assert_equal(b, static_cast<algebra*>(&c)->second()) ;
    }

    // given:
    //// Scalar coefficients: 0...(dimension - 1).
    //scalar   operator[]( unsigned int index ) const ;
    //scalar & operator[]( unsigned int index ) ;
    {
        for( unsigned int i = 0 ; i != dim ; ++i )
        {
            assert_equal(static_cast<algebra*>(&a)->operator[](i),
                         static_cast<const algebra*>(&a)->operator[](i)) ;
        }
    }

    //// (*this)[0]
    //scalar   realpart() const ;
    //scalar & realpart() ;
    {
        assert_equal(a[0], static_cast<const algebra*>(&a)->realpart()) ;
        assert_equal(a[0], static_cast<algebra*>(&a)->realpart()) ;
    }

    //// member operators
    //algebra & operator+=( const algebra & a ) ;
    //algebra & operator-=( const algebra & a ) ;
    //algebra & operator*=( const algebra & a ) ;
    //algebra & operator+=( scalar a ) ;
    //algebra & operator-=( scalar a ) ;
    //algebra & operator*=( scalar a ) ;
    //algebra & operator/=( scalar a ) ;
    {
        algebra c = a ;
        c += b ;
        assert_close(a + b, c) ;
    }
    {
        algebra c = a ;
        c -= b ;
        assert_close(a - b, c) ;
    }
    {
        algebra c = a ;
        c *= b ;
        assert_close(a*b, c) ;
    }
    {
        algebra c = a ;
        c += s ;
        assert_close(a + s, c) ;
    }
    {
        algebra c = a ;
        c -= s ;
        assert_close(a - s, c) ;
    }
    {
        algebra c = a ;
        c *= s ;
        assert_close(a*s, c) ;
    }
    {
        algebra c = a ;
        c /= s ;
        assert_close(a/s, c) ;
    }

    //template< typename subalg >
    //algebra<subalg>
    //operator+( const algebra<subalg> & a,
    //           const algebra<subalg> & b ) ;
    {
        algebra c ;
        for( unsigned int i = 0 ; i != dim ; ++i )
        {
            c[i] = a[i] + b[i] ;
        }
        assert_close(a + b, c) ;
    }
    
    //
    //template< typename subalg >
    //algebra<subalg>
    //operator+( typename traits<subalg>::scalar a,
    //           const algebra<subalg> & b ) ;
    assert_close(s + b, algebra(s) + b) ;
    
    //
    //template< typename subalg >
    //algebra<subalg>
    //operator+( const algebra<subalg> & a,
    //           typename traits<subalg>::scalar b ) ;
    assert_close(a + s, a + algebra(s)) ;

    //
    //template< typename subalg >
    //const algebra<subalg> &
    //operator+( const algebra<subalg> & a ) ;
    assert_equal(a, +a) ;

    //
    //template< typename subalg >
    //algebra<subalg>
    //operator-( const algebra<subalg> & a,
    //           const algebra<subalg> & b ) ;
    {
        algebra c ;
        for( unsigned int i = 0 ; i != dim ; ++i )
        {
            c[i] = a[i] - b[i] ;
        }
        assert_close(a - b, c) ;
    }
    
    //
    //template< typename subalg >
    //algebra<subalg>
    //operator-( typename traits<subalg>::scalar a,
    //           const algebra<subalg> & b ) ;
    assert_close(s - b, algebra(s) - b) ;
    
    //
    //template< typename subalg >
    //algebra<subalg>
    //operator-( const algebra<subalg> & a,
    //           typename traits<subalg>::scalar b ) ;
    assert_close(a - s, a - algebra(s)) ;

    //
    //template< typename subalg >
    //algebra<subalg>
    //operator-( const algebra<subalg> & a ) ;
    assert_close(-a, algebra(0) - a) ;

    //template< typename subalg >
    //algebra<subalg>
    //operator*( typename traits<subalg>::scalar a,
    //           const algebra<subalg> & b ) ;
    assert_close(s*a, algebra(s)*a) ;

    //
    //template< typename subalg >
    //algebra<subalg>
    //operator*( const algebra<subalg> & a,
    //           typename traits<subalg>::scalar b ) ;
    assert_close(a*s, a*algebra(s)) ;

    // given:
    //template< typename subalg >
    //algebra<subalg>
    //operator*( const algebra<subalg> & a,
    //           const algebra<subalg> & b ) ;

    //
    //template< typename subalg >
    //algebra<subalg>
    //operator/( const algebra<subalg> & a,
    //           typename traits<subalg>::scalar b ) ;
    assert_close(a/s, a*(scalar(1)/s)) ;
    
    //
    //template< typename subalg >
    //bool
    //operator==( const algebra<subalg> & a,
    //            const algebra<subalg> & b ) ;
    //
    //template< typename subalg >
    //bool
    //operator==( typename algebra<subalg>::scalar a,
    //            const algebra<subalg> & b ) ;
    //
    //template< typename subalg >
    //bool
    //operator==( const algebra<subalg> & a,
    //            typename algebra<subalg>::scalar b ) ;
    //
    //template< typename subalg >
    //bool
    //operator!=( const algebra<subalg> & a,
    //            const algebra<subalg> & b ) ;
    //
    //template< typename subalg >
    //bool
    //operator!=( typename algebra<subalg>::scalar a,
    //            const algebra<subalg> & b ) ;
    //
    //template< typename subalg >
    //bool
    //operator!=( const algebra<subalg> & a,
    //            typename algebra<subalg>::scalar b ) ;
    {
        algebra a = b ;
        assert_equal((a == b), true) ;
        assert_equal((a != b), false) ;
    }
    {
        algebra a = algebra(s) ;
        assert_equal((a == s), true) ;
        assert_equal((s == a), true) ;
        assert_equal((a != s), false) ;
        assert_equal((s != a), false) ;
    }
    
    //
    //// realpart(a) == a[0]
    //template< typename subalg >
    //typename traits<subalg>::scalar
    //realpart( const algebra<subalg> & a ) ;
    assert_equal(realpart(a), a[0]) ;
    
    //
    //// purepart(a) == a - realpart(a)
    //template< typename subalg >
    //algebra<subalg>
    //purepart( const algebra<subalg> & a ) ;
    assert_close(purepart(a), a - realpart(a)) ;
    
    //// conj(a) == a - 2*purepart(a)
    //template< typename subalg >
    //algebra<subalg>
    //conj( const algebra<subalg> & a ) ;
    assert_close(conj(a), a - scalar(2)*purepart(a)) ;
    
    //
    //// norm(a) == a*conj(a)
    //template< typename subalg >
    //typename traits<subalg>::scalar
    //norm( const algebra<subalg> & a ) ;
    assert_close(algebra(norm(a)), a*conj(a)) ;
    assert_close(norm(a), realpart(a*conj(a))) ;

    //
    //// abs(a) == sqrt(norm(a))
    //template< typename subalg >
    //typename traits<subalg>::scalar
    //abs( const algebra<subalg> & a ) ;
    assert_close(abs(a), sqrt(norm(a))) ;

    //
    //// a*inv(a) == 1
    //template< typename subalg >
    //algebra<subalg>
    //inv( const algebra<subalg> & a ) ;
    assert_close(a*inv(a), algebra(scalar(1))) ;
    
    //
    //// dot(a) == (b*conj(a) + a*conj(b))/2
    //
    //template< typename subalg >
    //typename traits<subalg>::scalar
    //dot( const algebra<subalg> & a,
    //     const algebra<subalg> & b ) ;
    assert_close(dot(a, b), (b*conj(a) + a*conj(b))/scalar(2)) ;

    //
    //template< typename subalg >
    //typename traits<subalg>::scalar
    //dot( typename traits<subalg>::scalar a,
    //     const algebra<subalg> & b ) ;
    assert_close(dot(a, s), dot(a, algebra(s))) ;

    //
    //template< typename subalg >
    //typename traits<subalg>::scalar
    //dot( const algebra<subalg> & a,
    //     typename traits<subalg>::scalar b ) ;
    assert_close(dot(s, b), dot(algebra(s), b)) ;

    //
    //// cross(a) == (b*conj(a) - a*conj(b))/2
    //
    //template< typename subalg >
    //algebra<subalg>
    //cross( const algebra<subalg> & a,
    //       const algebra<subalg> & b ) ;
    assert_close(cross(a, b), (b*conj(a) - a*conj(b))/scalar(2)) ;
    
    //
    //template< typename subalg >
    //algebra<subalg>
    //cross( typename traits<subalg>::scalar a,
    //       const algebra<subalg> & b ) ;
    assert_close(cross(s, b), cross(algebra(s), b)) ;
    
    //
    //template< typename subalg >
    //algebra<subalg>
    //cross( const algebra<subalg> & a,
    //       typename traits<subalg>::scalar b ) ;
    assert_close(cross(a, s), cross(a, algebra(s))) ;

    //// streams
    //
    //template< typename subalg, typename T_Char, typename T_Traits >
    //std::basic_ostream<T_Char, T_Traits> &
    //operator<<( std::basic_ostream<T_Char, T_Traits> & out,
    //            const algebra<subalg> & a ) ;
    //
    //template< typename subalg, typename T_Char, typename T_Traits >
    //std::basic_istream<T_Char, T_Traits> &
    //operator>>( std::basic_istream<T_Char, T_Traits> & in,
    //            algebra<subalg> & a ) ;
    {
        std::ostringstream out ;
        out.precision(std::numeric_limits<scalar>::digits10 + 1) ;
        out << a ;

        std::istringstream in(out.str()) ;
        algebra b ;
        in >> b ;
        assert_close(a, b) ;
    }
}

/////////////////////////////////////////////////////////
// series
/////////////////////////////////////////////////////////

template< typename algebra >
void series()
{
    using namespace dickson::scalar_func ;
    typedef typename dickson::traits<algebra>::scalar scalar ;

    // ensure is abs(a) less than pi
    algebra a = rand<algebra>() ;
    a /= abs(a) ;
    a *= scalar(3.14)*rand<scalar>() ;

    {
        assert_close(exp(log(a)), a) ;
        assert_close(log(exp(a)), a) ;
    }
    {
        algebra u = purepart(a) ;
        u /= abs(u) ;

        scalar r = rand<scalar>() ;
        
        assert_close(cos(a),  (exp(u*a) + exp(-u*a))/real(2)) ;
        assert_close(sin(a),  (exp(u*a) - exp(-u*a))*inv(u)/real(2)) ;
        assert_close(cosh(a), (exp(a) + exp(-a))/real(2)) ;
        assert_close(sinh(a), (exp(a) - exp(-a))/real(2)) ;
        
        assert_close(exp(u*r*a), cos(r*a) + u*sin(r*a)) ;
        assert_close(exp(r*a), cosh(r*a) + sinh(r*a)) ;
    }
    {
        algebra exp_a = algebra(scalar(0)) ;
        algebra cos_a = algebra(scalar(0)) ;
        algebra sin_a = algebra(scalar(0)) ;
        algebra cosh_a = algebra(scalar(0)) ;
        algebra sinh_a = algebra(scalar(0)) ;

        // log needs a smaller radius
        algebra b = rand<algebra>() ;
        b /= abs(b) ;
        b *= scalar(0.4)*rand<real>() ;
        algebra log_b_plus_1 = algebra(0) ;
        
        for( unsigned int n = 0 ; n <= 20 ; ++n )
        {
            const algebra t = pow(a, n)/factorial<scalar>(n) ;
            const algebra p = pow(b, n)/real(n) ;
            
            exp_a += t ;
            
            switch( n % 4 )
            {
            case 0:
                cos_a += t ;
                cosh_a += t ;
                if( n > 0 ) log_b_plus_1 -= p ;
                break ;
            case 1:
                sin_a += t ;
                sinh_a += t ;
                log_b_plus_1 += p ;
                break ;
            case 2:
                cos_a -= t ;
                cosh_a += t ;
                log_b_plus_1 -= p ;
                break ;
            case 3:
                sin_a -= t ;
                sinh_a += t ;
                log_b_plus_1 += p ;
                break ;
            }
        }
        
        assert_close(exp_a, exp(a)) ;
        assert_close(cos_a, cos(a)) ;
        assert_close(sin_a, sin(a)) ;
        assert_close(cosh_a, cosh(a)) ;
        assert_close(sinh_a, sinh(a)) ;
        assert_close(log_b_plus_1, log(b + real(1))) ;
    }    
    {
        assert_close(tan(a), sin(a)*inv(cos(a))) ;
        assert_close(tan(a), inv(cos(a))*sin(a)) ;
        assert_close(tanh(a), sinh(a)*inv(cosh(a))) ;
        assert_close(tanh(a), inv(cosh(a))*sinh(a)) ;
    }
    {
        const real r = rand<real>() ;
        assert_close(pow(a, r), exp(r*log(a))) ;
        assert_close(pow(abs(r), a), exp(a*log(abs(r)))) ;
    }
    {
        typedef unsigned int uint ;

        algebra an = algebra(scalar(1)) ;

        assert_close(pow(a, uint(0)), an) ;
        an *= a ;
        assert_close(pow(a, uint(1)), an) ;
        an *= a ;
        assert_close(pow(a, uint(2)), an) ;
        an *= a ;
        assert_close(pow(a, uint(3)), an) ;
        an *= a ;
        assert_close(pow(a, uint(4)), an) ;
        an *= a ;
        assert_close(pow(a, uint(5)), an) ;
    }
    {
        assert_close(a, sqrt(a)*sqrt(a)) ;
    }
}

/////////////////////////////////////////////////////////
// battery
/////////////////////////////////////////////////////////

template< typename algebra >
void battery()
{
    core<algebra>() ;
    series<algebra>() ;
}

/////////////////////////////////////////////////////////
// main
/////////////////////////////////////////////////////////

int main()
{
    std::srand(std::time(0)) ;
    std::cout.precision(std::numeric_limits<real>::digits10 + 1) ;
    
    for( int i = 0 ; i != 5 ; ++i )
    {
        battery<dickson::generate_algebra<real, 1>::type>() ;
        battery<dickson::generate_algebra<real, 2>::type>() ;
        battery<dickson::generate_algebra<real, 3>::type>() ;
        battery<dickson::generate_algebra<real, 4>::type>() ;
        battery<dickson::generate_algebra<real, 5>::type>() ;
        battery<dickson::generate_algebra<real, 6>::type>() ;
        battery<dickson::generate_algebra<real, 7>::type>() ;
        battery<dickson::generate_algebra<real, 8>::type>() ;
    }

    std::cout << g_assertion_count
              << " assertions succeeded."
              << std::endl ;
    
    return 0 ;
}

