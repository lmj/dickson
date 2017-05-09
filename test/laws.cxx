/////////////////////////////////////////////////////////////////////////
// 
// Copyright (C) 2006 James M. Lawrence.  All rights reserved.
// 
// software: dickson-1.0.3
// file: test/laws.cxx
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
// A test of Dickson algebra identities for real-valued and
// integer-valued scalars.
//
// For large algebras which exceed the stack, #define DICKSON_USE_HEAP
// (much slower).
// 
// Note that overflow and insufficient precision can produce bogus
// failure results.
// 

//#define DICKSON_USE_HEAP
#include "common.hxx"
#include <iostream>
#include <limits>
#include <cstdlib>
#include <ctime>

void section( const char* title ) ;

template< typename scalar, unsigned int log2_dimension >
void battery( scalar epsilon = scalar(0) ) ;

int main()
{
    // Choose precision.
    typedef long double real ;
    typedef long int integer ;

    // Choose epsilon comparison for reals.
    const real epsilon = 1e-6 ;

    //
    // Choose a large algebra for testing; dimension is 2^n.
    //
    // Depending on how high you set this, you may need to define
    // DICKSON_USE_HEAP to avoid blowing the stack.  (I could handle
    // up to 2^14 = 16384 dimensions with 96-bit long doubles.)
    //
    // Incidentally, n can be zero.  1-dimensional scalars (2^0) are
    // valid dickson algebras.
    //
    const unsigned int n = 5 ;
    
    section("real scalars") ;
    battery<real, n>(epsilon) ;

    section("integer scalars") ;
    battery<integer, n>() ;

    return 0 ;
}

///////////////////////////////////////////////////////////////
// test mode
///////////////////////////////////////////////////////////////

enum Test_Mode
{
    TEST_EQUAL,
    TEST_UNEQUAL
} ;

///////////////////////////////////////////////////////////////
// output spec
///////////////////////////////////////////////////////////////

void section( const char* title )
{
    std::cout
        << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
        << std::endl
        << title
        << std::endl
        << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
        << std::endl ;
}

void header( const char * name, Test_Mode test_mode )
{
    std::cout
        << "=========================================================="
        << std::endl 
        << name ;

    if( test_mode == TEST_UNEQUAL )
    {
        std::cout << " (verify failure with a larger algebra)" ;
    }

    std::cout
        << std::endl 
        << "=========================================================="
        << std::endl ;
}

void property( const char* name )
{
    std::cout
        << "-----------------------"
        << name
        << std::endl ;
}

template< typename T >
void puts( const T & a )
{
    std::cout << a << std::endl ;
}

//
// This has to be a macro since we are printing #x and #y.
//
// Depend on locally defined 'scalar', 'test_mode', and 'epsilon'
//

#define test_pair(x, y)                                 \
    {                                                   \
        puts(#x) ;                                      \
                                                        \
        if( test_mode == TEST_EQUAL )                   \
        {                                               \
            puts("    ==") ;                            \
        }                                               \
        else                                            \
        {                                               \
            puts("    !=") ;                            \
        }                                               \
                                                        \
        puts(#y) ;                                      \
                                                        \
        std::cout << "..." << std::flush ;              \
                                                        \
        if( within((x),(y),epsilon)                     \
            ==                                          \
            (test_mode == TEST_EQUAL) )                 \
        {                                               \
            puts("true.") ;                             \
        }                                               \
        else                                            \
        {                                               \
            puts("**********FALSE**********!!") ;       \
            puts(">>>>>>>>") ;                          \
            puts(x) ;                                   \
            puts(">>>>>>>>") ;                          \
            puts(y) ;                                   \
            puts(">>>>>>>>") ;                          \
            exit(1) ;                                   \
        }                                               \
        puts("") ;                                      \
    }                                                   \


/////////////////////////////////////////////////////////
// pull names for one-dimensional algebra test (double, etc)
/////////////////////////////////////////////////////////

using namespace dickson::scalar_func ;

///////////////////////////////////////////////////////////////
//
// dickson laws
//
///////////////////////////////////////////////////////////////

template< typename algebra >
void dickson_laws( typename dickson::traits<algebra>::scalar epsilon,
                   Test_Mode test_mode )
{
    typedef typename dickson::traits<algebra>::scalar scalar ;

    algebra a = rand<algebra>() ;
    algebra b = rand<algebra>() ;
    algebra c = rand<algebra>() ;

    const unsigned int m = 3 ;
    const unsigned int n = 5 ;
    
    header("dickson laws", test_mode) ;

    property("elements have scalar part and pure part") ;
    test_pair(a,
              dot(scalar(1), a) + cross(scalar(1), a)) ;

    property("conjugate negates pure part") ;
    test_pair(conj(a),
              dot(scalar(1), a) - cross(scalar(1), a)) ;

    property("dot is scalar part of b*conj(a)") ;
    test_pair(dot(a, b),
              (b*conj(a) + a*conj(b))/scalar(2)) ;

    property("cross is pure part of b*conj(a)") ;
    test_pair(cross(a, b),
              (b*conj(a) - a*conj(b))/scalar(2)) ;

    property("multiplication is dot minus cross, conj second term") ;
    test_pair(a*b,
              dot(a, conj(b)) - cross(a, conj(b))) ;

    property("cross is anti-commutative") ;
    test_pair(cross(a, b),
              -cross(b, a)) ;

    property("cross is nilpotent") ;
    test_pair(cross(a, a),
              scalar(0)) ;

    property("a*b - b*a") ;
    test_pair(a*b - b*a,
              cross(a, b) + cross(conj(a), conj(b))) ;
    test_pair(a*b - b*a,
              - cross(a, conj(b)) - cross(conj(a), b)) ;
    test_pair(a*b - b*a,
              b*conj(a) - conj(a)*b) ;
    test_pair(a*b - b*a,
              conj(b)*a - a*conj(b)) ;

    property("cross of pures is orthogonal to both") ;
    test_pair(dot(purepart(a), cross(purepart(a), purepart(b))),
              dot(purepart(b), cross(purepart(a), purepart(b)))) ;

    property("norm") ;
    test_pair(norm(a),
              norm(conj(a))) ;
    test_pair(norm(a),
              dot(a, a)) ;
    test_pair(norm(a),
              a*conj(a)) ;
    test_pair(norm(a),
              conj(a)*a) ;
    test_pair(norm(a),
              scalar(realpart(a)*realpart(a) + norm(purepart(a)))) ;

    property("conjugation of product") ;
    test_pair(conj(a*b),
              conj(b)*conj(a)) ;

    property("norm of product") ;
    test_pair(norm(a*b),
              (a*b)*(conj(b)*conj(a))) ;
    test_pair(norm(a*b),
              norm(a*conj(b))) ;
    test_pair(norm(a*b),
              norm(conj(a)*b)) ;
    test_pair(norm(a*b),
              norm(b*a)) ;

    property("norm of addition") ;
    test_pair(norm(a + b),
              scalar(norm(a) + norm(b) + scalar(2)*dot(a, b))) ;

    property("braid") ;
    test_pair(dot(a*b, c),
              dot(a, c*conj(b))) ;
    test_pair(dot(a*b, c),
              dot(b, conj(a)*c)) ;

    property("bi-multiplication is well-defined") ;
    test_pair((a*b)*a,
              a*(b*a)) ;
    test_pair((a*b)*conj(a),
              a*(b*conj(a))) ;
    test_pair((pow(a, m)*b)*pow(a, n),
              pow(a, m)*(b*pow(a, n))) ;

    property("associative when outside pure parts are parallel") ;
    {
        algebra c = rand<scalar>() + rand<scalar>()*purepart(a) ;
        test_pair((a*b)*c,
                  a*(b*c)) ;
    }

    property("consecutive powers may be swapped") ;
    test_pair(conj(a)*(a*b),
              a*(conj(a)*b)) ;
    test_pair((a*b)*conj(b),
              (a*conj(b))*b) ;
    test_pair(pow(a, m)*(pow(a, n)*b),
              pow(a, n)*(pow(a, m)*b)) ;
    test_pair((a*pow(b, m))*pow(b, n),
              (a*pow(b, n))*pow(b, m)) ;

    property("consecutive conjugate pairs may be reversed") ;
    test_pair(conj(a)*(a*b),
              (b*a)*conj(a)) ;
    test_pair((a*b)*conj(b),
              conj(b)*(b*a)) ;

    property("conjugate both terms in dot") ;
    test_pair(dot(a, b),
              dot(conj(a), conj(b))) ;
    
    property("distributive dot") ;
    test_pair(dot(a, b + c),
              scalar(dot(a, b) + dot(a, c))) ;

    property("distributive cross") ;
    test_pair(cross(a, b + c),
              cross(a, b) + cross(a, c)) ;

    property("dot dot") ;
    test_pair(dot(a, dot(b, c)),
              scalar((dot(a*b, c) + dot(a*c, b))/scalar(2))) ;

    property("dot cross") ;
    test_pair(dot(a, cross(b, c)),
              scalar((dot(a*b, c) - dot(a*c, b))/scalar(2))) ;

    property("a*b + b*a") ;
    test_pair(a*b + b*a,
              - cross(a, conj(b))
              + cross(conj(a), b)
              + scalar(2)*dot(a, conj(b))) ;
}

///////////////////////////////////////////////////////////////
//
// octonion laws
//
///////////////////////////////////////////////////////////////

template< typename algebra >
void octonion_laws_common( typename dickson::traits<algebra>::scalar epsilon,
                           Test_Mode test_mode )
{
    typedef typename dickson::traits<algebra>::scalar scalar ;

    algebra a = rand<algebra>() ;
    algebra b = rand<algebra>() ;
    algebra c = rand<algebra>() ;
    algebra d = rand<algebra>() ;

    const unsigned int m = 3 ;
    const unsigned int n = 5 ;
    
    header("octonion laws", test_mode) ;

    property("norm of product") ;
    test_pair(norm(a*b),
              scalar(norm(a)*norm(b))) ;

    property("consecutive powers are well-defined") ;
    test_pair(pow(a, m)*(pow(a, n)*b),
              pow(a, m + n)*b) ;
    test_pair((a*pow(b, m))*pow(b, n),
              a*pow(b, m + n)) ;

    property("dot times dot") ;
    test_pair(scalar(dot(a, b)*dot(c, d)),
              scalar((dot(a*c, b*d) + dot(a*d, b*c))/scalar(2))) ;

    property("moufang triplets") ;
    test_pair((a*b)*(c*a),
              a*(b*c)*a) ;
              
    property("cross dot") ;
    test_pair(cross(a, dot(b, c)),
              (cross(a*b, c) + cross(a*c, b))/scalar(2)) ;

    property("dot of squares") ;
    test_pair(dot(a*a, b*b),
              scalar(dot(a, b)*dot(a, b) +
                     dot(cross(a, b), cross(conj(a), conj(b))))) ;

    property("split a norm when outer pure parts parallel") ;
    {
        algebra c = rand<scalar>() + rand<scalar>()*purepart(a) ;
        test_pair(a*norm(b)*c,
                  (a*b)*(conj(b)*c)) ;
    }
}

template< typename algebra, bool is_integer >
struct select_octonion_laws ;

template< typename algebra >
struct select_octonion_laws<algebra, true>
{
    static void laws( typename dickson::traits<algebra>::scalar epsilon,
                      Test_Mode test_mode )
    {
        octonion_laws_common<algebra>(epsilon, test_mode) ;
    }
} ;

template< typename algebra >
struct select_octonion_laws<algebra, false>
{
    static void laws( typename dickson::traits<algebra>::scalar epsilon,
                      Test_Mode test_mode ) ;
} ;

template< typename algebra >
void
select_octonion_laws<algebra, false>::
laws( typename dickson::traits<algebra>::scalar epsilon,
      Test_Mode test_mode )
{
    octonion_laws_common<algebra>(epsilon, test_mode) ;
    
    typedef typename dickson::traits<algebra>::scalar scalar ;
    
    algebra a = rand<algebra>() ;
    algebra b = rand<algebra>() ;
    algebra c = rand<algebra>() ;
    
    property("inverse of product") ;
    test_pair(inv(a*b),
              inv(b)*inv(a)) ;
        
    property("moufang triplets") ;
    test_pair((a*b)*c,
              (a*inv(c))*(c*b*c)) ;
    test_pair(a*(b*c),
              (a*b*a)*(inv(a)*c)) ;
    
    property("conjugate terms in cross") ;
    test_pair(cross(a, b),
              -a*cross(conj(a), conj(b))*inv(a)) ;
    
    property("a*b + b*a") ;
    test_pair(a*b + b*a,
              a
              *
              (cross(a, b) -
               cross(conj(a), conj(b)) + 
               scalar(2)*dot(a, b))
              *
              inv(conj(a))) ;

    property("exponent similarity") ;
    test_pair(a*exp(b)*inv(a),
              exp(a*b*inv(a))) ;
}

template< typename algebra >
void octonion_laws( typename dickson::traits<algebra>::scalar epsilon,
                    Test_Mode test_mode )
{
    select_octonion_laws
    <
        algebra,
        std::numeric_limits
        <
            typename dickson::traits<algebra>::scalar
        >::is_integer
    >::laws(epsilon, test_mode) ;
}

///////////////////////////////////////////////////////////////
//
// quaternion laws
//
///////////////////////////////////////////////////////////////

template< typename algebra >
void quaternion_laws_common( typename dickson::traits<algebra>::scalar epsilon,
                             Test_Mode test_mode )
{
    typedef typename dickson::traits<algebra>::scalar scalar ;

    algebra a = rand<algebra>() ;
    algebra b = rand<algebra>() ;
    algebra c = rand<algebra>() ;
    algebra d = rand<algebra>() ;

    header("quaternion laws", test_mode) ;

    property("associative") ;
    test_pair(a*(b*c),
              (a*b)*c) ;

    property("cross braid") ;
    test_pair(cross(a*b, c),
              cross(a, c*conj(b))) ;

    property("cross cross") ;
    test_pair(cross(a, cross(b, c)),
              (cross(a*b, c) - cross(a*c, b))/scalar(2)) ;

    property("four crosses") ;
    test_pair(cross(cross(a, b), cross(c, d)),
              dot(a, c)*cross(b, d) - dot(a, d)*cross(b, c) +
              dot(b, d)*cross(a, c) - dot(b, c)*cross(a, d)) ;
}

template< typename algebra, bool is_integer >
struct select_quaternion_laws ;

template< typename algebra >
struct select_quaternion_laws<algebra, true>
{
    static void laws( typename dickson::traits<algebra>::scalar epsilon,
                      Test_Mode test_mode )
    {
        quaternion_laws_common<algebra>(epsilon, test_mode) ;
    }
} ;

template< typename algebra >
struct select_quaternion_laws<algebra, false>
{
    static void laws( typename dickson::traits<algebra>::scalar epsilon,
                      Test_Mode test_mode ) ;
} ;

template< typename algebra >
void
select_quaternion_laws<algebra, false>::
laws( typename dickson::traits<algebra>::scalar epsilon,
      Test_Mode test_mode )
{
    quaternion_laws_common<algebra>(epsilon, test_mode) ;
    
    typedef typename dickson::traits<algebra>::scalar scalar ;
    
    algebra a = rand<algebra>() ;
    algebra b = rand<algebra>() ;
    algebra c = rand<algebra>() ;
    
    property("cross braid") ;
    test_pair(cross(a*b, c),
              a*cross(b, conj(a)*c)*inv(a)) ;
}

template< typename algebra >
void quaternion_laws( typename dickson::traits<algebra>::scalar epsilon,
                      Test_Mode test_mode )
{
    select_quaternion_laws
    <
        algebra,
        std::numeric_limits
        <
            typename dickson::traits<algebra>::scalar
        >::is_integer
    >::laws(epsilon, test_mode) ;
}

///////////////////////////////////////////////////////////////
//
// complex laws
//
///////////////////////////////////////////////////////////////

template< typename algebra >
void complex_laws( typename dickson::traits<algebra>::scalar epsilon,
                   Test_Mode test_mode )
{
    typedef typename dickson::traits<algebra>::scalar scalar ;

    algebra a = rand<algebra>() ;
    algebra b = rand<algebra>() ;

    header("complex laws", test_mode) ;
    
    property("commutative") ;
    test_pair(a*b,
              b*a) ;
}

///////////////////////////////////////////////////////////////
//
// battery
//
///////////////////////////////////////////////////////////////

template< typename scalar, unsigned int log2_dimension >
void battery( scalar epsilon )
{
    std::srand(std::time(0)) ;
    std::cout.precision(std::numeric_limits<scalar>::digits10 + 1) ;

    typedef typename dickson::named_algebras<scalar>::complex complex ;
    typedef typename dickson::named_algebras<scalar>::quaternion quaternion ;
    typedef typename dickson::named_algebras<scalar>::octonion octonion ;
    typedef typename dickson::named_algebras<scalar>::sedenion sedenion ;
    
    typedef typename 
        dickson::generate_algebra<scalar, log2_dimension>::type
        big_algebra ;
    
    dickson_laws<big_algebra>(epsilon, TEST_EQUAL) ; 
    octonion_laws<octonion>(epsilon, TEST_EQUAL) ;
    quaternion_laws<quaternion>(epsilon, TEST_EQUAL) ;
    complex_laws<complex>(epsilon, TEST_EQUAL) ;
    
    octonion_laws<sedenion>(epsilon, TEST_UNEQUAL) ;
    quaternion_laws<octonion>(epsilon, TEST_UNEQUAL) ;
    complex_laws<quaternion>(epsilon, TEST_UNEQUAL) ;
}

