/////////////////////////////////////////////////////////////////////////
// 
// Copyright (C) 2006 James M. Lawrence.  All rights reserved.
// 
// software: dickson-1.0.3
// file: dickson/dickson.hxx
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
// "Dickson double algebras" ---
//    complex, quaternions, octonions, sedenions, etc
//
// For large algebras which exceed the stack, #define DICKSON_USE_HEAP
// (much slower).
//
// Reference:
//   "On Quaternions and Octonions," John H. Conway and Derek A. Smith.
//

#ifndef dickson_dickson_hxx
#define dickson_dickson_hxx

#include <cmath>
#include <istream>
#include <ostream>
#include <limits>
#include <cstdlib>
#include <cassert>
#include <new>

#define dickson_assert(x)

namespace dickson {

/////////////////////////////////////////////////////////
// algebra (fwd)
/////////////////////////////////////////////////////////

template< typename subalg >
class algebra ;

////////////////////////////////////////////////////////
//
// named_algebras (convenience)
//
////////////////////////////////////////////////////////

template< typename scalar >
struct named_algebras
{
    typedef algebra< scalar >      complex ;
    typedef algebra< complex >     quaternion ;
    typedef algebra< quaternion >  octonion ;
    typedef algebra< octonion >    sedenion ;
} ;

////////////////////////////////////////////////////////
//
// generate_algebra<scalar, log2_dimension>::type
//
// For example the octonion class (2^3 == 8 dimensions) is
//
//     generate_algebra<double, 3>::type 
//
////////////////////////////////////////////////////////

template< typename scalar, unsigned int log2_dimension >
struct generate_algebra
{
    typedef
    algebra
    <
        typename
        generate_algebra
        <
            scalar,
            log2_dimension - 1
        >::type
    >
    type ;
} ;

template< typename scalar >
struct generate_algebra< scalar, 0 >
{
    typedef scalar type ;
} ;

/////////////////////////////////////////////////////////
// traits (fwd)
/////////////////////////////////////////////////////////

template< typename subalg >
class traits ;

/////////////////////////////////////////////////////////
//
// algebra
//
/////////////////////////////////////////////////////////

template< typename T_subalg >
class algebra
{
public:
    typedef T_subalg subalg ;
    typedef typename traits<subalg>::scalar scalar ;

    static const unsigned int dimension = traits<algebra>::dimension ;
    
    // Construct (0, 0, 0, ..., 0).
    algebra() ;

    // Construct (s, 0, 0, ..., 0).
    explicit algebra( scalar s ) ;

    // Construct (a + b i) where 'i' is the appropriate unit element
    // orthogonal to a.
    algebra( const subalg & a, const subalg & b ) ;

    // The first subalgebra; corresponds to 'a' in algebra(a, b).
    const subalg & first() const ;
          subalg & first() ;

    // The second subalgebra; corresponds to 'b' in algebra(a, b).
    const subalg & second() const ;
          subalg & second() ;

    // Scalar coefficients: 0...(dimension - 1).
    scalar   operator[]( unsigned int index ) const ;
    scalar & operator[]( unsigned int index ) ;

    // (*this)[0]
    scalar   realpart() const ;
    scalar & realpart() ;

    // member operators
    algebra & operator+=( const algebra & a ) ;
    algebra & operator-=( const algebra & a ) ;
    algebra & operator*=( const algebra & a ) ;
    algebra & operator+=( scalar a ) ;
    algebra & operator-=( scalar a ) ;
    algebra & operator*=( scalar a ) ;
    algebra & operator/=( scalar a ) ;

#if !defined(DICKSON_USE_HEAP)

private:
    subalg m_first ;
    subalg m_second ;

#elif defined(DICKSON_USE_HEAP)

public:
    algebra( const algebra & ) ;
    algebra & operator=( const algebra & ) ;
    ~algebra() ;

private:
    template< typename TF_algebra >
    friend class algebra_helper ;

    bool m_toplevel ;
    char* m_data ;

#endif // defined(DICKSON_USE_HEAP)
} ;

/////////////////////////////////////////////////////////
//
// traits
//
/////////////////////////////////////////////////////////

template< typename T_scalar >
struct traits
{
    typedef T_scalar scalar ;

    static const unsigned int dimension = 1 ;
    static const unsigned int log2_dimension = 0 ;
} ;

template< typename T_subalg >
struct traits< algebra<T_subalg> >
{
    typedef T_subalg subalg ;
    typedef typename traits<subalg>::scalar scalar ;

    static const unsigned int dimension =
    2*traits<subalg>::dimension ;

    static const unsigned int log2_dimension =
    1 + traits<subalg>::log2_dimension ;
} ;

/////////////////////////////////////////////////////////
// non-member functions (fwd)
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
operator+( const algebra<subalg> & a,
           const algebra<subalg> & b ) ;

template< typename subalg >
algebra<subalg>
operator+( typename traits<subalg>::scalar a,
           const algebra<subalg> & b ) ;

template< typename subalg >
algebra<subalg>
operator+( const algebra<subalg> & a,
           typename traits<subalg>::scalar b ) ;

template< typename subalg >
const algebra<subalg> &
operator+( const algebra<subalg> & a ) ;

template< typename subalg >
algebra<subalg>
operator-( const algebra<subalg> & a,
           const algebra<subalg> & b ) ;

template< typename subalg >
algebra<subalg>
operator-( typename traits<subalg>::scalar a,
           const algebra<subalg> & b ) ;

template< typename subalg >
algebra<subalg>
operator-( const algebra<subalg> & a,
           typename traits<subalg>::scalar b ) ;

template< typename subalg >
algebra<subalg>
operator-( const algebra<subalg> & a ) ;

template< typename subalg >
algebra<subalg>
operator*( typename traits<subalg>::scalar a,
           const algebra<subalg> & b ) ;

template< typename subalg >
algebra<subalg>
operator*( const algebra<subalg> & a,
           typename traits<subalg>::scalar b ) ;

template< typename subalg >
algebra<subalg>
operator*( const algebra<subalg> & a,
           const algebra<subalg> & b ) ;

template< typename subalg >
algebra<subalg>
operator/( const algebra<subalg> & a,
           typename traits<subalg>::scalar b ) ;

template< typename subalg >
bool
operator==( const algebra<subalg> & a,
            const algebra<subalg> & b ) ;

template< typename subalg >
bool
operator==( typename algebra<subalg>::scalar a,
            const algebra<subalg> & b ) ;

template< typename subalg >
bool
operator==( const algebra<subalg> & a,
            typename algebra<subalg>::scalar b ) ;

template< typename subalg >
bool
operator!=( const algebra<subalg> & a,
            const algebra<subalg> & b ) ;

template< typename subalg >
bool
operator!=( typename algebra<subalg>::scalar a,
            const algebra<subalg> & b ) ;

template< typename subalg >
bool
operator!=( const algebra<subalg> & a,
            typename algebra<subalg>::scalar b ) ;

// realpart(a) == a[0]
template< typename subalg >
typename traits<subalg>::scalar
realpart( const algebra<subalg> & a ) ;

// purepart(a) == a - realpart(a)
template< typename subalg >
algebra<subalg>
purepart( const algebra<subalg> & a ) ;

// conj(a) == a - 2*purepart(a)
template< typename subalg >
algebra<subalg>
conj( const algebra<subalg> & a ) ;

// norm(a) == a*conj(a)
template< typename subalg >
typename traits<subalg>::scalar
norm( const algebra<subalg> & a ) ;

// abs(a) == sqrt(norm(a))
template< typename subalg >
typename traits<subalg>::scalar
abs( const algebra<subalg> & a ) ;

// a*inv(a) == 1
template< typename subalg >
algebra<subalg>
inv( const algebra<subalg> & a ) ;

// dot(a) == (b*conj(a) + a*conj(b))/2

template< typename subalg >
typename traits<subalg>::scalar
dot( const algebra<subalg> & a,
     const algebra<subalg> & b ) ;

template< typename subalg >
typename traits<subalg>::scalar
dot( typename traits<subalg>::scalar a,
     const algebra<subalg> & b ) ;

template< typename subalg >
typename traits<subalg>::scalar
dot( const algebra<subalg> & a,
     typename traits<subalg>::scalar b ) ;

// cross(a) == (b*conj(a) - a*conj(b))/2

template< typename subalg >
algebra<subalg>
cross( const algebra<subalg> & a,
       const algebra<subalg> & b ) ;

template< typename subalg >
algebra<subalg>
cross( typename traits<subalg>::scalar a,
       const algebra<subalg> & b ) ;

template< typename subalg >
algebra<subalg>
cross( const algebra<subalg> & a,
       typename traits<subalg>::scalar b ) ;

// power-series family defined as usual for Lie algebras

template< typename subalg >
algebra<subalg>
exp( const algebra<subalg> & a ) ;

template< typename subalg >
algebra<subalg>
log( const algebra<subalg> & a ) ;

template< typename subalg >
algebra<subalg>
cos( const algebra<subalg> & a ) ;

template< typename subalg >
algebra<subalg>
sin( const algebra<subalg> & a ) ;

template< typename subalg >
algebra<subalg>
tan( const algebra<subalg> & a ) ;

template< typename subalg >
algebra<subalg>
cosh( const algebra<subalg> & a ) ;

template< typename subalg >
algebra<subalg>
sinh( const algebra<subalg> & a ) ;

template< typename subalg >
algebra<subalg>
tanh( const algebra<subalg> & a ) ;

// note: expect pow(0,0) to return some random nonsensical value

template< typename subalg >
algebra<subalg>
pow( const algebra<subalg> & a, unsigned int n ) ;

template< typename subalg >
algebra<subalg>
pow( const algebra<subalg> & a, typename traits<subalg>::scalar gamma ) ;

template< typename subalg >
algebra<subalg>
pow( typename traits<subalg>::scalar gamma, const algebra<subalg> & a ) ;

// sqrt(a) == pow(a, 0.5)
template< typename subalg >
algebra<subalg>
sqrt( const algebra<subalg> & a ) ;

// streams

template< typename subalg, typename T_Char, typename T_Traits >
std::basic_ostream<T_Char, T_Traits> &
operator<<( std::basic_ostream<T_Char, T_Traits> & out,
            const algebra<subalg> & a ) ;

template< typename subalg, typename T_Char, typename T_Traits >
std::basic_istream<T_Char, T_Traits> &
operator>>( std::basic_istream<T_Char, T_Traits> & in,
            algebra<subalg> & a ) ;

/////////////////////////////////////////////////////////////////////////
//
// scalar_func
//
// A new scalar may be introduced by adding the appropriate functions
// to scalar_func (sin, cos, exp, etc.).
//
// The separate namespace also provides the convenience of a 'using'
// directive which pulls in only these functions which do not have a
// class signature.
//
/////////////////////////////////////////////////////////////////////////

namespace scalar_func {

template< typename scalar >
scalar
realpart( scalar a ) ;

template< typename scalar >
scalar
purepart( scalar a ) ;

template< typename scalar >
scalar
conj( scalar a ) ;

template< typename scalar >
scalar
norm( scalar a ) ;

template< typename scalar >
scalar
inv( const scalar a ) ;

template< typename scalar >
scalar
dot( scalar a, scalar b ) ;

template< typename scalar >
scalar
cross( scalar a, scalar b ) ;

template< typename scalar >
scalar
pow( scalar a, scalar n ) ;

template< typename scalar >
scalar
pow( scalar a, unsigned int n ) ;

using std::abs ;
using std::exp ;
using std::log ;
using std::cos ;
using std::sin ;
using std::cosh ;
using std::sinh ;
using std::pow ;
using std::sqrt ;
using std::atan2 ;

} // namespace scalar_func

// needed for compile-time recursion
using scalar_func::dot ;
using scalar_func::conj ;
using scalar_func::realpart ;
using scalar_func::norm ;

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//
// implementation
//
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////
// stack-based construction
/////////////////////////////////////////////////////////

#if !defined(DICKSON_USE_HEAP)

template< typename subalg >
inline
algebra<subalg>::
algebra()
    : m_first(),
      m_second()
{
}

template< typename subalg >
inline
algebra<subalg>::
algebra( const subalg & a,
         const subalg & b )
    : m_first(a),
      m_second(b)
{
}

template< typename subalg >
inline
algebra<subalg>::
algebra( scalar a )
    : m_first(subalg(a)),
      m_second()
{
}

template< typename subalg >
inline
const subalg &
algebra<subalg>::
first() const
{
    return m_first ;
}
    
template< typename subalg >
inline
subalg &
algebra<subalg>::
first()
{
    return m_first ;
}

template< typename subalg >
inline
const subalg &
algebra<subalg>::
second() const
{
    return m_second ;
}
    
template< typename subalg >
inline
subalg &
algebra<subalg>::
second()
{
    return m_second ;
}

template< typename subalg >
inline
typename algebra<subalg>::scalar
algebra<subalg>::
operator[]( unsigned int index ) const
{
    dickson_assert(index < dimension) ;
    return reinterpret_cast<const scalar*>(this)[index] ;
}

template< typename subalg >
inline
typename algebra<subalg>::scalar &
algebra<subalg>::
operator[]( unsigned int index )
{
    dickson_assert(index < dimension) ;
    return reinterpret_cast<scalar*>(this)[index] ;
}

template< typename subalg >
inline
typename algebra<subalg>::scalar
algebra<subalg>::
realpart() const
{
    return *reinterpret_cast<const scalar*>(this) ;
}

template< typename subalg >
inline
typename algebra<subalg>::scalar &
algebra<subalg>::
realpart()
{
    return *reinterpret_cast<scalar*>(this) ;
}

#endif // !defined(DICKSON_USE_HEAP)

/////////////////////////////////////////////////////////
// heap-based construction
/////////////////////////////////////////////////////////
//
//
// top-level construction: 
//
//           +--------+
//           |        |
// this ---->| m_data |
//           |        |
//           +--------+
//               |
//               |
//               | malloc()
//               |
//              \|/
//               *
//       +----------------+----------------+
//       |                |                |
//       |  first subalg  |  second subalg |
//       |                |                |
//       +----------------+----------------+
//
//
// nested construction:
//
//           +--------+----------------+----------------+
//           |        |                |                |
// this ---->| m_data |  first subalg  |  second subalg |
//           |        |                |                |
//           +--------+----------------+----------------+
//               |           *
//               |          /|\                             .
//               |           |
//               |           |
//               +-----------+
//
//
// It would be nice to do away with the top-level flag, but I don't
// see a good alternative.
//
// Why not just check if 'this' is adjacent to m_data?  Because then
// there is a possibility of leaking memory.
//
// Why not just set m_data to null for nested data?  Because then you
// have to continually check for null, fouling up critical inlining
// and making things slow.
//
// Why not just make a separate nested_algebra class with the
// appropriate conversion operators?  Because conversion operators do
// not work with templates.
//

#if defined(DICKSON_USE_HEAP)

template< typename scalar >
class algebra_helper
{
private:
    friend class algebra_helper< algebra<scalar> > ;
    
    static const unsigned int chunk_size = sizeof(scalar) ;

    static scalar at( const scalar* self, unsigned int index )
    {
        dickson_assert(index == 0) ;
        return *self ;
    }

    static scalar & at( scalar* self, unsigned int index )
    {
        dickson_assert(index == 0) ;
        return *self ;
    }

    static void construct_nested( scalar* self )
    {
        *self = scalar(0) ;
    }

    static void construct_nested( scalar* self, scalar a )
    {
        *self = a ;
    }

    static scalar realpart( const scalar* self )
    {
        return *self ;
    }

    static scalar & realpart( scalar* self )
    {
        return *self ;
    }
} ;

template< typename subalg >
class algebra_helper< algebra<subalg> >
{
private:
    friend class algebra<subalg> ;
    friend class algebra_helper< algebra< algebra<subalg> > > ;
    
    typedef typename traits<subalg>::scalar scalar ;

    static const unsigned int chunk_size =
    sizeof(algebra<subalg>) + 2*algebra_helper<subalg>::chunk_size ;
    
    static const subalg* first_ptr( const algebra<subalg>* self )
    {
        return reinterpret_cast<const subalg*>(self->m_data) ;
    }
    
    static subalg* first_ptr( algebra<subalg>* self )
    {
        return reinterpret_cast<subalg*>(self->m_data) ;
    }
    
    static const subalg* second_ptr( const algebra<subalg>* self )
    {
        return reinterpret_cast<const subalg*>(
            self->m_data
            +
            algebra_helper<subalg>::chunk_size) ;
    }

    static subalg* second_ptr( algebra<subalg>* self )
    {
        return reinterpret_cast<subalg*>(
            self->m_data
            +
            algebra_helper<subalg>::chunk_size) ;
    }

    static void allocate( algebra<subalg>* self )
    {
        self->m_data = static_cast<char*>(
            std::malloc(2*algebra_helper<subalg>::chunk_size)) ;

        self->m_toplevel = true ;
    }
    
    static void nest( algebra<subalg>* self )
    {
        self->m_data = 
            reinterpret_cast<char*>(self)
            +
            sizeof(algebra<subalg>) ;
        
        self->m_toplevel = false ;
    }
    
    static void construct( algebra<subalg>* self )
    {
        allocate(self) ;
        algebra_helper<subalg>::construct_nested(first_ptr(self)) ;
        algebra_helper<subalg>::construct_nested(second_ptr(self)) ;
    }

    static void construct( algebra<subalg>* self, scalar a )
    {
        allocate(self) ;
        algebra_helper<subalg>::construct_nested(first_ptr(self), a) ;
        algebra_helper<subalg>::construct_nested(second_ptr(self)) ;
    }

    static void construct( algebra<subalg>* self,
                           const subalg & a,
                           const subalg & b )
    {
        allocate(self) ;
        algebra_helper<subalg>::construct_nested(first_ptr(self),  a) ;
        algebra_helper<subalg>::construct_nested(second_ptr(self), b) ;
    }

    static void construct( algebra<subalg>* self,
                           const algebra<subalg> & a )
    {
        allocate(self) ;
        algebra_helper<subalg>::construct_nested(first_ptr(self),
                                                 a.first()) ;
        algebra_helper<subalg>::construct_nested(second_ptr(self),
                                                 a.second()) ;
    }

    static void construct_nested( algebra<subalg>* self )
    {
        nest(self) ;
        algebra_helper<subalg>::construct_nested(first_ptr(self)) ;
        algebra_helper<subalg>::construct_nested(second_ptr(self)) ;
    }

    static void construct_nested( algebra<subalg>* self, scalar a )
    {
        nest(self) ;
        algebra_helper<subalg>::construct_nested(first_ptr(self), a) ;
        algebra_helper<subalg>::construct_nested(second_ptr(self)) ;
    }

    static void construct_nested( algebra<subalg>* self,
                                  const subalg & a,
                                  const subalg & b )
    {
        nest(self) ;
        algebra_helper<subalg>::construct_nested(first_ptr(self),  a) ;
        algebra_helper<subalg>::construct_nested(second_ptr(self), b) ;
    }
    
    static void construct_nested( algebra<subalg>* self,
                                  const algebra<subalg> & a )
    {
        nest(self) ;
        algebra_helper<subalg>::construct_nested(first_ptr(self),
                                                 a.first()) ;
        algebra_helper<subalg>::construct_nested(second_ptr(self),
                                                 a.second()) ;
    }
    
    static void destruct( algebra<subalg>* self )
    {
        if( self->m_toplevel )
        {
            std::free(self->m_data) ;
        }
    }

    static scalar at( const algebra<subalg>* self, unsigned int index )
    {
        // strip the most significant bit
        const unsigned int subindex =
            ((algebra<subalg>::dimension - 1) >> 1)
            &
            index ;
        
        return 
            index == subindex
            ?
            algebra_helper<subalg>::at(first_ptr(self), subindex)
            :
            algebra_helper<subalg>::at(second_ptr(self), subindex) ;
    }

    static scalar & at( algebra<subalg>* self, unsigned int index )
    {
        // strip the most significant bit
        const unsigned int subindex =
            ((algebra<subalg>::dimension - 1) >> 1)
            &
            index ;
        
        return 
            index == subindex
            ?
            algebra_helper<subalg>::at(first_ptr(self), subindex)
            :
            algebra_helper<subalg>::at(second_ptr(self), subindex) ;
    }

    static scalar realpart( const algebra<subalg>* self )
    {
        return algebra_helper<subalg>::realpart(first_ptr(self)) ;
    }

    static scalar & realpart( algebra<subalg>* self )
    {
        return algebra_helper<subalg>::realpart(first_ptr(self)) ;
    }
} ;

template< typename subalg >
inline
algebra<subalg>::
algebra()
    // data initialized by algebra_helper
{
    algebra_helper<algebra>::construct(this) ;
}

template< typename subalg >
inline
algebra<subalg>::
algebra( const subalg & a,
         const subalg & b )
    // data initialized by algebra_helper
{
    algebra_helper<algebra>::construct(this, a, b) ;
}

template< typename subalg >
inline
algebra<subalg>::
algebra( scalar a )
    // data initialized by algebra_helper
{
    algebra_helper<algebra>::construct(this, a) ;
}

template< typename subalg >
inline
algebra<subalg>::
algebra( const algebra & a )
    // data initialized by algebra_helper
{
    algebra_helper<algebra>::construct(this, a) ;
}

template< typename subalg >
inline
algebra<subalg> &
algebra<subalg>::
operator=( const algebra & a )
{
    algebra::first() = a.first() ;
    algebra::second() = a.second() ;
    return *this ;
}

template< typename subalg >
inline
algebra<subalg>::
~algebra()
{
    algebra_helper<algebra>::destruct(this) ;
}

template< typename subalg >
inline
const subalg &
algebra<subalg>::
first() const
{
    return *algebra_helper<algebra>::first_ptr(this) ;
}

template< typename subalg >
inline
subalg &
algebra<subalg>::
first()
{
    return *algebra_helper<algebra>::first_ptr(this) ;
}

template< typename subalg >
inline
const subalg &
algebra<subalg>::
second() const
{
    return *algebra_helper<algebra>::second_ptr(this) ;
}
    
template< typename subalg >
inline
subalg &
algebra<subalg>::
second()
{
    return *algebra_helper<algebra>::second_ptr(this) ;
}

template< typename subalg >
inline
typename algebra<subalg>::scalar
algebra<subalg>::
operator[]( unsigned int index ) const
{
    dickson_assert(index < dimension) ;
    return algebra_helper<algebra>::at(this, index) ;
}

template< typename subalg >
inline
typename algebra<subalg>::scalar &
algebra<subalg>::
operator[]( unsigned int index )
{
    dickson_assert(index < dimension) ;
    return algebra_helper<algebra>::at(this, index) ;
}

template< typename subalg >
inline
typename algebra<subalg>::scalar
algebra<subalg>::
realpart() const
{
    return algebra_helper<algebra>::realpart(this) ;
}

template< typename subalg >
inline
typename algebra<subalg>::scalar &
algebra<subalg>::
realpart()
{
    return algebra_helper<algebra>::realpart(this) ;
}

#endif // defined(DICKSON_USE_HEAP)

/////////////////////////////////////////////////////////
// member operators
/////////////////////////////////////////////////////////

template< typename subalg >
inline
algebra<subalg> &
algebra<subalg>::
operator+=( const algebra & a )
{
    algebra::first() += a.first() ;
    algebra::second() += a.second() ;
    return *this ;
}

template< typename subalg >
inline
algebra<subalg> &
algebra<subalg>::
operator-=( const algebra & a )
{
    algebra::first() -= a.first() ;
    algebra::second() -= a.second() ;
    return *this ;
}

template< typename subalg >
inline
algebra<subalg> &
algebra<subalg>::
operator*=( const algebra & a )
{
    *this = (*this)*a ;
    return *this ;
}

template< typename subalg >
inline
algebra<subalg> &
algebra<subalg>::
operator+=( scalar a )
{
    algebra::realpart() += a ;
    return *this ;
}

template< typename subalg >
inline
algebra<subalg> &
algebra<subalg>::
operator-=( scalar a )
{
    algebra::realpart() -= a ;
    return *this ;
}

template< typename subalg >
inline
algebra<subalg> &
algebra<subalg>::
operator*=( scalar a )
{
    algebra::first() *= a ;
    algebra::second() *= a ;
    return *this ;
}

template< typename subalg >
inline
algebra<subalg> &
algebra<subalg>::
operator/=( scalar a )
{
    algebra::first() /= a ;
    algebra::second() /= a ;
    return *this ;
}

/////////////////////////////////////////////////////////
//
// dickson_private
//
/////////////////////////////////////////////////////////

namespace dickson_private {

/////////////////////////////////////////////////////////
// metaprogrammed conj()
/////////////////////////////////////////////////////////

#if 0

template< typename subalg, unsigned int index >
struct conj_helper ;

template< typename subalg, unsigned int index >
struct conj_helper<algebra<subalg>, index>
{
    static void func( algebra<subalg> & a )
    {
        typedef typename traits<subalg>::scalar scalar ;
    
        a[index] *= scalar(-1) ;
        conj_helper< algebra<subalg>, index - 1>::func(a) ;
    }
} ;

template< typename subalg >
struct conj_helper<algebra<subalg>, 0>
{
    static void func( algebra<subalg> & )
    {
    }
} ;

#endif 

/////////////////////////////////////////////////////////
// functions for real scalars only
/////////////////////////////////////////////////////////

template< typename subalg, bool is_integer >
struct select_function ;

template< typename subalg >
struct select_function<subalg, false>
{
    typedef typename traits<subalg>::scalar scalar ;

    static algebra<subalg> inv( const algebra<subalg> & a ) ;
    static algebra<subalg> exp( const algebra<subalg> & ) ;
    static algebra<subalg> log( const algebra<subalg> & ) ;
    static algebra<subalg> cos( const algebra<subalg> & ) ;
    static algebra<subalg> sin( const algebra<subalg> & ) ;
    static algebra<subalg> tan( const algebra<subalg> & ) ;
    static algebra<subalg> cosh( const algebra<subalg> & ) ;
    static algebra<subalg> sinh( const algebra<subalg> & ) ;
    static algebra<subalg> tanh( const algebra<subalg> & ) ;
    static algebra<subalg> pow( const algebra<subalg> & a, scalar gamma ) ;
    static algebra<subalg> pow( scalar gamma, const algebra<subalg> & a ) ;
    static algebra<subalg> sqrt( const algebra<subalg> & a ) ;

    static scalar abs( const algebra<subalg> & a )
    {
        return scalar_func::sqrt(norm(a)) ;
    }

    static scalar inv( scalar a )
    {
        return scalar(1)/a ;
    }

    static scalar pow( scalar a, scalar n )
    {
        return scalar_func::pow(a, n) ;
    }
} ;

/////////////////////////////////////////////////////////
// inv
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
select_function<subalg, false>::
inv( const algebra<subalg> & a )
{
    algebra<subalg> work = conj(a) ;
    work /= norm(a) ;
    return work ;
}

/////////////////////////////////////////////////////////
// exp
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
select_function<subalg, false>::
exp( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;
    
    algebra<subalg> work = a ;
    const scalar realpart_a = a.realpart() ;
    work.realpart() = scalar(0) ;
    
    if( work == scalar(0) )
    {
        return algebra<subalg>(scalar_func::exp(realpart_a)) ;
    }

    const scalar abs_purepart_a = abs(work) ;

    work *= 
        scalar_func::sin(abs_purepart_a)
        /
        abs_purepart_a ;

    work.realpart() = scalar_func::cos(abs_purepart_a) ;

    work *= scalar_func::exp(realpart_a) ;

    return work ;
}

/////////////////////////////////////////////////////////
// natural log
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
select_function<subalg, false>::
log( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;
    
    algebra<subalg> work = a ;
    const scalar realpart_a = a.realpart() ;
    work.realpart() = scalar(0) ;
    
    if( work == scalar(0) )
    {
        return algebra<subalg>(scalar_func::log(realpart_a)) ;
    }

    const scalar norm_purepart_a = norm(work) ;
    const scalar abs_a =
        scalar_func::sqrt(realpart_a*realpart_a + norm_purepart_a) ;
    const scalar abs_purepart_a = scalar_func::sqrt(norm_purepart_a) ;

    work *= 
        scalar_func::atan2(abs_purepart_a, realpart_a)
        /
        abs_purepart_a ;

    work.realpart() = scalar_func::log(abs_a) ;
    
    return work ;
}

/////////////////////////////////////////////////////////
// cos
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
select_function<subalg, false>::
cos( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;

    algebra<subalg> work = a ;
    const scalar realpart_a = a.realpart() ;
    work.realpart() = scalar(0) ;
    
    if( work == scalar(0) )
    {
        return algebra<subalg>(scalar_func::cos(realpart_a)) ;
    }

    const scalar abs_purepart_a = abs(work) ;

    work *= 
        -scalar_func::sin(realpart_a)
        *
        scalar_func::sinh(abs_purepart_a)
        /
        abs_purepart_a ;
    
    work.realpart() = 
        scalar_func::cos(realpart_a)
        *
        scalar_func::cosh(abs_purepart_a) ;
    
    return work ;
}

/////////////////////////////////////////////////////////
// sin
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
select_function<subalg, false>::
sin( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;

    algebra<subalg> work = a ;
    const scalar realpart_a = a.realpart() ;
    work.realpart() = scalar(0) ;
    
    if( work == scalar(0) )
    {
        return algebra<subalg>(scalar_func::sin(realpart_a)) ;
    }

    const scalar abs_purepart_a = abs(work) ;

    work *=
        scalar_func::cos(realpart_a)
        *
        scalar_func::sinh(abs_purepart_a)
        /
        abs_purepart_a ;

    work.realpart() =
        scalar_func::sin(realpart_a)
        *
        scalar_func::cosh(abs_purepart_a) ;
        
    return work ;
}

/////////////////////////////////////////////////////////
// tan
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
select_function<subalg, false>::
tan( const algebra<subalg> & a )
{
    return sin(a)*inv(cos(a)) ;
}

/////////////////////////////////////////////////////////
// cosh
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
select_function<subalg, false>::
cosh( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;

    algebra<subalg> work = a ;
    const scalar realpart_a = a.realpart() ;
    work.realpart() = scalar(0) ;
    
    if( work == scalar(0) )
    {
        return algebra<subalg>(scalar_func::cos(realpart_a)) ;
    }

    const scalar abs_purepart_a = abs(work) ;

    work *=
        scalar_func::sinh(realpart_a)
        *
        scalar_func::sin(abs_purepart_a)
        /
        abs_purepart_a ;
        
    work.realpart() =
        scalar_func::cosh(realpart_a)
        *
        scalar_func::cos(abs_purepart_a) ;

    return work ;
}

/////////////////////////////////////////////////////////
// sinh
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
select_function<subalg, false>::
sinh( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;

    algebra<subalg> work = a ;
    const scalar realpart_a = a.realpart() ;
    work.realpart() = scalar(0) ;
    
    if( work == scalar(0) )
    {
        return algebra<subalg>(scalar_func::cos(realpart_a)) ;
    }

    const scalar abs_purepart_a = abs(work) ;

    work *=
        scalar_func::cosh(realpart_a)
        *
        scalar_func::sin(abs_purepart_a)
        /
        abs_purepart_a ;

    work.realpart() =
        scalar_func::sinh(realpart_a)
        *
        scalar_func::cos(abs_purepart_a) ;
        
    return work ;
}

/////////////////////////////////////////////////////////
// tanh
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
select_function<subalg, false>::
tanh( const algebra<subalg> & a )
{
    return sinh(a)*inv(cosh(a)) ;
}

/////////////////////////////////////////////////////////
// pow
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
select_function<subalg, false>::
pow( const algebra<subalg> & a, scalar gamma )
{
    if( a == algebra<subalg>() )
    {
        return algebra<subalg>() ;
    }
    
    algebra<subalg> work = log(a) ;
    work *= gamma ;
    return exp(work) ;
}

template< typename subalg >
algebra<subalg>
select_function<subalg, false>::
pow( scalar gamma, const algebra<subalg> & a )
{
    if( gamma == scalar(0) )
    {
        return algebra<subalg>() ;
    }

    algebra<subalg> work = a ;
    work *= scalar_func::log(gamma) ;
    return exp(work) ;
}

/////////////////////////////////////////////////////////
// sqrt
/////////////////////////////////////////////////////////

template< typename subalg >
inline
algebra<subalg>
select_function<subalg, false>::
sqrt( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;

    return pow(a, scalar(0.5)) ;
}

} // namespace dickson_private

/////////////////////////////////////////////////////////
//
// non-memeber functions
//
/////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////
// operator+
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
operator+( const algebra<subalg> & a,
           const algebra<subalg> & b )
{
    algebra<subalg> work = a ;
    work += b ;
    return work ;
}

template< typename subalg >
algebra<subalg>
operator+( typename traits<subalg>::scalar a,
           const algebra<subalg> & b )
{
    algebra<subalg> work = b ;
    work.realpart() += a ;
    return work ;
}

template< typename subalg >
algebra<subalg>
operator+( const algebra<subalg> & a,
           typename traits<subalg>::scalar b )
{
    algebra<subalg> work = a ;
    work.realpart() += b ;
    return work ;
}

template< typename subalg >
inline
const algebra<subalg> &
operator+( const algebra<subalg> & a )
{
    return a ;
}

/////////////////////////////////////////////////////////
// operator-
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
operator-( const algebra<subalg> & a,
           const algebra<subalg> & b )
{
    algebra<subalg> work = a ;
    work -= b ;
    return work ;
}

template< typename subalg >
algebra<subalg>
operator-( typename traits<subalg>::scalar a,
           const algebra<subalg> & b )
{
    typedef typename traits<subalg>::scalar scalar ;

    algebra<subalg> work = b ;
    work *= scalar(-1) ;
    work.realpart() += a ;
    return work ;
}

template< typename subalg >
algebra<subalg>
operator-( const algebra<subalg> & a,
           typename traits<subalg>::scalar b )
{
    algebra<subalg> work = a ;
    work.realpart() -= b ;
    return work ;
}

template< typename subalg >
algebra<subalg>
operator-( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;

    algebra<subalg> work = a ;
    work *= scalar(-1) ;
    return work ;
}

/////////////////////////////////////////////////////////
// operator*
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
operator*( typename traits<subalg>::scalar a,
           const algebra<subalg> & b )
{
    // (assume commutative scalars)
    algebra<subalg> work = b ;
    work *= a ;
    return work ;
}

template< typename subalg >
algebra<subalg>
operator*( const algebra<subalg> & a,
           typename traits<subalg>::scalar b )
{
    algebra<subalg> work = a ;
    work *= b ;
    return work ;
}

#if defined(DICKSON_USE_LEFT_HANDED_MULT)
//
// multiplication with "i" on the left
//
template< typename subalg >
algebra<subalg>
operator*( const algebra<subalg> & a,
           const algebra<subalg> & b )
{
    return algebra<subalg>(a.first()*b.first()
                           -
                           b.second()*conj(a.second()),

                           b.first()*a.second()
                           +
                           conj(a.first())*b.second()) ;
}
#endif // defined(DICKSON_USE_LEFT_HANDED_MULT)

#if !defined(DICKSON_USE_LEFT_HANDED_MULT)
//
// multiplication with "i" on the right
//
template< typename subalg >
algebra<subalg>
operator*( const algebra<subalg> & a,
           const algebra<subalg> & b )
{
    return algebra<subalg>(a.first()*b.first()
                           -
                           conj(b.second())*a.second(),
                           
                           b.second()*a.first()
                           +
                           a.second()*conj(b.first())) ;
}
#endif // !defined(DICKSON_USE_LEFT_HANDED_MULT)

/////////////////////////////////////////////////////////
// operator/
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
operator/( const algebra<subalg> & a,
           typename traits<subalg>::scalar b )
{
    algebra<subalg> work = a ;
    work /= b ;
    return work ;
}

/////////////////////////////////////////////////////////
// operator==
/////////////////////////////////////////////////////////

template< typename subalg >
bool
operator==( const algebra<subalg> & a,
            const algebra<subalg> & b )
{
    return
        a.first() == b.first() &&
        a.second() == b.second() ;
}

template< typename subalg >
bool
operator==( typename algebra<subalg>::scalar a,
            const algebra<subalg> & b )
{
    return algebra<subalg>(a) == b ;
}

template< typename subalg >
bool
operator==( const algebra<subalg> & a,
            typename algebra<subalg>::scalar b )
{
    return a == algebra<subalg>(b) ;
}

/////////////////////////////////////////////////////////
// operator!=
/////////////////////////////////////////////////////////

template< typename subalg >
inline
bool
operator!=( const algebra<subalg> & a,
            const algebra<subalg> & b )
{
    return !(a == b) ;
}

template< typename subalg >
inline
bool
operator!=( typename algebra<subalg>::scalar a,
            const algebra<subalg> & b )
{
    return !(a == b) ;
}

template< typename subalg >
inline
bool
operator!=( const algebra<subalg> & a,
            typename algebra<subalg>::scalar b )
{
    return !(a == b) ;
}

/////////////////////////////////////////////////////////
// realpart
/////////////////////////////////////////////////////////

template< typename subalg >
inline
typename traits<subalg>::scalar
realpart( const algebra<subalg> & a )
{
    return a.realpart() ;
}

/////////////////////////////////////////////////////////
// purepart
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
purepart( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;

    algebra<subalg> work = a ;
    work.realpart() = scalar(0) ;
    return work ;
}

/////////////////////////////////////////////////////////
// conj
/////////////////////////////////////////////////////////

//
// iterative version
//
template< typename subalg >
algebra<subalg>
conj( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;

    algebra<subalg> work = a ;

    for( unsigned int i = 1 ; i != algebra<subalg>::dimension ; ++i )
    {
        work[i] *= scalar(-1) ;
    }

    return work ;
}

#if 0
//
// metaprogrammed version --- your compiler's maximum template depth
// must be >= the dimension of the algebra.
//
template< typename subalg >
algebra<subalg>
conj( const algebra<subalg> & a )
{
    algebra<subalg> work = a ;

    dickson_private::
    conj_helper
    <
        algebra<subalg>,
        traits< algebra<subalg> >::dimension - 1
    >::func(work) ;
    
    return work ;
}
#endif

#if 0
//
// multiplying version
//
template< typename subalg >
algebra<subalg>
conj( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;

    algebra<subalg> work = a ;
    work.realpart() *= scalar(-1) ;
    work *= scalar(-1) ;
    return work ;
}
#endif

#if 0
//
// elegant version
//
template< typename subalg >
inline
algebra<subalg>
conj( const algebra<subalg> & a )
{
    return algebra<subalg>(conj(a.first()), -a.second()) ;
}
#endif

/////////////////////////////////////////////////////////
// norm
/////////////////////////////////////////////////////////

template< typename subalg >
inline
typename traits<subalg>::scalar
norm( const algebra<subalg> & a )
{
    return norm(a.first()) + norm(a.second()) ;
}

/////////////////////////////////////////////////////////
// abs
/////////////////////////////////////////////////////////

template< typename subalg >
inline
typename traits<subalg>::scalar
abs( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;
    
    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        abs(a) ;
}

/////////////////////////////////////////////////////////
// inv
/////////////////////////////////////////////////////////

template< typename subalg >
inline
algebra<subalg>
inv( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;

    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        inv(a) ;
}

/////////////////////////////////////////////////////////
// dot
/////////////////////////////////////////////////////////

template< typename subalg >
inline
typename traits<subalg>::scalar
dot( const algebra<subalg> & a,
     const algebra<subalg> & b )
{
    return dot(a.first(), b.first()) + dot(a.second(), b.second()) ;
}

template< typename subalg >
typename traits<subalg>::scalar
dot( typename traits<subalg>::scalar a,
     const algebra<subalg> & b )
{
    return a*realpart(b) ;
}

template< typename subalg >
typename traits<subalg>::scalar
dot( const algebra<subalg> & a,
     typename traits<subalg>::scalar b )
{
    return realpart(a)*b ;
}

/////////////////////////////////////////////////////////
// cross
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
cross( const algebra<subalg> & a,
       const algebra<subalg> & b )
{
    typedef typename traits<subalg>::scalar scalar ;

    return (b*conj(a) - a*conj(b))/scalar(2) ;
}

template< typename subalg >
algebra<subalg>
cross( typename traits<subalg>::scalar a,
       const algebra<subalg> & b )
{
    return a*purepart(b) ;
}

template< typename subalg >
algebra<subalg>
cross( const algebra<subalg> & a,
       typename traits<subalg>::scalar b )
{
    return purepart(a)*(-b) ;
}

/////////////////////////////////////////////////////////
// exp
/////////////////////////////////////////////////////////

template< typename subalg >
inline
algebra<subalg>
exp( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;
    
    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        exp(a) ;
}

/////////////////////////////////////////////////////////
// natural log
/////////////////////////////////////////////////////////

template< typename subalg >
inline
algebra<subalg>
log( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;
    
    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        log(a) ;
}

/////////////////////////////////////////////////////////
// cos
/////////////////////////////////////////////////////////

template< typename subalg >
inline
algebra<subalg>
cos( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;
    
    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        cos(a) ;
}

/////////////////////////////////////////////////////////
// sin
/////////////////////////////////////////////////////////

template< typename subalg >
inline
algebra<subalg>
sin( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;
    
    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        sin(a) ;
}

/////////////////////////////////////////////////////////
// tan
/////////////////////////////////////////////////////////

template< typename subalg >
inline
algebra<subalg>
tan( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;
    
    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        tan(a) ;
}

/////////////////////////////////////////////////////////
// cosh
/////////////////////////////////////////////////////////

template< typename subalg >
inline
algebra<subalg>
cosh( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;
    
    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        cosh(a) ;
}

/////////////////////////////////////////////////////////
// sinh
/////////////////////////////////////////////////////////

template< typename subalg >
inline
algebra<subalg>
sinh( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;
    
    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        sinh(a) ;
}

/////////////////////////////////////////////////////////
// tanh
/////////////////////////////////////////////////////////

template< typename subalg >
inline
algebra<subalg>
tanh( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;
    
    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        tanh(a) ;
}

/////////////////////////////////////////////////////////
// pow
/////////////////////////////////////////////////////////

template< typename subalg >
algebra<subalg>
pow( const algebra<subalg> & a, unsigned int n )
{
    typedef typename traits<subalg>::scalar scalar ;

    algebra<subalg> b = (n % 2 == 0) ? algebra<subalg>(scalar(1)) : a ;
    algebra<subalg> c = a ;
    
    while( n >>= 1 )
    {
        c = c*c ;
        if( n % 2 == 1 )
        {
            b = b*c ;
        }
    }

    return b ;
}

template< typename subalg >
inline
algebra<subalg>
pow( const algebra<subalg> & a, typename traits<subalg>::scalar gamma )
{
    typedef typename traits<subalg>::scalar scalar ;

    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        pow(a, gamma) ;
}

template< typename subalg >
inline
algebra<subalg>
pow( typename traits<subalg>::scalar gamma, const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;

    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        pow(gamma, a) ;
}

/////////////////////////////////////////////////////////
// sqrt
/////////////////////////////////////////////////////////

template< typename subalg >
inline
algebra<subalg>
sqrt( const algebra<subalg> & a )
{
    typedef typename traits<subalg>::scalar scalar ;

    return
        dickson_private::
        select_function<subalg, std::numeric_limits<scalar>::is_integer>::
        sqrt(a) ;
}

/////////////////////////////////////////////////////////
// streams
////////////////////////////////////////////////////////

template< typename subalg, typename T_Char, typename T_Traits >
std::basic_ostream<T_Char, T_Traits> &
operator<<( std::basic_ostream<T_Char, T_Traits> & out,
            const algebra<subalg> & a )
{
    out << '(' ;

    for( unsigned int i = 0 ; i != algebra<subalg>::dimension - 1 ; ++i )
    {
        out << a[i] << ',' ;
    }
    out << a[algebra<subalg>::dimension - 1] ;
    
    out << ')' ;

    return out ;
}

template< typename subalg, typename T_Char, typename T_Traits >
std::basic_istream<T_Char, T_Traits> &
operator>>( std::basic_istream<T_Char, T_Traits> & in,
            algebra<subalg> & a )
{
    T_Char ch ;
    
    ch = 0 ;
    in >> ch ;
    if( ch != '(' )
    {
        in.setstate(std::ios_base::failbit) ;
        return in ;
    }

    algebra<subalg> work ;

    for( unsigned int i = 0 ; i != algebra<subalg>::dimension - 1 ; ++i )
    {
        ch = 0 ;
        in >> work[i] >> ch ;
        if( ch != ',' )
        {
            in.setstate(std::ios_base::failbit) ;
            return in ;
        }
    }
    in >> work[algebra<subalg>::dimension - 1] ;

    ch = 0 ;
    in >> ch ;
    if( ch != ')' )
    {
        in.setstate(std::ios_base::failbit) ;
        return in ;
    }

    a = work ;

    return in ;
}

/////////////////////////////////////////////////////////////////////////
//
// scalar_func
//
/////////////////////////////////////////////////////////////////////////

namespace scalar_func {

template< typename scalar >
inline
scalar
realpart( scalar a )
{
    return a ;
}

template< typename scalar >
inline
scalar
purepart( scalar a )
{
    return scalar(0) ;
}

template< typename scalar >
inline
scalar
conj( scalar a )
{
    return a ;
}

template< typename scalar >
inline
scalar
norm( scalar a )
{
    return a*a ;
}

template< typename scalar >
inline
scalar
inv( const scalar a )
{
    return
        dickson_private::
        select_function<scalar, std::numeric_limits<scalar>::is_integer>::
        inv(a) ;
}

template< typename scalar >
inline
scalar
dot( scalar a, scalar b )
{
    return a*b ;
}

template< typename scalar >
inline
scalar
cross( scalar a, scalar b )
{
    return scalar(0) ;
}

template< typename scalar >
inline
scalar
pow( scalar a, scalar n )
{
    return
        dickson_private::
        select_function<scalar, std::numeric_limits<scalar>::is_integer>::
        pow(a, n) ;
}

template< typename scalar >
scalar
pow( scalar a, unsigned int n )
{
    scalar b = (n % 2 == 0) ? scalar(1) : a ;
    scalar c = a ;
    
    while( n >>= 1 )
    {
        c = c*c ;
        if( n % 2 == 1 )
        {
            b = b*c ;
        }
    }

    return b ;
}

} // namespace scalar_func

} // namespace dickson

#endif

