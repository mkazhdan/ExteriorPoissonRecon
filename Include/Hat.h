/*
Copyright (c) 2023, Michael Kazhdan, Sing-Chun Lee, Marc Alexa, and Maximilian Kohlbrenner
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef HAT_INCLUDED
#define HAT_INCLUDED

#include <unordered_map>
#include <functional>
#include <Eigen/Sparse>
#include <Misha/Array.h>
#include <Misha/Geometry.h>
#include <Misha/Miscellany.h>
#include <Misha/Algebra.h>

namespace Hat
{
	// A SkewSymmetricMatrix< Real , Dim > is a type used to represent the alternating product of two Dim-dimensional vectors.
	// It supports being cast to a SquareMatrix< Real , Dim >.
	template< typename Real , unsigned int _Dim >
	struct SkewSymmetricMatrix : public InnerProductSpace< Real , SkewSymmetricMatrix< Real , _Dim > >
	{
		static const unsigned int Dim = ( _Dim * (_Dim-1) ) / 2;

		SkewSymmetricMatrix( void ){}
		SkewSymmetricMatrix( const SquareMatrix< Real , _Dim > m )
		{
			for( unsigned int idx=0 , i=0 ; i<_Dim ; i++ ) for( unsigned int j=0 ; j<i ; j++ ) _p[idx++] = ( m(i,j) - m(j,i) ) / 2.;
		}
		SquareMatrix< Real , _Dim > operator()( void ) const
		{
			SquareMatrix< Real , _Dim > m;
			for( unsigned int idx=0 , i=0 ; i<_Dim ; i++ ) for( unsigned int j=0 ; j<i ; j++ ) m(i,j) = _p[idx++];
			return m - m.transpose();
		}
		Real  operator[]( unsigned int idx ) const { return _p[idx]; }
		Real &operator[]( unsigned int idx ){ return _p[idx]; }

		void Add( const SkewSymmetricMatrix &skew ){ _p += skew._p; }
		void Scale( Real s ){ _p *= s; }
		Real InnerProduct( const SkewSymmetricMatrix &skew ) const { return Point< Real , Dim >::Dot( _p , skew._p ); }

	protected:
		Point< Real , Dim > _p;
	};

	template< unsigned int Dim > SkewSymmetricMatrix< double , Dim > Wedge( Point< double , Dim > v1 , Point< double , Dim > v2 );

	template< unsigned int Dim >
	struct Index : public Point< int , Dim >
	{
		Index( void ){}
		Index( Point< int , Dim > p ) : Point< int , Dim >( p ){}
		template< typename ... Indices > static Index Min( Index i , Indices ... is );
		template< typename ... Indices > static Index Max( Index i , Indices ... is );
		bool operator == ( Index i ) const;
		bool operator != ( Index i ) const;
		bool operator <  ( Index i ) const;
		bool operator <= ( Index i ) const;
		bool operator >  ( Index i ) const;
		bool operator >= ( Index i ) const;
	};

	template< unsigned int Dim >
	struct Range : public std::pair< Index< Dim > , Index< Dim > >
	{
		using std::pair< Index< Dim > , Index< Dim > >::first;
		using std::pair< Index< Dim > , Index< Dim > >::second;
		template< typename ... Ranges > static Range Intersect( Ranges ... rs );
		bool empty( void ) const;
		bool contains( Index< Dim > i ) const;
		size_t size( void ) const;

		template< typename IndexFunctor /* = std::function< void (Index<Dim>) > */ > void process( IndexFunctor f ) const;
	};


	// The base case (when Dim=1)
	struct _Basis
	{
		// The range of elements in the support of the i-th function
		static std::pair< int , int > FunctionSupport( int i );
		static std::pair< int , int > FunctionSupport( std::pair< int , int > r );

		// The range of functions supported on the i-th element
		static std::pair< int , int >  ElementSupport( int i );
		static std::pair< int , int >  ElementSupport( std::pair< int , int > r );

		// The integral of x^p * (1-x)^q over the interval
		static double Integral( unsigned int p , unsigned int q );
	};


	// The representation of the restriction of the product of hat functions and derivatives to an element
	// Setting P(x) = x , Q(x) = (1-x) , dP(x) = 1 , dQ(x) = -1, this represents the function:
	//		P^p * Q^q * dP^dp * dQ^dq * coefficient
	template< unsigned int Dim >
	struct ElementFunction
	{
		ElementFunction( void ) : _coefficient(0) {}
		ElementFunction( Index< Dim > e , Index< Dim > i );

		ElementFunction d( unsigned int d ) const;

		ElementFunction &operator *= ( const ElementFunction &f );
		ElementFunction &operator *= ( double s );
		ElementFunction  operator *  ( const ElementFunction &f ) const { ElementFunction _f = (*this) ; _f *= f ; return _f; }
		ElementFunction  operator *  ( double s ) const { ElementFunction f = (*this) ; f *= s ; return f; }
		ElementFunction  operator - ( void ) const { return (*this) * (-1.); }

		double integral( unsigned int res ) const;
	protected:
		double _coefficient;
		Index< Dim > _p , _q , _d;
	};
	template< unsigned int Dim > ElementFunction< Dim > operator * ( double s , const ElementFunction< Dim > &f ){ return f*s; }

	template< unsigned int Dim >
	struct ElementVector
	{
		ElementVector( void ){}
		ElementVector( const ElementFunction< Dim > &f );
		ElementFunction< Dim > &operator[]( unsigned int d ){ return _c[d]; }
		const ElementFunction< Dim > &operator[]( unsigned int d ) const { return _c[d]; }
		ElementVector &operator *= ( const ElementFunction< Dim > &f ){ for( unsigned int d=0 ; d<Dim ; d++ ) _c[d] *= f ; return *this; }
		ElementVector &operator *= ( double s ){ for( unsigned int d=0 ; d<Dim ; d++ ) _c[d] *= s ; return *this; }
		ElementVector  operator *  ( const ElementFunction< Dim > &f ) const { ElementVector v = (*this) ; v *= f ; return v; }
		ElementVector  operator *  ( double s ) const { ElementVector v = (*this) ; v *= s ; return v; }
		ElementVector  operator - ( void ) const { return (*this) * (-1.); }
	protected:
		ElementFunction< Dim > _c[Dim];
	};
	template< unsigned int Dim > ElementVector< Dim > operator * ( double s , const ElementVector< Dim > &v ){ return v*s; }
	template< unsigned int Dim > ElementVector< Dim > operator * ( const ElementFunction< Dim > f , const ElementVector< Dim > &v ){ return v*f; }

	template< unsigned int _Dim >
	struct ElementWedge
	{
		static const unsigned int Dim = ( _Dim * ( _Dim-1 ) ) / 2;
		ElementWedge( void ){}
		ElementWedge( const ElementVector< _Dim > &v1 , const ElementVector< _Dim > &v2 );
		ElementFunction< _Dim > &operator()( unsigned int i , unsigned int d ){ return _c[i][d]; }
		const ElementFunction< _Dim > &operator()( unsigned int i , unsigned int d ) const { return _c[i][d]; }
		ElementWedge &operator *= ( const ElementFunction< _Dim > &f ){ for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) _c[i][d] *= f ; return *this; }
		ElementWedge &operator *= ( double s ){ for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) _c[i][d] *= s ; return *this; }
		ElementWedge  operator *  ( const ElementFunction< _Dim > &f ) const { ElementWedge w = (*this) ; w *= f ; return w; }
		ElementWedge  operator *  ( double s ) const { ElementWedge w = (*this) ; w *= s ; return w; }
		ElementWedge  operator - ( void ) const { return (*this) * (-1.); }
	protected:
		ElementFunction< _Dim > _c[2][Dim];
	};

	template< unsigned int Dim > ElementWedge< Dim > operator * ( double s , const ElementWedge< Dim > &w ){ return w*s; }
	template< unsigned int Dim > ElementWedge< Dim > operator * ( const ElementFunction< Dim > &f , const ElementWedge< Dim > &w ){ return w*f; }
	template< unsigned int Dim > ElementWedge< Dim > operator ^ ( const ElementVector< Dim > &v1 , const ElementVector< Dim > &v2 ){ return ElementWedge< Dim >( v1 , v2 ); }


	template< unsigned int Dim > struct ScalarFunctions;
	template< unsigned int Dim > struct WedgeFunctions;


	template< unsigned int Dim >
	struct Basis
	{
		// The range of elements in the support of a hat basis function
		static Range< Dim > FunctionSupport( Index< Dim > i );
		static Range< Dim > FunctionSupport( Range< Dim > r );

		// The range of hat basis functions supported on an element
		static Range< Dim > ElementSupport( Index< Dim > i );
		static Range< Dim > ElementSupport( Range< Dim > r );

		// The integral of the product of two hat functions over an element
		static double ElementMass( unsigned int res , Index< Dim > e , Index< Dim > i , Index< Dim > j );

		// The integral of the dot-product of the gradients of two hat functions over an element
		static double ElementStiffness( unsigned int res , Index< Dim > e , Index< Dim > i , Index< Dim > j );

		// The integral of the dot-product of the wedge-product of the gradients of two pairs of hat functions over an element
		static double ElementGradientWedgeMass( unsigned int res , Index< Dim > e , std::pair< Index< Dim > , Index< Dim > > i , std::pair< Index< Dim > , Index< Dim > > j );

		// Estimates the integral over the element
		template< typename F /* = std::function< double ( Point< double , Dim > ) */ >
		static double Integral( unsigned int res , Index< Dim > e , F f , unsigned samplingRes );

		// Get the function associated with an index
		static std::function< double ( Point< double , Dim > ) > Scalar( unsigned int res , Index< Dim > i );

		// Get the gradient of the function associated with an index
		static std::function< Point< double , Dim > ( Point< double , Dim > ) > Gradient( unsigned int res , Index< Dim > i );
	protected:
		static int _RelativeIndex( Index< Dim > e , Index< Dim > i );
		friend struct ScalarFunctions< Dim >;
		friend struct WedgeFunctions< Dim >;
	};


	template< unsigned int Dim >
	struct ScalarFunctions
	{
		struct Stencil
		{
			friend ScalarFunctions;
			Stencil( unsigned int res );
			Stencil( const Stencil &stencil );
			Stencil( Stencil &&stencil );
			~Stencil( void );

			Stencil &operator = ( const Stencil &stencil );

			double operator()( Index< Dim > e , Index< Dim > f1 , Index< Dim > f2 ) const;
			double operator()( Index< Dim > f1 , Index< Dim > f2 ) const;

			unsigned int resolution( void ) const;
			Range< Dim > eRange( void ) const;
		protected:
			double **_values;
			unsigned int _res;
			Range< Dim > _eRange;
		};

		struct FullStencil
		{
			static constexpr unsigned int StencilNum( void ); // 3^Dim
			using Entry = std::pair< Index< Dim > , double >;
			using Row = std::vector< Entry >;

			FullStencil( const Stencil &stencil );
			const Row &row( Index< Dim > f ) const;
			static unsigned int StencilIndex( Index< Dim > f , unsigned int res );

		protected:
			unsigned int _res;
			std::vector< Row > _rows;
		};

		static Stencil      MassStencil( unsigned int r );
		static Stencil StiffnessStencil( unsigned int r );
		double dot( const Stencil &stencil , const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;
		double dot( const FullStencil &stencil , const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;

		ScalarFunctions( unsigned int resolution ) : _r(resolution) { for( unsigned int d=0 ; d<Dim ; d++ ) _fRange.first[d] = _eRange.first[d] = 0 , _fRange.second[d] = _r+1 , _eRange.second[d] = _r; }
		unsigned int resolution( void ) const { return _r; }
		size_t functionNum( void ) const;
		size_t elementNum( void ) const;

		Range< Dim > fRange( void ) const { return _fRange; }
		Range< Dim > eRange( void ) const { return _eRange; }

		size_t elementIndex( Index< Dim > i ) const;
		Index< Dim > elementIndex( size_t i ) const;
		size_t functionIndex( Index< Dim > i ) const;
		Index< Dim > functionIndex( size_t i ) const;

		Eigen::SparseMatrix< double > mass( void ) const;
		Eigen::SparseMatrix< double > stiffness( void ) const;
		Eigen::SparseMatrix< char > incidence( void ) const;
		template< typename SampleFunctor /* = std::function< Point< double , Dim > ( unsigned int idx ) > */ >
		Eigen::SparseMatrix< double > evaluation( SampleFunctor F , size_t sampleNum ) const;
		Eigen::SparseMatrix< double > evaluation( const std::vector< Point< double , Dim > > &samples ) const;
		template< typename SampleFunctor /* = std::function< Point< double , Dim > ( unsigned int idx ) > */ >
		Eigen::SparseMatrix< double > dEvaluation( SampleFunctor F , size_t sampleNum ) const;
		Eigen::SparseMatrix< double > dEvaluation( const std::vector< Point< double , Dim > > &samples ) const;
		Eigen::SparseMatrix< double > prolongation( void ) const;

		double value( const Eigen::VectorXd &x , Point< double , Dim > p ) const;
		Point< double , Dim > gradient( const Eigen::VectorXd &x , Point< double , Dim > p ) const;
	protected:
		Eigen::SparseMatrix< double > _systemMatrix( Stencil stencil ) const;
		unsigned int _r;
		Range< Dim > _fRange , _eRange;
	};

	template< unsigned int Dim >
	struct WedgeFunctions
	{
		struct Stencil
		{
			friend WedgeFunctions;
			Stencil( unsigned int res );
			Stencil( const Stencil &stencil );
			Stencil( Stencil &&stencil );
			~Stencil( void );

			Stencil &operator = ( const Stencil &stencil );

			double operator()( Index< Dim > e , Index< Dim > f1 , Index< Dim > f2 , Index< Dim > g1 , Index< Dim > g2 ) const;
			double operator()( Index< Dim > f1 , Index< Dim > f2 , Index< Dim > g1 , Index< Dim > g2 ) const;

			unsigned int resolution( void ) const;
			Range< Dim > eRange( void ) const;
		protected:
			double ****_values;
			unsigned int _res;
			Range< Dim > _eRange;
		};

		struct FullStencil
		{
			struct Entry
			{
				Index< Dim > g1 , g2;
				double value;
				Entry( Index< Dim > g1 , Index< Dim > g2 , double value ) : g1(g1) , g2(g2) , value(value) {}
			};
			struct Row
			{
				Index< Dim > f2;
				std::vector< Entry > entries;
			};

			FullStencil( const Stencil &stencil );
			const std::vector< Row > &rows( Index< Dim > f1 ) const;

		protected:
			unsigned int _res;
			std::vector< std::vector< Row > > _rows;
			unsigned int _index( Index< Dim > f1 ) const;
		};

		static Stencil MassStencil( unsigned int r );
		double dot( const Stencil &stencil , const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;
		double dot( const FullStencil &stencil , const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;

		WedgeFunctions( unsigned int resolution );
		size_t resolution( void ) const { return _r; }
		size_t functionNum( void ) const { return _indices.size(); }

		size_t index( std::pair< size_t , size_t > idx , bool &flip ) const;
		std::pair< size_t , size_t > index( size_t idx ) const;
		bool setIndex( std::pair< size_t , size_t > idx , size_t &i , bool &flip ) const;

		Eigen::VectorXd wedge( const Eigen::VectorXd &f1 , const Eigen::VectorXd &f2 ) const;
		// Returns the function corresponding to g -> df ^ dg
		Eigen::SparseMatrix< double > wedge( const Eigen::VectorXd &f ) const;
		Eigen::SparseMatrix< double > mass( void ) const;
		Eigen::SparseMatrix< double > prolongation( void ) const;

		Eigen::VectorXd dualVector( ConstPointer( SkewSymmetricMatrix< double , Dim > ) skew ) const;
		Eigen::SparseMatrix< double > dualMatrix( ConstPointer( SkewSymmetricMatrix< double , Dim > ) skew ) const;

	protected:
		struct _Hasher{ size_t operator()( const std::pair< size_t , size_t > &i ) const { return i.first * i.second; } };

		unsigned int _r;

		std::unordered_map< std::pair< size_t , size_t > , size_t , _Hasher > _indexMap;
		std::vector< std::pair< size_t , size_t > > _indices;
	};

	// An iterator allowing for triplets to be set by multiple threads, without write-on-write conflicts, but still behaving as a single container
	template< typename T >
	struct OMPTriplets
	{
		struct iterator
		{
			unsigned int thread , index;
			iterator &operator++( void ){ if( index+1==_ompTriplets[thread].size() ) _nextThread() , index=0 ; else index++ ; return *this; }
			iterator  operator++( int ){ iterator iter = *this ; if( index+1==_ompTriplets[thread].size() ) _nextThread() , index=0 ; else index++ ; return iter; }
			bool operator == ( iterator iter ) const { return thread==iter.thread && index==iter.index; };
			bool operator != ( iterator iter ) const { return thread!=iter.thread || index!=iter.index; };
			Eigen::Triplet< T > &operator * ( void ){ return _ompTriplets[thread][index]; }
			Eigen::Triplet< T > *operator -> ( void ){ return &( _ompTriplets[thread][index] ); }
			const Eigen::Triplet< T > &operator *  ( void ) const { return _ompTriplets[thread][index]; }
			const Eigen::Triplet< T > *operator -> ( void ) const { return &( _ompTriplets[thread][index] ); }
		protected:
			OMPTriplets &_ompTriplets;

			iterator( OMPTriplets &ompTriplets , bool begin ) : _ompTriplets( ompTriplets ) , index(0)
			{
				if( begin )
				{
					thread = -1;
					_nextThread();
				}
				else thread = (unsigned int)_ompTriplets.size();
			}

			void _nextThread( void )
			{
				do{ thread++; }
				while( thread<_ompTriplets.size() && _ompTriplets[thread].size()==0 );
			}

			friend OMPTriplets;
		};

		struct const_iterator
		{
			unsigned int thread , index;
			const_iterator &operator++( void ){ if( index+1==_ompTriplets[thread].size() ) _nextThread() , index=0 ; else index++ ; return *this; }
			const_iterator  operator++( int ){ const_iterator iter = *this ; if( index+1==_ompTriplets[thread].size() ) _nextThread() , index=0 ; else index++ ; return iter; }
			bool operator == ( const_iterator iter ) const { return thread==iter.thread && index==iter.index; };
			bool operator != ( const_iterator iter ) const { return thread!=iter.thread || index!=iter.index; };
			const Eigen::Triplet< T > &operator *  ( void ) const { return _ompTriplets[thread][index]; }
			const Eigen::Triplet< T > *operator -> ( void ) const { return &( _ompTriplets[thread][index] ); }
		protected:
			const OMPTriplets &_ompTriplets;

			const_iterator( const OMPTriplets &ompTriplets , bool begin ) : _ompTriplets( ompTriplets ) , index(0)
			{
				if( begin )
				{
					thread = -1;
					_nextThread();
				}
				else thread = (unsigned int)_ompTriplets.size();
			}

			void _nextThread( void )
			{
				do{ thread++; }
				while( thread<_ompTriplets.size() && _ompTriplets[thread].size()==0 );
			}

			friend OMPTriplets;
		};

		OMPTriplets( void ) : _triplets( omp_get_num_procs() ){}
		std::vector< std::vector< Eigen::Triplet< T > > > _triplets;
		std::vector< Eigen::Triplet< T > > &operator[] ( unsigned int thread ){ return _triplets[thread]; }
		const std::vector< Eigen::Triplet< T > > &operator[] ( unsigned int thread ) const { return _triplets[thread]; }
		void clear( void ) 
		{
			for( unsigned int i=0 ; i<_triplets.size() ; i++ )
			{
				std::vector< Eigen::Triplet< T > > t;
				std::swap( _triplets[i] , t );
			}
		}
		size_t size( void ) const { return _triplets.size(); }
		iterator begin( void ){ return iterator( *this , true ); }
		iterator end( void ){ return iterator( *this , false ); }
	};
#include "Hat.inl"
}
#endif // HAT_INCLUDED