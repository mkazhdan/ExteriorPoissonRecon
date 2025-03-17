/*
Copyright (c) 2025, Michael Kazhdan, Sing-Chun Lee, Marc Alexa, and Maximilian Kohlbrenner
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

///////////////////
// SquareMatrix //
///////////////////
template< typename Real , unsigned int Dim , bool Sym >
SquareMatrix< Real , Dim , Sym >::SquareMatrix( Point< Real , Dim > v1 , Point< Real , Dim > v2 )
{
	if( Sym ) for( unsigned int i=0 , idx=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<=i ; j++ , idx++ ) _p[idx] = ( v1[j]*v2[i] + v1[i]*v2[j] ) / 2.;
	else      for( unsigned int i=0 , idx=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j< i ; j++ , idx++ ) _p[idx] = ( v1[j]*v2[i] - v1[i]*v2[j] ) / 2.;
}

template< typename Real , unsigned int Dim , bool Sym >
SquareMatrix< Real , Dim , Sym >::SquareMatrix( const MishaK::SquareMatrix< Real , Dim > m )
{
	unsigned int idx=0;
	if( Sym ) for( unsigned int i=0 , idx=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<=i ; j++ , idx++ ) _p[idx] = ( m(i,j) + m(j,i) ) / 2.;
	else      for( unsigned int i=0 , idx=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j< i ; j++ , idx++ ) _p[idx] = ( m(i,j) - m(j,i) ) / 2.;
}

template< typename Real , unsigned int Dim , bool Sym >
MishaK::SquareMatrix< Real , Dim > SquareMatrix< Real , Dim , Sym >::operator()( void ) const
{
	MishaK::SquareMatrix< Real , Dim > m;
	if( Sym ) for( unsigned int i=0 , idx=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<=i ; j++ , idx++ ) m(i,j) = _p[idx] , m(j,i) =  _p[idx];
	else      for( unsigned int i=0 , idx=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j< i ; j++ , idx++ ) m(i,j) = _p[idx] , m(j,i) = -_p[idx];
	return m;
}

////////////////////////////
// SquareMatrix::UnitTest //
////////////////////////////
template< typename Real , unsigned int Dim , bool Sym >
void SquareMatrix< Real , Dim , Sym >::UnitTest::operator()( void ) const
{
	outerProduct();
	innerProduct();
}

template< typename Real , unsigned int Dim , bool Sym >
void SquareMatrix< Real , Dim , Sym >::UnitTest::outerProduct( void ) const
{
	double err1 = 0 , err2 = 0;
	for( unsigned int n=0 ; n<_N ; n++ )
	{
		// Check that the notions of inner-product are consistent
		{
			Point< double , Dim > p1 = RandomBallPoint< Real , Dim >() , p2 = RandomBallPoint< Real , Dim >();
			MishaK::SquareMatrix< double , Dim > m = OuterProduct(p1,p2);
			MishaK::SquareMatrix< double , Dim > v1 = SquareMatrix< double , Dim , Sym >( p1 , p2 )();
			MishaK::SquareMatrix< double , Dim > v2 = SquareMatrix< double , Dim , Sym >( m )();
			MishaK::SquareMatrix< double , Dim > v3 = Sym ? ( m + m.transpose() ) / 2 : ( m - m.transpose() ) / 2;
			err1 += (v1-v2).squareNorm();
			err2 += (v2-v3).squareNorm();
		}
	}
	std::cout << "\tOuter product error: " << err1 << " / " << err2 << std::endl;
}

template< typename Real , unsigned int Dim , bool Sym >
void SquareMatrix< Real , Dim , Sym >::UnitTest::innerProduct( void ) const
{
	double err1 = 0 , err2 = 0;
	for( unsigned int n=0 ; n<_N ; n++ )
	{
		// Check that the inner-products are consistent
		{
			SquareMatrix< double , Dim , Sym > m , n;
			for( unsigned int d=0 ; d<SquareMatrix< double , Dim , Sym >::Coefficients ; d++ ) m[d] = Random< double >()*2-1 , n[d] = Random< double >()*2-1;
			double v1 = SquareMatrix< double , Dim , Sym >::Dot( m , n );
			double v2 = MishaK::SquareMatrix< double , Dim >::Dot( m() , n() );
			err1 += (v1-v2) * (v1-v2);
		}
		{
			MishaK::SquareMatrix< double , Dim > m , n;
			for( unsigned int i=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) m(i,j) = Random< double >()*2 - 1 , n(i,j) = Random< double >()*2 - 1;
			if constexpr( Sym ) m = ( m + m.transpose() )/2 , n = ( n + n.transpose() )/2;
			else                m = ( m - m.transpose() )/2 , n = ( n - n.transpose() )/2;
			double v1 = SquareMatrix< double , Dim , Sym >::Dot( SquareMatrix< double , Dim , Sym >(m) , SquareMatrix< double , Dim , Sym >(n) );
			double v2 = MishaK::SquareMatrix< double , Dim >::Dot(m,n);
			err2 += (v1-v2) * (v1-v2);
		}
	}
	std::cout << "\tInner product error: " << err1 << " / " << err2 << std::endl;
}

////////////
// _Basis //
////////////
std::pair< int , int > _Basis::FunctionSupport( int i ){ return std::pair< int , int >( i-1 , i+1 ); }
std::pair< int , int > _Basis::FunctionSupport( std::pair< int , int > r ){ return std::pair< int , int >( r.first-1 , r.second+0 ); }
std::pair< int , int > _Basis:: ElementSupport( int i ){ return std::pair< int , int >( i+0 , i+2 ); }
std::pair< int , int > _Basis:: ElementSupport( std::pair< int , int > r ){ return std::pair< int , int >( r.first+0 , r.second+1 ); }

double _Basis::Integral( unsigned int p , unsigned int q )
{
	auto Prod =[]( unsigned int first , unsigned int last )
	{
		double f=1;
		for( unsigned int i=first ; i<=last ; i++ ) f *= i;
		return f;
	};
	if( p>q ) return Prod(1,q)/Prod(p+1,p+q+1);
	else      return Prod(1,p)/Prod(q+1,p+q+1);
}

/////////////////////
// ElementFunction //
/////////////////////
template< unsigned int Dim >
ElementFunction< Dim >::ElementFunction( Index< Dim > e , Index< Dim > i )
{
	if( Basis< Dim >::FunctionSupport( i ).contains( e ) )
	{
		_coefficient = 1;
		for( unsigned int d=0 ; d<Dim ; d++ ) 
			if( e[d]==i[d] ) _q[d] = 1;
			else             _p[d] = 1;
	}
	else _coefficient = 0;
}

template< unsigned int Dim >
ElementFunction< Dim > ElementFunction< Dim >::d( unsigned int d ) const
{
	if( _d[d] ) ERROR_OUT( "Cannot compute second derivative" );
	if( _p[d] && _q[d] ) ERROR_OUT( "Derivative is not an ElementFunction" );

	ElementFunction< Dim > f = (*this);
	if( _p[d] )
	{
		f._coefficient *= f._p[d];
		f._p[d]--;
		f._d[d] = 1;
	}
	else if( _q[d] )
	{
		f._coefficient *= -f._q[d];
		f._q[d]--;
		f._d[d] = 1;
	}
	return f;
}

template< unsigned int Dim >
ElementFunction< Dim > &ElementFunction< Dim >::operator *= ( const ElementFunction &f )
{
	_p += f._p;
	_q += f._q;
	_d += f._d;
	_coefficient *= f._coefficient;
	return *this;
}

template< unsigned int Dim >
ElementFunction< Dim > &ElementFunction< Dim >::operator *= ( double s )
{
	_coefficient *= s;
	return *this;
}

template< unsigned int Dim >
double ElementFunction< Dim >::integral( unsigned int res ) const
{
	double integral = _coefficient;
	for( unsigned int d=0 ; d<Dim ; d++ ) integral *= _Basis::Integral( _p[d] , _q[d] ) * pow( res , _d[d] ) / res;
	return integral;
}

template< unsigned int Dim >
ElementVector< Dim > ElementFunction< Dim >::gradient( void ) const
{
	ElementVector< Dim > grad;
	for( unsigned int d=0 ; d<Dim ; d++ ) grad[d] = this->d(d);
	return grad;
}

template< unsigned int Dim >
double ElementFunction< Dim >::operator()( Point< double , Dim > p ) const
{
	double v = _coefficient;
	for( unsigned int d=0 ; d<Dim ; d++ ) v *= pow( p[d] , _p[d] ) * pow( 1.-p[d] , _q[d] );
	return v;
}

///////////////////
// ElementVector //
///////////////////
template< unsigned int Dim >
Point< double , Dim > ElementVector< Dim >::operator()( Point< double , Dim > p ) const
{
	Point< double , Dim > v;
	for( unsigned int d=0 ; d<Dim ; d++ ) v[d] = _c[d]( p );
	return v;
}

/////////////////////////////////////////
// ElementProduct::IndexToCoefficients //
/////////////////////////////////////////
template< unsigned int Dim , bool Sym >
ElementProduct< Dim , Sym >::IndexToCoefficients::IndexToCoefficients( void )
{
	static std::mutex mut;
	std::lock_guard< std::mutex > guard( mut );
	{
		if( Sym ) for( unsigned int i=0 , idx=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<=i ; j++ , idx++ ) coefficients[idx].first = i , coefficients[idx].second = j;
		else      for( unsigned int i=0 , idx=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j< i ; j++ , idx++ ) coefficients[idx].first = i , coefficients[idx].second = j;
	}
}

////////////////////
// ElementProduct //
////////////////////
template< unsigned int Dim , bool Sym >
typename ElementProduct< Dim , Sym >::IndexToCoefficients ElementProduct< Dim , Sym >::_IndexToCoefficients;

template< unsigned int Dim , bool Sym >
ElementProduct< Dim , Sym >::ElementProduct( const ElementVector< Dim > &v1 , const ElementVector< Dim > &v2 )
{
	if( Sym ) for( unsigned int i=0 , idx=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<=i ; j++ , idx++ ) _c[0][idx] = v1[j]*v2[i]/2 , _c[1][idx] =  v1[i]*v2[j]/2;
	else      for( unsigned int i=0 , idx=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j< i ; j++ , idx++ ) _c[0][idx] = v1[j]*v2[i]/2 , _c[1][idx] = -v1[i]*v2[j]/2;
}

template< unsigned int Dim , bool Sym >
MishaK::SquareMatrix< double , Dim > ElementProduct< Dim , Sym >::operator()( Point< double , Dim > p ) const
{
	MishaK::SquareMatrix< double , Dim > v;
	if( Sym ) for( unsigned int i=0 , idx=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<=i ; j++ , idx++ ) v(i,j) = _c[0][idx](p) + _c[1][idx](p) , v(j,i) =  _c[1][idx](p) + _c[0][idx](p);
	else      for( unsigned int i=0 , idx=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j< i ; j++ , idx++ ) v(i,j) = _c[0][idx](p) + _c[1][idx](p) , v(j,i) = -_c[1][idx](p) - _c[0][idx](p);
	return v;
}

//////////////////////////////
// ElementProduct::UnitTest //
//////////////////////////////
template< unsigned int Dim , bool Sym >
void ElementProduct< Dim , Sym >::UnitTest::operator()( void ) const
{
	evaluation();
}

template< unsigned int Dim , bool Sym >
void ElementProduct< Dim , Sym >::UnitTest::evaluation( void ) const
{
	Index< Dim > E;

	double err = 0;
	auto f = [&]( Index< Dim > F1 , Index< Dim > F2 )
	{
		ElementVector< Dim > G1 = ElementFunction< Dim >( E , F1 ).gradient() , G2 = ElementFunction< Dim >( E , F2 ).gradient();
		ElementProduct< Dim , Sym > P( G1 , G2 );
		for( unsigned int n=0 ; n<_N ; n++ )
		{
			Point< double , Dim > p;
			for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = Random< double >();
			// (Anti)-symmetrize the outer product
			MishaK::SquareMatrix< double , Dim > v1 = SquareMatrix< double , Dim , Sym >( OuterProduct( G1(p) , G2(p) ) )();
			MishaK::SquareMatrix< double , Dim > v2 = P(p);
			err += ( v1 - v2 ).squareNorm();
		}
	};
	Basis< Dim >::ElementSupport( E ).template process< 2 >( f );

	std::cout << "\tEvaluation: " << err << std::endl;
}


///////////
// Basis //
///////////
template< unsigned int Dim >
Range< Dim > Basis< Dim >::FunctionSupport( Index< Dim > f )
{
	Range< Dim > r;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		std::pair< unsigned int , unsigned int > _r = _Basis::FunctionSupport( f[d] );
		r.first[d] = _r.first , r.second[d] = _r.second;
	}
	return r;
}

template< unsigned int Dim >
Range< Dim > Basis< Dim >::ElementSupport( Index< Dim > e )
{
	Range< Dim > r;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		std::pair< unsigned int , unsigned int > _r = _Basis::ElementSupport( e[d] );
		r.first[d] = _r.first , r.second[d] = _r.second;
	}
	return r;
}

template< unsigned int Dim >
Range< Dim > Basis< Dim >::FunctionSupport( Range< Dim > r )
{
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		std::pair< unsigned int , unsigned int > _r = _Basis::FunctionSupport( std::pair< int , int >( r.first[d] , r.second[d] ) );
		r.first[d] = _r.first , r.second[d] = _r.second;
	}
	return r;
}

template< unsigned int Dim >
Range< Dim > Basis< Dim >::ElementSupport( Range< Dim > r )
{
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		std::pair< unsigned int , unsigned int > _r = _Basis::ElementSupport( std::pair< int , int >( r.first[d] , r.second[d] ) );
		r.first[d] = _r.first , r.second[d] = _r.second;
	}
	return r;
}

template< unsigned int Dim >
double Basis< Dim >::ElementMass( unsigned int res , Index< Dim > e , Index< Dim > i , Index< Dim > j )
{
	return ( ElementFunction< Dim >( e , i ) * ElementFunction< Dim >( e , j ) ).integral( res );
}

template< unsigned int Dim >
double Basis< Dim >::ElementStiffness( unsigned int res , Index< Dim > e , Index< Dim > i , Index< Dim > j )
{
	ElementVector< Dim > v1( ElementFunction< Dim >( e , i ) ) , v2( ElementFunction< Dim >( e , j ) );
	double integral = 0;
	for( unsigned int d=0 ; d<Dim ; d++ ) integral += ( v1[d] * v2[d] ).integral( res );
	return integral;
}

template< unsigned int Dim >
template< bool Sym >
double Basis< Dim >::ElementGradientProductMass( unsigned int res , Index< Dim > e , std::pair< Index< Dim > , Index< Dim > > i , std::pair< Index< Dim > , Index< Dim > > j )
{
	double integral = 0;
	ElementFunction< Dim > f1 = ElementFunction< Dim >( e , i.first ) , f2 = ElementFunction< Dim >( e , i.second );
	ElementFunction< Dim > g1 = ElementFunction< Dim >( e , j.first ) , g2 = ElementFunction< Dim >( e , j.second );
	ElementProduct< Dim , Sym > p1( f1.gradient() , f2.gradient() );
	ElementProduct< Dim , Sym > p2( g1.gradient() , g2.gradient() );
	for( unsigned int d=0 ; d<ElementProduct< Dim , Sym >::Coefficients ; d++ ) for( unsigned int ii=0 ; ii<2 ; ii++ ) for( unsigned int jj=0 ; jj<2 ; jj++ )
		integral += ( p1(ii,d) * p2(jj,d) ).integral( res ) * ( ElementProduct< Dim , Sym >::IsDiagonal(d) ? 1 : 2 );
	return integral;
}

template< unsigned int Dim >
template< typename T , typename F /* = std::function< T ( Point< double , Dim > ) */ >
T Basis< Dim >::Integral( unsigned int res , Index< Dim > e , F f , unsigned samplingRes )
{
	T sum{};
	Range< Dim > sRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) sRange.second[d] = samplingRes;

	auto g = [&]( Index< Dim > i )
	{
		Point< double , Dim > p;
		for( unsigned int d=0 ; d<Dim; d++ ) p[d] = e[d] + ( i[d] + 0.5 ) / samplingRes;
		p /= res;
		sum += f(p);
	};
	sRange.process( g );
	for( unsigned int d=0 ; d<Dim ; d++ ) sum /= samplingRes * res;
	return sum;
}

template< unsigned int Dim >
template< typename T , typename F /* = std::function< T ( Point< double , Dim > ) */ >
T Basis< Dim >::Integral( F f , unsigned samplingRes )
{
	T sum{};
	Range< Dim > sRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) sRange.second[d] = samplingRes;

	auto g = [&]( Index< Dim > I )
	{
		Point< double , Dim > p;
		for( unsigned int d=0 ; d<Dim; d++ ) p[d] = ( I[d] + 0.5 ) / samplingRes;
		sum += f(p);
	};
	sRange.process( g );
	for( unsigned int d=0 ; d<Dim ; d++ ) sum /= samplingRes;
	return sum;
}

template< unsigned int Dim >
std::function< double ( Point< double , Dim > ) > Basis< Dim >::Scalar( unsigned int res , Index< Dim > idx )
{
	std::function< double ( Point< double , Dim > ) > s = [res,idx]( Point< double , Dim > p )
	{
		p *= res;
		double value = 1.;
		for( unsigned int d=0 ; d<Dim ; d++ )
		{
			int ip = (int)floor( p[d] );
			if( ip+1==idx[d] ) value *= p[d] - ip;
			else if( ip==idx[d] ) value *= ( 1. - ( p[d] - ip ) );
			else value = 0;
		}
		return value;
	};
	return s;
}

template< unsigned int Dim >
std::function< Point< double , Dim > ( Point< double , Dim > ) > Basis< Dim >::Gradient( unsigned int res , Index< Dim > idx )
{
	std::function< Point< double , Dim > ( Point< double , Dim > ) > g = [res,idx]( Point< double , Dim > p )
	{
		p *= res;
		Point< double , Dim > value;
		for( unsigned int d=0 ; d<Dim ; d++ ) value[d] = 1.;
		for( unsigned int d=0 ; d<Dim ; d++ )
		{
			int ip = (int)floor( p[d] );
			for( unsigned int _d=0 ; _d<Dim ; _d++ )
			{
				if( _d!=d )
				{
					if( ip+1==idx[d] ) value[_d] *= p[d] - ip;
					else if( ip==idx[d] ) value[_d] *= ( 1. - ( p[d] - ip ) );
					else value[_d] = 0;
				}
				else
				{
					if( ip+1==idx[d] ) value[_d] *= res;
					else if( ip==idx[d] ) value[_d] *= -(int)res;
					else value[_d] = 0;
				}
			}
		}
		return value;
	};
	return g;
}

template< unsigned int Dim >
template< unsigned int Radius , unsigned int D >
void Basis< Dim >::_RelativeIndex( Index< Dim > e , Index< Dim > f , unsigned int &idx )
{
	idx = idx * _Width< Radius >() + ( f[D]-e[D]+Radius );
	if constexpr( D+1<Dim ) _RelativeIndex< Radius , D+1 >( e , f , idx );
}

template< unsigned int Dim >
template< unsigned int Radius >
int Basis< Dim >::_RelativeIndex( Index< Dim > e , Index< Dim > f )
{
	unsigned int idx = 0;
	_RelativeIndex< Radius , 0 >( e , f , idx );
	return idx;
}

////////////////////////////
// Stencil< Dim , Width > //
////////////////////////////
template< typename T , unsigned int Width >
Stencil< T , Width >::Stencil( void ) : _values( NewPointer< T >( Width ) ){}

template< typename T , unsigned int Width >
Stencil< T , Width >::Stencil( const Stencil &stencil ) : Stencil() { for( unsigned int i=0 ; i<Width ; i++ ) _values[i] = stencil._values[i]; }

template< typename T , unsigned int Width >
Stencil< T , Width >::Stencil( Stencil &&stencil ) : Stencil() { std::swap( _values , stencil._values ); }

template< typename T , unsigned int Width >
Stencil< T , Width > &Stencil< T , Width >::operator = ( const Stencil &stencil )
{
	for( unsigned int i=0 ; i<Width ; i++ ) _values[i] = stencil._values[i];
	return *this;
}

template< typename T , unsigned int Width >
Stencil< T , Width > &Stencil< T , Width >::operator = ( Stencil &&stencil )
{
	for( unsigned int i=0 ; i<Width ; i++ ) _values[i] = std::move( stencil._values[i] );
	return *this;
}

template< typename T , unsigned int Width >
Stencil< T , Width >::~Stencil( void ){ if( _values ) DeletePointer( _values ); }

template< typename T , unsigned int Width >
void Stencil< T , Width >::Add( const Stencil &s ){ for( unsigned int i=0 ; i<Width ; i++ ) _values[i] += s._values[i]; }

template< typename T , unsigned int Width >
void Stencil< T , Width >::Scale( double s ){ for( unsigned int i=0 ; i<Width ; i++ ) _values[i] *= s; }

////////////////////////////////
// Stencil< Dim , Widths... > //
////////////////////////////////
template< typename T , unsigned int Width , unsigned int ... Widths >
Stencil< T , Width , Widths... >::Stencil( void ) : _values( NewPointer< SubStencil >( Width ) ){}

template< typename T , unsigned int Width , unsigned int ... Widths >
Stencil< T , Width , Widths... >::Stencil( const Stencil &stencil ) : Stencil() { for( unsigned int i=0 ; i<Width ; i++ ) _values[i] = stencil._values[i]; }

template< typename T , unsigned int Width , unsigned int ... Widths >
Stencil< T , Width , Widths... >::Stencil( Stencil &&stencil ) : Stencil() { std::swap( _values , stencil._values ); }

template< typename T , unsigned int Width , unsigned int ... Widths >
Stencil< T , Width , Widths... > &Stencil< T , Width , Widths... >::operator = ( const Stencil &stencil )
{
	for( unsigned int i=0 ; i<Width ; i++ ) _values[i] = stencil._values[i];
	return *this;
}

template< typename T , unsigned int Width , unsigned int ... Widths >
Stencil< T , Width , Widths... > &Stencil< T , Width , Widths... >::operator = ( Stencil &&stencil )
{
	for( unsigned int i=0 ; i<Width ; i++ ) _values[i] = std::move( stencil._values[i] );
	return *this;
}

template< typename T , unsigned int Width , unsigned int ... Widths >
Stencil< T , Width , Widths... >::~Stencil( void ){ if( _values ) DeletePointer( _values ); }

template< typename T , unsigned int Width , unsigned int ... Widths >
void Stencil< T , Width , Widths... >::Add( const Stencil &s ){ for( unsigned int i=0 ; i<Width ; i++ ) _values[i] += s._values[i]; }

template< typename T , unsigned int Width , unsigned int ... Widths >
void Stencil< T , Width , Widths... >::Scale( double s ){ for( unsigned int i=0 ; i<Width ; i++ ) _values[i] *= s; }
