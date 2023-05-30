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

template< unsigned int Dim >
SkewSymmetricMatrix< double , Dim > Wedge( Point< double , Dim > v1 , Point< double , Dim > v2 )
{
	SquareMatrix< double , Dim > skew;
	for( unsigned int i=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) skew(i,j) = v1[i] * v2[j] - v1[j] * v2[i];
	return Hat::SkewSymmetricMatrix< double , Dim >( skew );
}

///////////
// Index //
///////////
template< unsigned int Dim >
template< typename ... Indices >
Index< Dim > Index< Dim >::Min( Index i , Indices ... is )
{
	if constexpr( sizeof...(is)==0 ) return i;
	else
	{
		Index mn = Min( is ... );
		for( unsigned int d=0 ; d<Dim ; d++ ) mn[d] = std::min< int >( i[d] , mn[d] );
		return mn;
	}
}

template< unsigned int Dim >
template< typename ... Indices >
Index< Dim > Index< Dim >::Max( Index i , Indices ... is )
{
	if constexpr( sizeof...(is)==0 ) return i;
	else
	{
		Index mx = Max( is ... );
		for( unsigned int d=0 ; d<Dim ; d++ ) mx[d] = std::max< int >( i[d] , mx[d] );
		return mx;
	}
}

template< unsigned int Dim >
bool Index< Dim >::operator == ( Index i ) const
{
	for( unsigned int d=0 ; d<Dim ; d++ ) if( (*this)[d]!=i[d] ) return false;
	return true;
}

template< unsigned int Dim >
bool Index< Dim >::operator != ( Index i ) const { return !( (*this)==i ); }

template< unsigned int Dim >
bool Index< Dim >::operator < ( Index i ) const
{
	for( unsigned int d=0 ; d<Dim ; d++ )
		if( (*this)[d]<i[d] ) return true;
		else if( (*this)[d]>i[d] ) return false;
	return false;
}

template< unsigned int Dim >
bool Index< Dim >::operator <= ( Index i ) const
{
	for( unsigned int d=0 ; d<Dim ; d++ )
		if( (*this)[d]<i[d] ) return true;
		else if( (*this)[d]>i[d] ) return false;
	return true;
}

template< unsigned int Dim >
bool Index< Dim >::operator > ( Index i ) const { return !(*this<=i); }

template< unsigned int Dim >
bool Index< Dim >::operator >= ( Index i ) const { return !(*this<i); }

///////////
// Range //
///////////
template< unsigned int Dim >
template< typename ... Ranges > 
Range< Dim > Range< Dim >::Intersect( Ranges ... rs )
{
	Range r;
	r.first = Index< Dim >::Max( rs.first ... );
	r.second = Index< Dim >::Min( rs.second ... );
	return r;
}

template< unsigned int Dim >
bool Range< Dim >::empty( void ) const
{
	for( unsigned int d=0 ; d<Dim ; d++ ) if( first[d]>=second[d] ) return true;
	return false;
}

template< unsigned int Dim >
size_t Range< Dim >::size( void ) const
{
	size_t s = 1;
	for( unsigned int d=0 ; d<Dim ; d++ )
		if( first[d]>=second[d] ) return 0;
		else s *= second[d]-first[d];
	return s;
}

template< unsigned int Dim >
bool Range< Dim >::contains( Index< Dim > i ) const
{
	for( unsigned int d=0 ; d<Dim ; d++ ) if( i[d]<first[d] || i[d]>=second[d] ) return false;
	return true;
}

template< unsigned int Dim >
template< typename IndexFunctor > 
void Range< Dim >::process( IndexFunctor f ) const
{
	if constexpr( Dim==1 ) for( int i=first[0] ; i<second[0] ; i++ ) f( Index< Dim >(i) );
	else
	{
		Range< Dim-1 > _r;
		for( unsigned int d=0 ; d<Dim-1 ; d++ ) _r.first[d] = first[d+1] , _r.second[d] = second[d+1];

		Index< Dim > idx;
		auto _f = [&]( Index< Dim-1 > _idx )
		{
			for( unsigned int d=0 ; d<Dim-1 ; d++ ) idx[d+1] = _idx[d];
			f( idx );
		};
		for( int i=first[0] ; i<second[0] ; i++ )
		{
			idx[0] = i;
			_r.process( _f );
		}
	}
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

///////////////////
// ElementVector //
///////////////////
template< unsigned int Dim >
ElementVector< Dim >::ElementVector( const ElementFunction< Dim > &f ){ for( unsigned int d=0 ; d<Dim ; d++ ) _c[d] = f.d(d); }

///////////////////
// ElementVector //
///////////////////
template< unsigned int _Dim >
ElementWedge< _Dim >::ElementWedge( const ElementVector< _Dim > &v1 , const ElementVector< _Dim > &v2 )
{
	for( unsigned int d=0 , d1=0 ; d1<_Dim ; d1++ ) for( unsigned int d2=0 ; d2<d1 ; d2++ , d++ )
	{
		_c[0][d] =  v1[d1] * v2[d2];
		_c[1][d] = -v1[d2] * v2[d1];
	}
}

///////////
// Basis //
///////////
template< unsigned int Dim >
Range< Dim > Basis< Dim >::FunctionSupport( Index< Dim > i )
{
	Range< Dim > r;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		std::pair< unsigned int , unsigned int > _r = _Basis::FunctionSupport( i[d] );
		r.first[d] = _r.first , r.second[d] = _r.second;
	}
	return r;
}

template< unsigned int Dim >
Range< Dim > Basis< Dim >::ElementSupport( Index< Dim > i )
{
	Range< Dim > r;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		std::pair< unsigned int , unsigned int > _r = _Basis::ElementSupport( i[d] );
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
double Basis< Dim >::ElementGradientWedgeMass( unsigned int res , Index< Dim > e , std::pair< Index< Dim > , Index< Dim > > i , std::pair< Index< Dim > , Index< Dim > > j )
{
	double integral = 0;
	ElementFunction< Dim > f1 = ElementFunction< Dim >( e , i.first ) , f2 = ElementFunction< Dim >( e , i.second );
	ElementFunction< Dim > g1 = ElementFunction< Dim >( e , j.first ) , g2 = ElementFunction< Dim >( e , j.second );
	ElementWedge< Dim > w1 = ElementVector< Dim >( f1 ) ^ ElementVector< Dim >( f2 );
	ElementWedge< Dim > w2 = ElementVector< Dim >( g1 ) ^ ElementVector< Dim >( g2 );
	for( unsigned int d=0 ; d<ElementWedge< Dim >::Dim ; d++ ) for( unsigned int ii=0 ; ii<2 ; ii++ ) for( unsigned int jj=0 ; jj<2 ; jj++ ) integral += ( w1(ii,d) * w2(jj,d) ).integral( res );
	return integral;
}

template< unsigned int Dim >
template< typename F /* = std::function< double ( Point< double , Dim > ) */ >
double Basis< Dim >::Integral( unsigned int res , Index< Dim > e , F f , unsigned samplingRes )
{
	double sum = 0;
	Range< Dim > sRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) sRange.second[d] = samplingRes;

	auto g = [&]( Index< Dim > i )
	{
		Point< double , Dim > p;
		for( unsigned int d=0 ; d<3; d++ ) p[d] = e[d] + ( i[d] + 0.5 ) / samplingRes;
		p /= res;
		sum += f(p);
	};
	sRange.process( g );
	for( unsigned int d=0 ; d<Dim ; d++ ) sum /= samplingRes * res;
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
int Basis< Dim >::_RelativeIndex( Index< Dim > e , Index< Dim > i )
{
	unsigned int idx = 0;
	for( unsigned int d=0 ; d<Dim ; d++ ) if( e[d]!=i[d] ) idx |= 1<<d;
	return idx;
}

//////////////////////////////
// ScalarFunctions::Stencil //
//////////////////////////////
template< unsigned int Dim >
ScalarFunctions< Dim >::Stencil::Stencil( unsigned int res ) : _res(res)
{
	_values = new double *[1<<Dim];
	for( unsigned int i=0 ; i<(1<<Dim) ; i++ ) _values[i] = new double[1<<Dim];
	for( unsigned int d=0 ; d<Dim ; d++ ) _eRange.first[d] = 0 , _eRange.second[d] = _res;
}

template< unsigned int Dim >
ScalarFunctions< Dim >::Stencil::Stencil( const Stencil &stencil ) : _res(stencil._res) , _eRange(stencil._eRange)
{
	_values = new double *[1<<Dim];
	for( unsigned int i=0 ; i<(1<<Dim) ; i++ )
	{
		_values[i] = new double[1<<Dim];
		memcpy( _values[i] , stencil._values[i] , sizeof(double)*(1<<Dim) );
	}
}

template< unsigned int Dim >
ScalarFunctions< Dim >::Stencil::Stencil( Stencil &&stencil ) : _values(NULL)
{
	std::swap( _eRange , stencil._eRange );
	std::swap( _res , stencil._res );
	std::swap( _values , stencil._values );
}

template< unsigned int Dim >
typename ScalarFunctions< Dim >::Stencil &ScalarFunctions< Dim >::Stencil::operator =( const Stencil &stencil )
{
	_res = stencil._res;
	_eRange = stencil._eRange;
	for( unsigned int i=0 ; i<(1<<Dim) ; i++ ) memcpy( _values[i] , stencil._values[i] , sizeof(double)*(1<<Dim) );
	return *this;
}

template< unsigned int Dim >
ScalarFunctions< Dim >::Stencil::~Stencil( void )
{
	if( _values )
	{
		for( unsigned int i=0 ; i<(1<<Dim) ; i++ ) delete[] _values[i];
		delete[] _values;
	}
}

template< unsigned int Dim >
double ScalarFunctions< Dim >::Stencil::operator()( Index< Dim > e , Index< Dim > f1 , Index< Dim > f2 ) const
{
	return _values[ Basis< Dim >::_RelativeIndex( e , f1 ) ][ Basis< Dim >::_RelativeIndex( e , f2 ) ];
}

template< unsigned int Dim >
double ScalarFunctions< Dim >::Stencil::operator()( Index< Dim > f1 , Index< Dim > f2 ) const
{
	double value = 0;
	auto f = [&]( Index< Dim > e ){ value += this->operator()(e,f1,f2); };
	Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( f1 ) , Basis< Dim >::FunctionSupport( f2 ) , _eRange ).process( f );
	return value;
}

template< unsigned int Dim >
unsigned int ScalarFunctions< Dim >::Stencil::resolution( void ) const { return _res; }

template< unsigned int Dim >
Range< Dim > ScalarFunctions< Dim >::Stencil::eRange( void ) const { return _eRange; }

//////////////////////////////////
// ScalarFunctions::FullStencil //
//////////////////////////////////
template< unsigned int Dim >
constexpr unsigned int ScalarFunctions< Dim >::FullStencil::StencilNum( void )
{
	if constexpr( Dim==1 ) return 3;
	else return ScalarFunctions< Dim-1 >::FullStencil::StencilNum() * 3;
}


template< unsigned int Dim >
ScalarFunctions< Dim >::FullStencil::FullStencil( const Stencil &stencil ) : _res( stencil.resolution() ) , _rows( StencilNum() )
{
	Range< Dim > range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.first[d] = 0 , range.second[d] = 3;

	auto f = [&]( Index< Dim > f1 )
	{
		for( unsigned int d=0 ; d<Dim ; d++ )
			if     ( f1[d]==1 ) f1[d] = _res/2;
			else if( f1[d]==2 ) f1[d] = _res;

		Row _row;
		std::map< Index< Dim > , double > row;

		auto f = [&]( Index< Dim > e )
		{
			auto f = [&]( Index< Dim > f2 ){ row[ f2-f1 ] += stencil( e , f1 , f2 ); };
			Basis< Dim >::ElementSupport( e ).process(f);
		};
		Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( f1 ) , stencil.eRange() ).process( f );

		_row.reserve( row.size() );
		for( auto iter=row.begin() ; iter!=row.end() ; iter++ ) _row.push_back( Entry( iter->first , iter->second ) );

		_rows[ StencilIndex(f1,_res) ] = _row;
	};
	range.process( f );
}

template< unsigned int Dim >
unsigned int ScalarFunctions< Dim >::FullStencil::StencilIndex( Index< Dim > f , unsigned int res )
{
	unsigned int idx=0;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		if( f[d]<0 || f[d]>(int)res ) ERROR_OUT( "Bad function index: " , f , " / " , res );
		if     ( f[d]==0   ) idx = idx*3+0;
		else if( f[d]==res ) idx = idx*3+2;
		else                 idx = idx*3+1;
	}
	return idx;
}

template< unsigned int Dim >
const typename ScalarFunctions< Dim >::FullStencil::Row &ScalarFunctions< Dim >::FullStencil::row( Index< Dim > f ) const { return _rows[ StencilIndex( f , _res ) ]; }

/////////////////////
// ScalarFunctions //
/////////////////////
template< unsigned int Dim >
size_t ScalarFunctions< Dim >::functionNum( void ) const
{
	size_t sz = 1;
	for( unsigned int d=0 ; d<Dim ; d++ ) sz *= _r+1;
	return sz;
}

template< unsigned int Dim >
size_t ScalarFunctions< Dim >::elementNum( void ) const
{
	size_t sz = 1;
	for( unsigned int d=0 ; d<Dim ; d++ ) sz *= _r;
	return sz;
}

template< unsigned int Dim >
size_t ScalarFunctions< Dim >::elementIndex( Index< Dim > i ) const
{
	size_t idx=0;
	for( unsigned int d=0 ; d<Dim ; d++ ) idx = idx*_r + i[Dim-d-1];
	return idx;
}

template< unsigned int Dim >
Index< Dim > ScalarFunctions< Dim >::elementIndex( size_t i ) const
{
	Index< Dim > idx;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		idx[d] = i % _r;
		i /= _r;
	}
	return idx;
}

template< unsigned int Dim >
size_t ScalarFunctions< Dim >::functionIndex( Index< Dim > i ) const
{
	size_t idx=0;
	for( unsigned int d=0 ; d<Dim ; d++ ) idx = idx*(_r+1) + i[Dim-d-1];
	return idx;
}

template< unsigned int Dim >
Index< Dim > ScalarFunctions< Dim >::functionIndex( size_t i ) const
{
	Index< Dim > idx;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		idx[d] = i % (_r+1);
		i /= _r+1;
	}
	return idx;
}

template< unsigned int Dim >
typename ScalarFunctions< Dim >::Stencil ScalarFunctions< Dim >::MassStencil( unsigned int res )
{
	ScalarFunctions< Dim >::Stencil stencil(res);

	Index< Dim > e;
	Range< Dim > r = Basis< Dim >::ElementSupport( e );

	auto f = [&]( Index< Dim > i1 )
	{
		auto f = [&]( Index< Dim > i2 ){ stencil._values[ Basis< Dim >::_RelativeIndex(e,i1) ][ Basis< Dim >::_RelativeIndex(e,i2) ] = ( ElementFunction( e , i1 ) * ElementFunction( e , i2 ) ).integral(res); };
		r.process( f );
	};
	r.process( f );

	return stencil;
}

template< unsigned int Dim >
typename ScalarFunctions< Dim >::Stencil ScalarFunctions< Dim >::StiffnessStencil( unsigned int res )
{
	ScalarFunctions< Dim >::Stencil stencil(res);

	Index< Dim > e;
	Range< Dim > r = Basis< Dim >::ElementSupport( e );

	auto f = [&]( Index< Dim > i1 )
	{
		auto f = [&]( Index< Dim > i2 )
		{
			ElementVector< Dim > v1( ElementFunction< Dim >(e,i1) ) , v2( ElementFunction< Dim >(e,i2) );
			double integral = 0;
			for( unsigned int d=0 ; d<Dim ; d++ ) integral +=( v1[d] * v2[d] ).integral(res);
			stencil._values[ Basis< Dim >::_RelativeIndex(e,i1) ][ Basis< Dim >::_RelativeIndex(e,i2) ] = integral;
		};
		r.process( f );
	};
	r.process( f );

	return stencil;
}

template< unsigned int Dim >
double ScalarFunctions< Dim >::dot( const Stencil &stencil , const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const
{
	std::vector< double > integrals( omp_get_max_threads() , 0 );
#pragma omp parallel for
	for( int i=0 ; i<elementNum() ; i++ )
	{
		unsigned int t = omp_get_thread_num();
		Index< Dim > e = elementIndex( i );

		auto f = [&]( Index< Dim > f1 )
		{
			size_t i1 = functionIndex( f1 );
			auto f = [&]( Index< Dim > f2 )
			{
				size_t i2 = functionIndex( f1 );
				integrals[t] += stencil( e , f1 , f2 ) * x[i1] * y[i2];
			};
			Range< Dim >::Intersect( Basis< Dim >::ElementSupport( e ) ).process( f );
		};
		Range< Dim >::Intersect( Basis< Dim >::ElementSupport( e ) ).process( f );
	}

	double integral = 0;
	for( unsigned int i=0 ; i<integrals.size() ; i++ ) integral += integrals[i];
	return integral;
}

template< unsigned int Dim >
double ScalarFunctions< Dim >::dot( const FullStencil &stencil , const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const
{
	std::vector< double > integrals( omp_get_max_threads() , 0 );
#pragma omp parallel for
	for( int i=0 ; i<functionNum() ; i++ )
	{
		unsigned int t = omp_get_thread_num();
		Index< Dim > f1 = functionIndex( i );
		const typename Hat::ScalarFunctions< Dim >::FullStencil::Row &row = stencil.row( f1 );
		for( unsigned int j=0 ; j<row.size() ; j++ )
		{
			Index< Dim > f2 = f1 + row[j].first;
			integrals[t] += row[j].second * x[i] * y[ functionIndex(f2) ];
		}
	}

	double integral = 0;
	for( unsigned int i=0 ; i<integrals.size() ; i++ ) integral += integrals[i];
	return integral;
}

template< unsigned int Dim >
Eigen::SparseMatrix< char > ScalarFunctions< Dim >::incidence( void ) const
{
	Stencil stencil = MassStencil(_r);
	OMPTriplets< char > _triplets;

	FullStencil fullStencil( stencil );
#pragma omp parallel for
	for( int i=0 ; i<functionNum() ; i++ )
	{
		int t = omp_get_thread_num();
		size_t i1 = (size_t)i;
		Index< Dim > f1 = functionIndex( i1 );
		const typename Hat::ScalarFunctions< Dim >::FullStencil::Row &row = fullStencil.row( f1 );
		for( unsigned int j=0 ; j<row.size() ; j++ )
		{
			Index< Dim > f2 = f1 + row[j].first;
			size_t i2 = functionIndex( f2 );
			if( i1<=i2 ) _triplets[t].push_back( Eigen::Triplet< char >( (int)i1 , (int)i2 , 0 ) );
		}
	}

	Eigen::SparseMatrix< char > M( functionNum() , functionNum() );
	M.setFromTriplets( _triplets.begin() , _triplets.end() );
	_triplets.clear();
	Eigen::SparseMatrix< char > Mt = M.transpose();
	return M + Mt;
}

template< unsigned int Dim >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::_systemMatrix( Stencil stencil ) const
{
	OMPTriplets< double > _triplets;

	FullStencil fullStencil( stencil );
#pragma omp parallel for
	for( int i=0 ; i<functionNum() ; i++ )
	{
		int t = omp_get_thread_num();
		size_t i1 = (size_t)i;
		Index< Dim > f1 = functionIndex( i1 );
		const typename Hat::ScalarFunctions< Dim >::FullStencil::Row &row = fullStencil.row( f1 );
		for( unsigned int j=0 ; j<row.size() ; j++ )
		{
			Index< Dim > f2 = f1 + row[j].first;
			size_t i2 = functionIndex( f2 );
			if( i1<=i2 )
			{
				double integral = row[j].second;
				if( i1==i2 ) integral /= 2;
				_triplets[t].push_back( Eigen::Triplet< double >( (int)i1 , (int)i2 , integral ) );
			}
		}
	}

	Eigen::SparseMatrix< double > M( functionNum() , functionNum() );
	M.setFromTriplets( _triplets.begin() , _triplets.end() );
	_triplets.clear();
	Eigen::SparseMatrix< double > Mt = M.transpose();
	return M + Mt;
}

template< unsigned int Dim >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::mass( void ) const { return _systemMatrix( MassStencil(_r) ); }

template< unsigned int Dim >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::stiffness( void ) const { return _systemMatrix( StiffnessStencil(_r) ); }

template< unsigned int Dim >
template< typename SampleFunctor >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::evaluation( SampleFunctor F , size_t sampleNum ) const
{
	std::vector< Eigen::Triplet< double > > triplets;
	triplets.reserve( sampleNum*(1<<Dim) );
	for( unsigned int s=0 ; s<sampleNum ; s++ )
	{
		Point< double , Dim > p = F(s) * _r;
		Index< Dim > e;
		for( unsigned int d=0 ; d<Dim ; d++ ) e[d] = (int)floor( p[d] );
		for( unsigned int d=0 ; d<Dim ; d++ ) if( e[d]==_r ) e[d]--;
		for( unsigned int d=0 ; d<Dim ; d++ ) p[d] -= e[d];

		double values[2][Dim];
		for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) values[i][d] = i==0 ? (1.-p[d]) : p[d];
		auto f = [&]( Index< Dim > f )
		{
			double value = 1.;
			for( unsigned int d=0 ; d<Dim ; d++ ) value *= f[d]==e[d] ? values[0][d] : values[1][d];
			triplets.push_back( Eigen::Triplet< double >( s , (unsigned int)functionIndex( f ) , value ) );
		};
		Basis< Dim >::ElementSupport( e ).process( f );
	}

	Eigen::SparseMatrix< double > E( sampleNum , functionNum() );
	E.setFromTriplets( triplets.begin() , triplets.end() );
	return E;
}
template< unsigned int Dim >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::evaluation( const std::vector< Point< double , Dim > > &samples ) const
{
	return evaluation( [&]( unsigned int idx ){ return samples[idx]; } , samples.size() );
}

template< unsigned int Dim >
template< typename SampleFunctor >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::dEvaluation( SampleFunctor F , size_t sampleNum ) const
{
	std::vector< Eigen::Triplet< double > > triplets;
	triplets.reserve( Dim*sampleNum*(1<<Dim) );
	for( unsigned int s=0 ; s<sampleNum ; s++ )
	{
		Point< double , Dim > p = F(s) * _r;
		Index< Dim > e;
		for( unsigned int d=0 ; d<Dim ; d++ ) e[d] = (int)floor( p[d] );
		for( unsigned int d=0 ; d<Dim ; d++ ) if( e[d]==_r ) e[d]--;
		for( unsigned int d=0 ; d<Dim ; d++ ) p[d] -= e[d];

		double values[2][Dim];
		for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) values[i][d] = i==0 ? (1.-p[d]) : p[d];
		auto f = [&]( Index< Dim > f )
		{
			for( unsigned int dd=0 ; dd<Dim ; dd++ )
			{
				double dValue = 1<<_r;
				for( unsigned int d=0 ; d<Dim ; d++ )
					if( d==dd ) dValue *= f[d]==e[d] ? -1. : 1.;
					else        dValue *= f[d]==e[d] ? values[0][d] : values[1][d];
				triplets.push_back( Eigen::Triplet< double >( s*Dim+dd , (unsigned int)functionIndex( f ) , dValue ) );
			}
		};
		Basis< Dim >::ElementSupport( e ).process( f );
	}

	Eigen::SparseMatrix< double > E( sampleNum*Dim , functionNum() );
	E.setFromTriplets( triplets.begin() , triplets.end() );
	return E;
}
template< unsigned int Dim >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::dEvaluation( const std::vector< Point< double , Dim > > &samples ) const
{
	return dEvaluation( [&]( unsigned int idx ){ return samples[idx]; } , samples.size() );
}

template< unsigned int Dim >
double ScalarFunctions< Dim >::value( const Eigen::VectorXd &x , Point< double , Dim > p ) const
{
	p *= _r;
	Index< Dim > e;
	for( unsigned int d=0 ; d<Dim ; d++ ) e[d] = (int)floor( p[d] );
	for( unsigned int d=0 ; d<Dim ; d++ ) if( e[d]==_r ) e[d]--;
	for( unsigned int d=0 ; d<Dim ; d++ ) p[d] -= e[d];

	double values[2][Dim];
	for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) values[i][d] = i==0 ? (1.-p[d]) : p[d];

	double value = 0;
	auto f = [&]( Index< Dim > f )
	{
		double v = x[ functionIndex(f) ];
		for( unsigned int d=0 ; d<Dim ; d++ ) v *= f[d]==e[d] ? values[0][d] : values[1][d];
		value += v;
	};
	Basis< Dim >::ElementSupport( e ).process(f);
	return value;
}

template< unsigned int Dim >
Point< double , Dim > ScalarFunctions< Dim >::gradient( const Eigen::VectorXd &x , Point< double , Dim > p ) const
{
	p *= _r;
	Index< Dim > e;
	for( unsigned int d=0 ; d<Dim ; d++ ) e[d] = (int)floor( p[d] );
	for( unsigned int d=0 ; d<Dim ; d++ ) if( e[d]==_r ) e[d]--;
	for( unsigned int d=0 ; d<Dim ; d++ ) p[d] -= e[d];

	double values[2][Dim] , dValues[2][Dim];
	for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ )
	{
		values[i][d] = i==0 ? (1.-p[d]) : p[d];
		dValues[i][d] = i==0 ? -1./_r : 1./_r;
	}

	Point< double , Dim > gradient;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		auto f = [&]( Index< Dim > f )
		{
			double v = x[ functionIndex(f) ];
			for( unsigned int dd=0 ; dd<Dim ; dd++ )
				if( d!=dd ) v *= f[dd]==e[dd] ?  values[0][dd] :  values[1][dd];
				else        v *= f[dd]==e[dd] ? dValues[0][dd] : dValues[1][dd];
			gradient[d] += v;
		};
		Basis< Dim >::ElementSupport( e ).process(f);
	}
	return gradient;
}

template< unsigned int Dim >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::prolongation( void ) const
{
	if( _r&1 ) ERROR_OUT( "Expected even resolution: " , _r );
	ScalarFunctions coarseScalars( _r>>1 );

	std::vector< Eigen::Triplet< double > > triplets;
	{
		unsigned int x=1;
		for( unsigned int d=0 ; d<Dim ; d++ ) x *= 3;
		triplets.reserve( coarseScalars.functionNum() * x );
	}
	for( size_t c=0 ; c<coarseScalars.functionNum() ; c++ )
	{
		Index< Dim > coarse = coarseScalars.functionIndex(c);
		Index< Dim > fine;
		for( unsigned int d=0 ; d<Dim ; d++ ) fine[d] = coarse[d]*2;
		Range< Dim > r;
		for( unsigned int d=0 ; d<Dim ; d++ ) r.first[d] = fine[d]-1 , r.second[d] = fine[d]+2;
		auto f = [&]( Index< Dim > _fine )
		{
			double v = 1;
			for( unsigned int d=0 ; d<Dim ; d++ ) if( _fine[d]!=fine[d] ) v *= 0.5;
			triplets.push_back( Eigen::Triplet< double >( (unsigned int)functionIndex(_fine) , (unsigned int)coarseScalars.functionIndex(coarse) , v ) );
		};
		Range< Dim >::Intersect( r , _fRange ).process( f );
	}

	Eigen::SparseMatrix< double > P( functionNum() , coarseScalars.functionNum() );
	P.setFromTriplets( triplets.begin() , triplets.end() );
	return P;
};

/////////////////////////////
// WedgeFunctions::Stencil //
/////////////////////////////
template< unsigned int Dim >
WedgeFunctions< Dim >::Stencil::Stencil( unsigned int res ) : _res(res)
{
	_values = new double ***[1<<Dim];
	for( unsigned int i=0 ; i<(1<<Dim) ; i++ )
	{
		_values[i] = new double**[1<<Dim];
		for( unsigned int j=0 ; j<(1<<Dim) ; j++ )
		{
			_values[i][j] = new double*[1<<Dim];
			for( unsigned int k=0 ; k<(1<<Dim) ; k++ ) _values[i][j][k] = new double[1<<Dim];
		}
	}
	for( unsigned int d=0 ; d<Dim ; d++ ) _eRange.first[d] = 0 , _eRange.second[d] = _res;
}

template< unsigned int Dim >
WedgeFunctions< Dim >::Stencil::Stencil( const Stencil &stencil ) : _res(stencil._res) , _eRange(stencil._eRange)
{
	_values = new double ***[1<<Dim];
	for( unsigned int i=0 ; i<(1<<Dim) ; i++ )
	{
		_values[i] = new double**[1<<Dim];
		for( unsigned int j=0 ; j<(1<<Dim) ; j++ )
		{
			_values[i][j] = new double*[1<<Dim];
			for( unsigned int k=0 ; k<(1<<Dim) ; k++ )
			{
				_values[i][j][k] = new double[1<<Dim];
				memcpy( _values[i][j][k] , stencil._values[i][j][k] , sizeof(double)*(1<<Dim) );
			}
		}
	}
}

template< unsigned int Dim >
WedgeFunctions< Dim >::Stencil::Stencil( Stencil &&stencil ) : _values(NULL)
{
	std::swap( _res , stencil._res );
	std::swap( _eRange , stencil._eRange );
	std::swap( _values , stencil._values );
}

template< unsigned int Dim >
typename WedgeFunctions< Dim >::Stencil & WedgeFunctions< Dim >::Stencil::operator = ( const Stencil &stencil )
{
	_res = stencil._res;
	_eRange = stencil._eRange;
	for( unsigned int i=0 ; i<(1<<Dim) ; i++ ) for( unsigned int j=0 ; j<(1<<Dim) ; j++ ) for( unsigned int k=0 ; k<(1<<Dim) ; k++ ) memcpy( _values[i][j][k] , stencil._values[i][j][k] , sizeof(double)*(1<<Dim) );
}

template< unsigned int Dim >
WedgeFunctions< Dim >::Stencil::~Stencil( void )
{
	if( _values )
	{
		for( unsigned int i=0 ; i<(1<<Dim) ; i++ )
		{
			for( unsigned int j=0 ; j<(1<<Dim) ; j++ )
			{
				for( unsigned int k=0 ; k<(1<<Dim) ; k++ ) delete[] _values[i][j][k];
				delete[] _values[i][j];
			}
			delete[] _values[i];
		}
		delete[] _values;
	}
}

template< unsigned int Dim >
double WedgeFunctions< Dim >::Stencil::operator()( Index< Dim > e , Index< Dim > f1 , Index< Dim > f2 , Index< Dim > g1 , Index< Dim > g2 ) const
{
	return _values[ Basis< Dim >::_RelativeIndex(e,f1) ][ Basis< Dim >::_RelativeIndex(e,f2) ][ Basis< Dim >::_RelativeIndex(e,g1) ][ Basis< Dim >::_RelativeIndex(e,g2) ];
}

template< unsigned int Dim >
double WedgeFunctions< Dim >::Stencil::operator()( Index< Dim > f1 , Index< Dim > f2 , Index< Dim > g1 , Index< Dim > g2 ) const
{
	double value = 0;
	auto f = [&]( Index< Dim > e ){ value += this->operator()(e,f1,f2,g1,g2); };
	Range< Dim >::Intersect( _eRange , Basis< Dim >::FunctionSupport(f1) , Basis< Dim >::FunctionSupport(f2) , Basis< Dim >::FunctionSupport(g1) , Basis< Dim >::FunctionSupport(g2) ).process( f );
	return value;
}

template< unsigned int Dim >
unsigned int WedgeFunctions< Dim >::Stencil::resolution( void ) const { return _res; }

template< unsigned int Dim >
Range< Dim > WedgeFunctions< Dim >::Stencil::eRange( void ) const { return _eRange; }

/////////////////////////////////
// WedgeFunctions::FullStencil //
/////////////////////////////////
template< unsigned int Dim >
WedgeFunctions< Dim >::FullStencil::FullStencil( const Stencil &stencil ) : _res(stencil.resolution() ) , _rows( ScalarFunctions< Dim >::FullStencil::StencilNum() )
{
	Range< Dim > range , eRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.first[d] = 0 , range.second[d] = 3 , eRange.first[d] = 0 , eRange.second[d] = _res;

	struct WedgeIndex
	{
		Index< Dim > g1 , g2;

		WedgeIndex( Index< Dim > g1 , Index< Dim > g2 ) : g1(g1) , g2(g2) {}
		bool operator < ( WedgeIndex i ) const
		{
			if( g1<i.g1 ) return true;
			else if( g1>i.g1 ) return false;
			else return g2<i.g2;
		}
	};

	auto f = [&]( Index< Dim > f1 )
	{
		// Transform the stencil index to a function index
		for( unsigned int d=0 ; d<Dim ; d++ )
			if     ( f1[d]==1 ) f1[d] = _res/2;
			else if( f1[d]==2 ) f1[d] = _res;

		// Get the indices of all functions whose support overlaps the support of f1
		std::vector< Index< Dim > > f2s;
		auto f = [&]( Index< Dim > f2 ){ if( f2!=f1 ) f2s.push_back( f2 ); };
		Basis< Dim >::ElementSupport( Range< Dim >::Intersect( eRange , Basis< Dim >::FunctionSupport( f1 ) ) ).process(f);

		std::vector< Row > _row( f2s.size() );
		for( unsigned int i=0 ; i<f2s.size() ; i++ )
		{
			Index< Dim > f2 = f2s[i];
			_row[i].f2 = f2 - f1;

			// The range of elements supported on both functions
			std::map< WedgeIndex , double > row;

			auto f = [&]( Index< Dim > e )
			{
				auto f = [&]( Index< Dim > g1 )
				{
					auto f = [&]( Index< Dim > g2 )
					{
						if( g1<g2 ) row[ WedgeIndex( g1 , g2 ) ] += stencil( e , f1 , f2 , g1 , g2 );
					};
					Basis< Dim >::ElementSupport( e ).process(f);
				};
				Basis< Dim >::ElementSupport( e ).process(f);
			};
			Range< Dim >::Intersect( eRange , Basis< Dim >::FunctionSupport(f1) , Basis< Dim >::FunctionSupport(f2) ).process( f );
			_row[i].entries.reserve( row.size() );
			for( auto iter=row.begin() ; iter!=row.end() ; iter++ ) _row[i].entries.push_back( Entry( iter->first.g1-f1 , iter->first.g2-f1 , iter->second ) );
		}
		_rows[ ScalarFunctions< Dim >::FullStencil::StencilIndex(f1,_res) ] = _row;
	};
	range.process( f );
}

template< unsigned int Dim >
const std::vector< typename WedgeFunctions< Dim >::FullStencil::Row > &WedgeFunctions< Dim >::FullStencil::rows( Index< Dim > f1 ) const { return _rows[ ScalarFunctions< Dim >::FullStencil::StencilIndex( f1 , _res ) ]; }

////////////////////
// WedgeFunctions //
////////////////////
template< unsigned int Dim >
WedgeFunctions< Dim >::WedgeFunctions( unsigned int resolution ) : _r(resolution)
{
	ScalarFunctions< Dim > scalars( _r );
	Range< Dim > fRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) fRange.first[d] = 0 , fRange.second[d] = _r+1;

	for( size_t i=0 ; i<scalars.functionNum() ; i++ )
	{
		Index< Dim > idx = scalars.functionIndex(i);
		Range< Dim > r = Range< Dim >::Intersect( Basis< Dim >::ElementSupport( Basis< Dim >::FunctionSupport( idx ) ) , fRange );

		auto f = [&]( Index< Dim > _idx )
		{
			size_t _i = scalars.functionIndex( _idx );
			if( idx<_idx ) _indexMap[ std::make_pair(i,_i) ] = 0;
		};
		r.process( f );
	}
	_indices.resize( _indexMap.size() );

	size_t count=0;
	for( auto iter=_indexMap.begin() ; iter!=_indexMap.end() ; iter++ )
	{
		iter->second = count;
		_indices[ count ] = iter->first;
		count++;
	}
}

template< unsigned int Dim >
size_t WedgeFunctions< Dim >::index( std::pair< size_t , size_t > idx , bool &flip ) const
{
	auto iter = _indexMap.find( idx );
	if( iter!=_indexMap.end() ) flip = false;
	else
	{
		iter = _indexMap.find( std::make_pair( idx.second , idx.first ) );
		if( iter!=_indexMap.end() ) flip = true;
		else
		{
			ScalarFunctions< Dim > scalars( _r );
			ERROR_OUT( "Could not find index in map: " , idx.first , " , "  , idx.second , " -> " , scalars.functionIndex( idx.first ) , " , "  , scalars.functionIndex( idx.second ) );
		}
	}
	return iter->second;
}

template< unsigned int Dim >
bool WedgeFunctions< Dim >::setIndex( std::pair< size_t , size_t > idx , size_t &i , bool &flip ) const
{
	auto iter = _indexMap.find( idx );
	if( iter!=_indexMap.end() )
	{
		i = iter->second;
		flip = false;
		return true;
	}

	iter = _indexMap.find( std::make_pair( idx.second , idx.first ) );
	if( iter!=_indexMap.end() )
	{
		i = iter->second;
		flip = true;
		return true;
	}

	return false;
}

template< unsigned int Dim >
std::pair< size_t , size_t > WedgeFunctions< Dim >::index( size_t idx ) const
{
	return _indices[idx];
}

template< unsigned int Dim >
Eigen::VectorXd WedgeFunctions< Dim >::wedge( const Eigen::VectorXd &f1 , const Eigen::VectorXd &f2 ) const
{
	Eigen::VectorXd w( _indices.size() );
	for( unsigned int i=0 ; i<_indices.size() ; i++ ) w[i] = f1[ _indices[i].first ] * f2[ _indices[i].second ] - f1[ _indices[i].second ] * f2[ _indices[i].first ];
	return w;
}


template< unsigned int Dim >
Eigen::SparseMatrix< double > WedgeFunctions< Dim >::wedge( const Eigen::VectorXd &f ) const
{
	ScalarFunctions< Dim > scalars( _r );
	if( f.size()!=scalars.functionNum() ) ERROR_OUT( "Resolutions don't match: " , f.size() , " != " , scalars.functionNum() );

	OMPTriplets< double > _triplets;
	for( unsigned int t=0 ; t<_triplets.size() ; t++ ) _triplets[t].reserve( 2*_indices.size() );
#pragma omp parallel for
	for( int i=0 ; i<_indices.size() ; i++ )
	{
		int t = omp_get_thread_num();
		_triplets[t].push_back( Eigen::Triplet< double >( i , (unsigned int)_indices[i].first , -f[ _indices[i].second ] ) );
		_triplets[t].push_back( Eigen::Triplet< double >( i , (unsigned int)_indices[i].second ,  f[ _indices[i].first ] ) );
	}

	Eigen::SparseMatrix< double > W( functionNum() , scalars.functionNum() );
	W.setFromTriplets( _triplets.begin() , _triplets.end() );
	_triplets.clear();
	return W;
}

template< unsigned int Dim >
typename WedgeFunctions< Dim >::Stencil WedgeFunctions< Dim >::MassStencil( unsigned int res )
{
	Stencil stencil(res);

	Index< Dim > e;
	Range< Dim > r = Basis< Dim >::ElementSupport( e );
	auto f = [&]( Index< Dim > i1 )
	{
		auto f = [&]( Index< Dim > i2 )
		{
			auto f = [&]( Index< Dim > j1 )
			{
				auto f = [&]( Index< Dim > j2 )
				{
					stencil._values[ Basis< Dim >::_RelativeIndex(e,i1) ][ Basis< Dim >::_RelativeIndex(e,i2) ][ Basis< Dim >::_RelativeIndex(e,j1) ][ Basis< Dim >::_RelativeIndex(e,j2) ] = Basis< Dim >::ElementGradientWedgeMass( res , e , std::make_pair( i1 , i2 ) , std::make_pair( j1 , j2 ) );
				};
				r.process( f );
			};
			r.process( f );
		};
		r.process( f );
	};
	r.process( f );

	return stencil;
}

template< unsigned int Dim >
double WedgeFunctions< Dim >::dot( const Stencil &stencil , const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const
{
	ScalarFunctions< Dim > scalars( (unsigned int)resolution() );
	std::vector< double > integrals( omp_get_max_threads() , 0 );

#pragma omp parallel for
	for( int i=0 ; i<scalars.elementNum() ; i++ )
	{
		unsigned int t = omp_get_thread_num();
		Index< Dim > e = scalars.elementIndex( i );

		auto f = [&]( Index< Dim > f1 )
		{
			size_t i1 = scalars.functionIndex( f1 );
			auto f = [&]( Index< Dim > f2 )
			{
				size_t i2 = scalars.functionIndex( f2 );
				if( i1!=i2 )
				{
					bool flipI;
					size_t idxI = index( std::make_pair( i1 , i2 ) , flipI );
					if( !flipI )
					{
						auto f = [&]( Index< Dim > g1 )
						{
							size_t j1 = scalars.functionIndex( g1 );
							auto f = [&]( Index< Dim > g2 )
							{
								size_t j2 = scalars.functionIndex( g2 );
								if( j1!=j2 )
								{
									bool flipJ;
									size_t idxJ = index( std::make_pair( j1 , j2 ) , flipJ );
									if( !flipJ ) integrals[t] += stencil( e , f1 , f2 , g1 , g2 ) * x[idxI] * y[idxJ];
								}
							};
							Basis< Dim >::ElementSupport( e ).process( f );
						};
						Basis< Dim >::ElementSupport( e ).process( f );
					};
				}
			};
			Basis< Dim >::ElementSupport( e ).process( f );
		};
		Basis< Dim >::ElementSupport( e ).process( f );
	}

	double integral = 0;
	for( unsigned int i=0 ; i<integrals.size() ; i++ ) integral += integrals[i];

	return integral;
}

template< unsigned int Dim >
double WedgeFunctions< Dim >::dot( const FullStencil &stencil , const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const
{
	ScalarFunctions< Dim > scalars( (unsigned int)resolution() );
	std::vector< double > integrals( omp_get_max_threads() , 0 );

#pragma omp parallel for
	for( int i=0 ; i<scalars.functionNum() ; i++ )
	{
		unsigned int t = (unsigned int)omp_get_thread_num();
		size_t i1 = (size_t)i;
		Index< Dim > f1 = scalars.functionIndex( i );
		const std::vector< typename Hat::WedgeFunctions< Dim >::FullStencil::Row > &rows = stencil.rows( f1 );
		for( unsigned int j=0 ; j<rows.size() ; j++ )
		{
			Index< Dim > f2 = f1 + rows[j].f2;
			if( f1<f2 )
			{
				size_t i2 = scalars.functionIndex(f2);
				bool flipI;
				size_t idxI = index( std::make_pair(i1,i2) , flipI );
				for( unsigned int k=0 ; k<rows[j].entries.size() ; k++ )
				{
					Index< Dim > g1 = f1 + rows[j].entries[k].g1 , g2 = f1 + rows[j].entries[k].g2;
					size_t j1 = scalars.functionIndex(g1) , j2 = scalars.functionIndex(g2);
					bool flipJ;
					size_t idxJ = index( std::make_pair(j1,j2) , flipJ );
					integrals[t] += rows[j].entries[k].value * x[idxI] * y[idxJ];
				}
			}
		}
	}

	double integral = 0;
	for( unsigned int i=0 ; i<integrals.size() ; i++ ) integral += integrals[i];

	return integral;
}

template< unsigned int Dim >
Eigen::SparseMatrix< double > WedgeFunctions< Dim >::mass( void ) const
{
	ScalarFunctions< Dim > scalars( _r );
	Range< Dim > eRange , fRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) fRange.first[d] = eRange.first[d] = 0 , fRange.second[d] = _r+1 , eRange.second[d] = _r;

	auto RelativeIndex = []( Index< Dim > e , Index< Dim > i ){ return Basis< Dim >::_RelativeIndex( e , i ); };

	Stencil stencil = MassStencil( _r );

	OMPTriplets< double > _triplets;

	// Iterate over all wedge-basis functions
#pragma omp parallel for
	for( int i=0 ; i<functionNum() ; i++ )
	{
		int t = omp_get_thread_num();

		// The pair of functions associated with index i
		Index< Dim > i1 = scalars.functionIndex( _indices[i].first ) , i2 = scalars.functionIndex( _indices[i].second );

		auto f = [&]( Index< Dim > j1 )
		{
			auto f = [&]( Index< Dim > j2 )
			{
				if( j1==j2 ) return;


				bool flip;
				size_t j = index( std::pair< size_t , size_t >( scalars.functionIndex( j1 ) , scalars.functionIndex( j2 ) ) , flip );
				if( !flip && i<=j )
				{
					double integral = stencil( i1 , i2 , j1 , j2 );
					if( i==j ) integral /= 2;
					if( fabs(integral)>1e-10 ) _triplets[t].push_back( Eigen::Triplet< double >( i , (unsigned int)j , integral * ( flip ? -1. : 1. ) ) );
				}
			};
			// Iterate over functions supported on the intersection of the support of i1, i2, and j1;
			Range< Dim >::Intersect( fRange , Basis< Dim >::ElementSupport( Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( i1 ) , Basis< Dim >::FunctionSupport( i2 ) , Basis< Dim >::FunctionSupport( j1 ) ) ) ).process( f );
		};
		// Iterate over functions supported on the intersection of the support of i1 and i2;
		Range< Dim >::Intersect( fRange , Basis< Dim >::ElementSupport( Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( i1 ) , Basis< Dim >::FunctionSupport( i2 ) ) ) ).process( f );
	}

	Eigen::SparseMatrix< double > M( functionNum() , functionNum() );

	M.setFromTriplets( _triplets.begin() , _triplets.end() );
	_triplets.clear();

	Eigen::SparseMatrix< double > Mt = M.transpose();
	return M + Mt;
}

template< unsigned int Dim >
Eigen::SparseMatrix< double > WedgeFunctions< Dim >::prolongation( void ) const
{
	if( _r&1 ) ERROR_OUT( "Expected even resolution: " , _r );
	WedgeFunctions coarseWedges( _r>>1 );
	ScalarFunctions< Dim > scalars( _r );
	Eigen::SparseMatrix< double > sP = scalars.prolongation();

	OMPTriplets< double > _triplets;
#pragma omp parallel for
	for( int c=0 ; c<(int)coarseWedges.functionNum() ; c++ )
	{
		int t = omp_get_thread_num();
		size_t c1 = coarseWedges._indices[c].first , c2 = coarseWedges._indices[c].second;
		for( Eigen::InnerIterator it1(sP,c1) ; it1 ; ++it1 ) for( Eigen::InnerIterator it2(sP,c2) ; it2 ; ++it2 )
		{
			size_t f1 = (size_t)it1.row() , f2 = (size_t)it2.row();
			bool flip;
			size_t f;
			if( setIndex( std::make_pair( f1 , f2 ) , f , flip ) ) _triplets[t].push_back( Eigen::Triplet< double >( (unsigned int)f , (unsigned int)c , it1.value()*it2.value() * ( flip ? -1. : 1. ) ) );
		}
	}
	Eigen::SparseMatrix< double > P( functionNum() , coarseWedges.functionNum() );
	P.setFromTriplets( _triplets.begin() , _triplets.end() );
	_triplets.clear();

	return P;
}

template< unsigned int Dim >
Eigen::SparseMatrix< double > WedgeFunctions< Dim >::dualMatrix( ConstPointer( SkewSymmetricMatrix< double , Dim > ) skew ) const
{
	ScalarFunctions< Dim > scalars( _r );
	Range< Dim > eRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) eRange.first[d] = 0 , eRange.second[d] = _r;

	OMPTriplets< double > _triplets;

#pragma omp parallel for
	for( int i=0 ; i<(int)functionNum() ; i++ )
	{
		unsigned int t = (unsigned int)omp_get_thread_num();
		size_t i1 = _indices[i].first , i2 = _indices[i].second;
		Index< Dim > f1 = scalars.functionIndex( i1 ) , f2 = scalars.functionIndex( i2 );
		double dot = 0;
		auto f = [&]( Index< Dim > e )
		{
			const SkewSymmetricMatrix< double , Dim > &w = skew[ scalars.elementIndex(e) ];
			double integral = 0;
			ElementFunction< Dim > ef1 = ElementFunction< Dim >( e , f1 ) , ef2 = ElementFunction< Dim >( e , f2 );
			ElementWedge< Dim > wedge = ElementVector< Dim >( ef1 ) ^ ElementVector< Dim >( ef2 );
			for( unsigned int d=0 ; d<ElementWedge< Dim >::Dim ; d++ ) for( unsigned int ii=0 ; ii<2 ; ii++ ) dot += w[d] * wedge(ii,d).integral( _r );
		};
		// Iterate over all cells supported by the i-th function
		Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( f1 ) , Basis< Dim >::FunctionSupport( f2 ) , eRange ).process( f );

		if( dot )
		{
			_triplets[t].push_back( Eigen::Triplet< double >( (int)i2 , (int)i1 ,  dot ) );
			_triplets[t].push_back( Eigen::Triplet< double >( (int)i1 , (int)i2 , -dot ) );
		}
	}

	Eigen::SparseMatrix< double > d( scalars.functionNum() , scalars.functionNum() );
	d.setFromTriplets( _triplets.begin() , _triplets.end() );
//	_triplets.clear();
	return d;
}

template< unsigned int Dim >
Eigen::VectorXd WedgeFunctions< Dim >::dualVector( ConstPointer( SkewSymmetricMatrix< double , Dim > ) skew ) const
{
	ScalarFunctions< Dim > scalars( _r );
	Eigen::VectorXd dot( functionNum() );
	Range< Dim > eRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) eRange.first[d] = 0 , eRange.second[d] = _r;

	for( size_t i=0 ; i<functionNum() ; i++ )
	{
		dot[i] = 0;
		Index< Dim > i1 = scalars.functionIndex( _indices[i].first ) , i2 = scalars.functionIndex( _indices[i].second );

		auto f = [&]( Index< Dim > e )
		{
			const SkewSymmetricMatrix< double , Dim > &w = skew[ scalars.elementIndex(e) ];
			double integral = 0;
			ElementFunction< Dim > f1 = ElementFunction< Dim >( e , i1 ) , f2 = ElementFunction< Dim >( e , i2 );
			ElementWedge< Dim > wedge = ElementVector< Dim >( f1 ) ^ ElementVector< Dim >( f2 );
			for( unsigned int d=0 ; d<ElementWedge< Dim >::Dim ; d++ ) for( unsigned int ii=0 ; ii<2 ; ii++ ) dot[i] += w[d] * wedge(ii,d).integral( _r );
		};
		// Iterate over all cells supported by the i-th function
		Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( i1 ) , Basis< Dim >::FunctionSupport( i2 ) , eRange ).process( f );
	}
	return dot;
}
