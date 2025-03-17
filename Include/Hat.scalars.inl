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

/////////////////////////////////////////
// ScalarFunctions::IntegrationStencil //
/////////////////////////////////////////
template< unsigned int Dim >
template< typename T , unsigned int Rank , unsigned int Radius >
template< typename ... F >
T ScalarFunctions< Dim >::IntegrationStencil< T , Rank , Radius >::operator()( Index< Dim > e , F ... f ) const
{
	static_assert( sizeof...(f)==Rank , "[ERROR] Index count does not match rank" );
	Index< Rank > idx( Point< int , Rank >( Basis< Dim >::template _RelativeIndex< Radius >( e , f )... ) );
	return SquareStencil< T , Rank , StencilSize< Radius >() >::operator()( idx );
}

template< unsigned int Dim >
template< typename T , unsigned int Rank , unsigned int Radius >
template< typename ... F >
T ScalarFunctions< Dim >::IntegrationStencil< T , Rank , Radius >::operator()( Range< Dim > eRange , F ... f ) const
{
	static_assert( sizeof...(f)==Rank , "[ERROR] Index count does not match rank" );
	T value{};
	Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( f ).dilate( Radius ) ... , eRange ).process( [&]( Index< Dim > e ){ value += this->operator()( e , f... ); } );
	return value;
}

///////////////////////////////////////
// ScalarFunctions::StencilCaseTable //
///////////////////////////////////////
template< unsigned int Dim >
template< unsigned int Radius >
Hat::Range< Dim > ScalarFunctions< Dim >::StencilCaseTable< Radius >::Range( unsigned int res )
{
	Hat::Range< Dim > range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = std::min< unsigned int >( Width , res+1 );
	return range;
}

template< unsigned int Dim >
template< unsigned int Radius >
size_t ScalarFunctions< Dim >::StencilCaseTable< Radius >::Size( Hat::Index< Dim > C , unsigned int res )
{
	size_t sz = 1;
	if( res>Width ) for( unsigned int d=0 ; d<Dim ; d++ ) if( C[d]==Radius ) sz *= ( res+1 - 2*Radius );
	return sz;
}

template< unsigned int Dim >
template< unsigned int Radius >
size_t ScalarFunctions< Dim >::StencilCaseTable< Radius >::SubIndex( Hat::Index< Dim > I , unsigned int res )
{
	size_t idx = 0;
	if( res>Width )
	{
		Index< Dim > C = IndexToCase( I , res );
		unsigned int _subWidth = 1;
		for( unsigned int d=0 ; d<Dim ; d++ ) if( C[d]==Radius )
		{
			idx += ( I[d]-Radius ) * _subWidth;
			_subWidth *= res+1 - 2*Radius;
		}
	}
	return idx;
}

template< unsigned int Dim >
template< unsigned int Radius >
Index< Dim > ScalarFunctions< Dim >::StencilCaseTable< Radius >::CaseToIndex( Index< Dim > C , unsigned int res )
{
	if( res>Width )
	{
		for( unsigned int d=0 ; d<Dim ; d++ )
		{
			if     ( C[d]<(int)Radius );
			else if( C[d]>(int)Radius ) C[d] = C[d] - Radius + (res-Radius);
			else                        C[d] = res/2;
		}
	}
	return C;
}

template< unsigned int Dim >
template< unsigned int Radius >
Index< Dim > ScalarFunctions< Dim >::StencilCaseTable< Radius >::IndexToCase( Index< Dim > I , unsigned int res )
{
	if( res>Width )
	{
		for( unsigned int d=0 ; d<Dim ; d++ )
		{
			if     ( I[d]<(int)(    Radius) );
			else if( I[d]>(int)(res-Radius) ) I[d] = Radius + I[d] - (res-Radius);
			else                              I[d] = Radius;
		}
	}
	return I;
}

/////////////////////////////////
// ScalarFunctions::MatrixInfo //
/////////////////////////////////
template< unsigned int Dim >
template< unsigned int Radius , bool Sym >
ScalarFunctions< Dim >::MatrixInfo< Radius , Sym >::MatrixInfo( unsigned int res ) : _caseRange( StencilCaseTable< Radius >::Range( res ) ) , _res(res)
{
	for( unsigned int d=0 ; d<Dim ; d++ ) _fRange.second[d] = _res+1 , _off[d] = Radius , _infoRange.second[d] = 2*Radius+1;
	size_t offset = 0;

	auto f = [&]( Index< Dim > C )
	{
		// Get an index corresponding to the case
		Index< Dim > f1 = StencilCaseTable< Radius >::CaseToIndex( C , res );

		// The information for this case
		_Info &info = _info(C);
		info.sizes[0] = info.sizes[1] = 0;
		info.offset = offset;

		// First clear
		_infoRange.process( [&]( Index< Dim > I ){ info.stencil(I) = -1; } );

		// Then set
		{
			auto f = [&]( Index< Dim > f2 )
			{
				if( f1<f2 || ( f1==f2 && Sym ) ) info.stencil( f2 - f1 + _off ) = info.sizes[0]++;
				if( f1>f2 ) info.sizes[1]++;
			};
			Range< Dim >::Intersect( Range< Dim >( f1 ).dilate( Radius ) , _fRange ).process( f );
			offset += StencilCaseTable< Radius >::Size( C , _res ) * info.sizes[0];
		}
	};
	// Iterate over the different cases and set the associated information
	_caseRange.process( f );
}

template< unsigned int Dim >
template< unsigned int Radius , bool Sym >
size_t ScalarFunctions< Dim >::MatrixInfo< Radius , Sym >::entries( bool all ) const
{
	size_t e = 0;

	auto f = [&]( Index< Dim > C )
	{
		e += _info(C).sizes[0] * StencilCaseTable< Radius >::Size( C , _res );
		if( all ) e += _info(C).sizes[1] * StencilCaseTable< Radius >::Size( C , _res );
	};
	_caseRange.process( f );

	return e;
}

template< unsigned int Dim >
template< unsigned int Radius , bool Sym >
size_t ScalarFunctions< Dim >::MatrixInfo< Radius , Sym >::entry( Index< Dim > F1 , Index< Dim > F2 ) const
{
	const _Info &info = _info( StencilCaseTable< Radius >::IndexToCase( F1 , _res ) );
	Index< Dim > I = F2 - F1 + _off;
	if( !_infoRange.contains(I) ) ERROR_OUT( "Bad index pair: " , F1 , " : " , F2 );
	size_t e = info.stencil(I);
	if( e==-1 ) ERROR_OUT( "No entry: " , F1 , " : " , F2 );
	return info.offset + StencilCaseTable< Radius >::SubIndex( F1 , _res ) * info.sizes[0] + e;
}

template< unsigned int Dim >
template< unsigned int Radius , bool Sym >
template< typename F /* = std::function< void ( Index< Dim > , size_t ) > */ >
void ScalarFunctions< Dim >::MatrixInfo< Radius , Sym >::process( Index< Dim > f1 , F f ) const
{
	const _Info &info = _info( StencilCaseTable< Radius >::IndexToCase( f1 , _res ) );
	auto _f = [&]( Index< Dim > I )
	{
		size_t e = info.stencil(I);
		if( e!=-1 )
		{
			Index< Dim > f2 = f1 + I - _off;
			if( !_fRange.contains( f2 ) ) ERROR_OUT( "should not be happening" );
			f( f2 , info.offset + StencilCaseTable< Radius >::SubIndex( f1 , _res ) * info.sizes[0] + e );
		}
	};
	_infoRange.process( _f );
}

template< unsigned int Dim >
template< unsigned int Radius , bool Sym >
template< typename F /* = std::function< void ( Index< Dim > , size_t , bool ) > */ >
void ScalarFunctions< Dim >::MatrixInfo< Radius , Sym >::processAll( Index< Dim > f1 , F f ) const
{
	const _Info &info = _info( StencilCaseTable< Radius >::IndexToCase( f1 , _res ) );

	auto _f = [&]( Index< Dim > I )
	{
		size_t e = info.stencil(I);
		Index< Dim > f2 = f1 + I - _off;
		// Confirm that this is an in-bounds function index
		if( !_fRange.contains( f2 ) ) return;
		if( e!=-1 ) f( f2 , info.offset + StencilCaseTable< Radius >::SubIndex( f1 , _res ) * info.sizes[0] + e , true );
		else
		{
			// I = (f2-f1) + _off
			// Flipping the roles of f1 and f2
			// I <- (f1-f2) + _off = - ( I - _off ) + _off
			const _Info &info = _info( StencilCaseTable< Radius >::IndexToCase( f2 , _res ) );
			// The index of f1 relative to f2
			I = - ( I - _off ) + _off;
			e = info.stencil(I);
			if( e!=-1 ) f( f2 , info.offset + StencilCaseTable< Radius >::SubIndex( f2 , _res ) * info.sizes[0] + e , false );
		}
	};
	// Process _all_ functions within the prescribed radius
	_infoRange.process( _f );
}

template< unsigned int Dim >
template< unsigned int Radius , bool Sym >
size_t ScalarFunctions< Dim >::MatrixInfo< Radius , Sym >::entries( Index< Dim > I , bool all ) const
{
	if( all ) return _info( StencilCaseTable< Radius >::IndexToCase( I , _res ) ).sizes[0] + _info( StencilCaseTable< Radius >::IndexToCase( I , _res ) ).sizes[1];
	else      return _info( StencilCaseTable< Radius >::IndexToCase( I , _res ) ).sizes[0];
}

/////////////////////////////////////////////
// ScalarFunctions::FullIntegrationStencil //
/////////////////////////////////////////////
template< unsigned int Dim >
template< typename T , unsigned int Radius >
constexpr unsigned int ScalarFunctions< Dim >::FullIntegrationStencil< T , Radius >::StencilNum( void )
{
	constexpr unsigned int Size = 3 + 4*Radius;
	if constexpr( Dim==1 ) return Size;
	else return ScalarFunctions< Dim-1 >::template FullIntegrationStencil< T , Radius >::StencilNum() * Size;
}

template< unsigned int Dim >
template< typename T , unsigned int Radius>
ScalarFunctions< Dim >::FullIntegrationStencil< T , Radius >::FullIntegrationStencil( const IntegrationStencil< T , 2 , Radius > &stencil , unsigned int res ) : _res( res ) , _rows( StencilNum() )
{
	Range< Dim > eRange , fRange , range;
	for( unsigned int d=0 ; d<Dim ; d++ ) eRange.first[d] = 0 , eRange.second[d] = _res , fRange.first[d] = 0 , fRange.second[d] = _res+1 , range.first[d] = 0 , range.second[d] = 3+4*Radius;

	auto f = [&]( Index< Dim > f1 )
	{
		// [0,2*Radius]         <-> [0,2*Radius]
		// [_res/2]             <-> [2*Radius+1]
		// [_res-2*Radius,_res] <-> [2*Radius+2,4*Radius+2]
		for( unsigned int d=0 ; d<Dim ; d++ )
			if     ( f1[d]<=2*Radius   );
			else if( f1[d]==2*Radius+1 ) f1[d] = _res/2;
			else                         f1[d] = _res - 4*Radius-2 + f1[d];
		Row _row;
		std::map< Index< Dim > , T > row;

		auto f = [&]( Index< Dim > e )
		{
			auto f = [&]( Index< Dim > f2 ){ row[ f2-f1 ] += stencil( e , f1 , f2 ); };
			Range< Dim >::Intersect( Basis< Dim >::ElementSupport( e ).dilate( Radius ) , fRange ).process(f);
//			Basis< Dim >::ElementSupport( e ).dilate( Radius ).process(f);
		};
		Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( f1 ).dilate( Radius ) , eRange ).process( f );

		_row.reserve( row.size() );
		for( auto iter=row.begin() ; iter!=row.end() ; iter++ ) _row.push_back( Entry( iter->first , iter->second ) );

		_rows[ StencilIndex(f1,_res) ] = _row;
	};
	Range< Dim >::Intersect( range , fRange ).process( f );
}

template< unsigned int Dim >
template< typename T , unsigned int Radius >
unsigned int ScalarFunctions< Dim >::FullIntegrationStencil< T , Radius >::StencilIndex( Index< Dim > f , unsigned int res )
{
	unsigned int idx=0;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		if( f[d]<0 || f[d]>(int)res ) ERROR_OUT( "Bad function index: " , f , " / " , res );
		if     ( f[d]<=(int)(    2*Radius) ) idx = idx*(3+4*Radius) + f[d];
		else if( f[d]>=(int)(res-2*Radius) ) idx = idx*(3+4*Radius) + 2*Radius+2+f[d]-(res-2*Radius);
		else                                 idx = idx*(3+4*Radius) + 2*Radius+1;
	}
	return idx;
}

template< unsigned int Dim >
template< typename T , unsigned int Radius >
const typename ScalarFunctions< Dim >::template FullIntegrationStencil< T , Radius >::Row &ScalarFunctions< Dim >::FullIntegrationStencil< T , Radius >::row( Index< Dim > f ) const
{
	return _rows[ StencilIndex( f , _res ) ];
}

template< unsigned int Dim >
template< typename T >
typename ScalarFunctions< Dim >::template FullIntegrationStencil< T , 0 > ScalarFunctions< Dim >::Restrict( const FullIntegrationStencil< T , 0 > &stencilF )
{
	if( stencilF._res&1 ) ERROR_OUT( "Cannot restrict when reslution is odd" );
	FullIntegrationStencil< T , 0 > stencilC;
	stencilC._res = stencilF._res/2;
	stencilC._rows.resize( FullIntegrationStencil< T , 0 >::StencilNum() );


	Range< Dim > fRange , range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = 3 , fRange.second[d] = stencilC._res+1;

	auto f = [&]( Index< Dim > F_coarse1 )
		{
			// Transform the case-index to a coarse index
			// [0]      <-> [0]
			// [_res/2] <-> [1]
			// [_res]   <-> [2]
			for( unsigned int d=0 ; d<Dim ; d++ )
				if     ( F_coarse1[d]==0 );
				else if( F_coarse1[d]==1 ) F_coarse1[d] = stencilC._res/2;
				else                       F_coarse1[d] = stencilC._res;

			std::map< Index< Dim > , T > row;

			auto f = [&]( Index< Dim > F_fine1 , double weight1 )
				{
					// Iterate over all stencil entries supported on F_fine1
					const typename FullIntegrationStencil< T , 0 >::Row &rowF = stencilF.row( F_fine1 );
					for( unsigned int i=0 ; i<rowF.size() ; i++ )
					{
						Index< Dim > F_fine2 = F_fine1 + rowF[i].first;
						T value = rowF[i].second * weight1;

						auto f = [&]( Index< Dim > F_coarse2 , double weight2 ){ row[ F_coarse2-F_coarse1 ] += value * weight2; };
						// Iterate over all the coarser nodes F_fine2 restricts to
						Restriction::Process( F_fine2 , f );
					}
				};
			// Iterate over all the finer nodes F_coarse1 prolongs to
			Prolongation::Process( F_coarse1 , stencilC._res , f );

			typename FullIntegrationStencil< T , 0 >::Row &_row = stencilC._rows[ FullIntegrationStencil< T , 0 >::StencilIndex( F_coarse1 , stencilC._res ) ];

			_row.reserve( row.size() );
			for( auto iter=row.begin() ; iter!=row.end() ; iter++ ) _row.push_back( typename FullIntegrationStencil< T , 0 >::Entry( iter->first , iter->second ) );
		};
	// Iterate over the cases
	Range< Dim >::Intersect( range , fRange ).process( f );
	return stencilC;
}

//////////////////////////////////////////
// ScalarFunctions::ProlongationStencil //
//////////////////////////////////////////
template< unsigned int Dim >
ScalarFunctions< Dim >::ProlongationStencil::ProlongationStencil( void )
{
	Index< Dim > off;
	for( unsigned int d=0 ; d<Dim ; d++ ) _range.first[d] = -1 , _range.second[d] = 2 , off[d] = 1;

	auto f = [&]( Index< Dim > I )
	{
		double value = 1;
		for( unsigned int d=0 ; d<Dim; d++ ) if( I[d]!=0 ) value *= 0.5;
		SquareStencil< double , Dim , StencilCaseTable<1>::Width >::operator()( I+off ) = value;
	};
	_range.process( f );
}

//////////////////
// Prolongation //
//////////////////
template< unsigned int Dim >
template< typename ProlongationFunctor/*=std::function< void ( Index< Dim > F_fine , double weight ) >*/ >
void ScalarFunctions< Dim >::Prolongation::Process( Index< Dim > F_coarse , unsigned int coarseRes , ProlongationFunctor pFunctor )
{
	Range< Dim > fRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) fRange.second[d] = 2*coarseRes+1;
	auto f = [&]( Index< Dim > F_fine )
		{
			double scale = 1.;
			for( unsigned int d=0 ; d<Dim ; d++ ) if( F_fine[d]!=(F_coarse[d]*2) ) scale *= 0.5;
			pFunctor( F_fine , scale );
		};
	Range< Dim >::Intersect( Range< Dim >( F_coarse*2 ).dilate(1) , fRange ).process( f );
}

/////////////////
// Restriction //
/////////////////
template< unsigned int Dim >
template< typename RestrictionFunctor/*=std::function< void ( Index< Dim > F_coarse , double weight ) >*/ >
void ScalarFunctions< Dim >::Restriction::Process( Index< Dim > F_fine , RestrictionFunctor rFunctor )
{
	Index< Dim > F_coarse;
	_Process< 0 >( F_fine , F_coarse , 1. , rFunctor );
}

template< unsigned int Dim >
template< unsigned int D , typename RestrictionFunctor/*=std::function< void ( Index< Dim > F_coarse , double weight ) >*/ >
void ScalarFunctions< Dim >::Restriction::_Process( Index< Dim > F_fine , Index< Dim > F_coarse , double weight , RestrictionFunctor &rFunctor )
{
	if constexpr( D==Dim ) rFunctor( F_coarse , weight );
	else
	{
		F_coarse[D] = F_fine[D]>>1;
		if( F_fine[D]&1 )
		{
			weight *= 0.5;
			_Process< D+1 >( F_fine , F_coarse , weight , rFunctor );
			F_coarse[D]++;
			_Process< D+1 >( F_fine , F_coarse , weight , rFunctor );
		}
		else _Process< D+1 >( F_fine , F_coarse , weight , rFunctor );
	}
}

/////////////////////
// ScalarFunctions //
/////////////////////
template< unsigned int Dim >
size_t ScalarFunctions< Dim >::FunctionNum( unsigned int res )
{
	size_t sz = 1;
	for( unsigned int d=0 ; d<Dim ; d++ ) sz *= res+1;
	return sz;
}

template< unsigned int Dim >
size_t ScalarFunctions< Dim >::ElementNum( unsigned int res )
{
	size_t sz = 1;
	for( unsigned int d=0 ; d<Dim ; d++ ) sz *= res;
	return sz;
}

template< unsigned int Dim >
size_t ScalarFunctions< Dim >::ElementIndex( Index< Dim > i , unsigned int res )
{
	size_t idx=0;
	for( unsigned int d=0 ; d<Dim ; d++ ) idx = idx*res + i[Dim-d-1];
	return idx;
}

template< unsigned int Dim >
Index< Dim > ScalarFunctions< Dim >::ElementIndex( size_t i , unsigned int res )
{
	Index< Dim > idx;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		idx[d] = i % res;
		i /= res;
	}
	return idx;
}

template< unsigned int Dim >
size_t ScalarFunctions< Dim >::FunctionIndex( Index< Dim > i , unsigned int res )
{
	size_t idx=0;
	for( unsigned int d=0 ; d<Dim ; d++ ) idx = idx*(res+1) + i[Dim-d-1];
	return idx;
}

template< unsigned int Dim >
Index< Dim > ScalarFunctions< Dim >::FunctionIndex( size_t i , unsigned int res )
{
	Index< Dim > idx;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		idx[d] = i % (res+1);
		i /= res+1;
	}
	return idx;
}

template< unsigned int Dim >
Range< Dim > ScalarFunctions< Dim >::ElementRange( unsigned int res )
{
	Range< Dim > eRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) eRange.first[d] = 0 , eRange.second[d] = res;
	return eRange;
}

template< unsigned int Dim >
Range< Dim > ScalarFunctions< Dim >::FunctionRange( unsigned int res )
{
	Range< Dim > fRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) fRange.first[d] = 0 , fRange.second[d] = res+1;
	return fRange;
}



template< unsigned int Dim >
typename ScalarFunctions< Dim >::template IntegrationStencil< double , 2 , 0 > ScalarFunctions< Dim >::MassStencil( unsigned int res )
{
	ScalarFunctions< Dim >::IntegrationStencil< double , 2 , 0 > stencil;

	Index< Dim > e;
	Range< Dim > r = Basis< Dim >::ElementSupport( e );

	auto f = [&]( Index< Dim > i1 , Index< Dim > i2 )
	{
		stencil._values[ Basis< Dim >::template _RelativeIndex<0>(e,i1) ][ Basis< Dim >::template _RelativeIndex<0>(e,i2) ] = ( ElementFunction< Dim >( e , i1 ) * ElementFunction< Dim >( e , i2 ) ).integral(res);
	};
	r.template process< 2 >( f );

	return stencil;
}

template< unsigned int Dim >
typename ScalarFunctions< Dim >::template IntegrationStencil< double , 2 , 0 > ScalarFunctions< Dim >::StiffnessStencil( unsigned int res )
{
	ScalarFunctions< Dim >::IntegrationStencil< double , 2 , 0 > stencil;

	Index< Dim > e;

	auto f = [&]( Index< Dim > i1 , Index< Dim > i2 )
	{
		ElementVector< Dim > v1 = ElementFunction< Dim >(e,i1).gradient() , v2 = ElementFunction< Dim >(e,i2).gradient();
		double integral = 0;
		for( unsigned int d=0 ; d<Dim ; d++ ) integral +=( v1[d] * v2[d] ).integral(res);
		stencil._values[ Basis< Dim >::template _RelativeIndex<0>(e,i1) ][ Basis< Dim >::template _RelativeIndex<0>(e,i2) ] = integral;
	};
	Basis< Dim >::ElementSupport( e ).template process< 2 >( f );

	return stencil;
}


template< unsigned int Dim >
typename ScalarFunctions< Dim >::template IntegrationStencil< MishaK::SquareMatrix< double , Dim > , 2 , 0 > ScalarFunctions< Dim >::FullStiffnessStencil( unsigned int res )
{
	ScalarFunctions< Dim >::IntegrationStencil< MishaK::SquareMatrix< double , Dim > , 2 , 0 > stencil;

	Index< Dim > e;
	auto f = [&]( Index< Dim > i1 , Index< Dim > i2 )
	{
		ElementVector< Dim > v1 = ElementFunction< Dim >(e,i1).gradient() , v2 = ElementFunction< Dim >(e,i2).gradient();
		MishaK::SquareMatrix< double , Dim > &m = stencil._values[ Basis< Dim >::template _RelativeIndex<0>(e,i1) ][ Basis< Dim >::template _RelativeIndex<0>(e,i2) ];
		for( unsigned int d1=0 ; d1<Dim ; d1++ ) for( unsigned int d2=0 ; d2<Dim ; d2++ ) m(d2,d1) += ( v1[d1] * v2[d2] ).integral(res);
	};
	Basis< Dim >::ElementSupport( e ).template process< 2 >( f );

	return stencil;
}

template< unsigned int Dim >
typename ScalarFunctions< Dim >::template IntegrationStencil< double , 1 , 0 > ScalarFunctions< Dim >::ValueStencil( unsigned int res )
{
	ScalarFunctions< Dim >::IntegrationStencil< double , 1 , 0 > stencil;
	Index< Dim > e;

	auto f = [&]( Index< Dim > f )
	{
		ElementFunction< Dim > v(e,f);
		stencil._values[ Basis< Dim >::template _RelativeIndex<0>(e,f) ] = v.integral(res);
	};
	Basis< Dim >::ElementSupport( e ).process( f );
	return stencil;
}


template< unsigned int Dim >
typename ScalarFunctions< Dim >::template IntegrationStencil< Point< double , Dim > , 1 , 0 > ScalarFunctions< Dim >::PartialDerivativeStencil( unsigned int res )
{
	ScalarFunctions< Dim >::IntegrationStencil< Point< double , Dim > , 1 , 0 > stencil;
	Index< Dim > e;

	auto f = [&]( Index< Dim > f )
	{
		ElementVector< Dim > v = ElementFunction< Dim >(e,f).gradient();
		Point< double , Dim > &_v = stencil._values[ Basis< Dim >::template _RelativeIndex<0>(e,f) ];
		for( unsigned int d=0 ; d<Dim ; d++ ) _v[d] += v[d].integral(res);
	};
	Basis< Dim >::ElementSupport( e ).process( f );
	return stencil;
}

template< unsigned int Dim >
template< typename T , typename Indexer /* = Hat::BaseIndexer< Dim > */ >
T ScalarFunctions< Dim >::operator()( const Indexer & indexer , const IntegrationStencil< T , 2 , 0 > &stencil , const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	std::vector< double > integrals( ThreadPool::NumThreads() , 0 );
	ThreadPool::ParallelFor
	(
		0 , indexer.numElements() ,
		[&]( unsigned int t , size_t e )
		{
			Index< Dim > E = indexer.elementIndex( e );

			Window::IsotropicStaticWindow< size_t , Dim , 2 > neighbors = indexer.efNeighbors( e , t );

			for( unsigned int i1=0 ; i1<Window::IsotropicSize< Dim , 2 >() ; i1++ ) if( neighbors.neighbors.data[i1]!=-1 )
			{
				unsigned int f1 = neighbors.neighbors.data[i1];
				Index< Dim > F1 = indexer.functionIndex( f1 );
				for( unsigned int i2=0 ; i2<Window::IsotropicSize< Dim , 2 >() ; i2++ ) if( neighbors.neighbors.data[i2]!=-1 )
				{
					unsigned int f2 = neighbors.neighbors.data[i2];
					Index< Dim > F2 = indexer.functionIndex( f2 );
					integrals[t] += stencil( E , F1 , F2 ) * x[f1] * y[f2];
				}
			}
		}
	);


	T integral{};
	for( unsigned int i=0 ; i<integrals.size() ; i++ ) integral += integrals[i];
	return integral;
}

template< unsigned int Dim >
template< typename T , typename Indexer /* = Hat::BaseIndexer< Dim > */ >
T ScalarFunctions< Dim >::operator()( const Indexer & indexer , const FullIntegrationStencil< T , 0 > &stencil , const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	Index< Dim > Off;
	for( unsigned int d=0 ; d<Dim ; d++ ) Off[d] = 1;

	std::vector< double > integrals( ThreadPool::NumThreads() , 0 );
	ThreadPool::ParallelFor
	(
		0 , indexer.numFunctions() ,
		[&]( unsigned int t , size_t f1 )
		{
			Index< Dim > F1 = indexer.functionIndex( f1 );
			Window::IsotropicStaticWindow< size_t , Dim , 3 > neighbors = indexer.ffNeighbors( f1 , t );
			const typename Hat::ScalarFunctions< Dim >::FullIntegrationStencil< T , 0 >::Row &row = stencil.row( F1 );
			for( unsigned int j=0 ; j<row.size() ; j++ )
			{
				Index< Dim > I = row[j].first + Off;
				size_t f2 = neighbors( &I[0] );
				if( f2!=-1 ) integrals[t] += row[j].second * x[f1] * y[f2];
			}
		}
	);

	T integral{};
	for( unsigned int i=0 ; i<integrals.size() ; i++ ) integral += integrals[i];
	return integral;
}


template< unsigned int Dim >
template< typename T , typename Indexer /* = Hat::BaseIndexer< Dim > */ >
Eigen::Matrix< T , 1 , Eigen::Dynamic > ScalarFunctions< Dim >::operator()( const Indexer & indexer , const IntegrationStencil< T , 2 , 0 > &stencil , const Eigen::VectorXd &x ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	Eigen::Matrix< T , 1 , Eigen::Dynamic > v = Eigen::Matrix< T , 1 , Eigen::Dynamic >::Zero( indexer.numFunctions() );

	for( int e=0 ; e<indexer.numElements() ; e++ )
	{
		Index< Dim > E = indexer.elementIndex( e );

		Window::IsotropicStaticWindow< size_t , Dim , 2 > neighbors = indexer.efNeighbors( e , 0 );

		for( unsigned int i1=0 ; i1<Window::IsotropicSize< Dim , 2 >() ; i1++ ) if( neighbors.neighbors.data[i1]!=-1 )
		{
			unsigned int f1 = neighbors.neighbors.data[i1];
			Index< Dim > F1 = indexer.functionIndex( f1 );
			for( unsigned int i2=0 ; i2<Window::IsotropicSize< Dim , 2 >() ; i2++ ) if( neighbors.neighbors.data[i2]!=-1 )
			{
				unsigned int f2 = neighbors.neighbors.data[i2];
				Index< Dim > F2 = indexer.functionIndex( f2 );
				v[f2] += stencil( E , f1 , f2 ) * x[f1];
			}
		}
	}

	return v;
}

template< unsigned int Dim >
template< typename T , typename Indexer /* = Hat::BaseIndexer< Dim > */ >
Eigen::Matrix< T , 1 , Eigen::Dynamic > ScalarFunctions< Dim >::operator()( const Indexer & indexer , const FullIntegrationStencil< T , 0 > &stencil , const Eigen::VectorXd &x ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	Index< Dim > Off;
	for( unsigned int d=0 ; d<Dim ; d++ ) Off[d] = 1;

	Eigen::Matrix< T , 1 , Eigen::Dynamic > v = Eigen::Matrix< T , 1 , Eigen::Dynamic >::Zero( indexer.numFunctions() );

	for( int f1=0 ; f1<indexer.numFunctions() ; f1++ )
	{
		Index< Dim > F1 = indexer.functionIndex( f1 );

		Window::IsotropicStaticWindow< size_t , Dim , 3 > neighbors = indexer.ffNeighbors( f1 , 0 );

		const typename Hat::ScalarFunctions< Dim >::FullIntegrationStencil< T , 0 >::Row &row = stencil.row( F1 );

		for( unsigned int j=0 ; j<row.size() ; j++ )
		{
			Index< Dim > I = row[j].first + Off;
			size_t f2 = neighbors( &I[0] );
			if( f2!=-1 ) v[f2] += row[j].second * x[f1];
		}
	}
	return v;
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::systemMatrix( const Indexer & indexer , IntegrationStencil< double , 2 , 0 > stencil ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	static const unsigned int Radius = 0;
	Hat::Index< Dim > Off;
	for( unsigned int d=0 ; d<Dim ; d++ ) Off[d] = 1;

	// If a function is supported on a radius of (1+Radius) cells, then two functions have overlapping support if
	// they are in each other's (1+2*Radius)-ring
	MatrixInfo< 1+2*Radius , true > matrixInfo( _r );
	Eigen::SparseMatrix< double > M( indexer.numFunctions() , indexer.numFunctions() );
	Eigen::VectorXi rowSizes( (int)indexer.numFunctions() );
	ThreadPool::ParallelFor( 0 , indexer.numFunctions() , [&]( unsigned int , size_t i ){ Index< Dim > I = indexer.functionIndex( (size_t)i ) ; rowSizes[i] = (int)matrixInfo.entries(I,true); } );
	M.reserve( rowSizes );

	ThreadPool::ParallelFor
	(
		0 , indexer.numFunctions() ,
		[&]( unsigned int t , size_t f1 )
		{
			Window::IsotropicStaticWindow< size_t , Dim , 3 > neighbors = indexer.ffNeighbors( f1 , 0 );
			Index< Dim > F1 = indexer.functionIndex( f1 );

			auto Kernel = [&]( Index< Dim > F2 , size_t e , bool flip )
				{
					Hat::Index< Dim > I = F2-F1+Off;
					size_t f2 = neighbors( &I[0] );
					if( f2!=-1 )
					{
						Range< Dim > range = Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( F1 ) , Basis< Dim >::FunctionSupport( F2 ) , _eRange );
						M.insert( (int)f2 , (int)f1 ) = stencil( range , F1 , F2 );
					}
				};
			matrixInfo.processAll( F1 , Kernel );
		}
	);

	M.makeCompressed();
	return M;
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::systemMatrix( const Indexer & indexer , FullIntegrationStencil< double , 0 > stencil ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	static const unsigned int Radius = 0;
	Hat::Index< Dim > Off;
	for( unsigned int d=0 ; d<Dim ; d++ ) Off[d] = 1;

	Eigen::SparseMatrix< double > M( indexer.numFunctions() , indexer.numFunctions() );
	Eigen::VectorXi rowSizes( (int)indexer.numFunctions() );
	ThreadPool::ParallelFor( 0 , indexer.numFunctions() , [&]( unsigned int , size_t i ){ rowSizes[i] = (int)stencil.row( indexer.functionIndex( (size_t)i ) ).size(); } );

	M.reserve( rowSizes );

	ThreadPool::ParallelFor
	(
		0 , indexer.numFunctions() ,
		[&]( unsigned int , size_t f1 )
		{
			Window::IsotropicStaticWindow< size_t , Dim , 3 > neighbors = indexer.ffNeighbors( f1 , 0 );
			Index< Dim > F1 = indexer.functionIndex( f1 );

			const typename Hat::ScalarFunctions< Dim >::FullIntegrationStencil< double , Radius >::Row &row = stencil.row( F1 );
			for( unsigned int j=0 ; j<row.size() ; j++ )
			{
				Hat::Index< Dim > I = row[j].first + Off;
				size_t f2 = neighbors( &I[0] );
				if( f2!=-1 ) M.insert( (int)f2 , (int)f1 ) = row[j].second;
			}
		}
	);

	M.makeCompressed();
	return M;
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::boundarySystemMatrix( const Indexer & indexer , typename ScalarFunctions< Dim-1 >::template IntegrationStencil< double , 2 , 0 > stencil ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	static const unsigned int Radius = 0;
	Hat::Index< Dim > Off;
	for( unsigned int d=0 ; d<Dim ; d++ ) Off[d] = 1;

	if constexpr( Dim==1 )
	{
		Eigen::SparseMatrix< double > M( indexer.numFunctions() , indexer.numFunctions() );
		return M;
	}
	else
	{
		auto FullIndex = [&]( Index< Dim-1 > _F , unsigned int d , unsigned int offset )
			{
				Index< Dim > F;
				for( unsigned int _d=0 ; _d<d ; _d++ ) F[_d] = _F[_d];
				F[d] = offset;
				for( unsigned int _d=d+1 ; _d<Dim ; _d++ ) F[_d] = _F[_d-1];
				return F;
			};

		// Not parallelized, but handles multiple inserts of the same coefficients
		ScalarFunctions< Dim-1 > _scalars( _r );
		typename ScalarFunctions< Dim-1 >::template FullIntegrationStencil< double , Radius > fullStencil( stencil , _r );

		std::vector< Eigen::Triplet< double > > triplets;

		auto AddEntries = [&]( size_t f1 , Hat::Index< Dim > F1 , unsigned int d , unsigned int offset , Window::IsotropicStaticWindow< size_t , Dim , 3 > &neighbors )
			{
				Hat::Index< Dim > _F1;
				for( unsigned int dd=0 , _d=0 ; dd<Dim ; dd++ ) if( dd!=d ) _F1[_d++] = F1[dd];

				const typename Hat::ScalarFunctions< Dim-1 >::template FullIntegrationStencil< double , Radius >::Row &row = fullStencil.row( _F1 );
				for( unsigned int j=0 ; j<row.size() ; j++ )
				{
					Index< Dim-1 > _F2 = _F1 + row[j].first;
					Index< Dim > F2 = FullIndex( _F2 , d , offset );
					Index< Dim > I = F2 - F1 + Off;
					size_t f2 = neighbors( &I[0] );
					if( f2!=-1 )
					{
						if( f1<=f2 )
						{
							double integral = row[j].second;
							if( f1==f2 ) integral /= 2;
							for( unsigned int d=0 ; d<Dim ; d++ ) triplets.push_back(Eigen::Triplet< double >( (int)f1 , (int)f2 , integral));
						}
					}
				}
			};

		for( size_t f1=0 ; f1<indexer.numFunctions() ; f1++ )
		{
			Hat::Index< Dim > F1 = indexer.functionIndex( f1 );

			bool onBoundary = false;
			for( unsigned int d=0 ; d<Dim ; d++ ) onBoundary |= F1[d]==0 || F1[d]==_r;

			if( onBoundary )
			{
				Window::IsotropicStaticWindow< size_t , Dim , 3 > neighbors = indexer.ffNeighbors( f1 , 0 );

				for( unsigned int d=0 ; d<Dim ; d++ )
				{
					if( F1[d]== 0 ) AddEntries( f1 , F1 , d ,  0 , neighbors );
					if( F1[d]==_r ) AddEntries( f1 , F1 , d , _r , neighbors );
				}
			}
		}

		{
			Eigen::SparseMatrix< double > M( indexer.numFunctions() , indexer.numFunctions() );
			M.setFromTriplets(triplets.begin(), triplets.end());
			triplets.clear();
			Eigen::SparseMatrix< double > Mt = M.transpose();
			return M + Mt;
		}
	}
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ , typename SampleFunctor /* = std::function< Point< double , Dim > ( unsigned int idx ) > */ , typename WeightFunctor /* = std::function< double ( Point< double , Dim > ) */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::deltaMass( const Indexer & indexer , SampleFunctor && F , size_t sampleNum , WeightFunctor && wF ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	// The values of the corner functions at the samples
	struct FunctionValues
	{
		double values[1<<Dim];
		FunctionValues( Point< double , Dim > p )
		{
			Index< Dim > e;
			double _values[2][Dim];
			for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) _values[i][d] = i==0 ? (1.-p[d]) : p[d];

			auto f = [&]( Index< Dim > f )
				{
					double value = 1.;
					for( unsigned int d=0 ; d<Dim ; d++ ) value *= f[d]==e[d] ? _values[0][d] : _values[1][d];
					values[ Basis< Dim >::template _RelativeIndex<0>( e , f ) ] = value;
				};
			Basis< Dim >::ElementSupport( e ).process( f );
		}
	};

	MultiThreadedTriplets< double > triplets;

	OrderedSamples< Dim , double > samples( [&]( size_t i ){ Point< double , Dim > p = F(i) ; return std::pair< Point< double , Dim > , double >( p , wF(p) ); } , sampleNum , _r );

	ThreadPool::ParallelFor
		(
			0 , samples.size() ,
			[&]( unsigned int t , size_t i )
			{
				Hat::Index< Dim > E = samples[i].first;
				const std::vector< std::pair< Point< double , Dim > , double > > &subSamples = samples[i].second;

				Point< double , Dim > p;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = ( E[d] + 0.5 ) / _r;
				Window::IsotropicStaticWindow< size_t , Dim , 2 > nbrs = indexer.fNeighbors( p , t );
				// Get the values of each of the function corner functions at each of the samples
				std::vector< FunctionValues > functionValues;
				functionValues.reserve( subSamples.size() );
				for( unsigned int j=0 ; j<subSamples.size() ; j++ )
				{
					Point< double , Dim > p = subSamples[j].first * _r;
					for( unsigned int d=0 ; d<Dim ; d++ ) p[d] -= E[d];
					functionValues.emplace_back( p );
				}

				for( unsigned int i1=0 ; i1<Window::IsotropicSize< Dim , 2 >() ; i1++ )
				{
					size_t idx1 = nbrs.data[i1];
					if( idx1!=-1 )
					{
						Index< Dim > F1 = indexer.functionIndex( idx1 );
						for( unsigned int i2=0 ; i2<Window::IsotropicSize< Dim , 2 >() ; i2++ )
						{
							size_t idx2 = nbrs.data[i2];
							if( idx2!=-1 )
							{
								Index< Dim > F2 = indexer.functionIndex( idx2 );
								if( F1<=F2 )
								{
									double value = 0;
									for( unsigned int j=0 ; j<subSamples.size() ; j++ ) value += functionValues[j].values[ Basis< Dim >::template _RelativeIndex<0>(E,F1) ] * functionValues[j].values[ Basis< Dim >::template _RelativeIndex<0>(E,F2) ] * subSamples[j].second;
									if( F1==F2 ) value /= 2.;
									triplets[t].push_back( Eigen::Triplet< double >( (int)idx1 , (int)idx2 , value ) );
								}
							}
						}
					}
				}
			}
		);

	Eigen::SparseMatrix< double > E( indexer.numFunctions(), indexer.numFunctions());
	E.setFromTriplets( triplets.begin() , triplets.end() );
	Eigen::SparseMatrix< double > Et = E.transpose();
	return E + Et;
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ , typename SampleFunctor /* = std::function< Point< double , Dim > ( unsigned int idx ) > */ , typename WeightFunctor /* = std::function< double ( Point< double , Dim > ) */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::deltaMass( const Indexer & indexer , SampleFunctor && F , const OrderedSampler< Dim > &orderedSampler , WeightFunctor && wF ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );
	static_assert( std::is_convertible_v< SampleFunctor , std::function< Point< double , Dim > ( size_t ) > > , "[ERROR] SampleFunctor poorly formed" );
	static_assert( std::is_convertible_v< WeightFunctor , std::function< double ( Point< double , Dim > , unsigned int ) > > , "[ERROR] WeightFunctor poorly formed" );

	// The values of the corner functions at the samples
	struct FunctionValues
	{
		double values[1<<Dim];
		FunctionValues( Point< double , Dim > p )
		{
			Index< Dim > e;
			double _values[2][Dim];
			for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) _values[i][d] = i==0 ? (1.-p[d]) : p[d];

			auto f = [&]( Index< Dim > f )
				{
					double value = 1.;
					for( unsigned int d=0 ; d<Dim ; d++ ) value *= f[d]==e[d] ? _values[0][d] : _values[1][d];
					values[ Basis< Dim >::template _RelativeIndex<0>( e , f ) ] = value;
				};
			Basis< Dim >::ElementSupport( e ).process( f );
		}
	};

	MultiThreadedTriplets< double > triplets;

	ThreadPool::ParallelFor
	(
		0 , orderedSampler.size() ,
		[&]( unsigned int t , size_t i )
		{
			Hat::Index< Dim > E = orderedSampler[i].first;
			const std::vector< size_t > &subSampleIndices = orderedSampler[i].second;
			std::vector< std::pair< Point< double , Dim > , double > > subSamples( subSampleIndices.size() );
			for( size_t j=0 ; j<subSampleIndices.size() ; j++ )
			{
				subSamples[j].first = F( subSampleIndices[j] );
				subSamples[j].second = wF( subSamples[j].first , t );
			}

			Point< double , Dim > p;
			for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = ( E[d] + 0.5 ) / _r;
			Window::IsotropicStaticWindow< size_t , Dim , 2 > nbrs = indexer.fNeighbors( p , t );
			// Get the values of each of the function corner functions at each of the samples
			std::vector< FunctionValues > functionValues;
			functionValues.reserve( subSamples.size() );
			for( unsigned int j=0 ; j<subSamples.size() ; j++ )
			{
				Point< double , Dim > p = subSamples[j].first * _r;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] -= E[d];
				functionValues.emplace_back( p );
			}

			for( unsigned int i1=0 ; i1<Window::IsotropicSize< Dim , 2 >() ; i1++ )
			{
				size_t idx1 = nbrs.data[i1];
				if( idx1!=-1 )
				{
					Index< Dim > F1 = indexer.functionIndex( idx1 );
					for( unsigned int i2=0 ; i2<Window::IsotropicSize< Dim , 2 >() ; i2++ )
					{
						size_t idx2 = nbrs.data[i2];
						if( idx2!=-1 )
						{
							Index< Dim > F2 = indexer.functionIndex( idx2 );
							if( F1<=F2 )
							{
								double value = 0;
								for( unsigned int j=0 ; j<subSamples.size() ; j++ ) value += functionValues[j].values[ Basis< Dim >::template _RelativeIndex<0>(E,F1) ] * functionValues[j].values[ Basis< Dim >::template _RelativeIndex<0>(E,F2) ] * subSamples[j].second;
								if( F1==F2 ) value /= 2.;
								triplets[t].push_back( Eigen::Triplet< double >( (int)idx1 , (int)idx2 , value ) );
							}
						}
					}
				}
			}
		}
	);

	Eigen::SparseMatrix< double > E( indexer.numFunctions(), indexer.numFunctions());
	E.setFromTriplets( triplets.begin() , triplets.end() );
	Eigen::SparseMatrix< double > Et = E.transpose();
	return E + Et;
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ , typename WeightFunctor /* = std::function< double ( Point< double , Dim > ) */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::deltaMass( const Indexer & indexer , const std::vector< Point< double , Dim > > &samples , WeightFunctor && wF ) const
{
	return deltaMass( std::forward< Indexer >( indexer ) , [&]( unsigned int idx ){ return samples[idx]; } , samples.size() , std::forward< WeightFunctor >( wF ) );
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ , typename SampleFunctor /* = std::function< Point< double , Dim > ( unsigned int idx ) > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::deltaMass( const Indexer & indexer , SampleFunctor && F , size_t sampleNum ) const
{
	return deltaMass( std::forward< Indexer >( indexer ) , std::forward< SampleFunctor >( F ) , sampleNum , []( Point< double , Dim > ){ return 1.;} );
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ , typename SampleFunctor /* = std::function< Point< double , Dim > ( unsigned int idx , unsigned int ) > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::deltaMass( const Indexer & indexer , SampleFunctor && F , const OrderedSampler< Dim > &orderedSampler ) const
{
	return deltaMass( std::forward< Indexer >( indexer ) , std::forward< SampleFunctor >( F ) , orderedSampler , []( Point< double , Dim > , unsigned int ){ return 1.;} );
}


template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::deltaMass( const Indexer & indexer , const std::vector< Point< double , Dim > > &samples ) const
{
	return deltaMass( std::forward< Indexer >( indexer ) , samples , []( Point< double , Dim > ){ return 1.;} );
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::boundaryMass( const Indexer &indexer ) const
{
	if constexpr( Dim==0 )
	{
		ERROR_OUT( "Cannot compute boundary for dimension zero" );
		return Eigen::SparseMatrix< double >( indexer.numFunctions() , indexer.numFunctions() );
	}
	else
	{
		typename ScalarFunctions< Dim-1 >::template IntegrationStencil< double , 2 , 0 > stencil = ScalarFunctions< Dim-1 >::MassStencil( _r );
		return boundarySystemMatrix( indexer , stencil );
	}
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::boundaryStiffness( const Indexer &indexer ) const
{
	if constexpr( Dim==0 || Dim==1 )
	{
		ERROR_OUT( "Cannot compute boundary stiffness for dimension zero or one" );
		return Eigen::SparseMatrix< double >( indexer.numFunctions() , indexer.numFunctions() );
	}
	else
	{
		typename ScalarFunctions< Dim-1 >::template IntegrationStencil< double , 2 , 0 > stencil = ScalarFunctions< Dim-1 >::StiffnessStencil( _r );
		return boundarySystemMatrix( indexer , stencil );
	}
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ >
double ScalarFunctions< Dim >::value( const Indexer & indexer , const Eigen::VectorXd &x , Point< double , Dim > p , unsigned int thread ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	Window::IsotropicStaticWindow< size_t , Dim , 2 > neighbors = indexer.fNeighbors( p , thread );

	p *= _r;
	Index< Dim > E;
	for( unsigned int d=0 ; d<Dim ; d++ ) E[d] = (int)floor( p[d] );
	for( unsigned int d=0 ; d<Dim ; d++ ) if( E[d]==_r ) E[d]--;
	for( unsigned int d=0 ; d<Dim ; d++ ) p[d] -= E[d];

	double values[2][Dim];
	for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) values[i][d] = i==0 ? (1.-p[d]) : p[d];

	double value = 0;
	auto Kernel = [&]( Index< Dim > dF )
		{
			size_t f = neighbors( &dF[0] );
			if( f!=-1 )
			{
				double v = x[f];
				for( unsigned int d=0 ; d<Dim ; d++ ) v *= dF[d]==0 ? values[0][d] : values[1][d];
				value += v;
			}
		};
	Basis< Dim >::ElementSupport( Index< Dim >() ).process( Kernel );
	return value;
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ >
Point< double , Dim > ScalarFunctions< Dim >::gradient( const Indexer & indexer , const Eigen::VectorXd &x , Point< double , Dim > p , unsigned int thread ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	Window::IsotropicStaticWindow< size_t , Dim , 2 > neighbors = indexer.fNeighbors( p , thread );

	p *= _r;
	Index< Dim > E;
	for( unsigned int d=0 ; d<Dim ; d++ ) E[d] = (int)floor( p[d] );
	for( unsigned int d=0 ; d<Dim ; d++ ) if( E[d]==_r ) E[d]--;
	for( unsigned int d=0 ; d<Dim ; d++ ) p[d] -= E[d];

	double values[2][Dim] , dValues[2][Dim];
	for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ )
	{
		values[i][d] = i==0 ? (1.-p[d]) : p[d];
		dValues[i][d] = i==0 ? -(int)_r : (int)_r;
	}

	Point< double , Dim > gradient;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		auto Kernel = [&]( Index< Dim > dF )
			{
				size_t f = neighbors( &dF[0] );
				if( f!=-1 )
				{
					double v = x[f];
					for( unsigned int dd=0 ; dd<Dim ; dd++ )
						if( d!=dd ) v *= dF[dd]==0 ?  values[0][dd] :  values[1][dd];
						else        v *= dF[dd]==0 ? dValues[0][dd] : dValues[1][dd];
					gradient[d] += v;
				}
			};
		Basis< Dim >::ElementSupport( Index< Dim >() ).process( Kernel );
	}
	return gradient;
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::Indexer< Dim > */ , typename PiecewiseConstantScalarField /* = std::function< double ( Hat::Index< Dim > E ) > */ >
Eigen::VectorXd ScalarFunctions< Dim >::valueDual( const Indexer &indexer , PiecewiseConstantScalarField SF ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );
	static_assert( std::is_convertible_v< PiecewiseConstantScalarField , std::function< double ( Hat::Index< Dim > ) > > , "[ERROR] PiecewiseConstantScalarField poorly formed" );

	IntegrationStencil< double , 1 , 0 >  stencil = ValueStencil( _r );

	Eigen::VectorXd d( indexer.numFunctions() );

	ThreadPool::ParallelFor
		(
			0 , indexer.numFunctions() ,
			[&]( unsigned int t , size_t f )
			{
				Index< Dim > F = indexer.functionIndex( f );
				Window::IsotropicStaticWindow< size_t , Dim , 2 > neighbors = indexer.feNeighbors( f , t );

				double value = 0;
				for( unsigned int i=0 ; i<Window::IsotropicSize< Dim , 2 >() ; i++ )
				{
					size_t e = neighbors.data[i];
					if( e!=-1 )
					{
						Index< Dim > E = indexer.elementIndex( e );
						value += stencil( E , F ) * SF( E );
					}
				}
				d[f] = value;
			}
		);

	return d;
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::Indexer< Dim > */ >
Eigen::VectorXd ScalarFunctions< Dim >::valueDual( const Indexer & indexer , const Eigen::VectorXd &x ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	IntegrationStencil< double , 2 , 0 >  stencil = MassStencil( _r );

	Eigen::VectorXd d( indexer.numFunctions() );
	ThreadPool::ParallelFor
		(
			0 , indexer.numFunctions() ,
			[&]( unsigned int t , size_t f1 )
			{
				Index< Dim > F1 = indexer.functionIndex( f1 );
				Window::IsotropicStaticWindow< size_t , Dim , 2 > eNeighbors = indexer.feNeighbors( f1 , t );

				double value = 0;

				for( unsigned int i=0 ; i<Window::IsotropicSize< Dim , 2 >() ; i++ )
				{
					size_t e = eNeighbors.data[i];
					if( e!=-1 )
					{
						Index< Dim > E = indexer.elementIndex( e );
						Window::IsotropicStaticWindow< size_t , Dim , 2 > fNeighbors = indexer.efNeighbors( e , t );
						for( unsigned int i=0 ; i<Window::IsotropicSize< Dim , 2 >() ; i++ )
						{
							size_t f2 = fNeighbors.data[i];
							if( f2!=-1 )
							{
								Index< Dim > F2 = indexer.functionIndex( f2 );
								value += stencil( E , F1 , F2 ) * x[f2];
							}
						}
					}
				}

				d[f1] = value;
			}
		);

	return d;
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::Indexer< Dim > */ , typename PiecewiseConstantVectorField /* = std::function< Point< double , Dim > ( Hat::Index< Dim > E ) > */ >
Eigen::VectorXd ScalarFunctions< Dim >::gradientDual( const Indexer & indexer , PiecewiseConstantVectorField VF ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );
	static_assert( std::is_convertible_v< PiecewiseConstantVectorField , std::function< Point< double , Dim > ( Hat::Index< Dim > ) > > , "[ERROR] PiecewiseConstantVectorField poorly formed" );

	IntegrationStencil< Point< double , Dim > , 1 , 0 >  stencil = PartialDerivativeStencil( _r );

	Eigen::VectorXd d( indexer.numFunctions() );

	ThreadPool::ParallelFor
		(
			0 , indexer.numFunctions() ,
			[&]( unsigned int t , size_t f )
			{
				Index< Dim > F = functionIndex( f );

				Window::IsotropicStaticWindow< size_t , Dim , 2 > neighbors = indexer.feNeighbors( f , t );

				double value = 0;
				for( unsigned int i=0 ; i<Window::IsotropicSize< Dim , 2 >() ; i++ )
				{
					size_t e = neighbors.data[i];
					if( e!=-1 )
					{
						Hat::Index< Dim > E = indexer.elementIndex( e );
						value += Point< double , Dim >::Dot( stencil( E , F ) , VF( E ) );
					}
				}
				d[f] = value;
			}
		);

	return d;
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::Indexer< Dim > */ >
Eigen::VectorXd ScalarFunctions< Dim >::gradientDual( const Indexer & indexer , const Eigen::VectorXd &x ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	IntegrationStencil< double , 2 , 0 >  stencil = StiffnessStencil( _r );

	Eigen::VectorXd d( indexer.numFunctions() );

	ThreadPool::ParallelFor
		(
			0 , indexer.numFunctions() , 
			[&]( unsigned int t , size_t f1 )
			{
				Index< Dim > F1 = indexer.functionIndex( f1 );
				Window::IsotropicStaticWindow< size_t , Dim , 2 > eNeighbors = indexer.feNeighbors( f1 , t );

				double value = 0;

				for( unsigned int i=0 ; i<Window::IsotropicSize< Dim , 2 >() ; i++ )
				{
					size_t e = eNeighbors.data[i];
					if( e!=-1 )
					{
						Index< Dim > E = indexer.elementIndex( e );
						Window::IsotropicStaticWindow< size_t , Dim , 2 > fNeighbors = indexer.efNeighbors( e , t );
						for( unsigned int i=0 ; i<Window::IsotropicSize< Dim , 2 >() ; i++ )
						{
							size_t f2 = fNeighbors.data[i];
							if( f2!=-1 )
							{
								Index< Dim > F2 = indexer.functionIndex( f2 );
								value += stencil( E , F1 , F2 ) * x[f2];
							}
						}
					}
				}

				d[f1] = value;
			}
		);

	return d;
}

template< unsigned int Dim >
template< typename Data , typename Indexer /* = Hat::Indexer< Dim > */ , typename SampleFunctor /* = std::function< std::pair< Point< double , Dim > , Data > ( unsigned int idx ) > */ , typename WeightFunctor /* = std::function< double ( Point< double , Dim > ) */ >
std::vector< Data > ScalarFunctions< Dim >::deltaDual( const Indexer &indexer , SampleFunctor F , size_t sampleNum , WeightFunctor wF ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

	std::vector< Data > b( indexer.numFunctions() );
	for( unsigned int i=0 ; i<indexer.numFunctions() ; i++ ) b[i] = {};

	// The values of the corner functions at the samples
	struct FunctionValues
	{
		double values[1<<Dim];
		FunctionValues( void ){}

		FunctionValues( Point< double , Dim > p )
		{
			Index< Dim > e;
			double _values[2][Dim];
			for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) _values[i][d] = i==0 ? (1.-p[d]) : p[d];

			auto f = [&]( Index< Dim > f )
				{
					double value = 1.;
					for( unsigned int d=0 ; d<Dim ; d++ ) value *= f[d]==e[d] ? _values[0][d] : _values[1][d];
					values[ Basis< Dim >::template _RelativeIndex<0>( e , f ) ] = value;
				};
			Basis< Dim >::ElementSupport( e ).process( f );
		}
	};

	struct IndexedSample
	{
		Point< double , Dim > position;
		Data data;
		size_t index;
		double weight;
		IndexedSample( Point< double , Dim > p , Data d , double w , unsigned int res ) : position(p) , data(d) , weight(w) , index(0)
		{
			static const unsigned int NumBits = std::min< unsigned int >( ( sizeof(size_t)*8 ) / Dim , 16 );
			static const size_t MaxRes = (size_t)1<<NumBits;
			if( res>MaxRes ) ERROR_OUT( "Resolution too large: " , res , " >= " , MaxRes , " " , NumBits );
			Index< Dim > idx;
			for( unsigned int d=0 ; d<Dim ; d++ )
			{
				if( position[d]<0 || position[d]>=1 ) ERROR_OUT( "Sample out bounds: " , p );
				idx[d] = (unsigned int)floor( position[d] * res );
			}
			for( unsigned int i=0 ; i<NumBits ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ ) index |= ( (idx[d]>>i)&1 ) << ( i*Dim + d );
		}
	};

	std::vector< IndexedSample > samples;
	samples.reserve( sampleNum );
	for( unsigned int i=0 ; i<sampleNum ; i++ )
	{
		std::pair< Point< double , Dim > , Data > sample = F(i);
		double w = wF( sample.first );
		if( w>0 ) samples.emplace_back( sample.first , sample.second , w , _r );
	}
	std::sort( samples.begin() , samples.end() , []( const IndexedSample &i1 , const IndexedSample &i2 ){ return i1.index<i2.index; } );

	std::vector< FunctionValues > functionValues;
	unsigned int start=0;
	while( start<samples.size() )
	{
		// compute the [start,end) range for samples mapping to the same cell
		unsigned int end;
		for( end=start ; end<samples.size() && samples[end].index==samples[start].index ; end++ );

		auto SampleCell = [&]( Point< double , Dim > p )
			{
				Index< Dim > E;
				p *= _r;
				for( unsigned int d=0 ; d<Dim ; d++ ) E[d] = (int)floor( p[d] );
				for( unsigned int d=0 ; d<Dim ; d++ ) if( E[d]==_r ) E[d]--;
				return E;
			};

		// Get the cell using one of the samples (e.g. the first)
		Window::IsotropicStaticWindow< size_t , Dim , 2 > neighbors = indexer.fNeighbors( samples[start].position , 0 );
		Index< Dim > E = SampleCell( samples[start].position );

		// Get the values of each of the function corner functions at each of the samples
		functionValues.resize(0);
		functionValues.reserve( end-start );
		for( unsigned int i=start ; i<end ; i++ )
		{
			Point< double , Dim > p = samples[i].position * _r;
			for( unsigned int d=0 ; d<Dim ; d++ ) p[d] -= E[d];
			functionValues.emplace_back( p );
		}

		for( unsigned int i=0 ; i<Window::IsotropicSize< Dim , 2 >() ; i++ )
		{
			size_t f = neighbors.data[i];
			if( f!=-1 )
			{
				Hat::Index< Dim > F = indexer.functionIndex( f );
				Data data = {};
				for( unsigned int i=start ; i<end ; i++ ) data += samples[i].data * functionValues[i-start].values[ Basis< Dim >::template _RelativeIndex<0>( E , F ) ] * samples[i].weight;
				b[f] += data;
			}
		}

		start = end;
	}

	return b;
}

template< unsigned int Dim >
template< typename Data , typename Indexer /* = Hat::Indexer< Dim > */ , typename WeightFunctor /* = std::function< double ( Point< double , Dim > ) */ >
std::vector< Data > ScalarFunctions< Dim >::deltaDual( const Indexer &indexer , const std::vector< std::pair< Point< double , Dim > , Data > > &samples , WeightFunctor wF ) const
{
	return deltaDual< Data >( indexer , [&]( unsigned int idx ){ return samples[idx]; } , samples.size() , wF );
}

template< unsigned int Dim >
template< typename Data , typename Indexer /* = Hat::Indexer< Dim > */ , typename SampleFunctor /* = std::function< std::pair< Point< double , Dim > , Data > ( unsigned int idx ) > */ >
std::vector< Data > ScalarFunctions< Dim >::deltaDual( const Indexer &indexer , SampleFunctor F , size_t sampleNum ) const
{
	return deltaDual< Data >( indexer , F , sampleNum , []( Point< double , Dim > ){ return 1.;} );
}

template< unsigned int Dim >
template< typename Data , typename Indexer /* = Hat::Indexer< Dim > */ >
std::vector< Data > ScalarFunctions< Dim >::deltaDual( const Indexer &indexer , const std::vector< std::pair< Point< double , Dim > , Data > > &samples ) const
{
	return deltaDual< Data >( indexer , samples , []( Point< double , Dim > ){ return 1.;} );
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::BaseIndexer< Dim > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::mass( const Indexer & indexer ) const { return systemMatrix( indexer , MassStencil(_r) ); }

template< unsigned int Dim >
template< typename Indexer /* = Hat::Indexer< Dim > */ , typename MassFunctor /* = std::function< double ( Index< Dim > E ) > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::mass( const Indexer & indexer , MassFunctor && mf ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );
	static_assert( std::is_convertible_v< MassFunctor , std::function< double ( Index< Dim > ) > > , "[ERROR] Poorly formed MassFunctor" );

	Hat::Index< Dim > Off;
	for( unsigned int d=0 ; d<Dim ; d++ ) Off[d] = 1;
	Range< Dim > eRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) eRange.first[d] = 0 , eRange.second[d] = _r;
	IntegrationStencil< double , 2 , 0 > stencil = MassStencil( _r );

	// Determine which functions are supported on some cell
	std::vector< bool > supportedFunctions( indexer.numFunctions() , false );

	// [WARNING] Don't parallelize because of write-on-write conflict
	for( size_t e=0 ; e<indexer.numElements() ; e++ )
	{
		Hat::Index< Dim > E = indexer.elementIndex(e);
		if( mf( E )!=0 )
		{
			Window::IsotropicStaticWindow< size_t , Dim , 2 > neighbors = indexer.efNeighbors( e , 0 );
			for( unsigned int i=0 ; i<Window::IsotropicSize< Dim , 2 >() ; i++ )
			{
				size_t f = neighbors.data[i];
				if( f!=-1 ) supportedFunctions[f] = true;
			}
		}
	}

	Eigen::SparseMatrix< double > S( indexer.numFunctions() , indexer.numFunctions() );

	MatrixInfo< 1 , true > matrixInfo( _r );

	Eigen::VectorXi rowSizes( (int)indexer.numFunctions() );
	ThreadPool::ParallelFor
	(
		0 , indexer.numFunctions() ,
		[&]( unsigned int , size_t f )
		{
			Index< Dim > F = indexer.functionIndex(f);
			rowSizes[f] = supportedFunctions[f] ? (int)matrixInfo.entries(F,true) : 0;
		}
	);
	S.reserve( rowSizes );

	ThreadPool::ParallelFor
	(
		0 , indexer.numFunctions() ,
		[&]( unsigned int t , size_t f1 )
		{
			Index< Dim > F1 = indexer.functionIndex( f1 );

			if( supportedFunctions[f1] )
			{
				Window::IsotropicStaticWindow< size_t , Dim , 3 > neighbors = indexer.ffNeighbors( f1 , 0 );

				auto Kernel = [&]( Index< Dim > F2 , size_t e , bool flip )
					{
						Hat::Index< Dim > I = F2-F1+Off;
						size_t f2 = neighbors( &I[0] );
						double dot = 0;
						if( f2!=-1 )
						{
							auto Kernel = [&]( Index< Dim > E ){ dot += stencil( E , F1 , F2 ) * mf( E ); };
							// Iterate over all cells supported by the i-th function
							Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( F1 ) , Basis< Dim >::FunctionSupport( F2 ) , eRange ).process( Kernel );
						}

						if( dot ) S.insert( (int)f2 , (int)f1 ) = dot;
					};
				matrixInfo.processAll( F1 , Kernel );
			}
		}
	);

	S.makeCompressed();
	return S;
}

template< unsigned int Dim >
template< typename Indexer /* = Hat::Indexer< Dim > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::stiffness( const Indexer & indexer ) const { return systemMatrix( indexer , StiffnessStencil(_r) ); }

template< unsigned int Dim >
template< typename Indexer /* = Hat::Indexer< Dim > */ , typename InnerProductFunctor /* = std::function< MishaK::SquareMatrix< double , Dim > ( Index< Dim > E [ , unsigned int t [ ) > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::stiffness( const Indexer & indexer , InnerProductFunctor && ipf ) const
{
	static_assert( std::is_base_of_v< Hat::BaseIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );
	static_assert( std::is_convertible_v< InnerProductFunctor , std::function< MishaK::SquareMatrix< double , Dim > ( Index< Dim > ) > > || std::is_convertible_v< InnerProductFunctor , std::function< MishaK::SquareMatrix< double , Dim > ( Index< Dim > , unsigned int ) > > , "[ERROR] InnerProductFunctor poorly formed" );

	static const bool NeedsThreadID = std::is_convertible_v< InnerProductFunctor , std::function< MishaK::SquareMatrix< double , Dim > ( Index< Dim > , unsigned int ) > >;
	IntegrationStencil< MishaK::SquareMatrix< double, Dim >, 2, 0 > stencil = FullStiffnessStencil(_r);
	Hat::Index< Dim > Off;
	Range< Dim > eRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) eRange.first[d] = 0, eRange.second[d] = _r , Off[d] = 1;


	// Determine which functions are supported on a cell where the inner-product functor is non-zero
	std::vector< char > supportedFunctions( indexer.numFunctions() , 0 );
	ThreadPool::ParallelFor
	(
		0 , indexer.numElements() ,
		[&]( unsigned int t , size_t e )
		{
			Hat::Index< Dim > E = indexer.elementIndex(e);
			bool supported;
			if constexpr( NeedsThreadID ) supported = ipf( E ,t ).squareNorm() != 0;
			else                          supported = ipf( E ).squareNorm() != 0;
			if( supported )
			{
				Window::IsotropicStaticWindow< size_t , Dim , 2 > neighbors = indexer.efNeighbors( e , t );
				for( unsigned int i=0 ; i<Window::IsotropicSize< Dim , 2 >() ; i++ )
				{
					size_t f = neighbors.data[i];
					if( f!=-1 ) SetAtomic( supportedFunctions[f] , (char)1 );
				}
			}
		}
	);

	Eigen::SparseMatrix< double > S( indexer.numFunctions() , indexer.numFunctions() );

	MatrixInfo< 1 , true > matrixInfo(_r);

	Eigen::VectorXi rowSizes((int)indexer.numFunctions());
	ThreadPool::ParallelFor
	(
		0 , indexer.numFunctions() ,
		[&]( unsigned int , size_t f )
		{
			Index< Dim > F = indexer.functionIndex(f);
			rowSizes[f] = supportedFunctions[f] ? (int)matrixInfo.entries( F , true ) : 0;
		}
	);
	S.reserve( rowSizes );

	ThreadPool::ParallelFor
	(
		0 , indexer.numFunctions() ,
		[&]( unsigned int t , size_t f1 )
		{
			Index< Dim > F1 = indexer.functionIndex(f1);

			if( supportedFunctions[f1] )
			{
				Window::IsotropicStaticWindow< size_t , Dim , 3 > fNeighbors = indexer.ffNeighbors( f1 , t );
				Window::IsotropicStaticWindow< size_t , Dim , 2 > eNeighbors = indexer.feNeighbors( f1 , t );

				// Precompute the values on the elements supported by F1
				Window::IsotropicStaticWindow< MishaK::SquareMatrix< double , Dim > , Dim , 2 >  ips;
				{
					auto Kernel = [&]( Hat::Index< Dim > E )
						{
							Hat::Index< Dim > I = E - F1 + Off;
							if( eNeighbors( &I[0] )!=-1 )
								if constexpr( NeedsThreadID ) ips( &I[0] ) = ipf( E , t );
								else                          ips( &I[0] ) = ipf( E );
						};
					Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( F1 ) , eRange ).process( Kernel );
				}

				auto Kernel = [&]( Index< Dim > F2 , size_t e , bool flip )
					{
						Hat::Index< Dim > I = F2-F1+Off;
						size_t f2 = fNeighbors( &I[0] );
						if( f2!=-1 && supportedFunctions[f2] )
						{
							double dot = 0;
							auto f = [&]( Index< Dim > E )
								{
									Hat::Index< Dim > I = E - F1 + Off;
									MishaK::SquareMatrix< double, Dim > m = ips( &I[0] );
									unsigned int i1 = Basis< Dim >::template _RelativeIndex< 0 >( E , F1 );
									unsigned int i2 = Basis< Dim >::template _RelativeIndex< 0 >( E , F2 );
									dot += MishaK::SquareMatrix< double , Dim >::Dot( m , stencil[i1][i2] );
								};
							// Iterate over all cells supported by the i-th function
							Range< Dim >::Intersect( Basis< Dim >::FunctionSupport(F1) , Basis< Dim >::FunctionSupport(F2) , eRange ).process(f);

							if( dot ) S.insert( (int)f2 , (int)f1 ) = dot;
						}
					};
				matrixInfo.processAll( F1 , Kernel );
			}
		}
	);

	S.makeCompressed();
	return S;
}



template< unsigned int Dim >
template< typename ProlongationIndexer /* = Hat::BaseProlongationIndexer< Dim > */ >
Eigen::SparseMatrix< double > ScalarFunctions< Dim >::prolongation( const ProlongationIndexer &pIndexer , size_t numFineFunctions ) const
{
	static_assert( std::is_base_of_v< Hat::BaseProlongationIndexer< Dim > , ProlongationIndexer > , "[ERROR] ProlongationIndexer poorly formed" );

	Eigen::SparseMatrix< double > P( numFineFunctions , pIndexer.numFunctions() );
	Eigen::VectorXi rowSizes( (int)pIndexer.numFunctions() );
	ProlongationStencil pStencil;

	ThreadPool::ParallelFor
		(
			0 , pIndexer.numFunctions() ,
			[&]( unsigned int t , size_t f )
			{
				Index< Dim > F = pIndexer.functionIndex(f);
				Window::IsotropicStaticWindow< size_t , Dim , 3 > neighbors = pIndexer.ffChildNeighbors( f , t );
				unsigned int sz = 0;
				for( unsigned int i=0 ; i<Window::IsotropicSize< Dim , 3 >() ; i++ ) if( neighbors.data[i]!=-1 ) sz++;
				rowSizes[f] = sz;
			}
		);

	P.reserve( rowSizes );

	Range< Dim > range;
	Index< Dim > Off;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = 2*resolution()+1 , Off[d] = 1;

	ThreadPool::ParallelFor
		(
			0 , pIndexer.numFunctions() ,
			[&]( unsigned int t , size_t f1 )
			{
				Window::IsotropicStaticWindow< size_t , Dim , 3 > neighbors = pIndexer.ffChildNeighbors( f1 , t );
				Index< Dim > F1 = pIndexer.functionIndex(f1)*2;

				auto Kernel = [&]( Index< Dim > F2 )
					{
						Hat::Index< Dim > I = F2 - F1 + Off;
						size_t f2 = neighbors( &I[0] );
						if( f2!=-1 ) P.insert( (int)f2 , (int)f1 ) = pStencil( I );
					};
				Range< Dim >::Intersect( range , Range< Dim >(F1).dilate(1) ).process( Kernel );
			}
		);

	P.makeCompressed();

	return P;
};

template< unsigned int Dim >
template< typename ProlongationIndexer /* = Hat::BaseProlongationIndexer< Dim > */ >
Eigen::VectorXd ScalarFunctions< Dim >::prolongation( const Eigen::VectorXd &x , const ProlongationIndexer &pIndexer , size_t numFineFunctions ) const
{
	static_assert( std::is_base_of_v< Hat::BaseProlongationIndexer< Dim > , ProlongationIndexer > , "[ERROR] ProlongationIndexer poorly formed" );

	Eigen::VectorXd xP = Eigen::VectorXd::Zero( numFineFunctions );

	Range< Dim > range;
	Index< Dim > Off;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = 2*resolution()+1 , Off[d] = 1;

	ProlongationStencil pStencil;

	ThreadPool::ParallelFor
	(
		0 , pIndexer.numFunctions() ,
		[&]( unsigned int t , size_t f1 )
		{
			Window::IsotropicStaticWindow< size_t , Dim , 3 > neighbors = pIndexer.ffChildNeighbors( f1 , t );
			Index< Dim > F1 = pIndexer.functionIndex(f1)*2;

			auto Kernel = [&]( Index< Dim > F2 )
				{
					Hat::Index< Dim > I = F2 - F1 + Off;
					size_t f2 = neighbors( &I[0] );
					if( f2!=-1 ) AddAtomic( xP[f2] , x[f1] * pStencil(I) );
				};
			Range< Dim >::Intersect( range , Range< Dim >(F1).dilate(1) ).process( Kernel );
		}
	);

	return xP;
}


//////////////////////////////////////
// ScalarFunctions< Dim >::UnitTest //
//////////////////////////////////////
template< unsigned int Dim >
Eigen::VectorXd ScalarFunctions< Dim >::UnitTest::RandomCoefficients( const ScalarFunctions &scalars )
{
	Eigen::VectorXd v( scalars.functionNum() );
	for( unsigned int i=0 ; i<scalars.functionNum() ; i++ ) v[i] = Random< double >()*2 - 1;
	return v;
}


template< unsigned int Dim >
void ScalarFunctions< Dim >::UnitTest::operator()( void ) const
{
	linearEvaluation();
	prolongation();
	stencils();
	boundaryStencils();
	integrationStencils();
	massStiffnessAndDual();
	dot();
}

template< unsigned int Dim >
void ScalarFunctions< Dim >::UnitTest::linearEvaluation( void ) const
{
	Eigen::VectorXd x( _scalars.functionNum() );
	Range< Dim > fRange = _scalars.functionRange();

	double vErr = 0 , gErr = 0;
	for( unsigned int n=0 ; n<_N ; n++ )
	{
		Point< double , Dim > q = RandomSpherePoint< double , Dim >();
		double off = Random< double >();
		auto L = [&]( Point< double , Dim > p ) { return off + Point< double , Dim >::Dot( p , q ); };
		
		auto f = [&]( Index< Dim > F )
		{
			Point< double , Dim > p;
			for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = (double)F[d]/_scalars.resolution();
			x[ _scalars.functionIndex(F) ] = L(p);
		};
		fRange.process( f );

		for( unsigned i=0 ; i<_numSamples ; i++ )
		{
			Point< double , Dim > p;
			for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = Random< double >();
			{
				double v1 = L(p);
				double v2 = _scalars.value( x , p );
				vErr += (v1-v2)*(v1-v2);
			}
			{
				Point< double , Dim > v1 = q;
				Point< double , Dim > v2 = _scalars.gradient( x , p );
				gErr += (v1-v2).squareNorm();
			}
		}
	}
	std::cout << "\tLinear evaluation: " << vErr << " / " << gErr << std::endl;
}

template< unsigned int Dim >
void ScalarFunctions< Dim >::UnitTest::prolongation( void ) const
{
	ScalarFunctions< Dim > pScalars( _scalars.resolution()*2 );
	Range< Dim > fRange = _scalars.functionRange();
	Eigen::SparseMatrix< double > prolongation = pScalars.prolongation();

	double err = 0;
	for( unsigned int n=0 ; n<_N ; n++ )
	{
		Eigen::VectorXd x = RandomCoefficients( _scalars );
		Eigen::VectorXd Px = prolongation * x;

		for( unsigned i=0 ; i<_numSamples ; i++ )
		{
			Point< double , Dim > p;
			for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = Random< double >();
			double v1 = _scalars.value( x , p ) , v2 = pScalars.value( Px , p );
			err += (v1-v2)*(v1-v2);
		}
	}
	std::cout << "\tProlongation: " << err << std::endl;
}

template< unsigned int Dim >
void ScalarFunctions< Dim >::UnitTest::stencils( void ) const
{
	double symErr = 0 , aSymErr = 0;
	IntegrationStencil< double , 2 , 0 > stencil;

	for( unsigned int n=0 ; n<_N ; n++ )
	{
		auto StencilError = [&]( void )
		{
			FullIntegrationStencil< double , 0 > fullStencil( stencil , _scalars.resolution() );
			Eigen::VectorXd x = RandomCoefficients( _scalars );
			Eigen::VectorXd y = RandomCoefficients( _scalars );
			Eigen::SparseMatrix< double > M = _scalars.systemMatrix( stencil );
			Eigen::SparseMatrix< double > fullM = _scalars.systemMatrix( fullStencil );

			double values[] =
			{
				_scalars( stencil , x , y ) ,
				_scalars( fullStencil , x , y ) ,
				( M * x ).dot(y) ,
				( fullM * x ).dot(y) ,
				_scalars( stencil , x ).dot(y) ,
				_scalars( fullStencil , x ).dot(y)
			};
			return Miscellany::AverageAndStandardDeviation( values , sizeof(values)/sizeof(double) ).second;
		};
		{
			for( unsigned int i=0 ; i<StencilSize<0>() ; i++ ) for( unsigned int j=0 ; j<=i ; j++ ) stencil[i][j] = stencil[j][i] = Random< double >()*2-1;
			symErr += StencilError();
		}
		{
			for( unsigned int i=0 ; i<StencilSize<0>() ; i++ )
			{
				stencil[i][i] = 0;
				for( unsigned int j=0 ; j<=i ; j++ )
				{
					stencil[j][i] = Random< double >()*2-1;
					stencil[i][j] = -stencil[j][i];
				}
			}
			aSymErr += StencilError();
		}
	}

	std::cout << "\tStencils: " << symErr << " / " << aSymErr << std::endl;
}

template< unsigned int Dim >
void ScalarFunctions< Dim >::UnitTest::boundaryStencils( void ) const
{
	if constexpr( Dim>1 )
	{
		ScalarFunctions< Dim-1 > __scalars( _scalars.resolution() );

		auto GetIndex = [&]( Index< Dim-1 > _F , unsigned int d , bool start )
		{
			Index< Dim > F;
			for( unsigned int _d=0 ; _d<d ; _d++ ) F[_d] = _F[_d];
			F[d] = start ? 0 : _scalars.resolution();
			for( unsigned int _d=d+1 ; _d<Dim ; _d++ ) F[_d] = _F[_d-1];
			return F;
		};
		auto Slice = [&]( const Eigen::VectorXd &x , unsigned int d , bool start )
		{
			Eigen::VectorXd _x( __scalars.functionNum() );
			for( unsigned int i=0 ; i<__scalars.functionNum() ; i++ )
			{
				Index< Dim-1 > _F = __scalars.functionIndex(i);
				Index< Dim > F = GetIndex( _F , d , start );
				_x[i] = x[ _scalars.functionIndex(F) ];
			}
			return _x;
		};

		double mErr = 0 , sErr = 0;

		typename ScalarFunctions< Dim-1 >::template IntegrationStencil< double , 2 , 0 > mStencil = ScalarFunctions< Dim-1 >::     MassStencil( _scalars.resolution() );
		typename ScalarFunctions< Dim-1 >::template IntegrationStencil< double , 2 , 0 > sStencil = ScalarFunctions< Dim-1 >::StiffnessStencil( _scalars.resolution() );

		for( unsigned int n=0 ; n<_N ; n++ )
		{
			Eigen::VectorXd x = RandomCoefficients( _scalars );
			Eigen::VectorXd y = RandomCoefficients( _scalars );

			auto StencilError = [&]( typename ScalarFunctions< Dim-1 >::template IntegrationStencil< double , 2 , 0 > &stencil )
			{
				Eigen::SparseMatrix< double > M = _scalars.boundarySystemMatrix( stencil );
				double v1 = ( M * x ).dot(y);
				double v2 = 0;
				for( unsigned int d=0 ; d<Dim ; d++ ) for( unsigned int s=0 ; s<2 ; s++ )
				{
					Eigen::VectorXd _x = Slice( x , d , s==0 );
					Eigen::VectorXd _y = Slice( y , d , s==0 );
					v2 += __scalars( stencil , _x , _y );
				}
				return (v1-v2)*(v1-v2);
			};
			mErr += StencilError( mStencil );
			sErr += StencilError( sStencil );
		}
		std::cout << "\tBoundary stencils: " << mErr << " / " << sErr << std::endl;
	}
}

template< unsigned int Dim >
void ScalarFunctions< Dim >::UnitTest::integrationStencils( void ) const
{
	Range< Dim > sampleRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) sampleRange.second[d] = _numSamples;

	std::cout << "\tDiscrete integration:" << std::endl;
	Index< Dim > E;
	// Test single value integration
	{
		IntegrationStencil< double , 1 , 0 > stencil = ValueStencil( _scalars.resolution() );

		double err = 0;

		auto f = [&]( Index< Dim > F )
		{
			double v1 = stencil(E,F) , v2 = 0;

			auto f = [&]( Index< Dim > I )
			{
				Point< double , Dim > p;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = ( I[d] + 0.5 ) / _numSamples;
				v2 += Hat::ElementFunction< Dim >( E , F )( p );
			};
			sampleRange.process( f );
			for( unsigned int d=0 ; d<Dim ; d++ ) v2 /= _numSamples * _scalars.resolution();

			err += ( v1-v2 ) * ( v1-v2 );
		};
		Basis< Dim >::ElementSupport( E ).process( f );

		std::cout << "\t\t            Value: " << err << std::endl;
	}

	// Test single gradient integration
	{
		IntegrationStencil< Point< double , Dim > , 1 , 0 > stencil = PartialDerivativeStencil( _scalars.resolution() );

		double err = 0;

		auto f = [&]( Index< Dim > F )
		{
			Point< double, Dim > v1 = stencil(E,F) , v2;

			auto f = [&]( Index< Dim > I )
			{
				Point< double , Dim > p;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = ( I[d] + 0.5 ) / _numSamples;
				v2 += Hat::ElementFunction< Dim >( E , F ).gradient()( p ) * _scalars.resolution();
			};
			sampleRange.process( f );
			for( unsigned int d=0 ; d<Dim ; d++ ) v2 /= _numSamples * _scalars.resolution();

			err += ( v1-v2 ).squareNorm();
		};
		Basis< Dim >::ElementSupport( E ).process( f );

		std::cout << "\t\t          Partial: " << err << std::endl;
	}

	// Test pairwise value integration
	{
		IntegrationStencil< double , 2 , 0 > stencil = MassStencil( _scalars.resolution() );

		double err = 0;

		auto f = [&]( Index< Dim > F1 , Index< Dim > F2 )
		{
			double v1 = stencil(E,F1,F2) , v2 = 0;

			auto f = [&]( Index< Dim > I )
			{
				Point< double , Dim > p;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = ( I[d] + 0.5 ) / _numSamples;
				v2 += Hat::ElementFunction< Dim >( E , F1 )( p ) * Hat::ElementFunction< Dim >( E , F2 )( p );
			};
			sampleRange.process( f );
			for( unsigned int d=0 ; d<Dim ; d++ ) v2 /= _numSamples * _scalars.resolution();

			err += ( v1-v2 ) * ( v1-v2 );
		};
		Basis< Dim >::ElementSupport( E ).template process< 2 >( f );

		std::cout << "\t\t      Value-Value: " << err << std::endl;
	}

	// Test pairwise gradient integration
	{
		IntegrationStencil< double , 2 , 0 > stencil = StiffnessStencil( _scalars.resolution() );

		double err = 0;

		auto f = [&]( Index< Dim > F1 , Index< Dim > F2 )
		{
			double v1 = stencil(E,F1,F2) , v2 = 0;

			auto f = [&]( Index< Dim > I )
			{
				Point< double , Dim > p;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = ( I[d] + 0.5 ) / _numSamples;
				v2 += Point< double , Dim >::Dot( Hat::ElementFunction< Dim >( E , F1 ).gradient()( p ) , Hat::ElementFunction< Dim >( E , F2 ).gradient()( p ) ) * _scalars.resolution() * _scalars.resolution();
			};
			sampleRange.process( f );
			for( unsigned int d=0 ; d<Dim ; d++ ) v2 /= _numSamples * _scalars.resolution();

			err += ( v1-v2 ) * ( v1-v2 );
		};
		Basis< Dim >::ElementSupport( E ).template process< 2 >( f );

		std::cout << "\t\tGradient-Gradient: " << err << std::endl;
	}

	// Test full pairwise gradient integration
	{
		IntegrationStencil< MishaK::SquareMatrix< double , Dim > , 2 , 0 > stencil = FullStiffnessStencil( _scalars.resolution() );

		double err = 0;

		auto f = [&]( Index< Dim > F1 , Index< Dim > F2 )
		{
			MishaK::SquareMatrix< double , Dim > v1 = stencil(E,F1,F2) , v2;

			auto f = [&]( Index< Dim > I )
			{
				Point< double , Dim > p;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = ( I[d] + 0.5 ) / _numSamples;
				v2 += OuterProduct( Hat::ElementFunction< Dim >( E , F1 ).gradient()( p ) , Hat::ElementFunction< Dim >( E , F2 ).gradient()( p ) ) * _scalars.resolution() * _scalars.resolution();
			};
			sampleRange.process( f );
			for( unsigned int d=0 ; d<Dim ; d++ ) v2 /= _numSamples * _scalars.resolution();

			err += ( v1-v2 ).squareNorm();
		};
		Basis< Dim >::ElementSupport( E ).template process< 2 >( f );

		std::cout << "\t\t  Partial-Partial: " << err << std::endl;
	}

	// Test consistency between Stiffness and FullStiffness stencils
	{
		IntegrationStencil< double , 2 , 0 > stencil = StiffnessStencil( _scalars.resolution() );
		IntegrationStencil< MishaK::SquareMatrix< double , Dim > , 2 , 0 > fullStencil = FullStiffnessStencil( _scalars.resolution() );

		double err = 0;

		auto f = [&]( Index< Dim > F1 , Index< Dim > F2 )
		{
			double v1 = fullStencil(E,F1,F2).trace() , v2 = stencil(E,F1,F2);
			err += ( v1-v2 ) * ( v1-v2 );
		};
		Basis< Dim >::ElementSupport( E ).template process< 2 >( f );

		std::cout << "\tFull Stiffness -> Stiffness: " << err << std::endl;
	}
}

template< unsigned int Dim >
void ScalarFunctions< Dim >::UnitTest::massStiffnessAndDual( void ) const
{
	// Test mass discrete vs. analytic
	{
		double err = 0;
		std::vector< double > weights( _scalars.elementNum() );
		for( unsigned int n=0 ; n<_N ; n++ )
		{
			for( unsigned int i=0 ; i<_scalars.elementNum() ; i++ ) weights[i] = Random< double >();
			Eigen::SparseMatrix< double > mass = _scalars.mass( [&]( size_t i ){ return weights[i]; } );
			Eigen::VectorXd x = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars.resolution() );
			Eigen::VectorXd y = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars.resolution() );
			auto w = [&]( Point< double , Dim > p )
				{
					p *= _scalars.resolution();
					Index< Dim > E;
					for( unsigned int d=0 ; d<Dim ; d++ ) E[d] = (int)p[d];
					return weights[ _scalars.elementIndex( E ) ];
				};
			auto f = [&]( Point< double , Dim > p ){ return _scalars.value(x,p) * _scalars.value(y,p) * w(p); };
			double v1 = ( mass * x ).dot(y);
			double v2 = Basis< Dim >::template Integral< double >( f , _numSamples );
			err += (v1-v2) * (v1-v2);
		}
		std::cout << "\tMass matrix (disc): " << err << std::endl;
	}
	// Test the mass
	{
		Eigen::SparseMatrix< double > mass1 = _scalars.mass();
		Eigen::SparseMatrix< double > mass2 = _scalars.mass( []( size_t ){ return 1.; } );
		std::cout << "\t       Mass matrix: " << (mass1-mass2).squaredNorm() << std::endl;
	}

	// Test stiffness discrete vs. analytic
	{
		double err1 = 0 , err2 = 0;
		std::vector< MishaK::SquareMatrix< double , Dim > > weights( _scalars.elementNum() );
		for( unsigned int n=0 ; n<_N ; n++ )
		{
			for( unsigned int i=0 ; i<_scalars.elementNum() ; i++ )
				for( unsigned int d1=0 ; d1<Dim ; d1++ ) for( unsigned int d2=0 ; d2<Dim ; d2++ ) weights[i](d1,d2) = Random< double >();

			Eigen::SparseMatrix< double > stiffness = _scalars.stiffness( [&]( size_t i ){ return weights[i]; } );
			Eigen::VectorXd x = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars.resolution() );
			Eigen::VectorXd y = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars.resolution() );
			auto w = [&]( Point< double , Dim > p )
			{
				p *= _scalars.resolution();
				Index< Dim > E;
				for( unsigned int d=0 ; d<Dim ; d++ ) E[d] = (int)p[d];
				return weights[ _scalars.elementIndex( E ) ];
			};
			auto f1 = [&]( Point< double , Dim > p ){ return Point< double , Dim >::Dot( _scalars.gradient(x,p) , w(p) * _scalars.gradient(y,p) ); };
			auto f2 = [&]( Point< double , Dim > p ){ return MishaK::SquareMatrix< double , Dim >::Dot( w(p) , OuterProduct( _scalars.gradient(x,p) , _scalars.gradient(y,p) ) ); };
			double v1 = ( stiffness * x ).dot(y);
			double v2 = Basis< Dim >::template Integral< double >( f1 , _numSamples );
			double v3 = Basis< Dim >::template Integral< double >( f2 , _numSamples );
			err1 += (v1-v2) * (v1-v2);
			err2 += (v1-v3) * (v1-v3);
		}
		std::cout << "\tStiffness matrix (disc): " << err1 << " / " << err2 << std::endl;
	}

	// Test the stiffness
	{
		Eigen::SparseMatrix< double > stiffness1 = _scalars.stiffness();
		Eigen::SparseMatrix< double > stiffness2 = _scalars.stiffness( []( size_t ){ return MishaK::SquareMatrix< double , Dim >::Identity(); } );
		std::cout << "\t       Stiffness matrix: " << (stiffness1-stiffness2).squaredNorm() << std::endl;
	}

	// Test the value dual
	{
		double e = 0;
		std::vector< double > f1( _scalars.elementNum() );
		Eigen::VectorXd f2( _scalars.functionNum() );
		for( unsigned int n=0 ; n<_N ; n++ )
		{
			Eigen::VectorXd f2 = RandomCoefficients( _scalars );
			for( unsigned int i=0 ; i<_scalars.elementNum() ; i++ ) f1[i] = Random< double >()*2 - 1.;
			Eigen::VectorXd dual = _scalars.valueDual( [&]( size_t e ){ return f1[e]; } );
			auto F1 = [&]( Point< double , Dim > p )
			{
				p *= _scalars.resolution();
				Index< Dim > E;
				for( unsigned int d=0 ; d<Dim ; d++ ) E[d] = (int)floor( p[d] );
				return f1[ _scalars.elementIndex(E) ];
			};
			double dot1 = dual.dot( f2 );
			double dot2 = ScalarDotProduct( F1 , [&]( Point< double , Dim > p ){ return _scalars.value( f2 , p ); } , _numSamples );
			e += (dot1-dot2) * (dot1-dot2);
		}
		std::cout << "\t  Scalar dual: " << e << std::endl;
	}
	// Test the gradient dual
	{
		std::vector< Point< double , Dim > > v1( _scalars.elementNum() );
		double e = 0;
		for( unsigned int n=0 ; n<_N ; n++ )
		{
			for( unsigned int i=0 ; i<_scalars.elementNum() ; i++ ) v1[i] = RandomBallPoint< double , Dim >();
			Eigen::VectorXd f2 = RandomCoefficients( _scalars );
			Eigen::VectorXd dual = _scalars.gradientDual( [&]( size_t e ){ return v1[e]; } );
			auto V1 = [&]( Point< double , Dim > p )
			{
				p *= _scalars.resolution();
				Index< Dim > E;
				for( unsigned int d=0 ; d<Dim ; d++ ) E[d] = (int)floor( p[d] );
				return v1[ _scalars.elementIndex(E) ];
			};
			double dot1 = dual.dot( f2 );
			double dot2 = VectorDotProduct( V1 , [&]( Point< double , Dim > p ){ return _scalars.gradient( f2 , p ); } , _numSamples );
			e += (dot1-dot2) * (dot1-dot2);
		}
		std::cout << "\tGradient dual: " << e << std::endl;
	}
}

template< unsigned int Dim >
void ScalarFunctions< Dim >::UnitTest::dot( void ) const
{
	double massDeviation = 0 , stiffnessDeviation = 0;
	for( unsigned int n=0 ; n<_N ; n++ )
	{
		Eigen::VectorXd x = RandomCoefficients( _scalars );
		Eigen::VectorXd y = RandomCoefficients( _scalars );

		// Test the mass
		{
			IntegrationStencil< double , 2 , 0 > stencil = MassStencil( _scalars.resolution() );
			FullIntegrationStencil< double , 0 > fullStencil( stencil , _scalars.resolution() );
			Eigen::SparseMatrix< double > mass = _scalars.mass();
			Eigen::VectorXd _x = _scalars.valueDual( x );

			double values[] = { _scalars( stencil , x , y ) , _scalars( fullStencil , x , y ) , ( mass * x ).dot(y) , _scalars( stencil , x ).dot(y) , _scalars( fullStencil , x ).dot(y) };
			massDeviation += Miscellany::AverageAndStandardDeviation( values , sizeof(values)/sizeof(double) ).second;
		}
		// Test the stiffness
		{
			IntegrationStencil< double , 2 , 0 > stencil = StiffnessStencil( _scalars.resolution() );
			FullIntegrationStencil< double , 0 > fullStencil( stencil , _scalars.resolution() );
			Eigen::SparseMatrix< double > stiffness = _scalars.stiffness();
			Eigen::VectorXd _x = _scalars.gradientDual( x );

			double values[] = { _scalars( stencil , x , y ) , _scalars( fullStencil , x , y ) , ( stiffness * x ).dot(y) , _scalars( stencil , x ).dot(y) , _scalars( fullStencil , x ).dot(y) };
			stiffnessDeviation += Miscellany::AverageAndStandardDeviation( values , sizeof(values)/sizeof(double) ).second;
		}
	}
	massDeviation /= _N;
	stiffnessDeviation /= _N;
	std::cout << "\tMass/Stiffness dot deviation: " << massDeviation << " / " << stiffnessDeviation << std::endl;
}
