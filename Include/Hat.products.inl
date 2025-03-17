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

//////////////////////////////////////////////
// ProductFunctions::FullIntegrationStencil //
//////////////////////////////////////////////
template< unsigned int Dim , bool Sym >
template< typename T >
ProductFunctions< Dim , Sym >::FullIntegrationStencil< T >::FullIntegrationStencil( const IntegrationStencil< T , 2, 0 > &stencil , unsigned int res ) : _res( res ) , _rows( ScalarFunctions< Dim >::template FullIntegrationStencil< T , 0 >::StencilNum() )
{
	ScalarFunctions< Dim > scalars( _res );
	Range< Dim > range , eRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.first[d] = 0 , range.second[d] = 3 , eRange.first[d] = 0 , eRange.second[d] = _res;


	// Need to define ordering for std::map
	struct ProductIndex
	{
		Index< Dim > g1 , g2;

		ProductIndex( Index< Dim > g1 , Index< Dim > g2 ) : g1(g1) , g2(g2) {}
		bool operator < ( ProductIndex i ) const
		{
			if     ( g1<i.g1 ) return true;
			else if( g1>i.g1 ) return false;
			else return g2<i.g2;
		}
	};

	auto f = [&]( Index< Dim > f1 ) // Process the different boundary types
	{
		// Transform the stencil index to a function index
		for( unsigned int d=0 ; d<Dim ; d++ )
			if     ( f1[d]==1 ) f1[d] = _res/2;
			else if( f1[d]==2 ) f1[d] = _res;

		// Get the indices of all functions whose support overlaps the support of f1
		std::vector< Index< Dim > > f2s;
		auto f = [&]( Index< Dim > f2 ){ if( f1!=f2 || Sym ) f2s.push_back( f2 ); };
		Basis< Dim >::ElementSupport( Range< Dim >::Intersect( eRange , Basis< Dim >::FunctionSupport( f1 ) ) ).process(f);

		std::vector< Row > _row;
		_row.reserve( f2s.size() );
		for( unsigned int i=0 ; i<f2s.size() ; i++ )
		{
			Index< Dim > f2 = f2s[i];
			_row.emplace_back( f2-f1 , scalars.functionIndex(f2) - scalars.functionIndex(f1) );

			// The range of elements supported on both functions
			std::map< ProductIndex , T > row;

			auto f = [&]( Index< Dim > e )
			{
				auto f = [&]( Index< Dim > g1 , Index< Dim > g2 )
				{
					if( g1<=g2 && ( g1!=g2 || Sym ) ) row[ ProductIndex( g1 , g2 ) ] += stencil( e , f1 , f2 , g1 , g2 );
				};
				Basis< Dim >::ElementSupport( e ).template process< 2 >(f);
			};
			Range< Dim >::Intersect( eRange , Basis< Dim >::FunctionSupport(f1) , Basis< Dim >::FunctionSupport(f2) ).process( f );
			_row[i].entries.reserve( row.size() );
			for( auto iter=row.begin() ; iter!=row.end() ; iter++ ) _row[i].entries.push_back( Entry( iter->first.g1-f1 , iter->first.g2-f1 , scalars.functionIndex( iter->first.g1 ) - scalars.functionIndex(f1) , scalars.functionIndex( iter->first.g2 ) - scalars.functionIndex(f1) , iter->second ) );
		}
		_rows[ ScalarFunctions< Dim >::template FullIntegrationStencil< T , 0 >::StencilIndex(f1,_res) ] = _row;
	};
	range.process( f );

	// Move the case df2=0 to the front
	{
		if constexpr( Sym )
			for( unsigned int i=0 ; i<_rows.size() ; i++ ) for( unsigned int j=0 ; j<_rows[i].size() ; j++ ) if( _rows[i][j].df2==0 ) std::swap( _rows[i][j] , _rows[i][0] );
	}

	// Move the cases df2>0 before df2<0
	{
		auto Sort = []( const typename FullIntegrationStencil< T >::Row &r1 , const typename FullIntegrationStencil< T >::Row &r2 )
		{
			return ( r1.df2==0 && r2.df2!=0 ) || ( r1.df2>0 && r2.df2<0 );
		};
		for( unsigned int i=0 ; i<_rows.size() ; i++ ) std::sort( _rows[i].begin() , _rows[i].end() , Sort );
	}
}

template< unsigned int Dim , bool Sym >
template< typename T >
const std::vector< typename ProductFunctions< Dim , Sym >::template FullIntegrationStencil< T >::Row > &ProductFunctions< Dim , Sym >::FullIntegrationStencil< T >::rows( Index< Dim > f1 ) const { return _rows[ ScalarFunctions< Dim >::template FullIntegrationStencil< T , 0 >::StencilIndex( f1 , _res ) ]; }

//////////////////////
// ProductFunctions //
//////////////////////
template< unsigned int Dim , bool Sym >
ProductFunctions< Dim , Sym >::ProductFunctions( unsigned int resolution ) : _r(resolution) , _matrixInfo(resolution)
{
	ScalarFunctions< Dim > scalars( _r );
	Range< Dim > fRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) fRange.first[d] = 0 , fRange.second[d] = _r+1;
}

template< unsigned int Dim , bool Sym >
size_t ProductFunctions< Dim , Sym >::index( Index< Dim > F1 , Index< Dim > F2 , bool &flip ) const
{
	if( F1>F2 )
	{
		std::swap( F1 , F2 );
		flip = true;
	}
	else flip = false;
	for( unsigned int d=0 ; d<Dim ; d++ ) if( F1[d]<0 || F2[d]<0 || F1[d]>(int)_r || F2[d]>(int)_r || F1[d]-F2[d]<-1 || F1[d]-F2[d]>1 ) 
		THROW( "Bad index pair: " , F1 , " <-> " , F2 );
	size_t e = _matrixInfo.entry( F1 , F2 );
	if( e==-1 ) THROW( "Could not find index for: " , F1 , " <-> " , F2 );
	return e;
}

template< unsigned int Dim , bool Sym >
bool ProductFunctions< Dim , Sym >::setIndex( Index< Dim > F1 , Index< Dim > F2 , size_t &i , bool &flip ) const
{
	try{ i = index( F1 , F2 , flip ) ; return true; }
	catch( Exception & ){ return false; }
}

template< unsigned int Dim , bool Sym >
std::vector< std::pair< Index< Dim > , Index< Dim > > > ProductFunctions< Dim , Sym >::indices( void ) const
{
	std::vector< std::pair< Index< Dim > , Index< Dim > > > idx( functionNum() );
	ThreadPool::ParallelFor
		(
			0 , functionNum() , 
			[&]( unsigned int , size_t i )
			{
				Index< Dim > F1 = ScalarFunctions< Dim >::FunctionIndex( i , _r );
				_matrixInfo.process( F1 , [&]( Index< Dim > F2 , size_t e ){ idx[e] = std::make_pair( F1 , F2 ); } );
			}
		);
	return idx;
}

template< unsigned int Dim , bool Sym >
SquareMatrix< double , Dim , Sym > ProductFunctions< Dim , Sym  >::value( const Eigen::VectorXd &xy , Point< double , Dim > p ) const
{
	SquareMatrix< double , Dim , Sym > v;

	p *= _r;
	Index< Dim > E;
	for( unsigned int d=0 ; d<Dim ; d++ ) E[d] = (int)floor( p[d] );
	for( unsigned int d=0 ; d<Dim ; d++ ) if( E[d]==_r ) E[d]--;
	for( unsigned int d=0 ; d<Dim ; d++ ) p[d] -= E[d];

	auto f = [&]( Index< Dim > F1 )
	{
		ElementVector< Dim > g1 = ElementFunction< Dim >( E , F1 ).gradient() * _r;
		auto f = [&]( Index< Dim > F2 , size_t e )
		{
			ElementVector< Dim > g2 = ElementFunction< Dim >( E , F2 ).gradient() * _r;
			ElementProduct< Dim , Sym > g12( g1 , g2 );
			v += g12(p) * xy[e];
		};
		_matrixInfo.process( F1 , f );
	};
	Basis< Dim >::ElementSupport( E ).process( f );
	return v;
}

template< unsigned int Dim , bool Sym >
template< typename T , typename Indexer /* = Hat::BaseIndex< Dim > */ >
T ProductFunctions< Dim , Sym >::operator()( const Indexer &indexer , const FullIntegrationStencil< T > &stencil , const Eigen::VectorXd &x1 , const Eigen::VectorXd &x2 , const Eigen::VectorXd &y1 , const Eigen::VectorXd &y2 ) const
{
	ScalarFunctions< Dim > scalars( (unsigned int)resolution() );
	std::vector< T > integrals( ThreadPool::NumThreads() , T{} );
	Hat::Index< Dim > Off;
	for( unsigned int d=0 ; d<Dim ; d++ ) Off[d] = 1;

	ThreadPool::ParallelFor
	(
		0 , indexer.numFunctions() , 
		[&]( unsigned int t , size_t f1 )
		{
			Window::IsotropicStaticWindow< size_t , Dim , 3 > neighbors = indexer.ffNeighbors( f1 , t );

			Hat::Index< Dim > F1 = indexer.functionIndex( f1 );
			const std::vector< typename FullIntegrationStencil< T >::Row > &rows = stencil.rows( F1 );

			for( unsigned int j=0 ; j<rows.size() ; j++ )
			{
				size_t f2 = neighbors.data[ rows[j]._f2 ];
				if( f2!=-1 )
				{
					if( f1>f2 ) continue;
					double xValue = Coefficient( std::make_pair( x1[f1] , x1[f2] ) , std::make_pair( x2[f1] , x2[f2] ) , f1 , f2 );

					T value{};
					for( unsigned int k=0 ; k<rows[j].entries.size() ; k++ )
					{
						size_t g1 = neighbors.data[ rows[j].entries[k]._g1 ] , g2 = neighbors.data[ rows[j].entries[k]._g2 ];
						if( g1!=-1 && g2!=-1 ) value += rows[j].entries[k].value * Coefficient( std::make_pair( y1[g1] , y1[g2] ) , std::make_pair( y2[g1] , y2[g2] ) , g1 , g2 );
					}
					integrals[t] += value * xValue;
				}
			}
		}
	);

	T integral{};
	for( unsigned int i=0 ; i<integrals.size() ; i++ ) integral += integrals[i];
	return integral;
}

template< unsigned int Dim , bool Sym >
Eigen::SparseMatrix< double > ProductFunctions< Dim , Sym >::systemMatrix( IntegrationStencil< double , 2 , 0 > stencil ) const
{
	FullIntegrationStencil< double > fullStencil( stencil , _r );

	Eigen::SparseMatrix< double > M( functionNum() , functionNum() );

	Eigen::VectorXi rowSizes( (int)functionNum() );

	ThreadPool::ParallelFor
		(
			0 , ScalarFunctions< Dim >::FunctionNum(_r) , 
			[&]( unsigned int , size_t i )
			{
				Index< Dim > F1 = ScalarFunctions< Dim >::FunctionIndex( (size_t)i , _r );

				const std::vector< typename FullIntegrationStencil< double >::Row > &rows = fullStencil.rows( F1 );
				for( unsigned int i=0 ; i<rows.size() ; i++ )
				{
					Hat::Index< Dim > F2 = F1 + rows[i].f2;
					if( F1<=F2 )
					{
						bool flipF;
						rowSizes[ index( F1 , F2 , flipF ) ] = (int)rows[i].entries.size();
					}
				}
			}
		);

	M.reserve( rowSizes );

	ThreadPool::ParallelFor
		(
			0 , ScalarFunctions< Dim >::FunctionNum(_r) , 
			[&]( unsigned int , size_t i )
			{
				Index< Dim > F1 = ScalarFunctions< Dim >::FunctionIndex( (size_t)i , _r );

				const std::vector< typename FullIntegrationStencil< double >::Row > &rows = fullStencil.rows( F1 );
				for( unsigned int i=0 ; i<rows.size() ; i++ )
				{
					Hat::Index< Dim > F2 = F1 + rows[i].f2;
					if( F1<=F2 )
					{
						bool flipF;
						size_t f = index( F1 , F2 , flipF );

						const std::vector< typename FullIntegrationStencil< double >::Entry > &entries = rows[i].entries;
						for( unsigned int j=0 ; j<entries.size() ; j++ )
						{
							Hat::Index< Dim > G1 = F1 + entries[j].g1 , G2 = F1 + entries[j].g2;
							bool flipG;
							size_t g = index( G1 , G2 , flipG );
							if constexpr( Sym ) M.insert( (int)g , (int)f ) = entries[j].value;
							else                M.insert( (int)g , (int)f ) = flipF==flipG ? entries[j].value : -entries[j].value;
						}
					}
				}
			}
		);

	M.makeCompressed();

	return M;
}

template< unsigned int Dim , bool Sym >
typename ProductFunctions< Dim , Sym >::template IntegrationStencil< double , 2 , 0 > ProductFunctions< Dim , Sym >::MassStencil( unsigned int res )
{
	IntegrationStencil< double , 2 , 0 > stencil;

	Index< Dim > e;
	Range< Dim > r = Basis< Dim >::ElementSupport( e );
	auto f = [&]( Index< Dim > i1 , Index< Dim > i2 , Index< Dim > j1 , Index< Dim > j2 ) 
	{
		stencil._values[ Basis< Dim >::template _RelativeIndex<0>(e,i1) ][ Basis< Dim >::template _RelativeIndex<0>(e,i2) ][ Basis< Dim >::template _RelativeIndex<0>(e,j1) ][ Basis< Dim >::template _RelativeIndex<0>(e,j2) ] = Basis< Dim >::template ElementGradientProductMass< Sym >( res , e , std::make_pair( i1 , i2 ) , std::make_pair( j1 , j2 ) );
	};
	r.template process< 4 >( f );

	return stencil;
}

template< unsigned int Dim , bool Sym >
typename ProductFunctions< Dim , Sym >::template IntegrationStencil< SquareMatrix< double , Dim , Sym > , 1 , 0 > ProductFunctions< Dim , Sym >::ValueStencil( unsigned int r )
{
	IntegrationStencil< SquareMatrix< double , Dim , Sym > , 1 , 0 > stencil;
	{
		Index< Dim > E;
		auto f = [&]( Index< Dim > F1 , Index< Dim > F2 )
		{
			unsigned int i1 = Basis< Dim >::template _RelativeIndex< 0 >( E , F1 ) , i2 = Basis< Dim >::template _RelativeIndex< 0 >( E , F2 );
			ElementFunction< Dim > ef1 = ElementFunction< Dim >( E , F1 ) , ef2 = ElementFunction< Dim >( E , F2 );
			ElementProduct< Dim , Sym > product( ef1.gradient() , ef2.gradient() );
			for( unsigned int d=0 ; d<SquareMatrix< double , Dim , Sym >::Coefficients ; d++ )
				stencil[i1][i2][d] = ( product(0,d).integral(r) + product(1,d).integral(r) );
		};
		Basis< Dim >::ElementSupport( E ).template process< 2 >( f );
	}
	return stencil;
}

template< unsigned int Dim , bool Sym >
double ProductFunctions< Dim , Sym >::Coefficient( std::pair< double , double > f1 , std::pair< double , double > f2 , size_t i1 , size_t i2 )
{
	if constexpr( Sym )
	{
		if( i1!=i2 ) return   f1.first * f2.second + f1.second * f2.first    ;
		else         return ( f1.first * f2.second + f1.second * f2.first )/2;
	}
	else return f1.first * f2.second - f1.second * f2.first;
}

template< unsigned int Dim , bool Sym >
Eigen::VectorXd ProductFunctions< Dim , Sym >::product( const Eigen::VectorXd &f1 , const Eigen::VectorXd &f2 ) const
{
	Eigen::VectorXd p( _matrixInfo.entries( false ) );
	size_t fNum = ScalarFunctions< Dim >::FunctionNum( _r );
	for( size_t i1=0 ; i1<fNum ; i1++ )
	{
		Index< Dim > _f1 = ScalarFunctions< Dim >::FunctionIndex( i1 , _r );
		auto f = [&]( Index< Dim > _f2 , size_t e )
		{
			size_t i2 = ScalarFunctions< Dim >::FunctionIndex( _f2 , _r );
			p[e] = Coefficient( std::pair< double , double >( f1[(int)i1] , f1[(int)i2] ) , std::pair< double , double >( f2[(int)i1] , f2[(int)i2] ) , i1 , i2 );
		};
		_matrixInfo.process( _f1 , f );
	}
	return p;
}

template< unsigned int Dim , bool Sym >
Eigen::VectorXd ProductFunctions< Dim , Sym >::productTranspose( const Eigen::VectorXd &x , const Eigen::VectorXd &Z ) const
{
	size_t fNum = ScalarFunctions< Dim >::FunctionNum( _r );
	Eigen::VectorXd p( fNum );

	ThreadPool::ParallelFor
		(
			0 , fNum ,
			[&]( unsigned int , size_t i1 )
			{
				Index< Dim > F1 = ScalarFunctions< Dim >::FunctionIndex( i1 , _r );
				double value = 0;
				auto f = [&]( Index< Dim > F2 , size_t e , bool flip )
					{
						size_t i2 = ScalarFunctions< Dim >::FunctionIndex( F2 , _r );
						value += Coefficient( std::pair< double , double >( 1 , i1==i2 ? 1 : 0 ) , std::pair< double , double >( x[(int)i1] , x[(int)i2] ) , i1 , i2 ) * Z[e] * ( ( flip && !Sym ) ? -1 : 1 );
					};
				_matrixInfo.processAll( F1 , f );
				p[i1] = value;
			}
		);

	return p;
}

template< unsigned int Dim , bool Sym >
Eigen::SparseMatrix< double > ProductFunctions< Dim , Sym >::product( const Eigen::VectorXd &x ) const
{
	ScalarFunctions< Dim > scalars( _r );
	if( x.size()!=scalars.functionNum() ) ERROR_OUT( "Resolutions don't match: " , x.size() , " != " , scalars.functionNum() );

	std::vector< Eigen::Triplet< double > > _triplets;
	_triplets.reserve( _matrixInfo.entries( true ) );
	for( size_t i1=0 ; i1<scalars.functionNum() ; i1++ )
	{
		Index< Dim > f1 = scalars.functionIndex(i1);
		auto f = [&]( Index< Dim > f2 , size_t e )
		{
			size_t i2 = scalars.functionIndex(f2);
			if( Sym )
			{
				if( i1!=i2 )
				{
					_triplets.push_back( Eigen::Triplet< double >( (int)e , (int)i1 , x[i2] ) );
					_triplets.push_back( Eigen::Triplet< double >( (int)e , (int)i2 , x[i1] ) );
				}
				else _triplets.push_back( Eigen::Triplet< double >( (int)e , (unsigned)i2 , x[i1] ) );
			}
			else
			{
				_triplets.push_back( Eigen::Triplet< double >( (int)e , (int)i1 , -x[i2] ) );
				_triplets.push_back( Eigen::Triplet< double >( (int)e , (int)i2 ,  x[i1] ) );
			}
		};
		_matrixInfo.process( f1 , f );
	}

	Eigen::SparseMatrix< double > P( functionNum() , scalars.functionNum() );
	P.setFromTriplets( _triplets.begin() , _triplets.end() );
	_triplets.clear();
	return P;
}

template< unsigned int Dim , bool Sym >
Eigen::SparseMatrix< double > ProductFunctions< Dim , Sym >::prolongation( void ) const
{
	if( _r&1 ) ERROR_OUT( "Expected even resolution: " , _r );
	ProductFunctions coarseProducts( _r>>1 );
	ScalarFunctions< Dim > scalars( _r );
	Eigen::SparseMatrix< double > sP = scalars.prolongation();

#if 1
	std::vector< Eigen::Triplet< double > > triplets;
	ERROR_OUT( "Method not supported" );
#else
	std::vector< Eigen::Triplet< double > > triplets;
	ThreadPool::ParallelFor
		(
			0 , coarseProducts.functionNum() , 
			[&]( unsigned int t , size_t i )
			{
				size_t c1 = coarseProducts._indices[c].first , c2 = coarseProducts._indices[c].second;
				for( Eigen::InnerIterator it1(sP,c1) ; it1 ; ++it1 ) for( Eigen::InnerIterator it2(sP,c2) ; it2 ; ++it2 )
				{
					size_t f1 = (size_t)it1.row() , f2 = (size_t)it2.row();
					bool flip;
					size_t f;
					if( f1!=f2 || Sym )
					{
						size_t f = index( std::make_pair( f1 , f2 ) , flip );
						triplets[t].push_back( Eigen::Triplet< double >( (unsigned int)f , (unsigned int)c , it1.value()*it2.value() * ( !Sym && flip ? -1. : 1. ) ) );
					}
				}
			} ,
			1
		);
#endif
	Eigen::SparseMatrix< double > P( functionNum() , coarseProducts.functionNum() );
	P.setFromTriplets( triplets.begin() , triplets.end() );
	triplets.clear();

	return P;
}

template< unsigned int Dim , bool Sym >
template< typename PiecewiseConstantTensorField /* = std::function< SquareMatrix< double , Dim , Sym > ) ( size_t e ) > */ >
Eigen::VectorXd ProductFunctions< Dim , Sym >::valueDual( PiecewiseConstantTensorField T ) const
{
	Eigen::VectorXd dual( functionNum() );
	Range< Dim > eRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) eRange.first[d] = 0 , eRange.second[d] = _r;
	IntegrationStencil< SquareMatrix< double , Dim , Sym > , 1 , 0 > stencil = ValueStencil( _r );

	ThreadPool::ParallelFor
		(
			0 , ScalarFunctions< Dim >::FunctionNum(_r) , 
			[&]( unsigned int , size_t i1 )
			{
				Index< Dim > F1 = ScalarFunctions< Dim >::FunctionIndex(i1,_r);

				auto f = [&]( Index< Dim > F2 , size_t e )
					{
						dual[e] = 0;
						auto f = [&]( Index< Dim > E )
							{
								// The value of the tensor on the cell
								SquareMatrix< double , Dim , Sym > t = T( ScalarFunctions< Dim >::ElementIndex(E,_r) );

								// The dual associated to the product
								SquareMatrix< double , Dim , Sym > _t = stencil( E , F1 , F2 );
//								for( unsigned int d=0 ; d<SquareMatrix< double , Dim , Sym >::Coefficients ; d++ ) dual[e] += t[d] * _t[d];
								dual[e] += SquareMatrix< double , Dim , Sym >::Dot( t , _t );
							};
						// Iterate over all cells supported by the i-th function
						Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( F1 ) , Basis< Dim >::FunctionSupport( F2 ) , eRange ).process( f );
					};
				_matrixInfo.process( F1 , f );
			}
		);

	return dual;
}

template< unsigned int Dim , bool Sym >
Eigen::VectorXd ProductFunctions< Dim , Sym >::valueDual( const Eigen::VectorXd &xy ) const
{
	IntegrationStencil< double , 2 , 0 >  stencil = MassStencil( _r );
	Eigen::VectorXd dual( functionNum() );
	Range< Dim > eRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) eRange.first[d] = 0 , eRange.second[d] = _r;

	ThreadPool::ParallelFor
		(
			0 , ScalarFunctions< Dim >::FunctionNum(_r) , 
			[&]( unsigned int , size_t i1 )
			{
				Index< Dim > F1 = ScalarFunctions< Dim >::FunctionIndex(i1,_r);

				auto f = [&]( Index< Dim > F2 , size_t e )
					{
						dual[e] = 0;
						auto f = [&]( Index< Dim > E )
							{
								auto f = [&]( Index< Dim > G1 , Index< Dim > G2 )
									{
										if( ( Sym && G1<=G2 ) || (G1<G2) )
										{
											bool flip;
											size_t _e = index( G1 , G2 , flip );
											if( flip ) ERROR_OUT( "WTF!" );
											dual[e] += xy[_e] * stencil( E , F1 , F2 , G1 , G2 );
										}
									};
								// Iterate over all pairs of functions supported on the element
								Basis< Dim >::ElementSupport( E ).template process< 2 >( f );
							};
						// Iterate over all cells supported by the i-th function
						Range< Dim >::Intersect( Basis< Dim >::FunctionSupport( F1 ) , Basis< Dim >::FunctionSupport( F2 ) , eRange ).process( f );
					};
				_matrixInfo.process( F1 , f );
			}
		);

	return dual;
}

template< unsigned int Dim , bool Sym >
Eigen::SparseMatrix< double > ProductFunctions< Dim , Sym >::toMatrix( const Eigen::VectorXd &xy ) const
{
	ScalarFunctions< Dim > scalars( _r );
	Eigen::SparseMatrix< double > M( scalars.functionNum() , scalars.functionNum() );

	// Pre-allocate enough room to store the rows
	Eigen::VectorXi rowSizes( (int)scalars.functionNum() );
	ThreadPool::ParallelFor( 0 , scalars.functionNum() , [&]( unsigned int , size_t i ){ rowSizes[i] = (int)_matrixInfo.entries( scalars.functionIndex( (size_t)i ) , true ); } );
	M.reserve( rowSizes );

	ThreadPool::ParallelFor
		(
			0 , scalars.functionNum() , 
			[&]( unsigned int , size_t i )
			{
				size_t i1 = (size_t)i;
				Index< Dim > f1 = scalars.functionIndex( i1 );

				auto f = [&]( Index< Dim > f2 , size_t e , bool flip )
					{
						size_t i2 = scalars.functionIndex( f2 );
						if constexpr( Sym ) M.insert( (int)i2 , (int)i1 ) = i1==i2 ? xy[e] : xy[e];
						else                M.insert( (int)i2 , (int)i1 ) = flip ? xy[e] : -xy[e];
					};
				_matrixInfo.processAll( f1 , f );
			}
		);

	M.makeCompressed();
	return M;
}

template< unsigned int Dim , bool Sym >
Eigen::VectorXd ProductFunctions< Dim , Sym >::toVector( const Eigen::SparseMatrix< double > &XY ) const
{
	ScalarFunctions< Dim > scalars( _r );
	Eigen::VectorXd xy = Eigen::VectorXd::Zero( functionNum() );

	typename ScalarFunctions< Dim >::template MatrixInfo< 1 , Sym > matrixInfo( _r );

	ThreadPool::ParallelFor
		(
			0 , XY.outerSize() ,
			[&]( unsigned int , size_t i1 )
			{
				Index< Dim > f1 = scalars.functionIndex( i1 );
				for( Eigen::InnerIterator it(XY,(int)i1) ; it ; ++it )
				{
					size_t i2 = (size_t) it.row();
					Index< Dim > f2 = scalars.functionIndex(i2);
					if constexpr( Sym )
					{
						if( f1==f2 ) xy[ matrixInfo.entry(f1,f2) ] = it.value();
						else if( f1<f2 ) xy[ matrixInfo.entry(f1,f2) ] = it.value();
					}
					else if( f1<=f2 ) xy[ matrixInfo.entry(f1,f2) ] = it.value();
				}
			}
		);
	return xy;
}

/////////////////////////////////////////////
// ProductFunctions< Dim , Sym >::UnitTest //
/////////////////////////////////////////////
template< unsigned int Dim , bool Sym >
void ProductFunctions< Dim , Sym >::UnitTest::operator()( void ) const
{
	evaluation();
	product();
	integrationStencils();
	massAndDual();
	matrixVector();
}

template< unsigned int Dim , bool Sym >
void ProductFunctions< Dim , Sym >::UnitTest::evaluation( void ) const
{
	{
		double err = 0;

		Eigen::VectorXd x = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars );
		Eigen::VectorXd y = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars );
		Eigen::VectorXd xy = _products.product(x,y);
		for( unsigned int n=0 ; n<_N ; n++ )
		{
			Point< double , Dim > p;
			for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = Random< double >();
			SquareMatrix< double , Dim , Sym > v1( _scalars.gradient( x , p ) , _scalars.gradient( y , p ) );
			SquareMatrix< double , Dim , Sym > v2 = _products.value( xy , p );
			err += ( v1 - v2 ).squareNorm();
		}
		err = sqrt( err / _N );
		std::cout << "\tEvaluation error 1: " << err << std::endl;
	}

	{
		double err = 0;

		// W(i,j) = ( \nabla\phi_i . \nabla\phi_j )
		std::vector< std::pair< Index< Dim > , Index< Dim > > > indices = _products.indices();
		Eigen::VectorXd x = Eigen::VectorXd::Zero( _scalars.functionNum() );
		Eigen::VectorXd y = Eigen::VectorXd::Zero( _scalars.functionNum() );
		Eigen::VectorXd xy = Eigen::VectorXd::Zero( _products.functionNum() );
		for( unsigned int n=0 ; n<_N ; n++ )
		{
			unsigned int idx = rand()%_products.functionNum();
			size_t f1 = _scalars.functionIndex( indices[idx].first );
			size_t f2 = _scalars.functionIndex( indices[idx].second );

			x[f1] = y[f2] = xy[idx] = 1;
			for( unsigned int i=0 ; i<_numSamples ; i++ )
			{
				Point< double ,  Dim > p;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = Random< double >();
				SquareMatrix< double , Dim , Sym > v1 = _products.value( xy , p );
				SquareMatrix< double , Dim , Sym > v2( _scalars.gradient(x,p) , _scalars.gradient(y,p) );
				err += (v1-v2).squareNorm();
			}
			x[f1] = y[f2] = xy[idx] = 0;
		}
		std::cout << "\tEvaluation error 2: " << err << std::endl;
	}
}

template< unsigned int Dim , bool Sym >
void ProductFunctions< Dim , Sym >::UnitTest::product( void ) const
{
	std::vector< std::pair< Index< Dim > , Index< Dim > > > indices = _products.indices();
	double err1 = 0 , err2 = 0 , err3 = 0 , err4 = 0;
	Eigen::SparseMatrix< double > mass = _products.mass();
	Eigen::VectorXd z( _products.functionNum() );
	for( unsigned int n=0 ; n<_N ; n++ )
	{
		Eigen::VectorXd x = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars );
		Eigen::VectorXd y = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars );
		Eigen::SparseMatrix< double > X = _products.product(x);
		Eigen::SparseMatrix< double > Y = _products.product(y);
		for( unsigned int i=0 ; i<_products.functionNum() ; i++ ) z[i] = Random< double >()*2-1;

		{
			Eigen::VectorXd v1 = _products.product( x , y );
			Eigen::VectorXd v2 = X * y;
			err1 += (v1-v2).squaredNorm();
		}

		Eigen::VectorXd xy = _products.product( x , y );
		for( unsigned int i=0 ; i<_products.functionNum() ; i++ )
		{
			size_t i1 = _scalars.functionIndex( indices[i].first ) , i2 = _scalars.functionIndex( indices[i].second );
			double v1 = xy[i];
			double v2;
			if( i1==i2 ) v2 = Sym ? ( x[i1] * y[i2] + x[i2] * y[i1] )/2 : ( x[i1] * y[i2] - x[i2] * y[i1] )/2;
			else         v2 = Sym ? ( x[i1] * y[i2] + x[i2] * y[i1] ) : ( x[i1] * y[i2] - x[i2] * y[i1] );
			err2 += (v1-v2) * (v1-v2);
		}
		{
			double v1 = xy.dot(z);
			double v2 = ( _products.toMatrix( z )*x ).dot( y );
			err3 += (v1-v2)*(v1-v2);
		}

		{
			Eigen::VectorXd v1 = _products.productTranspose( x , z );
			Eigen::VectorXd v2 = X.transpose() * z;
			err4 += (v1-v2).squaredNorm();
		}
	}
	std::cout << "\tProduct error: " << err1 << " / " << err2 << " / " << err3 << " / " << err4 << std::endl;
}

template< unsigned int Dim , bool Sym >
void ProductFunctions< Dim , Sym >::UnitTest::integrationStencils( void ) const
{
	Range< Dim > sampleRange;
	for( unsigned int d=0 ; d<Dim ; d++ ) sampleRange.second[d] = _numSamples;

	std::cout << "\tDiscrete integration:" << std::endl;
	Index< Dim > E;
	// Test single value integration
	{
		IntegrationStencil< SquareMatrix< double , Dim , Sym > , 1 , 0 > stencil = ValueStencil( _scalars.resolution() );

		double err = 0;

		auto f = [&]( Index< Dim > F1 , Index< Dim > F2 )
		{
			if constexpr( Sym ){ if( F1> F2 ) return; }
			else               { if( F1>=F2 ) return; }

			MishaK::SquareMatrix< double , Dim > v1 = stencil( E , F1 , F2 )() , v2;

			auto f = [&]( Index< Dim > I )
			{
				Point< double , Dim > p;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = ( I[d] + 0.5 ) / _numSamples;
				v2 += Hat::ElementProduct< Dim , Sym >( Hat::ElementFunction< Dim >(E,F1).gradient() , Hat::ElementFunction< Dim >(E,F2).gradient() )( p ) * _scalars.resolution() * _scalars.resolution();
			};
			sampleRange.process( f );
			for( unsigned int d=0 ; d<Dim ; d++ ) v2 /= _numSamples * _scalars.resolution();
			err += ( v1-v2 ).squareNorm();
		};
		Basis< Dim >::ElementSupport( E ).template process< 2 >( f );

		std::cout << "\t\tValue: " << err << std::endl;
	}
}

template< unsigned int Dim , bool Sym >
void ProductFunctions< Dim , Sym >::UnitTest::massAndDual( void ) const
{
	// Test the value dual
	{
		Eigen::SparseMatrix< double > mass = _products.mass();
		IntegrationStencil< double , 2 , 0 > stencil = MassStencil( _scalars.resolution() );
		FullIntegrationStencil< double > fullStencil( stencil , _scalars.resolution() );

		double err = 0;
		for( unsigned int n=0 ; n<_N ; n++ )
		{
			Eigen::VectorXd x1 = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars );
			Eigen::VectorXd y1 = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars );
			Eigen::VectorXd x2 = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars );
			Eigen::VectorXd y2 = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars );
			Eigen::VectorXd xy1 = _products.product( x1 , y1 ) , xy2 = _products.product( x2 , y2 );
			Eigen::VectorXd _xy1 = _products.valueDual( xy1 );
			double v[] = { _xy1.dot( xy2 ) , ( mass * xy1 ).dot( xy2 ) , _products( stencil , xy1 , xy2 ) , _products( fullStencil , xy1 , xy2 ) , _products( fullStencil , x1 , y1 , x2 , y2 ) };
			std::pair< double , double > avg_and_dev = Miscellany::AverageAndStandardDeviation( v , sizeof(v)/sizeof(double) );
			err += avg_and_dev.second * avg_and_dev.second;
		}
		std::cout << "\tMass/dual: " << err << std::endl;
	}
	// Test the piecewise-constant value dual
	{
		std::vector< std::pair< Index< Dim > , Index< Dim > > > indices = _products.indices();
		std::vector< SquareMatrix< double , Dim , Sym > > T( _scalars.elementNum() );
		for( size_t i=0 ; i<_scalars.elementNum() ; i++ ) T[i] = SquareMatrix< double , Dim , Sym >( RandomBallPoint< double , Dim >() , RandomBallPoint< double , Dim >() );
		Eigen::SparseMatrix< double > stiffness = _scalars.stiffness( [&]( size_t e ){ return T[e](); } );
		Eigen::VectorXd dual = _products.valueDual( [&]( size_t e ){ return T[e]; } );

		Eigen::VectorXd x = Eigen::VectorXd::Zero( _scalars.functionNum() ) , y = Eigen::VectorXd::Zero( _scalars.functionNum() );
		double err1 = 0 , err2 = 0;
		for( unsigned int n=0 ; n<_N ; n++ )
		{
			size_t idx = rand() % _products.functionNum();
			size_t i1 = _scalars.functionIndex( indices[idx].first );
			size_t i2 = _scalars.functionIndex( indices[idx].second );
			x[i1] = 1 , y[i2] = 1;

			{
				double v1 = dual[idx];
				double v2 = ( stiffness * x ).dot(y);
				double v3 = Sym ? ( stiffness * y ).dot(x) : -( stiffness * y ).dot(x);
				err1 += (v1-v2)*(v1-v2);
				err2 += (v1-v3)*(v1-v3);
			}

			x[i1] = 0 , y[i2] = 0;
		}
		std::cout << "\tStiffness/mass: " << err1 << " / " << err2 << std::endl;
	}
	// Test the piecewise-constant value dual
	{
		double err = 0;
		std::vector< SquareMatrix< double , Dim , Sym > > T( _scalars.elementNum() );
		for( unsigned int n=0 ; n<_N ; n++ )
		{
			for( size_t i=0 ; i<_scalars.elementNum() ; i++ ) T[i] = SquareMatrix< double , Dim , Sym >( RandomBallPoint< double , Dim >() , RandomBallPoint< double , Dim >() );
			Eigen::SparseMatrix< double > v1 = _scalars.stiffness( [&]( size_t e ){ return T[e](); } );
			Eigen::SparseMatrix< double > v2 = _products.toMatrix( _products.valueDual( [&]( size_t e ){ return T[e]; } ) );
			err += (v1-v2).squaredNorm();
		}
		std::cout << "\tStiffness/mass: " << err << std::endl;
	}
}

template< unsigned int Dim , bool Sym >
void ProductFunctions< Dim , Sym >::UnitTest::matrixVector( void ) const
{
	double err1 = 0 , err2 = 0 , err3 = 0 , err4 = 0 , err5 = 0 , err6 = 0;

	Eigen::SparseMatrix< double > mass = _products.mass();
	Eigen::VectorXd Z( _products.functionNum() );
	for( unsigned int n=0 ; n<_N ; n++ )
	{
		Eigen::VectorXd x1 = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars );
		Eigen::VectorXd y1 = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars );
		Eigen::VectorXd x2 = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars );
		Eigen::VectorXd y2 = ScalarFunctions< Dim >::UnitTest::RandomCoefficients( _scalars );
		for( size_t i=0 ; i<_products.functionNum() ; i++ ) Z[i] = Random< double >()*2-1;

		Eigen::VectorXd xy1 = _products.product( x1 , y1 );
		Eigen::VectorXd xy2 = _products.product( x2 , y2 );
		Eigen::SparseMatrix< double > XY1 = _products.toMatrix( xy1 );
		Eigen::SparseMatrix< double > _XY1 = _products.toMatrix( _products.valueDual(xy1) );
		err1 += ( xy1 - _products.toVector( XY1 ) ).squaredNorm();
		{
			double v1 = xy1.dot(xy2) , v2 = (XY1*x2).dot(y2);
			err2 += (v1-v2) * (v1-v2);
		}
		{
			Eigen::SparseMatrix< double > A_Z = _products.toMatrix( Z );
			double v1 = xy1.dot( Z );
			double v2 = ( A_Z * x1 ).dot(y1);
			err3 += (v1-v2) * (v1-v2);
		}
		{
			double v1 = (mass*xy1).dot(xy2);
			double v2 = (_XY1*x2).dot(y2);
			err4 += (v1-v2) * (v1-v2);
		}
		{
			Eigen::SparseMatrix< double > X = _products.product(x1);
			{
				Eigen::VectorXd v1 = X.transpose() * Z;
				Eigen::VectorXd v2 = _products.toMatrix( Z ) * x1;
				err5 += (v1-v2).squaredNorm();
			}
		}
	}

	std::cout << "\t     Matrix <-> vector: " << err1 << std::endl;
	std::cout << "\t     Matrix/vector dot: " << err2 << std::endl;
	std::cout << "\t     Vector/matrix dot: " << err3 << std::endl;
	std::cout << "\tMatrix/vector dual dot: " << err4 << std::endl;
	std::cout << "\t                   ...: " << err5 << std::endl;
}
