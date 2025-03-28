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

#include <vector>
#include <unordered_map>
#include <Eigen/Eigenvalues>
#include "Misha/Miscellany.h"
#include "Misha/Geometry.h"
#include "Misha/RegularGrid.h"
#include "Misha/Algebra.h"
#include "Misha/MultiThreading.h"
#include <Misha/RegularTree.h>
#include "ExteriorPoissonSolvers.h"
#include "../Include/GridSamples.h"

namespace MishaK
{
	namespace ExteriorPoisson
	{
		// Note that the functions here take values in the unit cube whereas regular grids take values in the range defined by the resolution
		static const unsigned int CoDim = 2;

		// Computes the transformation mapping the point set into the unit cube
		template< unsigned int Dim , typename SampleFunctor /* = std::function< Point< double , Dim > > ( unsigned int idx ) */ >
		SquareMatrix< double , Dim+1 > ToUnitCube( SampleFunctor F , size_t sampleNum , double scale );

		// Computes the blade whose wedge product best fits the skew-symmetric matrix
		template< unsigned int Dim >
		std::pair< Point< double , Dim > , Point< double , Dim > > FitBlade( Hat::SkewSymmetricMatrix< double , Dim > S , unsigned int &steps , double eps );

		template< unsigned int Dim , typename Energy >
		struct Reconstructor
		{
			template< typename HermiteSampleFunctor /* = std::function< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > ( unsigned int ) >*/ , typename HierarchicaFunctionlIndexer /* = HierarchicalIndexer< Dim > */ >
			Reconstructor
			(
				const HierarchicaFunctionlIndexer & hierarchicalIndexer ,
				unsigned int depth , unsigned int minSolveDepth , unsigned int maxSolveDepth ,
				HermiteSampleFunctor F , const OrderedSampler< Dim > &orderedSampler ,
				const GridSamples::Estimator< Dim > &estimator ,
				std::function< double ( double ) > noiseToWeight ,
				double dWeight , double sWeight ,
				unsigned int kernelRadius ,
				unsigned int projectionIterations ,
				double projectionEpsilon ,
				bool verbose
			);

			void solve
				(
					bool singleLevel ,
					double iterationMultiplier ,
					unsigned int gsIterations ,
					bool separate ,
					bool verbose
					);

			const Eigen::VectorXd & x( unsigned int depth ) const { return _x[depth]; }
			const Eigen::VectorXd & y( unsigned int depth ) const { return _y[depth]; }
			const std::vector< Eigen::VectorXd > & x( void ) const { return _x; }
			const std::vector< Eigen::VectorXd > & y( void ) const { return _y; }

		protected:
			unsigned int _depth , _minSolveDepth , _maxSolveDepth;
			Hat::Range< Dim > _fRange;
			std::vector< Energy > _energies;
			std::vector< Eigen::VectorXd > _x , _y;
			std::vector< std::vector< std::vector< unsigned int > > > _mcIndices;
		};


		////////////////////
		// Implementation //
		////////////////////
		template< unsigned int Dim , typename SampleFunctor /* = std::function< Point< double , Dim > > ( unsigned int idx ) */ >
		SquareMatrix< double , Dim+1 > ToUnitCube( SampleFunctor F , size_t sampleNum , double scale )
		{
			if( !sampleNum ) MK_ERROR_OUT( "Expected some samples: " , sampleNum );

			SquareMatrix< double , Dim+1 > xForm = SquareMatrix< double , Dim+1 >::Identity();
			Point< double , Dim > bBox[2];
			bBox[0] = bBox[1] = F(0);
			for( unsigned int i=1 ; i<sampleNum ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ )
			{
				Point< double , Dim > p = F(i);
				bBox[0][d] = std::min< double >( bBox[0][d] , p[d] );
				bBox[1][d] = std::max< double >( bBox[1][d] , p[d] );
			}
			double _scale = bBox[1][0] - bBox[0][0];
			for( unsigned int d=1 ; d<Dim ; d++ ) _scale = std::max< double >( _scale , bBox[1][d] - bBox[0][d] );

			SquareMatrix< double , Dim+1 > t1 = SquareMatrix< double , Dim+1 >::Identity() , s = SquareMatrix< double , Dim+1 >::Identity() , t2 = SquareMatrix< double , Dim+1 >::Identity();
			for( unsigned int d=0 ; d<Dim ; d++ ) t2(Dim,d) = -( bBox[0][d] + bBox[1][d] ) / 2.;
			for( unsigned int d=0 ; d<Dim ; d++ ) s(d,d) = 1. / ( _scale * scale );
			for( unsigned int d=0 ; d<Dim ; d++ ) t1(Dim,d) = 0.5;
			xForm = t1 * s * t2;
			return xForm;
		}


		template< unsigned int Dim , typename Energy >
		template< typename HermiteSampleFunctor /* = std::function< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > ( unsigned int ) >*/ , typename HierarchicalIndexer /* = Hat::HierarchicalIndexer< Dim > */ >
		Reconstructor< Dim , Energy >::Reconstructor
		(
			const HierarchicalIndexer & hierarchicalIndexer ,
			unsigned int depth , unsigned int minSolveDepth , unsigned int maxSolveDepth ,
			HermiteSampleFunctor F , const OrderedSampler< Dim > &orderedSampler ,
			const GridSamples::Estimator< Dim > &estimator ,
			std::function< double ( double ) > noiseToWeight ,
			double dWeight , double sWeight ,
			unsigned int kernelRadius ,
			unsigned int projectionIterations ,
			double projectionEpsilon ,
			bool verbose
		) : _depth(depth) , _minSolveDepth(minSolveDepth) , _maxSolveDepth(maxSolveDepth)
		{
			using Indexer = typename HierarchicalIndexer::Indexer;
			static_assert( std::is_base_of_v< Hat::BaseHierarchicalIndexer< Dim > , HierarchicalIndexer > , "[ERROR] Poorly formed HierarchicalIndexer" );
			static_assert( std::is_base_of_v< SystemEnergy< Dim , Indexer > , Energy > , "[ERROR] Poorly formed Energy" );
			using namespace Hat;
			for( unsigned int d=0 ; d<Dim ; d++ ) _fRange.first[d] = 0 , _fRange.second[d] = (1<<_depth) + 1;

			_x.resize( _depth+1 );
			_y.resize( _depth+1 );
			_mcIndices.resize( _depth+1 );

			{
				for( unsigned int d=0 ; d<=_depth ; d++ )
				{
					auto indexer = hierarchicalIndexer[d];
					ScalarFunctions< Dim > scalars( 1<<d );

					_x[d].resize( indexer.numFunctions() );
					_y[d].resize( indexer.numFunctions() );

					_mcIndices[d].resize( 1<<Dim );

					for( unsigned int f=0 ; f<indexer.numActiveFunctions() ; f++ )
					{
						Index< Dim > F = indexer.functionIndex( f );
						unsigned mcIndex = 0;
						for( unsigned int d=0 ; d<Dim ; d++ ) if( F[d]&1 ) mcIndex |= (1<<d);
						_mcIndices[d][ mcIndex ].push_back( f );
					}
				}
			}

			Miscellany::PerformanceMeter pMeter;

			Hat::ScalarFunctions< Dim > scalars( 1<<_depth );
			Indexer indexer = hierarchicalIndexer[_depth];
			GridSamples::TreeSplatter< Dim , SkewSymmetricMatrix< double , Dim > > treeSplatter( _depth );

			// [BEGIN] Populate the skew-symmetric grid
			{
				auto _F = [&]( size_t idx , unsigned int thread )
					{
						std::pair< Point< double , Dim > , SkewSymmetricMatrix< double , Dim > > sample = F(idx);
						sample.second /= sqrt( sample.second.squareNorm() );
						// Have the magnitude of the samples contribution be proportional to its measure
						sample.second *= estimator.measure( sample.first , thread );
						sample.second *= 1<<(_depth*CoDim);
						//				sample.second /= (1<<(_depth*Dim));
						return sample;
					};


				{
					std::vector< std::tuple< Point< double , Dim > , SkewSymmetricMatrix< double , Dim > , double > > samples;
					{
						size_t sampleNum = 0;
						for( size_t i=0 ; i<orderedSampler.size() ; i++ ) sampleNum += orderedSampler[i].second.size();
						samples.resize( sampleNum );
						ThreadPool::ParallelFor
						(
							0 , orderedSampler.size() ,
							[&]( unsigned int t , size_t i )
							{
								const std::vector< size_t > &subSampleIndices = orderedSampler[i].second;
								for( size_t j=0 ; j<subSampleIndices.size() ; j++ )
								{
									std::pair< Point< double , Dim > , SkewSymmetricMatrix< double , Dim > > sample = _F( subSampleIndices[j] , t );
									std::get< 0 >( samples[ subSampleIndices[j] ] ) = sample.first;
									std::get< 1 >( samples[ subSampleIndices[j] ] ) = sample.second;
									std::get< 2 >( samples[ subSampleIndices[j] ] ) = estimator.depth( sample.first , t );
								}
							}
						);
					}
					auto DepthWeight = []( unsigned int depth ) { return (double)( 1<<(depth*Dim) ); };
					auto __F = [&]( size_t idx ){ return std::pair< Point< double , Dim > , SkewSymmetricMatrix< double , Dim > >( std::get< 0 >( samples[idx] ) , std::get< 1 >( samples[idx] ) ); };
					auto __D = [&]( size_t idx ){ return std::get< 2 >( samples[idx] ); };

					switch( kernelRadius )
					{
					case 0: treeSplatter.template addSamples< 0 >( __F , orderedSampler , __D , DepthWeight , false ) ; break;
					case 1: treeSplatter.template addSamples< 1 >( __F , orderedSampler , __D , DepthWeight , false ) ; break;
					case 2: treeSplatter.template addSamples< 2 >( __F , orderedSampler , __D , DepthWeight , false ) ; break;
					case 3: treeSplatter.template addSamples< 3 >( __F , orderedSampler , __D , DepthWeight , false ) ; break;
					case 4: treeSplatter.template addSamples< 4 >( __F , orderedSampler , __D , DepthWeight , false ) ; break;
					default: MK_ERROR_OUT( "Kernel radius must be in range [0,4]: " , kernelRadius );
					}
				}
			}
			// [END] Populate the ske-symmetric grid

			// Assume that we are given a d-dimensional manifold M living in \R^n, with a map w:M -> \wedge^{n-d}\R^n
			// We define an alternating product field on \R^n via a unit-norm kernel k as:
			//		W_l(p) = \int_M k(p-q) * w(q) dq
			// The square norm of W is:
			//		|| W_l ||^2 = \int_{\R^n} ( \int_M k_q(p) * w(q) ) dq )^2 dp
			// Transitioning to the next level effects this in two ways:
			//		1. w    -> 2^{n-d} * w:		to account for the fact that we want the gradients to be twice as large
			//		2. k(p) -> 2^n * k(p*2):	to account for the fact that the kernel gets twice as narrow but is still unit-norm
			// This gives:
			//		W_{l+1}(p) = \int_M 2^n * k( (p-q)*2 ) * w(q) * 2^{n-d}
			//		           = 2^{n-d} * 2^n * \int_M k( (p-q)*2 ) * w(q)
			// Assuming, for simplicity, that w(q) is constant
			//		           ~ 2^{n-d} * 2^n * \int_M k( (p-q)*2 )
			// Then the square norm is:
			//		|| W_{l+1} ||^2 ~ 2^{2n-2d} * 2^{2n} * \int_{\R^n} ( \int_M k( (p-q)*2 ) dq )^2 dp
			// Assuming, further, that k is the indicator function of a unit ball gives:
			//		                ~ 2^{2n-2d} * 2^{2n} * \int_{\R^n} ( |B_0.5^n(p) \cap M| )^2 dp
			//		                ~ 2^{2n-2d} * 2^{2n} * \int_{ p | d(p,M)<0.5 } |B_0.5^n(p) \cap M|^2 dp
			// Assuming that M is flat, we have:
			//		|B_0.5^n(p) \cap M| = |B_{sqrt( (0.5)^2 - d^2(p,M) )}}^d|
			//		                    ~ 2^{-d} 
			// This gives:
			//		|| W_{l+1} ||^2 ~ 2^{2n-2d} * 2^{2n} * \int_{ p | d(p,M)<0.5 } (2^{-d})^2
			//		                = 2^{2n-2d} * 2^{2n} * 2^{-2d} * \int_{ p |d(p,M)<0.5 }
			//		                ~ 2^{2n-2d} * 2^{2n-2d} * |M| * 2^{d-n}
			//		                ~ 2^{2n-2d} * 2^{n-d}
			// So that:
			//		|| W_{l+1} || ~ sqrt( 2^{2n-2d} * 2^{n-d} )
			double scale = pow( pow( 0.5 , _depth ) , 1.5 * CoDim );
			{
				auto F = [&]( unsigned int depth , Hat::Index< Dim > I , SkewSymmetricMatrix< double , Dim > &s ){ s *= scale; }; 
				treeSplatter.process( F );
			}

			if( verbose ) std::cout << pMeter( "Splatted" ) << std::endl;
			Eigen::SparseMatrix< double > R( indexer.numFunctions() , indexer.numFunctions() );

			if( sWeight>0 )
			{
				double measure = 0;
				std::atomic< size_t > sampleNum = 0;
				ThreadPool::ParallelFor
				(
					0 , orderedSampler.size() ,
					[&]( unsigned int t , size_t i )
					{
						const std::vector< size_t > &indices = orderedSampler[i].second;
						sampleNum += indices.size();
						double m = 0;
						for( size_t j=0 ; j<indices.size() ; j++ ) m += estimator.measure( F( indices[j] ).first , t );
						AddAtomic( measure , m );
					}
				);
				measure /= sampleNum;
				R = scalars.deltaMass( indexer , [&]( size_t idx ){ return F(idx).first; } , orderedSampler , [&]( Point< double , Dim > p , unsigned int t ){ return noiseToWeight( estimator.noise( p , t ) ); } ) * sWeight * measure;
			}

			if( verbose ) std::cout << pMeter( "Regularization" ) << std::endl;
			{
				unsigned int res = 1<<_depth;
				unsigned int _res[Dim];
				for( unsigned int d=0 ; d<Dim ; d++ ) _res[d] = res;
				auto F = [&]( Hat::Index< Dim > I , unsigned int t )
					{
						Point< double , Dim > p;
						for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = ( I[d] + 0.5 ) / res;
						SkewSymmetricMatrix< double , Dim > s = treeSplatter( p , t );

						// Possibly, project onto the space of blades
						if( projectionIterations )
						{
							std::pair< Point< double , Dim > , Point< double , Dim > > blade = FitBlade< Dim >( s , projectionIterations , projectionEpsilon );
							s = Hat::SkewSymmetricMatrix< double , Dim >( blade.first , blade.second );
						}
						return s;
					};
				_energies.emplace_back( hierarchicalIndexer[_depth] , 1<<_depth , F , dWeight , R );
			}
			if( verbose ) std::cout << pMeter( "Finest constraint" ) << std::endl;
			for( unsigned int d=0 ; d<_depth ; d++ ) _energies.emplace( _energies.begin() , _energies[0].restrict( hierarchicalIndexer[_depth-1-d] ) );
			if( verbose ) std::cout << pMeter( "Coarse constraints" ) << std::endl;
		}

		template< unsigned int Dim , typename Energy >
		void Reconstructor< Dim , Energy >::solve
			(
				bool singleLevel ,
				double iterationMultiplier ,
				unsigned int gsIterations ,
				bool separate ,
				bool verbose
				)
		{
			using namespace Hat;

			Miscellany::PerformanceMeter pMeter;

			auto Iters = [&]( unsigned int iters , unsigned int depth ){ return (int)ceil( pow( iterationMultiplier , _depth - depth ) * iters ); };

			auto MCIndices = [&]( unsigned int d ) -> const std::vector< std::vector< unsigned int > > & { return _mcIndices[d]; };

			{
				for( unsigned int i=0 ; i<_x[ _minSolveDepth ].size() ; i++ ) _x[ _minSolveDepth ][i] = Random< double >() , _y[ _minSolveDepth ][i] = Random< double >();
				double targetSquareNorm = _energies.back().squareNorm();
				// Rescale the initial guess so that it is on the order of the target skew-symmetric matrix field
				{
					double scale = _energies[ _minSolveDepth ].wedgeSquareNorm( _x[ _minSolveDepth ] , _y[ _minSolveDepth ] );
					scale = sqrt( sqrt( targetSquareNorm / scale ) );
					_x[ _minSolveDepth ] *= scale , _y[ _minSolveDepth ] *= scale;
				}
				if( singleLevel )
				{
					// Prolong from coarsest to finest and clear all the coarser levels
					for( unsigned int d=_minSolveDepth ; d<_maxSolveDepth ; d++ ) _x[d+1] = _energies[d].prolongation( _x[d] ) , _y[d+1] = _energies[d].prolongation( _y[d] );
				}
			}

			// Clear all levels except for the initial guess
			if( singleLevel ) for( unsigned int d=_minSolveDepth ; d<_maxSolveDepth ; d++ ) _x[d].setZero() , _y[d].setZero();
			else              for( unsigned int d=_minSolveDepth+1 ; d<=_maxSolveDepth ; d++ ) _x[d].setZero() , _y[d].setZero();

			double _err = _energies[ _maxSolveDepth ].squareNorm();

			if( singleLevel )
			{
				unsigned int d = _maxSolveDepth;
				if( separate ) ExteriorPoissonSolvers::SeparateGaussSeidelSolve< Dim >( d , MCIndices( d ) , _energies[d] , _x[d] , _y[d] , gsIterations , verbose );
				else           ExteriorPoissonSolvers::  PairedGaussSeidelSolve< Dim >( d , MCIndices( d ) , _energies[d] , _x[d] , _y[d] , gsIterations , verbose );
			}
			else
			{
				double err = 0;
				if( verbose ) err = _energies[ _maxSolveDepth ]( _x[ _maxSolveDepth ] , _y[ _maxSolveDepth ] );

				{
					// Fine-to-coarse
					for( unsigned int d=_maxSolveDepth ; d+1!=_minSolveDepth ; d-- )
					{
						if( d>_minSolveDepth ) _x[d].setZero() , _y[d].setZero();

						// Restrict to the next resolution
						if( d>_minSolveDepth ) _energies[d-1].update( _energies[d] , _x[d] , _y[d] );
					}

					// Coarse-to-fine
					for( unsigned int d=_minSolveDepth ; d<=_maxSolveDepth ; d++ )
					{
						// Solve at the current resolution
						{
							if( separate ) ExteriorPoissonSolvers::SeparateGaussSeidelSolve< Dim >( d , MCIndices( d ) , _energies[d] , _x[d] , _y[d] , Iters( gsIterations , d ) , verbose );
							else           ExteriorPoissonSolvers::  PairedGaussSeidelSolve< Dim >( d , MCIndices( d ) , _energies[d] , _x[d] , _y[d] , Iters( gsIterations , d ) , verbose );
						}

						// Prolong to the next resolution
						if( d<_maxSolveDepth ) _x[d+1] += _energies[d].prolongation( _x[d] ) , _y[d+1] += _energies[d].prolongation( _y[d] );
					}
				}
			}

			std::pair< double , double > errs = _energies[ _maxSolveDepth ].energies( _x[ _maxSolveDepth ] , _y[ _maxSolveDepth ] );
			{
				std::stringstream sStream1 , sStream2;
				sStream1 << "Error";
				sStream2 << _err << " -> " << errs.first + errs.second << " = " << errs.first << " + " << errs.second << " => " << sqrt( (errs.first+errs.second )/_err );
				std::cout << pMeter( sStream1.str() , sStream2.str() ) << std::endl;;
			}

			// Prolong from the finest solve depth to the finest depth
			for( unsigned int d=_maxSolveDepth ; d<_depth ; d++ ) _x[d+1] = _energies[d].prolongation( _x[d] ) , _y[d+1] = _energies[d].prolongation( _y[d] );


			// Remove the part from the higher frequency that was prolonged from the lower frequency
			for( unsigned int d=_depth ; d>0 ; d-- )
			{
				_x[d] -= _energies[d-1].prolongation( _x[d-1] );
				_y[d] -= _energies[d-1].prolongation( _y[d-1] );
			}
		}

		template< unsigned int Dim >
		std::pair< Point< double , Dim > , Point< double , Dim > > FitBlade( Hat::SkewSymmetricMatrix< double , Dim > S , unsigned int &steps , double eps )
		{
			if( S.squareNorm()==0 ) return std::make_pair( Point< double , Dim >() , Point< double , Dim >() );

			Point< double , Dim > v1 , v2;
			v1 = RandomSpherePoint< double , Dim >();
			v1 /= sqrt( v1.squareNorm() );

			auto Next = [&]( Point< double , Dim > v )
				{
					Point< double , Dim > _v = S() * v;
					_v -= Point< double , Dim >::Dot( v , _v ) * v;
					return _v / sqrt( _v.squareNorm() );
				};

			for( unsigned int s=0 ; s<steps ; s++ )
			{
				v2 = Next( v1 );

				Point< double , Dim > _v1 = -Next( v2 );

				if( ( v1 - _v1 ).squareNorm() < eps ){ steps = s ; break; }
				v1 = _v1;
			}
			double scale = sqrt( sqrt( S.squareNorm() ) );
			return std::make_pair( v1*scale , v2*scale );
		}
	}
}