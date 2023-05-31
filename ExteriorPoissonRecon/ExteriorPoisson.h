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

#ifndef EXTERIOR_POISSON_INCLUDED
#define EXTERIOR_POISSON_INCLUDED

#include <vector>
#include <unordered_map>
#include <Eigen/Eigenvalues>
#include "Misha/Miscellany.h"
#include "Misha/Geometry.h"
#include "Misha/RegularGrid.h"
#include "Misha/Algebra.h"
#include "Include/GridSamples.h"
#include "ExteriorPoissonSolvers.h"

namespace ExteriorPoisson
{
	// Note that the functions here take values in the unit cube whereas regular grids take values in the range defined by the resolution
	static const unsigned int CoDim = 2;

	// Computes the transformation mapping the point set into the unit cube
	template< unsigned int Dim , typename SampleFunctor /* = std::function< Point< double , Dim > > ( unsigned int idx ) */ >
	SquareMatrix< double , Dim+1 > ToUnitCube( SampleFunctor F , size_t sampleNum , double scale );

	// Outputs the header, running time, and peak memory usage
	void OutputPerformance( std::string header , const Miscellany::Timer &timer )
	{
		std::cout << header << ": " << timer.elapsed() << " (s), " << Miscellany::MemoryInfo::PeakMemoryUsageMB() << " (MB)" << std::endl;
	}

	// Computes the blade whose wedge product best fits the skew-symmetric matrix
	template< unsigned int Dim >
	std::pair< Point< double , Dim > , Point< double , Dim > > FitBlade( Hat::SkewSymmetricMatrix< double , Dim > S , unsigned int &steps , double eps );

	template< unsigned int Dim , typename Energy >
	struct Reconstructor
	{
		template< typename HermiteSampleFunctor /* = std::function< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > ( unsigned int ) >*/ >
		Reconstructor
		(
			unsigned int depth , unsigned int minSolveDepth , unsigned int maxSolveDepth ,
			HermiteSampleFunctor F , size_t sampleNum ,
			const GridSamples::Estimator< Dim > &estimator ,
			std::function< double ( double ) > noiseToWeight ,
			double dWeight , double sWeight ,
			unsigned int kernelRadius ,
			unsigned int projectionIterations ,
			double projectionEpsilon ,
			unsigned int verbosity
		);

		std::pair< Eigen::VectorXd , Eigen::VectorXd > solve
		(
			const Eigen::VectorXd initialGuess[] ,
			bool noCascadic , bool fullMultiGrid , bool singleLevel , double iterationMultiplier ,
			unsigned int vCycles ,
			unsigned int gsIterations ,
			bool separate ,
			unsigned int verbosity
		);
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
		if( !sampleNum ) ERROR_OUT( "Expected some samples: " , sampleNum );

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
	template< typename HermiteSampleFunctor /* = std::function< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > ( unsigned int ) >*/ >
	Reconstructor< Dim , Energy >::Reconstructor
	(
		unsigned int depth , unsigned int minSolveDepth , unsigned int maxSolveDepth ,
		HermiteSampleFunctor F , size_t sampleNum ,
		const GridSamples::Estimator< Dim > &estimator ,
		std::function< double ( double ) > noiseToWeight ,
		double dWeight , double sWeight ,
		unsigned int kernelRadius ,
		unsigned int projectionIterations ,
		double projectionEpsilon ,
		unsigned int verbosity
	) : _depth(depth) , _minSolveDepth(minSolveDepth) , _maxSolveDepth(maxSolveDepth)
	{
		using namespace Hat;
		static_assert( std::is_base_of< SystemEnergy< Dim > , Energy >::value , "[ERROR] Expected SystemEnergy< Dim >" );
		for( unsigned int d=0 ; d<Dim ; d++ ) _fRange.first[d] = 0 , _fRange.second[d] = (1<<_depth) + 1;

		_x.resize( _depth+1 );
		_y.resize( _depth+1 );
		_mcIndices.resize( _depth+1 );

		for( unsigned int d=0 ; d<=_depth ; d++ )
		{
			ScalarFunctions< Dim > scalars( 1<<d );

			_x[d].resize( scalars.functionNum() );
			_y[d].resize( scalars.functionNum() );

			_mcIndices[d].resize( 1<<Dim );

			for( unsigned int i=0 ; i<scalars.functionNum() ; i++ )
			{
				Index< Dim > idx = scalars.functionIndex( i );
				unsigned mcIndex = 0;
				for( unsigned int d=0 ; d<Dim ; d++ ) if( idx[d]&1 ) mcIndex |= (1<<d);
				_mcIndices[d][ mcIndex ].push_back( i );
			}
		}

		Miscellany::Timer timer;

		Hat::ScalarFunctions< Dim > scalars( 1<<_depth );
		RegularGrid< SkewSymmetricMatrix< double , Dim > , Dim > V;
		{
			unsigned int res[Dim];
			for( unsigned int d=0 ; d<Dim ; d++ ) res[d] = 1<<_depth;
			V.resize( res );
		}
		{
			auto _F = [&]( unsigned int idx )
			{
				std::pair< Point< double , Dim > , SkewSymmetricMatrix< double , Dim > > sample = F(idx);
				sample.second /= sqrt( sample.second.squareNorm() );
				// Have the magnitude of the samples contribution be proportional to its measure
				sample.second *= estimator.measure( sample.first );
				sample.second *= 1<<(_depth*CoDim);
				return sample;
			};
			std::vector< RegularGrid< SkewSymmetricMatrix< double , Dim > , Dim > > _Vs;
			{
				unsigned int threads = (unsigned int) omp_get_max_threads();
				std::vector< std::vector< GridSamples::Splatter< Dim > > > splatters( threads );
				std::vector< std::vector< RegularGrid< SkewSymmetricMatrix< double , Dim > , Dim > > > __Vs(threads);
				for( unsigned int t=0 ; t<threads ; t++ )
				{
					splatters[t].reserve( _depth+1 );
					__Vs[t].resize( _depth+1 );
					for( unsigned int depth=0 ; depth<=_depth ; depth++ )
					{
						splatters[t].emplace_back( depth , 2*kernelRadius+1 , true );
						unsigned int res[Dim];
						for( unsigned int d=0 ; d<Dim ; d++ ) res[d] = (1<<depth)+1;
						__Vs[t][depth].resize( res );
					}
				}
#pragma omp parallel for
				for( int s=0 ; s<(int)sampleNum ; s++ )
				{
					unsigned int t = (unsigned int)omp_get_thread_num();
					std::pair< Point< double , Dim > , SkewSymmetricMatrix< double , Dim > > sample = _F(s);
					double d = estimator.depth( sample.first );
					if( d>=_depth ) splatters[t][_depth]( __Vs[t][_depth]() , sample );
					else if( d<=0 ) splatters[t][_depth]( __Vs[t][0]() , sample );
					else
					{
						unsigned int d1 = (unsigned int)floor(d) , d2 = (unsigned int)floor(d)+1;
						double dd = d-d1;
						splatters[t][d1]( __Vs[t][d1]() , std::make_pair( sample.first , sample.second * (1.-dd) ) );
						splatters[t][d2]( __Vs[t][d2]() , std::make_pair( sample.first , sample.second *     dd  ) );
					}
				}

				_Vs.resize( _depth+1 );
				for( unsigned int depth=0 ; depth<=_depth ; depth++ )
				{
					unsigned int res[Dim];
					for( unsigned int d=0 ; d<Dim ; d++ ) res[d] = (1<<depth)+1;
					_Vs[depth].resize( res );
#pragma omp parallel for
					for( int i=0 ; i<(int)_Vs[depth].resolution() ; i++ ) for( unsigned int t=0 ; t<threads ; t++ ) _Vs[depth][i] += __Vs[t][depth][i];
				}
			}
			// Evaluate at the element centers
#pragma omp parallel for
			for( int i=0 ; i<(int)scalars.elementNum() ; i++ )
			{
				Index< Dim > e = scalars.elementIndex(i);
				Point< double , Dim > p;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = ( e[d] + 0.5 ) / (1<<_depth);

				SkewSymmetricMatrix< double , Dim > v;
				for( unsigned int depth=0 ; depth<=_depth ; depth++ ) v += _Vs[depth]( p * ( 1<<depth ) );
				V[i] = v;
				if( projectionIterations )
				{
					std::pair< Point< double , Dim > , Point< double , Dim > > blade = FitBlade< Dim >( V[i] , projectionIterations , projectionEpsilon );
					V[i] = Wedge( blade.first , blade.second );
				}
			}
		}
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
		for( unsigned int i=0 ; i<V.resolution() ; i++ ) V[i] *= scale;

		if( verbosity ) OutputPerformance( "Splatted" , timer );

		timer.reset();
		Eigen::SparseMatrix< double > R( scalars.functionNum() , scalars.functionNum() );

		if( sWeight>0 )
		{
			sWeight /= sampleNum;
			Eigen::SparseMatrix< double > E = scalars.evaluation( [&]( unsigned int idx ){ return F(idx).first; } , sampleNum );
			double measure = 0;
			Eigen::SparseMatrix< double > D( sampleNum , sampleNum );
			std::vector< Eigen::Triplet< double > > triplets( sampleNum );
			for( unsigned int i=0 ; i<sampleNum ; i++ ) measure += estimator.measure( F(i).first ) ,  triplets[i] = Eigen::Triplet< double >( i , i , noiseToWeight( estimator.noise( F(i).first ) ) );
			D.setFromTriplets( triplets.begin() , triplets.end() );
			R = E.transpose() * D * E * sWeight * measure;
		}
		if( verbosity ) OutputPerformance( "Set regularization energy" , timer );

		timer.reset();
		_energies.emplace_back( 1<<_depth , GetPointer( V() , V.resolution() ) , dWeight , R );
		for( unsigned int d=0 ; d<_depth ; d++ ) _energies.emplace( _energies.begin() , _energies[0].restrict() );
		if( verbosity ) OutputPerformance( "Set constraints" , timer );
		if( verbosity ) std::cout << "Target square norm: " << _energies.back().squareNorm() << std::endl;
	}

	template< unsigned int Dim , typename Energy >
	std::pair< Eigen::VectorXd , Eigen::VectorXd > Reconstructor< Dim , Energy >::solve
	(
		const Eigen::VectorXd initialGuess[] ,
		bool noCascadic , bool fullMultiGrid , bool singleLevel ,
		double iterationMultiplier ,
		unsigned int vCycles ,
		unsigned int gsIterations ,
		bool separate ,
		unsigned int verbosity
	)
	{
		using namespace Hat;

		Miscellany::Timer timer;

		auto Iters = [&]( unsigned int iters , unsigned int depth ){ return (int)ceil( pow( iterationMultiplier , _depth - depth ) * iters ); };

		auto MCIndices = [&]( unsigned int d ) -> const std::vector< std::vector< unsigned int > > & { return _mcIndices[d]; };

		bool initialGuessAtFinest = initialGuess || ( noCascadic && !fullMultiGrid ) || singleLevel;

		if( initialGuess ) _x.back() = initialGuess[0] , _y.back() = initialGuess[1];
		else
		{
			for( unsigned int i=0 ; i<_x[ _minSolveDepth ].size() ; i++ ) _x[ _minSolveDepth ][i] = Random< double >() , _y[ _minSolveDepth ][i] = Random< double >();
			double targetSquareNorm = _energies.back().squareNorm();
			// Rescale the initial guess so that it is on the order of the target skew-symmetric matrix field
			{
				Eigen::VectorXd w = _energies[ _minSolveDepth ].wedges.wedge( _x[ _minSolveDepth ] , _y[ _minSolveDepth ] );
				double scale = _energies[ _minSolveDepth ].wedgeSquareNorm(w);
				scale = sqrt( sqrt( targetSquareNorm / scale ) );
				_x[ _minSolveDepth ] *= scale , _y[ _minSolveDepth ] *= scale;
			}
			if( initialGuessAtFinest )
			{
				// Prolong from coarsest to finest and clear all the coarser levels
				for( unsigned int d=_minSolveDepth ; d<_maxSolveDepth ; d++ ) _x[d+1] = _energies[d+1].scalarProlongation( _x[d] ) , _y[d+1] = _energies[d+1].scalarProlongation( _y[d] );
			}
		}

		// Clear all levels except for the initial guess
		if( initialGuessAtFinest ) for( unsigned int d=_minSolveDepth ; d<_maxSolveDepth ; d++ ) _x[d].setZero() , _y[d].setZero();
		else                       for( unsigned int d=_minSolveDepth+1 ; d<=_maxSolveDepth ; d++ ) _x[d].setZero() , _y[d].setZero();

		double _err;
		if( initialGuess ) _err = _energies[ _maxSolveDepth ]( _x[ _maxSolveDepth ] , _y[ _maxSolveDepth ] );
		else               _err = _energies[ _maxSolveDepth ].squareNorm();

		if( singleLevel )
		{
			unsigned int d = _maxSolveDepth;
			if( separate ) ExteriorPoissonSolvers::SeparateGaussSeidelSolve< Dim >( d , MCIndices( d ) , _energies[d] , _x[d] , _y[d] , gsIterations , verbosity );
			else           ExteriorPoissonSolvers::  PairedGaussSeidelSolve< Dim >( d , MCIndices( d ) , _energies[d] , _x[d] , _y[d] , gsIterations , verbosity );
		}
		else
		{
			for( unsigned int v=0 ; v<vCycles ; v++ )
			{
				initialGuessAtFinest |= v!=0;

				Miscellany::Timer timer;
				double err = 0;
				if( verbosity>=1 ) err = _energies[ _maxSolveDepth ]( _x[ _maxSolveDepth ] , _y[ _maxSolveDepth ] );

				// Using a (half) W-cycle
				if( fullMultiGrid )
				{
					// Push the system to the coarser levels
					for( unsigned int d=_maxSolveDepth ; d+1!=_minSolveDepth ; d-- )
					{
						if( !initialGuessAtFinest )
						{
							if( d>_minSolveDepth ) _x[d].setZero() , _y[d].setZero();
						}
						else
						{
							if( d<_maxSolveDepth ) _x[d].setZero() , _y[d].setZero();
						}
						// Restrict to the next resolution
						if( d>_minSolveDepth ) _energies[d-1].update( _energies[d] , _x[d] , _y[d] );
					}

					for( unsigned int d=_minSolveDepth ; d<=_maxSolveDepth ; d++ )
					{
						// Solve to level d
						for( unsigned int dd=_minSolveDepth ; dd<=d ; dd++ )
						{
							// Solve at depth dd
							{
								if( separate ) ExteriorPoissonSolvers::SeparateGaussSeidelSolve< Dim >( dd , MCIndices( dd ) , _energies[dd] , _x[dd] , _y[dd] , Iters( gsIterations , dd ) , verbosity );
								else           ExteriorPoissonSolvers::  PairedGaussSeidelSolve< Dim >( dd , MCIndices( dd ) , _energies[dd] , _x[dd] , _y[dd] , Iters( gsIterations , dd ) , verbosity );
							}

							// Prolong to the next resolution
							if( dd<d ) _x[dd+1] += _energies[dd+1].scalarProlongation( _x[dd] ) , _y[dd+1] += _energies[d+1].scalarProlongation( _y[dd] );
						}
						// Restrict back down
						if( d<_maxSolveDepth )
						{
							// Push the system to the coarser levels
							for( unsigned int dd=d ; dd+1!=_minSolveDepth ; dd-- )
							{
								if( dd!=d ) _x[dd].setZero() , _y[dd].setZero();

								if( noCascadic )
								{
									if( separate ) ExteriorPoissonSolvers::SeparateGaussSeidelSolve< Dim >( dd , MCIndices( dd ) , _energies[dd] , _x[dd] , _y[dd] , Iters( gsIterations , dd ) , verbosity );
									else           ExteriorPoissonSolvers::  PairedGaussSeidelSolve< Dim >( dd , MCIndices( dd ) , _energies[dd] , _x[dd] , _y[dd] , Iters( gsIterations , dd ) , verbosity );
								}
								// Restrict to the next resolution
								if( dd>_minSolveDepth ) _energies[dd-1].update( _energies[dd] , _x[dd] , _y[dd] );
							}
						}
					}
				}
				else
				{
					// Fine-to-coarse
					for( unsigned int d=_maxSolveDepth ; d+1!=_minSolveDepth ; d-- )
					{
						if( !initialGuessAtFinest )
						{
							if( d>_minSolveDepth ) _x[d].setZero() , _y[d].setZero();
						}
						else
						{
							if( d<_maxSolveDepth ) _x[d].setZero() , _y[d].setZero();
						}
						// Solve at the current resolution
						if( noCascadic )
						{
							if( separate ) ExteriorPoissonSolvers::SeparateGaussSeidelSolve< Dim >( d , MCIndices( d ) , _energies[d] , _x[d] , _y[d] , Iters( gsIterations , d ) , verbosity );
							else           ExteriorPoissonSolvers::  PairedGaussSeidelSolve< Dim >( d , MCIndices( d ) , _energies[d] , _x[d] , _y[d] , Iters( gsIterations , d ) , verbosity );
						}
						// Restrict to the next resolution
						if( d>_minSolveDepth ) _energies[d-1].update( _energies[d] , _x[d] , _y[d] );
					}

					// Coarse-to-fine
					for( unsigned int d=_minSolveDepth ; d<=_maxSolveDepth ; d++ )
					{
						// Solve at the current resolution
						{
							if( separate ) ExteriorPoissonSolvers::SeparateGaussSeidelSolve< Dim >( d , MCIndices( d ) , _energies[d] , _x[d] , _y[d] , Iters( gsIterations , d ) , verbosity );
							else           ExteriorPoissonSolvers::  PairedGaussSeidelSolve< Dim >( d , MCIndices( d ) , _energies[d] , _x[d] , _y[d] , Iters( gsIterations , d ) , verbosity );
						}

						// Prolong to the next resolution
						if( d<_maxSolveDepth ) _x[d+1] += _energies[d+1].scalarProlongation( _x[d] ) , _y[d+1] += _energies[d+1].scalarProlongation( _y[d] );
					}
				}
				if( verbosity>=1 )
				{
					std::pair< double , double > _err = _energies[ _maxSolveDepth ].energies( _x[ _maxSolveDepth ] , _y[ _maxSolveDepth ] );
					{
						Miscellany::StreamFloatPrecision sfp( std::cout , 5 );
						std::cout << "Error[" << v << "]: " << err << " -> " << ( _err.first+_err.second ) << " = " << _err.first << " + " << _err.second;
					}
					{
						Miscellany::StreamFloatPrecision sfp( std::cout , 2 );
						std::cout << "\t" << std::setw(6) << timer.elapsed() << " (s)" << std::endl;
					}
				}
			}
		}

		OutputPerformance( "Solved" , timer );
		std::pair< double , double > errs = _energies[ _maxSolveDepth ].energies( _x[ _maxSolveDepth ] , _y[ _maxSolveDepth ] );
		std::cout << "Error: " << _err << " -> " << errs.first + errs.second << " = " << errs.first << " + " << errs.second << std::endl;

		// Prolong from the finest solve depth to the finest depth
		for( unsigned int d=_maxSolveDepth ; d<_depth ; d++ ) _x[d+1] = _energies[d+1].scalarProlongation( _x[d] ) , _y[d+1] = _energies[d+1].scalarProlongation( _y[d] );

		return std::make_pair( _x.back() , _y.back() );
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
#endif // EXTERIOR_POISSON_INCLUDED