#ifndef GRID_SAMPLES_INCLUDED

#include <vector>
#include <Eigen/Eigenvalues>
#include "Misha/Miscellany.h"
#include "Misha/Geometry.h"
#include "Misha/RegularGrid.h"
#include "Misha/Algebra.h"
#include "Hat.h"

//#define NEW_GRID_SAMPLES	// For multi-sampling, samplesPerCell is defined by depth

namespace GridSamples
{
	template< unsigned int Degree=(unsigned int)-1 >struct BSpline;

	// Run-time implementation of B-spline weight evaluation
	template<>
	struct BSpline< (unsigned int)-1 >
	{
		// Assuming x \in [0,1], evaluates the supported B-splines
		static void Evaluate( unsigned int degree , double x , double values[/*degree+1*/] );
	};

	// Templated implementation of B-spline weight evaluation
	template< unsigned int Degree >
	struct BSpline
	{
		// Assuming x \in [0,1], evaluates the supported B-splines
		static void Evaluate( double x , double values[Degree+1] );
	};


	// A pure abstract class that supports evaluating the measure, noise, and depth at a point in the unit cube.
	template< unsigned int Dim > struct Estimator
	{
		// The measure per unit volume at a given position
		virtual double measure( Point< double , Dim > ) const = 0;
		// The estimate of noise at a given position
		virtual double noise( Point< double , Dim > ) const = 0;
		// The measur eof the sample's depth at a given position
		virtual double depth( Point< double , Dim > ) const = 0;
	};

	// A struct describing the distribution of samples
	template< unsigned int Dim >
	struct DensityAndNoiseInfo : public VectorSpace< double , DensityAndNoiseInfo< Dim > >
	{
		double weightSum;
		Point< double , Dim > positionSum;
		SquareMatrix< double , Dim > covarianceSum;
		DensityAndNoiseInfo( void );
		DensityAndNoiseInfo( Point< double , Dim > p , double weight );
		//////////////////////////
		// Vector space methods //
		void Add  ( const DensityAndNoiseInfo &ni );
		void Scale( double s );
		//////////////////////////

		double noise( unsigned int coDim ) const;
		double measure( void ) const;
	};

	// A struct for splatting a sample into a grid
	template< unsigned int Dim >
	struct Splatter
	{
		Splatter( unsigned int depth , unsigned int bSplineDegree , bool unitKernel );
		~Splatter( void );

		template< typename DataType >
		void operator()( DataType *g , std::pair< Point< double , Dim > , DataType > sample );

		// Splats the values into the corners of the cell containg the sample using the B-spline weights
		template< typename DataType , typename SampleFunctor /* = std::function< std::pair< Point< double , Dim > , DataType > ( unsigned int idx ) */ >
		static RegularGrid< DataType , Dim > SplatSamples( unsigned int bSplineDegree , unsigned int depth , SampleFunctor F , size_t sampleNum , DataType zeroValue , bool unitKernel );

		// Splats the values into the corners of the cell containg the sample using multi-linear weights
		template< typename DataType , typename SampleFunctor /* = std::function< std::pair< Point< double , Dim > , DataType > ( unsigned int idx ) */ >
		static void SplatSamples( unsigned int bSplineDegree , RegularGrid< DataType , Dim > &g , SampleFunctor F , size_t sampleNum , bool unitKernel );

		// Splats the values into the corners of the cell containg the sample using multi-linear weights
		template< typename DataType , typename SampleFunctor /* = std::function< std::pair< Point< double , Dim > , DataType > ( unsigned int idx ) */ >
		static void SplatSamples_parallel( unsigned int bSplineDegree , RegularGrid< DataType , Dim > &g , SampleFunctor F , size_t sampleNum , DataType zeroValue , bool unitKernel );

	protected:
		double *_values;
		Hat::ScalarFunctions< Dim > _scalars;
		int _bSplineStart;
		unsigned int _bSplineDegree , _res;
		Hat::Range< Dim > _range;
		double _kernelScale;
	};

	// A struct for splatting a sample from a (Dim-CoDim)-dimensional manifold, living in Dim-dimensional Euclidean space into a grid
	template< unsigned int Dim , unsigned int CoDim >
	struct SingleEstimator : public RegularGrid< DensityAndNoiseInfo< Dim > , Dim > , Estimator< Dim >
	{
		SingleEstimator( unsigned int depth );
		unsigned int depth( void ) const;
		DensityAndNoiseInfo< Dim > operator()( Point< double , Dim > p ) const;
		double        measure( Point< double , Dim > p ) const;
		double          noise( Point< double , Dim > p ) const;
		double samplesPerCell( Point< double , Dim > p ) const;
		double          depth( Point< double , Dim > p ) const;

		template< typename SampleFunctor /* = std::function< Point< double , Dim > > ( unsigned int idx ) */ >
		static SingleEstimator Get( unsigned int kernelRadius , unsigned int depth , SampleFunctor F , size_t sampleNum );

		template< typename SampleFunctor /* = std::function< Point< double , Dim > > ( unsigned int idx ) */ >
		static void Set( SingleEstimator &estimator , unsigned int kernelRadius , SampleFunctor F , size_t sampleNum );

	protected:
		unsigned int _depth;
	};

	template< unsigned int Dim , unsigned int CoDim >
	struct MultiEstimator : public Estimator< Dim >
	{
		MultiEstimator( unsigned int depth , double samplesPerCell );
		double depth( Point< double , Dim > p ) const;
		DensityAndNoiseInfo< Dim > operator()( Point< double , Dim > p ) const;
		double        measure( Point< double , Dim > p ) const;
		double          noise( Point< double , Dim > p ) const;
#ifdef NEW_GRID_SAMPLES
#else // !NEW_GRID_SAMPLES
		template< bool Finest >
#endif // NEW_GRID_SAMPLES
		double samplesPerCell( Point< double , Dim > p ) const;

		template< typename SampleFunctor /* = std::function< Point< double , Dim > > ( unsigned int idx ) */ >
		static MultiEstimator Get( unsigned int kernelRadius , unsigned int depth , SampleFunctor F , size_t sampleNum , double samplesPerCell=pow(2.,Dim-CoDim) );
	protected:
		std::vector< SingleEstimator< Dim , CoDim > > _info;
		double _samplesPerCell;
	};

	/////////////
	// BSpline //
	/////////////
	template< unsigned int Degree >
	void BSpline< Degree >::Evaluate( double x , double values[Degree+1] )
	{
		if constexpr( Degree==0 ) values[0] = 1.;
		else
		{
			const double Scale = 1./Degree;
			BSpline< Degree-1 >::Evaluate( x , values+1 );
			values[0] = values[1] * (1.-x) * Scale;
			for( unsigned int i=1 ; i<Degree ; i++ )
			{
				double x1 = (x-i+Degree) , x2 = (-x+i+1);
				values[i] = ( values[i]*x1 + values[i+1]*x2 ) * Scale;
			}
			values[Degree] *= x * Scale;
		}
	}

	void BSpline< (unsigned int)-1 >::Evaluate( unsigned int degree , double x , double values[/*degree+1*/] )
	{
		if( degree==0 ) values[0] = 1.;
		else
		{
			const double Scale = 1./degree;
			Evaluate( degree-1 , x , values+1 );
			values[0] = values[1] * (1.-x) * Scale;
			for( unsigned int i=1 ; i<degree ; i++ )
			{
				double x1 = (x-i+degree) , x2 = (-x+i+1);
				values[i] = ( values[i]*x1 + values[i+1]*x2 ) * Scale;
			}
			values[degree] *= x * Scale;
		}
	}

	/////////////////////////
	// DensityAndNoiseInfo //
	/////////////////////////
	template< unsigned int Dim >
	DensityAndNoiseInfo< Dim >::DensityAndNoiseInfo( void ) : weightSum(0) {}

	template< unsigned int Dim >
	DensityAndNoiseInfo< Dim >::DensityAndNoiseInfo( Point< double , Dim > p , double weight ) : weightSum(weight) , positionSum(p*weight){ for( unsigned int i=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) covarianceSum(i,j) = p[i] * p[j] * weight; }

	template< unsigned int Dim >
	void DensityAndNoiseInfo< Dim >::Add  ( const DensityAndNoiseInfo &ni ){ weightSum += ni.weightSum , positionSum += ni.positionSum , covarianceSum += ni.covarianceSum; }
	
	template< unsigned int Dim >
	void DensityAndNoiseInfo< Dim >::Scale( double s ){ weightSum *= s , positionSum *= s , covarianceSum *= s; }

	template< unsigned int Dim >
	double DensityAndNoiseInfo< Dim >::noise( unsigned int coDim ) const
	{
		// Recall that if \sum_i w_i = 1 and c = \sum_i w_i * p_i
		//		C = \sum_i w_i * ( p_i - c ) * ( p_i - c )^t
		//		  = \sum_i w_i * ( p_i * p_i^t - p_i * c^t - c * p_i^t + c * c^t )
		//		  = \sum_i w_i * p_i * p_i^t - \sum_i w_i( p_i * c^t + c * p_i^t ) + ( \sum_i w_i ) * c * c^t
		//		  = \sum_i w_i * p_i * p_i^t - 2 * c * c^t + c * c^t
		//		  = \sum_i w_i * p_i * p_i^t - c * c^t
		Eigen::Matrix< double , Dim , Dim > C;

		for( unsigned int i=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) C(i,j) = covarianceSum(i,j)/weightSum - positionSum[i]/weightSum * positionSum[j]/weightSum;
		Eigen::SelfAdjointEigenSolver< Eigen::Matrix< double , Dim , Dim > > es( C );
		double sum = 0 , lowSum = 0;
		for( unsigned int d=0 ; d<Dim ; d++ )
		{
			sum += es.eigenvalues()[d];
			if( d<coDim ) lowSum += es.eigenvalues()[d];
		}
		return ( lowSum * Dim ) /( sum * coDim );
	}

	template< unsigned int Dim >
	double DensityAndNoiseInfo< Dim >::measure( void ) const{ return 1./weightSum; }

	//////////////
	// Splatter //
	//////////////
	template< unsigned int Dim >
	Splatter< Dim >::Splatter( unsigned int depth , unsigned int bSplineDegree , bool unitKernel ) : _res(1<<depth) , _scalars(1<<depth) , _bSplineDegree(bSplineDegree)
	{
		_values = new double[ Dim * ( _bSplineDegree+1 ) ];
		_bSplineStart = -(int)(_bSplineDegree>>1);
		for( unsigned int d=0 ; d<Dim ; d++ )
		{
			_range.first [d] = _bSplineStart;
			_range.second[d] = _bSplineStart+_bSplineDegree+1;
		}
		_kernelScale = unitKernel ? ( 1<<(depth*Dim) ) : 1.;
	}
	template< unsigned int Dim >
	Splatter< Dim >::~Splatter( void ){ delete _values; }

	template< unsigned int Dim >
	template< typename DataType >
	void Splatter< Dim >::operator()( DataType *g , std::pair< Point< double , Dim > , DataType > sample )
	{
		Hat::Index< Dim > e;
		sample.first *= _res;
		for( unsigned int d=0 ; d<Dim ; d++ )
		{
			e[d] = (int)floor( sample.first[d] );
			BSpline<>::Evaluate( _bSplineDegree , sample.first[d]-e[d] , _values+d*(_bSplineDegree+1) );
		}

		auto f = [&]( Hat::Index< Dim > _f )
		{
			Hat::Index< Dim > f = _f + e;
			double scale = _kernelScale;
			bool inRange = true;
			for( unsigned int d=0 ; d<Dim ; d++ )
				if( f[d]>=0 && f[d]<=(int)_res ) scale *= _values[ d * ( _bSplineDegree+1) + _f[d] - _bSplineStart ];
				else inRange = false;
			if( inRange ) g[ _scalars.functionIndex( f ) ] += scale * sample.second;
		};
		_range.process(f);
	}

	template< unsigned int Dim >
	template< typename DataType , typename SampleFunctor /* = std::function< std::pair< Point< double , Dim > , DataType > ( unsigned int idx ) */ >
	void Splatter< Dim >::SplatSamples( unsigned int bSplineDegree , RegularGrid< DataType , Dim > &g , SampleFunctor F , size_t sampleNum , bool unitKernel )
	{
		unsigned int res = g.res(0)-1;
		for( unsigned int d=1 ; d<Dim ; d++ ) if( g.res(d)-1!=res ) ERROR_OUT( "Not a cubical grid" );
		unsigned int depth = 0;
		while( (1u<<depth)<res ) depth++;
		if( res!=(1<<depth) ) ERROR_OUT( "Resolution is not a power of two: " , res );

		Splatter splatter( depth , bSplineDegree , unitKernel );
		for( unsigned int s=0 ; s<sampleNum ; s++ ) splatter( g() , F(s) );
	}

	template< unsigned int Dim >
	template< typename DataType , typename SampleFunctor /* = std::function< std::pair< Point< double , Dim > , DataType > ( unsigned int idx ) */ >
	void Splatter< Dim >::SplatSamples_parallel( unsigned int bSplineDegree , RegularGrid< DataType , Dim > &g , SampleFunctor F , size_t sampleNum , DataType zeroValue , bool unitKernel )
	{
		unsigned int res = g.res(0)-1;
		for( unsigned int d=1 ; d<Dim ; d++ ) if( g.res(d)-1!=res ) ERROR_OUT( "Not a cubical grid" );
		unsigned int depth = 0;
		while( (1u<<depth)<res ) depth++;
		if( res!=(1<<depth) ) ERROR_OUT( "Resolution is not a power of two: " , res );

		std::vector< RegularGrid< DataType , Dim > > gs( omp_get_max_threads() );
		std::vector< Splatter > splatters;
		splatters.reserve( omp_get_max_threads() );
		for( int t=0 ; t<omp_get_max_threads() ; t++ )
		{
			gs[t].resize( g.res() );
			splatters.emplace_back( depth , bSplineDegree , unitKernel );
		}

#pragma omp parallel for
		for( int s=0 ; s<(int)sampleNum ; s++ )
		{
			unsigned int t = omp_get_thread_num();
			splatters[t]( gs[t]() , F(s) );
		}
#pragma omp parallel for
		for( int i=0 ; i<(int)g.resolution() ; i++ ) for( unsigned int t=0 ; t<gs.size() ; t++ ) g[i] += gs[t][i];
	}

	template< unsigned int Dim >
	template< typename DataType , typename SampleFunctor /* = std::function< std::pair< Point< double , Dim > , DataType > ( unsigned int idx ) */ >
	RegularGrid< DataType , Dim > Splatter< Dim >::SplatSamples( unsigned int bSplineDegree , unsigned int depth , SampleFunctor F , size_t sampleNum , DataType zeroValue , bool unitKernel )
	{
		RegularGrid< DataType , Dim > g;
		{
			unsigned int res[Dim];
			for( unsigned int d=0 ; d<Dim ; d++ ) res[d] = (1<<depth) + 1;
			g.resize( res );
			for( size_t i=0 ; i<g.resolution() ; i++ ) g[i] = zeroValue;
		}
		SplatSamples( bSplineDegree , g , F , sampleNum , unitKernel );
		return g;
	}

	/////////////////////
	// SingleEstimator //
	/////////////////////
	template< unsigned int Dim , unsigned int CoDim >
	SingleEstimator< Dim , CoDim >::SingleEstimator( unsigned int depth ) : _depth(depth)
	{
		unsigned int res[Dim];
		for( unsigned int dd=0 ; dd<Dim ; dd++ ) res[dd] = (1<<depth) + 1;
		RegularGrid< DensityAndNoiseInfo< Dim > , Dim >::resize( res );
	}

	template< unsigned int Dim , unsigned int CoDim >
	unsigned int SingleEstimator< Dim , CoDim >::depth( void ) const { return _depth; }

	template< unsigned int Dim , unsigned int CoDim >
	DensityAndNoiseInfo< Dim > SingleEstimator< Dim , CoDim >::operator()( Point< double , Dim > p ) const { return RegularGrid< DensityAndNoiseInfo< Dim >  , Dim >::operator()( p*(1<<_depth) ); }

	template< unsigned int Dim , unsigned int CoDim >
	double SingleEstimator< Dim , CoDim >::measure( Point< double , Dim > p ) const { return this->operator()( p ).measure(); }

	template< unsigned int Dim , unsigned int CoDim >
	double SingleEstimator< Dim , CoDim >::noise( Point< double , Dim > p ) const { return this->operator()( p ).noise( CoDim ); }

	template< unsigned int Dim , unsigned int CoDim >
	double SingleEstimator< Dim , CoDim >::depth( Point< double , Dim > p ) const { return _depth; }

	template< unsigned int Dim , unsigned int CoDim >
	double SingleEstimator< Dim , CoDim >::samplesPerCell( Point< double , Dim > p ) const { return this->operator()( p ).weightSum / ( 1<<(_depth*(Dim-CoDim) ) ); }

	template< unsigned int Dim , unsigned int CoDim >
	template< typename SampleFunctor /* = std::function< Point< double , Dim > > ( unsigned int idx ) */ >
	SingleEstimator< Dim , CoDim > SingleEstimator< Dim , CoDim >::Get( unsigned int kernelRadius , unsigned int depth , SampleFunctor F , size_t sampleNum )
	{
		SingleEstimator estimator( depth );
		Set( estimator , kernelRadius , F , sampleNum );
		return estimator;
	}

	template< unsigned int Dim , unsigned int CoDim >
	template< typename SampleFunctor /* = std::function< Point< double , Dim > > ( unsigned int idx ) */ >
	void SingleEstimator< Dim , CoDim >::Set( SingleEstimator< Dim , CoDim > &estimator , unsigned int kernelRadius , SampleFunctor F , size_t sampleNum )
	{
		double scale;
		{
			unsigned int bSplineDegree = 2*kernelRadius + 1;
			double *_values = new double[bSplineDegree+1];
			BSpline<>::Evaluate( bSplineDegree , 0.5 , _values );
			scale = pow( 1./_values[kernelRadius+1] , CoDim );
			delete[] _values;
		}
		auto _F = [&]( unsigned int idx ){ return std::make_pair( F(idx) , DensityAndNoiseInfo< Dim >( F(idx) , scale / pow( pow( 0.5 , estimator.depth() ) , Dim - CoDim ) ) ); };
		Splatter< Dim >::template SplatSamples_parallel< DensityAndNoiseInfo< Dim > >( 2*kernelRadius+1 , estimator , _F , sampleNum , DensityAndNoiseInfo< Dim >() , false );
	}

	////////////////////
	// MultiEstimator //
	////////////////////
	template< unsigned int Dim , unsigned int CoDim >
	MultiEstimator< Dim , CoDim >::MultiEstimator( unsigned int depth , double samplesPerCell ) : _samplesPerCell(samplesPerCell)
	{
		_info.reserve( depth+1 );
		for( unsigned int d=0 ; d<=depth ; d++ ) _info.emplace_back( d );
	}

	template< unsigned int Dim , unsigned int CoDim >
	double MultiEstimator< Dim , CoDim >::depth( Point< double , Dim > p ) const
	{
		// Find the finest depth at which the number of samples per node exceeds _samplesPerCell
		unsigned int d1 = (unsigned int)_info.size()-1 , d2 = (unsigned int)_info.size();
		while( d2 && _info[d1].samplesPerCell(p)<_samplesPerCell ) d1-- , d2--;
		if( d2==0 || d2==_info.size() ) // Extremal case
		{
			unsigned int d;
			if     ( d2==0 ) d=0;
			else if( d2==_info.size() ) d = (unsigned int)_info.size()-1;
			double samplesPerCell = _info[d].samplesPerCell(p);
			// Generically the assumption is that increase the depth by should decrease the number of samples per cell by a factor of 2^{Dim-CoDim}
			// Solve for the depth such that:
			//		_samplesPerCell = samplesPerCell / 2^{depth-d}
			//		2^{depth-d} = samplesPerCell/_samplesPerCell
			//		depth - d = log( samplesPerCell/_samplesPerCell ) / log(2)
			//		depth = log( samplesPerCell/_samplesPerCell ) / log(2) + d
			return d + log( samplesPerCell/_samplesPerCell ) / log(2.);
		}
		else // Interior case
		{
			// Here we assume that the change in the samples per cell is simply samplesPerCell[d2]/samplesPerCell[d1]
			// Solve for the depth such that:
			//		_samplesPerCell = samplesPerCell[d1] * pow( samplesPerCell[d2]/samplesPerCell[d1] , d-d1 )
			//		pow( samplesPerCell[d2]/samplesPerCell[d1] , d-d1 ) = _samplesPerCell / samplesPerCell[d1]
			//		d-d1 = log( _samplesPerCell / samplesPerCell[d1] ) / log( samplesPerCell[d2]/samplesPerCell[d1] )
			//		d = log( _samplesPerCell / samplesPerCell[d1] ) / log( samplesPerCell[d2]/samplesPerCell[d1] ) + d1
			double samplesPerCell1 = _info[d1].samplesPerCell(p) , samplesPerCell2 = _info[d2].samplesPerCell(p);
			if( !samplesPerCell2 ) return d1;
			else return d1 + log( _samplesPerCell/samplesPerCell1 ) / log( samplesPerCell2/samplesPerCell1 );
		}
	}

	template< unsigned int Dim , unsigned int CoDim >
	DensityAndNoiseInfo< Dim > MultiEstimator< Dim , CoDim >::operator()( Point< double , Dim > p ) const
	{
		double depth = this->depth( p );
		if( depth>=_info.size()-1 ) return _info.back()( p );
		else if( depth<0 ) return _info[0](p);
		else
		{
			int d1 = (int)floor(depth) , d2 = (int)floor(depth)+1;
			double d = depth-d1;
			return _info[d1](p) * (1.-d) + _info[d2](p) * d;
		}
	}

	template< unsigned int Dim , unsigned int CoDim >
	double MultiEstimator< Dim , CoDim >::       measure( Point< double , Dim > p ) const { return this->operator()( p ).measure(); }

	template< unsigned int Dim , unsigned int CoDim >
	double MultiEstimator< Dim , CoDim >::         noise( Point< double , Dim > p ) const { return this->operator()( p ).noise( CoDim ); }

#ifdef NEW_GRID_SAMPLES
	template< unsigned int Dim , unsigned int CoDim >
	double MultiEstimator< Dim , CoDim >::samplesPerCell( Point< double , Dim > p ) const
	{
		// If at depth d, the cell has _samplesPerCell samples per cell then,
		// assuming that the number samples per cell decreases by a factor of 2^{Dim-CoDim} per level,
		// the number of samples per cell at depth _depth is _samplesPerCell * pow( 0.5 , (Dim-CoDim) * (_depth-d)
		return _samplesPerCell * pow( 0.5 , (Dim-CoDim) * ( (int)_info.size()-1 - depth(p) ) );
	}
#else // !NEW_GRID_SAMPLES
	template< unsigned int Dim , unsigned int CoDim >
	template< bool Finest >
	double MultiEstimator< Dim , CoDim >::samplesPerCell( Point< double , Dim > p ) const
	{
		if constexpr( Finest ) return _info.back().samplesPerCell(p);
		else return this->operator()( p ).weightSum / ( 1<<((_info.size()-1)*(Dim-CoDim) ) );
	}
#endif // NEW_GRID_SAMPLES

	template< unsigned int Dim , unsigned int CoDim >
	template< typename SampleFunctor /* = std::function< Point< double , Dim > > ( unsigned int idx ) */ >
	MultiEstimator< Dim , CoDim > MultiEstimator< Dim , CoDim >::Get( unsigned int kernelRadius , unsigned int depth , SampleFunctor F , size_t sampleNum , double samplesPerCell )
	{
		MultiEstimator estimator( depth , samplesPerCell );
		for( unsigned int d=0 ; d<=depth ; d++ ) SingleEstimator< Dim , CoDim >::Set( estimator._info[d] , kernelRadius , F , sampleNum );
		return estimator;
	}
}

#endif // GRID_SAMPLES_INCLUDED