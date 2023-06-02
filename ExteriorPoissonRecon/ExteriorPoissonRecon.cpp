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

#include <stdio.h>
#include <stdlib.h>
#if !defined( __clang__ )
#include <omp.h>
#endif // !__clang__
#include <vector>
#include <unordered_map>

#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Geometry.h"
#include "Misha/MarchingSimplices.h"
#include "ExteriorPoisson.h"
#include "Include/Samples.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"

static const unsigned int ProjectionIterations = 1000;
static const double ProjectionEpslion = 1e-10;
static const unsigned int CoDim = 2;

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineParameter< double > NoiseToWeightSigma( "nSigma" , 0.25 );
Misha::CmdLineParameter< unsigned int > Depth( "depth" , 5 ) , Seed( "seed" , 0 ) , TensorKernelRadius( "tRadius" , 1 ) , DensityKernelRadius( "dRadius" , 1 ) , Verbosity( "verbose" , 0 ) , MinSolveDepth( "minSolveDepth" , 0 ) , MaxSolveDepth( "maxSolveDepth" );
Misha::CmdLineParameter< double > Scale( "scale" , 1.1 ) , ScreeningWeight( "sWeight" , 50. ) , DirichletWeight( "dWeight" , 0.003125 ) , Jitter( "jitter" , 0 );
Misha::CmdLineReadable SingleLevel( "singleLevel" ) , FullMultiGrid( "fmg" ) , NoCascadic( "noCascadic" ) , ReCenter( "reCenter" ) , Project( "project" );
Misha::CmdLineParameter< unsigned int > GSIterations( "gsIters" , 10 ) , VCycles( "vCycles" , 1 );
Misha::CmdLineParameter< double > IterationMultiplier( "iMult" , 2. );
Misha::CmdLineReadable ForceHierarchical( "hierarchical" ) , Separate( "separate" );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&Depth ,
	&Seed ,
	&NoiseToWeightSigma ,
	&DensityKernelRadius ,
	&TensorKernelRadius ,
	&VCycles ,
	&GSIterations ,
	&Verbosity ,
	&ScreeningWeight ,
	&DirichletWeight ,
	&SingleLevel ,
	&Scale ,
	&IterationMultiplier ,
	&MinSolveDepth ,
	&MaxSolveDepth ,
	&NoCascadic ,
	&FullMultiGrid ,
	&Separate ,
	&ForceHierarchical ,
	&ReCenter ,
	&Project ,
	&Jitter ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input points and frames>\n" , In.name.c_str() );
	printf( "\t[--%s <output voxel grid header>]\n" , Out.name.c_str() );
	printf( "\t[--%s <grid depth>=%d]\n" , Depth.name.c_str() , Depth.value );
	printf( "\t[--%s <minimum hierarchy solver depth>=%d]\n" , MinSolveDepth.name.c_str() , MinSolveDepth.value );
	printf( "\t[--%s <maximum hierarchy solver depth>]\n" , MaxSolveDepth.name.c_str() );
	printf( "\t[--%s <noise-to-weight sigma>=%f]\n" , NoiseToWeightSigma.name.c_str() , NoiseToWeightSigma.value );
	printf( "\t[--%s <density kernel radius>=%d]\n" , DensityKernelRadius.name.c_str() , DensityKernelRadius.value );
	printf( "\t[--%s <tensor kernel radius>=%d]\n" , TensorKernelRadius.name.c_str() , TensorKernelRadius.value );
	printf( "\t[--%s <jitter magnitude>=%f]\n" , Jitter.name.c_str() , Jitter.value );
	printf( "\t[--%s <v-cycles>=%d]\n" , VCycles.name.c_str() , VCycles.value );
	printf( "\t[--%s <gs iterations>=%d]\n" , GSIterations.name.c_str() , GSIterations.value );
	printf( "\t[--%s <random seed>=%d]\n" , Seed.name.c_str() , Seed.value );
	printf( "\t[--%s <screening weight>=%g]\n" , ScreeningWeight.name.c_str() , ScreeningWeight.value );
	printf( "\t[--%s <dirichlet weight>=%g]\n" , DirichletWeight.name.c_str() , DirichletWeight.value );
	printf( "\t[--%s <bounding box scale>=%f]\n" , Scale.name.c_str() , Scale.value );
	printf( "\t[--%s <iteration multiplier>=%f]\n" , IterationMultiplier.name.c_str() , IterationMultiplier.value );
	printf( "\t[--%s]\n" , Separate.name.c_str() );
	printf( "\t[--%s]\n" , ForceHierarchical.name.c_str() );
	printf( "\t[--%s]\n" , SingleLevel.name.c_str() );
	printf( "\t[--%s]\n" , NoCascadic.name.c_str() );
	printf( "\t[--%s]\n" , FullMultiGrid.name.c_str() );
	printf( "\t[--%s]\n" , ReCenter.name.c_str() );
	printf( "\t[--%s <verbosity>=%d]\n" , Verbosity.name.c_str() , Verbosity.value );
}

template< unsigned int Dim , typename Energy >
void Execute( std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > hermiteData );

template< unsigned int Dim >
std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > ReadSamples( std::string fileName );

template< unsigned int Dim >
bool ValidSampleFile( std::string fileName );

int main( int argc , char* argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	if( !ScreeningWeight.value && !ReCenter.set )
	{
		WARN( "Forcing re-centering when there is no screening" );
		ReCenter.set = true;
	}

	bool useSingleCycleCascadic = VCycles.value==1 &&  !NoCascadic.set && !ForceHierarchical.set;

	Miscellany::Timer timer;

	unsigned int dim;
	if     ( ValidSampleFile<4>( In.value ) ) dim = 4;
	else if( ValidSampleFile<3>( In.value ) ) dim = 3;
	else ERROR_OUT( "Invalid sample file" );


	if( useSingleCycleCascadic )
	{
		switch( dim )
		{
			case 3: Execute< 3 , SingleCycleCascadicSystemEnergy< 3 > >( ReadSamples< 3 >( In.value ) ) ; break;
			case 4: Execute< 4 , SingleCycleCascadicSystemEnergy< 4 > >( ReadSamples< 4 >( In.value ) ) ; break;
			default: ERROR_OUT( "Only dimensions 3 and 4 supported" );
		}
	}
	else
	{
		switch( dim )
		{
			case 3: Execute< 3 , HierarchicalSystemEnergy< 3 > >( ReadSamples< 3 >( In.value ) ) ; break;
			case 4: Execute< 4 , HierarchicalSystemEnergy< 4 > >( ReadSamples< 4 >( In.value ) ) ; break;
			default: ERROR_OUT( "Only dimensions 3 and 4 supported" );
		}
	}

	if( Verbosity.value && Verbosity.value!=-1 ) ExteriorPoisson::OutputPerformance( "Performance" , timer );

	return EXIT_SUCCESS;
}

template< unsigned int Dim , typename Energy >
void Execute( std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > hermiteData )
{
	using namespace Hat;

	Miscellany::Timer timer;

	if( !MaxSolveDepth.set ) MaxSolveDepth.value = Depth.value;
	MaxSolveDepth.value = std::min< int >( MaxSolveDepth.value , Depth.value );
	MinSolveDepth.value = std::min< int >( MinSolveDepth.value , Depth.value );

	// Map the samples into the unit cube
	SquareMatrix< double , Dim+1 > unitCubeToWorld;
	{
		SquareMatrix< double , Dim+1 > worldToUnitCube = ExteriorPoisson::ToUnitCube< Dim >( [&]( unsigned int idx ){ return hermiteData[idx].first; } , hermiteData.size() , Scale.value );
		for( unsigned int i=0 ; i<hermiteData.size() ; i++ ) hermiteData[i].first = worldToUnitCube( hermiteData[i].first );
		unitCubeToWorld = worldToUnitCube.inverse();
	}

	if( Jitter.value>0 )
	{
		SquareMatrix< double , Dim+1 > translation = SquareMatrix< double , Dim+1 >::Identity();
		for( unsigned int d=0 ; d<Dim ; d++ ) translation( Dim , d ) = ( Random< double >() * 2. - 1. ) * Jitter.value;
		for( unsigned int i=0 ; i<hermiteData.size() ; i++ ) hermiteData[i].first = translation( hermiteData[i].first );
		unitCubeToWorld = unitCubeToWorld * translation.inverse();
	}

	// Splat the samples in
	timer.reset();
	GridSamples::MultiEstimator< Dim , CoDim > estimator = GridSamples::MultiEstimator< Dim , CoDim >::Get( DensityKernelRadius.value , Depth.value , [&]( unsigned int idx ){ return hermiteData[idx].first; } , hermiteData.size() );
	if( Verbosity.value && Verbosity.value!=-1 ) ExteriorPoisson::OutputPerformance( "Got density" , timer );
	if( Verbosity.value && Verbosity.value!=-1 ) 
	{
		unsigned int threads = (unsigned int)omp_get_max_threads();
		std::vector< double > measureSum(threads,0.) , depthSum(threads,0.) , noiseSum(threads,0.);
#pragma omp parallel for
		for( int i=0 ; i<(int)hermiteData.size() ; i++ )
		{
			unsigned int t = (unsigned int)omp_get_thread_num();
			measureSum[t] += estimator.measure( hermiteData[i].first );
			depthSum[t] += estimator.depth( hermiteData[i].first );
			noiseSum[t] += estimator.noise( hermiteData[i].first );
		}
		double measure = 0 , averageDepth = 0 , averageNoise = 0;
		for( unsigned int t=0 ; t<threads ; t++ ) measure += measureSum[t] , averageDepth += depthSum[t] , averageNoise += noiseSum[t];
		averageDepth /= hermiteData.size() , averageNoise /= hermiteData.size();
		std::cout << "Total measure / average depth / average noise: " << measure << " / " << averageDepth << " / " << averageNoise << std::endl; 
	}

	srand( Seed.value );

	auto NoiseToWeight = []( double n ){ return exp( -n*n/(2.*NoiseToWeightSigma.value*NoiseToWeightSigma.value) ); };
	timer.reset();
	ExteriorPoisson::Reconstructor< Dim , Energy > recon
	(
		Depth.value , MinSolveDepth.value , MaxSolveDepth.value ,
		[&]( unsigned int idx ){ return hermiteData[idx]; } , hermiteData.size() ,
		estimator ,
		NoiseToWeight ,
		DirichletWeight.value , ScreeningWeight.value ,
		TensorKernelRadius.value ,
		Project.set ? ProjectionIterations : 0 , ProjectionEpslion ,
		Verbosity.value
	);
	if( Verbosity.value && Verbosity.value!=-1 ) ExteriorPoisson::OutputPerformance( "Initialized reconstructor" , timer );

	std::pair< Eigen::VectorXd , Eigen::VectorXd > solution = recon.solve
	(
		NULL ,
		NoCascadic.set , FullMultiGrid.set , SingleLevel.set , IterationMultiplier.value , 
		VCycles.value , GSIterations.value , Separate.set , 
		Verbosity.value
	);

	if( ReCenter.set )
	{
		ScalarFunctions< Dim > scalars( 1<<Depth.value );
		double xErr = 0 , yErr = 0 , xAvg = 0 , yAvg = 0;
		for( unsigned int i=0 ; i<hermiteData.size() ; i++ )
		{
			Point< double , Dim > p = hermiteData[i].first;
			double xVal = scalars.value( solution.first , p ) , yVal = scalars.value( solution.second , p );
			xAvg += xVal , yAvg += yVal , xErr += xVal * xVal , yErr += yVal * yVal;
		}
		xAvg /= hermiteData.size() , yAvg /= hermiteData.size() , xErr /= hermiteData.size() , yErr /= hermiteData.size();
		if( Verbosity.value>2 )
		{
			std::cout << "Average: " << xAvg << " / " << yAvg << std::endl;
			std::cout << "Interpolation Error: " << xErr << " / " << yErr << std::endl;
		}
		for( unsigned int i=0 ; i<solution.first.size() ; i++ ) solution.first[i] -= xAvg;
		for( unsigned int i=0 ; i<solution.second.size() ; i++ ) solution.second[i] -= yAvg;
	}

	if( Out.set )
	{
		Hat::ScalarFunctions< Dim > scalars( 1<<Depth.value );

		XForm< double , Dim+1 > voxelToWorld;
		{
			SquareMatrix< double , Dim+1 > gridToCube = SquareMatrix< double , Dim+1 >::Identity();
			for( unsigned int d=0 ; d<Dim ; d++ ) gridToCube(d,d) = 1./(1<<Depth.value);
			voxelToWorld = unitCubeToWorld * gridToCube;
		}

		RegularGrid< Point< double , 2 > , Dim > solutionGrid;
		unsigned int res[Dim];
		for( unsigned int d=0 ; d<Dim ; d++ ) res[d] = (1<<Depth.value) + 1;
		solutionGrid.resize( res );
		for( unsigned int i=0 ; i<solution.first .size() ; i++ ) solutionGrid[i][0] = solution.first[i] , solutionGrid[i][1] = solution.second[i];
		solutionGrid.write( Out.value + std::string( ".grid" ) , voxelToWorld );

		{
			Hat::ScalarFunctions< Dim > scalars( 1<<Depth.value );
			RegularGrid< double , Dim > samplesPerCell;
			unsigned int res[Dim];
			for( unsigned int d=0 ; d<Dim ; d++ ) res[d] = (1<<Depth.value) + 1;
			samplesPerCell.resize( res );
			for( size_t i=0 ; i<scalars.functionNum() ; i++ )
			{
				Hat::Index< Dim > idx = scalars.functionIndex( i );
				Point< double , Dim > p;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = idx[d];
				p /= 1<<Depth.value;
#ifdef NEW_GRID_SAMPLES
				samplesPerCell[i] = estimator.samplesPerCell(p);
#else // !NEW_GRID_SAMPLES
				samplesPerCell[i] = estimator.template samplesPerCell< true >(p);
#endif // NEW_GRID_SAMPLES
				//				samplesPerCell[i] = measureAndNoiseEstimator.depth(p) / Depth.value;
			}
			samplesPerCell.write( Out.value + std::string( ".density.grid" ) , voxelToWorld );
		}
	}
}

// read from a text file
template< unsigned int Dim >
std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > ReadSamples( std::string fileName )
{
	std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > hermiteData;
	if( Misha::GetFileExtension( fileName )==std::string( "ply" ) )
	{
		using Factory = VertexFactory::Factory< double, VertexFactory::PositionFactory< double, Dim >, VertexFactory::MatrixFactory< double , Dim , Dim > >;
		Factory factory = Factory( VertexFactory::PositionFactory< double , Dim >() , VertexFactory::MatrixFactory< double , Dim , Dim >( "skew" ) );
		std::vector< typename Factory::VertexType > vertices;
		std::vector< SimplexIndex< Dim-CoDim, unsigned int > > simplexIndices;
		int fileType;
		PLY::ReadSimplices( fileName , factory , vertices , simplexIndices , NULL , fileType );
		hermiteData.resize( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ )
		{
			hermiteData[i].first = vertices[i].template get<0>();
			hermiteData[i].second = Hat::SkewSymmetricMatrix< double , Dim >( vertices[i].template get<1>() );
		}
	}
	else ERROR_OUT( "Only .ply files supported" );

	return hermiteData;
}

template< unsigned int Dim >
bool ValidSampleFile( std::string fileName )
{
	VertexFactory::PositionFactory< double , Dim > pFactory;
	int fileType;
	std::vector< PlyProperty > plyProperties = PLY::ReadVertexHeader( fileName , fileType );
	return pFactory.plyValidReadProperties( plyProperties );
}