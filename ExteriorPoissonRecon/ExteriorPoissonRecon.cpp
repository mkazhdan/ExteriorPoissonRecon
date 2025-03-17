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
#include <vector>
#include <unordered_map>

#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Geometry.h"
#include "Misha/MarchingSimplices.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include "Misha/MultiThreading.h"
#include "Misha/Streams.h"
#include "Misha/EigenStreams.h"
#include "Include/Samples.h"
#include "ExteriorPoisson.h"

using namespace MishaK;

static const unsigned int ProjectionIterations = 1000;
static const double ProjectionEpslion = 1e-10;
static const unsigned int CoDim = 2;

CmdLineParameter< std::string >
	In( "in" ) ,
	Out( "out" );

CmdLineParameter< double >
	NoiseToWeightSigma( "nSigma" , 0.25 ) ,
	Scale( "scale" , 1.1 ) ,
	ScreeningWeight( "sWeight" , 100. ) ,
	DirichletWeight( "dWeight" , 0.00625 ) ,
	IterationMultiplier( "iMult" , 2. ) ,
	Jitter( "jitter" , 0 );

CmdLineParameter< unsigned int >
	Depth( "depth" , 5 ) ,
	Seed( "seed" , 0 ) ,
	TensorKernelRadius( "tRadius" , 1 ) ,
	DensityKernelRadius( "dRadius" , 1 ) ,
	AdaptiveRadius( "aRadius" , 3 ) ,
	MinSolveDepth( "minSolveDepth" , 0 ) ,
	GSIterations( "gsIters" , 10 ) ,
	MaxSolveDepth( "maxSolveDepth" );

CmdLineReadable
	Project( "project" ) ,
	SingleLevel( "singleLevel" ) ,
	Separate( "separate" ) ,
	Verbose( "verbose" );

CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&Depth ,
	&Seed ,
	&NoiseToWeightSigma ,
	&DensityKernelRadius ,
	&TensorKernelRadius ,
	&AdaptiveRadius ,
	&GSIterations ,
	&Verbose ,
	&ScreeningWeight ,
	&DirichletWeight ,
	&Scale ,
	&IterationMultiplier ,
	&MinSolveDepth ,
	&MaxSolveDepth ,
	&Separate ,
	&SingleLevel ,
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
	printf( "\t[--%s <adaptive radius>=%d]\n" , AdaptiveRadius.name.c_str() , AdaptiveRadius.value );
	printf( "\t[--%s <jitter magnitude>=%f]\n" , Jitter.name.c_str() , Jitter.value );
	printf( "\t[--%s <gs iterations>=%d]\n" , GSIterations.name.c_str() , GSIterations.value );
	printf( "\t[--%s <random seed>=%d]\n" , Seed.name.c_str() , Seed.value );
	printf( "\t[--%s <screening weight>=%g]\n" , ScreeningWeight.name.c_str() , ScreeningWeight.value );
	printf( "\t[--%s <dirichlet weight>=%g]\n" , DirichletWeight.name.c_str() , DirichletWeight.value );
	printf( "\t[--%s <bounding box scale>=%f]\n" , Scale.name.c_str() , Scale.value );
	printf( "\t[--%s <iteration multiplier>=%f]\n" , IterationMultiplier.name.c_str() , IterationMultiplier.value );
	printf( "\t[--%s]\n" , Separate.name.c_str() );
	printf( "\t[--%s]\n" , SingleLevel.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

template< unsigned int Dim >
void Execute( std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > hermiteData );

template< unsigned int Dim , typename Energy , typename HierarchicalIndexer >
void Execute( std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > & hermiteData , ::SquareMatrix< double , Dim+1 > unitCubeToWorld , const HierarchicalIndexer & hierarchicalIndexer );

template< unsigned int Dim >
std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > ReadSamples( std::string fileName );

template< unsigned int Dim >
bool ValidSampleFile( std::string fileName );

int main( int argc , char* argv[] )
{
	Miscellany::PerformanceMeter::Width = 20;

	CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	Miscellany::PerformanceMeter pMeter;

	if     ( ValidSampleFile< 4 >( In.value ) ) Execute< 4 >( ReadSamples< 4 >( In.value ) );
	else if( ValidSampleFile< 3 >( In.value ) ) Execute< 3 >( ReadSamples< 3 >( In.value ) );
	else ERROR_OUT( "Invalid sample file" );

	if( Verbose.set ) std::cout << pMeter( "Executed" ) << std::endl;

	return EXIT_SUCCESS;
}

template< unsigned int Dim >
void Execute( std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > hermiteData )
{
	// Map the samples into the unit cube
	::SquareMatrix< double , Dim+1 > unitCubeToWorld;
	{
		::SquareMatrix< double , Dim+1 > worldToUnitCube = ExteriorPoisson::ToUnitCube< Dim >( [&]( size_t idx ){ return hermiteData[idx].first; } , hermiteData.size() , Scale.value );
		for( unsigned int i=0 ; i<hermiteData.size() ; i++ ) hermiteData[i].first = worldToUnitCube( hermiteData[i].first );
		unitCubeToWorld = worldToUnitCube.inverse();
	}

	if( Jitter.value>0 )
	{
		::SquareMatrix< double , Dim+1 > translation = ::SquareMatrix< double , Dim+1 >::Identity();
		for( unsigned int d=0 ; d<Dim ; d++ ) translation( Dim , d ) = ( Random< double >() * 2. - 1. ) * Jitter.value;
		for( unsigned int i=0 ; i<hermiteData.size() ; i++ ) hermiteData[i].first = translation( hermiteData[i].first );
		unitCubeToWorld = unitCubeToWorld * translation.inverse();
	}

	Miscellany::PerformanceMeter *pMeter = new Miscellany::PerformanceMeter();
	if( AdaptiveRadius.value==-1 )
	{
		using HierarchicalIndexer = Hat::HierarchicalRegularIndexer< Dim >;
		HierarchicalIndexer hierarchicalIndexer( Depth.value );
		if( Verbose.set ) std::cout << (*pMeter)( "Indexer" ) << std::endl;
		delete pMeter;
		return Execute< Dim , SingleCycleCascadicSystemEnergy< Dim , typename HierarchicalIndexer::Indexer > >( hermiteData , unitCubeToWorld , hierarchicalIndexer );
	}
	else
	{
		using HierarchicalIndexer = Hat::HierarchicalAdaptedIndexer< Dim >;
		HierarchicalIndexer hierarchicalIndexer( hermiteData.size() , [&]( size_t i ){ return hermiteData[i].first; } , AdaptiveRadius.value , Depth.value );
		if( Verbose.set ) std::cout << (*pMeter)( "Indexer" ) << std::endl;
		delete pMeter;
		return Execute< Dim , SingleCycleCascadicSystemEnergy< Dim , typename HierarchicalIndexer::Indexer > >( hermiteData , unitCubeToWorld , hierarchicalIndexer );
	}
}

template< unsigned int Dim , typename Energy , typename Indexer >
void Execute( std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > &  hermiteData , ::SquareMatrix< double , Dim+1 > unitCubeToWorld , const Indexer & hierarchicalIndexer )
{
	using namespace Hat;

	Miscellany::PerformanceMeter pMeter;

	if( !MaxSolveDepth.set ) MaxSolveDepth.value = Depth.value;
	MaxSolveDepth.value = std::min< int >( MaxSolveDepth.value , Depth.value );
	MinSolveDepth.value = std::min< int >( MinSolveDepth.value , Depth.value );

	// Order the samples
	OrderedSampler< Dim > orderedSampler( [&]( size_t idx ){ return hermiteData[idx].first; } , hermiteData.size() , 1<<Depth.value );
	if( Verbose.set ) std::cout << pMeter( "Ordered samples" ) << std::endl;

	// Splat the samples in
	GridSamples::TreeEstimator< Dim , 2 > estimator( DensityKernelRadius.value , Depth.value , [&]( size_t idx ){ return hermiteData[idx].first; } , orderedSampler , true );
	if( Verbose.set ) std::cout << pMeter( "Density" ) << std::endl;

	srand( Seed.value );

	auto NoiseToWeight = []( double n ){ return exp( -n*n/(2.*NoiseToWeightSigma.value*NoiseToWeightSigma.value) ); };
	pMeter.reset();
	ExteriorPoisson::Reconstructor< Dim , Energy > recon
	(
		hierarchicalIndexer ,
		Depth.value , MinSolveDepth.value , MaxSolveDepth.value ,
		[&]( size_t idx ){ return hermiteData[idx]; } , orderedSampler ,
		estimator ,
		NoiseToWeight ,
		DirichletWeight.value , ScreeningWeight.value ,
		TensorKernelRadius.value ,
		Project.set ? ProjectionIterations : 0 , ProjectionEpslion ,
		Verbose.set
	);
	if( Verbose.set ) std::cout << pMeter( "Reconstructor" ) << std::endl;

	recon.solve
	(
		SingleLevel.set ,
		IterationMultiplier.value , 
		GSIterations.value ,
		Separate.set , 
		Verbose.set
	);


	if( Out.set )
	{
		Hat::ScalarFunctions< Dim > scalars( 1<<Depth.value );

		XForm< double , Dim+1 > voxelToWorld;
		{
			::SquareMatrix< double , Dim+1 > gridToCube = ::SquareMatrix< double , Dim+1 >::Identity();
			for( unsigned int d=0 ; d<Dim ; d++ ) gridToCube(d,d) = 1./(1<<Depth.value);
			voxelToWorld = unitCubeToWorld * gridToCube;
		}

		// Output the grid of solution coefficients
		{
			FileStream fStream( Out.value + std::string( ".tree" ) , false , true );
			unsigned int dim = Dim;
			size_t hierarchicalIndexerHashCode = typeid( hierarchicalIndexer ).hash_code();
			fStream.write( dim );
			fStream.write( hierarchicalIndexerHashCode );
			fStream.write( voxelToWorld );
			fStream.write( recon.x() );
			fStream.write( recon.y() );
			hierarchicalIndexer.write( fStream , false );
		}

		// Output the density estimator
		{
			FileStream fStream( Out.value + std::string( ".density.tree" ) , false , true );
			fStream.write( voxelToWorld );
			estimator.write( fStream , false );
		}
	}
}

// read from a text file
template< unsigned int Dim >
std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > ReadSamples( std::string fileName )
{
	std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > hermiteData;
	if( GetFileExtension( fileName )==std::string( "ply" ) )
	{
		using Factory = VertexFactory::Factory< double, VertexFactory::PositionFactory< double, Dim >, VertexFactory::MatrixFactory< double , Dim , Dim > >;
		Factory factory = Factory( VertexFactory::PositionFactory< double , Dim >() , VertexFactory::MatrixFactory< double , Dim , Dim >( "skew" ) );
		std::vector< typename Factory::VertexType > vertices;
		std::vector< SimplexIndex< Dim-CoDim, unsigned int > > simplexIndices;
		int fileType = PLY::ReadSimplices( fileName , factory , vertices , simplexIndices , NULL );
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