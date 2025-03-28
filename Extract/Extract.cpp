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
#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Geometry.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include "Misha/MarchingSimplices.h"
#include "Misha/RegularGrid.h"
#include "Misha/EigenStreams.h"
#include "Include/GridSamples.h"
#include "Include/Hat.h"
#include "Include/IsoTree.h"

using namespace MishaK;

static const unsigned int CoDim = 2;

CmdLineParameter< std::string >
	In( "in" ) ,
	InDensity( "density" ) ,
	Out( "out" ) ,
	OutLevelSets( "levelSets" ) ,
	OutGridHeader( "grid" );

CmdLineParameter< double >
	TrimDensity( "trimDensity" , 0. ) ;

CmdLineReadable
	Verbose( "verbose" );

CmdLineReadable* params[] =
{
	&In ,
	&InDensity ,
	&Out ,
	&OutLevelSets ,
	&OutGridHeader ,
	&TrimDensity ,
	&Verbose ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input grid>\n" , In.name.c_str() );
	printf( "\t[--%s <density grid>]\n" , InDensity.name.c_str() );
	printf( "\t[--%s <output mesh>]\n" , Out.name.c_str() );
	printf( "\t[--%s <output level sets>]\n" , OutLevelSets.name.c_str() );
	printf( "\t[--%s <output grid header>]\n" , OutGridHeader.name.c_str() );
	printf( "\t[--%s <trim value>=%f]\n" , TrimDensity.name.c_str() , TrimDensity.value );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

Point< double , 3 > SurfaceColor( unsigned int idx )
{
	switch( idx )
	{
		case 0: return Point< double , 3 >( 223. , 41. , 53. );
		case 1: return Point< double , 3 >( 55. , 114. , 255. );
		default: MK_ERROR_OUT( "Only two surfaces supported" );
	}
	return Point< double , 3 >();
}

template< unsigned int Dim , unsigned int CoDim , typename HierarchicalIndexer >
int Execute( void )
{
	Miscellany::PerformanceMeter pMeter;

	XForm< double , Dim+1 > voxelToWorld;
	GridSamples::TreeEstimator< Dim , CoDim > density;
	IsoTree< Dim , CoDim > *isoTree = nullptr;

	{
		Miscellany::PerformanceMeter pMeter;

		unsigned int dim;
		size_t hierarchicalIndexerHashCode;
		std::vector< Eigen::VectorXd > v[CoDim];

		FileStream fStream( In.value , true , false );
		fStream.read( dim );
		fStream.read( hierarchicalIndexerHashCode );
		fStream.read( voxelToWorld );
		for( unsigned int c=0 ; c<CoDim ; c++ ) fStream.read( v[c] );
		HierarchicalIndexer hierarchicalIndexer( fStream );

		using Node = typename IsoTree< Dim , CoDim >::Node;
		auto CoefficientFunctor = [&]( unsigned int d , size_t f )
			{
				Point< double , CoDim > value;
				for( unsigned int c=0 ; c<CoDim ; c++ ) value[c] = v[c][d][f];
				return value;
			};
		isoTree = new IsoTree< Dim , CoDim >( hierarchicalIndexer , CoefficientFunctor );
		if( Verbose.set ) std::cout << pMeter( "Iso-tree" ) << std::endl;

		isoTree->refineZeroCrossing( !OutLevelSets.set );
		if( Verbose.set ) std::cout << pMeter( "Dilated" ) << std::endl;

		if( OutGridHeader.set )
		{
			Eigen::VectorXd _v[CoDim];
			for( unsigned int c=0 ; c<CoDim ; c++ ) _v[c] = hierarchicalIndexer.regularCoefficients( v[c] );
			RegularGrid< Dim , double > grid;
			grid.resize( ( 1<< hierarchicalIndexer.maxDepth() ) + 1 );
			for( unsigned int c=0 ; c<CoDim ; c++ )
			{
				for( size_t i=0 ; i<grid.resolution() ; i++ ) grid[i] = _v[c][i];
				std::string fileName = OutGridHeader.value + std::string( "." ) + std::to_string( c ) + std::string( ".grid" );
				grid.write( fileName , voxelToWorld );
			}
		}
	}

	unsigned int depth = isoTree->spaceRoot().maxDepth();

	double scale = 1.;
	for( unsigned int d=0 ; d<Dim ; d++ ) scale *= voxelToWorld(d,d);
	scale = pow( scale , 1./Dim );

	if( InDensity.set )
	{
		XForm< double , Dim+1 > xForm;
		FileStream stream( InDensity.value , true , false );
		stream.read( xForm );
		density.read( stream );
	}
	pMeter.reset();
	MarchingSimplices::SimplicialMesh< Dim , unsigned int , Point< double , Dim > > sMesh;
	{
		pMeter.reset();
		std::vector< const typename IsoTree< Dim , CoDim >::Node * > levelSetNodes = isoTree->levelSetNodes( !OutLevelSets.set );
		if( Verbose.set ) std::cout << pMeter( "Level-set nodes" ) << std::endl;

		auto LevelSetNodeIndex = [&]( size_t i )
			{
				Point< unsigned int , Dim > I;
				Hat::Index< Dim > _I = levelSetNodes[i]->offset();
				for( unsigned int d=0 ; d<Dim ; d++ ) I[d] = (unsigned int)_I[d];
				return I;
			};
		sMesh = MarchingSimplices::RegularSubGridTriangulation< Dim >( levelSetNodes.size() , LevelSetNodeIndex , true , false );
	}
	if( Verbose.set ) std::cout << pMeter( "Triangulation" ) << std::endl;

	std::vector< Point< double , CoDim > > values( sMesh.vertices.size() );
	{
		unsigned int res = 1<<depth;
		auto SimplexElement = [&]( size_t i )
			{
				Point< double , Dim > p;
				for( unsigned int d=0 ; d<=Dim ; d++ ) p += sMesh.vertices[ sMesh.simplexIndices[i][d] ];
				p /= Dim+1;
				Hat::Index< Dim > E;
				for( unsigned int d=0 ; d<Dim ; d++ ) E[d] = std::max< int >( 0 , std::min< int >( res-1 , (int)floor( p[d] ) ) );
				return E;
			};
		std::vector< Hat::Index< Dim > > vertexElements( sMesh.vertices.size() );
		for( unsigned int s=0 ; s<sMesh.simplexIndices.size() ; s++ )
		{
			Hat::Index< Dim > E = SimplexElement(s);
			for( unsigned int d=0 ; d<=Dim ; d++ ) vertexElements[ sMesh.simplexIndices[s][d] ] = E;
		}
		ThreadPool::ParallelFor
		(
			0 , sMesh.vertices.size() ,
			[&]( unsigned int t , size_t v )
			{
				Point< double , Dim > p = sMesh.vertices[v];
				Hat::Index< Dim > F;
				for( unsigned int d=0 ; d<Dim ; d++ ) F[d] = (int)floor( p[d]+0.5 );
				values[v] = isoTree->functionValue( F );
			}
		);
	}
	if( Verbose.set ) std::cout << pMeter( "Sampled" ) << std::endl;

	MarchingSimplices::SimplicialMesh< Dim-CoDim , unsigned int , Point< double , Dim > > levelSet;
	levelSet = MarchingSimplices::LevelSet< double >( sMesh , [&]( unsigned int idx ){ return values[idx]; } , Point< double , CoDim >() );
	if( Verbose.set ) std::cout << pMeter( "Level set" ) << std::endl;

	MarchingSimplices::Orient( levelSet.simplexIndices );
	if( Verbose.set ) std::cout << pMeter( "Oriented" ) << std::endl;

	if( InDensity.set && TrimDensity.value>0 )
	{
		unsigned int res = 1<<depth;
		{
			std::vector< SimplexIndex< Dim-CoDim , unsigned int > > temp;
			std::vector< std::vector< size_t > > components = MarchingSimplices::FaceAdjacentConnectedComponents( levelSet.simplexIndices );
			for( unsigned int i=0 ; i<components.size() ; i++ )
			{
				double maxDensity = 0;
				for( unsigned int j=0 ; j<components[i].size() ; j++ ) for( unsigned int v=0 ; v<=(Dim-CoDim) ; v++ )
					maxDensity = std::max< double >( maxDensity , density.samplesPerCell( levelSet.vertices[ levelSet.simplexIndices[ components[i][j] ][v] ] / res , 0 ) );
				if( Verbose.set ) std::cout << "Max Density[ " << i << " ] " << maxDensity << std::endl;

				if( maxDensity>TrimDensity.value )
				{
					temp.reserve( temp.size() + components[i].size() );
					for( unsigned int j=0 ; j<components[i].size() ; j++ ) temp.push_back( levelSet.simplexIndices[ components[i][j] ] );
				}
			}
			levelSet.simplexIndices = temp;
		}
	}

	std::vector< Point< double , Dim > > frame[2];
	{
		unsigned int res = 1<<depth;
		Hat::ScalarFunctions< Dim > scalars( res );

		frame[0].resize( levelSet.vertices.size() );
		frame[1].resize( levelSet.vertices.size() );

		Hat::HierarchicalRegularIndexer< Dim > hierarchicalIndexer( depth );
		typename Hat::HierarchicalRegularIndexer< Dim >::Indexer indexer = hierarchicalIndexer[depth];

		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ )
		{
			Point< Point< double , Dim > , CoDim , double > grad = isoTree->gradient( levelSet.vertices[i] / res , 0 , false ) / res / scale;
			frame[0][i] = grad[0];
			frame[1][i] = grad[1];
		}
	}
	delete isoTree;

	for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) levelSet.vertices[i] = voxelToWorld( levelSet.vertices[i] );

	if( Out.set )
	{
		using Factory = VertexFactory::Factory< double , VertexFactory::PositionFactory< double , Dim > , VertexFactory::MatrixFactory< double , CoDim , Dim > >;
		using Vertex = typename Factory::VertexType;
		Factory factory( VertexFactory::PositionFactory< double , Dim >() , VertexFactory::MatrixFactory< double , CoDim , Dim >( "n" ) );
		std::vector< std::vector< int > > simplices( levelSet.simplexIndices.size() );
		for( unsigned int i=0 ; i<levelSet.simplexIndices.size() ; i++ )
		{
			simplices[i].resize( Dim-CoDim+1 );
			for( unsigned int j=0 ; j<Dim-CoDim+1 ; j++ ) simplices[i][j] = (int)levelSet.simplexIndices[i][j];
		}
		std::vector< Vertex > vertices( levelSet.vertices.size() );
		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ )
		{
			vertices[i].template get<0>() = levelSet.vertices[i];
			for( unsigned int c=0 ; c<CoDim ; c++ ) for( unsigned int d=0 ; d<Dim ; d++ ) vertices[i].template get<1>()(c,d) = frame[c][i][d];
		}
		PLY::WritePolygons< Factory , int >( Out.value , factory , vertices , simplices , PLY_BINARY_NATIVE );
	}

	// Output the level sets
	if( OutLevelSets.set )
	{
		pMeter.reset();

		using Factory = VertexFactory::Factory< double , VertexFactory::PositionFactory< double , Dim > , VertexFactory::RGBColorFactory< double > >;
		using Vertex = typename Factory::VertexType;

		std::vector< Vertex > vertices;
		std::vector< SimplexIndex< Dim-1 , int > > simplices;

		for( unsigned int c=0 ; c<CoDim ; c++ )
		{
			pMeter.reset();
			MarchingSimplices::SimplicialMesh< Dim-1 , unsigned int , Point< double , Dim > > levelSet;
			levelSet = MarchingSimplices::LevelSet< double >( sMesh , [&]( unsigned int idx ){ return values[idx][c]; } , Point< double , 1 >() );
			for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) levelSet.vertices[i] = voxelToWorld( levelSet.vertices[i] );
			if( Verbose.set ) std::cout << pMeter( "Level-set" ) << std::endl;

			auto AddLevelSet = []( const MarchingSimplices::SimplicialMesh< Dim-1 , unsigned int , Point< double , Dim > > &levelSet , std::vector< Vertex > &vertices , std::vector< SimplexIndex< Dim-1 , int > > &simplices , Point< double , 3 > color )
				{
					size_t sz = vertices.size();
					vertices.reserve( sz + levelSet.vertices.size() );
					simplices.reserve( simplices.size() + levelSet.simplexIndices.size() );

					Vertex v;
					v.template get<1>() = color;

					for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ )
					{
						v.template get<0>() = levelSet.vertices[i];
						vertices.push_back( v );
					}

					for( unsigned int i=0 ; i<levelSet.simplexIndices.size() ; i++ )
					{
						SimplexIndex< Dim-1 , int > si;
						for( unsigned int j=0 ; j<=Dim-1 ; j++ ) si[j] = (int)( levelSet.simplexIndices[i][j] + (unsigned int)sz );
						simplices.push_back( si );
					}
				};

			AddLevelSet( levelSet , vertices , simplices , SurfaceColor(c) );
		}

		Factory factory;
		PLY::WriteSimplices( OutLevelSets.value , factory , vertices , simplices , PLY_BINARY_NATIVE );
	}
	return EXIT_SUCCESS;
}

int main( int argc , char* argv[] )
{
	CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	unsigned int dim;
	size_t hierarchicalIndexerHashCode;
	{
		FileStream fStream( In.value , true , false );
		fStream.read( dim );
		fStream.read( hierarchicalIndexerHashCode );
	}
	
	if( dim==3 )
	{
		if     ( hierarchicalIndexerHashCode==typeid( Hat::HierarchicalRegularIndexer< 3 > ).hash_code() ) return Execute< 3 , CoDim , Hat::HierarchicalRegularIndexer< 3 > >();
		else if( hierarchicalIndexerHashCode==typeid( Hat::HierarchicalAdaptedIndexer< 3 > ).hash_code() ) return Execute< 3 , CoDim , Hat::HierarchicalAdaptedIndexer< 3 > >();
		else MK_ERROR_OUT( "Unrecognized hierarchical indexer hash: "  , hierarchicalIndexerHashCode );
	}
	else if( dim==4 )
	{
		if     ( hierarchicalIndexerHashCode==typeid( Hat::HierarchicalRegularIndexer< 4 > ).hash_code() ) return Execute< 4 , CoDim , Hat::HierarchicalRegularIndexer< 4 > >();
		else if( hierarchicalIndexerHashCode==typeid( Hat::HierarchicalAdaptedIndexer< 4 > ).hash_code() ) return Execute< 4 , CoDim , Hat::HierarchicalAdaptedIndexer< 4 > >();
		else MK_ERROR_OUT( "Unrecognized hierarchical indexer hash: "  , hierarchicalIndexerHashCode );
	}
	else MK_ERROR_OUT( "Only dimensions 3 and 4 supported: " , dim );
	return EXIT_FAILURE;

	return EXIT_SUCCESS;
}
