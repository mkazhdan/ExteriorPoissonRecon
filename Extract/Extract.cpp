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
#include "Include/Hat.h"

#define OUTPUT_GRADIENTS
static const unsigned int CoDim = 2;

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" ) , LevelSets( "levelSets" ) , Density( "density" );
Misha::CmdLineParameter< double > TrimDensity( "trimDensity" , 0. ) ;
Misha::CmdLineReadable Verbose( "verbose" );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&LevelSets ,
	&Density ,
	&TrimDensity ,
	&Verbose ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input grid>\n" , In.name.c_str() );
	printf( "\t[--%s <density grid>]\n" , Density.name.c_str() );
	printf( "\t[--%s <output mesh>]\n" , Out.name.c_str() );
	printf( "\t[--%s <output level sets>]\n" , LevelSets.name.c_str() );
	printf( "\t[--%s <trim value>=%f]\n" , TrimDensity.name.c_str() , TrimDensity.value );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

#ifdef OUTPUT_GRADIENTS
const std::string GradientNames1[] = { "nx1" , "ny1" , "nz1" , "nw1" };
const std::string GradientNames2[] = { "nx2" , "ny2" , "nz2" , "nw2" };
#endif // OUTPUT_GRADIENTS

template< typename Data , unsigned int Dim >
unsigned int GridDepth( const RegularGrid< Data , Dim > &grid )
{
	unsigned int depth = 0;
	for( unsigned int d=0 ; d<Dim ; d++ )
	{
		unsigned int _depth = 0;
		while( ( (1<<_depth)+1u )<grid.res(d) ) _depth++;
		if( ( (1<<_depth)+1 )!=grid.res(d) ) ERROR_OUT( "Resolution not a power of two: " , _depth );
		if( d==0 ) depth = _depth;
		else if( depth!=_depth ) ERROR_OUT( "Grid resolutions vary" );
	}
	return depth;
}

template< unsigned int Dim >
void Execute( void )
{
	Miscellany::Timer timer;

	RegularGrid< Point< double , CoDim > , Dim > grid;
	RegularGrid< double , Dim > density;
	XForm< double , Dim+1 > voxelToWorld;

	grid.read( In.value , voxelToWorld );
	unsigned int depth = GridDepth( grid );
	double scale = 1.;
	for( unsigned int d=0 ; d<Dim ; d++ ) scale *= voxelToWorld(d,d);
	scale = pow( scale , 1./Dim );

	double densityScale = 1. / ( 1<<depth );
	if( Density.set )
	{
		XForm< double , Dim+1 > xForm;
		density.read( Density.value , xForm );
		if( depth!=GridDepth( density ) ) ERROR_OUT( "Density resolution does not match: " , depth , " != " , GridDepth( density ) );
		densityScale *= 1<<GridDepth( density );
	}

	timer.reset();
	MarchingSimplices::SimplicialMesh< Dim , unsigned int , Point< double , Dim > > sMesh = MarchingSimplices::RegularGridTriangulation< Dim >( 1<<depth );
	if( Verbose.set ) std::cout << "Got ambient space triangulation: " << timer.elapsed() << " (s)" << std::endl;

	timer.reset();
	std::vector< Point< double , CoDim > > values( sMesh.vertices.size() );
	for( unsigned int i=0 ; i<sMesh.vertices.size() ; i++ ) values[i] = grid( sMesh.vertices[i] );
	if( Verbose.set ) std::cout << "Sampled implicit function: " << timer.elapsed() << " (s)" << std::endl;

	timer.reset();
	MarchingSimplices::SimplicialMesh< Dim-CoDim , unsigned int , Point< double , Dim > > levelSet;
	levelSet = MarchingSimplices::LevelSet< double >( sMesh , [&]( unsigned int idx ){ return values[idx]; } , Point< double , CoDim >() );
	if( Verbose.set ) std::cout << "Level set extracted: " << timer.elapsed() << " (s)" << std::endl;

	{
		Miscellany::Timer timer;
		MarchingSimplices::Orient( levelSet.simplexIndices );
		std::cout << "Oriented: " << timer.elapsed() << " (s)" << std::endl;
	}

	if( Density.set && TrimDensity.value>0 )
	{
		std::vector< SimplexIndex< Dim-CoDim , unsigned int > > temp;
		std::vector< std::vector< size_t > > components = MarchingSimplices::FaceAdjacentConnectedComponents( levelSet.simplexIndices );
		for( unsigned int i=0 ; i<components.size() ; i++ )
		{
			double maxDensity = 0;
			for( unsigned int j=0 ; j<components[i].size() ; j++ ) for( unsigned int v=0 ; v<=(Dim-CoDim) ; v++ )
				maxDensity = std::max< double >( maxDensity , density( levelSet.vertices[ levelSet.simplexIndices[ components[i][j] ][v] ]*densityScale ) );
			if( Verbose.set ) std::cout << "Max Density[ " << i << " ] " << maxDensity << std::endl;

			if( maxDensity>TrimDensity.value )
			{
				temp.reserve( temp.size() + components[i].size() );
				for( unsigned int j=0 ; j<components[i].size() ; j++ ) temp.push_back( levelSet.simplexIndices[ components[i][j] ] );
			}
		}
		levelSet.simplexIndices = temp;
	}

#ifdef OUTPUT_GRADIENTS
	std::vector< Point< double , Dim > > frame[2];
	{
		Hat::ScalarFunctions< Dim > scalars( 1<<depth );

		frame[0].resize( levelSet.vertices.size() );
		frame[1].resize( levelSet.vertices.size() );

		Eigen::VectorXd x( grid.resolution() ) , y( grid.resolution() );
		for( unsigned int i=0 ; i<grid.resolution() ; i++ )
		{
			Point< double , CoDim > v = grid[i];
			x[i] = v[0] , y[i] = v[1];
		}
		double r = (double)(1<<depth);
		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ )
		{
			frame[0][i] = scalars.gradient( x , levelSet.vertices[i] / r ) / r / scale;
			frame[1][i] = scalars.gradient( y , levelSet.vertices[i] / r ) / r / scale;
		}
	}
#endif // OUTPUT_GRADIENTS

	for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) levelSet.vertices[i] = voxelToWorld( levelSet.vertices[i] );

	if( Out.set )
	{
#ifdef OUTPUT_GRADIENTS
		using Factory = VertexFactory::Factory< double , VertexFactory::PositionFactory< double , Dim > , VertexFactory::StaticFactory< double , Dim > , VertexFactory::StaticFactory< double , Dim > >;
		using Vertex = typename Factory::VertexType;
		VertexFactory::PositionFactory< double , Dim > pFactory;
		VertexFactory::StaticFactory< double , Dim > gFactory1( GradientNames1 );
		VertexFactory::StaticFactory< double , Dim > gFactory2( GradientNames2 );
		Factory factory( pFactory , gFactory1 , gFactory2 );
#else // !OUTPUT_GRADIENTS
		using Factory = VertexFactory::PositionFactory< double , Dim >;
		Factory factory;
#endif // OUTPUT_GRADIENTS
		std::vector< std::vector< int > > simplices( levelSet.simplexIndices.size() );
		for( unsigned int i=0 ; i<levelSet.simplexIndices.size() ; i++ )
		{
			simplices[i].resize( Dim-CoDim+1 );
			for( unsigned int j=0 ; j<Dim-CoDim+1 ; j++ ) simplices[i][j] = (int)levelSet.simplexIndices[i][j];
		}
#ifdef OUTPUT_GRADIENTS
		std::vector< Vertex > vertices( levelSet.vertices.size() );
		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ )
		{
			vertices[i].template get<0>() = levelSet.vertices[i];
			vertices[i].template get<1>() = frame[0][i];
			vertices[i].template get<2>() = frame[1][i];
		}
		PLY::WritePolygons< Factory , int >( Out.value , factory , vertices , simplices , PLY_BINARY_NATIVE );
#else // !OUTPUT_GRADIENTS
		PLY::WritePolygons< Factory , int >( Out.value , factory , levelSet.vertices , simplices , PLY_BINARY_NATIVE );
#endif // OUTPUT_GRADIENTS
	}

	// Output the level sets
	if( LevelSets.set )
	{
		timer.reset();
		std::vector< double > xValues( sMesh.vertices.size() ) , yValues( sMesh.vertices.size() );
		for( unsigned int i=0 ; i<sMesh.vertices.size() ; i++ )
		{
			Point< double , CoDim > v = grid( sMesh.vertices[i] );
			xValues[i] = v[0];
			yValues[i] = v[1];
		}
		if( Verbose.set ) std::cout << "Sampled implicit function: " << timer.elapsed() << " (s)" << std::endl;

		timer.reset();
		MarchingSimplices::SimplicialMesh< Dim-1 , unsigned int , Point< double , Dim > > xLevelSet , yLevelSet;
		xLevelSet = MarchingSimplices::LevelSet< double >( sMesh , [&]( unsigned int idx ){ return xValues[idx]; } , Point< double , 1 >() );
		yLevelSet = MarchingSimplices::LevelSet< double >( sMesh , [&]( unsigned int idx ){ return yValues[idx]; } , Point< double , 1 >() );
		if( Verbose.set ) std::cout << "Level set extracted: " << timer.elapsed() << " (s)" << std::endl;

		for( unsigned int i=0 ; i<xLevelSet.vertices.size() ; i++ ) xLevelSet.vertices[i] = voxelToWorld( xLevelSet.vertices[i] );
		for( unsigned int i=0 ; i<yLevelSet.vertices.size() ; i++ ) yLevelSet.vertices[i] = voxelToWorld( yLevelSet.vertices[i] );

		using Factory = VertexFactory::Factory< double , VertexFactory::PositionFactory< double , Dim > , VertexFactory::RGBColorFactory< double > >;
		using Vertex = typename Factory::VertexType;

		std::vector< Vertex > vertices;
		std::vector< std::vector< int > > simplices;
		vertices.reserve( xLevelSet.vertices.size() + yLevelSet.vertices.size() );
		simplices.reserve( xLevelSet.simplexIndices.size() + yLevelSet.simplexIndices.size() );

		Vertex v;

		v.template get<1>() = Point< double , 3 >( 223. , 41. , 53. );
		for( unsigned int i=0 ; i<xLevelSet.vertices.size() ; i++ )
		{
			v.template get<0>() = xLevelSet.vertices[i];
			vertices.push_back( v );
		}
		v.template get<1>() = Point< double , 3 >( 55. , 114. , 255. );
		for( unsigned int i=0 ; i<yLevelSet.vertices.size() ; i++ )
		{
			v.template get<0>() = yLevelSet.vertices[i];
			vertices.push_back( v );
		}

		for( unsigned int i=0 ; i<xLevelSet.simplexIndices.size() ; i++ )
		{
			SimplexIndex< Dim-1 , unsigned int > si = xLevelSet.simplexIndices[i];
			std::vector< int > _si(Dim);
			for( unsigned int j=0 ; j<=Dim-1 ; j++ ) _si[j] += (int)si[j];
			simplices.push_back( _si );
		}
		for( unsigned int i=0 ; i<yLevelSet.simplexIndices.size() ; i++ )
		{
			SimplexIndex< Dim-1 , unsigned int > si = yLevelSet.simplexIndices[i];
			std::vector< int > _si(Dim);
			for( unsigned int j=0 ; j<=Dim-1 ; j++ ) _si[j] += (int)( si[j] + (unsigned int)xLevelSet.vertices.size() );
			simplices.push_back( _si );
		}

		Factory factory;
		PLY::WritePolygons< Factory , int >( LevelSets.value , factory , vertices , simplices , PLY_BINARY_NATIVE );
	}
}

int main( int argc , char* argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}


	unsigned int dim;
	if( !RegularGrid< double , 1 >::ReadDimension( In.value , dim ) ) ERROR_OUT( "Failed to read dimension: " , In.value );
	switch( dim )
	{
		case 4: Execute< 4 >() ; break;
		case 3: Execute< 3 >() ; break;
		default: ERROR_OUT( "Only dimensions 3 and 4 supported" );
	}

	return EXIT_SUCCESS;
}
