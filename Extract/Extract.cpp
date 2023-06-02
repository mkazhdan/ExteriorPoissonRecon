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

#if 0
template< unsigned int Dim , unsigned int CoDim >
std::vector< std::string > GradientNames( void )
{
	const std::string _names[]  = { "nx" , "ny" , "nz" , "nw" };
	std::vector< std::string > names( Dim * CoDim );
	for( unsigned int c=0 ; c<CoDim ; c++ ) for( unsigned int d=0 ; d<Dim ; d++ ) names[d+c*Dim] = _names[d] + std::to_string( c+1 );
	return names;
}
#endif

Point< double , 3 > SurfaceColor( unsigned int idx )
{
	switch( idx )
	{
		case 0: return Point< double , 3 >( 223. , 41. , 53. );
		case 1: return Point< double , 3 >( 55. , 114. , 255. );
		default: ERROR_OUT( "Only two surfaces supported" );
	}
	return Point< double , 3 >();
}

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

template< unsigned int Dim , unsigned int CoDim >
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
		if( Verbose.set ) std::cout << "Oriented: " << timer.elapsed() << " (s)" << std::endl;
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

	std::vector< Point< double , Dim > > frame[CoDim];
	{
		Hat::ScalarFunctions< Dim > scalars( 1<<depth );

		double r = (double)(1<<depth);
		for( unsigned int c=0 ; c<CoDim ; c++ )
		{
			frame[c].resize( levelSet.vertices.size() );
			Eigen::VectorXd x( grid.resolution() );
			for( unsigned int i=0 ; i<grid.resolution() ; i++ ) x[i] = grid[i][c];
			for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) frame[c][i] = scalars.gradient( x , levelSet.vertices[i] / r ) / r / scale;
		}
	}

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
	if( LevelSets.set )
	{
		timer.reset();

		using Factory = VertexFactory::Factory< double , VertexFactory::PositionFactory< double , Dim > , VertexFactory::RGBColorFactory< double > >;
		using Vertex = typename Factory::VertexType;

		std::vector< Vertex > vertices;
		std::vector< std::vector< int > > simplices;

		for( unsigned int c=0 ; c<CoDim ; c++ )
		{
			timer.reset();
			std::vector< double > values( sMesh.vertices.size() );
			for( unsigned int i=0 ; i<sMesh.vertices.size() ; i++ )
			{
				Point< double , CoDim > v = grid( sMesh.vertices[i] );
				values[i] = v[c];
			}
			if( Verbose.set ) std::cout << "Sampled implicit function: " << timer.elapsed() << " (s)" << std::endl;

			timer.reset();
			MarchingSimplices::SimplicialMesh< Dim-1 , unsigned int , Point< double , Dim > > levelSet;
			levelSet = MarchingSimplices::LevelSet< double >( sMesh , [&]( unsigned int idx ){ return values[idx]; } , Point< double , 1 >() );
			if( Verbose.set ) std::cout << "Level set extracted: " << timer.elapsed() << " (s)" << std::endl;

			for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) levelSet.vertices[i] = voxelToWorld( levelSet.vertices[i] );

			size_t vSize = vertices.size();
			vertices.reserve( vertices.size() + levelSet.vertices.size() );
			simplices.reserve( simplices.size() + levelSet.simplexIndices.size() );

			Vertex v;
			v.template get<1>() = SurfaceColor( c );
			for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ )
			{
				v.template get<0>() = levelSet.vertices[i];
				vertices.push_back( v );
			}
			for( unsigned int i=0 ; i<levelSet.simplexIndices.size() ; i++ )
			{
				SimplexIndex< Dim-1 , unsigned int > si = levelSet.simplexIndices[i];
				std::vector< int > _si(Dim);
				for( unsigned int j=0 ; j<=Dim-1 ; j++ ) _si[j] = (int)( si[j]  + (unsigned int)vSize );
				simplices.push_back( _si );
			}
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
		case 4: Execute< 4 , CoDim >() ; break;
		case 3: Execute< 3 , CoDim >() ; break;
		default: ERROR_OUT( "Only dimensions 3 and 4 supported" );
	}

	return EXIT_SUCCESS;
}
