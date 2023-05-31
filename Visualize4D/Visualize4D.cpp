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
#include <map>
#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Geometry.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include "Misha/MarchingSimplices.h"
#include "Misha/RegularGrid.h"
#include "Include/Hat.h"

static const unsigned int Dim = 4;
static const unsigned int CoDim = 2;


Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" ) , LevelSets( "levelSets" ) , Density( "density" );
Misha::CmdLineParameter< double > TrimDensity( "trimDensity" , 0. ) ;
Misha::CmdLineReadable NoOrient( "noOrient" ) , Verbose( "verbose" );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&LevelSets ,
	&Density ,
	&TrimDensity ,
	&NoOrient ,
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
	printf( "\t[--%s]\n" , NoOrient.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
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

// Orients the triangles explicitly
template< typename Index >
void Orient( std::vector< SimplexIndex< 2 , Index > > &triangles );

int main( int argc , char* argv[] )
{
	Miscellany::Timer timer;

	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	RegularGrid< Point< double , CoDim > , Dim > grid;
	RegularGrid< double , Dim > density;
	XForm< double , Dim+1 > voxelToWorld;

	grid.read( In.value , voxelToWorld );
	unsigned int depth = GridDepth( grid );

	if( Density.set )
	{
		XForm< double , Dim+1 > xForm;
		density.read( Density.value , xForm );
		if( depth!=GridDepth( density ) ) ERROR_OUT( "Density resolution does not match: " , depth , " != " , GridDepth( density ) );
	}

	timer.reset();
	MarchingSimplices::SimplicialMesh< Dim , unsigned int , Point< double , Dim > > sMesh = MarchingSimplices::RegularGridTriangulation< Dim >( 1<<depth );
	if( Verbose.set ) std::cout << "Got ambient space triangulation: " << timer.elapsed() << " (s)" << std::endl;

	timer.reset();
	std::vector< Point< double , CoDim > > values( sMesh.vertices.size() );
	for( unsigned int i=0 ; i<sMesh.vertices.size() ; i++ ) values[i] = grid( sMesh.vertices[i] );
	if( Verbose.set ) std::cout << "Sampled implicit function: " << timer.elapsed() << " (s)" << std::endl;

	timer.reset();
	MarchingSimplices::SimplicialMesh< Dim-CoDim , unsigned int , Point< double , Dim > > levelSet = MarchingSimplices::LevelSet< double >( sMesh , [&]( unsigned int idx ){ return values[idx]; } , Point< double , 2 >() );
	if( Verbose.set ) std::cout << "Level set extracted: " << timer.elapsed() << " (s)" << std::endl;

	if( !NoOrient.set )
	{
		Miscellany::Timer timer;
		MarchingSimplices::Orient( levelSet.simplexIndices );
		std::cout << "Oriented: " << timer.elapsed() << " (s)" << std::endl;
	}

	if( Density.set && TrimDensity.value>0 )
	{
		std::vector< SimplexIndex< Dim-CoDim , unsigned int > > temp;
		std::vector< std::vector< size_t > > components = MarchingSimplices::FaceAdjacentConnectedComponents( levelSet.simplexIndices );
		std::cout << "Components: " << components.size() << std::endl;
		for( unsigned int i=0 ; i<components.size() ; i++ )
		{
			double maxDensity = 0;
			for( unsigned int j=0 ; j<components[i].size() ; j++ ) for( unsigned int v=0 ; v<=(Dim-CoDim) ; v++ )
				maxDensity = std::max< double >( maxDensity , density( levelSet.vertices[ levelSet.simplexIndices[ components[i][j] ][v] ] ) );

			if( Verbose.set ) std::cout << "Max Density[ " << i << " ] " << maxDensity << std::endl;

			if( maxDensity>TrimDensity.value )
			{
				temp.reserve( temp.size() + components[i].size() );
				for( unsigned int j=0 ; j<components[i].size() ; j++ ) temp.push_back( levelSet.simplexIndices[ components[i][j] ] );
			}
		}
		levelSet.simplexIndices = temp;
	}

	for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) levelSet.vertices[i] = voxelToWorld( levelSet.vertices[i] );

	auto AbsoluteValueToColor = []( double v )
	{
		v = std::min< double >( 1. , std::max< double >( 0. , v ) );
		Point< double , 3 > hsv(0.,1.,1.) , rgb;
		hsv[0] = (1-v/3) * 2.*M_PI;
		Miscellany::HSVtoRGB( &hsv[0] , &rgb[0] );
		return rgb * 255.;
	};

	auto ValueToColor = []( double v )
	{
		Point< double , 3 > rgb;
		Miscellany::TurboValueToColor( v , &rgb[0] );
		return rgb * 255.;
	};

	if( Out.set )
	{
		using Factory = VertexFactory::PositionFactory< double , Dim >;
		using Vertex = Factory::VertexType;
		Factory factory;
		PLY::WriteTriangles< Factory , unsigned int >( Out.value , factory , levelSet.vertices , levelSet.simplexIndices , PLY_BINARY_NATIVE );
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
		using Vertex = Factory::VertexType;

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

	return EXIT_SUCCESS;
}

template< typename Index >
void Orient( std::vector< SimplexIndex< 2 , Index > > &triangles )
{
	using EdgeIndex = MarchingSimplices::MultiIndex< 2 , Index >;
	using EdgeMap = typename EdgeIndex::map;
	struct EdgeData
	{
		Index t1 , t2;
		EdgeData( void ) : t1(-1) , t2(-1) {}
	};

	EdgeMap eMap;

	for( unsigned int i=0 ; i<triangles.size() ; i++ ) for( unsigned int j=0 ; j<3 ; j++ )
	{
		Index idx[] = { triangles[i][(j+1)%3] , triangles[i][(j+2)%3] };
		eMap[ EdgeIndex( idx ) ] = 0;
	}

	Index count = 0;
	for( auto iter=eMap.begin() ; iter!=eMap.end() ; iter++ ) iter->second = count++;

	std::vector< EdgeData > eData( count );

	for( unsigned int i=0 ; i<triangles.size() ; i++ ) for( unsigned int j=0 ; j<3 ; j++ )
	{
		Index idx[] = { triangles[i][(j+1)%3] , triangles[i][(j+2)%3] };
		auto iter = eMap.find( EdgeIndex(idx)  );
		if( iter==eMap.end() ) ERROR_OUT( "Could not find edge" );
		if     ( eData[ iter->second ].t1==-1 ) eData[ iter->second ].t1 = i;
		else if( eData[ iter->second ].t2==-1 ) eData[ iter->second ].t2 = i;
		else ERROR_OUT( "Both triangles occupied" );
	}

	std::vector< bool > oriented( triangles.size() , false );
	auto NeedsOrienting = [&]( void )
	{
		for( unsigned int i=0 ; i<oriented.size() ; i++ ) if( !oriented[i] ) return (Index)i;
		return (Index)-1;
	};

	Index idx;
	while( (idx=NeedsOrienting())!=-1 )
	{
		std::vector< Index > queue;
		queue.push_back( idx );

		while( queue.size() )
		{
			Index t = queue.back();
			queue.pop_back();

			for( unsigned int e=0 ; e<3 ; e++ )
			{
				Index idx[] = { triangles[t][(e+1)%3] , triangles[t][(e+2)%3] };
				auto iter = eMap.find( EdgeIndex(idx)  );
				if( iter==eMap.end() ) ERROR_OUT( "Could not find edge" );

				Index _t = t==eData[ iter->second ].t1  ? eData[ iter->second ].t2 : eData[ iter->second ].t1;
				if( _t==-1 || oriented[_t] ) continue;

				int sign1=0 , sign2=0;
				for( unsigned int j=0 ; j<3 ; j++ )
				{
					if     ( triangles[ t][(j+1)%3]==iter->first[0] && triangles[ t][(j+2)%3]==iter->first[1] ) sign1 =  1;
					else if( triangles[ t][(j+1)%3]==iter->first[1] && triangles[ t][(j+2)%3]==iter->first[0] ) sign1 = -1;
					if     ( triangles[_t][(j+1)%3]==iter->first[0] && triangles[_t][(j+2)%3]==iter->first[1] ) sign2 =  1;
					else if( triangles[_t][(j+1)%3]==iter->first[1] && triangles[_t][(j+2)%3]==iter->first[0] ) sign2 = -1;
				}
				if( !sign1 || !sign2 ) ERROR_OUT( "Could not find orientations" );
				if( sign1==sign2 ) std::swap( triangles[_t][0] , triangles[_t][1] );
				oriented[_t] = true;
				queue.push_back( _t );
			}
		}
	}
}
