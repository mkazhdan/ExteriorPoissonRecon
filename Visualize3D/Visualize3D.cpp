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

static const unsigned int Dim = 3;
static const unsigned int CoDim = 2;

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" ) , LevelSets( "levelSets" ) , Density( "density" );
Misha::CmdLineParameter< double > TrimDensity( "trimDensity" , 0. ) , TubularRadius( "tubular" , 0. );
Misha::CmdLineReadable NoOrient( "noOrient" ) , Verbose( "verbose" );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&LevelSets ,
	&Density ,
	&TubularRadius ,
	&TrimDensity ,
	&NoOrient ,
	&Verbose ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input grid>\n" , In.name.c_str() );
	printf( "\t[--%s <input density>]\n" , Density.name.c_str() );
	printf( "\t[--%s <output geometry>]\n" , Out.name.c_str() );
	printf( "\t[--%s <trimming density>=%f]\n" , TrimDensity.name.c_str() , TrimDensity.value );
	printf( "\t[--%s <tubular radius (in voxel units)>=%f]\n" , TubularRadius.name.c_str() , TubularRadius.value );
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
	double scale = pow( voxelToWorld(0,0) * voxelToWorld(1,1) * voxelToWorld(2,2) , 1./3 );

	double densityScale = 1. / ( 1<<depth );
	if( Density.set )
	{
		XForm< double , Dim+1 > xForm;
		density.read( Density.value , xForm );
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

	if( !NoOrient.set )
	{
		Miscellany::Timer timer;
		auto Wedge = [] ( Point< double , Dim > v1 , Point< double , Dim > v2 )
		{
			SquareMatrix< double , Dim > S;
			for( unsigned int i=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) S(i,j) = v1[i]*v2[j] - v1[j] * v2[i];
			return S;
		};

		Point< double , Dim > frame[Dim];
		for( unsigned int i=0 ; i<levelSet.simplexIndices.size() ; i++ )
		{
			auto Orthogonalize = []( Point< double , Dim > frame[] , unsigned int i )
			{
				for( unsigned j=0 ; j<i ; j++ ) frame[i] -= Point< double , Dim >::Dot( frame[i] , frame[j] ) * frame[j];
			};

			SquareMatrix< double , Dim-CoDim > D;
			for( unsigned int j=0 ; j<(Dim-CoDim) ; j++ ) for( unsigned int k=0 ; k<(Dim-CoDim) ; k++ )
			{
				Point< double , Dim > v1 = levelSet.vertices[ levelSet.simplexIndices[i][j+1] ] - levelSet.vertices[ levelSet.simplexIndices[i][0] ];
				Point< double , Dim > v2 = levelSet.vertices[ levelSet.simplexIndices[i][k+1] ] - levelSet.vertices[ levelSet.simplexIndices[i][0] ];
				D(j,k) = Point< double , Dim >::Dot( v1 , v2 );
			}

			// If the simplex is singular, skip
			if( std::abs( D.determinant() )<1e-12 ){ WARN_ONCE( "Vanishing determinant: " , i ); }
			else
			{
				// Get the orthonormal frame spanning the simplex
				for( unsigned int j=0 ; j<(Dim-CoDim) ; j++ )
				{
					frame[j] = levelSet.vertices[ levelSet.simplexIndices[i][j+1] ] - levelSet.vertices[ levelSet.simplexIndices[i][0] ];
					Orthogonalize( frame , j );
					frame[j] /= sqrt( frame[j].squareNorm() );
				}

				// Get an orthonormal frame spanning the simplex normal
				for( unsigned int j=Dim-CoDim ; j<Dim-1 ; j++ )
				{
					while( true )
					{
						frame[j] = RandomBallPoint< double , Dim >();
						Orthogonalize( frame , j );
						if( frame[j].squareNorm()>1e-10 )
						{
							frame[j] /= sqrt( frame[j].squareNorm() );
							break;
						}
					}
				}
				frame[Dim-1] = Point< double , Dim >::CrossProduct( frame );
				SquareMatrix< double , Dim > skew1 = Wedge( frame[Dim-2] , frame[Dim-1] );
				Point< double , Dim > center;
				for( unsigned int j=0 ; j<=(Dim-CoDim) ; j++ ) center += levelSet.vertices[ levelSet.simplexIndices[i][j] ];
				center /= (Dim-CoDim+1);
				Point< double , Dim > d1 , d2;
				for( unsigned int d=0 ; d<Dim ; d++ )
				{
					Point< double , 2 > partial = grid.partial( d , center );
					d1[d] = partial[0];
					d2[d] = partial[1];
				}

				SquareMatrix< double , Dim > skew2 = Wedge( d1 , d2 );
				if( SquareMatrix< double , Dim >::Dot( skew1 , skew2 )<0 ) std::swap( levelSet.simplexIndices[i][0] , levelSet.simplexIndices[i][1] );
			}
		}
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
		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ )
		{
			frame[0][i] = scalars.gradient( x , levelSet.vertices[i] / (double)(1<<depth ) );
			frame[1][i] = scalars.gradient( y , levelSet.vertices[i] / (double)(1<<depth ) );
			frame[0][i] /= sqrt( frame[0][i].squareNorm() );
			frame[1][i] /= sqrt( frame[1][i].squareNorm() );
		}
	}

	for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) levelSet.vertices[i] = voxelToWorld( levelSet.vertices[i] );

	// Output the mesh
	if( Out.set )
	{
		if( TubularRadius.value>0 )
		{
			static const unsigned int AngularSamples = 12;
			std::vector< Point< double , Dim > > vertices( levelSet.vertices.size() * AngularSamples );
			for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ )
			{
				Point< double , Dim > n[] = { frame[0][i] + frame[1][i] , frame[0][i] - frame[1][i] };
				n[0] /= sqrt( n[0].squareNorm() );
				n[1] /= sqrt( n[1].squareNorm() );
				for( unsigned int j=0 ; j<AngularSamples ; j++ )
				{
					double theta = ( 2. * M_PI * j ) / AngularSamples;
					vertices[ i*AngularSamples + j ] = levelSet.vertices[i] + ( n[0] * cos( theta ) + n[1] * sin( theta ) ) * TubularRadius.value * scale;
				}
			}
			std::vector< std::vector< int > > quads( levelSet.simplexIndices.size() * AngularSamples );
			for( unsigned int i=0 ; i<levelSet.simplexIndices.size() ; i++ )
				for( unsigned int j=0 ; j<AngularSamples ; j++ )
				{
					quads[ i*AngularSamples + j ].resize( 4 );
					quads[ i*AngularSamples + j ][0] = levelSet.simplexIndices[i][0] * AngularSamples + ( j + 0 ) % AngularSamples;
					quads[ i*AngularSamples + j ][1] = levelSet.simplexIndices[i][0] * AngularSamples + ( j + 1 ) % AngularSamples;
					quads[ i*AngularSamples + j ][2] = levelSet.simplexIndices[i][1] * AngularSamples + ( j + 1 ) % AngularSamples;
					quads[ i*AngularSamples + j ][3] = levelSet.simplexIndices[i][1] * AngularSamples + ( j + 0 ) % AngularSamples;
				}

			using Factory = VertexFactory::PositionFactory< double , Dim >;
			Factory factory;
			PLY::WritePolygons< Factory , int >( Out.value , factory , vertices , quads , PLY_BINARY_NATIVE );
		}
		else
		{
			std::vector< std::vector< int > > edges( levelSet.simplexIndices.size() );
			for( unsigned int i=0 ; i<levelSet.simplexIndices.size() ; i++ )
			{
				edges[i].resize( 2 );
				edges[i][0] = levelSet.simplexIndices[i][0];
				edges[i][1] = levelSet.simplexIndices[i][1];
			}
			using Factory = VertexFactory::PositionFactory< double , Dim >;
			Factory factory;
			PLY::WritePolygons< Factory , int >( Out.value , factory , levelSet.vertices , edges , PLY_BINARY_NATIVE );
		}
	}

	// Output the level sets
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

		if( LevelSets.set )
		{
			using Factory = VertexFactory::Factory< double , VertexFactory::PositionFactory< double , 3 > , VertexFactory::RGBColorFactory< double > >;
			using Vertex = Factory::VertexType;

			std::vector< Vertex > vertices;
			std::vector< SimplexIndex< Dim-1 , int > > simplices;
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
				SimplexIndex< Dim-1 , int > _si;
				for( unsigned int j=0 ; j<=Dim-1 ; j++ ) _si[j] += (int)si[j];
				simplices.push_back( _si );
			}
			for( unsigned int i=0 ; i<yLevelSet.simplexIndices.size() ; i++ )
			{
				SimplexIndex< Dim-1 , unsigned int > si = yLevelSet.simplexIndices[i];
				SimplexIndex< Dim-1 , int > _si;
				for( unsigned int j=0 ; j<=Dim-1 ; j++ ) _si[j] += (int)( si[j] + (unsigned int)xLevelSet.vertices.size() );
				simplices.push_back( _si );
			}

			Factory factory;
			PLY::WriteTriangles< Factory , int >( LevelSets.value , factory , vertices , simplices , PLY_BINARY_NATIVE );
		}
	}


	return EXIT_SUCCESS;
}