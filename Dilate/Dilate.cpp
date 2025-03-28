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

using namespace MishaK;

static const unsigned int Dim = 3;
static const unsigned int CoDim = 2;

CmdLineParameter< std::string > In( "in" ) , Out( "out" );
CmdLineParameter< unsigned int > AngularResolution( "res" , 12 );
CmdLineParameter< double > TubularRadius( "radius" , 1./64 );
CmdLineReadable Verbose( "verbose" );

CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&AngularResolution ,
	&TubularRadius ,
	&Verbose ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input curve>\n" , In.name.c_str() );
	printf( "\t[--%s <output tube>]\n" , Out.name.c_str() );
	printf( "\t[--%s <tubular radius (in units of diameter)>=%f]\n" , TubularRadius.name.c_str() , TubularRadius.value );
	printf( "\t[--%s <angular resolution>=%d]\n" , AngularResolution.name.c_str() , AngularResolution.value );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

const std::string GradientNames1[] = { "nx1" , "ny1" , "nz1" };
const std::string GradientNames2[] = { "nx2" , "ny2" , "nz2" };


int main( int argc , char* argv[] )
{
	Miscellany::Timer timer;

	CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	std::vector< Point< double , Dim > > frame[CoDim];

	MarchingSimplices::SimplicialMesh< Dim-CoDim , unsigned int , Point< double , Dim > > levelSet;

	{
		using Factory = VertexFactory::Factory< double , VertexFactory::PositionFactory< double , Dim > , VertexFactory::MatrixFactory< double , CoDim , Dim > >;
		using Vertex = typename Factory::VertexType;
		Factory factory( VertexFactory::PositionFactory< double , Dim >() , VertexFactory::MatrixFactory< double , CoDim , Dim >( "n" ) );

		std::vector< std::vector< int > > edges;
		std::vector< Vertex > vertices;

		int file_type = PLY::ReadPolygons< Factory >( In.value , factory , vertices , edges , NULL );

		levelSet.vertices.resize( vertices.size() );
		for( unsigned int c=0 ; c<CoDim ; c++ ) frame[c].resize( vertices.size() );
		levelSet.simplexIndices.resize( edges.size() );

		for( unsigned int i=0 ; i<vertices.size() ; i++ )
		{
			levelSet.vertices[i] = vertices[i].template get<0>();
			for( unsigned int c=0 ; c<CoDim ; c++ ) for( unsigned int d=0 ; d<Dim ; d++ ) frame[c][i][d] = vertices[i].template get<1>()(c,d);
		}

		for( unsigned int i=0 ; i<edges.size() ; i++ )
		{
			if( edges[i].size()!=Dim-CoDim+1 ) MK_ERROR_OUT( "Expected edge: " , edges[i].size() );
			for( unsigned int j=0 ; j<Dim-CoDim+1 ; j++ ) levelSet.simplexIndices[i][j] = edges[i][j];
		}
	}

	double radius=0;

	{
		Point< double , Dim > center;
		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) center += levelSet.vertices[i];
		center /= (double)levelSet.vertices.size();

		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) radius = std::max< double >( radius , Point< double , Dim >::SquareNorm( center - levelSet.vertices[i] ) );
		radius = sqrt( radius );
	}

	timer.reset();



	for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) frame[0][i] /= sqrt( frame[0][i].squareNorm() ) , frame[1][i] /= sqrt( frame[1][i].squareNorm() );

	// Output the mesh
	if( Out.set && TubularRadius.value>0 && AngularResolution.value>0 )
	{
		std::vector< Point< double , Dim > > vertices( levelSet.vertices.size() * AngularResolution.value );
		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ )
		{
			Point< double , Dim > n[] = { frame[0][i] + frame[1][i] , frame[0][i] - frame[1][i] };
			n[0] /= sqrt( n[0].squareNorm() );
			n[1] /= sqrt( n[1].squareNorm() );
			for( unsigned int j=0 ; j<AngularResolution.value ; j++ )
			{
				double theta = ( 2. * M_PI * j ) / AngularResolution.value;
				vertices[ i*AngularResolution.value + j ] = levelSet.vertices[i] + ( n[0] * cos( theta ) + n[1] * sin( theta ) ) * TubularRadius.value;
			}
		}
		std::vector< std::vector< int > > quads( levelSet.simplexIndices.size() * AngularResolution.value );
		for( unsigned int i=0 ; i<levelSet.simplexIndices.size() ; i++ )
			for( unsigned int j=0 ; j<AngularResolution.value ; j++ )
			{
				quads[ i*AngularResolution.value + j ].resize( 4 );
				quads[ i*AngularResolution.value + j ][0] = levelSet.simplexIndices[i][0] * AngularResolution.value + ( j + 0 ) % AngularResolution.value;
				quads[ i*AngularResolution.value + j ][1] = levelSet.simplexIndices[i][0] * AngularResolution.value + ( j + 1 ) % AngularResolution.value;
				quads[ i*AngularResolution.value + j ][2] = levelSet.simplexIndices[i][1] * AngularResolution.value + ( j + 1 ) % AngularResolution.value;
				quads[ i*AngularResolution.value + j ][3] = levelSet.simplexIndices[i][1] * AngularResolution.value + ( j + 0 ) % AngularResolution.value;
			}

		using Factory = VertexFactory::PositionFactory< double , Dim >;
		Factory factory;
		PLY::WritePolygons< Factory , int >( Out.value , factory , vertices , quads , PLY_BINARY_NATIVE );
	}

	return EXIT_SUCCESS;
}
