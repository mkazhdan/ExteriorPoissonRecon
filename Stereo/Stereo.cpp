/*
Copyright (c) 2025, Michael Kazhdan, Sing-Chun Lee, Marc Alexa, and Maximilian Kohlbrenner
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

using namespace MishaK;

static const unsigned int Dim = 4;
static const unsigned int CoDim = 2;


const double DefaultStereographic[] = { 0 , 0 , 0 , 1 };

CmdLineParameter< std::string > In( "in" ) , Out( "out" );
CmdLineParameterArray< double , Dim > Stereographic( "stereo" , DefaultStereographic );
CmdLineReadable Normalize( "normalize" ) , Verbose( "verbose" );

CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&Normalize ,
	&Stereographic ,
	&Verbose ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input 4D mesh>\n" , In.name.c_str() );
	printf( "\t[--%s <output 3D mesh>]\n" , Out.name.c_str() );
	printf( "\t[--%s <stereographic projection axis>]\n" , Stereographic.name.c_str() );
	printf( "\t[--%s]\n" , Normalize.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}


// A mapping taking the sphere to the plane through the origin perpendicular to a
Point< double , 4 > StereographicProjection( Point< double , 4 > p , Point< double , 4 > a );

int main( int argc , char* argv[] )
{
	Miscellany::Timer timer;

	CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	MarchingSimplices::SimplicialMesh< Dim-CoDim , unsigned int , Point< double , Dim > > levelSet;
	int file_type;
	{
		using Factory = VertexFactory::PositionFactory< double , Dim >;
		using Vertex = Factory::VertexType;
		Factory factory;
		bool *propertiesFlag = NULL;
		file_type = PLY::ReadTriangles< Factory , unsigned int >( In.value , factory , levelSet.vertices , levelSet.simplexIndices , propertiesFlag );
	}

	if( Normalize.set )
	{
		auto Center = [&]( void )
		{
			Point< double , Dim > center;
			for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) center += levelSet.vertices[i];
			return center / (double)levelSet.vertices.size();
		};

		auto Radius = [&]( void )
		{
			double radius = 0;
			for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) radius += levelSet.vertices[i].squareNorm();
			return sqrt( radius / (double)levelSet.vertices.size() );
		};

		Point< double , Dim > center = Center();
		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) levelSet.vertices[i] -= center;
		if( Verbose.set ) std::cout << "Center: " << center << std::endl;

		double radius = Radius();
		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) levelSet.vertices[i] /= radius;
		if( Verbose.set ) std::cout << "Radius: " << radius << std::endl;
	}


	Point< double , Dim > stereographicAxis;
	for( unsigned int d=0 ; d<Dim ; d++ ) stereographicAxis[d] = Stereographic.values[d];

	if( stereographicAxis.squareNorm() )
	{
		std::vector< Point< double , Dim-1 > > vertices( levelSet.vertices.size() );
		Point< double,  4 > frame[Dim];

		// The first vector is the perpendicular to the plane onto which the stereographic projection maps the sphere
		frame[0] = stereographicAxis / sqrt( stereographicAxis.squareNorm() );

		// Compute the remaining frames by using coordinate axis and Gram-Schmidt (initially checking that frame[0] is not itself a coordinate axis
		unsigned int count = 0;
		for( unsigned int d=0 ; d<Dim ; d++ ) if( frame[0][d] ) count++;
		if( count==1 )
		{
			unsigned int dir = -1;
			for( unsigned int d=0 ; d<Dim ; d++ ) if( frame[0][d] ) dir=d;
			if( dir==-1 ) ERROR_OUT( "Could not determine direction" );
			for( unsigned int d=1 ; d<Dim ; d++ ) frame[d][(dir+d)%Dim] = 1;
		}
		else
		{
			for( unsigned int d=1 ; d<Dim ; d++ )
			{
				frame[d][d] = 1.;
				for( unsigned int _d=0 ; _d<d ; _d++ ) frame[d] -= Point< double , 4 >::Dot( frame[_d] , frame[d] ) * frame[_d];
				frame[d] /= sqrt( frame[d].squareNorm() );
			}
		}
		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ )
		{
			Point< double , Dim > v = StereographicProjection( levelSet.vertices[i] , frame[0] );
			for( unsigned int d=0 ; d<3 ; d++ ) vertices[i][d] = Point< double , Dim >::Dot( v , frame[d+1] );
		}

		if( Out.set )
		{
			std::vector< SimplexIndex< Dim-CoDim , int > > simplices( levelSet.simplexIndices.size() );
			for( unsigned int i=0 ; i<levelSet.simplexIndices.size() ; i++ ) for( unsigned int j=0 ; j<=(Dim-CoDim) ; j++ ) simplices[i][j] = (int)levelSet.simplexIndices[i][j];
			{
				using Factory = VertexFactory::PositionFactory< double , Dim-1 >;
				using Vertex = Factory::VertexType;
				Factory factory;

				PLY::WriteTriangles< Factory , int >( Out.value , factory , vertices , simplices , PLY_BINARY_NATIVE );
			}
		}
	}
	else ERROR_OUT( "Vanishing stereographic axis: " , stereographicAxis );


	return EXIT_SUCCESS;
}

// From https://en.wikipedia.org/wiki/M%C3%B6bius_transformation
Point< double , Dim > StereographicProjection( Point< double , Dim > p , Point< double , Dim > a , Point< double , Dim > b , double alpha , bool denom )
{
	return b + alpha * ( p - a ) / ( denom ? (p-a).squareNorm() : 1. );
}
Point< double , Dim > StereographicProjection( Point< double , Dim > p , Point< double , Dim > a )
{
	// a -> infinity
	// p = -a -> ( 0 , 0 , 0 , 0 )
	//	 => b + alpha * ( p - a ) / || p - a ||^2 -> ( 0 , 0 , 0 , 0 )
	//	<=> b - 2 * alpha * a / 4 -> ( 0 , 0 , 0 , 0 )
	//	<=> b = alpha * a / 2
	// Given p s.t. < p , a > = 0
	//	< F(p) , a > = < b + alpha * ( p - a ) / || p - a ||^2 , a >
	//	             = alpha * < a/2 + ( p - a ) / || p - a ||^2 , a >
	//	             = alpha * ( 0.5 + < ( p - a ) / || p - a ||^2 , a >
	//	             = alpha * ( 0.5 - 1. / || p - a ||^2 )
	//	             = alpha * ( 0.5 - 1. / ( || p ||^2 + || a ||^2 )
	//	             = 0 if || p ||^2 = 0
	a /= sqrt( a.squareNorm() );
	return StereographicProjection( p , a , a/2 , 1. , true );
}