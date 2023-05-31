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

#define NEW_VISUALIZE_4D_CODE

static const unsigned int Dim = 4;
static const unsigned int CoDim = 2;

enum
{
	CURVATURE_NONE ,
	CURVATURE_INTEGRATED_GAUSSIAN ,
	CURVATURE_INTEGRATED_MEAN ,
	CURVATURE_POINTWISE_GAUSSIAN ,
	CURVATURE_POINTWISE_MEAN ,
	CURVATURE_COUNT
};
const std::string CurvatureNames[] =
{
	std::string( "none" ) ,
	std::string( "gaussian (integrated)" ) ,
	std::string( "mean (integrated)" ) ,
	std::string( "gaussian (pointwise)" ) ,
	std::string( "mean (pointwise)" )
};

const double DefaultStereographic[] = { 1 , 2 , 3 , 4 };

Misha::CmdLineParameter< std::string > In( "in" ) , OutHeader( "out" ) , Density( "density" );
Misha::CmdLineParameter< double > TrimDensity( "trimDensity" , 0. ) , CurvatureThreshold( "cThreshold" , 90. );
Misha::CmdLineParameterArray< double , Dim > Stereographic( "stereo" , DefaultStereographic );
Misha::CmdLineParameter< int > CurvatureType( "cType" , CURVATURE_INTEGRATED_GAUSSIAN );
Misha::CmdLineReadable Normalize( "normalize" ) , NoOrient( "noOrient" ) , Verbose( "verbose" ) , VoxelCoordinates( "vCoordinates" ) , AbsoluteCurvature( "abs" );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&OutHeader ,
	&Density ,
	&TrimDensity ,
	&Normalize ,
	&Stereographic ,
	&CurvatureThreshold ,
	&NoOrient ,
	&AbsoluteCurvature ,
	&CurvatureType ,
	&VoxelCoordinates ,
	&Verbose ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input grid>\n" , In.name.c_str() );
	printf( "\t[--%s <output header>]\n" , OutHeader.name.c_str() );
	printf( "\t[--%s <density grid>]\n" , Density.name.c_str() );
	printf( "\t[--%s <trim value>=%f]\n" , TrimDensity.name.c_str() , TrimDensity.value );
	printf( "\t[--%s <stereographic projection axis>]\n" , Stereographic.name.c_str() );
	printf( "\t[--%s <curvature threshold (in degrees)>=%f]\n" , CurvatureThreshold.name.c_str() , CurvatureThreshold.value );
	printf( "\t[--%s <curvature type>=%d]\n" , CurvatureType.name.c_str() , CurvatureType.value );
	for( unsigned int i=0 ; i<CURVATURE_COUNT ; i++ ) printf( "\t\t%d] %s\n" , i , CurvatureNames[i].c_str() );
	printf( "\t[--%s]\n" , Normalize.name.c_str() );
	printf( "\t[--%s]\n" , NoOrient.name.c_str() );
	printf( "\t[--%s]\n" , AbsoluteCurvature.name.c_str() );
	printf( "\t[--%s]\n" , VoxelCoordinates.name.c_str() );
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

// A mapping taking the sphere to the plane through the origin perpendicular to a
Point< double , 4 > StereographicProjection( Point< double , 4 > p , Point< double , 4 > a );

// Comptues the Gaussian curvature (angle deficit) at each vertex
template< unsigned int Dim , typename Index >
std::vector< double > GaussianCurvature( const std::vector< Point< double , Dim > > &vertices , const std::vector< SimplexIndex< 2 , Index > > &triangles );

// Comptues the mean curvature (edge-weighted average of dihedral angles of incident edges) at each vertex
template< unsigned int Dim , typename Index >
std::vector< double > MeanCurvature( const std::vector< Point< double , Dim > > &vertices , const std::vector< SimplexIndex< 2 , Index > > &triangles );

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
	if( VoxelCoordinates.set ) voxelToWorld = XForm< double , Dim+1 >::Identity();
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

	std::vector< double > curvatureValues;

	switch( CurvatureType.value )
	{
		case CURVATURE_NONE:
			break;
		case CURVATURE_INTEGRATED_GAUSSIAN:
		case CURVATURE_POINTWISE_GAUSSIAN:
			curvatureValues = GaussianCurvature( levelSet.vertices , levelSet.simplexIndices );
			break;
		case CURVATURE_INTEGRATED_MEAN:
		case CURVATURE_POINTWISE_MEAN:
			curvatureValues = MeanCurvature( levelSet.vertices , levelSet.simplexIndices ) ; break;
		default: ERROR_OUT( "Unrecognized curvature type" );
	}
	for( unsigned int i=0 ; i<curvatureValues.size() ; i++ ) curvatureValues[i] *= 180. / M_PI;
	if( CurvatureType.value==CURVATURE_POINTWISE_GAUSSIAN || CurvatureType.value==CURVATURE_POINTWISE_MEAN )
	{
		std::vector< double > areas( levelSet.vertices.size() , 0 );
		for( unsigned int i=0 ; i<levelSet.simplexIndices.size() ; i++ )
		{
			Simplex< double , Dim , 2 > tri;
			for( unsigned int k=0 ; k<=2 ; k++ ) tri[k] = levelSet.vertices[ levelSet.simplexIndices[i][k] ];
			double area = tri.measure() / 3.;
			for( unsigned int k=0 ; k<=2 ; k++ ) areas[ levelSet.simplexIndices[i][k] ] += area;
		}
		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) curvatureValues[i] /= areas[i];
	}

	if( CurvatureThreshold.value<=0 )
	{
		double dev = 0;
		for( unsigned int i=0 ; i<curvatureValues.size() ; i++ ) dev += curvatureValues[i] * curvatureValues[i];
		dev = sqrt( dev/curvatureValues.size() );
		for( unsigned int i=0 ; i<curvatureValues.size() ; i++ ) curvatureValues[i] /= dev;
		CurvatureThreshold.value = 1.;
		std::cout << "Deviation: " << dev << std::endl;
	}

	if( CurvatureType.value!=CURVATURE_NONE )
	{

		double avg = 0 , dev = 0 , max = 0;
		for( unsigned int i=0 ; i<curvatureValues.size() ; i++ ) avg += curvatureValues[i] , max = std::max< double >( max , std::abs( curvatureValues[i] ) );
		avg /= curvatureValues.size();
		for( unsigned int i=0 ; i<curvatureValues.size() ; i++ ) dev += curvatureValues[i] * curvatureValues[i];
		dev = sqrt( dev / curvatureValues.size() );
		std::cout << "Curvature average: "   << avg << std::endl;
		std::cout << "Curvature deviation: " << dev << std::endl;
		std::cout << "Curvature maximum: "   << max << std::endl;
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

		double radius = Radius();
		for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ ) levelSet.vertices[i] /= radius;
	}

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

	if( OutHeader.set )
	{
		std::string fileName = OutHeader.value + ".4D.ply";
		if( CurvatureType.value==CURVATURE_NONE )
		{
			using Factory = VertexFactory::PositionFactory< double , Dim >;
			using Vertex = Factory::VertexType;
			Factory factory;
			PLY::WriteTriangles< Factory , unsigned int >( fileName , factory , levelSet.vertices , levelSet.simplexIndices , PLY_BINARY_NATIVE );
		}
		else
		{
			using Factory = VertexFactory::Factory< double , VertexFactory::PositionFactory< double , Dim > , VertexFactory::RGBColorFactory< double > >;
			using Vertex = Factory::VertexType;
			std::vector< Vertex > vertices( levelSet.vertices.size() );
			for( unsigned int i=0 ; i<levelSet.vertices.size() ; i++ )
			{
				vertices[i].template get<0>() = levelSet.vertices[i];
				vertices[i].template get<1>() = AbsoluteCurvature.set ? AbsoluteValueToColor( std::abs( curvatureValues[i] ) / CurvatureThreshold.value ) : ValueToColor( curvatureValues[i] / CurvatureThreshold.value );
			}

			Factory factory;

			PLY::WriteTriangles< Factory , unsigned int >( fileName , factory , vertices , levelSet.simplexIndices , PLY_BINARY_NATIVE );
		}
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

		if( OutHeader.set )
		{
			std::string fileName = OutHeader.value + ".stereo.ply";
			std::vector< SimplexIndex< Dim-CoDim , int > > simplices( levelSet.simplexIndices.size() );
			for( unsigned int i=0 ; i<levelSet.simplexIndices.size() ; i++ ) for( unsigned int j=0 ; j<=(Dim-CoDim) ; j++ ) simplices[i][j] = (int)levelSet.simplexIndices[i][j];
			if( CurvatureType.value==CURVATURE_NONE )
			{
				using Factory = VertexFactory::PositionFactory< double , Dim-1 >;
				using Vertex = Factory::VertexType;
				Factory factory;

				PLY::WriteTriangles< Factory , int >( fileName , factory , vertices , simplices , PLY_BINARY_NATIVE );
			}
			else
			{
				using Factory = VertexFactory::Factory< double , VertexFactory::PositionFactory< double , Dim-1 > , VertexFactory::RGBColorFactory< double > >;
				using Vertex = Factory::VertexType;
				std::vector< Vertex > _vertices( vertices.size() );
				for( unsigned int i=0 ; i<vertices.size() ; i++ )
				{
					_vertices[i].template get<0>() = vertices[i];
					_vertices[i].template get<1>() = AbsoluteCurvature.set ? AbsoluteValueToColor( std::abs( curvatureValues[i] ) / CurvatureThreshold.value ) : ValueToColor( curvatureValues[i] / CurvatureThreshold.value );
				}

				Factory factory;

				PLY::WriteTriangles< Factory , int >( fileName , factory , _vertices , simplices , PLY_BINARY_NATIVE );
			}
		}
	}


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

// Computes the Gaussian curvature (angle deficit) at each vertex
template< unsigned int D , typename Index >
std::vector< double > GaussianCurvature( const std::vector< Point< double , D > > &vertices , const std::vector< SimplexIndex< 2 , Index > > &triangles )
{
	std::vector< double > gCurvatures( vertices.size() , 0 );
	for( unsigned int i=0 ; i<triangles.size() ; i++ )
	{
		for( unsigned int j=0 ; j<3 ; j++ )
		{
			Point< double , Dim > v[3];
			for( unsigned int k=0 ; k<3 ; k++ ) v[k] = vertices[ triangles[i][(j+k)%3] ];

			Point< double , Dim > d[] = { v[1]-v[0] , v[2]-v[0] };
			d[0] /= sqrt( d[0].squareNorm() ) , d[1] /= sqrt( d[1].squareNorm() );
			gCurvatures[ triangles[i][j] ] += acos( Point< double , Dim >::Dot( d[0] , d[1] ) );
		}
	}
	for( unsigned int i=0 ; i<gCurvatures.size() ; i++ ) gCurvatures[i] = ( 2. * M_PI - gCurvatures[i] );
	return gCurvatures;
}

template< unsigned int Dim , typename Index >
std::vector< double > MeanCurvature( const std::vector< Point< double , Dim > > &vertices , const std::vector< SimplexIndex< 2 , Index > > &_triangles )
{
	using EdgeIndex = MarchingSimplices::MultiIndex< 2 , Index >;
	using EdgeMap = typename EdgeIndex::map;
	struct EdgeData
	{
		Index t1 , t2;
		double area;
		EdgeData( void ) : t1(-1) , t2(-1) , area(0) {}
	};

	std::vector< SimplexIndex< 2 , Index > > triangles = _triangles;
	Orient( triangles );

	auto Frame = [&]( Index i )
	{
		Point< double , Dim > d[] = { vertices[ triangles[i][1] ] - vertices[ triangles[i][0] ] , vertices[ triangles[i][2] ] - vertices[ triangles[i][0] ] };
		return Hat::Wedge( d[0] , d[1] );
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

	for( unsigned int i=0 ; i<triangles.size() ; i++ )
	{
		Simplex< double , Dim , 2 > tri;
		for( unsigned int k=0 ; k<=2 ; k++ ) tri[k] = vertices[ triangles[i][k] ];
		double area = tri.measure()/3.;
		for( unsigned int j=0 ; j<3 ; j++ )
		{
			Index idx[] = { triangles[i][(j+1)%3] , triangles[i][(j+2)%3] };
			auto iter = eMap.find( EdgeIndex(idx)  );
			if( iter==eMap.end() ) ERROR_OUT( "Could not find edge" );
			if     ( eData[ iter->second ].t1==-1 ) eData[ iter->second ].t1 = i , eData[ iter->second ].area += area;
			else if( eData[ iter->second ].t2==-1 ) eData[ iter->second ].t2 = i , eData[ iter->second ].area += area;
			else ERROR_OUT( "Both triangles occupied" );
		}
	}

	std::vector< double > mCurvatures( vertices.size() , 0 ) , weights( vertices.size() , 0 );
	for( auto iter=eMap.begin() ; iter!=eMap.end() ; iter++ )
	{
		Index e = iter->second;
		unsigned int t1 = eData[e].t1 , t2 = eData[e].t2;
		if( t2==-1 ) WARN( "Boundary edge: " , t1 );
		else
		{
			int sign1=0 , sign2=0;
			for( unsigned int j=0 ; j<3 ; j++ )
			{
				if     ( triangles[t1][(j+1)%3]==iter->first[0] && triangles[t1][(j+2)%3]==iter->first[1] ) sign1 =  1;
				else if( triangles[t1][(j+1)%3]==iter->first[1] && triangles[t1][(j+2)%3]==iter->first[0] ) sign1 = -1;
				if     ( triangles[t2][(j+1)%3]==iter->first[0] && triangles[t2][(j+2)%3]==iter->first[1] ) sign2 =  1;
				else if( triangles[t2][(j+1)%3]==iter->first[1] && triangles[t2][(j+2)%3]==iter->first[0] ) sign2 = -1;
			}
			if( !sign1 || !sign2 ) ERROR_OUT( "Could not find orientations" );
			Hat::SkewSymmetricMatrix< double , Dim > skew1 = Frame( eData[e].t1 ) , skew2 = Frame( eData[e].t2 );
			double n1 = skew1.squareNorm() , n2 = skew2.squareNorm();
			if( !n1 || !n2 ) WARN( "Vanishing frame" );
			else
			{
				n1 = sqrt(n1) , n2 = sqrt(n2);
				skew1 /= n1 , skew2 /= n2;
				if( sign1 * sign2>0 ) ERROR_OUT( "Bad orientation" );
				double theta , dot = Hat::SkewSymmetricMatrix< double , Dim >::Dot( skew1 , skew2 );
				if     ( dot>= 1. ) theta = 0;
				else if( dot<=-1. ) theta = M_PI;
				else                theta = acos( dot );
				Point< double , Dim > c = ( vertices[ triangles[t1][0] ] + vertices[ triangles[t1][1] ] + vertices[ triangles[t1][2] ] ) / 3.;
				Hat::SkewSymmetricMatrix< double , Dim > skew = Hat::Wedge( vertices[ iter->first[0] ]-c , vertices[ iter->first[1] ]-c );
				Hat::SkewSymmetricMatrix< double , Dim > v = skew1 - skew2;

				if( Hat::SkewSymmetricMatrix< double , Dim >::Dot( v , skew )<0 ) theta = -theta;

				weights[ iter->first[0] ]++;
				weights[ iter->first[1] ]++;
				mCurvatures[ iter->first[0] ] += theta;
				mCurvatures[ iter->first[1] ] += theta;
			}
		}
	}
	for( unsigned int i=0 ; i<mCurvatures.size() ; i++ ) mCurvatures[i] /= weights[i];

	return mCurvatures;
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
