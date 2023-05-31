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


const double DefaultStereographic[] = { 0 , 0 , 0 , 1 };

Misha::CmdLineParameter< std::string > In( "in" ) , Out( "out" );
Misha::CmdLineParameterArray< double , Dim > Stereographic( "stereo" , DefaultStereographic );
Misha::CmdLineReadable Normalize( "normalize" ) , NoOrient( "noOrient" ) , Verbose( "verbose" );

Misha::CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&Normalize ,
	&Stereographic ,
	&NoOrient ,
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
	printf( "\t[--%s]\n" , NoOrient.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}


// A mapping taking the sphere to the plane through the origin perpendicular to a
Point< double , 4 > StereographicProjection( Point< double , 4 > p , Point< double , 4 > a );

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

	MarchingSimplices::SimplicialMesh< Dim-CoDim , unsigned int , Point< double , Dim > > levelSet;
	int file_type;
	{
		using Factory = VertexFactory::PositionFactory< double , Dim >;
		using Vertex = Factory::VertexType;
		Factory factory;
		bool *propertiesFlag = NULL;
		PLY::ReadTriangles< Factory , unsigned int >( In.value , factory , levelSet.vertices , levelSet.simplexIndices , propertiesFlag , file_type );
	}

	if( !NoOrient.set )
	{
		Miscellany::Timer timer;
		MarchingSimplices::Orient( levelSet.simplexIndices );
		std::cout << "Oriented: " << timer.elapsed() << " (s)" << std::endl;
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
