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
#if !defined( __clang__ )
#include <omp.h>
#endif // !__clang__
#include <vector>
#include <unordered_map>

#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Geometry.h"
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include "Include/Samples.h"

Misha::CmdLineParameter< std::string > Out( "out" );
Misha::CmdLineParameter< std::string > SampleType( "type" );
Misha::CmdLineParameter< unsigned int > SampleResolution( "res" , 1024 );
Misha::CmdLineParameter< double > AngularNoise( "aNoise" , 0. ) , PositionalNoise( "pNoise" , 0. ) , SpotNoise( "sNoise" , 0. );
Misha::CmdLineReadable RegularSample( "regular" );

Misha::CmdLineReadable* params[] =
{
	&SampleType ,
	&SampleResolution ,
	&Out ,
	&RegularSample ,
	&AngularNoise ,
	&PositionalNoise ,
	&SpotNoise ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <sample type>\n" , SampleType.name.c_str() );
	for( unsigned int i=0 ; i<Samples::  Curve::TypeCount ; i++ ) printf( "\t\t%s\n" , Samples::  Curve::Descriptions[i].c_str() );
	for( unsigned int i=0 ; i<Samples::Surface::TypeCount ; i++ ) printf( "\t\t%s\n" , Samples::Surface::Descriptions[i].c_str() );
	printf( "\t[--%s <sample resolution>=%d]\n" , SampleResolution.name.c_str() , SampleResolution.value );
	printf( "\t[--%s <angular noise sigma (in units of degrees)>=%f]\n" , AngularNoise.name.c_str() , AngularNoise.value );
	printf( "\t[--%s <positional noise sigma (in normalized units)>=%f]\n" , PositionalNoise.name.c_str() , PositionalNoise.value );
	printf( "\t[--%s <spot noise (fraction of input samples)>=%f]\n" , SpotNoise.name.c_str() , SpotNoise.value );
	printf( "\t[--%s <output file>]\n" , Out.name.c_str() );
	printf( "\t[--%s]\n" , RegularSample.name.c_str() );
}

template< unsigned int Dim >
std::pair< Point< double , Dim > , Point< double , Dim > > FitBlade( Hat::SkewSymmetricMatrix< double , Dim > S , unsigned int &steps , double eps=1e-10 )
{
	if( S.squareNorm()==0 ) return std::make_pair( Point< double , Dim >() , Point< double , Dim >() );

	Point< double , Dim > v1 , v2;
	v1 = RandomSpherePoint< double , Dim >();
	v1 /= sqrt( v1.squareNorm() );

	auto Next = [&]( Point< double , Dim > v )
	{
		Point< double , Dim > _v = S() * v;
		_v -= Point< double , Dim >::Dot( v , _v ) * v;
		return _v / sqrt( _v.squareNorm() );
	};

	for( unsigned int s=0 ; s<steps ; s++ )
	{
		v2 = Next( v1 );

		Point< double , Dim > _v1 = -Next( v2 );

		if( ( v1 - _v1 ).squareNorm() < eps ){ steps = s ; break; }
		v1 = _v1;
	}
	double scale = sqrt( sqrt( S.squareNorm() ) );
	return std::make_pair( v1*scale , v2*scale );
}

template< unsigned int Dim , typename SampleFunctor /* = std::function< Point< double , Dim > > ( unsigned int idx ) */ >
SquareMatrix< double , Dim+1 > ToUnitCube( SampleFunctor F , size_t sampleNum )
{
	if( !sampleNum ) ERROR_OUT( "Expected some samples: " , sampleNum );

	SquareMatrix< double , Dim+1 > xForm = SquareMatrix< double , Dim+1 >::Identity();
	Point< double , Dim > bBox[2];
	bBox[0] = bBox[1] = F(0);
	for( unsigned int i=1 ; i<sampleNum ; i++ ) for( unsigned int d=0 ; d<Dim ; d++ )
	{
		Point< double , Dim > p = F(i);
		bBox[0][d] = std::min< double >( bBox[0][d] , p[d] );
		bBox[1][d] = std::max< double >( bBox[1][d] , p[d] );
	}

	double _scale = 0;
	{
		// use shape diameter
		Point< double, Dim > center;
		for( int i=0 ; i<sampleNum ; i++ ) center += F(i);
		center /= (double)sampleNum;
		for( int i=0 ; i<sampleNum ; i++ ) _scale = std::max< double >( _scale , ( center-F(i) ).squareNorm() );
		_scale = sqrt( _scale ) * 2;
	}

	SquareMatrix< double , Dim+1 > t1 = SquareMatrix< double , Dim+1 >::Identity() , s = SquareMatrix< double , Dim+1 >::Identity() , t2 = SquareMatrix< double , Dim+1 >::Identity();
	for( unsigned int d=0 ; d<Dim ; d++ ) t2(Dim,d) = -( bBox[0][d] + bBox[1][d] ) / 2.;
	for( unsigned int d=0 ; d<Dim ; d++ ) s(d,d) = 1. / ( _scale );
	for( unsigned int d=0 ; d<Dim ; d++ ) t1(Dim,d) = 0.5;
	xForm = t1 * s * t2;
	return xForm;
}

template< unsigned int Dim >
void AddNoise( std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > >& hermiteData )
{
	// Add positional noise
	if( PositionalNoise.value>0 )
	{
		std::default_random_engine generator;
		std::normal_distribution< double > distribution( 0. , PositionalNoise.value );
		for( unsigned int i=0 ; i<hermiteData.size() ; i++ )
		{
			Point< double , Dim > delta;
			while( true )
			{
				bool inBounds = true;
				for( unsigned int d=0 ; d<Dim && inBounds ; d++ )
				{
					delta[d] = distribution( generator );
					inBounds &= hermiteData[i].first[d] + delta[d]>0 && hermiteData[i].first[d] + delta[d]<1;
				}
				if( inBounds ) break;
			}
			hermiteData[i].first += delta;
		}
	}

	// Add angular noise
	if( AngularNoise.value>0 )
	{
		std::default_random_engine generator;
		std::normal_distribution< double > distribution( 0. , AngularNoise.value * M_PI/180. );

		for( unsigned int i=0 ; i<hermiteData.size() ; i++ )
		{
			SquareMatrix< double , Dim > R = Samples::RandomRotation< Dim >( [&]( void ){ return distribution( generator ); } );
			hermiteData[i].second = R.transpose() * hermiteData[i].second() * R;
		}
	}

	// Add spot noise
	if( SpotNoise.value )
	{
		std::default_random_engine generator;
		std::uniform_real_distribution< double > pDistribution( 0. , 1. );
		std::uniform_real_distribution< double > sDistribution( -1. , 1. );

		unsigned int samples = (unsigned int)ceil( hermiteData.size() * SpotNoise.value );
		hermiteData.reserve( hermiteData.size() + samples );
		std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > sample;
		for( unsigned int i=0 ; i<samples ; i++ )
		{
			for( unsigned int d=0 ; d<Dim ; d++ ) sample.first[d] = pDistribution( generator );
			SquareMatrix< double , Dim > skew;
			for( unsigned int j=0 ; j<Dim ; j++ ) for( unsigned int k=0 ; k<j ; k++ )
			{
				skew(j,k) = sDistribution( generator );
				skew(k,j) = -skew(j,k);
			}
			sample.second = Hat::SkewSymmetricMatrix< double , Dim >( skew );
			hermiteData.push_back( sample );
		}
	}

}

template< unsigned int Dim >
void Execute( std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > hermiteData )
{
	// Add noise (after remapping to the unit cube)
	{
		SquareMatrix< double , Dim+1 > worldToUnitCube = ToUnitCube< Dim >( [&]( unsigned int idx ){ return hermiteData[idx].first; } , hermiteData.size() );
		SquareMatrix< double , Dim+1 > unitCubeToWorld  = worldToUnitCube.inverse();

		for( unsigned int i=0 ; i<hermiteData.size() ; i++ ) hermiteData[i].first = worldToUnitCube( hermiteData[i].first );
		AddNoise( hermiteData );
		for( unsigned int i=0 ; i<hermiteData.size() ; i++ ) hermiteData[i].first = unitCubeToWorld( hermiteData[i].first );
	}

	// generate skew symmetric matrix header
	using Factory = VertexFactory::Factory< double , VertexFactory::PositionFactory< double, Dim > , VertexFactory::MatrixFactory< double , Dim , Dim > >;
	Factory factory = Factory( VertexFactory::PositionFactory< double, Dim >(), VertexFactory::MatrixFactory< double , Dim , Dim >( "skew" ) );
	std::vector< typename Factory::VertexType > vertices( hermiteData.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i].template get<0>() = hermiteData[i].first , vertices[i].template get<1>() = hermiteData[i].second();
	if( Out.set ) PLY::WritePoints< Factory >( Out.value , factory , vertices , PLY_BINARY_NATIVE );
}

int main( int argc , char* argv[] )
{
	Misha::CmdLineParse( argc-1 , argv+1 , params );
	if( !SampleType.set )
	{
		ShowUsage( argv[0] );
		return EXIT_SUCCESS;
	}

	try{ Execute< 3 >( Samples::Curve::GetSamples( SampleResolution.value , SampleType.value , !RegularSample.set ) ); }
	catch( const Samples::InvalidSamplerException & )
	{
		try{ Execute< 4 >( Samples::Surface::GetSamples( SampleResolution.value , SampleType.value , !RegularSample.set ) ); }
		catch( const Samples::InvalidSamplerException & ){ ERROR_OUT( "Could not parse sample type: " , SampleType.value ); }
	}

	return EXIT_SUCCESS;
}