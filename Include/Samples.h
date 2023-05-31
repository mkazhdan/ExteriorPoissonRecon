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

#ifndef CURVE_SAMPLES_INCLUDED
#define CURVE_SAMPLES_INCLUDED

#include <vector>
#include <random>
#include "Include/Hat.h"
#include "Misha/Miscellany.h"
#include "Misha/Geometry.h"
#include "Misha/AutoDiff/AutoDiff.h"

namespace Samples
{
	struct InvalidSamplerException : public std::exception { const char *what( void ){ return "Bad sampler type"; } };

	template< unsigned int Dim >
	Point< double , Dim > toPoint( AutoDiff::Tensor< UIntPack::Pack< Dim > > t )
	{
		Point< double , Dim > p;
		for( unsigned int i=0 ; i<Dim ; i++ ) p[i] = t[i];
		return p;
	}

	template< unsigned int Dim , typename AngleFunctor /* = std::function< double ( void ) > */ >
	SquareMatrix< double , Dim > RandomRotation( AngleFunctor af )
	{
		SquareMatrix< double , Dim > R = SquareMatrix< double , Dim >::Identity();
		for( unsigned int d=0 ; d<Dim/2 ; d++ )
		{
			double theta = af();
			double ct = cos(theta) , st = sin(theta);
			R( 2*d+0 , 2*d+0 ) =  ct;
			R( 2*d+0 , 2*d+1 ) =  st;
			R( 2*d+1 , 2*d+0 ) = -st;
			R( 2*d+1 , 2*d+1 ) =  ct;
		}
		SquareMatrix< double , Dim > _R = RandomRotationMatrix< double , Dim >();
		return _R.transpose() * R * _R;
	}

	// Curve in 3D
	namespace Curve
	{
		// Functionality for transforming tangent to normals
		Hat::SkewSymmetricMatrix< double , 3 > Dual( Point< double , 3 > t )
		{
			SquareMatrix< double , 3 > n;
			for( unsigned int d=0 ; d<3 ; d++ ) n( (d+1)%3 , (d+2)%3 ) = t[d] , n( (d+2)%3 , (d+1)%3 ) = -t[d];
			return Hat::SkewSymmetricMatrix< double , 3 >( n );
		}

		// Functionality for transforming normal to tangent
		Point< double , 3 > Dual( const Hat::SkewSymmetricMatrix< double , 3 > &n )
		{
			Point< double , 3 > t;
			SquareMatrix< double , 3 > _n =n();
			for( unsigned int d=0 ; d<3 ; d++ ) t[d] = ( _n( (d+1)%3 , (d+2)%3 ) - _n( (d+2)%3 , (d+1)%3 ) ) / 2;
			return t;
		}

		// A virtual class representing a type that takes in a 1D parameter and returns the associated sample (in 3D)
		struct Sampler
		{
			virtual std::pair< Point< double , 3 > , Hat::SkewSymmetricMatrix< double , 3 > > operator()( double s ) const = 0;
		};


		// Returns a sample from two linked circles, one in the xy-plane, the other in the xz-plane
		struct LinkSampler : public Sampler
		{
			std::pair< Point< double , 3 > , Hat::SkewSymmetricMatrix< double , 3 > > operator()( double s ) const
			{
				double radius = 2./3;
				static const double X = 2.*M_PI;
				AutoDiff::Identity< UIntPack::Pack<> > S;
				AutoDiff::Constant< UIntPack::Pack<> , UIntPack::Pack<> > zero;
				auto conditional = [&]( AutoDiff::Tensor< UIntPack::Pack<> > s ){ return (double)s<0.5; };
				auto F1 = AutoDiff::Conditional( conditional , AutoDiff::Cos( X * 2 * S ) * radius - radius/2 , AutoDiff::Cos( X * 2 * (S-0.5) ) * radius + radius/2 );
				auto F2 = AutoDiff::Conditional( conditional , AutoDiff::Sin( X * 2 * S ) * radius , zero );
				auto F3 = AutoDiff::Conditional( conditional , zero , AutoDiff::Sin( X * 2 * (S-0.5) ) * radius );
				auto F = AutoDiff::Concatenation( F1 , F2 , F3 );
				auto dF = F.d();
				return std::make_pair( toPoint( F(s) ) , Dual( toPoint( dF(s) ) ) );
			}
		};

		struct TrefoilAlgebraicSampler : public Sampler
		{
			std::pair< Point< double , 3 > , Hat::SkewSymmetricMatrix< double , 3 > > operator()( double s ) const
			{
                // Samples from a trefoil for which we know an algebraic formulation: https://math.stackexchange.com/questions/148916/trefoil-knot-as-an-algebraic-curve
				static const double X = 2.*M_PI;
				AutoDiff::Identity< UIntPack::Pack<> > S;
				auto Denom = 2. - AutoDiff::Sin( 2. * X * S ) - AutoDiff::Cos( 2. * X * S );
				auto F1 = ( AutoDiff::Cos( 2. * X * S ) - AutoDiff::Sin( 2. * X * S ) ) / Denom;
				auto F2 = AutoDiff::Cos( 3. * X * S ) / Denom;
				auto F3 = AutoDiff::Sin( 3. * X * S ) / Denom;
				auto F = AutoDiff::Concatenation( F1 , F2 , F3 );
				auto dF = F.d();
				return std::make_pair( toPoint( F(s) ) , Dual( toPoint( dF(s) ) ) );
			}
		};

		// Returns a sample from a circle in the xy-plane
		struct CircleSampler : public Sampler
		{
			std::pair< Point< double , 3 > , Hat::SkewSymmetricMatrix< double , 3 > > operator()( double s ) const
			{
				static const double X = 2. * M_PI;
				AutoDiff::Identity< UIntPack::Pack<> > S;
				auto F1 = AutoDiff::Cos( X * S );
				auto F2 = AutoDiff::Sin( X * S );
				auto F = AutoDiff::Concatenation( F1 , AutoDiff::Constant< UIntPack::Pack<> , UIntPack::Pack<> >() , F2 );
				auto dF = F.d();
				return std::make_pair( toPoint( F(s) ) , Dual( toPoint( dF(s) ) ) );
			}
		};

		// Returns a sample from a line-segment aligned with the x-axis
		struct LineSegmentSampler : public Sampler
		{
			std::pair< Point< double , 3 > , Hat::SkewSymmetricMatrix< double , 3 > > operator()( double s ) const
			{
				AutoDiff::Identity< UIntPack::Pack<> > S;
				auto F = AutoDiff::Concatenation( 2. * S - 1. , AutoDiff::Constant< UIntPack::Pack<> , UIntPack::Pack<> >() , AutoDiff::Constant< UIntPack::Pack<> , UIntPack::Pack<> >() );
				auto dF = F.d();
				return std::make_pair( toPoint( F(s) ) , Dual( toPoint( dF(s) ) ) );
			}
		};

		struct BorromeanRingsSampler : public Sampler
		{
			BorromeanRingsSampler( double r=0.25 ) : _r(r)
			{
				_h = 0.5*(sqrt(7)-1)*_r;

				// angles in [0,1]
				_alpha = sin(_h/(2.*_r)) / (2.*M_PI);
				_theta = 1. - 2.*acos(_h/(2.*_r)) / (2.*M_PI);
				_length = 4 * _alpha + 2*_theta;

				double centers[4][3] =
				{
					{   _r+_h ,  0. , 0. }, 
					{      0. ,  _h , 0. }, 
					{-(_r+_h) ,  0. , 0. },
					{      0. , -_h , 0. }
				};
				memcpy( _centers , centers , sizeof( double )*4*3 );
				double borders[4] = { 2.*_alpha , 2.*_alpha+_theta , 4.*_alpha+_theta , 4.*_alpha + 2.*_theta };
				memcpy( _borders , borders , sizeof( double )*4 );
			}
			std::pair< Point< double , 3 > , Hat::SkewSymmetricMatrix< double , 3 > > operator()( double s ) const
			{
				// which of the 3 knots
				int k=0;
				if     ( 1./3<=s && s<2./3 ) k=1 , s -= 1./3.;
				else if( 2./3<=s           ) k=2 , s -= 2./3.;
				AutoDiff::Identity< UIntPack::Pack<> > S;
				auto L = S * 3. * _length;

				Point< double , 3 > p , t;
				double l = _length * s * 3;

				if( l<_borders[0] )
				{
					auto Angle = 0.5 + _alpha - L;
					auto F1 = _centers[0][0] + AutoDiff::Cos( 2. * M_PI * Angle ) * _r;
					auto F2 = _centers[0][1] + AutoDiff::Sin( 2. * M_PI * Angle ) * _r;
					auto dF1 = F1.d();
					auto dF2 = F2.d();

					p[ _orientations[k][0] ] = (double)F1(s);
					p[ _orientations[k][1] ] = (double)F2(s);
					t[ _orientations[k][0] ] = (double)dF1(s);
					t[ _orientations[k][1] ] = (double)dF2(s);
				}
				else if( _borders[0]<=l && l<_borders[1] )
				{
					auto Angle = ( 0.75 + (1.-_theta)/2. - _borders[0] ) + L;
					auto F1 = _centers[1][0] + AutoDiff::Cos( 2. * M_PI * Angle ) * _r;
					auto F2 = _centers[1][1] + AutoDiff::Sin( 2. * M_PI * Angle ) * _r;
					auto dF1 = F1.d();
					auto dF2 = F2.d();

					p[ _orientations[k][0] ] = (double)F1(s);
					p[ _orientations[k][1] ] = (double)F2(s);
					t[ _orientations[k][0] ] = (double)dF1(s);
					t[ _orientations[k][1] ] = (double)dF2(s);
				}
				else if( _borders[1]<=l && l<_borders[2] )
				{
					auto Angle = _alpha - ( L - _borders[1] );
					auto F1 = _centers[2][0] + AutoDiff::Cos( 2. * M_PI * Angle ) * _r;
					auto F2 = _centers[2][1] + AutoDiff::Sin( 2. * M_PI * Angle ) * _r;
					auto dF1 = F1.d();
					auto dF2 = F2.d();

					p[ _orientations[k][0] ] = (double)F1(s);
					p[ _orientations[k][1] ] = (double)F2(s);
					t[ _orientations[k][0] ] = (double)dF1(s);
					t[ _orientations[k][1] ] = (double)dF2(s);
				}
				else // ( _borders[2]<=l && l<_borders[3] )
				{
					auto Angle = ( 0.75 -_theta/2. - _borders[2] ) + L;
					auto F1 = _centers[3][0] + AutoDiff::Cos( 2. * M_PI * Angle ) * _r;
					auto F2 = _centers[3][1] + AutoDiff::Sin( 2. * M_PI * Angle ) * _r;
					auto dF1 = F1.d();
					auto dF2 = F2.d();

					p[ _orientations[k][0] ] = (double)F1(s);
					p[ _orientations[k][1] ] = (double)F2(s);
					t[ _orientations[k][0] ] = (double)dF1(s);
					t[ _orientations[k][1] ] = (double)dF2(s);
				}
				// Assuming that p \in [-_r-_h,_r+_h]
				p /= ( _r + _h );
				t /= ( _r + _h );
				return std::make_pair( p , Dual(t) );
			}
		protected:
			double _r , _h , _alpha , _theta , _borders[4] , _length , _centers[4][3];
			const int _orientations[3][2] = { {0,1} , {1,2} , {2,0} };
		};


		// Returns a sample from a spiral with a prescribed number of rotations and height
		struct SpiralSampler : public Sampler
		{
			SpiralSampler( int r ) : _r(r) {}
			std::pair< Point< double , 3 > , Hat::SkewSymmetricMatrix< double , 3 > > operator()( double s ) const
			{
				static const double X = 2. * M_PI * _r;
				AutoDiff::Identity< UIntPack::Pack<> > S;
				auto F = AutoDiff::Concatenation( AutoDiff::Cos( X * S ) , -1. + S * 2. , AutoDiff::Sin( X * S ) );
				auto dF = F.d();
				return std::make_pair( toPoint( F(s) ) , Dual( toPoint( dF(s) ) ) );
			}
		protected:
			int _r;
		};

		struct TorusKnotSampler : public Sampler
		{
			TorusKnotSampler( int p , int q , double innerRadius=2./3 , double outerRadius=1./3 ) : _p(p) , _q(q) , _ir(innerRadius) , _or(outerRadius)
			{
				unsigned int max = std::min< unsigned int >( std::abs(_p) , std::abs(_q) );
				if( max==0 ) _gcd = 1;
				else
				{
					for( _gcd=max ; _gcd>=1 ; _gcd-- ) if( (_p%_gcd)==0 && (_q%_gcd)==0 ) break;
					_p /= _gcd , _q /= _gcd;
				}
			}
			std::pair< Point< double, 3 > , Hat::SkewSymmetricMatrix< double , 3 > > operator()( double s ) const
			{
				// ref: https://en.wikipedia.org/wiki/Torus_knot
				// (2, -3) should be the Trefoil knot
				// Point: ( r * cos( p * phi ) , r * sin( p * phi ) , -sin( q * phi ) ) , where r = cos( q * phi ) + 2 and phi \in (0, 2\pi) 
				// r' = -q * sin( q * phi )
				// Tangent: (-r * p * sin( p * phi ) + r' * cos( p * phi ) , r * p * cos( p * phi ) + r' * sin( q * phi ) , -q cos( q * phi ) )
				AutoDiff::Identity< UIntPack::Pack<> > S;
				int is = (int)floor(s*_gcd);
				double theta_offset = _q ? ( 2. * M_PI * is ) / ( _gcd * _q ) : 0.;
				auto Theta = 2.*M_PI*( S*(double)_gcd - (double)is );
				auto SinP = AutoDiff::Sin( Theta * (double)_p + theta_offset );
				auto CosP = AutoDiff::Cos( Theta * (double)_p + theta_offset );
				auto SinQ = AutoDiff::Sin( Theta * (double)_q );
				auto CosQ = AutoDiff::Cos( Theta * (double)_q );
				auto R = ( _or * CosQ + _ir );
				auto F = AutoDiff::Concatenation( R * CosP , R * SinP , -SinQ * _or );
				auto dF = F.d();
				return std::make_pair( toPoint( F(s) ) , Dual( toPoint( dF(s) ) ) );
			}
		protected:
			int _p , _q , _gcd;
			double _ir , _or;
		};

		struct BiTorusKnotSampler : public Sampler
		{
			BiTorusKnotSampler( int p1 , int q1 , int p2 , int q2 , double innerRadius1=0.6 , double outerRadius1=0.4 , double innerRadius2=0.6 , double outerRadius2=0.1 )
				: _tk1( p1 , q1 , innerRadius1 , outerRadius1 ) , _tk2( p2 , q2 , innerRadius2 , outerRadius2 ) {}
			std::pair< Point< double, 3 > , Hat::SkewSymmetricMatrix< double , 3 > > operator()( double s ) const
			{
				std::pair< Point< double, 3 > , Hat::SkewSymmetricMatrix< double , 3 > > sample = s<0.5 ? _tk1(s*2.) : _tk2(s*2.-1.);
				sample.second * 2;
				return sample;
			}
		protected:
			TorusKnotSampler _tk1 , _tk2;
		};

		struct SphericalCurveSampler : public Sampler
		{
			SphericalCurveSampler( int type , int n , double l ) : _type(type), _n(n), _l(l) { if (l < 0.) l = 0.; if (l > M_PI * 2) l = 2 * M_PI; }
			std::pair< Point< double, 3 >, Hat::SkewSymmetricMatrix< double, 3 > > operator()(double s) const
			{
				// in the form: ref: https://mathcurve.com/courbes3d.gb/spheric/spheric.shtml
				// x = cos f * cos g
				// y = sin f * cos g
				// z = sin g
				if (_type == 0) {
					// x = cos( (2 * M_PI - l) * cos (s)) * cos ( s / n + l * sin 2(s) )
					// y = cos( (2 * M_PI - l) * cos (s)) * sin ( s / n + l * sin 2(s) )
					// z = sin( (2 * M_PI - l) * cos (s))
					// s \in [0, 2n \pi]
					AutoDiff::Identity< UIntPack::Pack<> > S;
					auto CosTerm1 = AutoDiff::Cos( ( 2. * M_PI - _l ) * AutoDiff::Cos(      S * 2. * (double)_n * M_PI ) );
					auto CosTerm2 = AutoDiff::Cos( S * 2. * M_PI + _l * AutoDiff::Sin( 2. * S * 2. * (double)_n * M_PI ) );
					auto SinTerm1 = AutoDiff::Sin( ( 2. * M_PI - _l ) * AutoDiff::Cos(      S * 2. * (double)_n * M_PI ) );
					auto SinTerm2 = AutoDiff::Sin( S * 2. * M_PI + _l * AutoDiff::Sin( 2. * S * 2. * (double)_n * M_PI ) );
					auto F = AutoDiff::Concatenation( CosTerm1 * CosTerm2 , CosTerm1 * SinTerm2 , SinTerm1 );
					auto dF = F.d();
					return std::make_pair( toPoint( F(s) ) , Dual( toPoint( dF(s) )  ) );
				}
				else
				{
					// x = cos(l * cos ( s / n ) ) * cos( s )
					// y = cos(l * cos ( s / n ) ) * sin( s )
					// z = sin(l * cos ( s / n ) ) 
					// s \in [0, 2\pi]
					AutoDiff::Identity< UIntPack::Pack<> > S;
					auto CosTerm1 = AutoDiff::Cos( _l * AutoDiff::Cos( S * 2. * (double)_n * M_PI ) );
					auto CosTerm2 = AutoDiff::Cos( S * 2. * M_PI );
					auto SinTerm1 = AutoDiff::Sin( _l * AutoDiff::Cos( S * 2. * (double)_n * M_PI ) );
					auto SinTerm2 = AutoDiff::Sin( S * 2. * M_PI );
					auto F = AutoDiff::Concatenation( CosTerm1 * CosTerm2 , CosTerm1 * SinTerm2 , SinTerm1 );
					auto dF = F.d();
					return std::make_pair( toPoint( F(s) ) , Dual( toPoint( dF(s) ) ) );
				}
			}
		protected:
			int _type;
			int _n;
			double _l;
		};

		enum Type
		{
			LINE_SEGMENT = 0,
			CIRCLE ,
			LINK ,
			SPIRAL ,
			TORUS_KNOT ,
			BI_TORUS_KNOT ,
			SPHERICAL_CURVE ,
			BORROMEAN_RINGS ,
			TREFOIL_ALGEBRAIC ,
			COUNT
		};
		static const unsigned int TypeCount = Type::COUNT;
		static const std::string Descriptions[] =
		{ 
			std::string( "line_segment" ) , 
			std::string( "circle" ) , 
			std::string( "link" ) , 
			std::string( "spiral:<r>" ) , 
			std::string( "torus_knot:<p>:<q>" ) ,
			std::string( "bi_torus_knot:<op>:<oq>:<ip>:<iq>" ) ,
			std::string( "spherical:<type>:<n>:<l>" ) ,
			std::string( "borromean_rings" ),
			std::string( "trefoil_algebraic" )
		};

		std::string Header( Type type )
		{
			std::istringstream ss( Descriptions[type] );
			std::string header;
			if( !std::getline( ss , header , ':' ) ) ERROR_OUT( "Could not read header: " , Descriptions[type] );
			return header;
		};

		std::vector< std::pair< Point< double , 3 > , Hat::SkewSymmetricMatrix< double , 3 > > > GetSamples( unsigned int count , std::string str , bool random=false )
		{
			std::vector< std::pair< Point< double , 3 > , Hat::SkewSymmetricMatrix< double , 3 > > > samples(count);

			Sampler *sampler = NULL;

			std::istringstream ss( str );
			std::string header;
			if( !std::getline( ss , header , ':' ) ) ERROR_OUT( "Could not read header: " , str );

			if     ( header==Header(LINE_SEGMENT)      ) sampler = new LineSegmentSampler();
			else if( header==Header(CIRCLE)            ) sampler = new CircleSampler();
			else if( header==Header(LINK)              ) sampler = new LinkSampler();
			else if( header==Header(BORROMEAN_RINGS )  ) sampler = new BorromeanRingsSampler();
			else if( header==Header(TREFOIL_ALGEBRAIC) ) sampler = new TrefoilAlgebraicSampler();
			else if( header==Header(SPIRAL) )
			{
				std::string r;
				if( !std::getline( ss , r , ':' ) ) ERROR_OUT( "Could not read spiral rotations: " , str );
				sampler = new SpiralSampler( std::stoi( r ) );
			}
			else if( header==Header(TORUS_KNOT) )
			{
				std::string p , q;
				if( !std::getline( ss , p , ':' ) ) ERROR_OUT( "Could not read torus knot p: " , str );
				if( !std::getline( ss , q , ':' ) ) ERROR_OUT( "Could not read torus knot q: " , str );
				sampler = new TorusKnotSampler( std::stoi(p) , std::stoi(q) );
			}
			else if( header==Header(BI_TORUS_KNOT) )
			{
				std::string op , oq , ip , iq;
				if( !std::getline( ss , op , ':' ) ) ERROR_OUT( "Could not read bi torus knot op: " , str );
				if( !std::getline( ss , oq , ':' ) ) ERROR_OUT( "Could not read bi torus knot oq: " , str );
				if( !std::getline( ss , ip , ':' ) ) ERROR_OUT( "Could not read bi torus knot ip: " , str );
				if( !std::getline( ss , iq , ':' ) ) ERROR_OUT( "Could not read bi torus knot iq: " , str );
				sampler = new BiTorusKnotSampler( std::stoi(op) , std::stoi(oq) , std::stoi(ip) , std::stoi(iq) );
			}
			else if( header==Header(SPHERICAL_CURVE) )
			{
				std::string type , n , l;
				if( !std::getline( ss , type , ':' ) ) ERROR_OUT( "Could not read spherical curve type: " , str );
				if( !std::getline( ss , n    , ':' ) ) ERROR_OUT( "Could not read spherical curve n: " , str );
				if( !std::getline( ss , l    , ':' ) ) ERROR_OUT( "Could not read spherical curve l: " , str );
				sampler = new SphericalCurveSampler( std::stoi(type) , std::stoi(n) , std::stod(l) );
			}
			if( !sampler ) throw InvalidSamplerException();


			if( random )
			{
				std::default_random_engine generator;
				std::uniform_real_distribution< double > distribution( 0. , 1. );
				for( unsigned int i=0 ; i<count ; i++ ) samples[i] = sampler->operator()( distribution( generator ) );
			}
			else for( unsigned int i=0 ; i<count ; i++ ) samples[i] = sampler->operator()( (i+0.5)/count );
			return samples;
		}
	}

	// Surface in 4D
	namespace Surface
	{
		// Functionality for transforming tangents to normals
		Hat::SkewSymmetricMatrix< double , 4 > Dual( Hat::SkewSymmetricMatrix< double , 4 > t )
		{
			// In 4D, ref: https://en.wikipedia.org/wiki/Hodge_star_operator <https://en.wikipedia.org/wiki/Hodge_star_operator> 

			SquareMatrix< double , 4 > n , _t = t();
			n( 0 , 0 ) = 0           , n( 1 , 0 ) = _t( 2 , 3 ) , n( 2 , 0 ) = _t( 3 , 1 ) , n( 3 , 0 ) = _t( 1 , 2 );
			n( 0 , 1 ) = _t( 3 , 2 ) , n( 1 , 1 ) = 0           , n( 2 , 1 ) = _t( 0 , 3 ) , n( 3 , 1 ) = _t( 2 , 0 );
			n( 0 , 2 ) = _t( 1 , 3 ) , n( 1 , 2 ) = _t( 3 , 0 ) , n( 2 , 2 ) = 0           , n( 3 , 2 ) = _t( 0 , 1 );
			n( 0 , 3 ) = _t( 2 , 1 ) , n( 1 , 3 ) = _t( 0 , 2 ) , n( 2 , 3 ) = _t( 1 , 0 ) , n( 3 , 3 ) = 0          ;
			return Hat::SkewSymmetricMatrix< double , 4 >( n );
		}

		// A virtual class representing a type that takes in a 2D parameter and returns the associated sample (in 4D)
		struct Sampler
		{
			virtual std::pair< Point< double , 4 > , Hat::SkewSymmetricMatrix< double , 4 > > operator()( double s , double t ) const = 0;
		};

		struct CliffordTorusSampler : public Sampler
		{
			std::pair< Point< double , 4 >,Hat::SkewSymmetricMatrix< double , 4 > > operator()( double s , double t ) const
			{
				AutoDiff::Identity< UIntPack::Pack<> > S , T;
				static const double X = 2. * M_PI;
				auto F1 = AutoDiff::Cos( X*S ) / sqrt( 2. );
				auto F2 = AutoDiff::Sin( X*S ) / sqrt( 2. );
				auto F3 = AutoDiff::Cos( X*T ) / sqrt( 2. );
				auto F4 = AutoDiff::Sin( X*T ) / sqrt( 2. );
				auto dF1 = F1.d();
				auto dF2 = F2.d();
				auto dF3 = F3.d();
				auto dF4 = F4.d();
				return std::make_pair
				(
					Point< double , 4 >( F1(s) , F2(s) , F3(t) , F4(t) ) ,
					Hat::Wedge( Point< double , 4 >( dF1(s) , dF2(s) , 0 , 0 ) , Point< double , 4 >( 0 , 0 , dF3(t) , dF4(t) ) )
				);
			}
		};

		struct MobiusTubeSampler : public Sampler
		{
			std::pair< Point< double , 4 > , Hat::SkewSymmetricMatrix< double , 4 > > operator()( double s , double t ) const
			{
				// ref: https://en.wikipedia.org/wiki/Klein_bottle
				double theta = 2. * M_PI;
				double phi = 2. * M_PI;
				double R = 3./5;
				double r = 2./5;
				Point< double , 4 > tangent1( 
					-r * theta * sin( theta * s ) * cos( phi * t ) , 
					-r * theta * sin( theta * s ) * sin( phi * t ) , 
					r * theta * cos( theta * s ) * cos( phi * t * 0.5 ) , 
					r * theta * cos( theta * s ) * sin( phi * t * 0.5 ) 
				);
				Point< double , 4 > tangent2( 
					-( R + r * cos( theta * s ) ) * phi * sin( phi * t ) , 
					( R + r * cos( theta * s ) ) * phi * cos( phi * t ) , 
					-r * sin( theta * s ) * phi * 0.5 * sin( phi * t * 0.5 ) , 
					r * sin( theta * s ) * phi * 0.5 * cos( phi * t * 0.5 ) 
				);
				tangent1 /= Point< double, 4 >::Length( tangent1 );
				tangent2 /= Point< double, 4 >::Length( tangent2 );
				return std::make_pair
				(
					Point< double , 4 >( 
						( R + r * cos( theta * s ) ) * cos( phi * t ) , 
						( R + r * cos( theta * s ) ) * sin( phi * t ) , 
						r * sin( theta * s ) * cos( phi * t * 0.5 ) , 
						r * sin( theta * s ) * sin( phi * t * 0.5 ) 
						) ,
					Dual( Hat::Wedge( tangent1 , tangent2 ) )
				);
			}
		};

		struct HelixSurfaceSampler : public Sampler
		{
			HelixSurfaceSampler( int type , double d ) : _type(type) , _d(d) { if (_d <= 0. || _d >= 1.) throw InvalidSamplerException(); }
			std::pair< Point< double , 4 > , Hat::SkewSymmetricMatrix< double , 4 > > operator()( double s , double v ) const
			{
				// ref: https://arxiv.org/pdf/1705.10090.pdf
				SquareMatrix< double , 4 > J1, J2, J3;
				J1 *= 0, J2 *= 0, J3 *= 0;
				J1( 1 , 0 ) = -1; J1( 0 , 1 ) =  1; J1( 3 , 2 ) = -1; J1( 2 , 3 ) =  1;
				J2( 3 , 0 ) = -1; J2( 2 , 1 ) = -1; J2( 1 , 2 ) =  1; J2( 0 , 3 ) =  1;
				J3( 2 , 0 ) = -1; J3( 3 , 1 ) =  1; J3( 0 , 2 ) =  1; J3( 1 , 3 ) = -1;
				double _s(0), _v(0), xi0(0), xi1(0), xi2(0), xi3(0), nu(0), epsilon(0);
				double ds(0), dv(0), dxi0(0), dxi1(0), dxi2(0), dxi3(0);
				double lambda = pow(2 * epsilon * abs(nu) * _d / (1 - _d * _d), 2) - nu * nu;
				double invD = 1. / sqrt(1 + _d * _d);

				// Note: Q is a function in v while b is a function in s
				// parameters
				switch (_type) {
				case 0:
					_s = s * 8. * M_PI - 4. * M_PI; // s \in [-4\pi , 4\pi]
					ds = 8. * M_PI;

					_v = v * 4. * M_PI - 2. * M_PI; // v \in [-2\pi , 2\pi]
					dv = 4. * M_PI;
					xi0 = M_PI * 0.5;
					xi1 = M_PI * 0.25;
					xi2 = xi3 = _v;
					dxi2 = dxi3 = dv; // w.r.t. v
					nu = 4.; 
					epsilon = 2.;
					break;
				case 1:
					_s = s * 4. * M_PI - 2. * M_PI; // s \in [-2\pi , 2\pi]
					ds = 4. * M_PI;

					_v = v * 4. - 2.; // v \in [-2 , 2]
					dv = 4.;
					xi0 = M_PI * 0.5;
					xi1 = M_PI * 0.25;
					xi2 = xi3 = exp( _v );
					dxi2 = dxi3 = exp( _v ) * dv; // w.r.t. v
					nu = 2.;
					epsilon = 1.;
					break;
				case 2:
					_s = s * 8. * M_PI - 4. * M_PI; // s \in [-4\pi , 4\pi]
					ds = 8. * M_PI;

					_v = v * 4. * M_PI - 2. * M_PI; // v \in [-2\pi , 2\pi]
					dv = 4. * M_PI;

					xi0 = M_PI * 0.5;
					xi1 = invD;
					nu = sqrt(5.);
					xi2 = sqrt(lambda + nu * nu) / _d * _v;
					xi3 = sqrt(lambda + nu * nu) * _d * _v;
					dxi2 = sqrt(lambda + nu * nu) / _d * dv;
					dxi3 = sqrt(lambda + nu * nu) * _d * dv;
					epsilon = 1.;
					break;
				case 3:
					_s = s * 8. * M_PI - 4. * M_PI; // s \in [-4\pi , 4\pi]
					ds = 8. * M_PI;

					_v = v * 4. * M_PI - 2. * M_PI; // v \in [-2\pi , 2\pi]
					dv = 4. * M_PI;

					xi0 = M_PI * 0.5;
					xi1 = invD;
					nu = 2;
					xi2 = sqrt(lambda + nu * nu) / _d * _v;
					xi3 = sqrt(lambda + nu * nu) * _d * _v;
					dxi2 = sqrt(lambda + nu * nu) / _d * dv;
					dxi3 = sqrt(lambda + nu * nu) * _d * dv;
					epsilon = 1.;
					break;
				default: throw InvalidSamplerException();
				}

				double cos0 = cos( xi0 ), sin0 = sin( xi0 );
				double cos1 = cos( xi1 ), sin1 = sin( xi1 );
				double cos2 = cos( xi2 ), sin2 = sin( xi2 );
				double cos3 = cos( xi3 ), sin3 = sin( xi3 );
				// Taking derivative w.r.t. v
				double dcos0 = -sin0 * dxi0, dsin0 = cos0 * dxi0;
				double dcos1 = -sin1 * dxi1, dsin1 = cos1 * dxi1;
				double dcos2 = -sin2 * dxi2, dsin2 = cos2 * dxi2;
				double dcos3 = -sin3 * dxi3, dsin3 = cos3 * dxi3;

				Point< double , 4 > r1( cos1 * cos2 , -cos1 * sin2 , sin1 * cos3 , -sin1 * sin3 );
				Point< double , 4 > r2 = J1 * r1;
				Point< double , 4 > r3 = J2 * r1 * cos0 + J3 * r1 * sin0;
				Point< double , 4 > r4 = J2 * r1 * sin0 - J3 * r1 * cos0;
				// Taking derivative w.r.t. v
				Point< double , 4 > dr1( dcos1 * cos2 + cos1 * dcos2 , -dcos1 * sin2 + -cos1 * dsin2 , dsin1 * cos3 + sin1 * dcos3 , -dsin1 * sin3 + -sin1 * dsin3 );
				Point< double , 4 > dr2 = J1 * dr1;
				Point< double , 4 > dr3 = J2 * ( dr1 * cos0 + r1 * dcos0 ) + J3 * ( dr1 * sin0 + r1 * dsin0 );
				Point< double , 4 > dr4 = J2 * ( dr1 * sin0 + r1 * dsin0 ) - J3 * ( dr1 * cos0 + r1 * dcos0 );

				SquareMatrix< double , 4 > Q;
				Q( 0 , 0 ) = r1[0]; Q( 1 , 0 ) = r1[1]; Q( 2 , 0 ) = r1[2]; Q( 3 , 0 ) = r1[3];
				Q( 0 , 1 ) = r2[0]; Q( 1 , 1 ) = r2[1]; Q( 2 , 1 ) = r2[2]; Q( 3 , 1 ) = r2[3];
				Q( 0 , 2 ) = r3[0]; Q( 1 , 2 ) = r3[1]; Q( 2 , 2 ) = r3[2]; Q( 3 , 2 ) = r3[3];
				Q( 0 , 3 ) = r4[0]; Q( 1 , 3 ) = r4[1]; Q( 2 , 3 ) = r4[2]; Q( 3 , 3 ) = r4[3];
				// Taking derivative w.r.t. v
				SquareMatrix< double , 4 > dQ;
				dQ( 0 , 0 ) = dr1[0]; dQ( 1 , 0 ) = dr1[1]; dQ( 2 , 0 ) = dr1[2]; dQ( 3 , 0 ) = dr1[3];
				dQ( 0 , 1 ) = dr2[0]; dQ( 1 , 1 ) = dr2[1]; dQ( 2 , 1 ) = dr2[2]; dQ( 3 , 1 ) = dr2[3];
				dQ( 0 , 2 ) = dr3[0]; dQ( 1 , 2 ) = dr3[1]; dQ( 2 , 2 ) = dr3[2]; dQ( 3 , 2 ) = dr3[3];
				dQ( 0 , 3 ) = dr4[0]; dQ( 1 , 3 ) = dr4[1]; dQ( 2 , 3 ) = dr4[2]; dQ( 3 , 3 ) = dr4[3];

				Point< double , 4 > b( _d * cos( _s / _d ) , lambda * _d * sin( _s / _d ) , cos( _d * _s ) , lambda * sin( _d * _s ) );
				b *= invD;
				// Taking derivative w.r.t. s
				Point< double , 4 > db( -sin( _s / _d ) , lambda * cos( _s / _d ) , -_d * sin( _d * _s ) , lambda * _d * cos( _d * _s ) ); 
				db *= ds * invD;

				Point< double , 4 > tangent1 = Q * db; // w.r.t. s
				Point< double , 4 > tangent2 = dQ * b; // w.r.t. t
				tangent1 /= Point< double, 4 >::Length( tangent1 );
				tangent2 /= Point< double, 4 >::Length( tangent2 );
				return std::make_pair
				(
					Q * b , 
					Dual( Hat::Wedge<4>( tangent1 , tangent2 ) )
				);
			}
		protected:
			int _type;
			double _d;
		};

		struct HopfTorusSampler : public Sampler
		{
			HopfTorusSampler( int type , int n , double l ) : curveSampler( type , n , l ) {}
			std::pair< Point< double, 4 >, Hat::SkewSymmetricMatrix< double, 4 > > operator()( double s , double t ) const
			{
				// use spherical curve and the Hopf fibration to map it to the 3-sphere
				// ref: https://en.wikipedia.org/wiki/Hopf_fibration
				// given (a, b, c) on the 2-sphere and t \in [0, 2\pi]
				// x = ( 1 + c ) * cos( t )
				// y = a * sin( t ) - b * cos( t )
				// z = a * cos( t ) + b * sin( t )
				// w = ( 1 + c ) * sin( t )
				// normalized by factor = 1 / sqrt(2(1 + c))
				const auto& [cp, cn] = curveSampler(s);
				const auto ct = Curve::Dual(cn);
				double _t = t * 2 * M_PI;
				double dt = 2 * M_PI;
				double cost = cos(_t), sint = sin(_t);
				double dcosdt = -sin(_t), dsindt = cos(_t);
				double a = cp[0], b = cp[1], c = cp[2];
				double dads = ct[0], dbds = ct[1], dcds = ct[2];
				double factor = 1. / sqrt(2 * (1 + c));
				double dfactords = -dcds / (pow(2 * (1 + c), 1.5));
				double x = ( 1 + c ) * cost * factor;
				double y = (a * sint - b * cost) * factor;
				double z = (a * cost + b * sint) * factor;
				double w = (1 + c) * sint * factor;
				double dxdt = (1 + c) * dcosdt * factor;
				double dydt = (a * dsindt - b * dcosdt) * factor;
				double dzdt = (a * dcosdt + b * dsindt) * factor;
				double dwdt = (1 + c) * dsindt * factor;
				double dxds = dcds * cost * factor + (1 + c) * cost * dfactords;
				double dyds = (dads * sint - dbds * cost) * factor + (a * sint - b * cost) * dfactords;
				double dzds = (dads * cost + dbds * sint) * factor + (a * cost + b * sint) * dfactords;
				double dwds = dcds * sint * factor + (1 + c) * sint * dfactords;
				return std::make_pair
				(
					Point< double , 4 >( x , y , z , w ) ,
					Dual( Hat::Wedge<4>( Point< double , 4 >( dxdt , dydt , dzdt , dwdt ) , Point< double , 4 >( dxds , dyds , dzds , dwds ) ) )
				);
			}
		protected:
			Curve::SphericalCurveSampler curveSampler;
		};

		struct KYSampler : public Sampler
		{
			KYSampler(int m, double xi) : m(m), xi(xi) { 
				tau = 1. / m * exp(-m * m / 4 / M_PI + xi) ;
			}
			std::pair< Point< double , 4 >,Hat::SkewSymmetricMatrix< double , 4 > > operator()( double s , double t ) const
			{
				double x = 2. * M_PI;
				s = sqrt(2) * s;
				t = sqrt(2) * t;
				double val = sqrt(s * s + t * t);
				double u;
				if (val < tau) {
					if( val>1./m ) throw std::out_of_range( "out of range" );
					u = -tau * acosh(val / tau);
				}
				else {
					if ((std::max)(s, t) > M_PI / sqrt(2) / m) throw std::out_of_range( "out of range" );
					double psi = m * val;
					if (psi < 1) psi = 1;
					if (psi > 2) psi = 0;
					u = -tau * psi * acosh(val / tau) + tau * (1 - psi) * acosh(1. / m / tau);
				}

				return std::make_pair
				(
					Point< double , 4 >( sin(u + M_PI / 4) * cos(x*s) , sin(u + M_PI / 4) * sin(x*s) , cos(u + M_PI / 4) * cos(x*t) , cos(u + M_PI / 4) * sin(x*t) ),
					Hat::Wedge( Point< double , 4 >(sin(u + M_PI / 4) * cos(x*s) , sin(u + M_PI / 4) * sin(x*s) , 0 , 0 ) , Point< double , 4 >( 0 , 0 , cos(u + M_PI / 4) * cos(x*t) , cos(u + M_PI / 4) * sin(x*t) ) )
				);
			}
		protected:
			int m;
			double xi;
			double tau;
		};

		enum Type
		{
			CLIFFORD_TORUS = 0,
			MOBIUS_TUBE ,
			HELIX_SURFACE ,		// Helix surface in 3-sphere: ref: https://arxiv.org/pdf/1705.10090.pdf 
			LAWSON_SURFACE ,	// Lawson minmal surfaces in 3-sphere ref: https://math.jhu.edu/~js/Math748/lawson.s3.pdf, we will use \xi(m,k)
			HOPF_TORUS ,		// using curve on 2-sphere and the Hopf fibration
			KY ,				// Kapouleas-Yang minimal surfaces in S3
			COUNT
		};
		static const unsigned int TypeCount = Type::COUNT;
		static const std::string Descriptions[] =
		{ 
			std::string( "clifford_torus" ) , 
			std::string( "mobius_tube" ) , 
			std::string( "helix_surface:<type>:<d>" ) ,
			std::string( "lawson_surface:<m>:<k>" ),
			std::string( "hopf_torus:<n>:<l>" ) ,
			std::string( "ky:<m>:<u>" )
		};

		std::string Header( Type type )
		{
			std::istringstream ss( Descriptions[type] );
			std::string header;
			if( !std::getline( ss , header , ':' ) ) ERROR_OUT( "Could not read header: " , Descriptions[type] );
			return header;
		};

		std::vector< std::pair< Point< double , 4 > , Hat::SkewSymmetricMatrix< double , 4 > > > GetSamples( unsigned int count , std::string str , bool random=false )
		{
			std::vector< std::pair< Point< double , 4 > , Hat::SkewSymmetricMatrix< double , 4 > > > samples(count*count);

			Sampler *sampler = NULL;

			std::istringstream ss( str );
			std::string header;
			if( !std::getline( ss , header , ':' ) ) ERROR_OUT( "Could not read header: " , str );
			if     ( header==Header(CLIFFORD_TORUS) ) sampler = new CliffordTorusSampler();
			else if( header==Header(MOBIUS_TUBE)    ) sampler = new MobiusTubeSampler();
			else if( header==Header(HELIX_SURFACE)  ) 
			{
				std::string type, d;
				if ( !std::getline(ss, type, ':') ) ERROR_OUT("Could not read helix sufrace type: ", str);
				if ( !std::getline(ss,    d, ':') ) ERROR_OUT("Could not read helix sufrace param: ", str);
				sampler = new HelixSurfaceSampler( std::stoi(type) , std::stod(d) );
			}
			else if (header == Header(LAWSON_SURFACE))
			{
				std::string m, k;
				if ( !std::getline(ss, m, ':') ) ERROR_OUT("Could not read Lawson sufrace m: ", str);
				if ( !std::getline(ss, k, ':') ) ERROR_OUT("Could not read Lawson sufrace k: ", str);
				unsigned int _m = std::stoi(m), _k = std::stoi(k);

				// ref: https://math.jhu.edu/~js/Math748/lawson.s3.pdf
				// sample the fundamental polygon
				double phi0 = 0.;
				double phi1 = phi0 + M_PI / (_k + 1);
				Point< double, 4> P1(0, 0, cos(phi0), sin(phi0));
				Point< double, 4> P2(0, 0, cos(phi1), sin(phi1));

				double theta0 = 0.;
				double theta1 = theta0 + M_PI / (_m + 1);
				Point< double, 4> Q1(cos(theta0), sin(theta0), 0, 0);
				Point< double, 4> Q2(cos(theta1), sin(theta1), 0, 0);

				std::vector< std::pair< Point< double, 4 >, Hat::SkewSymmetricMatrix< double, 4 > > > patch(count * count);
				// the polygon is P1 Q1 P2 Q2;
				// samples P1 Q1 Q2


			}
			else if( header==Header(HOPF_TORUS) )
			{
				std::string n, l;
				if( !std::getline( ss , n , ':' ) ) ERROR_OUT( "Could not read Hopf torus nodes: " , str );
				if( !std::getline( ss , l , ':' ) ) ERROR_OUT( "Could not read Hopf torus amplitude: " , str );
				sampler = new HopfTorusSampler( 1 , std::stoi(n) , std::stod(l) );
			}
			else if( header==Header(KY)  ) 
			{
				std::string m, xi;
				if ( !std::getline(ss, m, ':') ) ERROR_OUT("Could not read ky surface m: ", str);
				if ( !std::getline(ss, xi, ':') ) ERROR_OUT("Could not read ky surface xi: ", str);
				sampler = new KYSampler( std::stoi(m) , std::stod(xi) );
				samples.clear();
				samples.reserve(count * count);
				for (unsigned int i = 0; i < count; i++) for (unsigned int j = 0; j < count; j++) {
					try {
						auto s = sampler->operator()((i + 0.5) / count, (j + 0.5) / count);
						samples.push_back(s);
					}
					catch (std::exception& ) {}
				}
				std::cout << "# of samples: " << samples.size() << std::endl;
				delete sampler;
				return samples;
			}
			if( !sampler ) throw InvalidSamplerException();


			if( random )
			{
				std::default_random_engine generator;
				std::uniform_real_distribution< double > distribution( 0. , 1. );
				for( unsigned int i=0 ; i<count*count ; i++ ) samples[i] = sampler->operator()( distribution( generator ) , distribution( generator ) );
			}
			else for( unsigned int i=0 ; i<count ; i++ ) for( unsigned int j=0 ; j<count ; j++ ) samples[i*count+j] = sampler->operator()( (i+0.5)/count , (j+0.5)/count );
			return samples;
		}
	}
}
#endif // CURVE_SAMPLES_INCLUDED
