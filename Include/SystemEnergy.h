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

#ifndef SYSTEM_ENERGY_INCLUDED
#define SYSTEM_ENERGY_INCLUDED

#include "Misha/Polynomial.h"
#include "Hat.h"

namespace MishaK
{
	// A structure describing the system energy
	template< unsigned int Dim , typename _Indexer /* = Hat::BaseIndexer< Dim > */ >
	struct SystemEnergy
	{
		using Indexer = _Indexer;

		struct QuadraticApproximation
		{
			Eigen::SparseMatrix< double > q;
			Eigen::VectorXd l;
			double c;

			QuadraticApproximation( void ) : c(0){}
			QuadraticApproximation( const Eigen::SparseMatrix< double > &q , const Eigen::VectorXd &l , double c ) : q(q) , l(l) , c(c){}
			QuadraticApproximation( Eigen::SparseMatrix< double > &&q , const Eigen::VectorXd &&l , double c ) : q( std::move(q) ) , l( std::move(l) ) , c(c){}

			double operator()( const Eigen::VectorXd &x ) const { return ( q * x ).dot( x ) + l.dot(x) + c; }
		};

		Indexer indexer;
		Hat::ScalarFunctions< Dim > scalars;					// The scalar system
		Hat::AlternatingProductFunctions< Dim > wedges;			// The wedge system

		SystemEnergy( const Indexer &indexer , unsigned int r ) : indexer(indexer) , wedges(r) , scalars(r){}

		// Computes the energy components given by the alternating-field fitting and regularization energy terms
		double operator()( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;

		// Given the current estimate of the solution (x,y) and the two directions (_x,_y), computes the value of s minimizing E(x+s*_x,y+s*_y)
		double stepSize( const Eigen::VectorXd &x , const Eigen::VectorXd &y , const Eigen::VectorXd &_x , const Eigen::VectorXd &_y ) const;

		// Given the current estimate of the solution (x,y) and the two directions (_x,_y), computes the values of (s,t) minimizing P(s,t) = E(x+s*_x,y+t*_y) using Newton iterations
		Point< double , 2 > newtonUpdate( const Eigen::VectorXd &x , const Eigen::VectorXd &y , const Eigen::VectorXd &dx , const Eigen::VectorXd &dy , unsigned int steps ) const;

		// Computes the differential of the energy with respect to the two function's coefficients
		virtual std::pair< Eigen::VectorXd , Eigen::VectorXd > d( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;

		// Given the current estimate of the solution (x,y) computes the quadratic polynomial P(s,t) = E(x+s*e[idx],y+t*e[idx])
		virtual Polynomial::Polynomial2D< 2 > quadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int idx , bool setConstantTerm , unsigned int thread ) const;

		// Given the current estimate of the solution (x,y) and the two directions (_x,_y), computes the bi-quadratic polynomial P(s,t) = E(x+s*_x,y+t*_y)
		virtual Polynomial::Polynomial2D< 4 > biQuadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , const Eigen::VectorXd &_x , const Eigen::VectorXd &_y ) const;

		// Given the current estimate of the solution (x,y), gives the quadratic and linear approximations to the energy in x
		virtual QuadraticApproximation quadraticApproximation1( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;

		// Given the current estimate of the solution (x,y), gives the quadratic and linear approximations to the energy in y
		virtual QuadraticApproximation quadraticApproximation2( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;

		//////////////////////////
		// Pure virtual methods //
		//////////////////////////

		// Separately computes the energy components given by the alternating-field fitting (first) and regularization (second) energy terms
		virtual std::pair< double , double > energies( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const = 0;

		// Returns the square norm of the target alternating product
		virtual double squareNorm( void ) const = 0;

		// Returns the square norm of the wedge
		virtual double wedgeSquareNorm( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const = 0;
		virtual double wedgeSquareNorm( const Eigen::VectorXd &w ) const = 0;

		// Returns the prolongation of the scalar function
		virtual Eigen::VectorXd prolongation( const Eigen::VectorXd &coarse ) const = 0;

		////////////////
		// Unit tests //
		////////////////
		// Test everything
		void test( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int numTests , double eps=1e-4 ) const;

		// Test the computation of the differential
		void testD( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int numTests , double eps=1e-4 ) const;

		// Test the computation of the bi-quadratic polynomial
		void testBiQuadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int numTests , double eps=1e-4 ) const;

		// Test the quadratic approximation in the first variable
		void testQuadraticApproximation1( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int numTests , double eps=1e-4 ) const;

		// Test the quadratic approximation in the second variable
		void testQuadraticApproximation2( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int numTests , double eps=1e-4 ) const;

		// Test the quadratic approximation in the paired variables
		void testQuadraticApproximation( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int numTests ) const;

	protected:
		SystemEnergy( void ) : SystemEnergy(1){}
	};

	template< unsigned int Dim , typename Indexer /* = Hat::BaseIndexer< Dim > */ >
	struct SingleCycleCascadicSystemEnergy : public SystemEnergy< Dim , Indexer >
	{
		using SystemEnergy< Dim , Indexer >::scalars;
		using SystemEnergy< Dim , Indexer >::wedges;
		using SystemEnergy< Dim , Indexer >::indexer;

		// Construct the system energy
		template< typename SkewField /* = std::function< SkewSymmetriMatrix< double , Dim > ( size_t ) > */ >
		SingleCycleCascadicSystemEnergy( const Indexer &indexer , unsigned int r , SkewField && skew , double sWeight );
		template< typename SkewField /* = std::function< SkewSymmetriMatrix< double , Dim > ( size_t ) > */ >
		SingleCycleCascadicSystemEnergy( const Indexer &indexer , unsigned int r , SkewField && skew , double sWeight , Eigen::SparseMatrix< double > R );

		// Restrict from the finer resolution
		SingleCycleCascadicSystemEnergy restrict( const Indexer &coarseProlongationIndexer ) const;

		// Update the system given the finer solution
		void update( const SingleCycleCascadicSystemEnergy &finer , const Eigen::VectorXd &x , const Eigen::VectorXd &y );

		// Separately computes the energy components given by the alternating-field fitting (first) and regularization (second) energy terms
		std::pair< double , double > energies( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;

		// Given the current estimate of the solution (x,y) computes the quadratic polynomial P(s,t) = E(x+s*e[idx],y+t*e[idx])
		Polynomial::Polynomial2D< 2 > quadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int idx , bool setConstantTerm , unsigned int thread ) const;

		// Given the current estimate of the solution (x,y) and the two directions (_x,_y), computes the bi-quadratic polynomial P(s,t) = E(x+s*_x,y+t*_y)
		Polynomial::Polynomial2D< 4 > biQuadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , const Eigen::VectorXd &_x , const Eigen::VectorXd &_y ) const;

		// Returns the square norm of the target alternating product
		double squareNorm( void ) const { return _c; }

		// Returns the square norm of the wedge
		double wedgeSquareNorm( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const { return wedges( indexer , _wMStencil , x , y , x , y ); }

		double wedgeSquareNorm( const Eigen::VectorXd &w ) const { return wedges( _wMStencil , w , w ); }

		// Returns the prolongation of the scalar function
		Eigen::VectorXd prolongation( const Eigen::VectorXd &coarse ) const { return _sP * coarse; }

		// Returns the scalar prolongation matrix
		const Eigen::SparseMatrix< double > &prolongation( void ) const { return _sP; }

	protected:
		SingleCycleCascadicSystemEnergy( const Indexer &indexer , unsigned int r );

		typename Hat::AlternatingProductFunctions< Dim >::template FullIntegrationStencil< double > _wMStencil;	// The full wedge mass matrix stencil
		typename Hat::ScalarFunctions< Dim >::template FullIntegrationStencil< double , 0 >  _sSStencil;	// The full scalar stiffness matrix stencil
		Eigen::SparseMatrix< double > _sP;								// The scalar prolongation matrix (from the coarse resolution into this one)
		Eigen::SparseMatrix< double > _R;								// The symmetric matrix defining the regularization energy
		Eigen::SparseMatrix< double > _B;								// The integral of the target alternating-form field against the alternating-form basis
		double _c;														// The constant terms coming from fitting
		double _sWeight;
	};

	template< unsigned int Dim , typename Indexer /* = Hat::BaseIndexer< Dim > */ >
	struct HierarchicalSystemEnergy : public SystemEnergy< Dim , Indexer >
	{
		using SystemEnergy< Dim , Indexer >::scalars;
		using SystemEnergy< Dim , Indexer >::wedges;
		using SystemEnergy< Dim , Indexer >::indexer;

		// Construct the system energy
		HierarchicalSystemEnergy( const Indexer &indexer , unsigned int r , ConstPointer( Hat::SkewSymmetricMatrix< double , Dim > ) skew , double sWeight );
		HierarchicalSystemEnergy( const Indexer &indexer , unsigned int r , ConstPointer( Hat::SkewSymmetricMatrix< double , Dim > ) skew , double sWeight , Eigen::SparseMatrix< double > R );

		// Construct the system energy by restricting from the finer resolution
		HierarchicalSystemEnergy restrict( const Indexer &coarseProlongationIndexer ) const;

		// Update the system given the finer solution
		void update( const HierarchicalSystemEnergy &finer , const Eigen::VectorXd &x , const Eigen::VectorXd &y );

		// Separately computes the energy components given by the alternating-field fitting (first) and regularization (second) energy terms
		std::pair< double , double > energies( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;

		// Given the current estimate of the solution (x,y) computes the quadratic polynomial P(s,t) = E(x+s*e[idx],y+t*e[idx])
		Polynomial::Polynomial2D< 2 > quadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int idx , bool setConstantTerm , unsigned int thread ) const;

		// Returns the square norm of the target alternating product
		double squareNorm( void ) const { return _c[0]; }

		// Returns the square norm of the wedge
		double wedgeSquareNorm( const Eigen::VectorXd &w ) const { return w.dot( _M[0][0] * w ); }

		// Returns the prolongation of the scalar function
		Eigen::VectorXd prolongation( const Eigen::VectorXd &coarse ) const { return _sP * coarse; }

		// Returns the scalar prolongation matrix
		const Eigen::SparseMatrix< double > &prolongation( void ) const { return _sP; }

		// Returns the wedge prolongation matrix
		const Eigen::SparseMatrix< double > &wedgeProlongation( void ) const { return _wP; }

		// Given the current estimate of the solution (x,y) and the two directions (_x,_y), computes the bi-quadratic polynomial P(s,t) = E(x+s*_x,y+t*_y)
		Polynomial::Polynomial2D< 4 > biQuadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , const Eigen::VectorXd &_x , const Eigen::VectorXd &_y ) const;

		// Given the current estimate of the solution (x,y), gives the quadratic and linear approximations to the energy in x
		typename SystemEnergy< Dim , Indexer >::QuadraticApproximation quadraticApproximation1( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;

		// Given the current estimate of the solution (x,y), gives the quadratic and linear approximations to the energy in y
		typename SystemEnergy< Dim , Indexer >::QuadraticApproximation quadraticApproximation2( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;

		// Computes the differential of the energy with respect to the two function's coefficients (faster than the base implementation)
		std::pair< Eigen::VectorXd , Eigen::VectorXd > d( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;
	protected:
		HierarchicalSystemEnergy( const Indexer &indexer , unsigned int r );

		bool _initFromFiner;
		Eigen::SparseMatrix< char > _sI;				// Scalar incidence matrix
		Eigen::SparseMatrix< double > _sP , _wP;		// The scalar- and wedge-product prolongation matrix (from the coarse resolution into this one)
		Eigen::SparseMatrix< double > _M[3][3];			// The wedge-product mass-matrices
		Eigen::SparseMatrix< double > _R;				// The symmetric matrix defining the regularization energy
		Eigen::VectorXd _w[3] , _b[3] , _r[2];			// The integral of the target alternating-form field against the alternating-form basis
		double _c[2];									// The constant terms coming from fitting (c[0]) and regularization (c[1])
	};

	//////////////////
	// SystemEnergy //
	//////////////////

	template< unsigned int Dim , typename Indexer >
	double SystemEnergy< Dim , Indexer >::operator()( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const
	{
		std::pair< double , double > errs = energies( x , y );
		return errs.first + errs.second;
	}

	template< unsigned int Dim , typename Indexer >
	std::pair< Eigen::VectorXd , Eigen::VectorXd > SystemEnergy< Dim , Indexer >::d( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const
	{
		QuadraticApproximation q1 = quadraticApproximation1( x , y );
		QuadraticApproximation q2 = quadraticApproximation2( x , y );
		return std::make_pair( 2. * q1.q * x + q1.l , 2. * q2.q * y + q2.l );
	}

	template< unsigned int Dim , typename Indexer >
	double SystemEnergy< Dim , Indexer >::stepSize( const Eigen::VectorXd &x , const Eigen::VectorXd &y , const Eigen::VectorXd &dx , const Eigen::VectorXd &dy ) const
	{
		Matrix< double , 2 , 2 > M;
		M(0,0) = M(0,1) = 1.;

		Polynomial::Polynomial1D< 4 > Q = biQuadraticFit( x , y , dx , dy ).template pullBack< 2 >( M );
		Polynomial::Polynomial1D< 3 > dQ = Q.d(0);

		double roots[3];
		unsigned int rNum = Polynomial::Roots( dQ , roots );
		if( !rNum ) ERROR_OUT( "Expected a root for an odd-degree polynomial" );
		double s = roots[0];
		for( unsigned int i=1 ; i<rNum ; i++ ) if( Q( roots[i] )<Q(s) ) s = roots[i];
		return s;
	}

	template< unsigned int Dim , typename Indexer >
	Point< double , 2 > SystemEnergy< Dim , Indexer >::newtonUpdate( const Eigen::VectorXd &x , const Eigen::VectorXd &y , const Eigen::VectorXd &dx , const Eigen::VectorXd &dy , unsigned int steps ) const
	{
		Point< double , 2 > c;
		Polynomial::Polynomial2D< 4 > Q2 = biQuadraticFit( x , y , dx , dy );
		for( unsigned int step=0 ; step<steps ; step++ )
		{
			Point< double , 2 > d = - Q2.hessian( c ).inverse() * Q2.gradient( c );
			if( d.squareNorm()<1e-20 ) break;

			// Compute the 1D polynomial Q1(s) = cQ2( p * s );
			Polynomial::Polynomial1D< 4 > Q1;
			{
				Matrix< double , 2 , 2 > S;
				S(0,0) = d[0];
				S(0,1) = d[1];
				S(1,0) = c[0];
				S(1,1) = c[1];
				Q1 = Q2.template pullBack< 2 >( S );
			}
			Polynomial::Polynomial1D< 3 > dQ1 = Q1.d(0);

			double roots[3];
			unsigned int rNum = Polynomial::Roots( dQ1 , roots );
			if( !rNum ) ERROR_OUT( "Expected a root for an odd-degree polynomial: " , Q1 );
			double s = roots[0];
			for( unsigned int i=1 ; i<rNum ; i++ ) if( Q1( roots[i] )<Q1(s) ) s = roots[i];
			c += d * s;
		}
		return c;
	}

	template< unsigned int Dim , typename Indexer >
	Polynomial::Polynomial2D< 2 > SystemEnergy< Dim , Indexer >::quadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int idx , bool setConstantTerm , unsigned int thread ) const
	{
		Eigen::VectorXd d;
		d.setZero( x.size() );
		d[idx] = 1.;
		return biQuadraticFit( x , y , d , d );
	}

	template< unsigned int Dim , typename Indexer >
	Polynomial::Polynomial2D< 4 > SystemEnergy< Dim , Indexer >::biQuadraticFit( const Eigen::VectorXd & , const Eigen::VectorXd & , const Eigen::VectorXd & , const Eigen::VectorXd & ) const
	{
		ERROR_OUT( "Method not supported" );
		return Polynomial::Polynomial2D< 4 >();
	}

	template< unsigned int Dim , typename Indexer >
	typename SystemEnergy< Dim , Indexer >::QuadraticApproximation SystemEnergy< Dim , Indexer >::quadraticApproximation1( const Eigen::VectorXd & , const Eigen::VectorXd & ) const
	{
		ERROR_OUT( "Method not supported" );
		return typename SystemEnergy< Dim , Indexer >::QuadraticApproximation();
	}

	template< unsigned int Dim , typename Indexer >
	typename SystemEnergy< Dim , Indexer >::QuadraticApproximation SystemEnergy< Dim , Indexer >::quadraticApproximation2( const Eigen::VectorXd & , const Eigen::VectorXd & ) const
	{
		ERROR_OUT( "Method not supported" );
		return typename SystemEnergy< Dim , Indexer >::QuadraticApproximation();
	}

	template< unsigned int Dim , typename Indexer >
	void SystemEnergy< Dim , Indexer >::test( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int numTests , double eps ) const
	{
		testD( x , y , numTests , eps );
		testBiQuadraticFit( x , y , numTests , eps );
		testQuadraticApproximation1( x , y , numTests , eps );
		testQuadraticApproximation2( x , y , numTests , eps );
		testQuadraticApproximation( x , y , numTests );
	}

	template< unsigned int Dim , typename Indexer >
	void SystemEnergy< Dim , Indexer >::testD( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int numTests , double eps ) const
	{
		std::cout << "Finite difference vs. Gradient (eps = " << eps << "):" << std::endl;

		std::pair< Eigen::VectorXd , Eigen::VectorXd > diff = d(x,y);

		for( unsigned int i=0 ; i<numTests ; i++ )
		{
			Eigen::VectorXd _x = x;
			Eigen::VectorXd _y = y;
			for( unsigned int i=0 ; i<x.size() ; i++ ) _x[i] += ( Random< double >()*2.-1. ) * eps;
			for( unsigned int i=0 ; i<y.size() ; i++ ) _y[i] += ( Random< double >()*2.-1. ) * eps;
			double d1 = ( operator()( _x , _y ) - operator()( x , y ) ) / eps , d2 = ( diff.first.dot(_x-x) + diff.second.dot(_y-y) ) / eps;
			std::cout << "\t" << fabs( d1 - d2 ) << " <- " << d1 << " / " << d2 << std::endl;
		}
	}

	template< unsigned int Dim , typename Indexer >
	void SystemEnergy< Dim , Indexer >::testBiQuadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int numTests , double eps ) const
	{
		std::cout << "Energy vs. polynomial 2D:" << std::endl;

		for( unsigned int i=0 ; i<numTests ; i++ )
		{
			Eigen::VectorXd _x( x.size() ) , _y( y.size() );
			for( unsigned int i=0 ; i<x.size() ; i++ ) _x[i] = ( Random< double >()*2.-1. ) * eps;
			for( unsigned int i=0 ; i<y.size() ; i++ ) _y[i] = ( Random< double >()*2.-1. ) * eps;
			Polynomial::Polynomial2D< 4 > Q = biQuadraticFit( x , y , _x , _y );

			Eigen::VectorXd __x = x + _x , __y = y + _y;
			double v1 = operator()( __x , __y ) , v2 = Q(1.,1.);
			std::cout << "\t" << fabs( v1 - v2 ) << " <- " << v1 << " / " << v2 << std::endl;
		}
	}

	template< unsigned int Dim , typename Indexer >
	void SystemEnergy< Dim , Indexer >::testQuadraticApproximation1( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int numTests , double eps ) const
	{
		std::cout << "Energy vs. quadratic approximation 1:" << std::endl;
		QuadraticApproximation q1 = quadraticApproximation1( x , y );

		Eigen::VectorXd _x( x.size() );
		for( unsigned int i=0 ; i<numTests ; i++ )
		{
			for( unsigned int j=0 ; j<x.size() ; j++ ) _x[j] = x[j] + ( Random< double >()*2. - 1. ) * eps;
			double v1 = operator()( _x , y ) , v2 = q1( _x );
			std::cout << "\t" << fabs( v1 - v2 ) << " <- " << v1 << " / " << v2 << std::endl;
		}
	}

	template< unsigned int Dim , typename Indexer >
	void SystemEnergy< Dim , Indexer >::testQuadraticApproximation2( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int numTests , double eps ) const
	{
		std::cout << "Energy vs. quadratic approximation 2:" << std::endl;
		QuadraticApproximation q2 = quadraticApproximation2( x , y );

		Eigen::VectorXd _y( y.size() );
		for( unsigned int i=0 ; i<numTests ; i++ )
		{
			for( unsigned int j=0 ; j<y.size() ; j++ ) _y[j] = y[j] + ( Random< double >()*2. - 1. ) * eps;
			double v1 = operator()( x , _y ) , v2 = q2( _y );
			std::cout << "\t" << fabs( v1 - v2 ) << " <- " << v1 << " / " << v2 << std::endl;
		}
	}

	template< unsigned int Dim , typename Indexer >
	void SystemEnergy< Dim , Indexer >::testQuadraticApproximation( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int numTests ) const
	{
		std::cout << "Bi-quadratic approximation vs. quadratic approximation:" << std::endl;

		Eigen::VectorXd d;
		d.setZero( x.size() );
		for( unsigned int i=0 ; i<numTests ; i++ )
		{
			unsigned int idx = rand() % x.size();
			d[idx] = 1.;

			Polynomial::Polynomial2D< 4 > Q4 = biQuadraticFit( x , y , d , d );
			Polynomial::Polynomial2D< 2 > Q2 = quadraticFit( x , y , idx , false , 0 );

			std::cout << "\t" << ( Q4 - Q2 ) << std::endl;

			d[idx] = 0.;
		}
	}


	/////////////////////////////////////
	// SingleCycleCascadicSystemEnergy //
	/////////////////////////////////////
	template< unsigned int Dim , typename Indexer >
	SingleCycleCascadicSystemEnergy< Dim , Indexer >::SingleCycleCascadicSystemEnergy( const Indexer &indexer , unsigned int r )
		: SystemEnergy< Dim , Indexer >( indexer , r ) , _wMStencil( Hat::AlternatingProductFunctions< Dim >::MassStencil(r) , r ) , _sSStencil( Hat::ScalarFunctions< Dim >::StiffnessStencil(r) , r ) , _c(0) , _sWeight(0)
	{}

	template< unsigned int Dim , typename Indexer >
	template< typename SkewField /* = std::function< SkewSymmetriMatrix< double , Dim > ( Hat::Index< Dim > ) > */ >
	SingleCycleCascadicSystemEnergy< Dim , Indexer >::SingleCycleCascadicSystemEnergy( const Indexer &indexer , unsigned int r , SkewField && skew , double sWeight )
		: SingleCycleCascadicSystemEnergy< Dim , Indexer >( indexer , r , std::forward< SkewField >( skew ) , sWeight , Eigen::SparseMatrix< double >( indexer.numFunctions() , indexer.numFunctions() ) )
	{}

	template< unsigned int Dim , typename Indexer >
	template< typename SkewField /* = std::function< SkewSymmetriMatrix< double , Dim > ( Hat::Index< Dim > ) > */ >
	SingleCycleCascadicSystemEnergy< Dim , Indexer >::SingleCycleCascadicSystemEnergy( const Indexer &indexer , unsigned int r , SkewField && skew , double sWeight , Eigen::SparseMatrix< double > R )
		: SingleCycleCascadicSystemEnergy< Dim , Indexer >( indexer , r )
	{
		_sWeight = sWeight;
		_B = scalars.stiffness( indexer , [&]( Hat::Index< Dim > E , unsigned int t ){ return skew( E , t )(); } );
		_R = R;
		ThreadPool::ParallelFor
		(
			0 , indexer.numElements() ,
			[&]( unsigned int t , size_t e ){ AddAtomic( _c , skew( indexer.elementIndex(e) , t ).squareNorm() ); }
		);

		_c /= scalars.elementNum();
	}

	template< unsigned int Dim , typename Indexer >
	SingleCycleCascadicSystemEnergy< Dim , Indexer > SingleCycleCascadicSystemEnergy< Dim , Indexer >::restrict( const Indexer &coarseProlongationIndexer ) const
	{
		static_assert( std::is_base_of_v< Hat::BaseProlongationIndexer< Dim > , Indexer > , "[ERROR] Indexer poorly formed" );

		if( scalars.resolution()&1 ) ERROR_OUT( "Expected even resolution: " , scalars.resolution() );
		SingleCycleCascadicSystemEnergy coarser( coarseProlongationIndexer , scalars.resolution() / 2 );
		coarser._sP = coarser.scalars.prolongation( coarseProlongationIndexer , indexer.numFunctions() );
		coarser._R = coarser._sP.transpose() * _R * coarser._sP;
		coarser._B = coarser._sP.transpose() * _B * coarser._sP;
		coarser._c = _c;
		coarser._sWeight = _sWeight;
		return coarser;
	}

	template< unsigned int Dim , typename Indexer >
	void SingleCycleCascadicSystemEnergy< Dim , Indexer >::update( const SingleCycleCascadicSystemEnergy &finer , const Eigen::VectorXd &x , const Eigen::VectorXd &y )
	{
		if( finer.scalars.resolution()!=scalars.resolution()*2 ) ERROR_OUT( "Resolutions don't match: " , finer.scalars.resolution() , " != 2 * " , scalars.resolution() );
		if( !( x.isZero() && y.isZero() ) ) ERROR_OUT( "Non-trivial restriction not supported" );
	}

	template< unsigned int Dim , typename Indexer >
	std::pair< double , double > SingleCycleCascadicSystemEnergy< Dim , Indexer >::energies( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const
	{
		// E(x,y) = || dx ^ dy - b ||^2  + x^t * R * x + y^t * R * y 
		//        = || dx ^ dy ||^2 - 2 < dx ^ dy , b > + x^t * R * x + y^t * R * y
		std::pair< double , double > e;

		e.first = wedges( indexer , _wMStencil , x , y , x , y ) - 2. * y.dot( _B * x ) + _c;
		e.second = ( _R * x ).dot( x ) + ( _R * y ).dot( y );
		if( _sWeight ) e.second += ( scalars( indexer , _sSStencil , x , x ) + scalars( indexer , _sSStencil , y , y ) ) * _sWeight;
		return e;
	}

	template< unsigned int Dim , typename Indexer >
	Polynomial::Polynomial2D< 4 > SingleCycleCascadicSystemEnergy< Dim , Indexer >::biQuadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , const Eigen::VectorXd &_x , const Eigen::VectorXd &_y ) const
	{
		Polynomial::Polynomial2D< 4 > Q2;
		Eigen::VectorXd xy = wedges.product( x , y ) , _x_y = wedges.product( _x , _y );
		Eigen::VectorXd _xy = wedges.product( _x , y ) , x_y = wedges.product( x , _y );

		Q2.coefficient(2,2) +=     wedges( _wMStencil , _x_y , _x_y );
		Q2.coefficient(2,1) += 2 * wedges( _wMStencil , _x_y , _xy );
		Q2.coefficient(1,2) += 2 * wedges( _wMStencil , _x_y , xy );
		Q2.coefficient(2,0) +=     wedges( _wMStencil , _xy , _xy );
		Q2.coefficient(0,2) +=     wedges( _wMStencil , x_y , x_y );
		Q2.coefficient(1,1) += 2 * wedges( _wMStencil , _xy , x_y ) + 2 * wedges( _wMStencil , _x_y , xy );
		Q2.coefficient(1,0) += 2 * wedges( _wMStencil , _xy , xy );
		Q2.coefficient(0,1) += 2 * wedges( _wMStencil , x_y , xy );
		Q2.coefficient(0,0) +=     wedges( _wMStencil , xy , xy );

		Q2.coefficient(2,0) +=      ( _R * _x ).dot( _x );
		Q2.coefficient(1,0) += 2. * ( _R * x ).dot( _x );
		Q2.coefficient(0,0) +=      ( _R * x ).dot( x );

		Q2.coefficient(0,2) +=      ( _R * _y ).dot( _y );
		Q2.coefficient(0,1) += 2. * ( _R * y ).dot( _y );
		Q2.coefficient(0,0) +=      ( _R * y ).dot( y );

		if( _sWeight )
		{
			Q2.coefficient(2,0) +=      scalars( indexer , _sSStencil , _x , _x );
			Q2.coefficient(1,0) += 2. * scalars( indexer , _sSStencil , x , _x );
			Q2.coefficient(0,0) +=      scalars( indexer , _sSStencil , x , x );

			Q2.coefficient(0,2) +=      scalars( indexer , _sSStencil , _y , _y );
			Q2.coefficient(0,1) += 2. * scalars( indexer , _sSStencil , y , _y );
			Q2.coefficient(0,0) +=      scalars( indexer , _sSStencil , y , y );
		}

		Q2.coefficient(1,1) += - 2. * _y.dot( _B * _x );
		Q2.coefficient(0,1) += - 2. * _y.dot( _B * x );
		Q2.coefficient(1,0) += - 2. * y.dot( _B * _x );
		Q2.coefficient(0,0) += - 2. * y.dot( _B * x );

		Q2.coefficient(0,0) += _c;

		return Q2;
	}

	template< unsigned int Dim , typename Indexer >
	Polynomial::Polynomial2D< 2 > SingleCycleCascadicSystemEnergy< Dim , Indexer >::quadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int idx , bool setConstantTerm , unsigned int thread ) const
	{
		Window::IsotropicStaticWindow< size_t , Dim , 3 > neighbors = indexer.ffNeighbors( idx , thread );
		Hat::Index< Dim > Off;
		for( unsigned int d=0 ; d<Dim ; d++ ) Off[d] = 1;

		Polynomial::Polynomial2D< 2 > Q2;

		Hat::Index< Dim > F1 = indexer.functionIndex( idx );

		//	 E(a,b) = || (x + a e_i) ^ ( y + b e_i ) - w ||^2 + ...
		//	        = || a e_i^y + b x ^ e_i + x^y - w ||^2 + ...
		//	        = a^2 || e_i^y ||^2 + b^2 || x^e_i ||^2 + 2 a b < e_i^y , x^e_i > + 2a < e_i^y  , x^y - w > + 2b < x^e_i , x^y  - w > + ...
		//	        = a^2 || e_i^y ||^2 + b^2 || x^e_i ||^2 - 2 a b < e_i^x , e_i^y > + 2a < e_i^y  , x^y - w > + 2b < x^e_i , x^y  - w > + ...
		// Considering the individual terms:
		//	|| e_i^y ||^2 = || \sum_j y[j] e_i^e_j ||^2 = \sum_{j,k} y[j] y[k] < e_i^e_j , e_i^e_k >
		//	|| x^e_i ||^2 = || \sum_j x[j] e_j^e_i ||^2 = \sum_{j,k} x[j] x[k] < e_i^e_j , e_i^e_k >
		//	< e_i^x , e_i^y > = \sum_{j,k} x[j] y[k] < e_i^e_j , e_i^e_k>
		//	< e_i^y , x^y - w > =   \sum_{j,k,l} y[j] x[k] y[l] < e_i^e_j , e_k^e_l > - \sum_j y[j] < e_i^e_j , w >
		//	< x^e_i , x^y - w > = - \sum_{j,k,l} x[j] x[k] y[l] < e_i^e_j , e_k^e_l > + \sum_j x[j] < e_i^e_j , w >
		for( Eigen::InnerIterator it(_R,idx) ; it ; ++it )
		{
			if( it.row()==idx )
			{
				Q2.coefficient(2,0) += it.value();
				Q2.coefficient(0,2) += it.value();
			}
			Q2.coefficient(1,0) += 2. * it.value() * x[ it.row() ];
			Q2.coefficient(0,1) += 2. * it.value() * y[ it.row() ];
		}

		if( _sWeight )
		{
			const typename Hat::ScalarFunctions< Dim >::template FullIntegrationStencil< double , 1 >::Row &row = _sSStencil.row( F1 );
			for( unsigned int i=0 ; i<row.size() ; i++ )
			{
				double value = row[i].second * _sWeight;
				Hat::Index< Dim > I = row[i].first + Off;
				size_t _idx = neighbors( &I[0] );
				if( _idx!=-1 )
				{
					if( idx==_idx )
					{
						Q2.coefficient(2,0) += value;
						Q2.coefficient(0,2) += value;
					}
					Q2.coefficient(1,0) += 2. * value * x[_idx];
					Q2.coefficient(0,1) += 2. * value * y[_idx];
				}
			}
		}

		for( Eigen::InnerIterator it(_B,idx) ; it ; ++it ) if( it.row()!=it.col() )
		{
			Q2.coefficient(1,0) -= 2. * it.value() * y[ it.row() ];
			Q2.coefficient(0,1) += 2. * it.value() * x[ it.row() ];
		}

		// The index 
		size_t i1 = idx;

		const std::vector< typename Hat::AlternatingProductFunctions< Dim >::template FullIntegrationStencil< double >::Row > &rows = _wMStencil.rows( F1 );
		for( unsigned int i=0 ; i<rows.size() ; i++ )
		{
			Hat::Index< Dim > I = rows[i].f2 + Off;
			size_t i2 = neighbors( &I[0] );
			if( i2==-1 ) continue;
			const std::vector< typename Hat::AlternatingProductFunctions< Dim >::template FullIntegrationStencil< double >::Entry > &entries = rows[i].entries;
			for( unsigned int j=0 ; j<entries.size() ; j++ )
			{
				Hat::Index< Dim > J1 = entries[j].g1 + Off;
				Hat::Index< Dim > J2 = entries[j].g2 + Off;
				size_t j1 = neighbors( &J1[0] );
				size_t j2 = neighbors( &J2[0] );

				if( j1==-1 || j2==-1 ) continue;

				double xy = ( x[j1] * y[j2] - x[j2] * y[j1] );
				double s = entries[j].value;

				//	< e_i^y , x^y > = \sum_{j,k,l} y[j] x[k] y[l] < e_i^e_j , e_k^e_l >
				Q2.coefficient(1,0) += 2 * s * y[i2] * xy;
				Q2.coefficient(0,1) -= 2 * s * x[i2] * xy;

				if( j1==i1 )
				{
					//	|| e_i ^ y ||^2 = || \sum_j y_j e_i^e_j ||^2 = \sum_{j,k} y_j y_k < e_i^e_j , e_i^e_k >
					Q2.coefficient(2,0) +=     s * y[i2] * y[j2];

					//	|| x ^ e_i ||^2 = || \sum_j x_j e_j^e_i ||^2 = \sum_{j,k} x_j x_k < e_j^e_i , e_k^e_i > = \sum_{j,k} x_j x_k < e_i^e_j , e_i^e_k >
					Q2.coefficient(0,2) +=     s * x[i2] * x[j2];

					//	< e_i ^ y , x ^ e_i > = \sum_{j,k} x_j y_k < e_i^e_j , e_i^e_k>
					Q2.coefficient(1,1) -= 2 * s * x[i2] * y[j2];
				}
				if( j2==i1 )
				{
					Q2.coefficient(2,0) -=     s * y[i2] * y[j1];
					Q2.coefficient(0,2) -=     s * x[i2] * x[j1];
					Q2.coefficient(1,1) += 2 * s * x[i2] * y[j1];
				}
			}
		}
		return Q2;
	}

	//////////////////////////////
	// HierarchicalSystemEnergy //
	//////////////////////////////
	template< unsigned int Dim , typename Indexer >
	HierarchicalSystemEnergy< Dim , Indexer >::HierarchicalSystemEnergy( const Indexer &indexer , unsigned int r ) : SystemEnergy< Dim , Indexer >( indexer , r ) , _initFromFiner( false )
	{
		if( !(r&1) ) _wP = wedges.prolongation();
		_sI = scalars.incidence();
		_M[0][0] = wedges.mass();
		_c[0] = 0;
	}

	template< unsigned int Dim , typename Indexer >
	HierarchicalSystemEnergy< Dim , Indexer >::HierarchicalSystemEnergy( const Indexer &indexer , unsigned int r , ConstPointer( Hat::SkewSymmetricMatrix< double , Dim > ) skew , double sWeight )
		: HierarchicalSystemEnergy< Dim , Indexer >( indexer , r , skew , sWeight , Eigen::SparseMatrix< double >( indexer.numFunctions() , indexer.numFunctions() ) )
	{}

	template< unsigned int Dim , typename Indexer >
	HierarchicalSystemEnergy< Dim , Indexer >::HierarchicalSystemEnergy( const Indexer &indexer , unsigned int r , ConstPointer( Hat::SkewSymmetricMatrix< double , Dim > ) skew , double sWeight , Eigen::SparseMatrix< double > R )
		: HierarchicalSystemEnergy< Dim , Indexer >( indexer , r )
	{
		_R = R;
		if( sWeight>0 ) _R += scalars.stiffness() * sWeight;
		_b[0] = wedges.dualVector(skew);
		for( unsigned int i=0 ; i<indexer.elementNum() ; i++ ) _c[0] += skew[i].squareNorm();
		_c[0] /= scalars.elementNum();
	}

	template< unsigned int Dim , typename Indexer >
	HierarchicalSystemEnergy< Dim , Indexer > HierarchicalSystemEnergy< Dim , Indexer >::restrict( const Indexer &coarseIndexer ) const
	{
		if( scalars.resolution()&1 ) ERROR_OUT( "Expected even resolution: " , scalars.resolution() );
		HierarchicalSystemEnergy coarser( coarseIndexer , scalars.resolution() / 2 );
		coarser._sP = coarser.scalars.prolongation();
		coarser._R = coarser._sP.transpose() * _R * coarser._sP;
		coarser._b[0] = _wP.transpose() * _b[0];
		coarser._c[0] = _c[0];
		return coarser;
	}

	template< unsigned int Dim , typename Indexer >
	void HierarchicalSystemEnergy< Dim , Indexer >::update( const HierarchicalSystemEnergy &finer , const Eigen::VectorXd &x , const Eigen::VectorXd &y )
	{
		if( finer.scalars.resolution()!=scalars.resolution()*2 ) ERROR_OUT( "Resolutions don't match: " , finer.scalars.resolution() , " != 2 * " , scalars.resolution() );
		if( x.size()!=finer.indexer.numFunctions() ) ERROR_OUT( "Unexpected x resolution: " , x.size() , " != " , finer.indexer.numFunctions() );
		if( y.size()!=finer.indexer.numFunctions() ) ERROR_OUT( "Unexpected y resolution: " , y.size() , " != " , finer.indexer.numFunctions() );

		bool hasXY = !( x.isZero() && y.isZero() );
		_initFromFiner = finer._initFromFiner || hasXY;

		const Eigen::SparseMatrix< double > &sP = finer._sP , &wP = finer._wP;
		Eigen::SparseMatrix< double > sP_t = sP.transpose() , wP_t = wP.transpose();

		_c[0] = finer._c[0];
		_c[1] = finer._c[1];

		for( unsigned int i=1 ; i<3 ; i++ ) for( unsigned int j=0 ; j<=i ; j++ ) _M[i][j].setZero();
		for( unsigned int i=1 ; i<3 ; i++ ) _b[i].setZero();
		for( unsigned int i=0 ; i<3 ; i++ ) _w[i].setZero();
		for( unsigned int i=0 ; i<2 ; i++ ) _r[i].setZero();

		if( finer._initFromFiner )
		{
			for( unsigned int i=1 ; i<3 ; i++ ) for( unsigned int j=0 ; j<=i ; j++ ) _M[i][j] += sP_t * finer._M[i][j] * ( j==0 ? wP : sP );

			_b[1] += sP_t * finer._b[1];
			_b[2] += sP_t * finer._b[2];

			_w[0] += wP_t * finer._w[0];
			_w[1] += sP_t * finer._w[1];
			_w[2] += sP_t * finer._w[2];

			_r[0] += sP_t * finer._r[0];
			_r[1] += sP_t * finer._r[1];
		}

		if( hasXY )
		{
			Eigen::VectorXd xy = finer.wedges.wedge(x,y);
			Eigen::SparseMatrix< double > P[2] , P_t[2];
			P[0] =  finer.wedges.wedge( x ) * sP;
			P[1] = -finer.wedges.wedge( y ) * sP;
			for( unsigned int i=0 ; i<2 ; i++ ) P_t[i] = P[i].transpose();

			_w[0] += wP_t * finer._M[0][0] * xy;
			for( unsigned int i=1 ; i<3 ; i++ ) _b[i] += P_t[i-1] * finer._b[0];
			for( unsigned int i=1 ; i<3 ; i++ ) _w[i] += P_t[i-1] * finer._M[0][0] * xy;
			for( unsigned int i=1 ; i<3 ; i++ ) for( unsigned int j=0 ; j<=i ; j++ ) _M[i][j] += P_t[i-1] * finer._M[0][0] * ( j==0 ? wP : P[j-1] );

			_r[0] += sP_t * finer._R * x;
			_r[1] += sP_t * finer._R * y;

			_c[0] += ( finer._M[0][0] * xy ).dot( xy ) - 2. * xy.dot( finer._b[0] );
			_c[1] += ( finer._R * x ).dot( x ) + ( finer._R * y ).dot( y );

			if( finer._initFromFiner )
			{
				for( unsigned int i=1 ; i<3 ; i++ ) for( unsigned int j=1 ; j<=i ; j++ ) _M[i][j] += P_t[i-1] * finer._M[0][j] * sP + sP_t * finer._M[i][0] * P[j-1];

				_w[0] +=   wP_t * ( finer._M[0][1] * y + finer._M[0][2] * x );
				_w[1] += P_t[0] * ( finer._M[0][2] * x + finer._M[0][1] * y + finer._w[0] ) + sP_t * ( finer._M[1][0] * xy + finer._M[1][2] * x + finer._M[1][1] * y );
				_w[2] += P_t[1] * ( finer._M[0][2] * x + finer._M[0][1] * y + finer._w[0] ) + sP_t * ( finer._M[2][0] * xy + finer._M[2][2] * x + finer._M[2][1] * y );

				_c[0] +=    ( finer._M[2][2] * x ).dot( x );
				_c[0] +=    ( finer._M[1][1] * y ).dot( y );
				_c[0] += 2.*( finer._M[2][1] * y ).dot( x );

				_c[0] += 2.*( finer._M[0][2] * x ).dot( xy );
				_c[0] += 2.*( finer._M[0][1] * y ).dot( xy );

				_c[0] += 2.*xy.dot( finer._w[0] );

				_c[0] += 2.*x.dot( finer._w[2] );
				_c[0] += 2.*y.dot( finer._w[1] );

				_c[0] -= 2.*finer._b[2].dot( x );
				_c[0] -= 2.*finer._b[1].dot( y );

				_c[1] += 2.*x.dot( finer._r[0] ) + 2.*y.dot( finer._r[1] );
			}
		}

		if( _initFromFiner ) for( unsigned int i=0 ; i<3 ; i++ ) for( unsigned int j=0 ; j<i ; j++ ) _M[j][i] = _M[i][j].transpose();
	}

	template< unsigned int Dim , typename Indexer >
	std::pair< double , double > HierarchicalSystemEnergy< Dim , Indexer >::energies( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const
	{
		std::pair< double , double > e;

		Eigen::VectorXd xy = wedges.product( x , y );

		e.first = ( _M[0][0] * xy ).dot( xy ) - 2. * xy.dot( _b[0] ) + _c[0];
		e.second = ( _R * x ).dot( x ) + ( _R * y ).dot( y ) + _c[1];

		if( _initFromFiner )
		{
			e.first += 2. * ( _M[0][2] * x ).dot( xy );
			e.first += 2. * ( _M[0][1] * y ).dot( xy );
			e.first +=      ( _M[2][2] * x ).dot( x );
			e.first +=      ( _M[1][1] * y ).dot( y );
			e.first += 2. * xy.dot( _w[0] ) + 2. * ( _M[1][2] * x ).dot( y );
			e.first += 2. * x.dot( _w[2] ) - 2. * x.dot( _b[2] ) + 2. * y.dot( _w[1] ) - 2. * y.dot( _b[1] );
			e.second += 2. * x.dot( _r[0] ) + 2. * y.dot( _r[1] );
		}
		return e;
	}

	template< unsigned int Dim , typename Indexer >
	Polynomial::Polynomial2D< 4 > HierarchicalSystemEnergy< Dim , Indexer >::biQuadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , const Eigen::VectorXd &_x , const Eigen::VectorXd &_y ) const
	{
		Polynomial::Polynomial2D< 4 > Q2;
		Eigen::VectorXd xy = wedges.wedge( x , y ) , _x_y = wedges.wedge( _x , _y );
		Eigen::VectorXd _xy = wedges.wedge( _x , y ) , x_y = wedges.wedge( x , _y );

		Q2.coefficient(2,2) +=     ( _M[0][0] * _x_y ).dot( _x_y );
		Q2.coefficient(2,1) += 2 * ( _M[0][0] * _x_y ).dot( _xy );
		Q2.coefficient(1,2) += 2 * ( _M[0][0] * _x_y ).dot( x_y );
		Q2.coefficient(2,0) +=     ( _M[0][0] * _xy ).dot( _xy );
		Q2.coefficient(0,2) +=     ( _M[0][0] * x_y ).dot( x_y );
		Q2.coefficient(1,1) += 2 * ( _M[0][0] * _xy ).dot( x_y ) + 2 * ( _M[0][0] * _x_y ).dot( xy );
		Q2.coefficient(1,0) += 2 * ( _M[0][0] * _xy ).dot( xy );
		Q2.coefficient(0,1) += 2 * ( _M[0][0] * x_y ).dot( xy );
		Q2.coefficient(0,0) +=     ( _M[0][0] * xy ).dot( xy );

		Q2.coefficient(2,0) +=      ( _R * _x ).dot( _x );
		Q2.coefficient(1,0) += 2. * ( _R * x ).dot( _x );
		Q2.coefficient(0,0) +=      ( _R * x ).dot( x );

		Q2.coefficient(0,2) +=      ( _R * _y ).dot( _y );
		Q2.coefficient(0,1) += 2. * ( _R * y ).dot( _y );
		Q2.coefficient(0,0) +=      ( _R * y ).dot( y );

		Q2.coefficient(1,1) += - 2. * _x_y.dot( _b[0] );
		Q2.coefficient(0,1) += - 2. * x_y.dot( _b[0] );
		Q2.coefficient(1,0) += - 2. * _xy.dot( _b[0] );
		Q2.coefficient(0,0) += - 2. * xy.dot( _b[0] );

		Q2.coefficient(0,0) += _c[0] + _c[1];

		if( _initFromFiner )
		{
			Q2.coefficient(2,1) += 2. * ( _M[0][2] * _x ).dot( _x_y );
			Q2.coefficient(2,0) += 2. * ( _M[0][2] * _x ).dot( _xy );
			Q2.coefficient(1,1) += 2. * ( _M[0][2] * _x ).dot( x_y ) + 2. * ( _M[0][2] * x ).dot( _x_y );
			Q2.coefficient(1,0) += 2. * ( _M[0][2] * _x ).dot( xy ) + 2. * ( _M[0][2] * x ).dot( _xy );
			Q2.coefficient(0,1) += 2. * ( _M[0][2] *  x ).dot( x_y );
			Q2.coefficient(0,0) += 2. * ( _M[0][2] *  x ).dot( xy );

			Q2.coefficient(1,2) += 2. * ( _M[0][1] * _y ).dot( _x_y );
			Q2.coefficient(0,2) += 2. * ( _M[0][1] * _y ).dot( x_y );
			Q2.coefficient(1,1) += 2. * ( _M[0][1] * _y ).dot( _xy ) + 2. * ( _M[0][1] * y ).dot( _x_y );
			Q2.coefficient(0,1) += 2. * ( _M[0][1] * _y ).dot( xy ) + 2. * ( _M[0][1] * y ).dot( x_y );
			Q2.coefficient(1,0) += 2. * ( _M[0][1] *  y ).dot( _xy );
			Q2.coefficient(0,0) += 2. * ( _M[0][1] *  y ).dot( xy );

			Q2.coefficient(2,0) +=      ( _M[2][2] * _x ).dot( _x );
			Q2.coefficient(1,0) += 2. * ( _M[2][2] * x ).dot( _x );
			Q2.coefficient(0,0) +=      ( _M[2][2] * x ).dot( x );

			Q2.coefficient(0,2) +=      ( _M[1][1] * _y ).dot( _y );
			Q2.coefficient(0,1) += 2. * ( _M[1][1] * y ).dot( _y );
			Q2.coefficient(0,0) +=      ( _M[1][1] * y ).dot( y );

			Q2.coefficient(1,1) += 2. * _x_y.dot( _w[0] ) + 2. * ( _M[1][2] * _x ).dot( _y );
			Q2.coefficient(0,1) += 2. * x_y.dot( _w[0] ) + 2. * ( _M[1][2] * x ).dot( _y );
			Q2.coefficient(1,0) += 2. * _xy.dot( _w[0] ) + 2. * ( _M[1][2] * _x ).dot( y );
			Q2.coefficient(0,0) += 2. * xy.dot( _w[0] ) + 2. * ( _M[1][2] * x ).dot( y );

			Q2.coefficient(1,0) += 2. * _x.dot( _w[2] ) - 2. * _x.dot( _b[2] ) + 2. * _x.dot( _r[0] );
			Q2.coefficient(0,0) += 2. * x.dot( _w[2] ) - 2. * x.dot( _b[2] ) + 2. * x.dot( _r[0] );

			Q2.coefficient(0,1) += 2. * _y.dot( _w[1] ) - 2. * _y.dot( _b[1] ) + 2. * _y.dot( _r[1] );
			Q2.coefficient(0,0) += 2. * y.dot( _w[1] ) - 2. * y.dot( _b[1] ) + 2. * y.dot( _r[1] );
		}

		return Q2;
	}

	template< unsigned int Dim , typename Indexer >
	Polynomial::Polynomial2D< 2 > HierarchicalSystemEnergy< Dim , Indexer >::quadraticFit( const Eigen::VectorXd &x , const Eigen::VectorXd &y , unsigned int idx , bool setConstantTerm , unsigned int thread ) const
	{
		Polynomial::Polynomial2D< 2 > Q2;
		// Iterate over all other basis functions whose support overlaps the support of function idx=it.col()
		for( Eigen::InnerIterator it(_sI,idx) ; it ; ++it ) if( it.row()!=it.col() )
		{
			bool flip;
			// idx1: index of a wedge basis function one of whose components is idx
			size_t idx1 = wedges.index( std::make_pair( it.col() , it.row() ) , flip );
			for( Eigen::InnerIterator _it(_M[0][0],idx1) ; _it ; ++_it )
			{
				// The other wedge basis function whose support overlaps the support of function idx1
				std::pair< size_t , size_t > idx2 = wedges.index( _it.row() );
				double xy = x[ idx2.first ] * y[ idx2.second ] - x[ idx2.second ] * y[ idx2.first ];
				Q2.coefficient(1,0) += 2 * _it.value() * y[ it.row() ] * xy * ( flip ? -1. : 1. );
				Q2.coefficient(0,1) -= 2 * _it.value() * x[ it.row() ] * xy * ( flip ? -1. : 1. );

				if( idx2.first==idx )
				{
					Q2.coefficient(2,0) +=     _it.value() * y[ it.row() ] * y[ idx2.second ] * ( flip ? -1. : 1. );
					Q2.coefficient(0,2) +=     _it.value() * x[ it.row() ] * x[ idx2.second ] * ( flip ? -1. : 1. );
					Q2.coefficient(1,1) -= 2 * _it.value() * y[ it.row() ] * x[ idx2.second ] * ( flip ? -1. : 1. );
				}
				if( idx2.second==idx )
				{
					Q2.coefficient(2,0) -=     _it.value() * y[ it.row() ] * y[ idx2.first ] * ( flip ? -1. : 1. );
					Q2.coefficient(0,2) -=     _it.value() * x[ it.row() ] * x[ idx2.first ] * ( flip ? -1. : 1. );
					Q2.coefficient(1,1) += 2 * _it.value() * y[ it.row() ] * x[ idx2.first ] * ( flip ? -1. : 1. );
				}
			}
			Q2.coefficient(1,0) -= 2. * y[ it.row() ] * _b[0][idx1] * ( flip ? -1. : 1. );
			Q2.coefficient(0,1) += 2. * x[ it.row() ] * _b[0][idx1] * ( flip ? -1. : 1. );
		}

		for( Eigen::InnerIterator it(_R,idx) ; it ; ++it )
		{
			if( it.row()==idx )
			{
				Q2.coefficient(2,0) += it.value();
				Q2.coefficient(0,2) += it.value();
			}
			Q2.coefficient(1,0) += 2. * it.value() * x[ it.row() ];
			Q2.coefficient(0,1) += 2. * it.value() * y[ it.row() ];
		}

		if( setConstantTerm )
		{
			for( unsigned int i=0 ; i<wedges.functionNum() ; i++ )
			{
				std::pair< size_t , size_t > idx = wedges.index(i);
				double xy = x[ idx.first ] * y[ idx.second ] - x[ idx.second ] * y[ idx.first ];
				Q2.coefficient(0,0) += -2. * xy * _b[0][i];
			}
			Q2.coefficient(0,0) += _c[0] + _c[1];
		}

		if( _initFromFiner )
		{
			for( Eigen::InnerIterator _it(_M[0][2],idx) ; _it ; ++_it )
			{
				std::pair< size_t , size_t > idx2 = wedges.index( _it.row() );
				double xy = x[ idx2.first ] * y[ idx2.second ] - x[ idx2.second ] * y[ idx2.first ];
				Q2.coefficient(1,0) += 2. * _it.value() * xy;

				if( idx2.first==idx )
				{
					Q2.coefficient(2,0) += 2. * _it.value() * y[ idx2.second ];
					Q2.coefficient(1,1) -= 2. * _it.value() * x[ idx2.second ];
				}
				if( idx2.second==idx )
				{
					Q2.coefficient(2,0) -= 2. * _it.value() * y[ idx2.first ];
					Q2.coefficient(1,1) += 2. * _it.value() * x[ idx2.first ];
				}
			}

			for( Eigen::InnerIterator _it(_M[0][1],idx) ; _it ; ++_it )
			{
				std::pair< size_t , size_t > idx2 = wedges.index( _it.row() );
				double xy = x[ idx2.first ] * y[ idx2.second ] - x[ idx2.second ] * y[ idx2.first ];
				Q2.coefficient(0,1) += 2. * _it.value() * xy;

				if( idx2.first==idx )
				{
					Q2.coefficient(0,2) -= 2. * _it.value() * x[ idx2.second ];
					Q2.coefficient(1,1) += 2. * _it.value() * y[ idx2.second ];
				}
				if( idx2.second==idx )
				{
					Q2.coefficient(0,2) += 2. * _it.value() * x[ idx2.first ];
					Q2.coefficient(1,1) -= 2. * _it.value() * y[ idx2.first ];
				}
			}

			for( Eigen::InnerIterator _it(_M[2][2],idx) ; _it ; ++_it )
			{
				Q2.coefficient(1,0) += 2. * _it.value() * x[ _it.row() ];
				if( _it.row()==_it.col() ) Q2.coefficient(2,0) += _it.value();
			}

			for( Eigen::InnerIterator _it(_M[1][1],idx) ; _it ; ++_it )
			{
				Q2.coefficient(0,1) += 2. * _it.value() * y[ _it.row() ];
				if( _it.row()==_it.col() ) Q2.coefficient(0,2) += _it.value();
			}

			for( Eigen::InnerIterator _it(_M[2][1],idx) ; _it ; ++_it )
			{
				Q2.coefficient(0,1) += 2. * _it.value() * x[ _it.row() ];
				if( _it.row()==_it.col() ) Q2.coefficient(1,1) += _it.value();
			}
			for( Eigen::InnerIterator _it(_M[1][2],idx) ; _it ; ++_it )
			{
				Q2.coefficient(1,0) += 2. * _it.value() * y[ _it.row() ];
				if( _it.row()==_it.col() ) Q2.coefficient(1,1) += _it.value();
			}

			for( Eigen::InnerIterator it(_sI,idx) ; it ; ++it ) if( it.row()!=it.col() )
			{
				bool flip;
				size_t idx1 = wedges.index( std::make_pair( it.col() , it.row() ) , flip );
				for( Eigen::InnerIterator _it(_M[2][0],idx1) ; _it ; ++_it )
				{
					Q2.coefficient(1,0) += 2. * _it.value() * y[ it.row() ] * x[ _it.row() ] * ( flip ? -1. : 1. );
					Q2.coefficient(0,1) -= 2. * _it.value() * x[ it.row() ] * x[ _it.row() ] * ( flip ? -1. : 1. );
				}
				for( Eigen::InnerIterator _it(_M[1][0],idx1) ; _it ; ++_it )
				{
					Q2.coefficient(0,1) -= 2. * _it.value() * x[ it.row() ] * y[ _it.row() ] * ( flip ? -1. : 1. );
					Q2.coefficient(1,0) += 2. * _it.value() * y[ it.row() ] * y[ _it.row() ] * ( flip ? -1. : 1. );
				}

				Q2.coefficient(1,0) += 2. * y[ it.row() ] * _w[0][idx1] * ( flip ? -1. : 1. );
				Q2.coefficient(0,1) -= 2. * x[ it.row() ] * _w[0][idx1] * ( flip ? -1. : 1. );
			}

			Q2.coefficient(1,0) += 2. * _w[2][idx] - 2. * _b[2][idx] + 2. * _r[0][idx];
			Q2.coefficient(0,1) += 2. * _w[1][idx] - 2. * _b[1][idx] + 2. * _r[1][idx];

			if( setConstantTerm )
			{
				Eigen::VectorXd xy = wedges.wedge( x , y );
				Q2.coefficient(0,0) += 2. * ( _M[0][2] * x ).dot( xy );
				Q2.coefficient(0,0) += 2. * ( _M[0][1] * y ).dot( xy );
				Q2.coefficient(0,0) +=      ( _M[2][2] * x ).dot( x );
				Q2.coefficient(0,0) +=      ( _M[1][1] * y ).dot( y );
				Q2.coefficient(0,0) += 2. * xy.dot( _w[0] ) + 2. * ( _M[1][2] * x ).dot( y );
				Q2.coefficient(0,0) += 2. * x.dot( _w[2] ) - 2. * x.dot( _b[2] ) + 2. * x.dot( _r[0] );
				Q2.coefficient(0,0) += 2. * y.dot( _w[1] ) - 2. * y.dot( _b[1] ) + 2. * y.dot( _r[1] );
			}
		}

		return Q2;
	}

	template< unsigned int Dim , typename Indexer >
	typename SystemEnergy< Dim , Indexer >::QuadraticApproximation HierarchicalSystemEnergy< Dim , Indexer >::quadraticApproximation1( const Eigen::VectorXd & , const Eigen::VectorXd &y ) const
	{
		Eigen::SparseMatrix< double > Ly = -wedges.wedge( y );
		Eigen::SparseMatrix< double > Ly_t = Ly.transpose();
		typename SystemEnergy< Dim , Indexer >::QuadraticApproximation q;

		if( _initFromFiner )
		{
			q.q = Ly_t * _M[0][0] * Ly + Ly_t * _M[0][2] + _M[2][0] * Ly + _M[2][2] + _R;
			q.l = 2. * ( Ly_t * _M[0][1] * y + Ly_t * _w[0] - Ly_t * _b[0] + _M[2][1] * y + _w[2] - _b[2] + _r[0] );
			q.c = ( _M[1][1] * y ).dot( y ) + ( _R * y ).dot( y ) + 2. * y.dot( _w[1] - _b[1] + _r[1] ) + _c[0] + _c[1];
		}
		else
		{
			q.q = Ly_t * _M[0][0] * Ly + _R;
			q.l = 2. * ( - Ly_t * _b[0] );
			q.c = ( _R * y ).dot( y ) + _c[0] + _c[1];
		}
		return q;
	}

	template< unsigned int Dim , typename Indexer >
	typename SystemEnergy< Dim , Indexer >::QuadraticApproximation HierarchicalSystemEnergy< Dim , Indexer >::quadraticApproximation2( const Eigen::VectorXd &x , const Eigen::VectorXd & ) const
	{
		Eigen::SparseMatrix< double > Lx = wedges.wedge( x );
		Eigen::SparseMatrix< double > Lx_t = Lx.transpose();
		typename SystemEnergy< Dim , Indexer >::QuadraticApproximation q;

		if( _initFromFiner )
		{
			q.q = Lx_t * _M[0][0] * Lx  + Lx_t * _M[0][1] + _M[1][0] * Lx + _M[1][1] + _R;
			q.l = 2. * ( Lx_t * _M[0][2] * x + Lx_t * _w[0] - Lx_t * _b[0] + _M[1][2] * x + _w[1] - _b[1] + _r[1] );
			q.c = ( _M[2][2] * x ).dot( x ) + ( _R * x ).dot( x ) + 2. * x.dot( _w[2] - _b[2] + _r[0] ) + _c[0] + _c[1];
		}
		else
		{
			q.q = Lx_t * _M[0][0] * Lx + _R;
			q.l = 2. * ( - Lx_t * _b[0] );
			q.c = ( _R * x ).dot( x ) + _c[0] + _c[1];
		}
		return q;
	}

	template< unsigned int Dim , typename Indexer >
	std::pair< Eigen::VectorXd , Eigen::VectorXd > HierarchicalSystemEnergy< Dim , Indexer >::d( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const
	{
		std::pair< Eigen::VectorXd , Eigen::VectorXd > d;
		Eigen::VectorXd xy = wedges.wedge(x,y);
		Eigen::SparseMatrix< double > Ly_t = -wedges.wedge(y).transpose() , Lx_t = wedges.wedge(x).transpose();

		if( _initFromFiner )
		{
			d.first  = 2. * ( Ly_t * ( _M[0][0] * xy ) + 2. * _M[2][0] * xy + _M[2][2] * x + _R * x ) + 2. * ( Ly_t * _M[0][1] * y + Ly_t * _w[0] - Ly_t * _b[0] + _M[2][1] * y + _w[2] - _b[2] + _r[0] );
			d.second = 2. * ( Lx_t * ( _M[0][0] * xy ) + 2. * _M[1][0] * xy + _M[1][1] * y + _R * y ) + 2. * ( Lx_t * _M[0][2] * x + Lx_t * _w[0] - Lx_t * _b[0] + _M[1][2] * x + _w[1] - _b[1] + _r[1] );
		}
		else
		{
			d.first  = 2. * ( Ly_t * ( _M[0][0] * xy ) + _R * x ) + 2. * ( - Ly_t * _b[0] );
			d.second = 2. * ( Lx_t * ( _M[0][0] * xy ) + _R * y ) + 2. * ( - Lx_t * _b[0] );
		}
		return d;
	}
}
#endif // SYSTEM_ENERGY_INCLUDED