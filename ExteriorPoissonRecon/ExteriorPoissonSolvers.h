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

#ifndef POISSON_CURVE_RECON_INCLUDED
#define POISSON_CURVE_RECON_INCLUDED


#include <vector>
#include <iomanip>
#include "Misha/Geometry.h"
#include "Misha/Miscellany.h"
#include "Misha/Poly34.h"
#include "Misha/Polynomial.h"
#include "Misha/MultiThreading.h"
#include "../Include/Hat.h"
#include "../Include/MG.h"
#include "../Include/SystemEnergy.h"

namespace MishaK
{
	namespace ExteriorPoissonSolvers
	{
		enum SolverType
		{
			DIRECT ,
			SINGLE_LEVEL ,
			MULTIGRID
		};
		static const unsigned int SolverTypeCount = 3;
		static const std::string SolverTypeNames[] = { std::string( "direct" ) , std::string( "single-level" ) , std::string( "multigrid" ) };

		enum UpdateType
		{
			ALTERNATING ,
			PAIRED_GAUSS_SEIDEL ,
			GRADIENT_DESCENT
		};
		static const unsigned int UpdateTypeCount = 3;
		static const std::string UpdateTypeNames[] = { std::string( "alternating" ) , std::string( "paired Gauss-Seidel" ) , std::string( "gradient-descent" ) };

		// Given the current estimate of the solution (x,y), updates the estimate of the solution by alternately locking one of the veriables and solving for the minimizer of the quadratic equation in the second
		template< unsigned int Dim , typename SystemEnergyType , typename ProlongationFunctor , typename MCIndicesFunctor >
		void AlternatingSolve( unsigned int depth , SolverType sType , ProlongationFunctor P , MCIndicesFunctor MCI , const SystemEnergyType &energy , Eigen::VectorXd &x , Eigen::VectorXd &y , unsigned int iters , unsigned int vCycles , unsigned int gsIters , bool verbose );

		// Given the current estimate of the solution (x,y), updates the estimate of the solution by performing gradient descent on the energy
		template< unsigned int Dim , typename SystemEnergyType , typename ProlongationFunctor , typename MCIndicesFunctor >
		void GradientDescent( unsigned int depth , bool useMG , ProlongationFunctor P , MCIndicesFunctor MCI , const SystemEnergyType &energy , const Eigen::SparseMatrix< double > &sM , Eigen::VectorXd &x , Eigen::VectorXd &y , unsigned int iters , unsigned int vCycles , unsigned int gsIters , unsigned int newtonIters , bool verbose );

		// Given the current estimate of the solution (x,y), updates the estimate of the solution by performing Gauss-Seidel relaxation on pairs of variables at a time
		template< unsigned int Dim , typename SystemEnergyType >
		void PairedGaussSeidelSolve( unsigned int depth , const std::vector< std::vector< unsigned int > > &mcIndices , const SystemEnergyType &energy , Eigen::VectorXd &x , Eigen::VectorXd &y , unsigned int iters , bool verbose );

		// Given the current estimate of the solution (x,y), updates the estimate of the solution by performing Gauss-Seidel relaxation on one variable a time
		template< unsigned int Dim , typename SystemEnergyType >
		void SeparateGaussSeidelSolve( unsigned int depth , const std::vector< std::vector< unsigned int > > &mcIndices , const SystemEnergyType &energy , Eigen::VectorXd &x , Eigen::VectorXd &y , unsigned int iters , bool verbose );

		/////////////
		// Solvers //
		/////////////

		template< unsigned int Dim , typename SystemEnergyType , typename ProlongationFunctor , typename MCIndicesFunctor >
		void AlternatingSolve( unsigned int depth , SolverType sType , ProlongationFunctor P , MCIndicesFunctor MCI , const SystemEnergyType &energy , Eigen::VectorXd &x , Eigen::VectorXd &y , unsigned int iters , unsigned int vCycles , unsigned int gsIters , bool verbose )
		{
			static_assert( std::is_base_of< SystemEnergy< Dim , typename SystemEnergyType::Indexer > , SystemEnergyType >::value , "[ERROR] Poorly formed SystemEnergyType" );
			Miscellany::PerformanceMeter pMeter;
			Hat::ScalarFunctions< Dim > scalars( 1<<depth );

			double err = 0.;
			if( verbose ) err = energy( x , y );
			for( unsigned int iter=0 ; iter<iters ; iter++ )
			{
				double _err = 0;
				if( verbose ) _err = energy( x , y );
				// Fixing x, minimize
				{
					typename SystemEnergy< Dim , typename SystemEnergyType::Indexer >::QuadraticApproximation q = energy.quadraticApproximation2( x , y );
					q.q *= 2.;
					q.l = -q.l;
					double _err;
					if( verbose ) _err = ( q.q * y - q.l ).squaredNorm();
					if( sType==SINGLE_LEVEL )
					{
						GSRelaxer gsRelaxer( q.q );
						for( unsigned int i=0 ; i<gsIters*vCycles ; i++ ) gsRelaxer( q.l , y , MCI(depth) );
					}
					else
					{
						MGSolver mgSolver( q.q , P , sType==MULTIGRID ? depth : 0 );
						for( unsigned int v=0 ; v<vCycles ; v++ ) mgSolver.vCycle( q.l , y , gsIters , MCI );
					}
					if( verbose ) std::cout << "\t\t\t" << _err << " -> " << ( q.q * y - q.l ).squaredNorm() << std::endl;
				}
				// Fixing y, minimize
				{
					typename SystemEnergy< Dim , typename SystemEnergyType::Indexer >::QuadraticApproximation q = energy.quadraticApproximation1( x , y );
					q.q *= 2.;
					q.l = -q.l;
					double _err;
					if( verbose ) _err = ( q.q * x - q.l ).squaredNorm();
					if( sType==SINGLE_LEVEL )
					{
						GSRelaxer gsRelaxer( q.q );
						for( unsigned int i=0 ; i<gsIters*vCycles ; i++ ) gsRelaxer( q.l , x , MCI(depth) );
					}
					else
					{
						MGSolver mgSolver( q.q , P , sType==MULTIGRID ? depth : 0 );
						for( unsigned int v=0 ; v<vCycles ; v++ ) mgSolver.vCycle( q.l , x , gsIters , MCI );
					}
					if( verbose ) std::cout << "\t\t\t" << _err << " -> " << ( q.q * x - q.l ).squaredNorm() << std::endl;
				}
				if( verbose )
				{
					std::pair< double , double > err = energy.energies( x , y );
					std::cout << "\t\tError[" << (1<<depth) << "][" << iter << "]: " << _err << " -> " << err.first << " + " << err.second << std::endl;
				}
			}
			if( verbose )
			{
				std::pair< double , double > _err = energy.energies( x , y );
				std::stringstream sStream1 , sStream2;
				{
					Miscellany::StreamFloatPrecision sfp( sStream2 , 5 );
					sStream1 << "Error[" << std::setw(2) << (1<<depth) << "]";
					sStream2 << std::setw(9) << err << " -> " << ( _err.first+_err.second ) << " = " << _err.first << " + " << _err.second;
				}
				std::cout << pMeter( sStream1.str() , sStream2.str() ) << std::endl;
			}
		}

		template< unsigned int Dim , typename SystemEnergyType >
		void PairedGaussSeidelSolve( unsigned int depth , const std::vector< std::vector< unsigned int > > &mcIndices , const SystemEnergyType &energy , Eigen::VectorXd &x , Eigen::VectorXd &y , unsigned int iters , bool verbose )
		{
			static_assert( std::is_base_of< SystemEnergy< Dim , typename SystemEnergyType::Indexer > , SystemEnergyType >::value , "[ERROR] Poorly formed SystemEnergyType" );
			Miscellany::PerformanceMeter pMeter;

			double err = 0.;
			if( verbose ) err = energy( x , y );
			for( unsigned int iter=0 ; iter<iters ; iter++ )
			{
				for( unsigned int i=0 ; i<mcIndices.size() ; i++ )
					ThreadPool::ParallelFor
					(
						0 , mcIndices[i].size() ,
						[&]( unsigned int t , size_t j )
				{
					unsigned int idx = mcIndices[i][j];
					Polynomial::Polynomial2D< 2 > Q = energy.quadraticFit( x , y , idx , false , t );

					SquareMatrix< double , 2 > H = Q.hessian( Point< double , 2 >() );
					Point< double , 2 > grad = Q.gradient( Point< double , 2 >() );
					if( fabs( H.determinant() )>1e-20 )
					{
						Point< double , 2 > p = - H.inverse() * grad / 2.;
						x[idx] += p[0] , y[idx] += p[1];
					}
				}
					);
			}
			if( verbose )
			{
				std::pair< double , double > _err = energy.energies( x , y );
				std::stringstream sStream1 , sStream2;
				{
					Miscellany::StreamFloatPrecision sfp( sStream2 , 5 );
					sStream1 << "Error[" << std::setw(2) << (1<<depth) << "]";
					sStream2 << std::setw(9) << err << " -> " << ( _err.first+_err.second ) << " = " << _err.first << " + " << _err.second;
				}
				std::cout << pMeter( sStream1.str() , sStream2.str() ) << std::endl;
			}
		}

		template< unsigned int Dim , typename SystemEnergyType >
		void SeparateGaussSeidelSolve( unsigned int depth , const std::vector< std::vector< unsigned int > > &mcIndices , const SystemEnergyType &energy , Eigen::VectorXd &x , Eigen::VectorXd &y , unsigned int iters , bool verbose )
		{
			static_assert( std::is_base_of< SystemEnergy< Dim , typename SystemEnergyType::Indexer > , SystemEnergyType >::value , "[ERROR] Poorly formed SystemEnergyType" );
			Miscellany::PerformanceMeter pMeter;

			double err = 0.;
			if( verbose ) err = energy( x , y );
			for( unsigned int iter=0 ; iter<iters ; iter++ )
			{
				for( unsigned int i=0 ; i<mcIndices.size() ; i++ )
					ThreadPool::ParallelFor
					(
						0 , mcIndices[i].size() ,
						[&]( unsigned int t , size_t j )
				{
					unsigned int idx = mcIndices[i][j];

					Polynomial::Polynomial2D< 2 > Q = energy.quadraticFit( x , y , idx , false , t );

					SquareMatrix< double , 2 > H = Q.hessian( Point< double , 2 >() );
					Point< double , 2 > grad = Q.gradient( Point< double , 2 >() );
					// E(p) = p^t * H * p + p^t * grad
					// Fixing p[1], this gives:
					// E(p[0]) = p[0]^2 * H(0,0) + ( H(1,0) + H(0,1) ) * p[0] * p[1] + p[0] * grad[0]
					//  => 0 = 2 * H(0,0) * p[0] + ( H(1,0) + H(0,1) ) * p[1] + grad[0] 
					//  => p[0] = - [ ( H(1,0) + H(0,1) ) * p[1] + grad[0]  ] / ( 2 * H(0,0) )

					// Fixing p[0], this gives:
					// E(p[1]) = p[1]^2 * H(1,1) + ( H(1,0) + H(0,1) ) * p[0] * p[1] + p[1] * grad[1]
					//  => 0 = 2 * H(1,1) * p[1] + ( H(1,0) + H(0,1) ) * p[0] + grad[1] 
					//  => p[1] = - [ ( H(1,0) + H(0,1) ) * p[0] + grad[1]  ] / ( 2 * H(1,1) )

					Point< double , 2 > p;
					if( fabs( H(0,0)>1e-20 ) ) p[0] = - ( ( H(1,0) + H(0,1) )*p[1] + grad[0] ) / ( 2 * H(0,0) );
					if( fabs( H(1,1)>1e-20 ) ) p[1] = - ( ( H(1,0) + H(0,1) )*p[0] + grad[1] ) / ( 2 * H(1,1) );

					x[idx] += p[0] , y[idx] += p[1];
				}
					);
			}
			if( verbose )
			{
				std::pair< double , double > _err = energy.energies( x , y );
				std::stringstream sStream1 , sStream2;
				{
					Miscellany::StreamFloatPrecision sfp( sStream2 , 5 );
					sStream1 << "Error[" << std::setw(2) << (1<<depth) << "]";
					sStream2 << std::setw(9) << err << " -> " << ( _err.first+_err.second ) << " = " << _err.first << " + " << _err.second;
				}
				std::cout << pMeter( sStream1.str() , sStream2.str() ) << std::endl;
			}
		}

		template< unsigned int Dim , typename SystemEnergyType , typename ProlongationFunctor , typename MCIndicesFunctor >
		void GradientDescent( unsigned int depth , bool useMG , ProlongationFunctor P , MCIndicesFunctor MCI , const SystemEnergyType &energy , const Eigen::SparseMatrix< double > &sM , Eigen::VectorXd &x , Eigen::VectorXd &y , unsigned int iters , unsigned int vCycles , unsigned int gsIters , unsigned int newtonIters , bool verbose )
		{
			static_assert( std::is_base_of< SystemEnergy< Dim , typename SystemEnergyType::Indexer > , SystemEnergyType >::value , "[ERROR] Poorly formed SystemEnergyType" );

			Miscellany::PerformanceMeter pMeter;
			Hat::ScalarFunctions< Dim > scalars( 1<<depth );
			MGSolver mgSolver( sM , P , useMG ? depth : 0 );

			double err = 0;
			if( verbose ) err = energy( x , y );
			for( unsigned int iter=0 ; iter<iters ; iter++ )
			{
				// E(x,y) = || dx ^ dy - b ||^2 + x^t * R * x + y^t * R * y
				// Taking the partial derivative w.r.t. x:
				//		d_xE = 2 < dx ^ dy - b , d_x ( dx ^ dy - b ) > + 2 * R * x
				//		     = 2 < dx ^ dy - b , d_x ( dx ^ dy ) > + 2 * R * x
				//		     = 2 < dx ^ dy - b , d_x Ly( x ) > + 2 * R * x
				//		     = 2 Ly^t * ( wM ( dx ^ dy ) - b ) + 2 * R * x
				double _err = 0;
				if( verbose ) _err = energy( x , y );
				std::pair< Eigen::VectorXd , Eigen::VectorXd > d = energy.d( x , y );
				Eigen::VectorXd nabla_x , nabla_y;
				nabla_x.setZero( d.first .size() );
				nabla_y.setZero( d.second.size() );
				double _errX = 0 , _errY = 0;
				if( verbose ) _errX = ( sM * nabla_x - d.first ).squaredNorm() , _errY = ( sM * nabla_y - d.second ).squaredNorm();
				for( unsigned int v=0 ; v<vCycles ; v++ ) mgSolver.vCycle( d.first  , nabla_x , gsIters , MCI );
				for( unsigned int v=0 ; v<vCycles ; v++ ) mgSolver.vCycle( d.second , nabla_y , gsIters , MCI );
				if( verbose ) std::cout << "\t\t\t" << _errX << " , " << _errY << " -> " << ( sM * nabla_x - d.first ).squaredNorm() << " , " << ( sM * nabla_y - d.second ).squaredNorm() << std::endl;

				// Rescale so that the coefficients of the polynomial do not explode
				double scale = sqrt( nabla_x.squaredNorm() + nabla_y.squaredNorm() );
				nabla_x /= scale , nabla_y /= scale;

				if( newtonIters )
				{
					Point< double , 2 > p = energy.newtonUpdate( x , y , nabla_x , nabla_y , newtonIters );
					x += nabla_x * p[0];
					y += nabla_y * p[1];
				}
				else
				{
					double s = energy.stepSize( x , y , nabla_x , nabla_y );
					x += s * nabla_x;
					y += s * nabla_y;
				}

				if( verbose )
				{
					std::pair< double , double > err = energy.energies( x , y );
					std::cout << "\t\tError[" << (1<<depth) << "][" << iter << "]: " << _err << " -> " << err.first << " + " << err.second << std::endl;
				}
			}
			if( verbose )
			{
				std::pair< double , double > _err = energy.energies( x , y );
				std::stringstream sStream1 , sStream2;
				{
					Miscellany::StreamFloatPrecision sfp( sStream2 , 5 );
					sStream1 << "Error[" << std::setw(2) << (1<<depth) << "]";
					sStream2 << std::setw(9) << err << " -> " << ( _err.first+_err.second ) << " = " << _err.first << " + " << _err.second;
				}
				std::cout << pMeter( sStream1.str() , sStream2.str() ) << std::endl;
			}
		}
	}
}
#endif // POISSON_CURVE_RECON_INCLUDED