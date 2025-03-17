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

//tex:
// The space is spanned by $\{\omega_I\}$ where $I=(i,j)$ with $i\leq j$ (resp. $i < j$)
// index pairs of functions with overlapping support and $\omega_I = \frac{\nabla\phi_i\otimes\nabla\phi_j\pm\nabla\phi_j\otimes\nabla\phi_i}2$. //
// Given coefficients $X\in\mathbb{R}^N$, we denote by the associated function by $$\omega_X\equiv\sum_{I=1}^N x_I\cdot\omega_I.$$ 
template< unsigned int Dim , bool Sym >
struct ProductFunctions
{
	template< typename T , unsigned int Rank , unsigned int Radius >
	using IntegrationStencil = typename ScalarFunctions< Dim >::template IntegrationStencil< T , 2*Rank , Radius >;

	template< typename T >
	struct FullIntegrationStencil
	{
		struct Entry
		{
			Index< Dim > g1 , g2;	// All pairs of overlapping functions (except g1=g2 if not symmetric)
			long long dg1 , dg2;	// The relative (1D) index offsets
			unsigned int _g1, _g2; 
			T value;
			Entry(Index< Dim > g1, Index< Dim > g2, size_t dg1, size_t dg2, T value) : g1(g1), g2(g2), dg1(dg1), dg2(dg2), value(value) 
			{
				for (unsigned int d = 0; d < Dim; d++) g1[d]++, g2[d]++;
				_g1 = (unsigned int)ScalarFunctions< Dim >::FunctionIndex(g1, 2);
				_g2 = (unsigned int)ScalarFunctions< Dim >::FunctionIndex(g2, 2);
			}
		};

		struct Row
		{
			Index< Dim > f2;	// All functions with overlapping support (except f1=f2 if not symmetric)
			long long df2;		// The relative (1D) index offset
			std::vector< Entry > entries;
		};

		FullIntegrationStencil( const IntegrationStencil< T , 2 , 0 > &stencil , unsigned int res );
		const std::vector< Row > &rows( Index< Dim > f1 ) const;

	protected:
		unsigned int _res;
		std::vector< std::vector< Row > > _rows;
	};

	// Function for combining a values from a pair of scalar functions to get the coefficient of the product
	static double Coefficient( std::pair< double , double > f1 , std::pair< double , double > f2 , size_t i1 , size_t i2 );


	ProductFunctions( unsigned int resolution );
	size_t resolution( void ) const { return _r; }

	// The number of product functions
	size_t functionNum( void ) const { return _matrixInfo.entries(false); }

	// Returns the index associated to a pair vertices and sets the flag if the orientation is reversed
	// Throws an exception if the pair does not index a basis
	size_t index( Index< Dim > F1 , Index< Dim > F2 , bool &flip ) const;
	size_t index( std::pair< size_t , size_t > idx , bool &flip ) const { return index( ScalarFunctions< Dim >::FunctionIndex(idx.first,_r) , ScalarFunctions< Dim >::FunctionIndex(idx.second,_r) , flip ); }

	// Sets the index associated to a pair vertices and sets the flag if the orientation is reversed
	// Returns false if the pair does not index a basis
	bool setIndex( Index< Dim > F1 , Index< Dim > F2 , size_t &i , bool &flip ) const;
	bool setIndex( std::pair< size_t , size_t > idx , size_t &i , bool &flip ) const { return setIndex( ScalarFunctions< Dim >::FunctionIndex(idx.first,_r) , ScalarFunctions< Dim >::FunctionIndex(idx.second,_r) , i , flip ); }

	// Returns the index pairs associated with the functions
	std::vector< std::pair< Index< Dim > , Index< Dim > > > indices( void ) const;

	//////////////
	// Products //
	//////////////

#	//tex: Computes the product coefficients, should satisfy $$\omega_{x \times y} = \nabla\phi_x\times\nabla\phi_y.$$
	Eigen::VectorXd product( const Eigen::VectorXd &x , const Eigen::VectorXd &y ) const;

	//tex: Computes the transpose of the product operator associated to $x\in\mathbb{R}^n$ applied to $Z\in\mathbb{R}^N$:
	// $$(x\times\cdot)^\top\cdot Z.$$
	Eigen::VectorXd productTranspose( const Eigen::VectorXd &x , const Eigen::VectorXd &Z ) const;

	// Returns the function corresponding to g -> df * dg
	Eigen::SparseMatrix< double > product( const Eigen::VectorXd &f ) const;


	////////////////////////////////////////////
	// Evaluation and Monte-Carlo integration //
	////////////////////////////////////////////

	// Evaluation of a function's value, represented in the product function basis
	SquareMatrix< double , Dim , Sym > value( const Eigen::VectorXd &xy , Point< double , Dim > p ) const;

	// Estimates the dot-product of two scalar functions
	template< typename F1 /* = std::function< SquareMatrix< double , Dim , Sym > ( Point< double , Dim > ) > */ , typename F2 /* = std::function< SquareMatrix< double , Dim , Sym > ( Point< double , Dim > ) > */ >
	static double TensorDotProduct( F1 f1 , F2 f2 , unsigned int samplesPerDimension ){ return Basis< Dim >::template Integral< double >( [&]( Point< double , Dim > p ){ return SquareMatrix< double , Dim , Sym >::Dot( f1(p) , f2(p) ); } , samplesPerDimension ); }


	//////////////////////////////////////////
	// Single function integration stencils //
	//////////////////////////////////////////

	static IntegrationStencil< SquareMatrix< double , Dim , Sym > , 1 , 0 > ValueStencil( unsigned int r );

	////////////////////////////////////////////
	// Pair of functions integration stencils //
	////////////////////////////////////////////

	static IntegrationStencil< double , 2 , 0 > MassStencil( unsigned int r );

	////////////////////////
	// Stencil evaluation //
	////////////////////////

	// Evaluates the stencil on a pair of functions
	template< typename T >
	T operator()( const     IntegrationStencil< T , 2 , 0 > &stencil , const Eigen::VectorXd &xy1 , const Eigen::VectorXd &xy2 ) const;

	// Evaluates the stencil on a pair of functions
	template< typename T >
	T operator()( const FullIntegrationStencil< T > &stencil , const Eigen::VectorXd &xy1 , const Eigen::VectorXd &xy2 ) const;

	// Evaluates the stencil on the pair of functions generated as the product of pairs of functions
	template< typename T , typename Indexer /* = Hat::BaseIndex< Dim > */ >
	T operator()( const Indexer &indexer , const FullIntegrationStencil< T > &stencil , const Eigen::VectorXd &x1 , const Eigen::VectorXd &x2 , const Eigen::VectorXd &y1 , const Eigen::VectorXd &y2 ) const;

	// Transforms a stencil into a linear operator
	Eigen::SparseMatrix< double > systemMatrix( IntegrationStencil< double , 2 , 0 > stencil ) const;


	//////////////////////////
	// Dual representations //
	//////////////////////////

	// Integrates the basis functions against a piecewise constant tensor field
	template< typename PiecewiseConstantTensorField /* = std::function< SquareMatrix< double , Dim , Sym > ) ( size_t e ) > */ >
	Eigen::VectorXd valueDual( PiecewiseConstantTensorField T ) const;

	// Integrates the basis functions against the tensor field with the given coefficients
	Eigen::VectorXd valueDual( const Eigen::VectorXd &xy ) const;

	//////////////////////////////
	// Bilinear representations //
	//////////////////////////////
	
	// Returns the mass matrix for the product functions
	Eigen::SparseMatrix< double > mass( void ) const { return systemMatrix( MassStencil( _r ) ); }




	Eigen::SparseMatrix< double > prolongation( void ) const;

	//tex: Given a vector $Z\in\mathbb{R}^N$, if we denote by $\mathbf{A}_Z\in\mathbb{R}^{n\times n}$ the representation of $Z$ as a matrix,
	// then $y^\top\cdot\mathbf{A}_Z\cdot x = Z^\top\cdot(x\times y)$.
	Eigen::SparseMatrix< double > toMatrix( const Eigen::VectorXd &xy ) const;
	Eigen::VectorXd toVector( const Eigen::SparseMatrix< double > &XY ) const;

	struct UnitTest
	{
	protected:
		ScalarFunctions< Dim > _scalars;
		ProductFunctions< Dim , Sym > _products;
		unsigned int _N , _numSamples;
	public:
		UnitTest( unsigned int resolution , unsigned int N , unsigned int numSamples ) : _scalars(resolution) , _products(resolution) , _N(N) , _numSamples(numSamples){}
		void operator()( void ) const;

		// Check that the product of the evaluation of the gradients matches the evaluation of the products
		// ScalarFunctions< Dim >::gradient -> ProductFunctions< Dim , Sym >::value + ProductFunctions< Dim , Sym >::product
		void evaluation( void ) const;

		// Check that the vector-vector a vector-matrix representations of the product agree
		// ProductFunctions< Dim , Sym >::product <-> ProductFunctions< Dim , Sym >::product
		void product( void ) const;

		// Check that the stencil values match what would be obtained through discrete/monte-carlo integration
		// ProductFunctions::ValueStencil <-> ElementFunction::gradient
		void integrationStencils( void ) const;

		// The evaluations of mass should agree with the dual
		// ProductFunctions::mass <-> ProductFunctions::dual <->  ProductFunctions::operator() <-> ProductFunctions::operator() <-> ProductFunctions::operator()
		void massAndDual( void ) const;

		// The matrix to vector and vector matrix should be consistent
		void matrixVector( void ) const;
	};

protected:

	unsigned int _r;

	typename ScalarFunctions< Dim >::template MatrixInfo< 1 , Sym > _matrixInfo;
};

#include "Hat.products.inl"
