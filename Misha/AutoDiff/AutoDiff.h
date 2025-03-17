/*
Copyright (c) 2020, Michael Kazhdan
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
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMIRTED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef AUTO_DIFF_INCLUDED
#define AUTO_DIFF_INCLUDED

#define NEW_AUTO_DIFF_CODE

#include <iostream>
#include "Tensors.h"
#include "Misha/Exceptions.h"

namespace MishaK
{
	namespace AutoDiff
	{
		template< typename T > bool constexpr IsScalar( void );

		////////////////////////////
		////////////////////////////
		//// Short declarations ////
		////////////////////////////
		////////////////////////////

		//////////////
		// Function //
		//////////////
		// A general class describing a function taking in a tensor and outputting a tensor
		// (Parameters OutPack and InPack are assumed to be of variadic UIntPack<...> type giving the
		// dimensions along the different ranks of the input/output tensors.)
		template< typename OutPack , typename InPack , typename F > struct Function;

		/////////////////////
		// Basic Functions //
		/////////////////////

		// A Function that is constant in its input
		// InPack: a UIntPacK< ... > describing the dimensions of the input tensor
		// OutPack: a UIntPacK< ... > describing the dimensions of the output tensor
		template< typename OutPack , typename InPack > struct Constant;

		// A Function that is linear in its input
		// InPack: a UIntPacK< ... > describing the dimensions of the input tensor
		// OutPack: a UIntPacK< ... > describing the dimensions of the output tensor
		template< typename OutPack , typename InPack > struct Linear;

		// The identity Function
		// InOutPack: a UIntPack< ... > describing the dimensions of the input/output tensor
		template< typename InOutPack > struct Identity;


		////////////////////////////////////////
		// Functions that eat other Functions //
		////////////////////////////////////////

		// A function returning the negation of a Function
		// Output type derives from Function< F::OutPack , F::InPack >
		template< typename F > auto operator - ( const F &f );

		// A function returning the sum of two Function's (and specializations)
		// Assumes:
		//		F1::InPack==F2::InPack
		//		F1::OutPack==F2::OutPack
		// Output type derives from Function< F1::OutPack , F1::InPack >
		template< typename F1 , typename F2 > auto operator + ( const F1 &f1 , const F2 &f2 );

		// A function returning the difference of two Function's (and specializations)
		// Assumes:
		//		F1::InPack==F2::InPack
		//		F1::OutPack==F2::OutPack
		// Output type derives from Function< F1::OutPack , F1::InPack >
		template< typename F1 , typename F2 > auto operator - ( const F1 &f1 , const F2 &f2 );

		// A function returning the outer product of a Function with a Function (and specialization)
		// Assumes:
		//		F1::InPack==F2::InPack
		// Output type derives from Function< Concatenation< F1::OutPack , F2::OutPack > , F1::InPack >
		template< typename F1 , typename F2 > auto operator * ( const F1 &f , const F2 &f2 );

		// A function returning the quotient of a Function by a scalar-valued Function
		// Assumes:
		//		F1::InPack==F2::InPack
		//		F2::OutPack==UIntPack<>
		// Output type derives from Function< F1::OutPack , F1::InPack >
		template< typename F1 , typename F2 > auto operator / ( const F1 &f1 , const F2 &f2 );

		// A function returning the permutation a Function's output
		// Assumes:
		//		PermutationPack is a permutation
		//		PermutationPack::Size==F::OutPack::Size
		// Output type derives from Function< Permutation< F::OutPack , PermutationPack > >
		template< typename PermutationPack , typename F > auto Permutation( const F &f );

		// A function returning the contracted outer product of two Function's
		// Assumes:
		//		F1::InPack==F2::InPack
		//		F1::OutPack==Concatenation< Pack1 , Pack2 >
		//		F2::OutPack==Concatenation< Pack2 , Pack3 >
		//		Pack2::Size = I
		// Output type derives from Function< Concatenation< Pack1 , Pack3 > , F1::InPack >
		template< unsigned int I , typename F1 , typename F2 > auto ContractedOuterProduct( const F1 &f1 , const F2 &f2 );

		// A function returning the contraction of a Function
		// Assumes:
		//		F::OutPack::Size>=2
		// Output type derives from Function< F::OutPack::Remove(I1,I2) , F::InPack >
		template< unsigned int I1 , unsigned int I2 , typename F > auto Contraction( const F &f );

		// A function returning the composition two Function's
		// Assumes:
		//		F2::InPack==F1::OutPack
		// Output type derives from Function< F2::OutPack , F1::InPack >
		template< typename F1 , typename F2 > auto Composition( const F1 &f1 , const F2 &f2 );

#ifdef NEW_AUTO_DIFF_CODE
		// A function returning the concatenation of multiple functions
		// Assumes:
		//		F::InPack==Fs::InPack
		//		F::OutPack==Fs::OutPack
		// Output type derives from Function< ParameterPack::Concatenation< ParameterPack::UIntPack< sizeof ... (Fs) + 1 , F::OutPack > , F::InPack >
		template< typename F , typename ... Fs > auto Concatenation( const F &f , const Fs & ... fs );

		// A function returning the first function if the condition is met and the second otherwise
		// Assumes:
		//		F1::InPack==F2::InPack
		//		F1::OutPack==F2::OutPack
		//		ConditionFunctor = std::function< bool ( Tensory< F1::InPack ) >
		// Output type derives from Function< F1::OutPack , F1::InPack >
		template< typename ConditionFunctor , typename F1 , typename F2 > auto Conditional( ConditionFunctor c , const F1 &f1 , const F2 &f2 );
#endif // NEW_AUTO_DIFF_CODE

		// A function extracting a sub-tensor whose first I coefficients are given by indices
		// Assumes:
		//		I<=F::OutPack::Size
		//		indices[i]<F::OutPack[i]
		// Output type derives from Function< ParameterPack::Partition< I , F::OutPack >::Second , F::InPack >
		template< unsigned int I , typename F > auto Extract( const unsigned int indices[/*I*/] , const F &f );

		// A function extracting a single coefficient of the output
		// Assumes:
		//		sizeof...(indices)==F::OutPack::Size
		//		indices[i]<F::OutPack[i]
		// Output type derives from Function< UIntPack<> , F::InPack >
		template< typename F > auto Coefficient( const unsigned int indices[/*F::OutPack::Size*/] , const F &f );

		// A function returning the transpose of the output of a Function
		// Output type derives from Function< F::OutPack::Transpose , F::InPack >
		template< typename F > auto Transpose( const F &f );

		// A function returning the square norm of the output of a Function
		// Output type derives from Function< UIntPack<> , F::InPack >
		template< typename F > auto SquareNorm( const F &f );

		// A function returning the determinant of the output of a Function (assumed to return a square 2-tensor)
		// Assumes:
		//		F::OutPack==UIntPack< Dim , Dim >
		// Output type derives from Function< UIntPack<> , F::InPack >
		template< typename F > auto Determinant( const F &f );

		// A function returning the cofactor of the output of a Function (assumed to return a square 2-tensor)
		// Assumes:
		//		F::OutPack==UIntPack< Dim , Dim >
		// Output type derives from Function< UIntPack<> , F::InPack >
		template< typename F > auto Cofactor( const F &f );

		// A function returning the adjugate of the output of a function (assumed to return a square 2-tensor)
		// Assumes:
		//		F::OutPack==UIntPack< Dim , Dim >
		// Output type derives from Function< UIntPack<> , F::InPack >
		template< typename F > auto Adjugate( const F &f );

		// A function returning the inverse of the output of a function (assumed to return a square 2-tensor)
		// Assumes:
		//		F::OutPack==UIntPack< Dim , Dim >
		// Output type derives from Function< UIntPack<> , F::InPack >
		template< typename F > auto Inverse( const F &f );

		// A function returning the cross-product of the columns of the output of a Function (assumed to return a 2-tensor with one more row than columns)
		// Assumes:
		//		F::OutPack==UIntPack< Dim , Dim-1 >
		// Output type derives from Function< UIntPack<Dim> , F::InPack >
		template< typename F > auto CrossProduct( const F &f );

		// Some common transcendental Function's
		// Assumes:
		//		F::OutPack==UIntPack<>
		// Output type derives from Function< UIntPack<> , F::InPack >
		template< typename F > auto Pow( const F &f , double e );
		template< typename F > auto Exp( const F &f );
		template< typename F > auto Log( const F &f );
		template< typename F > auto Sin( const F &f );
		template< typename F > auto Cos( const F &f );
		template< typename F > auto Sinh( const F &f );
		template< typename F > auto Cosh( const F &f );

		// A function returning the discrete derivative of a Function at a particular input
		template< typename F > auto DiscreteDerivative( const F &f , const Tensor< typename F::InPack > &in , double eps );

		////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////
		//// Full declarations (including helper functions) ////
		////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////

		template< typename T >
		bool constexpr IsScalar( void )
		{
			if constexpr( std::is_arithmetic_v< T > || std::is_base_of< Tensor< ParameterPack::UIntPack<> > , T >::value ) return true;
			else return false;
		}

		// Some combinatoric functions that will be useful
		template< unsigned int D > struct Factorial      { static const unsigned int Value = Factorial< D-1 >::Value * D; };
		template<>                 struct Factorial< 0 > { static const unsigned int Value = 1; };

		template< unsigned int D , unsigned int K > struct Choose;
		template< unsigned int D , unsigned int K > struct Choose         { static const unsigned int Value = ( Choose< D-1 , K-1 >::Value * D ) / K; };
		template< unsigned int D >                  struct Choose< D , 0 >{ static const unsigned int Value = 1; };

		template< bool B , bool ... Bs >
		constexpr bool AND( void )
		{
			if constexpr( sizeof...(Bs)==0 ) return B;
			else return B && AND< Bs... >();
		}
		template< bool B , bool ... Bs >
		constexpr bool OR( void )
		{
			if constexpr( sizeof...(Bs)==0 ) return B;
			else return B || OR< Bs... >();
		}

		// A class for describing a function
		template< unsigned int ... OutDims , unsigned int ... InDims , typename F >
		struct Function< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > , F >
		{
			// Type of the tensor taken as input
			typedef ParameterPack::UIntPack< InDims ... > InPack;
			// Type of the tensor returned as output
			typedef ParameterPack::UIntPack< OutDims ... > OutPack;

			// [NOTE] There are three ways the function call operator can be invoked with a single argument:
			//	1. [Composition] Where the argument is a function, implying composition
			//	2. [Evaluation] Where the argument is a scalar and the input is a zero-order tensor
			//	3. [Evaluation] Where the argument is a tensor of the same type as the input
			template< typename V > auto operator()( const V &v ) const;
		};

		// A class for describing a constant function
		template< unsigned int ... OutDims , unsigned int ... InDims >
		struct Constant< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > > : public Function< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > , Constant< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > > >
		{
			typedef Function< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > , Constant > _Function;

			Constant( void ){}
			Constant( const Tensor< ParameterPack::UIntPack< OutDims ... > > &c ) : _c(c){}

			auto value( const Tensor< ParameterPack::UIntPack< InDims ... > > &t ) const;
			auto d( void ) const;
			template< unsigned int ... _OutDims , unsigned int ... _InDims >
			friend std::ostream &operator << ( std::ostream &os , const Constant< ParameterPack::UIntPack< _OutDims ... > , ParameterPack::UIntPack< _InDims ... > > &c );
		protected:
			const Tensor< ParameterPack::UIntPack< OutDims ... > > _c;
		};

		// A class for describing a linear function
		template< unsigned int ... OutDims , unsigned int ... InDims >
		struct Linear< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > > : public Function< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > , Linear< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > > >
		{
			typedef Function< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > , Linear > _Function;

			// A constructor generating the zero linear function
			Linear( void ){}

			// A constructor generating a linear function with the prescribed tensor taking input tensors to output tensors
			Linear( const Tensor< ParameterPack::Concatenation< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > > > &l ) : _l(l){}

			// A constructor generating a simple Linear function with one in the entry {out,in} and zero for all other entries.
			Linear( std::initializer_list< unsigned int > out , std::initializer_list< unsigned int > in );

			// A constructor generating a simple Linear function with one in the entry {_OutDims,_InDims} and zero for all other entries.
			template< unsigned int ... _OutDims , unsigned int ... _InDims > Linear( ParameterPack::UIntPack< _OutDims ... > , ParameterPack::UIntPack< _InDims ... > );

			auto value( const Tensor< ParameterPack::UIntPack< InDims ... > > &t ) const;
			auto d( void ) const;
			template< unsigned int ... _OutDims , unsigned int ... _InDims >
			friend std::ostream &operator << ( std::ostream &os , const Linear< ParameterPack::UIntPack< _OutDims ... > , ParameterPack::UIntPack< _InDims ... > > &l );

		protected:
			Tensor< ParameterPack::Concatenation< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > > > _l;
		};

		template< unsigned int ... Dims >
		struct Identity< ParameterPack::UIntPack< Dims ... > > : public Function< ParameterPack::UIntPack< Dims ... > , ParameterPack::UIntPack< Dims ... > , Identity< ParameterPack::UIntPack< Dims ... > > >
		{
			typedef Function< ParameterPack::UIntPack< Dims ... > , ParameterPack::UIntPack< Dims ... > , Identity > _Function;

			Identity( void ) : _identity( Tensor< ParameterPack::UIntPack< Dims ... > >::Identity() ) {}
			auto value( const Tensor< ParameterPack::UIntPack< Dims ... > > &t ) const;
			auto d( void ) const;
			template< unsigned int ... _Dims >
			friend std::ostream &operator << ( std::ostream &os , const Identity< ParameterPack::UIntPack< _Dims ... > > &id );
		protected:
			Tensor< ParameterPack::UIntPack< Dims ... , Dims ... > > _identity;
		};

		// A class for describing the product of a function with a scalar
		template< typename F >
		struct _Scale : public Function< typename F::OutPack , typename F::InPack , _Scale< F > >
		{
			typedef Function< typename F::OutPack , typename F::InPack , _Scale > _Function;
			template< typename _F > friend struct _Scale;

			_Scale( const F &f , double s ) : _f(f) , _s(s) {}

			auto value( const Tensor< typename _Function::InPack > &t ) const;
			auto d( void ) const;
			template< typename _F > friend std::ostream &operator << ( std::ostream &os , const _Scale< _F > &scale );
			const F &f( void ) const { return f; }
		protected:
			template< typename _F > friend auto operator * ( const _Scale< _F > & , double );
			template< typename _F > friend auto operator * ( double , const _Scale< _F > & );
			const F _f;
			const double _s;
		};

		// A class for describing the sum of two or more functions (with the same order input and the same order output)
		template< typename ... Fs > struct _Add;

		template< typename F , typename ... Fs >
		struct _Add< F , Fs ... > : public Function< typename F::OutPack , typename F::InPack , _Add< F , Fs ... > >
		{
			static_assert( AND< ParameterPack::Comparison< typename F::OutPack , typename Fs::OutPack >::Equal ... >() , "[ERROR] Output types differ" );
			static_assert( AND< ParameterPack::Comparison< typename F:: InPack , typename Fs:: InPack >::Equal ... >() , "[ERROR] Input types differ" );

			typedef Function< typename F::OutPack , typename F::InPack , _Add > _Function;

			_Add( const std::tuple< F , Fs ... > f ) : _f(f) {}

			auto value( const Tensor< typename _Function::InPack > &t ) const;
			auto d( void ) const;
			template< typename _F , typename ... _Fs >
			friend std::ostream &operator << ( std::ostream &os , const _Add< _F , _Fs... > &_Add );
			const std::tuple< F , Fs ... > &f_tuple( void ) const { return _f; }
		protected:
			template< unsigned int I > void _toStream( std::ostream &os ) const;
			template< unsigned int I > auto _d( void ) const;
			template< unsigned int I > auto _value( const Tensor< typename _Function::InPack > &t ) const;

			const std::tuple< F , Fs... > _f;
		};

		// A class that permutes the dimensions of the output tensor
		template< typename PermutationPack , typename F > struct _Permutation;

		template< unsigned int ... PermutationIndices , typename F >
		struct _Permutation< ParameterPack::UIntPack< PermutationIndices ... > , F > : public Function< ParameterPack::Permutation< typename F::OutPack , ParameterPack::UIntPack< PermutationIndices ... > > , typename F::InPack , _Permutation< ParameterPack::UIntPack< PermutationIndices ... > , F > >
		{
			typedef ParameterPack::UIntPack< PermutationIndices ... > PermutationPack;
			static_assert( PermutationPack::Size==F::OutPack::Size , "[ERROR] Sizes don't match" );
			typedef Function< ParameterPack::Permutation< typename F::OutPack , PermutationPack > , typename F::InPack , _Permutation< PermutationPack , F > > _Function;

			_Permutation( const F &f ) : _f(f) {}
			auto value( const Tensor< typename _Function::InPack > &t ) const;
			auto d( void ) const;
			template< unsigned int ... _PermutationIndices , typename _F > friend std::ostream &operator << ( std::ostream & , const _Permutation< ParameterPack::UIntPack< _PermutationIndices ... > , _F > & );

		protected:
			const F _f;
		};

		// A class for describing the product of two functions (with the same order input)
		template< unsigned int I , typename F1 , typename F2 >
		struct _ContractedOuterProduct : public Function< ParameterPack::Concatenation< typename ParameterPack::Partition< F1::OutPack::Size-I , typename F1::OutPack >::First , typename ParameterPack::Partition< I , typename F2::OutPack >::Second > , typename F1::InPack , _ContractedOuterProduct< I , F1 , F2 > >
		{
			static_assert( ParameterPack::Comparison< typename F2::InPack , typename F2::InPack >::Equal , "[ERROR] Input types differ" );
			typedef typename F1::OutPack OutPack1;
			typedef typename F2::OutPack OutPack2;

			typedef Function< ParameterPack::Concatenation< typename ParameterPack::Partition< F1::OutPack::Size-I , typename F1::OutPack >::First , typename ParameterPack::Partition< I , typename F2::OutPack >::Second > , typename F1::InPack , _ContractedOuterProduct > _Function;

			_ContractedOuterProduct( const F1 &f1 , const F2 &f2 ) : _f1(f1) , _f2(f2) {}

			auto value( const Tensor< typename _Function::InPack > &t ) const;
			auto d( void ) const;
			template< unsigned int _I , typename _F1 , typename _F2 > friend std::ostream &operator << ( std::ostream &os , const _ContractedOuterProduct< _I , _F1 , _F2 > &op );

		protected:
			const F1 _f1;
			const F2 _f2;
		};

		// A class that contracts along two dimensions of the output
		template< unsigned int I1 , unsigned int I2 , typename F >
		struct _Contraction : public Function< typename ParameterPack::Selection< I1 , typename ParameterPack::Selection< I2 , typename F::OutPack >::Complement >::Complement , typename F::InPack , _Contraction< I1 , I2 , F > >
		{
			typedef Function< typename ParameterPack::Selection< I1 , typename ParameterPack::Selection< I2 , typename F::OutPack >::Complement >::Complement , typename F::InPack , _Contraction > _Function;

			_Contraction( const F &f ) : _f(f) {}

			auto value( const Tensor< typename _Function::InPack > &t ) const;
			auto d( void ) const;
			template< unsigned int _I1 , unsigned int _I2 , typename _F > friend std::ostream &operator << ( std::ostream &os , const _Contraction< _I1 , _I2 , _F > &op );

		protected:
			const F _f;
		};

		// A class for describing the composition of two functions
		template< typename F1 , typename F2 >
		struct _Composition : public Function< typename F1::OutPack , typename F2::InPack , _Composition< F1 , F2 > >
		{
			static_assert( ParameterPack::Comparison< typename F1::InPack , typename F2::OutPack >::Equal , "[ERROR] Input/Output types differ" );

			typedef Function< typename F1::OutPack , typename F2::InPack , _Composition > _Function;

			_Composition( const F1 &f1 , const F2 &f2 ) : _f1(f1) , _f2(f2) {}

			auto value( const Tensor< typename _Function::InPack > &t ) const;
			auto d( void ) const;
			template< typename _F1 , typename _F2 >
			friend std::ostream &operator << ( std::ostream &os , const _Composition< _F1 , _F2 > &composition );
		protected:
			const F1 _f1;
			const F2 _f2;
		};

#ifdef NEW_AUTO_DIFF_CODE
		// A function returning the concatenation of the individual functions
		template< typename F , typename ... Fs >
		struct _Concatenation : public Function< ParameterPack::Concatenation< ParameterPack::UIntPack< sizeof...(Fs)+1 > , typename F::OutPack > , typename F::InPack , _Concatenation< F , Fs ... > >
		{
			static_assert( AND< ParameterPack::Comparison< typename F::OutPack , typename Fs::OutPack >::Equal ... >() , "[ERROR] Output types differ" );
			static_assert( AND< ParameterPack::Comparison< typename F:: InPack , typename Fs:: InPack >::Equal ... >() , "[ERROR] Input types differ" );

			typedef Function< ParameterPack::Concatenation< ParameterPack::UIntPack< sizeof...(Fs)+1 > , typename F::OutPack > , typename F::InPack , _Concatenation > _Function;

			_Concatenation( const std::tuple< F , Fs ... > f ) : _f(f) {}

			auto value( const Tensor< typename _Function::InPack > &t ) const;
			auto d( void ) const;
			template< typename _F , typename ... _Fs >
			friend std::ostream &operator << ( std::ostream &os , const _Concatenation< _F , _Fs ... > &concatenation );
		protected:
			template< unsigned int I > void _toStream( std::ostream &os ) const;
			template< unsigned int I , typename ... DFs > auto _d( DFs ... dFs ) const;
			template< unsigned int I > void _setValue( AutoDiff::Tensor< typename _Function::OutPack > &v , const Tensor< typename _Function::InPack > &t ) const;
			const std::tuple< F , Fs... > _f;
		};

		// A function returning the first function if the condition is met and the second otherwise
		template< typename ConditionFunctor , typename F1 , typename F2 >
		struct _Conditional : public Function< typename F1::OutPack , typename F1::InPack , _Conditional< ConditionFunctor , F1 , F2 > >
		{
			static_assert( ParameterPack::Comparison< typename F1::InPack , typename F1::InPack >::Equal , "[ERROR] Input types differ" );
			static_assert( ParameterPack::Comparison< typename F1::OutPack , typename F1::OutPack >::Equal , "[ERROR] Output types differ" );

			typedef Function< typename F1::OutPack , typename F1::InPack , _Conditional > _Function;

			_Conditional( ConditionFunctor c , const F1 &f1 , const F2 &f2 ) : _c(c) , _f1(f1) , _f2(f2) {}

			auto value( const Tensor< typename _Function::InPack > &t ) const;
			auto d( void ) const;
			template< typename _ConditionFunctor , typename _F1 , typename _F2 >
			friend std::ostream &operator << ( std::ostream &os , const _Conditional< _ConditionFunctor , _F1 , _F2 > &conditional );
		protected:
			const ConditionFunctor _c;
			const F1 _f1;
			const F2 _f2;
		};
#endif // NEW_AUTO_DIFF_CODE

		// A class for extracting a sub-tensor from the output
		template< unsigned int I , typename F >
		struct _Extract : public Function< typename ParameterPack::Partition< I , typename F::OutPack >::Second , typename F::InPack , _Extract< I , F > >
		{
			typedef Function< typename ParameterPack::Partition< I , typename F::OutPack >::Second , typename F::InPack , _Extract > _Function;

			_Extract( const unsigned int indices[/*I*/] , const F &f ) : _f(f) { memcpy( _indices , indices , sizeof(unsigned int)*I ); }

			auto value( const Tensor< typename _Function::InPack > &t ) const;
			auto d( void ) const;
			template< unsigned int _I , typename _F > friend std::ostream &operator << ( std::ostream &os , const _Extract< _I , _F > &ex );

		protected:
			const F _f;
			unsigned int _indices[I];
		};

		////////////////////////
		////////////////////////
		//// Implementation ////
		////////////////////////
		////////////////////////

		//////////////
		// Function //
		//////////////
		template< unsigned int ... OutDims , unsigned int ... InDims , typename F >
		template< typename V >
		auto Function< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > , F >::operator()( const V &v ) const
		{
			if constexpr( std::is_base_of< Tensor< InPack > , V >::value ) return static_cast< const F & >( *this ).value( v );
			else if constexpr( std::is_arithmetic_v< V > && InPack::Size==0 ) return static_cast< const F & >( *this ).value( Tensor< ParameterPack::UIntPack<> >( v ) );
			else return Composition< F , V >( static_cast< const F & >( *this ) , v );
		}

		//////////////
		// Constant //
		//////////////
		// A class for reprenting a constant function
		template< unsigned int ... OutDims , unsigned int ... InDims >
		auto Constant< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > >::value( const Tensor< ParameterPack::UIntPack< InDims ... > > &t ) const { return _c; }

		template< unsigned int ... OutDims , unsigned int ... InDims >
		auto Constant< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > >::d( void ) const { return Constant< ParameterPack::Concatenation< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > > , ParameterPack::UIntPack< InDims ... > >(); }

		template< unsigned int ... OutDims , unsigned int ... InDims >
		std::ostream &operator << ( std::ostream &os , const Constant< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > > &c ){ return os << c._c; }

		////////////
		// Linear //
		////////////
		// A constructor generating a simple Linear function with one in the entry {out,in} and zero for all other entries.
		template< unsigned int ... OutDims , unsigned int ... InDims >
		Linear< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > >::Linear( std::initializer_list< unsigned int > out , std::initializer_list< unsigned int > in )
		{
			if( out.size()!=sizeof...(OutDims) || in.size()!=sizeof...(InDims) ) ERROR_OUT( "Output dimensions don't match" );
			if constexpr( sizeof...(OutDims)+sizeof...(InDims)==0 ) _l = 1.;
			else
			{
				unsigned int idx[ sizeof...(OutDims)+sizeof...(InDims) ];
				unsigned int c = 0;
				for( auto it=out.begin() ; it!=out.end() ; it++ ) idx[c++] = *it;
				for( auto it= in.begin() ; it!= in.end() ; it++ ) idx[c++] = *it;
				_l(idx) = 1;
			}
		}

		// A constructor generating a simple Linear function with one in the entry {_OutDims,_InDims} and zero for all other entries.
		template< unsigned int ... OutDims , unsigned int ... InDims >
		template< unsigned int ... _OutDims , unsigned int ... _InDims >
		Linear< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > >::Linear( ParameterPack::UIntPack< _OutDims ... > , ParameterPack::UIntPack< _InDims ... > )
		{
			static_assert( sizeof...(_OutDims)==sizeof...(OutDims) && sizeof...(_InDims)==sizeof...(InDims) , "[ERROR] Size mismatch" );
			unsigned int idx[] = { _OutDims ... , _InDims ... };
			_l(idx) = 1;
		}

		template< unsigned int ... OutDims , unsigned int ... InDims >
		auto Linear< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > >::value( const Tensor< ParameterPack::UIntPack< InDims ... > > &t ) const
		{
			return _l.template contractedOuterProduct< sizeof ... ( InDims ) >( t );
		}

		template< unsigned int ... OutDims , unsigned int ... InDims >
		auto Linear< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > >::d( void ) const { return Constant< ParameterPack::Concatenation< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > > , ParameterPack::UIntPack< InDims ... > >(_l); }

		template< unsigned int ... OutDims , unsigned int ... InDims >
		std::ostream &operator << ( std::ostream &os , const Linear< ParameterPack::UIntPack< OutDims ... > , ParameterPack::UIntPack< InDims ... > > &l ){ return os << l._l; }

		//////////////
		// Identity //
		//////////////
		template< unsigned int ... Dims >
		auto Identity< ParameterPack::UIntPack< Dims ... > >::value( const Tensor< ParameterPack::UIntPack< Dims ... > > &t ) const { return t; }

		template< unsigned int ... Dims >
#ifdef NEW_AUTO_DIFF_CODE
		auto Identity< ParameterPack::UIntPack< Dims ... > >::d( void ) const { return Constant< ParameterPack::UIntPack< Dims ... > , ParameterPack::UIntPack< Dims ... > >( _identity ); }
#else // !NEW_AUTO_DIFF_CODE
		auto Identity< ParameterPack::UIntPack< Dims ... > >::d( void ) const { return Constant( _identity ); }
#endif // NEW_AUTO_DIFF_CODE

		template< unsigned int ... Dims >
		std::ostream &operator << ( std::ostream &os , const Identity< ParameterPack::UIntPack< Dims ... > > &id )
		{
			return os << "Id_{" << ParameterPack::UIntPack< Dims ... >() << "}";
		}

		////////////
		// _Scale //
		////////////
		template< typename F >
		auto _Scale< F >::value( const Tensor< typename _Function::InPack > &t ) const { return _f(t) * _s; }

		template< typename F >
		auto _Scale< F >::d( void ) const { return _Scale< decltype( std::declval< F >().d() ) >( _f.d() , _s ); }

		template< typename F >
		std::ostream &operator << ( std::ostream &os , const _Scale< F > &scale )
		{
			if( scale._s==-1 ) return os << "( -" << scale._f << " )";
			else               return os << "( " << scale._s << " * " << scale._f << " )";
		}

		//////////
		// _Add //
		//////////
		template< typename F , typename ... Fs >
		auto _Add< F , Fs ... >::value( const Tensor< typename _Function::InPack > &t ) const { return this->template _value<0>(t); }

		template< typename F , typename ... Fs >
		auto _Add< F , Fs ... >::d( void ) const { return this->template _d<0>(); }

		template< typename F , typename ... Fs >
		std::ostream &operator << ( std::ostream &os , const _Add< F , Fs ... > &add )
		{
			os << "( ";
			add.template _toStream<0>( os );
			os << " )";
			return os;
		}

		template< typename F , typename ... Fs >
		template< unsigned int I >
		void _Add< F , Fs ... >::_toStream( std::ostream &os ) const
		{
			if constexpr( I==0 ) os << std::get<I>( _f );
			else os << " + " << std::get<I>( _f );
			if constexpr( I<sizeof...(Fs) ) this->template _toStream< I+1 >( os );
		}

		template< typename F , typename ... Fs >
		template< unsigned int I >
		auto _Add< F , Fs ... >::_d( void ) const
		{
			if constexpr( I==sizeof...(Fs) ) return std::get<I>(_f).d();
			else return std::get<I>(_f).d() + this->template _d<I+1>();
		}

		template< typename F , typename ... Fs >
		template< unsigned int I >
		auto _Add< F , Fs ... >::_value( const Tensor< typename _Function::InPack > &t ) const
		{
			if constexpr( I==sizeof...(Fs) ) return std::get<I>(_f)(t);
			else return std::get<I>(_f)(t) + this->template _value<I+1>(t);
		}

		//////////////
		// Negation //
		//////////////
		template< typename F >
		auto operator - ( const F &f ){ return f * -1.; }

		//////////////
		// Addition //
		//////////////
		template< typename F1 , typename F2 >
		auto operator + ( const F1 &f1 , const F2 &f2 ){ return _Add< F1 , F2 >( std::make_tuple(f1,f2) ); }

		template< typename F1 , typename F2 , typename ... Fs >
		auto operator + ( const F1 &f , const _Add< F2 , Fs ... > &add ){ return _Add< F1 , F2 , Fs ... >( std::tuple_cat( std::make_tuple(f) , add.f_tuple() ) ); }

		template< typename F1 , typename ... Fs , typename F2 >
		auto operator + ( const _Add< F1 , Fs ... > &add , const F2 &f ){ return _Add< F1 , Fs ... , F2 >( std::tuple_cat( add.f_tuple() , std::make_tuple(f) ) ); }

		template< typename F1 , typename ... F1s , typename F2 , typename ... F2s >
		auto operator + ( const _Add< F1 , F1s ... > &add1 , const _Add< F2 , F2s ... > &add2 ){ return _Add< F1 , F1s ... , F2 , F2s ... >( std::tuple_cat( add1.f_tuple() , add2.f_tuple() ) ); }

		template< typename F >
		auto operator + ( const F &f , const Tensor< typename F::OutPack > &t ){ return f + Constant< typename F::OutPack , typename F::InPack >( t ); }

		template< typename F >
		auto operator + ( const Tensor< typename F::OutPack > &t , const F &f ){ return Constant< typename F::OutPack , typename F::InPack >( t ) + f; }

		template< typename F >
		auto operator + ( const F &f , double s )
		{
			static_assert( ParameterPack::Comparison< typename F::OutPack , ParameterPack::UIntPack<> >::Equal , "[ERROR] function must be scalar-valued" );
			return f + Constant< typename F::OutPack , typename F::InPack >( s );
		}

		template< typename F >
		auto operator + ( double s , const F &f )
		{
			static_assert( ParameterPack::Comparison< typename F::OutPack , ParameterPack::UIntPack<> >::Equal , "[ERROR] function must be scalar-valued" );
			return Constant< typename F::OutPack , typename F::InPack >( s ) + f;
		}

		/////////////////
		// Subtraction //
		/////////////////
		template< typename F1 , typename F2 >
		auto operator - ( const F1 &f1 , const F2 &f2 ){ return f1 + ( -f2 ); }

		template< typename F >
		auto operator - ( const F &f , const Tensor< typename F::OutPack > &t ){ return f + Constant< typename F::OutPack , typename F::InPack >( -t ); }

		template< typename F >
		auto operator - ( const Tensor< typename F::OutPack > &t , const F &f ){ Constant< typename F::OutPack , typename F::InPack >( t ) - f; }

		template< typename F >
		auto operator - ( const F &f , double s )
		{
			static_assert( ParameterPack::Comparison< typename F::OutPack , ParameterPack::UIntPack<> >::Equal , "[ERROR] function must be scalar-valued" );
			return f + Constant< typename F::OutPack , typename F::InPack >( -s );
		}

		template< typename F >
		auto operator - ( double s , const F &f )
		{
			static_assert( ParameterPack::Comparison< typename F::OutPack , ParameterPack::UIntPack<> >::Equal , "[ERROR] function must be scalar-valued" );
			return Constant< typename F::OutPack , typename F::InPack >( s ) - f;
		}

		////////////////////
		// Scalar product //
		////////////////////
		template< typename F >
		auto operator * ( const F &f , const double &s ){ return _Scale< F >(f,s); }

		template< typename F >
		auto operator * ( const _Scale< F > &f , double s ){ return _Scale< F >( f._f , f._s*s ); }

		template< typename F >
		auto operator * ( const double &s , const F &f ){ return f*s; }

		/////////////////////
		// Scalar quotient //
		/////////////////////
		template< typename F >
		auto operator / ( const F &f , const double &s ){ return f * (1./s); }

		////////////////////
		// Tensor product //
		////////////////////
		template< typename F1 , typename F2 >
		auto operator * ( const F1 &f1 , const F2 &f2 ){ return _ContractedOuterProduct< 0 , F1 , F2 >(f1,f2); }

		/////////////////////
		// Tensor quotient //
		/////////////////////
		template< typename F1 , typename F2 >
		auto operator / ( const F1 &f1 , const F2 &f2 )
		{
			static_assert( ParameterPack::Comparison< typename F2::OutPack , ParameterPack::UIntPack<> >::Equal , "[ERROR] denominator is not scalar-valued" );

			if constexpr( std::is_same< F1 , double >::value ) return Constant< ParameterPack::UIntPack<> , typename F2::InPack >( f1 );
			else if constexpr( std::is_same< F1 , Tensor< ParameterPack::UIntPack<> > >::value ) return Constant< ParameterPack::UIntPack<> , typename F2::InPack >( f1 );
			else
			{
				static_assert( ParameterPack::Comparison< typename F1::InPack , typename F2::InPack >::Equal , "[ERROR] Input types differ" );
				return f1 * Pow( f2 , -1. );
			}
		}

		/////////////////
		// Permutation //
		/////////////////
		template< unsigned int ... PermutationIndices , typename F >
		auto _Permutation< ParameterPack::UIntPack< PermutationIndices ... > , F >::value( const Tensor< typename _Function::InPack > &t ) const
		{
			return _f( t ).permute( PermutationPack() );
		}

		template< unsigned int ... PermutationIndices , typename F >
		auto _Permutation< ParameterPack::UIntPack< PermutationIndices ... > , F >::d( void ) const
		{
			return Permutation< ParameterPack::Concatenation< PermutationPack , ParameterPack::SequentialPack< unsigned int , _Function::InPack::Size , PermutationPack::Size > > >( _f.d() ); 
		}

		template< unsigned int ... PermutationIndices , typename F >
		std::ostream &operator << ( std::ostream &os , const _Permutation< ParameterPack::UIntPack< PermutationIndices ... > , F > &p )
		{
			return os << "Permutation" << _Permutation< ParameterPack::UIntPack< PermutationIndices ... > , F >::PermutationPack() << "( " << p._f << " )";
		}

		template< typename PermutationPack , typename F >
		auto Permutation( const F &f ){ return _Permutation< PermutationPack , F >( f ); }

		////////////////////////////
		// ContractedOuterProduct //
		////////////////////////////
		template< unsigned int I , typename F1 , typename F2 >
		auto _ContractedOuterProduct< I , F1 , F2 >::value( const Tensor< typename _Function::InPack > &t ) const
		{
			return _f1(t).template contractedOuterProduct< I >( _f2(t) );
		}

		template< unsigned int I , typename F1 , typename F2 >
		auto _ContractedOuterProduct< I , F1 , F2 >::d( void ) const
		{
			// [Note] As a function is differentiated, the input type is appended to the output type

			// Get the PermutationPack permuting _f1.d() so that the InPack is on the left
			typedef ParameterPack::SequentialPack< unsigned int , F1::InPack::Size , F1::OutPack::Size > PreP1Pack;
			typedef ParameterPack::SequentialPack< unsigned int , F1::OutPack::Size , 0 > PreP2Pack;
			typedef ParameterPack::Concatenation< PreP1Pack , PreP2Pack > PrePermutationPack;

			// Get the PermutationPack permuting _f1.d() * f2 so that the InPack is on the right
#ifdef NEW_AUTO_DIFF_CODE
			using namespace ParameterPack;
			typedef ParameterPack::Concatenation< typename Partition< F1::OutPack::Size-I , typename F1::OutPack >::First , typename ParameterPack::Partition< I , typename F2::OutPack >::Second > ContractionPack;
#else // !NEW_AUTO_DIFF_CODE
			typedef ParameterPack::Concatenation< typename ParameterPack::Partition< F1::OutPack::Size-I , typename F1::OutPack >::First , typename ParameterPack::Partition< I , typename F2::OutPack >::Second > ContractionPack;
#endif // NEW_AUTO_DIFF_CODE
			typedef ParameterPack::SequentialPack< unsigned int , ContractionPack::Size , F1::InPack::Size > PostP1Pack;
			typedef ParameterPack::SequentialPack< unsigned int , F1::InPack::Size , 0 > PostP2Pack;
			typedef ParameterPack::Concatenation< PostP1Pack , PostP2Pack > PostPermutationPack;

			return Permutation< PostPermutationPack >( ContractedOuterProduct< I >( Permutation< PrePermutationPack >( _f1.d() ) , _f2 ) ) + ContractedOuterProduct< I >( _f1 , _f2.d() );
		}

		template< unsigned int I , typename F1 , typename F2 >
		std::ostream &operator << ( std::ostream &os , const _ContractedOuterProduct< I , F1 , F2 > &op )
		{
			if( I==0 ) return os << "( " << op._f1 << " * " << op._f2 << " )";
			else       return os << "( " << op._f1 << " *_" << I << " " << op._f2 << " )";
		}

		template< unsigned int I , typename F1 , typename F2 >
		auto ContractedOuterProduct( const F1 &f1 , const F2 &f2 ){ return _ContractedOuterProduct< I , F1 , F2 >( f1 , f2 ); }

		//////////////////
		// _Contraction //
		//////////////////
		template< unsigned int I1 , unsigned int I2 , typename F >
		auto _Contraction< I1 , I2 , F >::value( const Tensor< typename _Function::InPack > &t ) const { return _f(t).template contract< I1 , I2 >(); }

		template< unsigned int I1 , unsigned int I2 , typename F >
		auto _Contraction< I1 , I2 , F >::d( void ) const { return Contraction< I1 , I2 >( _f.d() ); }

		template< unsigned int I1 , unsigned int I2 , typename F >
		std::ostream &operator << ( std::ostream &os , const _Contraction< I1 , I2 , F > &contraction )
		{
			return os << "Contraction_{" << I1 << " , " << I2 << "}( "<< contraction._f << ")";
		}

		/////////////////
		// Contraction //
		/////////////////
		template< unsigned I1 , unsigned I2 , typename F >
		auto Contraction( const F &f )
		{
			if constexpr( I1<I2 ) return _Contraction< I1 , I2 , F >( f );
			else                  return _Contraction< I2 , I1 , F >( f );
		}

		/////////////////
		// Composition //
		/////////////////
		template< typename F1 , typename F2 >
		auto _Composition< F1 , F2 >::value( const Tensor< typename _Function::InPack > &t ) const { return _f1( _f2(t) ); }

		template< typename F1 , typename F2 >
		auto _Composition< F1 , F2 >::d( void ) const { return ContractedOuterProduct< F1::InPack::Size >( Composition( _f1.d() , _f2 ) , _f2.d() ); }

		template< typename F1 , typename F2 >
		std::ostream &operator << ( std::ostream &os , const _Composition< F1 , F2 > &composition ){ return os << composition._f1 << "( " << composition._f2 << " )"; }

		template< typename F1 , typename F2 > auto Composition( const F1 &f1 , const F2 &f2 ){ return _Composition< F1 , F2 >( f1 , f2 ); }

#ifdef NEW_AUTO_DIFF_CODE
		///////////////////
		// Concatenation //
		///////////////////
		template< typename F , typename ... Fs >
		auto _Concatenation< F , Fs ... >::value( const Tensor< typename _Function::InPack > &t ) const
		{
			AutoDiff::Tensor< typename _Function::OutPack > out;
			_setValue< 0 >( out , t );
			return out;
		}
		template< typename F , typename ... Fs >
		template< unsigned int I >
		void _Concatenation< F , Fs ... >::_setValue( AutoDiff::Tensor< typename _Function::OutPack > &v , const Tensor< typename _Function::InPack > &t ) const
		{
			if constexpr( I==(sizeof...(Fs)+1) ) return;
			else
			{
				static const unsigned int Size = AutoDiff::Tensor< typename F::OutPack >::Size;
				AutoDiff::Tensor< typename F::OutPack > out = std::get< I >( _f )(t);
				for( unsigned int i=0 ; i<Size ; i++ ) v[I*Size+i] = out[i];
				_setValue<I+1>( v , t );
			}
		}

		template< typename F , typename ... Fs >
		auto _Concatenation< F , Fs ... >::d( void ) const { return _d<0>(); }

		template< typename F , typename ... Fs >
		template< unsigned int I , typename ... DFs >
		auto _Concatenation< F , Fs ... >::_d( DFs ... dFs ) const
		{
			if constexpr( I==sizeof...(Fs)+1 ) return _Concatenation< DFs ... >( std::make_tuple( dFs ... ) );
			else return _d<I+1>( dFs ... , std::get<I>( _f ).d() );
		}

		template< typename F , typename ... Fs >
		std::ostream &operator << ( std::ostream &os , const _Concatenation< F , Fs ... > &concatenation )
		{
			os << " { ";
			concatenation.template _toStream<0>( os );
			os << " } ";
			return os;
		}

		template< typename F , typename ... Fs >
		template< unsigned int I > void _Concatenation< F , Fs ... >::_toStream( std::ostream &os ) const
		{
			if constexpr( I==(sizeof...(Fs)+1) ) return;
			if( I ) os << " , ";
			os << _f.template get<I>();
			_toStream<I+1>( os );
		}

		template< typename F , typename ... Fs > auto Concatenation( const F &f , const Fs &...fs ){ return _Concatenation< F , Fs... >( std::make_tuple( f , fs ... ) ); }

		/////////////////
		// Conditional //
		/////////////////
		template< typename CFunctor , typename F1 , typename F2 >
		auto _Conditional< CFunctor , F1 , F2 >::value( const Tensor< typename _Function::InPack > &t ) const { return _c(t) ? _f1(t) : _f2(t); }

		template< typename CFunctor , typename F1 , typename F2 >
		auto _Conditional< CFunctor , F1 , F2 >::d( void ) const { return Conditional( _c , _f1.d() , _f2.d() ); }

		template< typename CFunctor , typename F1 , typename F2 >
		std::ostream &operator << ( std::ostream &os , const _Conditional< CFunctor , F1 , F2 > &conditional ){ return os << conditional._f1 << " or " << conditional._f2; }

		template< typename CFunctor , typename F1 , typename F2 > auto Conditional( CFunctor c , const F1 &f1 , const F2 &f2 ){ return _Conditional< CFunctor , F1 , F2 >( c , f1 , f2 ); }
#endif // NEW_AUTO_DIFF_CODE

		//////////////
		// _Extract //
		//////////////
		template< unsigned int I , typename F >
		auto _Extract< I , F >::value( const Tensor< typename _Function::InPack > &t ) const { return _f(t).template extract<I>( _indices ); }

		template< unsigned int I , typename F >
		auto _Extract< I , F >::d( void ) const { return Extract< I >( _indices , _f.d() ); }

		template< unsigned int I , typename F >
		std::ostream &operator << ( std::ostream &os , const _Extract< I , F > &extract )
		{
			if constexpr( I==F::OutPack::Size ) os << "Coefficient_{";
			else                                os << "Extract_{";
			for( unsigned int i=0 ; i<I ; i++ )
				if( i ) os << "," << extract._indices[i];
				else    os <<        extract._indices[i];
			return os << "}( " << extract._f << " )";
		}

		/////////////
		// Extract //
		/////////////
		template< unsigned int I , typename F >
		auto Extract( const unsigned int indices[] , const F &f ){ return _Extract< I , F >( indices , f ); }

		/////////////////
		// Coefficient //
		/////////////////
		template< typename F >
		auto Coefficient( const unsigned int indices[/*F::OutPack::Size*/] , const F &f ){ return Extract< F::OutPack::Size >( indices , f ); }

		///////////////
		// Transpose //
		///////////////
		template< typename F >
		auto Transpose( const F &f ){ return Permutation< typename ParameterPack::SequentialPack< unsigned int , F::OutPack::Size >::Transpose >( f ); }

		////////////////
		// SquareNorm //
		////////////////
		template< typename F >
		auto SquareNorm( const F &f ){ return ContractedOuterProduct< F::OutPack::Size >( f , f ); }

		/////////////////
		// Determinant //
		/////////////////
		template< unsigned int R , unsigned int C ,typename F >
		auto _SubMatrix( const F &f )
		{
			typedef typename F::OutPack OutPack;
			static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
			static_assert( OutPack::template Get<0>()==OutPack::template Get<1>() , "[ERROR] Output 2-tensor must be square" );
			static const unsigned int Dim = OutPack::template Get<0>();

			// Create the (Dim-1)x(Dim-1) matrix
			Tensor< ParameterPack::UIntPack< Dim-1 , Dim-1 , Dim , Dim > > l;
			unsigned int index[4];
			for( unsigned int i=0 ; i<Dim-1 ; i++ ) for( unsigned int j=0 ; j<Dim-1 ; j++ )
			{
				index[0] = i;
				index[1] = j;
				index[2] = i<R ? i : i+1;
				index[3] = j<C ? j : j+1;
				l( index ) = 1;
			}
			Linear< ParameterPack::UIntPack< Dim-1 , Dim-1 > , ParameterPack::UIntPack< Dim , Dim > > L( l );
			return L(f);
		}

		template< unsigned int C ,typename F >
		auto _Determinant0( const F &f )
		{
			typedef typename F::OutPack OutPack;
			static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
			static_assert( OutPack::template Get<0>()==OutPack::template Get<1>() , "[ERROR] Output 2-tensor must be square" );
			static const unsigned int Dim = OutPack::template Get<0>();

			auto _det = Determinant( _SubMatrix<0,C>(f) ) * Coefficient( ParameterPack::UIntPack<0,C>::Values , f );

			if constexpr( C==0 ) return _det;
			else
			{
				if constexpr( C&1 ) return _Determinant0< C-1 >( f ) - _det;
				else                return _Determinant0< C-1 >( f ) + _det;
			}
		}

		template< typename F >
		auto Determinant( const F &f )
		{
			typedef typename F::OutPack OutPack;
			static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
			static_assert( OutPack::template Get<0>()==OutPack::template Get<1>() , "[ERROR] Output 2-tensor must be square" );
			static const unsigned int Dim = OutPack::template Get<0>();
			if      constexpr( Dim==1 ) return Coefficient( ParameterPack::UIntPack<0,0>::Values , f );
			else if constexpr( Dim==2 ) return Coefficient( ParameterPack::UIntPack<0,0>::Values , f ) * Coefficient( ParameterPack::UIntPack<1,1>::Values , f ) - Coefficient( ParameterPack::UIntPack<1,0>::Values , f ) * Coefficient( ParameterPack::UIntPack<0,1>::Values , f );
			else return _Determinant0< Dim-1 >( f );
		}

		//////////////
		// Cofactor //
		//////////////

		template< unsigned int R , unsigned int C , typename F >
		auto _Cofactor( const F &f )
		{
			typedef typename F::InPack InPack;
			typedef typename F::OutPack OutPack;
			static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
			static_assert( OutPack::template Get<0>()==OutPack::template Get<1>() , "[ERROR] Output 2-tensor must be square" );
			static const unsigned int Dim = OutPack::template Get<0>();

			// Create the Dim x Dim matrix
			Tensor< ParameterPack::UIntPack< Dim , Dim > > e;
			{
				unsigned int index[2];
				index[0] = R;
				index[1] = C;
				e( index ) = 1;
			}

			auto _det = Constant< ParameterPack::UIntPack< Dim , Dim > , InPack >( e ) * Determinant( _SubMatrix< R , C >( f ) );

			if constexpr( R==0 && C==0 ) return _det;
			else if constexpr( C==0 )
			{
				if constexpr( (R+C)&1 ) return _Cofactor< R-1 , Dim-1 >(f) - _det;
				else                    return _Cofactor< R-1 , Dim-1 >(f) + _det;
			}
			else
			{
				if constexpr( (R+C)&1 ) return _Cofactor< R , C-1 >( f ) - _det;
				else                    return _Cofactor< R , C-1 >( f ) + _det;
			}
		}

		template< typename F >
		auto Cofactor( const F& f )
		{
			typedef typename F::OutPack OutPack;
			static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
			static_assert( OutPack::template Get<0>()==OutPack::template Get<1>() , "[ERROR] Output 2-tensor must be square" );
			static const unsigned int Dim = OutPack::template Get<0>();
			return _Cofactor< Dim-1 , Dim-1 >( f );
		}

		//////////////
		// Adjugate //
		//////////////
		template< typename F >
		auto Adjugate( const F &f ){ return Transpose( Cofactor( f ) ); }

		/////////////
		// Inverse //
		/////////////
		template< typename F >
		auto Inverse( const F &f ){ return Adjugate( f ) / Determinant( f ); }

		//////////////////
		// CrossProduct //
		//////////////////
		template< unsigned int R , typename F >
		auto _CrossProduct( const F &f )
		{
			typedef typename F::InPack InPack;
			typedef typename F::OutPack OutPack;
			static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
			static_assert( OutPack::template Get<0>()==OutPack::template Get<1>()+1 , "[ERROR] Output 2-tensor must have one more row than column" );
			static const unsigned int Dim = OutPack::template Get<0>();

			// Create the Dim-dimensional vector
			Tensor< ParameterPack::UIntPack< Dim > > c;
			{
				unsigned int index[1];
				index[0] = R;
				c( index ) = 1;
			}

			// Create the (Dim-1)x(Dim-1) matrix
			Tensor< ParameterPack::UIntPack< Dim-1 , Dim-1 , Dim , Dim-1 > > l;
			{
				unsigned int index[4];
				for( unsigned int i=0 ; i<Dim-1 ; i++ ) for( unsigned int j=0 ; j<Dim-1 ; j++ )
				{
					index[0] = i;
					index[1] = j;
					index[2] = i<R ? i : i+1;
					index[3] = j;
					l( index ) = 1;
				}
			}
			Constant< ParameterPack::UIntPack< Dim > , InPack > C( c );
			Linear< ParameterPack::UIntPack< Dim-1 , Dim-1 > , ParameterPack::UIntPack< Dim , Dim-1 > > L( l );

			auto _det = C * Determinant( L( f ) );

			if constexpr( R==0 ) return _det;
			else
			{
				if constexpr( R&1 ) return _CrossProduct< R-1 >( f ) - _det;
				else                return _CrossProduct< R-1 >( f ) + _det;
			}
		}

		template< typename F >
		auto CrossProduct( const F& f )
		{
			typedef typename F::OutPack OutPack;
			static_assert( OutPack::Size==2 , "[ERROR] Output must be a 2-tensor" );
			static_assert( OutPack::template Get<0>()==OutPack::template Get<1>()+1 , "[ERROR] Output 2-tensor must have one more row than column" );
			static const unsigned int Dim = OutPack::template Get<0>();
			return _CrossProduct< Dim-1 >( f );
		}

		////////////////////////
		// DiscreteDerivative //
		////////////////////////
		template< typename F >
		auto DiscreteDerivative( const F &f , const Tensor< typename F::InPack > &in , double eps )
		{
			typedef Tensor< typename F::InPack > InTensor;
			typedef Tensor< typename F::OutPack > OutTensor;

			if constexpr( F::InPack::Size==0 ) return ( F( in + InTensor( eps ) ) - F( in - InTensor( eps ) ) ) / ( 2. * eps );
			else
			{
				Tensor< ParameterPack::Concatenation< typename F::OutPack , typename F::InPack > > out;

				unsigned int index[ F::InPack::Size ];
				MultiDimensionalArray::Loop< F::InPack::Size >::Run
				(
					ParameterPack::IsotropicUIntPack< F::InPack::Size >::Values , F::InPack::Values ,
					[&]( int d , int i ){ index[d] = i; } ,
					[&]( void )
					{
						InTensor in1 = in , in2 = in;
						in1( index ) += eps;
						in2( index ) -= eps;
						OutTensor _out = ( f( in1 ) - f( in2 ) ) / (2.*eps);
						if constexpr( F::OutPack::Size==0 ) out( index ) = (double)_out;
						else
						{
							unsigned int _index[ F::OutPack::Size + F::InPack::Size ];
							for( unsigned int i=0 ; i<F::InPack::Size ; i++ ) _index[i+F::OutPack::Size] = index[i];
							MultiDimensionalArray::Loop< F::OutPack::Size >::Run
							(
								ParameterPack::IsotropicUIntPack< F::OutPack::Size >::Values , F::OutPack::Values ,
								[&]( int d , int i ){ _index[d] = i; } ,
								[&]( const double &v ){ out( _index ) = v; } ,
								_out
							);
						}
					}
				);
				return out;
			}
		}

#include "AutoDiff.Transcendental.inc"
	}
}
#endif // AUTO_DIFF_INCLUDED
