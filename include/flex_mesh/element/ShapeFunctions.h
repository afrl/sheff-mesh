#pragma once

#include "flex_mesh/common/Common.h"
#include "flex_mesh/common/GeometryTypes.h"

#include <vector>

namespace flex_mesh
{
	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	class ShapeFunctions
	{
	public:

		static void values(const ElementType type,
			const Coordinate& parametric_coordinate,
			ElementVector& shape_func_values);

		static void derivatives(const ElementType type,
			const Coordinate& parametric_coordinate,
			ElementVectorDeriv& shape_func_derivatives);

	private:
		// 1D
		static void shapeLine2(const Coordinate& p, ElementVector& S);
		static void shapeDerivLine2(const Coordinate& p, ElementVectorDeriv& SDeriv);
		static void shapeLine3(const Coordinate& p, ElementVector& S);
		static void shapeDerivLine3(const Coordinate& p, ElementVectorDeriv& SDeriv);
		static void shapeLine4(const Coordinate& xi, ElementVector& S);
		static void shapeDerivLine4(const Coordinate& p, ElementVectorDeriv& SDeriv);
		static void shapeLine5(const Coordinate& p, ElementVector& S);
		static void shapeDerivLine5(const Coordinate& p, ElementVectorDeriv& SDeriv);
		// 2D
		static void shapeTri3(const Coordinate& p, ElementVector& S);
		static void shapeDerivTri3(const Coordinate& p, ElementVectorDeriv& SDeriv);
		static void shapeTri6(const Coordinate& p, ElementVector& S);
		static void shapeDerivTri6(const Coordinate& p, ElementVectorDeriv& SDeriv);
		static void shapeQuad4(const Coordinate& p, ElementVector& S);
		static void shapeDerivQuad4(const Coordinate& p, ElementVectorDeriv& SDeriv);
		static void shapeQuad8(const Coordinate& p, ElementVector& S);
		static void shapeDerivQuad8(const Coordinate& p, ElementVectorDeriv& SDeriv);
		// 3D
		static void shapePyr5(const Coordinate& p, ElementVector& S);
		static void shapeDerivPyr5(const Coordinate& p, ElementVectorDeriv& SDeriv);
		static void shapePyr13(const Coordinate& p, ElementVector& S);
		static void shapeDerivPyr13(const Coordinate& p, ElementVectorDeriv& SDeriv);
		static void shapeHex8(const Coordinate& p, ElementVector& S);
		static void shapeDerivHex8(const Coordinate& p, ElementVectorDeriv& SDeriv);
		static void shapeHex20(const Coordinate& p, ElementVector& S);
		static void shapeDerivHex20(const Coordinate& p, ElementVectorDeriv& SDeriv);
		static void shapeTet4(const Coordinate& p, ElementVector& S);
		static void shapeDerivTet4(const Coordinate& p, ElementVectorDeriv& SDeriv);
		static void shapeTet10(const Coordinate& p, ElementVector& S);
		static void shapeDerivTet10(const Coordinate& p, ElementVectorDeriv& SDeriv);
		static void shapeWedge15(const Coordinate& p, ElementVector& S);
		static void shapeDerivWedge15(const Coordinate& p, ElementVectorDeriv& SDeriv);
	};

	/////////////////
	// Definitions //
	/////////////////


	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	void flex_mesh::ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::values(const ElementType type, const Coordinate & parametric_coordinate, ElementVector & shape_func_values)
	{
		const auto num_nodes = get_num_nodes_for_type(type);
		shape_func_values.resize(num_nodes);

		if (type == ElementType::L1D2N)
			shapeLine2(parametric_coordinate, shape_func_values);
		else if (type == ElementType::L1D3N)
			shapeLine3(parametric_coordinate, shape_func_values);
		/*
		else if (type == ElementType::L1D4N)
			shapeLine4(parametric_coordinate, shape_func_values);
		else if (type == ElementType::L1D5N)
			shapeLine5(parametric_coordinate, shape_func_values);
		*/
		else if (type == ElementType::T2D3N)
			shapeTri3(parametric_coordinate, shape_func_values);
		else if (type == ElementType::Q2D4N)
			shapeQuad4(parametric_coordinate, shape_func_values);
		else if (type == ElementType::T2D6N)
			shapeTri6(parametric_coordinate, shape_func_values);
		else if (type == ElementType::Q2D8N)
			shapeQuad8(parametric_coordinate, shape_func_values);
		else if (type == ElementType::T3D4N)
			shapeTet4(parametric_coordinate, shape_func_values);
		else if (type == ElementType::T3D10N)
			shapeTet10(parametric_coordinate, shape_func_values);
		else if (type == ElementType::H3D8N)
			shapeHex8(parametric_coordinate, shape_func_values);
		else if (type == ElementType::H3D20N)
			shapeHex20(parametric_coordinate, shape_func_values);
		else if (type == ElementType::P3D5N)
			shapePyr5(parametric_coordinate, shape_func_values);
		else if (type == ElementType::P3D13N)
			shapePyr13(parametric_coordinate, shape_func_values);
		else if (type == ElementType::W3D15N)
			shapeWedge15(parametric_coordinate, shape_func_values);
		else
			throw _UTIL_MSG_EXCEPT("Cannot find interpolation functions for element.");
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	void flex_mesh::ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::derivatives(const ElementType type, const Coordinate & parametric_coordinate, ElementVectorDeriv& shape_func_derivatives)
	{
		const auto num_nodes = get_num_nodes_for_type(type);
		shape_func_derivatives.resize(num_nodes, 3);
		shape_func_derivatives.setZero();

		if (type == ElementType::L1D2N)
			shapeDerivLine2(parametric_coordinate, shape_func_derivatives);
		else if (type == ElementType::L1D3N)
			shapeDerivLine3(parametric_coordinate, shape_func_derivatives);
		/*
		else if (type == ElementType::L1D4N)
			shapeDerivLine4(parametric_coordinate, shape_func_derivatives);
		else if (type == ElementType::L1D5N)
			shapeDerivLine5(parametric_coordinate, shape_func_derivatives);
		*/
		else if (type == ElementType::T2D3N)
			shapeDerivTri3(parametric_coordinate, shape_func_derivatives);
		else if (type == ElementType::Q2D4N)
			shapeDerivQuad4(parametric_coordinate, shape_func_derivatives);
		else if (type == ElementType::T2D6N)
			shapeDerivTri6(parametric_coordinate, shape_func_derivatives);
		else if (type == ElementType::Q2D8N)
			shapeDerivQuad8(parametric_coordinate, shape_func_derivatives);
		else if (type == ElementType::T3D4N)
			shapeDerivTet4(parametric_coordinate, shape_func_derivatives);
		else if (type == ElementType::T3D10N)
			shapeDerivTet10(parametric_coordinate, shape_func_derivatives);
		else if (type == ElementType::H3D8N)
			shapeDerivHex8(parametric_coordinate, shape_func_derivatives);
		else if (type == ElementType::H3D20N)
			shapeDerivHex20(parametric_coordinate, shape_func_derivatives);
		else if (type == ElementType::P3D5N)
			shapeDerivPyr5(parametric_coordinate, shape_func_derivatives);
		else if (type == ElementType::P3D13N)
			shapeDerivPyr13(parametric_coordinate, shape_func_derivatives);
		else if (type == ElementType::W3D15N)
			shapeDerivWedge15(parametric_coordinate, shape_func_derivatives);
		else
			throw _UTIL_MSG_EXCEPT("Cannot find interpolation functions for element.");
	}

#ifdef __MSC_VER
#pragma region 1D Shape Function Methods
#endif 

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeLine2(const Coordinate & xi, ElementVector & S)
	{
		S(0) = 1.0 / 2.0 - xi[0] / 2;
		S(1) = 1.0 / 2.0 + xi[0] / 2;
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivLine2(const Coordinate & xi, ElementVectorDeriv & SDeriv)
	{
		// Derivatives w.r.t. xi (column 0)
		SDeriv(0, 0) = -1.0 / 2.0;
		SDeriv(1, 0) = 1.0 / 2.0;
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeLine3(const Coordinate & xi, ElementVector & S)
	{
		const double t1 = xi[0] * xi[0];
		S(0) = -xi[0] / 2 + t1 / 2;
		S(1) = 1.0 - t1;
		S(2) = xi[0] / 2 + t1 / 2;
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivLine3(const Coordinate & xi, ElementVectorDeriv & SDeriv)
	{
		// Derivatives w.r.t. xi (column 0)
		SDeriv(0, 0) = -1.0 / 2.0 + xi[0];
		SDeriv(1, 0) = -2.0*xi[0];
		SDeriv(2, 0) = 1.0 / 2.0 + xi[0];
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeLine4(const Coordinate & xi, ElementVector & S)
	{
		const double t1 = xi[0] * xi[0];
		const double t2 = t1 * xi[0];
		S(0) = -1.0 / 16.0 + xi[0] / 16 + 9.0 / 16.0*t1 - 9.0 / 16.0*t2;
		S(1) = 9.0 / 16.0 - 27.0 / 16.0*xi[0] - 9.0 / 16.0*t1 + 27.0 / 16.0*t2;
		S(2) = 9.0 / 16.0 + 27.0 / 16.0*xi[0] - 9.0 / 16.0*t1 - 27.0 / 16.0*t2;
		S(3) = -1.0 / 16.0 - xi[0] / 16 + 9.0 / 16.0*t1 + 9.0 / 16.0*t2;
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivLine4(const Coordinate & xi, ElementVectorDeriv & SDeriv)
	{
		const double t1 = xi[0] * xi[0];
		// Derivatives w.r.t. xi (column 0)
		SDeriv(0, 0) = 1.0 / 16.0 + 9.0 / 8.0*xi[0] - 27.0 / 16.0*t1;
		SDeriv(1, 0) = -27.0 / 16.0 - 9.0 / 8.0*xi[0] + 81.0 / 16.0*t1;
		SDeriv(2, 0) = 27.0 / 16.0 - 9.0 / 8.0*xi[0] - 81.0 / 16.0*t1;
		SDeriv(3, 0) = -1.0 / 16.0 + 9.0 / 8.0*xi[0] + 27.0 / 16.0*t1;
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeLine5(const Coordinate & xi, ElementVector & S)
	{
		const double t1 = xi[0] * xi[0];
		const double t2 = t1 * xi[0];
		const double t3 = t1 * t1;
		S(0) = xi[0] / 6 - t1 / 6 - 2.0 / 3.0*t2 + 2.0 / 3.0*t3;
		S(1) = -4.0 / 3.0*xi[0] + 8.0 / 3.0*t1 + 4.0 / 3.0*t2 - 8.0 / 3.0*t3;
		S(2) = 1.0 - 5.0*t1 + 4.0*t3;
		S(3) = 4.0 / 3.0*xi[0] + 8.0 / 3.0*t1 - 4.0 / 3.0*t2 - 8.0 / 3.0*t3;
		S(4) = -xi[0] / 6 - t1 / 6 + 2.0 / 3.0*t2 + 2.0 / 3.0*t3;
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivLine5(const Coordinate & xi, ElementVectorDeriv & SDeriv)
	{
		const double t1 = xi[0] * xi[0];
		const double t2 = t1 * xi[0];
		// Derivatives w.r.t. xi (column 0)
		SDeriv(0, 0) = 1.0 / 6.0 - xi[0] / 3 - 2.0*t1 + 8.0 / 3.0*t2;
		SDeriv(1, 0) = -4.0 / 3.0 + 16.0 / 3.0*xi[0] + 4.0*t1 - 32.0 / 3.0*t2;
		SDeriv(2, 0) = -10.0*xi[0] + 16.0*t2;
		SDeriv(3, 0) = 4.0 / 3.0 + 16.0 / 3.0*xi[0] - 4.0*t1 - 32.0 / 3.0*t2;
		SDeriv(4, 0) = -1.0 / 6.0 - xi[0] / 3 + 2.0*t1 + 8.0 / 3.0*t2;
	}

#ifdef __MSC_VER
#pragma endregion

#pragma region 2D Shape Function Methods
#endif 

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeTri3(const Coordinate & p, ElementVector & S)
	{
		const double xi = p[0];
		const double eta = p[1];
		S(0) = 1 - xi - eta;
		S(1) = xi;
		S(2) = eta;
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivTri3(const Coordinate & p, ElementVectorDeriv & SDeriv)
	{
		// Derivatives w.r.t. xi (column 0)
		SDeriv(0, 0) = -1;
		SDeriv(1, 0) = 1;
		SDeriv(2, 0) = 0;
		// Derivatives w.r.t. eta (column 1)
		SDeriv(0, 1) = -1;
		SDeriv(1, 1) = 0;
		SDeriv(2, 1) = 1;
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeTri6(const Coordinate & p, ElementVector & S)
	{
		const double xi = p[0];
		const double eta = p[1];
		const double t1 = xi * xi;
		const double t2 = eta * eta;
		const double t3 = eta * xi;
		S(0) = 1 - 3 * xi - 3 * eta + 4 * t3 + 2 * t1 + 2 * t2;
		S(1) = 4 * xi - 4 * t3 - 4 * t1;
		S(2) = -xi + 2 * t1;
		S(3) = 4 * t3;
		S(4) = -eta + 2 * t2;
		S(5) = 4 * eta - 4 * t3 - 4 * t2;
	}


	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivTri6(const Coordinate & p, ElementVectorDeriv & SDeriv)
	{
		const double xi = p[0];
		const double eta = p[1];
		// Derivatives w.r.t. xi (column 0)
		SDeriv(0, 0) = -3 + 4 * eta + 4 * xi;
		SDeriv(1, 0) = 4 - 4 * eta - 8 * xi;
		SDeriv(2, 0) = -1 + 4 * xi;
		SDeriv(3, 0) = 4 * eta;
		SDeriv(4, 0) = 0;
		SDeriv(5, 0) = -4 * eta;
		// Derivatives w.r.t. eta (column 1)
		SDeriv(0, 1) = -3 + 4 * eta + 4 * xi;
		SDeriv(1, 1) = -4 * xi;
		SDeriv(2, 1) = 0;
		SDeriv(3, 1) = 4 * xi;
		SDeriv(4, 1) = -1 + 4 * eta;
		SDeriv(5, 1) = 4 - 4 * xi - 8 * eta;
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeQuad4(const Coordinate & p, ElementVector & S)
	{
		const double xi = p[0];
		const double eta = p[1];
		S(0) = .25 * (1 - xi) * (1 - eta);
		S(1) = .25 * (1 + xi) * (1 - eta);
		S(2) = .25 * (1 + xi) * (1 + eta);
		S(3) = .25 * (1 - xi) * (1 + eta);
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivQuad4(const Coordinate & p, ElementVectorDeriv & SDeriv)
	{
		const double xi = p[0];
		const double eta = p[1];

		// Derivatives w.r.t. xi (column 0)
		SDeriv(0, 0) = -.25*(1 - eta);
		SDeriv(1, 0) = .25*(1 - eta);
		SDeriv(2, 0) = .25*(1 + eta);
		SDeriv(3, 0) = -.25*(1 + eta);
		// Derivatives w.r.t. eta (column 1)
		SDeriv(0, 1) = -.25*(1 - xi);
		SDeriv(1, 1) = -.25*(1 + xi);
		SDeriv(2, 1) = .25*(1 + xi);
		SDeriv(3, 1) = .25*(1 - xi);
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeQuad8(const Coordinate & p, ElementVector & S)
	{
		const double xi = p[0];
		const double eta = p[1];
		S(0) = (-1 + eta)*(1 + eta - eta * xi - xi * xi) / 4;
		S(1) = (-1 + eta)*(-1 + xi * xi) / 2;
		S(2) = (-1 + eta)*(1 + eta + eta * xi - xi * xi) / 4;
		S(3) = (1 - eta * eta)*(1 + xi) / 2;
		S(4) = (1 + eta)*(-1 + eta + eta * xi + xi * xi) / 4;
		S(5) = (1 + eta)*(1 - xi * xi) / 2;
		S(6) = (1 + eta)*(-1 + eta - eta * xi + xi * xi) / 4;
		S(7) = (-1 + eta * eta)*(-1 + xi) / 2;
	}


	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivQuad8(const Coordinate & p, ElementVectorDeriv & SDeriv)
	{
		const double xi = p[0];
		const double eta = p[1];

		// Derivatives w.r.t. xi (column 0)
		SDeriv(0, 0) = -((1 - eta)*(1 - xi)) / 4 - (1 - eta)*(-1 - eta - xi) / 4;
		SDeriv(1, 0) = -((1 - eta)*xi);
		SDeriv(2, 0) = (1 - eta)*(1 + xi) / 4 + (1 - eta)*(-1 - eta + xi) / 4;
		SDeriv(3, 0) = (1 - eta * eta) / 2;
		SDeriv(4, 0) = (1 + eta)*(1 + xi) / 4 + (1 + eta)*(-1 + eta + xi) / 4;
		SDeriv(5, 0) = -((1 + eta)*xi);
		SDeriv(6, 0) = -((1 + eta)*(1 - xi)) / 4 - (1 + eta)*(-1 + eta - xi) / 4;
		SDeriv(7, 0) = -(1 - eta * eta) / 2;
		// Derivatives w.r.t. eta (column 1)
		SDeriv(0, 1) = -((1 - eta)*(1 - xi)) / 4 - (1 - xi)*(-1 - eta - xi) / 4;
		SDeriv(1, 1) = -(1 - xi * xi) / 2;
		SDeriv(2, 1) = -((1 - eta)*(1 + xi)) / 4 - (1 + xi)*(-1 - eta + xi) / 4;
		SDeriv(3, 1) = -(eta*(1 + xi));
		SDeriv(4, 1) = (1 + eta)*(1 + xi) / 4 + (1 + xi)*(-1 + eta + xi) / 4;
		SDeriv(5, 1) = (1 - xi * xi) / 2;
		SDeriv(6, 1) = (1 + eta)*(1 - xi) / 4 + (1 - xi)*(-1 + eta - xi) / 4;
		SDeriv(7, 1) = -(eta*(1 - xi));
	}

#ifdef __MSC_VER
#pragma endregion

#pragma region 3D Shape Function Methods
#endif 

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapePyr5(const Coordinate & p, ElementVector & S)
	{
		// Natural coordinates of nodes
		const std::array<double, 5> xi_n = { -1.0, 1.0, 1.0, -1.0, 0.0 };
		const std::array<double, 5> eta_n = { -1.0, -1.0, 1.0, 1.0, 0.0 };
		const std::array<double, 5> mu_n = { -1.0, -1.0, -1.0, -1.0, 1.0 };

		const double xi = p[0];
		const double eta = p[1];
		const double mu = p[2];
		// Shape function values
		for (auto i = 0; i < 4; i++)
			S(i) = 1.0 / 8.0*(1.0 + xi * xi_n[i])*(1.0 + eta * eta_n[i])*(1.0 + mu * mu_n[i]);
		S(4) = 1.0 / 2.0 * (1.0 + mu);
	}


	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivPyr5(const Coordinate & p, ElementVectorDeriv & SDeriv)
	{
		const double xi = p[0];
		const double eta = p[1];
		const double mu = p[2];

		// Derivatives
		SDeriv(0, 0) = -1.0 / 8.0 * (1.0 - eta)*(1.0 - mu);
		SDeriv(0, 1) = -1.0 / 8.0 * (1.0 - xi)*(1.0 - mu);
		SDeriv(0, 2) = -1.0 / 8.0 * (1.0 - xi)*(1.0 - eta);
		SDeriv(1, 0) = 1.0 / 8.0 * (1.0 - eta)*(1.0 - mu);
		SDeriv(1, 1) = -1.0 / 8.0 * (1.0 + xi)*(1.0 - mu);
		SDeriv(1, 2) = -1.0 / 8.0 * (1.0 + xi)*(1.0 - eta);
		SDeriv(2, 0) = 1.0 / 8.0 * (1.0 + eta)*(1.0 - mu);
		SDeriv(2, 1) = 1.0 / 8.0 * (1.0 + xi)*(1.0 - mu);
		SDeriv(2, 2) = -1.0 / 8.0 * (1.0 + xi)*(1.0 + eta);
		SDeriv(3, 0) = -1.0 / 8.0 * (1.0 + eta)*(1.0 - mu);
		SDeriv(3, 1) = 1.0 / 8.0 * (1.0 - xi)*(1.0 - mu);
		SDeriv(3, 2) = -1.0 / 8.0 * (1.0 - xi)*(1.0 + eta);
		SDeriv(4, 0) = 0.0;
		SDeriv(4, 1) = 0.0;
		SDeriv(4, 2) = 0.5;
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapePyr13(const Coordinate & p, ElementVector & S)
	{

		const double xi = p[0];
		const double eta = p[1];
		const double mu = p[2];

		// Shape function values
		// Base nodes
		S(0) = -1.0 / 16.0*(1.0 - xi)*(1.0 - eta)*(1.0 - mu)*(4.0 + 3.0*xi + 3.0*eta + 2.0*xi*eta + 2.0*mu + xi * mu + eta * mu + 2.0*xi*eta*mu);
		S(1) = -1.0 / 16.0*(1.0 + xi)*(1.0 - eta)*(1.0 - mu)*(4.0 - 3.0*xi + 3.0*eta - 2.0*xi*eta + 2.0*mu - xi * mu + eta * mu - 2.0*xi*eta*mu);
		S(2) = -1.0 / 16.0*(1.0 + xi)*(1.0 + eta)*(1.0 - mu)*(4.0 - 3.0*xi - 3.0*eta + 2.0*xi*eta + 2.0*mu - xi * mu - eta * mu + 2.0*xi*eta*mu);
		S(3) = -1.0 / 16.0*(1.0 - xi)*(1.0 + eta)*(1.0 - mu)*(4.0 + 3.0*xi - 3.0*eta - 2.0*xi*eta + 2.0*mu + xi * mu - eta * mu - 2.0*xi*eta*mu);
		// Tip
		S(4) = 1.0 / 2.0*mu*(1.0 + mu);
		// Mid side nodes
		S(5) = 1.0 / 8.0*mu*(1.0 - xi * xi)*(1.0 - eta)*(1.0 - mu)*(2.0 + eta + eta * mu);
		S(6) = 1.0 / 8.0*mu*(1.0 + xi)*(1.0 - eta * eta)*(1.0 - mu)*(2.0 - xi - xi * mu);
		S(7) = 1.0 / 8.0*mu*(1.0 - xi * xi)*(1.0 + eta)*(1.0 - mu)*(2.0 - eta - eta * mu);
		S(8) = 1.0 / 8.0*mu*(1.0 - xi)*(1.0 - eta * eta)*(1.0 - mu)*(2.0 + xi + xi * mu);
		S(9) = 1.0 / 4.0*(1.0 - xi)*(1.0 - eta)*(1.0 - mu * mu);
		S(10) = 1.0 / 4.0*(1.0 + xi)*(1.0 - eta)*(1.0 - mu * mu);
		S(11) = 1.0 / 4.0*(1.0 + xi)*(1.0 + eta)*(1.0 - mu * mu);
		S(12) = 1.0 / 4.0*(1.0 - xi)*(1.0 + eta)*(1.0 - mu * mu);
	}


	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivPyr13(const Coordinate & p, ElementVectorDeriv & SDeriv)
	{

		const double xi = p[0];
		const double eta = p[1];
		const double mu = p[2];

		// Derivatives
		// w.r.t xi
		// Base nodes
		SDeriv(0, 0) = 1.0 / 16.0*(1.0 - eta)*(1.0 - mu)*(1.0 + 6.0*xi + eta + 4.0*xi*eta + mu + 2.0*xi*mu - eta * mu + 4.0*xi*eta*mu);
		SDeriv(1, 0) = -1.0 / 16.0*(1.0 - eta)*(1.0 - mu)*(1.0 - 6.0*xi + eta - 4.0*xi*eta + mu - 2.0*xi*mu - eta * mu - 4.0*xi*eta*mu);
		SDeriv(2, 0) = -1.0 / 16.0*(1.0 + eta)*(1.0 - mu)*(1.0 - 6.0*xi - eta + 4.0*xi*eta + mu - 2.0*xi*mu + eta * mu + 4.0*xi*eta*mu);
		SDeriv(3, 0) = 1.0 / 16.0*(1.0 + eta)*(1.0 - mu)*(1.0 + 6.0*xi - eta - 4.0*xi*eta + mu + 2.0*xi*mu + eta * mu - 4.0*xi*eta*mu);
		// Tip
		SDeriv(4, 0) = 0.0;
		// Mid side nodes
		SDeriv(5, 0) = -1.0 / 4.0*xi*(1.0 - eta)*(1.0 - mu)*(2.0 + eta + eta * mu);
		SDeriv(6, 0) = 1.0 / 8.0*(1.0 - eta * eta)*(1.0 - mu)*(1.0 - 2.0*xi - mu - 2.0*xi*mu);
		SDeriv(7, 0) = -1.0 / 4.0*xi*(1.0 + eta)*(1.0 - mu)*(2.0 - eta - eta * mu);
		SDeriv(8, 0) = -1.0 / 8.0*(1.0 - eta * eta)*(1.0 - mu)*(1.0 + 2.0*xi - mu + 2.0*xi*mu);
		SDeriv(9, 0) = -1.0 / 4.0*(1.0 - eta)*(1.0 - mu * mu);
		SDeriv(10, 0) = 1.0 / 4.0*(1.0 - eta)*(1.0 - mu * mu);
		SDeriv(11, 0) = 1.0 / 4.0*(1.0 + eta)*(1.0 - mu * mu);
		SDeriv(12, 0) = -1.0 / 4.0*(1.0 + eta)*(1.0 - mu * mu);
		// w.r.t eta
		// Base nodes
		SDeriv(0, 1) = 1.0 / 16.0*(1.0 - xi)*(1.0 - mu)*(1.0 + xi + 6.0*eta + 4.0*xi*eta + mu - xi * mu + 2.0*eta*mu + 4.0*xi*eta*mu);
		SDeriv(1, 1) = 1.0 / 16.0*(1.0 + xi)*(1.0 - mu)*(1.0 - xi + 6.0*eta - 4.0*xi*eta + mu + xi * mu + 2.0*eta*mu - 4.0*xi*eta*mu);
		SDeriv(2, 1) = -1.0 / 16.0*(1.0 + xi)*(1.0 - mu)*(1.0 - xi - 6.0*eta + 4.0*xi*eta + mu + xi * mu - 2.0*eta*mu + 4.0*xi*eta*mu);
		SDeriv(3, 1) = -1.0 / 16.0*(1.0 - xi)*(1.0 - mu)*(1.0 + xi - 6.0*eta - 4.0*xi*eta + mu - xi * mu - 2.0*eta*mu - 4.0*xi*eta*mu);
		// Tip
		SDeriv(4, 1) = 0.0;
		// Mid side nodes
		SDeriv(5, 1) = -1.0 / 8.0*(1.0 - xi * xi)*(1.0 - mu)*(1.0 + 2.0*eta - mu + 2.0*eta*mu);
		SDeriv(6, 1) = -1.0 / 4.0*eta*(1.0 + xi)*(1.0 - mu)*(2.0 - xi - xi * mu);
		SDeriv(7, 1) = 1.0 / 8.0*(1.0 - xi * xi)*(1.0 - mu)*(1.0 - 2.0*eta - mu - 2.0*eta*mu);
		SDeriv(8, 1) = -1.0 / 4.0*eta*(1.0 - xi)*(1.0 - mu)*(2.0 + xi + xi * mu);
		SDeriv(9, 1) = -1.0 / 4.0*(1.0 - xi)*(1.0 - mu * mu);
		SDeriv(10, 1) = -1.0 / 4.0*(1.0 + xi)*(1.0 - mu * mu);
		SDeriv(11, 1) = 1.0 / 4.0*(1.0 + xi)*(1.0 - mu * mu);
		SDeriv(12, 1) = 1.0 / 4.0*(1.0 - xi)*(1.0 - mu * mu);
		// w.r.t mu
		// Base nodes
		SDeriv(0, 2) = 1.0 / 8.0*(1.0 - xi)*(1.0 - eta)*(1.0 + xi + eta + 2.0*mu + xi * mu + eta * mu + 2.0*xi*eta*mu);
		SDeriv(1, 2) = 1.0 / 8.0*(1.0 + xi)*(1.0 - eta)*(1.0 - xi + eta + 2.0*mu - xi * mu + eta * mu - 2.0*xi*eta*mu);
		SDeriv(2, 2) = 1.0 / 8.0*(1.0 + xi)*(1.0 + eta)*(1.0 - xi - eta + 2.0*mu - xi * mu - eta * mu + 2.0*xi*eta*mu);
		SDeriv(3, 2) = 1.0 / 8.0*(1.0 - xi)*(1.0 + eta)*(1.0 + xi - eta + 2.0*mu + xi * mu - eta * mu - 2.0*xi*eta*mu);
		// Tip
		SDeriv(4, 2) = 0.5 + mu;
		// Mid side nodes
		SDeriv(5, 2) = -1.0 / 4.0*(1.0 - xi * xi)*(1.0 - eta)*(1.0 + eta * mu);
		SDeriv(6, 2) = -1.0 / 4.0*(1.0 + xi)*(1.0 - eta * eta)*(1.0 - xi * mu);
		SDeriv(7, 2) = -1.0 / 4.0*(1.0 - xi * xi)*(1.0 + eta)*(1.0 - eta * mu);
		SDeriv(8, 2) = -1.0 / 4.0*(1.0 - xi)*(1.0 - eta * eta)*(1.0 + xi * mu);
		SDeriv(9, 2) = -0.5*(1.0 - xi)*(1.0 - eta)*mu;
		SDeriv(10, 2) = -0.5*(1.0 + xi)*(1.0 - eta)*mu;
		SDeriv(11, 2) = -0.5*(1.0 + xi)*(1.0 + eta)*mu;
		SDeriv(12, 2) = -0.5*(1.0 - xi)*(1.0 + eta)*mu;
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeHex8(const Coordinate & p, ElementVector & S)
	{

		const double xi = p[0];
		const double eta = p[1];
		const double zeta = p[2];
		{
			const double t1 = xi * eta;
			const double t2 = xi * zeta;
			const double t3 = eta * zeta;
			const double t4 = t1 * zeta;
			S(0) = 1.0 / 8.0 - xi / 8 - eta / 8 + t1 / 8 - zeta / 8 + t2 / 8 + t3 / 8 - t4 / 8;
			S(1) = 1.0 / 8.0 - xi / 8 + eta / 8 - t1 / 8 - zeta / 8 + t2 / 8 - t3 / 8 + t4 / 8;
			S(2) = 1.0 / 8.0 - xi / 8 + eta / 8 - t1 / 8 + zeta / 8 - t2 / 8 + t3 / 8 - t4 / 8;
			S(3) = 1.0 / 8.0 - xi / 8 - eta / 8 + t1 / 8 + zeta / 8 - t2 / 8 - t3 / 8 + t4 / 8;
			S(4) = 1.0 / 8.0 + xi / 8 - eta / 8 - t1 / 8 - zeta / 8 - t2 / 8 + t3 / 8 + t4 / 8;
			S(5) = 1.0 / 8.0 + xi / 8 + eta / 8 + t1 / 8 - zeta / 8 - t2 / 8 - t3 / 8 - t4 / 8;
			S(6) = 1.0 / 8.0 + xi / 8 + eta / 8 + t1 / 8 + zeta / 8 + t2 / 8 + t3 / 8 + t4 / 8;
			S(7) = 1.0 / 8.0 + xi / 8 - eta / 8 - t1 / 8 + zeta / 8 + t2 / 8 - t3 / 8 - t4 / 8;
		}
	}


	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivHex8(const Coordinate & p, ElementVectorDeriv & SDeriv)
	{

		const double xi = p[0];
		const double eta = p[1];
		const double zeta = p[2];
		// Derivatives w.r.t. xi (column 0)
		{
			const double t1 = eta * zeta;
			SDeriv(0, 0) = -1.0 / 8.0 + eta / 8 + zeta / 8 - t1 / 8;
			SDeriv(1, 0) = -1.0 / 8.0 - eta / 8 + zeta / 8 + t1 / 8;
			SDeriv(2, 0) = -1.0 / 8.0 - eta / 8 - zeta / 8 - t1 / 8;
			SDeriv(3, 0) = -1.0 / 8.0 + eta / 8 - zeta / 8 + t1 / 8;
			SDeriv(4, 0) = 1.0 / 8.0 - eta / 8 - zeta / 8 + t1 / 8;
			SDeriv(5, 0) = 1.0 / 8.0 + eta / 8 - zeta / 8 - t1 / 8;
			SDeriv(6, 0) = 1.0 / 8.0 + eta / 8 + zeta / 8 + t1 / 8;
			SDeriv(7, 0) = 1.0 / 8.0 - eta / 8 + zeta / 8 - t1 / 8;
		}
		// Derivatives w.r.t. eta (column 1)
		{
			const double t1 = xi * zeta;
			SDeriv(0, 1) = -1.0 / 8.0 + xi / 8 + zeta / 8 - t1 / 8;
			SDeriv(1, 1) = 1.0 / 8.0 - xi / 8 - zeta / 8 + t1 / 8;
			SDeriv(2, 1) = 1.0 / 8.0 - xi / 8 + zeta / 8 - t1 / 8;
			SDeriv(3, 1) = -1.0 / 8.0 + xi / 8 - zeta / 8 + t1 / 8;
			SDeriv(4, 1) = -1.0 / 8.0 - xi / 8 + zeta / 8 + t1 / 8;
			SDeriv(5, 1) = 1.0 / 8.0 + xi / 8 - zeta / 8 - t1 / 8;
			SDeriv(6, 1) = 1.0 / 8.0 + xi / 8 + zeta / 8 + t1 / 8;
			SDeriv(7, 1) = -1.0 / 8.0 - xi / 8 - zeta / 8 - t1 / 8;
		}
		// Derivatives w.r.t. zeta (column 2)
		{
			const double t1 = xi * eta;
			SDeriv(0, 2) = -1.0 / 8.0 + xi / 8 + eta / 8 - t1 / 8;
			SDeriv(1, 2) = -1.0 / 8.0 + xi / 8 - eta / 8 + t1 / 8;
			SDeriv(2, 2) = 1.0 / 8.0 - xi / 8 + eta / 8 - t1 / 8;
			SDeriv(3, 2) = 1.0 / 8.0 - xi / 8 - eta / 8 + t1 / 8;
			SDeriv(4, 2) = -1.0 / 8.0 - xi / 8 + eta / 8 + t1 / 8;
			SDeriv(5, 2) = -1.0 / 8.0 - xi / 8 - eta / 8 - t1 / 8;
			SDeriv(6, 2) = 1.0 / 8.0 + xi / 8 + eta / 8 + t1 / 8;
			SDeriv(7, 2) = 1.0 / 8.0 + xi / 8 - eta / 8 - t1 / 8;
		}
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeHex20(const Coordinate & p, ElementVector & S)
	{
		const double xi = p[0];
		const double eta = p[1];
		const double zeta = p[2];

		// S
		{
			const double t1 = xi * xi;
			const double t2 = eta * eta;
			const double t3 = zeta * zeta;
			const double t4 = xi * t2;
			const double t5 = t1 * eta;
			const double t6 = eta * t3;
			const double t7 = t2 * zeta;
			const double t8 = xi * t3;
			const double t9 = t1 * zeta;
			const double t10 = xi * eta;
			const double t11 = t10 * t3;
			const double t12 = t10 * zeta;
			const double t13 = t4 * zeta;
			const double t14 = t5 * zeta;
			const double t15 = -1.0 / 4.0 + t1 / 8 + eta / 8 + t2 / 8 + t3 / 8 + xi / 8 + zeta / 8 - t4 / 8 - t5 / 8 - t6 / 8 - t7 / 8 - t8 / 8 -
				t9 / 8 + t11 / 8 - t12 / 8 + t13 / 8 + t14 / 8;
			const double t16 = xi * zeta;
			const double t18 = -1.0 / 4.0 + t1 / 8 - eta / 8 + t2 / 8 + t3 / 8 + xi / 8 + zeta / 8 - t4 / 8 + t5 / 8 + t6 / 8 - t7 / 8 - t8 / 8 -
				t9 / 8 - t11 / 8 + t12 / 8 + t13 / 8 - t14 / 8;
			const double t20 = -1.0 / 4.0 + t1 / 8 - eta / 8 + t2 / 8 + t3 / 8 + xi / 8 - zeta / 8 - t4 / 8 + t5 / 8 + t6 / 8 + t7 / 8 - t8 / 8 +
				t9 / 8 - t11 / 8 - t12 / 8 - t13 / 8 + t14 / 8;
			const double t22 = -1.0 / 4.0 + t1 / 8 + eta / 8 + t2 / 8 + t3 / 8 + xi / 8 - zeta / 8 - t4 / 8 - t5 / 8 - t6 / 8 + t7 / 8 - t8 / 8 +
				t9 / 8 + t11 / 8 + t12 / 8 - t13 / 8 - t14 / 8;
			const double t24 = eta * zeta;
			const double t29 = -1.0 / 4.0 + t1 / 8 + eta / 8 + t2 / 8 + t3 / 8 - xi / 8 + zeta / 8 + t4 / 8 - t5 / 8 - t6 / 8 - t7 / 8 + t8 / 8 -
				t9 / 8 - t11 / 8 + t12 / 8 - t13 / 8 + t14 / 8;
			const double t31 = -1.0 / 4.0 + t1 / 8 - eta / 8 + t2 / 8 + t3 / 8 - xi / 8 + zeta / 8 + t4 / 8 + t5 / 8 + t6 / 8 - t7 / 8 + t8 / 8 -
				t9 / 8 + t11 / 8 - t12 / 8 - t13 / 8 - t14 / 8;
			const double t33 = -1.0 / 4.0 + t1 / 8 - eta / 8 + t2 / 8 + t3 / 8 - xi / 8 - zeta / 8 + t4 / 8 + t5 / 8 + t6 / 8 + t7 / 8 + t8 / 8 +
				t9 / 8 + t11 / 8 + t12 / 8 + t13 / 8 + t14 / 8;
			const double t35 = -1.0 / 4.0 + t1 / 8 + eta / 8 + t2 / 8 + t3 / 8 - xi / 8 - zeta / 8 + t4 / 8 - t5 / 8 - t6 / 8 + t7 / 8 + t8 / 8 +
				t9 / 8 - t11 / 8 - t12 / 8 + t13 / 8 - t14 / 8;
			S(0) = t15;
			S(1) = 1.0 / 4.0 - t2 / 4 - xi / 4 - zeta / 4 + t4 / 4 + t7 / 4 + t16 / 4 - t13 / 4;
			S(2) = t18;
			S(3) = 1.0 / 4.0 + eta / 4 - t3 / 4 - xi / 4 - t10 / 4 - t6 / 4 + t8 / 4 + t11 / 4;
			S(4) = t20;
			S(5) = 1.0 / 4.0 - t2 / 4 - xi / 4 + zeta / 4 + t4 / 4 - t7 / 4 - t16 / 4 + t13 / 4;
			S(6) = t22;
			S(7) = 1.0 / 4.0 - eta / 4 - t3 / 4 - xi / 4 + t10 / 4 + t6 / 4 + t8 / 4 - t11 / 4;
			S(8) = 1.0 / 4.0 - t1 / 4 - eta / 4 - zeta / 4 + t5 / 4 + t24 / 4 + t9 / 4 - t14 / 4;
			S(9) = 1.0 / 4.0 - t1 / 4 + eta / 4 - zeta / 4 - t5 / 4 - t24 / 4 + t9 / 4 + t14 / 4;
			S(10) = 1.0 / 4.0 - t1 / 4 + eta / 4 + zeta / 4 - t5 / 4 + t24 / 4 - t9 / 4 - t14 / 4;
			S(11) = 1.0 / 4.0 - t1 / 4 - eta / 4 + zeta / 4 + t5 / 4 - t24 / 4 - t9 / 4 + t14 / 4;
			S(12) = t29;
			S(13) = 1.0 / 4.0 - t2 / 4 + xi / 4 - zeta / 4 - t4 / 4 + t7 / 4 - t16 / 4 + t13 / 4;
			S(14) = t31;
			S(15) = 1.0 / 4.0 + eta / 4 - t3 / 4 + xi / 4 + t10 / 4 - t6 / 4 - t8 / 4 - t11 / 4;
			S(16) = t33;
			S(17) = 1.0 / 4.0 - t2 / 4 + xi / 4 + zeta / 4 - t4 / 4 - t7 / 4 + t16 / 4 - t13 / 4;
			S(18) = t35;
			S(19) = 1.0 / 4.0 - eta / 4 - t3 / 4 + xi / 4 - t10 / 4 + t6 / 4 - t8 / 4 + t11 / 4;
		}
	}


	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivHex20(const Coordinate & p, ElementVectorDeriv & SDeriv)
	{
		const double xi = p[0];
		const double eta = p[1];
		const double zeta = p[2];

		// Derivatives w.r.t. xi (column 0)
		{
			const double t1 = eta * eta;
			const double t2 = xi * eta;
			const double t3 = zeta * zeta;
			const double t4 = xi * zeta;
			const double t5 = eta * t3;
			const double t6 = eta * zeta;
			const double t7 = t1 * zeta;
			const double t8 = t2 * zeta;
			SDeriv(0, 0) = xi / 4 + 1.0 / 8.0 - t1 / 8 - t2 / 4 - t3 / 8 - t4 / 4 + t5 / 8 - t6 / 8 + t7 / 8 + t8 / 4;
			SDeriv(1, 0) = -1.0 / 4.0 + t1 / 4 + zeta / 4 - t7 / 4;
			SDeriv(2, 0) = xi / 4 + 1.0 / 8.0 - t1 / 8 + t2 / 4 - t3 / 8 - t4 / 4 - t5 / 8 + t6 / 8 + t7 / 8 - t8 / 4;
			SDeriv(3, 0) = -1.0 / 4.0 - eta / 4 + t3 / 4 + t5 / 4;
			SDeriv(4, 0) = xi / 4 + 1.0 / 8.0 - t1 / 8 + t2 / 4 - t3 / 8 + t4 / 4 - t5 / 8 - t6 / 8 - t7 / 8 + t8 / 4;
			SDeriv(5, 0) = -1.0 / 4.0 + t1 / 4 - zeta / 4 + t7 / 4;
			SDeriv(6, 0) = xi / 4 + 1.0 / 8.0 - t1 / 8 - t2 / 4 - t3 / 8 + t4 / 4 + t5 / 8 + t6 / 8 - t7 / 8 - t8 / 4;
			SDeriv(7, 0) = -1.0 / 4.0 + eta / 4 + t3 / 4 - t5 / 4;
			SDeriv(8, 0) = -xi / 2 + t2 / 2 + t4 / 2 - t8 / 2;
			SDeriv(9, 0) = -xi / 2 - t2 / 2 + t4 / 2 + t8 / 2;
			SDeriv(10, 0) = -xi / 2 - t2 / 2 - t4 / 2 - t8 / 2;
			SDeriv(11, 0) = -xi / 2 + t2 / 2 - t4 / 2 + t8 / 2;
			SDeriv(12, 0) = xi / 4 - 1.0 / 8.0 + t1 / 8 - t2 / 4 + t3 / 8 - t4 / 4 - t5 / 8 + t6 / 8 - t7 / 8 + t8 / 4;
			SDeriv(13, 0) = 1.0 / 4.0 - t1 / 4 - zeta / 4 + t7 / 4;
			SDeriv(14, 0) = xi / 4 - 1.0 / 8.0 + t1 / 8 + t2 / 4 + t3 / 8 - t4 / 4 + t5 / 8 - t6 / 8 - t7 / 8 - t8 / 4;
			SDeriv(15, 0) = 1.0 / 4.0 + eta / 4 - t3 / 4 - t5 / 4;
			SDeriv(16, 0) = xi / 4 - 1.0 / 8.0 + t1 / 8 + t2 / 4 + t3 / 8 + t4 / 4 + t5 / 8 + t6 / 8 + t7 / 8 + t8 / 4;
			SDeriv(17, 0) = 1.0 / 4.0 - t1 / 4 + zeta / 4 - t7 / 4;
			SDeriv(18, 0) = xi / 4 - 1.0 / 8.0 + t1 / 8 - t2 / 4 + t3 / 8 + t4 / 4 - t5 / 8 - t6 / 8 + t7 / 8 - t8 / 4;
			SDeriv(19, 0) = 1.0 / 4.0 - eta / 4 - t3 / 4 + t5 / 4;
		}
		// Derivatives w.r.t. eta (column 1)
		{
			const double t1 = xi * eta;
			const double t2 = xi * xi;
			const double t3 = zeta * zeta;
			const double t4 = eta * zeta;
			const double t5 = xi * t3;
			const double t6 = xi * zeta;
			const double t7 = t1 * zeta;
			const double t8 = t2 * zeta;
			SDeriv(0, 1) = 1.0 / 8.0 + eta / 4 - t1 / 4 - t2 / 8 - t3 / 8 - t4 / 4 + t5 / 8 - t6 / 8 + t7 / 4 + t8 / 8;
			SDeriv(1, 1) = -eta / 2 + t1 / 2 + t4 / 2 - t7 / 2;
			SDeriv(2, 1) = -1.0 / 8.0 + eta / 4 - t1 / 4 + t2 / 8 + t3 / 8 - t4 / 4 - t5 / 8 + t6 / 8 + t7 / 4 - t8 / 8;
			SDeriv(3, 1) = 1.0 / 4.0 - xi / 4 - t3 / 4 + t5 / 4;
			SDeriv(4, 1) = -1.0 / 8.0 + eta / 4 - t1 / 4 + t2 / 8 + t3 / 8 + t4 / 4 - t5 / 8 - t6 / 8 - t7 / 4 + t8 / 8;
			SDeriv(5, 1) = -eta / 2 + t1 / 2 - t4 / 2 + t7 / 2;
			SDeriv(6, 1) = 1.0 / 8.0 + eta / 4 - t1 / 4 - t2 / 8 - t3 / 8 + t4 / 4 + t5 / 8 + t6 / 8 - t7 / 4 - t8 / 8;
			SDeriv(7, 1) = -1.0 / 4.0 + xi / 4 + t3 / 4 - t5 / 4;
			SDeriv(8, 1) = -1.0 / 4.0 + t2 / 4 + zeta / 4 - t8 / 4;
			SDeriv(9, 1) = 1.0 / 4.0 - t2 / 4 - zeta / 4 + t8 / 4;
			SDeriv(10, 1) = 1.0 / 4.0 - t2 / 4 + zeta / 4 - t8 / 4;
			SDeriv(11, 1) = -1.0 / 4.0 + t2 / 4 - zeta / 4 + t8 / 4;
			SDeriv(12, 1) = 1.0 / 8.0 + eta / 4 + t1 / 4 - t2 / 8 - t3 / 8 - t4 / 4 - t5 / 8 + t6 / 8 - t7 / 4 + t8 / 8;
			SDeriv(13, 1) = -eta / 2 - t1 / 2 + t4 / 2 + t7 / 2;
			SDeriv(14, 1) = -1.0 / 8.0 + eta / 4 + t1 / 4 + t2 / 8 + t3 / 8 - t4 / 4 + t5 / 8 - t6 / 8 - t7 / 4 - t8 / 8;
			SDeriv(15, 1) = 1.0 / 4.0 + xi / 4 - t3 / 4 - t5 / 4;
			SDeriv(16, 1) = -1.0 / 8.0 + eta / 4 + t1 / 4 + t2 / 8 + t3 / 8 + t4 / 4 + t5 / 8 + t6 / 8 + t7 / 4 + t8 / 8;
			SDeriv(17, 1) = -eta / 2 - t1 / 2 - t4 / 2 - t7 / 2;
			SDeriv(18, 1) = 1.0 / 8.0 + eta / 4 + t1 / 4 - t2 / 8 - t3 / 8 + t4 / 4 - t5 / 8 - t6 / 8 + t7 / 4 - t8 / 8;
			SDeriv(19, 1) = -1.0 / 4.0 - xi / 4 + t3 / 4 + t5 / 4;
		}
		// Derivatives w.r.t. zeta (column 2)
		{
			const double t1 = eta * zeta;
			const double t2 = eta * eta;
			const double t3 = xi * zeta;
			const double t4 = xi * xi;
			const double t5 = xi * eta;
			const double t6 = t5 * zeta;
			const double t7 = xi * t2;
			const double t8 = t4 * eta;
			SDeriv(0, 2) = zeta / 4 + 1.0 / 8.0 - t1 / 4 - t2 / 8 - t3 / 4 - t4 / 8 + t6 / 4 - t5 / 8 + t7 / 8 + t8 / 8;
			SDeriv(1, 2) = -1.0 / 4.0 + t2 / 4 + xi / 4 - t7 / 4;
			SDeriv(2, 2) = zeta / 4 + 1.0 / 8.0 + t1 / 4 - t2 / 8 - t3 / 4 - t4 / 8 - t6 / 4 + t5 / 8 + t7 / 8 - t8 / 8;
			SDeriv(3, 2) = -zeta / 2 - t1 / 2 + t3 / 2 + t6 / 2;
			SDeriv(4, 2) = zeta / 4 - 1.0 / 8.0 + t1 / 4 + t2 / 8 - t3 / 4 + t4 / 8 - t6 / 4 - t5 / 8 - t7 / 8 + t8 / 8;
			SDeriv(5, 2) = 1.0 / 4.0 - t2 / 4 - xi / 4 + t7 / 4;
			SDeriv(6, 2) = zeta / 4 - 1.0 / 8.0 - t1 / 4 + t2 / 8 - t3 / 4 + t4 / 8 + t6 / 4 + t5 / 8 - t7 / 8 - t8 / 8;
			SDeriv(7, 2) = -zeta / 2 + t1 / 2 + t3 / 2 - t6 / 2;
			SDeriv(8, 2) = -1.0 / 4.0 + eta / 4 + t4 / 4 - t8 / 4;
			SDeriv(9, 2) = -1.0 / 4.0 - eta / 4 + t4 / 4 + t8 / 4;
			SDeriv(10, 2) = 1.0 / 4.0 + eta / 4 - t4 / 4 - t8 / 4;
			SDeriv(11, 2) = 1.0 / 4.0 - eta / 4 - t4 / 4 + t8 / 4;
			SDeriv(12, 2) = zeta / 4 + 1.0 / 8.0 - t1 / 4 - t2 / 8 + t3 / 4 - t4 / 8 - t6 / 4 + t5 / 8 - t7 / 8 + t8 / 8;
			SDeriv(13, 2) = -1.0 / 4.0 + t2 / 4 - xi / 4 + t7 / 4;
			SDeriv(14, 2) = zeta / 4 + 1.0 / 8.0 + t1 / 4 - t2 / 8 + t3 / 4 - t4 / 8 + t6 / 4 - t5 / 8 - t7 / 8 - t8 / 8;
			SDeriv(15, 2) = -zeta / 2 - t1 / 2 - t3 / 2 - t6 / 2;
			SDeriv(16, 2) = zeta / 4 - 1.0 / 8.0 + t1 / 4 + t2 / 8 + t3 / 4 + t4 / 8 + t6 / 4 + t5 / 8 + t7 / 8 + t8 / 8;
			SDeriv(17, 2) = 1.0 / 4.0 - t2 / 4 + xi / 4 - t7 / 4;
			SDeriv(18, 2) = zeta / 4 - 1.0 / 8.0 - t1 / 4 + t2 / 8 + t3 / 4 + t4 / 8 - t6 / 4 - t5 / 8 + t7 / 8 - t8 / 8;
			SDeriv(19, 2) = -zeta / 2 + t1 / 2 - t3 / 2 + t6 / 2;
		}
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeTet4(const Coordinate & p, ElementVector & S)
	{
		const double xi = p[0];
		const double eta = p[1];
		const double zeta = p[2];

		S(0) = 1 - xi - eta - zeta;
		S(1) = xi;
		S(2) = eta;
		S(3) = zeta;
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivTet4(const Coordinate & p, ElementVectorDeriv & SDeriv)
	{
		// Derivatives w.r.t. xi (column 0)
		{
			SDeriv(0, 0) = -1;
			SDeriv(1, 0) = 1;
			SDeriv(2, 0) = 0;
			SDeriv(3, 0) = 0;
		}
		// Derivatives w.r.t. eta (column 1)
		{
			SDeriv(0, 0) = -1;
			SDeriv(1, 0) = 0;
			SDeriv(2, 0) = 1;
			SDeriv(3, 0) = 0;
		}
		// Derivatives w.r.t. zeta (column 2)
		{
			SDeriv(0, 0) = -1;
			SDeriv(1, 0) = 0;
			SDeriv(2, 0) = 0;
			SDeriv(3, 0) = 1;
		}
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeTet10(const Coordinate & p, ElementVector & S)
	{
		const double xi = p[0];
		const double eta = p[1];
		const double zeta = p[2];
		{
			const double t4 = eta * zeta;
			const double t5 = 0.4e1 * t4;
			const double t6 = zeta * xi;
			const double t7 = 0.4e1 * t6;
			const double t8 = xi * eta;
			const double t9 = 0.4e1 * t8;
			const double t10 = xi * xi;
			const double t11 = 0.2e1 * t10;
			const double t12 = eta * eta;
			const double t13 = 0.2e1 * t12;
			const double t14 = zeta * zeta;
			const double t15 = 0.2e1 * t14;
			S(0) = 0.1e1 - 0.3e1 * xi - 0.3e1 * eta - 0.3e1 * zeta + t5 + t7 + t9 + t11 + t13 + t15;
			S(1) = -xi + t11;
			S(2) = -eta + t13;
			S(3) = -zeta + t15;
			S(4) = 0.4e1 * xi - 0.4e1 * t6 - 0.4e1 * t8 - 0.4e1 * t10;
			S(5) = t9;
			S(6) = 0.4e1 * eta - 0.4e1 * t4 - 0.4e1 * t8 - 0.4e1 * t12;
			S(7) = 0.4e1 * zeta - 0.4e1 * t4 - 0.4e1 * t6 - 0.4e1 * t14;
			S(8) = t7;
			S(9) = t5;
		}
	}


	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivTet10(const Coordinate & p, ElementVectorDeriv & SDeriv)
	{
		const double xi = p[0];
		const double eta = p[1];
		const double zeta = p[2];
		// Derivatives w.r.t. xi (column 0)
		{
			const double t1 = 0.4e1 * zeta;
			const double t2 = 0.4e1 * eta;
			const double t3 = 0.4e1 * xi;
			SDeriv(0, 0) = -0.3e1 + t1 + t2 + t3;
			SDeriv(1, 0) = -0.1e1 + t3;
			SDeriv(2, 0) = 0.0e0;
			SDeriv(3, 0) = 0.0e0;
			SDeriv(4, 0) = 0.4e1 - t1 - t2 - 0.8e1 * xi;
			SDeriv(5, 0) = t2;
			SDeriv(6, 0) = -t2;
			SDeriv(7, 0) = -t1;
			SDeriv(8, 0) = t1;
			SDeriv(9, 0) = 0.0e0;
		}
		// Derivatives w.r.t. eta (column 1)
		{
			const double t1 = 0.4e1 * zeta;
			const double t2 = 0.4e1 * eta;
			const double t3 = 0.4e1 * xi;
			SDeriv(0, 1) = -0.3e1 + t1 + t2 + t3;
			SDeriv(1, 1) = 0.0e0;
			SDeriv(2, 1) = -0.1e1 + t2;
			SDeriv(3, 1) = 0.0e0;
			SDeriv(4, 1) = -t3;
			SDeriv(5, 1) = t3;
			SDeriv(6, 1) = 0.4e1 - t1 - t3 - 0.8e1 * eta;
			SDeriv(7, 1) = -t1;
			SDeriv(8, 1) = 0.0e0;
			SDeriv(9, 1) = t1;
		}
		// Derivatives w.r.t. zeta (column 2)
		{
			const double t1 = 0.4e1 * zeta;
			const double t2 = 0.4e1 * eta;
			const double t3 = 0.4e1 * xi;
			SDeriv(0, 2) = -0.3e1 + t1 + t2 + t3;
			SDeriv(1, 2) = 0.0e0;
			SDeriv(2, 2) = 0.0e0;
			SDeriv(3, 2) = -0.1e1 + t1;
			SDeriv(4, 2) = -t3;
			SDeriv(5, 2) = 0.0e0;
			SDeriv(6, 2) = -t2;
			SDeriv(7, 2) = 0.4e1 - t2 - t3 - 0.8e1 * zeta;
			SDeriv(8, 2) = t3;
			SDeriv(9, 2) = t2;
		}
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeWedge15(const Coordinate & p, ElementVector & S)
	{
		const double xi = p[0];
		const double eta = p[1];
		const double zeta = p[2];
		{
			const double t1 = zeta / 0.2e1;
			const double t2 = eta * zeta;
			const double t3 = 0.3e1 / 0.2e1 * t2;
			const double t4 = zeta * xi;
			const double t5 = 0.3e1 / 0.2e1 * t4;
			const double t6 = xi * eta;
			const double t7 = 0.2e1 * t6;
			const double t8 = xi * xi;
			const double t9 = eta * eta;
			const double t10 = zeta * zeta;
			const double t11 = t10 / 0.2e1;
			const double t12 = t8 * zeta;
			const double t13 = t9 * zeta;
			const double t14 = t10 * xi;
			const double t15 = t14 / 0.2e1;
			const double t16 = t10 * eta;
			const double t17 = t16 / 0.2e1;
			const double t18 = t6 * zeta;
			const double t19 = 0.2e1 * t18;
			const double t20 = -xi - eta - t1 + t3 + t5 + t7 + t8 + t9 + t11 - t12 - t13 - t15 - t17 - t19;
			const double t22 = t4 / 0.2e1;
			const double t25 = t2 / 0.2e1;
			const double t31 = -xi - eta + t1 - t3 - t5 + t7 + t8 + t9 + t11 + t12 + t13 - t15 - t17 + t19;
			S(0) = t20;
			S(1) = 0.2e1 * xi - 0.2e1 * t4 - 0.2e1 * t6 - 0.2e1 * t8 + 0.2e1 * t12 + 0.2e1 * t18;
			S(2) = -xi + t22 + t8 - t12 + t15;
			S(3) = 0.2e1 * t6 - 0.2e1 * t18;
			S(4) = -eta + t25 + t9 - t13 + t17;
			S(5) = 0.2e1 * eta - 0.2e1 * t2 - 0.2e1 * t6 - 0.2e1 * t9 + 0.2e1 * t13 + 0.2e1 * t18;
			S(6) = 0.1e1 - xi - eta - t10 + t14 + t16;
			S(7) = xi - t14;
			S(8) = eta - t16;
			S(9) = t31;
			S(10) = 0.2e1 * xi + 0.2e1 * t4 - 0.2e1 * t6 - 0.2e1 * t8 - 0.2e1 * t12 - 0.2e1 * t18;
			S(11) = -xi - t22 + t8 + t12 + t15;
			S(12) = 0.2e1 * t6 + 0.2e1 * t18;
			S(13) = -eta - t25 + t9 + t13 + t17;
			S(14) = 0.2e1 * eta + 0.2e1 * t2 - 0.2e1 * t6 - 0.2e1 * t9 - 0.2e1 * t13 - 0.2e1 * t18;
		}
	}

	template<typename Coordinate, typename ElementVector, typename ElementVectorDeriv>
	inline void ShapeFunctions<Coordinate, ElementVector, ElementVectorDeriv>::shapeDerivWedge15(const Coordinate & p, ElementVectorDeriv & SDeriv)
	{
		const double xi = p[0];
		const double eta = p[1];
		const double zeta = p[2];
		// Derivatives w.r.t. xi (column 0)
		{
			const double t1 = 0.3e1 / 0.2e1 * zeta;
			const double t2 = 0.2e1 * eta;
			const double t3 = 0.2e1 * xi;
			const double t4 = zeta * xi;
			const double t5 = 0.2e1 * t4;
			const double t6 = zeta * zeta;
			const double t7 = t6 / 0.2e1;
			const double t8 = eta * zeta;
			const double t9 = 0.2e1 * t8;
			const double t11 = 0.2e1 * zeta;
			const double t12 = 0.4e1 * xi;
			const double t13 = 0.4e1 * t4;
			const double t15 = zeta / 0.2e1;
			const double t17 = eta - t8;
			const double t18 = -0.1e1 + t6;
			const double t22 = eta + t8;
			SDeriv(0, 0) = -0.1e1 + t1 + t2 + t3 - t5 - t7 - t9;
			SDeriv(1, 0) = 0.2e1 - t11 - t2 - t12 + t13 + t9;
			SDeriv(2, 0) = -0.1e1 + t15 + t3 - t5 + t7;
			SDeriv(3, 0) = 0.2e1 * t17;
			SDeriv(4, 0) = 0.0e0;
			SDeriv(5, 0) = -0.2e1 * t17;
			SDeriv(6, 0) = t18;
			SDeriv(7, 0) = -t18;
			SDeriv(8, 0) = 0.0e0;
			SDeriv(9, 0) = -0.1e1 - t1 + t2 + t3 + t5 - t7 + t9;
			SDeriv(10, 0) = 0.2e1 + t11 - t2 - t12 - t13 - t9;
			SDeriv(11, 0) = -0.1e1 - t15 + t3 + t5 + t7;
			SDeriv(12, 0) = 0.2e1 * t22;
			SDeriv(13, 0) = 0.0e0;
			SDeriv(14, 0) = -0.2e1 * t22;
		}
		// Derivatives w.r.t. eta (column 1)
		{
			const double t1 = 0.3e1 / 0.2e1 * zeta;
			const double t2 = 0.2e1 * eta;
			const double t3 = 0.2e1 * xi;
			const double t4 = zeta * xi;
			const double t5 = 0.2e1 * t4;
			const double t6 = zeta * zeta;
			const double t7 = t6 / 0.2e1;
			const double t8 = eta * zeta;
			const double t9 = 0.2e1 * t8;
			const double t11 = -xi + t4;
			const double t12 = zeta / 0.2e1;
			const double t14 = 0.2e1 * zeta;
			const double t15 = 0.4e1 * eta;
			const double t16 = 0.4e1 * t8;
			const double t18 = -0.1e1 + t6;
			const double t20 = -xi - t4;
			SDeriv(0, 1) = -0.1e1 + t1 + t2 + t3 - t5 - t7 - t9;
			SDeriv(1, 1) = 0.2e1 * t11;
			SDeriv(2, 1) = 0.0e0;
			SDeriv(3, 1) = -0.2e1 * t11;
			SDeriv(4, 1) = -0.1e1 + t12 + t2 - t9 + t7;
			SDeriv(5, 1) = 0.2e1 - t14 - t3 - t15 + t16 + t5;
			SDeriv(6, 1) = t18;
			SDeriv(7, 1) = 0.0e0;
			SDeriv(8, 1) = -t18;
			SDeriv(9, 1) = -0.1e1 - t1 + t2 + t3 + t5 - t7 + t9;
			SDeriv(10, 1) = 0.2e1 * t20;
			SDeriv(11, 1) = 0.0e0;
			SDeriv(12, 1) = -0.2e1 * t20;
			SDeriv(13, 1) = -0.1e1 - t12 + t2 + t9 + t7;
			SDeriv(14, 1) = 0.2e1 + t14 - t3 - t15 - t16 - t5;
		}
		// Derivatives w.r.t. zeta (column 2)
		{
			const double t1 = 0.3e1 / 0.2e1 * eta;
			const double t2 = 0.3e1 / 0.2e1 * xi;
			const double t3 = xi * xi;
			const double t4 = eta * eta;
			const double t5 = zeta * xi;
			const double t6 = eta * zeta;
			const double t7 = xi * eta;
			const double t8 = 0.2e1 * t7;
			const double t10 = -xi + t3 + t7;
			const double t11 = xi / 0.2e1;
			const double t13 = eta / 0.2e1;
			const double t15 = -eta + t4 + t7;
			SDeriv(0, 2) = -0.1e1 / 0.2e1 + t1 + t2 + zeta - t3 - t4 - t5 - t6 - t8;
			SDeriv(1, 2) = 0.2e1 * t10;
			SDeriv(2, 2) = t11 - t3 + t5;
			SDeriv(3, 2) = -t8;
			SDeriv(4, 2) = t13 - t4 + t6;
			SDeriv(5, 2) = 0.2e1 * t15;
			SDeriv(6, 2) = -0.2e1 * zeta + 0.2e1 * t5 + 0.2e1 * t6;
			SDeriv(7, 2) = -0.2e1 * t5;
			SDeriv(8, 2) = -0.2e1 * t6;
			SDeriv(9, 2) = 0.1e1 / 0.2e1 - t1 - t2 + zeta + t3 + t4 - t5 - t6 + t8;
			SDeriv(10, 2) = -0.2e1 * t10;
			SDeriv(11, 2) = -t11 + t3 + t5;
			SDeriv(12, 2) = t8;
			SDeriv(13, 2) = -t13 + t4 + t6;
			SDeriv(14, 2) = -0.2e1 * t15;
		}
	}

#ifdef __MSC_VER
#pragma endregion
#endif 

}
