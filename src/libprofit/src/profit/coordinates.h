/**
 * Coordinate classes definition
 *
 * ICRAR - International Centre for Radio Astronomy Research
 * (c) UWA - The University of Western Australia, 2021
 * Copyright by UWA (in the framework of the ICRAR)
 * All rights reserved
 *
 * Contributed by Rodrigo Tobar
 *
 * This file is part of libprofit.
 *
 * libprofit is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libprofit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libprofit.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PROFIT_COORDINATES_H
#define PROFIT_COORDINATES_H

#include <cassert>
#include <ostream>

#include "profit/common.h"

namespace profit {

/** An (x, y) pair in a 2-dimensional surface.
 *
 * Comparison between these objects can be done with the <, <=, ==, !=, > and >=
 * operators, but users should note that there is no way to order values based
 * on these operators (that is, objects of this type are by themselves
 * non-sortable).
 */
template <typename T, typename Derived>
class _2d_coordinate {

public:

	_2d_coordinate() = default;
	~_2d_coordinate() = default;
	_2d_coordinate(T x, T y) : x(x), y(y) {}

	T x {0};
	T y {0};

	/// greater or equal comparison across both dimensions
	bool operator>=(const _2d_coordinate &other) const
	{
		return x >= other.x && y >= other.y;
	}

	/// greater than comparison across both dimensions
	bool operator>(const _2d_coordinate &other) const
	{
		return x > other.x && y > other.y;
	}

	/// less or equal comparison across both dimensions
	bool operator<=(const _2d_coordinate &other) const
	{
		return x <= other.x && y <= other.y;
	}

	/// less than comparison across both dimensions
	bool operator<(const _2d_coordinate &other) const
	{
		return x < other.x && y < other.y;
	}

	bool operator==(const _2d_coordinate &other) const {
		return x == other.x && y == other.y;
	}

	bool operator!=(const _2d_coordinate &other) const {
		return x != other.x || y != other.y;
	}

	_2d_coordinate &operator+=(const _2d_coordinate &other)
	{
		x += other.x;
		y += other.y;
		return *this;
	}

	Derived operator+(const _2d_coordinate &other) const
	{
		Derived sum(*this);
		sum += other;
		return sum;
	}

	_2d_coordinate &operator+=(T n)
	{
		x += n;
		y += n;
		return *this;
	}

	Derived operator+(T n) const
	{
		Derived sum(*this);
		sum += n;
		return sum;
	}

	_2d_coordinate &operator-=(const _2d_coordinate &other)
	{
		x -= other.x;
		y -= other.y;
		return *this;
	}

	Derived operator-(const _2d_coordinate &other) const
	{
		Derived sum(*this);
		sum -= other;
		return sum;
	}

	_2d_coordinate &operator*=(T f)
	{
		x *= f;
		y *= f;
		return *this;
	}

	Derived operator*(T f) const
	{
		Derived mul(*this);
		mul *= f;
		return mul;
	}

	_2d_coordinate &operator/=(T f)
	{
		x /= f;
		y /= f;
		return *this;
	}

	Derived operator/(T f) const
	{
		Derived div(*this);
		div /= f;
		return div;
	}
};

/// Stream output operator
template <typename CharT, typename T, typename Derived>
std::basic_ostream<CharT> &operator<<(std::basic_ostream<CharT> &os, const _2d_coordinate<T, Derived> &coord)
{
	os << '[' << coord.x << ", " << coord.y << ']';
	return os;
}

/// Allow subtraction of different numeric scalars
template <typename T, typename Derived>
inline
Derived operator-(T x, const _2d_coordinate<T, Derived> &other)
{
	return Derived(x, x) - other;
}

template <typename T, typename Derived>
inline Derived operator-(const _2d_coordinate<T, Derived> &other, T x)
{
	return other - Derived(x, x);
}


/// Element-wise max() function for _2d_coordinate objects
template <typename T, typename Derived>
inline Derived max(const _2d_coordinate<T, Derived> a, const _2d_coordinate<T, Derived> b)
{
	return Derived{std::max(a.x, b.x), std::max(a.x, b.x)};
}

/// A 2-d coordinate in a discrete surface starting at (0, 0)
class discrete_2d_coordinate
  : public _2d_coordinate<unsigned int, discrete_2d_coordinate>
{
public:
	using _2d_coordinate::_2d_coordinate;

	template <typename T, typename Derived>
	discrete_2d_coordinate(const _2d_coordinate<T, Derived> &other)
	  : _2d_coordinate(other)
	{}

	discrete_2d_coordinate(int x, int y)
	{
		assert(x >= 0);
		assert(y >= 0);
		this->x = x;
		this->y = y;
	}

	explicit operator bool() const
	{
		return x > 0 && y > 0;
	}

	discrete_2d_coordinate &operator%=(unsigned int i)
	{
		x %= i;
		y %= i;
		return *this;
	}

	discrete_2d_coordinate operator%(unsigned int i) const
	{
		discrete_2d_coordinate mod(*this);
		mod %= i;
		return mod;
	}
};


/// Allow subtraction of different numeric scalars
inline discrete_2d_coordinate operator-(int x, const discrete_2d_coordinate &other)
{
	return discrete_2d_coordinate(x, x) - other;
}

inline discrete_2d_coordinate operator-(const discrete_2d_coordinate &other, int x)
{
	return other - discrete_2d_coordinate(x, x);
}

/// A 2-d coordinate in a continuous surface
template <typename Derived>
class continuous_2d_coordinate
  : public _2d_coordinate<double, Derived>
{
public:
	using Base = _2d_coordinate<double, Derived>;

	using Base::Base;
	continuous_2d_coordinate(int x, int y) : Base(x, y) {}

	template <typename U, typename Derived2>
	continuous_2d_coordinate(const _2d_coordinate<U, Derived2> &other)
	  : Base(other)
	{}
};

}  // namespace profit

#endif /* PROFIT_COORDINATES_H */