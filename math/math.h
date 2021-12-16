/*
  ==============================================================================

    math.h
    Created: 28 Mar 2021 6:06:23pm
    Author:  Julian Vanasse

        Improved vector operation templates.
 
    This namespace exists to make C++ programming as close as possible to the
    ease of MATLAB or python programming (don't laugh).

    Requires Boost and kiss_fft libraries.
 
 
  ==============================================================================
*/

#pragma once

#include <string>
#include <complex>
#include <iostream>
#include <algorithm>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/assignment.hpp>

#include "../kiss_fft/kiss_fft.h"
#include "../kiss_fft/_kiss_fft_guts.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace eers {
    namespace math {
        template <typename T> boost::numeric::ublas::vector<T> add(boost::numeric::ublas::vector<T> x, T a);

        template <typename T> boost::numeric::ublas::vector<T> max(boost::numeric::ublas::vector<T> x, T a);
        template <typename T> boost::numeric::ublas::vector<T> max(boost::numeric::ublas::vector<T> x, boost::numeric::ublas::vector<T> y);
        template <typename T> boost::numeric::ublas::vector<T> min(boost::numeric::ublas::vector<T> x, T a);
        template <typename T> boost::numeric::ublas::vector<T> min(boost::numeric::ublas::vector<T> x, boost::numeric::ublas::vector<T> y);

        template <typename T> T mean(boost::numeric::ublas::vector<T> x);
        template <typename T> T sum(boost::numeric::ublas::vector<T> x);

        template <typename T> boost::numeric::ublas::vector<T> real(boost::numeric::ublas::vector<std::complex<T> > z);
        template <typename T> boost::numeric::ublas::vector<T> imag(boost::numeric::ublas::vector<std::complex<T> > z);
        template <typename T> boost::numeric::ublas::vector<T> abs(boost::numeric::ublas::vector<std::complex<T> > z);
        template <typename T> boost::numeric::ublas::vector<T> abs(boost::numeric::ublas::vector<T> z);
        template <typename T> boost::numeric::ublas::vector<T> angle(boost::numeric::ublas::vector<std::complex<T> > z);

        template <typename T> boost::numeric::ublas::vector<T> pow(boost::numeric::ublas::vector<T> x, T p);
        template <typename T> boost::numeric::ublas::vector<T> pow(T b, boost::numeric::ublas::vector<T> p);
        template <typename T> boost::numeric::ublas::vector<T> log10(boost::numeric::ublas::vector<T> x);
        template <typename T> boost::numeric::ublas::vector<T> logb(boost::numeric::ublas::vector<T> x, T b);
        template <typename T> T dB2lin(T x);
        template <typename T> boost::numeric::ublas::vector<T> dB2lin(boost::numeric::ublas::vector<T> x);

        template <typename T> boost::numeric::ublas::vector<T> sin(boost::numeric::ublas::vector<T> t);
        template <typename T> boost::numeric::ublas::vector<T> cos(boost::numeric::ublas::vector<T> t);
        template <typename T> boost::numeric::ublas::vector<T> wrap_to_pi(boost::numeric::ublas::vector<T> t);

        template <typename T> boost::numeric::ublas::vector<T> zp(boost::numeric::ublas::vector<T> x, int len);
        template <typename T> boost::numeric::ublas::vector<std::complex<T> > zp(boost::numeric::ublas::vector<std::complex<T> > z, int len);

        template <typename T> boost::numeric::ublas::vector<std::complex<T> > complex(boost::numeric::ublas::vector<T> r);
        template <typename T> boost::numeric::ublas::vector<std::complex<T> > complex(boost::numeric::ublas::vector<T> r, boost::numeric::ublas::vector<T> i);
        template <typename T> boost::numeric::ublas::vector<std::complex<T> > pol2cart(boost::numeric::ublas::vector<T> r, boost::numeric::ublas::vector<T> i);

        template <typename T> boost::numeric::ublas::vector<std::complex<T> > fft(boost::numeric::ublas::vector<T> r, kiss_fft_cfg dir);
        template <typename T> boost::numeric::ublas::vector<std::complex<T> > fft(boost::numeric::ublas::vector<std::complex<T> > z, kiss_fft_cfg dir);

        template <typename T> boost::numeric::ublas::vector<T> poly(boost::numeric::ublas::vector<T> roots);
        template <typename T> T polyval(boost::numeric::ublas::vector<T> p, T x);

        template <typename T> bool has_nan(boost::numeric::ublas::vector<T> x);
        template <typename T> bool has_inf(boost::numeric::ublas::vector<T> x);

        template <typename T> boost::numeric::ublas::vector<T> linspace(T start, T end, int num);
        template <typename T> boost::numeric::ublas::vector<T> logspace(T start, T end, int num);

        template <typename T> T rescale(T value, T source_min, T source_max, T dest_min, T dest_max);
    };
};
