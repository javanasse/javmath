/*
  ==============================================================================

    VectorOperations2.cpp
    Created: 28 Mar 2021 6:06:23pm
    Author:  Julian Vanasse

  ==============================================================================
*/

#include "math.h"

template <typename T>
boost::numeric::ublas::vector<T> eers::math::add(boost::numeric::ublas::vector<T> x, T a)
{
    boost::numeric::ublas::vector<T> y (x.size());
    for (int n = 0; n < x.size(); n++)
    {
        y(n) = x(n) + a;
    }
    return y;
}

template <typename T>
boost::numeric::ublas::vector<T> eers::math::max(boost::numeric::ublas::vector<T> x, T a)
{
    boost::numeric::ublas::vector<T> y (x.size());
    for (int n = 0; n < x.size(); n++)
    {
        if (x(n) > a)
            y(n) = x(n);
        else
            x(n) = a;
    }
    return y;
}

template <typename T>
boost::numeric::ublas::vector<T> eers::math::max(boost::numeric::ublas::vector<T> x, boost::numeric::ublas::vector<T> y)
{
    boost::numeric::ublas::vector<T> w (x.size());
    for (int n = 0; n < x.size(); n++)
    {
        if (x(n) > y(n))
            w(n) = x(n);
        else
            w(n) = y(n);
    }
    return y;
}

template <typename T>
boost::numeric::ublas::vector<T> eers::math::min(boost::numeric::ublas::vector<T> x, T a)
{
    boost::numeric::ublas::vector<T> y (x.size());
    for (int n = 0; n < x.size(); n++)
    {
        if (x(n) < a)
            y(n) = x(n);
        else
            x(n) = a;
    }
    return y;
}

template <typename T>
boost::numeric::ublas::vector<T> eers::math::min(boost::numeric::ublas::vector<T> x, boost::numeric::ublas::vector<T> y)
{
    boost::numeric::ublas::vector<T> w (x.size());
    for (int n = 0; n < x.size(); n++)
    {
        if (x(n) < y(n))
            w(n) = x(n);
        else
            w(n) = y(n);
    }
    return y;
}

template <typename T>
T eers::math::mean(boost::numeric::ublas::vector<T> x)
{
    T mu = T();
    for (int n = 0; n < x.size(); n++)
    {
        mu += x(n);
    }
    mu /= static_cast<T>(x.size());
    return mu;
}

template <typename T>
T eers::math::sum(boost::numeric::ublas::vector<T> x)
{
    T s = T();
    for (int n = 0; n < x.size(); n++)
    {
        s += x(n);
    }
    return s;
}

template <typename T>
boost::numeric::ublas::vector<T> eers::math::real(boost::numeric::ublas::vector<std::complex<T> > z)
{
    boost::numeric::ublas::vector<T> r (z.size());
    for (int n = 0; n < z.size(); n++)
    {
        r(n) = z(n).real();
    }
    return r;
}

template <typename T>
boost::numeric::ublas::vector<T> eers::math::imag(boost::numeric::ublas::vector<std::complex<T> > z)
{
    boost::numeric::ublas::vector<T> i (z.size());
    for (int n = 0; n < z.size(); n++)
    {
        i(n) = z(n).imag();
    }
    return i;
}

template <typename T>
boost::numeric::ublas::vector<T> eers::math::abs(boost::numeric::ublas::vector<std::complex<T> > z)
{
    boost::numeric::ublas::vector<T> v (z.size());
    for (int n = 0; n < z.size(); n++)
    {
        v(n) = abs(z(n));
    }
    return v;
}

template <typename T>
boost::numeric::ublas::vector<T> eers::math::abs(boost::numeric::ublas::vector<T> x)
{
    boost::numeric::ublas::vector<T> y (x.size());
    for (int n = 0; n < x.size(); n++)
    {
        y(n) = std::abs(x(n));
    }
    return x;
}

template <typename T>
boost::numeric::ublas::vector<T> eers::math::angle(boost::numeric::ublas::vector<std::complex<T> > z)
{
    boost::numeric::ublas::vector<T> v (z.size());
    for (int n = 0; n < z.size(); n++)
    {
        v(n) = std::arg(z(n));
    }
    return v;
}

template <typename T>
boost::numeric::ublas::vector<T> eers::math::pow(boost::numeric::ublas::vector<T> x, T p)
{
    boost::numeric::ublas::vector<T> y (x.size());
    for (int n = 0; n < x.size(); n++)
    {
        y(n) = std::pow(x(n), p);
    }
    return y;
}

template<typename T>
boost::numeric::ublas::vector<T> eers::math::pow(T b, boost::numeric::ublas::vector<T> p)
{
    boost::numeric::ublas::vector<T> y(p.size(), T());
    for (int n = 0; n < p.size(); n++)
    {
        y(n) = std::pow(b, p(n));
    }

    return y;
}

template <typename T>
boost::numeric::ublas::vector<T> eers::math::log10(boost::numeric::ublas::vector<T> x)
{
    boost::numeric::ublas::vector<T> y (x.size());
    for (int n = 0; n < x.size(); n++)
    {
        y(n) = std::log10(x(n));
    }
    return y;
}

template<typename T>
boost::numeric::ublas::vector<T> eers::math::logb(boost::numeric::ublas::vector<T> x, T b)
{
    boost::numeric::ublas::vector<T> y(x.size());
    for (int n = 0; n < x.size(); n++)
    {
        y(n) = std::log(x(n)) / std::log(b);
    }
    return y;
}

template<typename T>
T eers::math::dB2lin(T x)
{
    return std::pow(T(10.0), x / T(20.0));
}

template<typename T>
boost::numeric::ublas::vector<T> eers::math::dB2lin(boost::numeric::ublas::vector<T> x)
{
    boost::numeric::ublas::vector<T> y (x.size(), T(0.0));
    for (int n = 0; n < y.size(); n++)
    {
        y(n) = eers::math::dB2lin(x(n));
    }
    return y;
}

template <typename T> boost::numeric::ublas::vector<T> eers::math::sin(boost::numeric::ublas::vector<T> t)
{
    boost::numeric::ublas::vector<T> y (t.size());
    for (int n = 0; n < t.size(); n++)
    {
        y(n) = std::sin(t(n));
    }
    return y;
}

template <typename T> boost::numeric::ublas::vector<T> eers::math::cos(boost::numeric::ublas::vector<T> t)
{
    boost::numeric::ublas::vector<T> y (t.size());
    for (int n = 0; n < t.size(); n++)
    {
        y(n) = std::cos(t(n));
    }
    return y;
}

template <typename T> boost::numeric::ublas::vector<T> eers::math::wrap_to_pi(boost::numeric::ublas::vector<T> t)
{
    boost::numeric::ublas::vector<T> y (t.size());
    for (int n = 0; n < t.size(); n++)
    {
        y(n) = t(n) - (2.0 * M_PI * std::floor((t(n) + M_PI) / (2.0 * M_PI)));
    }
    return y;
}

template <typename T>
boost::numeric::ublas::vector<T> eers::math::zp(boost::numeric::ublas::vector<T> x, int len)
{
    /* zero-pad */
    boost::numeric::ublas::vector<T> y (x.size() + len, 0.0);
    for (int n = 0; n < x.size(); n++)
    {
        y(n) = x(n);
    }
    return y;
}

template <typename T> boost::numeric::ublas::vector<std::complex<T> > eers::math::zp(boost::numeric::ublas::vector<std::complex<T> > z, int len)
{
    /* zero-pad */
    boost::numeric::ublas::vector<std::complex<T> > y (z.size() + len, 0.0);
    for (int n = 0; n < z.size(); n++)
    {
        y(n) = z(n);
    }
    return y;
}

template <typename T>
boost::numeric::ublas::vector<std::complex<T> > eers::math::complex(boost::numeric::ublas::vector<T> r)
{
    /* takes real vector and outputs complex vector */
    boost::numeric::ublas::vector<std::complex<T> > z (r.size());
    for (int n = 0; n < z.size(); n++)
    {
        z(n) = r(n);
    }
    return z;
}

template <typename T>
boost::numeric::ublas::vector<std::complex<T> > eers::math::complex(boost::numeric::ublas::vector<T> r, boost::numeric::ublas::vector<T> i)
{
    /* takes real vector and outputs complex vector */
    boost::numeric::ublas::vector<std::complex<T> > z (r.size());
    std::complex<T> u;
    for (int n = 0; n < z.size(); n++)
    {
        u = std::complex<T>(r(n), i(n));
        z(n) = u;
    }
    return z;
}

template <typename T>
boost::numeric::ublas::vector<std::complex<T> > eers::math::pol2cart(boost::numeric::ublas::vector<T> r, boost::numeric::ublas::vector<T> p)
{
    /* convert magnitude and phase vectors to cartesian complex vector */
    
    if (r.size() != p.size())
        throw std::length_error("r and p must be same size");
    
    boost::numeric::ublas::vector<std::complex<T> > z (r.size(), T());
    for (int n = 0; n < z.size(); n++)
    {
        T re = r(n) * std::cos(p(n));
        T im = r(n) * std::sin(p(n));
        z(n) = std::complex<T> (re, im);
    }
    return z;
}

template <typename T>
boost::numeric::ublas::vector<std::complex<T> > eers::math::fft(boost::numeric::ublas::vector<T> r, kiss_fft_cfg dir)
{
    /* wrapper for kiss_fft */
    
    if (dir->nfft < r.size())
        throw std::length_error("vector size must be <= n_fft");
    
    // result
    boost::numeric::ublas::vector<std::complex<T> > result (dir->nfft);
    // use kiss_fft_cpx complex type
    boost::numeric::ublas::vector<kiss_fft_cpx> buffer (dir->nfft);
    // copy r into real position of buffer
    for (int n = 0; n < r.size(); n++)
    {
        buffer(n).r = r(n);
        buffer(n).i = 0.0;
    }
    for (int n = r.size(); n < buffer.size(); n++)
    {
        buffer(n).r = 0.0;
        buffer(n).i = 0.0;
    }
    
    // transform
    kiss_fft(dir, buffer.data().begin(), buffer.data().begin());
    
    // copy to result
    for (int n = 0; n < buffer.size(); n++)
    {
        result(n) = std::complex<T>(static_cast<T>(buffer(n).r), static_cast<T>(buffer(n).i));
    }
    
    // rescale if inverse
    if (dir->inverse)
    {
        result *= (1.0 / static_cast<T>(dir->nfft));
    }
    
    return result;
}

template <typename T>
boost::numeric::ublas::vector<std::complex<T> > eers::math::fft(boost::numeric::ublas::vector<std::complex<T> > z, kiss_fft_cfg dir)
{
    /* wrapper for kiss_fft */
    
    if (dir->nfft < z.size())
        throw std::length_error("vector size must be <= n_fft");
    
    // result
    boost::numeric::ublas::vector<std::complex<T> > result (dir->nfft);
    // use kiss_fft_cpx complex type
    boost::numeric::ublas::vector<kiss_fft_cpx> buffer (dir->nfft);
    // copy r into real position of buffer
    for (int n = 0; n < z.size(); n++)
    {
        buffer(n).r = z(n).real();
        buffer(n).i = z(n).imag();
    }
    for (int n = z.size(); n < buffer.size(); n++)
    {
        buffer(n).r = 0.0;
        buffer(n).i = 0.0;
    }
    
    // transform
    kiss_fft(dir, buffer.data().begin(), buffer.data().begin());
    
    // copy to result
    for (int n = 0; n < buffer.size(); n++)
    {
        result(n) = std::complex<T>(static_cast<T>(buffer(n).r), static_cast<T>(buffer(n).i));
    }
    
    // rescale if inverse
    if (dir->inverse)
    {
        result *= (1.0 / static_cast<T>(dir->nfft));
    }
    
    return result;
}

///
/// Expands roots of polynomial to the polynomial coefficients.
///
template <typename T> boost::numeric::ublas::vector<T> eers::math::poly(boost::numeric::ublas::vector<T> roots)
{
    std::vector<T> result;
    result.push_back((T)1.0);
    for (size_t i = 0; i < roots.size(); ++i)
    {
        std::vector<T> temp;
        temp.resize(result.size());
        std::copy(result.begin(), result.end(), temp.begin());

        for (auto& item : temp)
            item *= (-roots[i]);

        result.push_back(0);
        for (size_t j = 1; j < result.size(); ++j)
        {
            result[j] += temp[j - 1];
        }
    }
    boost::numeric::ublas::vector<T> coefficients(result.size());
    std::copy(result.begin(), result.end(), coefficients.data().begin());

    return coefficients;
}

///
/// Evaluates polynomial p at value x.
/// Polynomial coefficients is stored in descending order as follows:
/// p(0) is the highest order coefficient, p(end) is the constant term.
///
template <typename T> T eers::math::polyval(boost::numeric::ublas::vector<T> p, T x)
{
    T acc = T();
    //for (int n = p.size()-1; n >= 0; n--)
    for (int n = 0; n < p.size(); n++)
    {
        acc = (acc * x) + p(n);
    }
    return acc;
}

template <typename T> bool eers::math::has_nan(boost::numeric::ublas::vector<T> x)
{
    /* Return true if x contains nan value(s) */
    
    for (int n = 0; n < x.size(); n++)
    {
        if (std::isnan(x(n)))
            return true;
    }
    return false;
}

template <typename T> bool eers::math::has_inf(boost::numeric::ublas::vector<T> x)
{
    for (int n = 0; n < x.size(); n++)
    {
        if (std::isinf(x(n)))
            return true;
    }
    return false;
}

template<typename T>
boost::numeric::ublas::vector<T> eers::math::linspace(T start, T end, int num)
{
    boost::numeric::ublas::vector<T> result(num, T());
    if (num == 0)
    {
        return result;
    }
    if (num == 1)
    {
        result(0) = start;
        return result;
    }
    if (num < 0)
    {
        throw std::invalid_argument("num must be > 0");
    }

    T delta = (end - start) / static_cast<T>(num);
    for (int n = 0; n < num; n++)
    {
        T tn = static_cast<T>(n);
        result(n) = start + (tn * delta);
    }

    return result;
}


/*
    eers::math::logspace(T start, T end, int num)

    returns a vector of logarithmically spaced points from START to END with
    a base of 10.

    START and END must be real and greater than zero.
*/
template<typename T>
boost::numeric::ublas::vector<T> eers::math::logspace(T start, T end, int num)
{
    if (start < 0 || end < 0)
    {
        throw std::invalid_argument("start and end must be positive.");
    }

    // take log of bounds
    start = std::log10(start);
    end = std::log10(end);

    // space points 
    boost::numeric::ublas::vector<T> line = eers::math::linspace(start, end, num);
    
    // take power 10 and return
    return eers::math::pow(T(10.0), line);
}

/*

    rescale a value in one interval to a different interval. 

    This is a shift and scale transformation, useful for graphics.

*/
template<typename T>
T eers::math::rescale(T value, T source_min, T source_max, T dest_min, T dest_max)
{
    return (((value - source_min) * (dest_max - dest_min)) / (source_max - source_min)) + dest_min;
}


template boost::numeric::ublas::vector<float>                   eers::math::add(boost::numeric::ublas::vector<float> x, float a);
template boost::numeric::ublas::vector<std::complex<float> >    eers::math::add(boost::numeric::ublas::vector<std::complex<float> > x, std::complex<float>  a);
template boost::numeric::ublas::vector<float>                   eers::math::max(boost::numeric::ublas::vector<float> x, float a);
template boost::numeric::ublas::vector<float>                   eers::math::max(boost::numeric::ublas::vector<float> x, boost::numeric::ublas::vector<float> y);
template boost::numeric::ublas::vector<float>                   eers::math::min(boost::numeric::ublas::vector<float> x, float a);
template boost::numeric::ublas::vector<float>                   eers::math::min(boost::numeric::ublas::vector<float> x, boost::numeric::ublas::vector<float> y);
template float                                                  eers::math::mean(boost::numeric::ublas::vector<float> x);
template float                                                  eers::math::sum(boost::numeric::ublas::vector<float> x);
template boost::numeric::ublas::vector<float>                   eers::math::real(boost::numeric::ublas::vector<std::complex<float> > z);
template boost::numeric::ublas::vector<float>                   eers::math::imag(boost::numeric::ublas::vector<std::complex<float> > z);
template boost::numeric::ublas::vector<float>                   eers::math::abs(boost::numeric::ublas::vector<std::complex<float> > z);
template boost::numeric::ublas::vector<float>                   eers::math::abs(boost::numeric::ublas::vector<float> x);
template boost::numeric::ublas::vector<float>                   eers::math::angle(boost::numeric::ublas::vector<std::complex<float> > z);
template boost::numeric::ublas::vector<float>                   eers::math::pow(boost::numeric::ublas::vector<float> x, float p);
template boost::numeric::ublas::vector<float>                   eers::math::pow(float b, boost::numeric::ublas::vector<float> p);
template boost::numeric::ublas::vector<float>                   eers::math::log10(boost::numeric::ublas::vector<float> x);
template boost::numeric::ublas::vector<float>                   eers::math::logb(boost::numeric::ublas::vector<float> x, float b);
template float                                                  eers::math::dB2lin(float x);
template boost::numeric::ublas::vector<float>                   eers::math::dB2lin(boost::numeric::ublas::vector<float> x);
template boost::numeric::ublas::vector<float>                   eers::math::sin(boost::numeric::ublas::vector<float> t);
template boost::numeric::ublas::vector<float>                   eers::math::cos(boost::numeric::ublas::vector<float> t);
template boost::numeric::ublas::vector<float>                   eers::math::wrap_to_pi(boost::numeric::ublas::vector<float> t);
template boost::numeric::ublas::vector<float>                   eers::math::zp(boost::numeric::ublas::vector<float> x, int len);
template boost::numeric::ublas::vector<std::complex<float> >    eers::math::zp(boost::numeric::ublas::vector<std::complex<float> > z, int len);
template boost::numeric::ublas::vector<std::complex<float> >    eers::math::complex(boost::numeric::ublas::vector<float> r);
template boost::numeric::ublas::vector<std::complex<float> >    eers::math::complex(boost::numeric::ublas::vector<float> r, boost::numeric::ublas::vector<float> i);
template boost::numeric::ublas::vector<std::complex<float> >    eers::math::pol2cart(boost::numeric::ublas::vector<float> r, boost::numeric::ublas::vector<float> p);
template boost::numeric::ublas::vector<std::complex<float> >    eers::math::fft(boost::numeric::ublas::vector<float> r, kiss_fft_cfg dir);
template boost::numeric::ublas::vector<std::complex<float> >    eers::math::fft(boost::numeric::ublas::vector<std::complex<float> > z, kiss_fft_cfg dir);
template boost::numeric::ublas::vector<float>                   eers::math::poly(boost::numeric::ublas::vector<float> roots);
template boost::numeric::ublas::vector<std::complex<float> >    eers::math::poly(boost::numeric::ublas::vector<std::complex<float> > roots);
template float                                                  eers::math::polyval(boost::numeric::ublas::vector<float> p, float x);
template std::complex<float>                                    eers::math::polyval(boost::numeric::ublas::vector<std::complex<float> > p, std::complex<float> x);
template bool                                                   eers::math::has_nan(boost::numeric::ublas::vector<float> x);
template bool                                                   eers::math::has_inf(boost::numeric::ublas::vector<float> x);
template boost::numeric::ublas::vector<float>                   eers::math::linspace(float start, float end, int num);
template boost::numeric::ublas::vector<float>                   eers::math::logspace(float start, float end, int num);
template float                                                  eers::math::rescale(float value, float source_min, float source_max, float dest_min, float dest_max);

template boost::numeric::ublas::vector<double>                  eers::math::add(boost::numeric::ublas::vector<double> x, double a);
template boost::numeric::ublas::vector<std::complex<double> >   eers::math::add(boost::numeric::ublas::vector<std::complex<double> > x, std::complex<double>  a);
template boost::numeric::ublas::vector<double>                  eers::math::max(boost::numeric::ublas::vector<double> x, double a);
template boost::numeric::ublas::vector<double>                  eers::math::max(boost::numeric::ublas::vector<double> x, boost::numeric::ublas::vector<double> y);
template boost::numeric::ublas::vector<double>                  eers::math::min(boost::numeric::ublas::vector<double> x, double a);
template boost::numeric::ublas::vector<double>                  eers::math::min(boost::numeric::ublas::vector<double> x, boost::numeric::ublas::vector<double> y);
template double                                                 eers::math::mean(boost::numeric::ublas::vector<double> x);
template double                                                 eers::math::sum(boost::numeric::ublas::vector<double> x);
template boost::numeric::ublas::vector<double>                  eers::math::real(boost::numeric::ublas::vector<std::complex<double> > z);
template boost::numeric::ublas::vector<double>                  eers::math::imag(boost::numeric::ublas::vector<std::complex<double> > z);
template boost::numeric::ublas::vector<double>                  eers::math::abs(boost::numeric::ublas::vector<std::complex<double> > z);
template boost::numeric::ublas::vector<double>                  eers::math::abs(boost::numeric::ublas::vector<double> x);
template boost::numeric::ublas::vector<double>                  eers::math::angle(boost::numeric::ublas::vector<std::complex<double> > z);
template boost::numeric::ublas::vector<double>                  eers::math::pow(boost::numeric::ublas::vector<double> x, double p);
template boost::numeric::ublas::vector<double>                  eers::math::pow(double b, boost::numeric::ublas::vector<double> p);
template boost::numeric::ublas::vector<double>                  eers::math::log10(boost::numeric::ublas::vector<double> x);
template boost::numeric::ublas::vector<double>                  eers::math::logb(boost::numeric::ublas::vector<double> x, double b);
template double                                                 eers::math::dB2lin(double x);
template boost::numeric::ublas::vector<double>                  eers::math::dB2lin(boost::numeric::ublas::vector<double> x);
template boost::numeric::ublas::vector<double>                  eers::math::sin(boost::numeric::ublas::vector<double> t);
template boost::numeric::ublas::vector<double>                  eers::math::cos(boost::numeric::ublas::vector<double> t);
template boost::numeric::ublas::vector<double>                  eers::math::wrap_to_pi(boost::numeric::ublas::vector<double> t);
template boost::numeric::ublas::vector<double>                  eers::math::zp(boost::numeric::ublas::vector<double> x, int len);
template boost::numeric::ublas::vector<std::complex<double> >   eers::math::zp(boost::numeric::ublas::vector<std::complex<double> > z, int len);
template boost::numeric::ublas::vector<std::complex<double> >   eers::math::complex(boost::numeric::ublas::vector<double> r);
template boost::numeric::ublas::vector<std::complex<double> >   eers::math::complex(boost::numeric::ublas::vector<double> r, boost::numeric::ublas::vector<double> i);
template boost::numeric::ublas::vector<std::complex<double> >   eers::math::pol2cart(boost::numeric::ublas::vector<double> r, boost::numeric::ublas::vector<double> p);
template boost::numeric::ublas::vector<std::complex<double> >   eers::math::fft(boost::numeric::ublas::vector<double> r, kiss_fft_cfg dir);
template boost::numeric::ublas::vector<std::complex<double> >   eers::math::fft(boost::numeric::ublas::vector<std::complex<double> > z, kiss_fft_cfg dir);
template boost::numeric::ublas::vector<double>                  eers::math::poly(boost::numeric::ublas::vector<double> roots);
template boost::numeric::ublas::vector<std::complex<double> >   eers::math::poly(boost::numeric::ublas::vector<std::complex<double> > roots);
template double                                                 eers::math::polyval(boost::numeric::ublas::vector<double> p, double x);
template std::complex<double>                                   eers::math::polyval(boost::numeric::ublas::vector<std::complex<double> > p, std::complex<double> x);
template bool                                                   eers::math::has_nan(boost::numeric::ublas::vector<double> x);
template bool                                                   eers::math::has_inf(boost::numeric::ublas::vector<double> x);
template boost::numeric::ublas::vector<double>                  eers::math::linspace(double start, double end, int num);
template boost::numeric::ublas::vector<double>                  eers::math::logspace(double start, double end, int num);
template double                                                 eers::math::rescale(double value, double source_min, double source_max, double dest_min, double dest_max);