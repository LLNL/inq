/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__OBSERVABLES__SPECTRUM
#define INQ__OBSERVABLES__SPECTRUM

/*
 Copyright (C) 2019 Xavier Andrade, Alfredo Correa.

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <inq_config.h>

#include <magnitude/energy.hpp>
#include <math/array.hpp>
#include <math/complex.hpp>

namespace inq {
namespace observables {

template <typename TimeType, typename TimeSeriesType, typename RetElementType = decltype(exp(complex{0.0, 1.0})*TimeSeriesType{}[0])>
math::array<RetElementType, 1> spectrum(quantity<magnitude::energy> maxw, quantity<magnitude::energy> dw, TimeType const & time, TimeSeriesType const & time_series) {

	CALI_CXX_MARK_FUNCTION;

  assert(time.size() == time_series.size());

  long ntime = time.size();
  long nfreq = maxw/dw + 1;
  
  math::array<RetElementType, 1> freq_series(nfreq);

  assert(freq_series.size() == nfreq);

  gpu::run(nfreq,
           [fse = begin(freq_series), tim = begin(time), tse = begin(time_series), ntime, dw] GPU_LAMBDA (auto ifreq){
             
             double ww = dw.in_atomic_units()*ifreq;
             
             RetElementType sum = 0.5*(tim[1] - tim[0])*tse[0];
             for(long itime = 1; itime < ntime - 1; itime++){
               assert(tim[itime] > tim[itime - 1]);
               sum += 0.5*(tim[itime + 1] - tim[itime - 1])*exp(complex{0.0, 1.0}*ww*tim[itime])*tse[itime];
             }
             sum += 0.5*(tim[ntime - 1] - tim[ntime - 2])*exp(complex{0.0, 1.0}*ww*tim[ntime - 1])*tse[ntime - 1];
             
             fse[ifreq] = sum;
           }
           );

  return freq_series;
}


}
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifdef INQ_OBSERVABLES_SPECTRUM_UNIT_TEST
#undef INQ_OBSERVABLES_SPECTRUM_UNIT_TEST

#include <catch2/catch_all.hpp>
#include <ions/unit_cell.hpp>
#include <fstream>

using namespace inq;
using namespace magnitude;
using namespace Catch::literals;

TEST_CASE("observables::spectrum", "[observables::spectrum]") {

  int ntime = 1000;
  double dtime = 0.1;

  math::array<double, 1> time(ntime);
  math::array<double, 1> tseries(ntime);
  double freq1 = 10.0;
  double freq2 = 6.39;
  double amp1 = 2.0;
  double amp2 = -1.5;
  
  for(int itime = 0; itime < ntime; itime++){
    time[itime] = dtime*itime;
    tseries[itime] = amp1*cos(time[itime]*freq1) + amp2*sin(time[itime]*freq2);
  }

  auto maxw = 20.0_Ha;
  auto dw = 0.1_Ha;
  
  auto fseries = observables::spectrum(maxw, dw, time, tseries);

  auto nfreq = fseries.size();

  CHECK(nfreq == 201);

  std::ofstream file("spectrum.dat");
  
  for(int ifreq = 0; ifreq < nfreq; ifreq++){
    file << ifreq*dw.in_atomic_units() << '\t' << real(fseries[ifreq]) << '\t' << imag(fseries[ifreq]) << std::endl;
  }
  
}

#endif
#endif
