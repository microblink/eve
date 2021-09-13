#include <eve/function/cyl_bessel_jy.hpp>
#include <eve/wide.hpp>
#include <eve/constant/inf.hpp>
#include <eve/constant/minf.hpp>
#include <eve/constant/nan.hpp>
#include <iostream>

using wide_ft = eve::wide<double, eve::fixed<8>>;

int main()
{

  wide_ft x  = {0.5, 0.0, 0.1, 1.0, 3.0, 0.5, 21.5, 10.0};
  wide_ft nu = {0.0, 1.0, 2.0, 3.0, 0.5, 1.5,  2.5,  3.5};

  auto [j, jp, n, np] = eve::cyl_bessel_jy(nu, x);
  std::cout << "---- simd" << '\n'
            << "<- nu  = " << nu << '\n'
            << "<- x   = " << x << '\n'
            << "   with: [j, jp, n, np] = eve::cyl_bessel_jy(nu, x);" << '\n'
            << "-> j   = " << j  << '\n'
            << "-> jp  = " << jp << '\n'
            << "-> n   = " << n  << '\n'
            << "-> np  = " << np << '\n';

  double snu = 1.5, sx = 0.5;
  auto [sj, sjp, sn, snp] = eve::cyl_bessel_jy(snu, sx);
  std::cout << "---- scalar" << '\n'
            << "<- snu  = " << snu << '\n'
            << "<- sx   = " << sx << '\n'
            << "   with: [j, jp, n, np] = eve::cyl_bessel_jy(snu, sx);" << '\n'
            << "-> sj   = " << sj  << '\n'
            << "-> sjp  = " << sjp << '\n'
            << "-> sn   = " << sn  << '\n'
            << "-> snp  = " << snp << '\n';
  return 0;
}
