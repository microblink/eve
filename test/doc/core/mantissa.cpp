#include <eve/function/mantissa.hpp>
#include <eve/constant/inf.hpp>
#include <eve/constant/minf.hpp>
#include <eve/constant/nan.hpp>
#include <eve/constant/mindenormal.hpp>
#include <eve/wide.hpp>

using wide_ft = eve::wide<float, eve::fixed<8>>;

int main()
{
  wide_ft pf = {-1.0f, 0.0f, 367.0f, -1005600.0f, eve::mindenormal(eve::as<float>()),
               eve::inf(eve::as<float>()), eve::minf(eve::as<float>()), eve::nan(eve::as<float>())};

  std::cout << "---- simd" << '\n'
            << "<- pf           = " << pf << '\n'
            << "-> mantissa(pf) = " << eve::mantissa(pf) << '\n';

  float        xf = -32768.0f;

  std::cout << "---- scalar" << '\n'
            << "<- xf           = " << xf << '\n'
            << "-> mantissa(xf) = " << eve::mantissa(xf) << '\n';
  return 0;
}
