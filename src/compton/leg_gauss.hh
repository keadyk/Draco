#include "print_handler.hh"
#include <sstream>
#include <vector>

#ifndef leg_gauss_hh
#define leg_gauss_hh

namespace leg_gauss {
// get the Gauss-Legendre quadrature points for a given (even) order:
static std::vector<double> get_gauss_legendre_pts(const size_t np) {
  std::vector<double> pts;
  switch (np) {
  case 2:
    pts.assign(1, 0.5773502691896257645091488);
    break;
  case 4:
    pts.assign(2, 0.0);
    pts[0] = 0.3399810435848562648026658;
    pts[1] = 0.8611363115940525752239465;
    break;
  case 6:
    pts.assign(3, 0.0);
    pts[0] = 0.2386191860831969086305017;
    pts[1] = 0.6612093864662645136613996;
    pts[2] = 0.9324695142031520278123016;
    break;
  case 8:
    pts.assign(4, 0.0);
    pts[0] = 0.1834346424956498049394761;
    pts[1] = 0.5255324099163289858177390;
    pts[2] = 0.7966664774136267395915539;
    pts[3] = 0.9602898564975362316835609;
    break;
  case 10:
    pts.assign(5, 0.0);
    pts[0] = 0.1488743389816312108848260;
    pts[1] = 0.4333953941292471907992659;
    pts[2] = 0.6794095682990244062343274;
    pts[3] = 0.8650633666889845107320967;
    pts[4] = 0.9739065285171717200779640;
    break;
  case 12:
    pts.assign(6, 0.0);
    pts[0] = 0.1252334085114689154724414;
    pts[1] = 0.3678314989981801937526915;
    pts[2] = 0.5873179542866174472967024;
    pts[3] = 0.7699026741943046870368938;
    pts[4] = 0.9041172563704748566784659;
    pts[5] = 0.9815606342467192506905491;
    break;
  case 14:
    pts.assign(7, 0.0);
    pts[0] = 0.1080549487073436620662447;
    pts[1] = 0.3191123689278897604356718;
    pts[2] = 0.5152486363581540919652907;
    pts[3] = 0.6872929048116854701480198;
    pts[4] = 0.8272013150697649931897947;
    pts[5] = 0.9284348836635735173363911;
    pts[6] = 0.9862838086968123388415973;
    break;
  default:
    std::ostringstream msg;
    msg << "Legendre points cannot be generated for np=" << np << std::endl;
    print::error(msg.str());
  }
  return pts;
}

// get the Gauss-Legendre quadrature weights for a given (even) order:
static std::vector<double> get_gauss_legendre_wts(const size_t np) {
  std::vector<double> wts;
  switch (np) {
  case 2:
    wts.assign(1, 1.0000000000000000000000000);
    break;
  case 4:
    wts.assign(2, 0.0);
    wts[0] = 0.6521451548625461426269361;
    wts[1] = 0.3478548451374538573730639;
    break;
  case 6:
    wts.assign(3, 0.0);
    wts[0] = 0.4679139345726910473898703;
    wts[1] = 0.3607615730481386075698335;
    wts[2] = 0.1713244923791703450402961;
    break;
  case 8:
    wts.assign(4, 0.0);
    wts[0] = 0.3626837833783619829651504;
    wts[1] = 0.3137066458778872873379622;
    wts[2] = 0.2223810344533744705443560;
    wts[3] = 0.1012285362903762591525314;
    break;
  case 10:
    wts.assign(5, 0.0);
    wts[0] = 0.2955242247147528701738930;
    wts[1] = 0.2692667193099963550912269;
    wts[2] = 0.2190863625159820439955349;
    wts[3] = 0.1494513491505805931457763;
    wts[4] = 0.0666713443086881375935688;
    break;
  case 12:
    wts.assign(6, 0.0);
    wts[0] = 0.2491470458134027850005624;
    wts[1] = 0.2334925365383548087608499;
    wts[2] = 0.2031674267230659217490645;
    wts[3] = 0.1600783285433462263346525;
    wts[4] = 0.1069393259953184309602547;
    wts[5] = 0.0471753363865118271946160;
    break;
  case 14:
    wts.assign(7, 0.0);
    wts[0] = 0.2152638534631577901958764;
    wts[1] = 0.2051984637212956039659241;
    wts[2] = 0.1855383974779378137417166;
    wts[3] = 0.1572031671581935345696019;
    wts[4] = 0.1215185706879031846894148;
    wts[5] = 0.0801580871597602098056333;
    wts[6] = 0.0351194603317518630318329;
    break;
  default:
    std::ostringstream msg;
    msg << "Legendre weights cannot be generated for np=" << np << std::endl;
    print::error(msg.str());
  }
  return wts;
}
}

#endif
