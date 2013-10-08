#include "implHypParab.hpp"
namespace Geometry
{
  double implHypParab (const bilinerarPatchPoints& patch, double const x, double const y, double const z)
  {
    Point3D const Pu = patch.P10 - patch.P00;
    Point3D const Pv = patch.P01 - patch.P00;
    Point3D const Puv = patch.P11 - patch.P10 + patch.P00 - patch.P01;

    const double a0 (patch.P00.x);
    const double a1 (Pu.x);
    const double a2 (Pv.x);
    const double a3 (Puv.x);
    const double b0 (patch.P00.y);
    const double b1 (Pu.y);
    const double b2 (Pv.y);
    const double b3 (Puv.y);
    const double c0 (patch.P00.z);;
    const double c1 (Pu.z);
    const double c2 (Pv.z);
    const double c3 (Puv.z);


    const double bs0 = b0 * b0;
    const double bs1 = b1 * b1;
    const double bs2 = b2 * b2;
    const double bs3 = b3 * b3;
    const double as0 = a0 * a0;
    const double as1 = a1 * a1;
    const double as2 = a2 * a2;
    const double as3 = a3 * a3;
    const double cs0 = c0 * c0;
    const double cs1 = c1 * c1;
    const double cs2 = c2 * c2;
    const double cs3 = c3 * c3;
    const double xs = x * x;
    const double ys = y * y;
    const double zs = z * z;


    return zs * (a1 * a2 * bs3 + as3 * b1 * b2 - a1 * a3 * b2 * b3 - a2 * a3 * b1 * b3) +
           ys * (a1 * a2 * cs3 + as3 * c1 * c2 - a1 * a3 * c2 * c3 - a2 * a3 * c1 * c3) +
           xs * (b1 * b2 * cs3 + bs3 * c1 * c2 - b1 * b3 * c2 * c3 - b2 * b3 * c1 * c3) +
           y * (x * (-a1 * b2 * cs3 - a2 * b1 * cs3 + a1 * b3 * c2 * c3 + a2 * b3 * c1 * c3 + a3 * b1 * c2 * c3 + a3 * b2 * c1 * c3 - a3 * b3 * c1 * c2 * 2.0) +
                as1 * b3 * cs2 + as2 * b3 * cs1 + a0 * a1 * b2 * cs3 + a0 * a2 * b1 * cs3 -
                a1 * a2 * b0 * cs3 * 2.0 - a1 * a3 * b1 * cs2 - a2 * a3 * b2 * cs1 - as3 * b0 * c1 * c2 * 2.0 +
                as3 * b1 * c0 * c2 + as3 * b2 * c0 * c1 - as2 * b1 * c1 * c3 - as1 * b2 * c2 * c3 - a0 * a1 * b3 * c2 * c3 -
                a0 * a2 * b3 * c1 * c3 - a0 * a3 * b1 * c2 * c3 - a0 * a3 * b2 * c1 * c3 + a0 * a3 * b3 * c1 * c2 * 2.0 +
                a1 * a2 * b1 * c2 * c3 + a1 * a2 * b2 * c1 * c3 + a1 * a2 * b3 * c0 * c3 * 2.0 - a1 * a2 * b3 * c1 * c2 * 2.0 +
                a1 * a3 * b0 * c2 * c3 * 2.0 - a1 * a3 * b2 * c0 * c3 + a1 * a3 * b2 * c1 * c2 - a1 * a3 * b3 * c0 * c2 +
                a2 * a3 * b0 * c1 * c3 * 2.0 - a2 * a3 * b1 * c0 * c3 + a2 * a3 * b1 * c1 * c2 - a2 * a3 * b3 * c0 * c1) -
           x * (-a3 * bs1 * cs2 - a3 * bs2 * cs1 + a0 * b1 * b2 * cs3 * 2.0 - a1 * b0 * b2 * cs3 -
                a2 * b0 * b1 * cs3 + a1 * b1 * b3 * cs2 + a2 * b2 * b3 * cs1 + a0 * bs3 * c1 * c2 * 2.0 -
                a1 * bs3 * c0 * c2 - a2 * bs3 * c0 * c1 + a1 * bs2 * c1 * c3 + a2 * bs1 * c2 * c3 -
                a0 * b1 * b3 * c2 * c3 * 2.0 - a0 * b2 * b3 * c1 * c3 * 2.0 + a1 * b0 * b3 * c2 * c3 -
                a1 * b1 * b2 * c2 * c3 + a1 * b2 * b3 * c0 * c3 - a1 * b2 * b3 * c1 * c2 + a2 * b0 * b3 * c1 * c3 -
                a2 * b1 * b2 * c1 * c3 + a2 * b1 * b3 * c0 * c3 - a2 * b1 * b3 * c1 * c2 + a3 * b0 * b1 * c2 * c3 +
                a3 * b0 * b2 * c1 * c3 - a3 * b0 * b3 * c1 * c2 * 2.0 - a3 * b1 * b2 * c0 * c3 * 2.0 +
                a3 * b1 * b2 * c1 * c2 * 2.0 + a3 * b1 * b3 * c0 * c2 + a3 * b2 * b3 * c0 * c1) +
           z * (y * (-as3 * b1 * c2 - as3 * b2 * c1 - a1 * a2 * b3 * c3 * 2.0 + a1 * a3 * b2 * c3 +
                     a1 * a3 * b3 * c2 + a2 * a3 * b1 * c3 + a2 * a3 * b3 * c1) +
                x * (-a1 * bs3 * c2 - a2 * bs3 * c1 + a1 * b2 * b3 * c3 + a2 * b1 * b3 * c3 -
                     a3 * b1 * b2 * c3 * 2.0 + a3 * b1 * b3 * c2 + a3 * b2 * b3 * c1) +
                as1 * bs2 * c3 + as2 * bs1 * c3 + a0 * a1 * bs3 * c2 +
                a0 * a2 * bs3 * c1 - a1 * a2 * bs3 * c0 * 2.0 - a1 * a3 * bs2 * c1 -
                a2 * a3 * bs1 * c2 + as3 * b0 * b1 * c2 + as3 * b0 * b2 * c1 - as3 * b1 * b2 * c0 * 2.0 -
                as2 * b1 * b3 * c1 - as1 * b2 * b3 * c2 - a0 * a1 * b2 * b3 * c3 - a0 * a2 * b1 * b3 * c3 +
                a0 * a3 * b1 * b2 * c3 * 2.0 - a0 * a3 * b1 * b3 * c2 - a0 * a3 * b2 * b3 * c1 +
                a1 * a2 * b0 * b3 * c3 * 2.0 - a1 * a2 * b1 * b2 * c3 * 2.0 + a1 * a2 * b1 * b3 * c2 +
                a1 * a2 * b2 * b3 * c1 - a1 * a3 * b0 * b2 * c3 - a1 * a3 * b0 * b3 * c2 + a1 * a3 * b1 * b2 * c2 +
                a1 * a3 * b2 * b3 * c0 * 2.0 - a2 * a3 * b0 * b1 * c3 - a2 * a3 * b0 * b3 * c1 + a2 * a3 * b1 * b2 * c1 +
                a2 * a3 * b1 * b3 * c0 * 2.0) -
           a0 * a3 * (bs1 * cs2 + bs2 * cs1) +
           a1 * a2 * (bs0 * cs3 + bs3 * cs0) + as0 * b1 * b2 * cs3 - as1 * b0 * b3 * cs2 - as2 * b0 * b3 * cs1 +
           as3 * b1 * b2 * cs0 + as0 * bs3 * c1 * c2 - as1 * bs2 * c0 * c3 - as2 * bs1 * c0 * c3 +
           as3 * bs0 * c1 * c2 - a0 * a1 * b0 * b2 * cs3 - a0 * a2 * b0 * b1 * cs3 + a0 * a1 * b1 * b3 * cs2 +
           a1 * a3 * b0 * b1 * cs2 + a0 * a2 * b2 * b3 * cs1 + a2 * a3 * b0 * b2 * cs1 - a1 * a3 * b2 * b3 * cs0 -
           a2 * a3 * b1 * b3 * cs0 - a0 * a1 * bs3 * c0 * c2 - a0 * a2 * bs3 * c0 * c1 + a0 * a1 * bs2 * c1 * c3 +
           a1 * a3 * bs2 * c0 * c1 + a0 * a2 * bs1 * c2 * c3 + a2 * a3 * bs1 * c0 * c2 - a1 * a3 * bs0 * c2 * c3 -
           a2 * a3 * bs0 * c1 * c3 - as3 * b0 * b1 * c0 * c2 - as3 * b0 * b2 * c0 * c1 + as2 * b0 * b1 * c1 * c3 +
           as2 * b1 * b3 * c0 * c1 + as1 * b0 * b2 * c2 * c3 + as1 * b2 * b3 * c0 * c2 - as0 * b1 * b3 * c2 * c3 -
           as0 * b2 * b3 * c1 * c3 + a0 * a1 * b0 * b3 * c2 * c3 - a0 * a1 * b1 * b2 * c2 * c3 + a0 * a1 * b2 * b3 * c0 * c3 -
           a0 * a1 * b2 * b3 * c1 * c2 + a0 * a2 * b0 * b3 * c1 * c3 - a0 * a2 * b1 * b2 * c1 * c3 + a0 * a2 * b1 * b3 * c0 * c3 -
           a0 * a2 * b1 * b3 * c1 * c2 + a0 * a3 * b0 * b1 * c2 * c3 + a0 * a3 * b0 * b2 * c1 * c3 -
           a0 * a3 * b0 * b3 * c1 * c2 * 2.0 - a0 * a3 * b1 * b2 * c0 * c3 * 2.0 +
           a0 * a3 * b1 * b2 * c1 * c2 * 2.0 + a0 * a3 * b1 * b3 * c0 * c2 + a0 * a3 * b2 * b3 * c0 * c1 -
           a1 * a2 * b0 * b1 * c2 * c3 - a1 * a2 * b0 * b2 * c1 * c3 - a1 * a2 * b0 * b3 * c0 * c3 * 2.0 +
           a1 * a2 * b0 * b3 * c1 * c2 * 2.0 + a1 * a2 * b1 * b2 * c0 * c3 * 2.0 - a1 * a2 * b1 * b3 * c0 * c2 -
           a1 * a2 * b2 * b3 * c0 * c1 + a1 * a3 * b0 * b2 * c0 * c3 - a1 * a3 * b0 * b2 * c1 * c2 +
           a1 * a3 * b0 * b3 * c0 * c2 - a1 * a3 * b1 * b2 * c0 * c2 + a2 * a3 * b0 * b1 * c0 * c3 -
           a2 * a3 * b0 * b1 * c1 * c2 + a2 * a3 * b0 * b3 * c0 * c1 - a2 * a3 * b1 * b2 * c0 * c1;
  }
}
