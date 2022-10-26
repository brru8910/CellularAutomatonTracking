#ifndef _LocalTrackFinderBR_Utl_h_
#define _LocalTrackFinderBR_Utl_h_

namespace iputl {

struct xyz {
  double x, y, z;

  xyz(double _x, double _y, double _z)
  : x {_x}, y {_y}, z {_z}
  { }

  xyz
  operator * (double a) const noexcept
  { return {a*x, a*y, a*z}; }

  xyz
  operator / (double a) const noexcept
  { return *this * (1/a); }

  double
  operator * (const xyz &other) const noexcept
  { return x*other.x + y*other.y + z*other.z; }

  double
  mag2() const noexcept
  { return (*this) * *this; }

  double
  mag() const noexcept
  { return sqrt(mag2()); }

  xyz
  project_xz() const noexcept
  { return {x, 0, z}; }

  xyz
  project_yz() const noexcept
  { return {0, y, z}; }

  double
  angle(const xyz &other) const noexcept
  { return acos( (*this)*other / (this->mag() * other.mag()) ); }
}; // struct xyz


struct line {
  double x0, y0, z0, dx, dy, dz;

  xyz
  direction() const noexcept
  { return {dx, dy, dz}; }
}; // struct line


static double
fit_with_line(const std::vector<xyz> &pts, const std::vector<xyz> &ers, line &l)
{
  if (pts.size() != ers.size())
    abort();

  l.z0 = pts[0].z;
  l.dz = 1;

  // x-part
  double xls2 = 0, xs2 = 0, l2s2 = 0, ls2 = 0, s2 = 0;
  const int n = pts.size();
  for (int i = 0; i < n; ++i)
  {
    const double li = pts[i].z - l.z0;
    const double xi = pts[i].x;
    const double si2 = ers[i].x;
    xls2 += xi*li/si2;
    xs2 += xi/si2;
    l2s2 += li*li/si2;
    ls2 += li/si2;
    s2 += 1/si2;
  }
  l.x0 = (xls2 - xs2*l2s2/ls2)/(ls2 - s2*l2s2/ls2);
  l.dx = (xs2 - l.x0*s2)/ls2;

  // y-part
  double yls2 = 0, ys2 = 0;
  l2s2 = ls2 = s2 = 0;
  for (int i = 0; i < n; ++i)
  {
    const double li = pts[i].z - l.z0;
    const double yi = pts[i].y;
    const double si2 = ers[i].y;
    yls2 += yi*li/si2;
    ys2 += yi/si2;
    l2s2 += li*li/si2;
    ls2 += li/si2;
    s2 += 1/si2;
  }
  l.y0 = (yls2 - ys2*l2s2/ls2)/(ls2 - s2*l2s2/ls2);
  l.dy = (ys2 - l.y0*s2)/ls2;

  // chi2
  double chi2 = 0;
  for (int i = 0; i < n; ++i)
  {
    const double li = pts[i].z - l.z0;
    const double x = l.x0 + li*l.dx;
    const double y = l.y0 + li*l.dy;

    const double xi = pts[i].x;
    const double sxi = ers[i].x;
    const double yi = pts[i].y;
    const double syi = ers[i].y;

    const double chi2x = (x - xi)*(x - xi)/(sxi*sxi);
    const double chi2y = (y - yi)*(y - yi)/(syi*syi);
    chi2 += chi2x + chi2y;
  }

  // normalize direction
  xyz linedir = l.direction();
  linedir = linedir / linedir.mag();
  l.dx = linedir.x;
  l.dy = linedir.y;
  l.dz = linedir.z;

  return chi2;
} // fit_with_line()

} // namespace iputl

#endif
