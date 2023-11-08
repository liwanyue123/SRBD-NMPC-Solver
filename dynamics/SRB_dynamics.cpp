
#include "SRB_dynamics.h"
#include "orientation_tool.h"

SRBDynamic::SRBDynamic()
{
  pf_.resize(6, 1);
  pf_.setZero();
  Rf_.resize(6, 3);
  Rf_.block<3, 3>(0, 0).setIdentity();
  Rf_.block<3, 3>(3, 0).setIdentity();
  m_ = 0.0;
  mu_ = 0.5;
  Lfx_ = 0.05;
  Lfz_ = 0.05;
  fmax_ = 1000;
  fmin_ = 0;
  dt_ = 0.0;
  Lbody_.resize(3, 3);
  Lbody_.diagonal() << 1, 1, 1;
  Nx_ = 12;
  Nu_ = 12;
  Ng_ = 24;
}

SRBDynamic::~SRBDynamic() {}

void SRBDynamic::SetFoot(const Vec &pr, const Vec &pl, const Mat &R0, const Mat &R1)
{
  pf_.block<3, 1>(0, 0) = pr.block<3, 1>(0, 0);
  pf_.block<3, 1>(3, 0) = pl.block<3, 1>(0, 0);
  Rf_.block<3, 3>(0, 0) = R0;
  Rf_.block<3, 3>(3, 0) = R1;
}

void SRBDynamic::SetMass(scale_t m)
{
  m_ = m;
}

void SRBDynamic::SetMPCdt(scale_t dt)
{
  dt_ = dt;
}

void SRBDynamic::SetInertia(const Mat &L)
{
  Lbody_ = L.inverse();
}

SRBDynamic::Vec SRBDynamic::GetFoot(int N)
{
  if (N == 0)
  {
    return pf_.block<3, 1>(0, 0);
  }
  else
  {
    return pf_.block<3, 1>(3, 0);
  }
}

SRBDynamic::Mat SRBDynamic::GetFootR(int N)
{
  if (N == 0)
  {
    return Rf_.block<3, 3>(0, 0);
  }
  else
  {
    return Rf_.block<3, 3>(3, 0);
  }
}

void SRBDynamic::GetDynamic(const Vec &x, const Vec &x_next, const Vec &u, Mat *pA, Mat *pB, Mat *pb, Mat *pf)
{
  // rk中的k表示下一时刻状态
  // rp中的p表示当前时刻状态
  Val rk = x_next.block<3, 1>(0, 0);
  Val wk = x_next.block<3, 1>(3, 0);
  Val xk = x_next.block<3, 1>(6, 0);
  Val vk = x_next.block<3, 1>(9, 0);
  Val rp = x.block<3, 1>(0, 0);
  Val wp = x.block<3, 1>(3, 0);
  Val xp = x.block<3, 1>(6, 0);
  Val vp = x.block<3, 1>(9, 0);
  Ual uk = u.block<12, 1>(0, 0);
  Jal Rk = expm(rk);
  Jal Rp = expm(rp);

  // Jal Rkpt = Rk*Rp.transpose();
  Val rkp = logm<scale_t>(Rk * (Rp.transpose()));
  Jal Jl = jl<scale_t>(rkp);
  Jal Jlt = jlt<scale_t>(rkp);
  Val pu = skew<scale_t>(GetFoot(0) - xp) * uk.block<3, 1>(0, 0) * dt_ +
           skew<scale_t>(GetFoot(1) - xp) * uk.block<3, 1>(6, 0) * dt_ + uk.block<3, 1>(3, 0) * dt_ + uk.block<3, 1>(9, 0) * dt_;

  if (pA)
  {
    *pA = Mat::Zero(Nx_, Nx_);
    (*pA).block<3, 3>(0, 0) = Jl * Jlt.transpose();
    (*pA).block<3, 3>(0, 3) = Jl * dt_;
    (*pA).block<3, 3>(3, 0) = Rp.transpose() * (skew<scale_t>(Lbody_ * Rp * pu) - Lbody_ * skew<scale_t>(Rp * pu));
    (*pA).block<3, 3>(3, 3).diagonal() << 1, 1, 1;
    (*pA).block<3, 3>(3, 6) = Rp.transpose() * Lbody_ * Rp * skew<scale_t>((uk.block<3, 1>(0, 0) + uk.block<3, 1>(6, 0)) * dt_);
    (*pA).block<3, 3>(6, 6).diagonal() << 1, 1, 1;
    (*pA).block<3, 3>(6, 9).diagonal() << dt_, dt_, dt_;
    (*pA).block<3, 3>(9, 9).diagonal() << 1, 1, 1;
  }

  if (pB)
  {
    (*pB) = Mat::Zero(Nx_, Nu_);
    (*pB).block<3, 3>(3, 0) = Rp * Lbody_ * Rp.transpose() * skew<scale_t>(GetFoot(0) - xp) * dt_;
    (*pB).block<3, 3>(3, 3) = Rp * Lbody_ * Rp.transpose() * dt_;
    (*pB).block<3, 3>(3, 6) = Rp * Lbody_ * Rp.transpose() * skew<scale_t>(GetFoot(1) - xp) * dt_;
    (*pB).block<3, 3>(3, 9) = Rp * Lbody_ * Rp.transpose() * dt_;
    (*pB).block<3, 3>(6, 0).diagonal() << 0.5 * dt_ * dt_ / m_, 0.5 * dt_ * dt_ / m_, 0.5 * dt_ * dt_ / m_;
    (*pB).block<3, 3>(6, 6).diagonal() << 0.5 * dt_ * dt_ / m_, 0.5 * dt_ * dt_ / m_, 0.5 * dt_ * dt_ / m_;
    (*pB).block<3, 3>(9, 0).diagonal() << dt_ / m_, dt_ / m_, dt_ / m_;
    (*pB).block<3, 3>(9, 6).diagonal() << dt_ / m_, dt_ / m_, dt_ / m_;
  }
  // 动力学推导：
  // 由f = 0，可以得到f收敛的离散递推公式：
  // d_{f}/d_{x_next} * dx_next + d_{f}/d_{x} * dx - f = 0 (1)
  // 由（1）可得
  // A =  ( d_{f} / d_{x_next} )^-1 * ( d_{f}/d_{x} )
  // B = ( d_{f} / d_{x_next} )^-1 * ( d_{f}/d_{u} )
  // b = ( d_{f} / d_{x_next} )^-1 * -f

  // 动力学方程：
  // f(1:3) = logm(Rk*Rp^T) - wp*dt 世界系轴角变化 = 世界系角速度
  // f(4:6) = wk - wp - Rp*Lbody*Rp^T*(cross(pr,F_r) + cross(pl,F_l) + Tau_r + Tau_l)*dt 角速度变化 = 躯干合力矩/躯干惯量
  // f(7:9) = xk - xp - vp*dt - 0.5*(F_r + F_l)/m*dt^2 - 0.5*g*dt^2 世界系位置变化
  // f(10:12) = vk - vp - (F_r + F_l)/m*dt - g*dt 世界系速度变化

  Vec f;

  f.resize(12);
  f.setZero();
  f.block<3, 1>(0, 0) = logm<scale_t>(Rk * (Rp.transpose())) - wp * dt_;

  f.block<3, 1>(3, 0) = wk - wp - Rp.transpose() * Lbody_ * Rp * skew<scale_t>(GetFoot(0) - xp) * uk.block<3, 1>(0, 0) * dt_ -
                        Rp.transpose() * Lbody_ * Rp * skew<scale_t>(GetFoot(1) - xp) * uk.block<3, 1>(6, 0) * dt_ -
                        Rp.transpose() * Lbody_ * Rp * uk.block<3, 1>(3, 0) * dt_ -
                        Rp.transpose() * Lbody_ * Rp * uk.block<3, 1>(9, 0) * dt_;
  // f.block<3, 1>(3, 0) = wk - wp - Rp*Lbody_*Rp.transpose()*skew<scale_t>(GetFoot(0)-xp)*uk.block<3, 1>(0, 0)*dt_ -
  //                       Rp*Lbody_*Rp.transpose()*skew<scale_t>(GetFoot(1)-xp)*uk.block<3, 1>(6, 0)*dt_ -
  //                       Rp*Lbody_*Rp.transpose()*uk.block<3, 1>(3, 0)*dt_ -
  //                       Rp*Lbody_*Rp.transpose()*uk.block<3, 1>(9, 0)*dt_ + Rp*Lbody_*Rp.transpose()*skew(wp)*Rp*Mbody_*Rp.transpose()*wp*dt_;

  f.block<3, 1>(6, 0) = xk - xp - vp * dt_ - 0.5 * (uk.block<3, 1>(0, 0) + uk.block<3, 1>(6, 0)) * dt_ * dt_ / m_ -
                        0.5 * Val(0.0, 0.0, -9.8) * dt_ * dt_;
  f.block<3, 1>(9, 0) = vk - vp - (uk.block<3, 1>(0, 0) + uk.block<3, 1>(6, 0)) * dt_ / m_ - Val(0.0, 0.0, -9.8) * dt_;

  if (pb)
  {
    *pb = -f;
    (*pb).block<3, 1>(0, 0) = Jl * ((*pb).block<3, 1>(0, 0));
  }

  if (pf)
  {
    *pf = f;
  }
}

void SRBDynamic::GetContinuousDynamic(const Vec &x, const Vec &u, const Mat &j_x, const Mat &j_u, Mat *dx, Mat *j_dx, Mat *j_du)
{
  Val r = x.block<3, 1>(0, 0);
  Val l = x.block<3, 1>(3, 0);
  Val p = x.block<3, 1>(6, 0);
  Val v = x.block<3, 1>(9, 0);
  Ual uf = u.block<12, 1>(0, 0);
  Jal R = expm(r);

  Jal Jlt = jlt<scale_t>(r);
  Val w = R * Lbody_ * R.transpose() * l;

  if (dx)
  {
    dx->resize(12, 1);
    (*dx).block<3, 1>(0, 0) = Jlt * w;

    (*dx).block<3, 1>(3, 0) = uf.block<3, 1>(3, 0) + uf.block<3, 1>(9, 0) +
                              skew<scale_t>(GetFoot(0) - p) * uf.block<3, 1>(0, 0) +
                              skew<scale_t>(GetFoot(1) - p) * uf.block<3, 1>(6, 0);

    (*dx).block<3, 1>(6, 0) = v;

    (*dx).block<3, 1>(9, 0) = (uf.block<3, 1>(0, 0) + uf.block<3, 1>(6, 0)) / m_ + Val(0.0, 0.0, -9.8);
  }
  if (!j_dx && !j_du)
  {
    return;
  }

  Jal djlt_x;
  Jal djlt_y;
  Jal djlt_z;

  djlt<scale_t>(r, &djlt_x, &djlt_y, &djlt_z);
  Jal djlt_w;
  djlt_w.block<3, 1>(0, 0) = djlt_x * w;
  djlt_w.block<3, 1>(0, 1) = djlt_y * w;
  djlt_w.block<3, 1>(0, 2) = djlt_z * w;

  Jal Jl = jl<scale_t>(r);

  Mat j_fx = Mat::Zero(12, 12);
  j_fx.block<3, 3>(0, 0) = djlt_w + Jlt * (R * Lbody_ * R.transpose() * skew(l) - skew(w)) * Jl;
  j_fx.block<3, 3>(0, 3) = Jlt * (R * Lbody_ * R.transpose());
  j_fx.block<3, 3>(3, 6) = skew<scale_t>(uf.block<3, 1>(0, 0) + uf.block<3, 1>(6, 0));
  j_fx.block<3, 3>(6, 9).diagonal() << 1.0, 1.0, 1.0;

  if (j_dx)
  {
    *j_dx = j_fx * j_x;
  }

  if (j_du)
  {
    Mat j_fu = Mat::Zero(12, 12);
    j_fu.block<3, 3>(3, 0) = skew<scale_t>(GetFoot(0) - p);
    j_fu.block<3, 3>(3, 3).diagonal() << 1.0, 1.0, 1.0;
    j_fu.block<3, 3>(3, 6) = skew<scale_t>(GetFoot(1) - p);
    j_fu.block<3, 3>(3, 9).diagonal() << 1.0, 1.0, 1.0;

    j_fu.block<3, 3>(9, 0).diagonal() << 1.0 / m_, 1.0 / m_, 1.0 / m_;
    j_fu.block<3, 3>(9, 6).diagonal() << 1.0 / m_, 1.0 / m_, 1.0 / m_;

    *j_du = j_fx * j_u + j_fu;
  }
}

void SRBDynamic::GetShootingDynamic(const Vec &x, const Vec &x_next, const Vec &u, Mat *pA, Mat *pB, Mat *pb, Mat *pf)
{
  // rk中的k表示下一时刻状态
  // rp中的p表示当前时刻状态
  Val rk = x_next.block<3, 1>(0, 0);
  Val lk = x_next.block<3, 1>(3, 0);
  Val pk = x_next.block<3, 1>(6, 0);
  Val vk = x_next.block<3, 1>(9, 0);
  Ual uk = u.block<12, 1>(0, 0);
  Jal Rk = expm(rk);

  // ode
  Mat k1 = Mat::Zero(12, 1);
  Mat k2 = Mat::Zero(12, 1);
  Mat k3 = Mat::Zero(12, 1);
  Mat k4 = Mat::Zero(12, 1);

  Mat k1_x = Mat::Zero(12, 12);
  Mat k2_x = Mat::Zero(12, 12);
  Mat k3_x = Mat::Zero(12, 12);
  Mat k4_x = Mat::Zero(12, 12);

  Mat k1_u = Mat::Zero(12, 12);
  Mat k2_u = Mat::Zero(12, 12);
  Mat k3_u = Mat::Zero(12, 12);
  Mat k4_u = Mat::Zero(12, 12);

  //

  Mat j_x0 = Mat::Identity(12, 12);
  Mat j_u0 = Mat::Zero(12, 12);
  GetContinuousDynamic(x, uk, j_x0, j_u0, &k1, &k1_x, &k1_u);
  GetContinuousDynamic(x + 0.5 * dt_ * k1, uk, j_x0 + 0.5 * dt_ * k1_x, j_u0 + 0.5 * dt_ * k1_u, &k2, &k2_x, &k2_u);
  GetContinuousDynamic(x + 0.5 * dt_ * k2, uk, j_x0 + 0.5 * dt_ * k2_x, j_u0 + 0.5 * dt_ * k2_u, &k3, &k3_x, &k3_u);
  GetContinuousDynamic(x + dt_ * k3, uk, j_x0 + dt_ * k3_x, j_u0 + dt_ * k3_u, &k4, &k4_x, &k4_u);

  Vec x_get = x + (dt_ / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
  Mat j_x = j_x0 + dt_ * k1_x; //(dt_/6.0)*(k1_x + 2.0*k2_x + 2.0*k3_x + k4_x);
  Mat j_u = j_u0 + dt_ * k1_u; //(dt_/6.0)*(k1_u + 2.0*k2_u + 2.0*k3_u + k4_u);

  Val rg = x_get.block<3, 1>(0, 0);
  Val lg = x_get.block<3, 1>(3, 0);
  Val pg = x_get.block<3, 1>(6, 0);
  Val vg = x_get.block<3, 1>(9, 0);
  Jal Rg = expm(rg);

  Vec f;

  f.resize(12);
  f.setZero();

  f.block<3, 1>(0, 0) = rk - rg; // logm<scale_t>(Rk*(Rg.transpose()));
  f.block<3, 1>(3, 0) = lk - lg;
  f.block<3, 1>(6, 0) = pk - pg;
  f.block<3, 1>(9, 0) = vk - vg;

  if (pA)
  {
    *pA = Mat::Zero(Nx_, Nx_);
    *pA = j_x;
  }

  if (pB)
  {
    *pB = Mat::Zero(Nx_, Nu_);
    *pB = j_u;
  }

  // 动力学推导：
  // 由f = 0，可以得到f收敛的离散递推公式：
  // d_{f}/d_{x_next} * dx_next + d_{f}/d_{x} * dx - f = 0 (1)
  // 由（1）可得
  // A =  ( d_{f} / d_{x_next} )^-1 * ( d_{f}/d_{x} )
  // B = ( d_{f} / d_{x_next} )^-1 * ( d_{f}/d_{u} )
  // b = ( d_{f} / d_{x_next} )^-1 * -f

  // 动力学方程：
  // f(1:3) = logm(Rk*Rp^T) - wp*dt 世界系轴角变化 = 世界系角速度
  // f(4:6) = wk - wp - Rp*Lbody*Rp^T*(cross(pr,F_r) + cross(pl,F_l) + Tau_r + Tau_l)*dt 角速度变化 = 躯干合力矩/躯干惯量
  // f(7:9) = xk - xp - vp*dt - 0.5*(F_r + F_l)/m*dt^2 - 0.5*g*dt^2 世界系位置变化
  // f(10:12) = vk - vp - (F_r + F_l)/m*dt - g*dt 世界系速度变化

  if (pb)
  {
    *pb = -f;
    // (*pb).block<3, 1>(0, 0) = Jl*((*pb).block<3, 1>(0, 0));
  }

  if (pf)
  {
    *pf = f;
  }
}

void SRBDynamic::GetConstrain(const Vec &u, Mat &Ac, Mat &f)
{
  // 两腿的摩擦锥约束
  Ac = Mat::Zero(Ng_, Nu_);
  Mat b = Mat::Zero(Ng_, 1);
  for (int leg = 0; leg < 2; ++leg)
  {
    Ac.block<12, 6>(12 * leg, 6 * leg) << -1, 0, mu_, 0, 0, 0,
        0, -1, mu_, 0, 0, 0,
        1, 0, mu_, 0, 0, 0,
        0, 1, mu_, 0, 0, 0,
        0, 0, -1, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        Lfx_ * Rf_.block<3, 1>(3 * leg, 2).transpose(), -Rf_.block<3, 1>(3 * leg, 1).transpose(),
        Lfx_ * Rf_.block<3, 1>(3 * leg, 2).transpose(), Rf_.block<3, 1>(3 * leg, 1).transpose(),
        Lfz_ * Rf_.block<3, 1>(3 * leg, 2).transpose(), -Rf_.block<3, 1>(3 * leg, 2).transpose(),
        Lfz_ * Rf_.block<3, 1>(3 * leg, 2).transpose(), Rf_.block<3, 1>(3 * leg, 2).transpose(),
        0, 0, 0, -Rf_.block<3, 1>(3 * leg, 0).transpose(),
        0, 0, 0, Rf_.block<3, 1>(3 * leg, 0).transpose();
    b(12 * leg + 4) = fmax_;
    b(12 * leg + 5) = -fmin_;
  }
  f = Ac * u + b;
}

void SRBDynamic::Barrier(scale_t value, scale_t mu, scale_t theta, scale_t *b, scale_t *db, scale_t *ddb)
{
  if (value > theta)
  {
    if (b)
    {
      *b = -mu * log(value);
    }
    if (db)
    {
      *db = -mu / value;
    }
    if (ddb)
    {
      *ddb = mu / (value * value);
    }
  }
  else
  {
    if (b)
    {
      *b = 0.5 * mu * (((value - 2.0 * theta) / theta) * ((value - 2.0 * theta) / theta) - 1.0) - mu * log(theta);
      ;
    }
    if (db)
    {
      *db = mu * (value - 2.0 * theta) / (theta * theta);
    }
    if (ddb)
    {
      *ddb = mu / (theta * theta);
    }
  }
}
