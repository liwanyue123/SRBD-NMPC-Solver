#ifndef FLOW_TOOL_H_
#define FLOW_TOOL_H_

#include <iostream>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Core>

/*! @brief 打印eigen矩阵
 */
template<typename T = int>
void print(const Eigen::MatrixXd& H, const std::string& s = "") {
  std::cout << s << std::endl << H << std::endl;
}

/*! @brief X轴旋转矩阵
 */
template<typename T>
Eigen::Matrix<T, 3, 3> rotx(T val) {
  Eigen::Matrix<T, 3, 3> R;
  R << 1, 0, 0,
       0,  cos(val), -sin(val),
       0,  sin(val),  cos(val);
  return R;
}

/*! @brief Y轴旋转矩阵
 */
template<typename T>
Eigen::Matrix<T, 3, 3> roty(T val) {
  Eigen::Matrix<T, 3, 3> R;
  R << cos(val),  0, sin(val),
              0,  1,        0,
      -sin(val),  0, cos(val);
  return R;
}

/*! @brief Z轴旋转矩阵
 */
template<typename T>
Eigen::Matrix<T, 3, 3> rotz(T val) {
  Eigen::Matrix<T, 3, 3> R;
  R << cos(val), -sin(val), 0,
       sin(val),  cos(val), 0,
              0,         0, 1;
  return R;
}

/*! @brief 向量转反对称矩阵
 */
template<typename T>
Eigen::Matrix<T, 3, 3> skew(const Eigen::Matrix<T, 3, 1>& val) {
  Eigen::Matrix<T, 3, 3> R;
  R << 0, -val(2), val(1),
       val(2),  0, -val(0),
       -val(1),         val(0), 0;
  return R;
}

/*! @brief 反对称矩阵转向量
 */
template<typename T>
Eigen::Matrix<T, 3, 1> skewt(const Eigen::Matrix<T, 3, 3>& val) {
  return Eigen::Matrix<T, 3, 1>(-val(1, 2), val(0, 2), -val(0, 1));
}

/*! @brief so3 -> SO3
 */
template<typename T>
Eigen::Matrix<T, 3, 3> expm(const Eigen::Matrix<T, 3, 1>& val) {
  T theta = sqrt(val(0)*val(0) + val(1)*val(1) + val(2)*val(2));
  T h = 1e-10;
  if (theta < h) {
    theta = h;
  }
  Eigen::Matrix<T, 3, 3> V = skew(val);
  return Eigen::Matrix<T, 3, 3>::Identity() + (sin(theta)/theta)*V + ((1.0-cos(theta))/(theta*theta))*(V*V);
}

/*! @brief SO3 -> so3
 */
template<typename T>
Eigen::Matrix<T, 3, 1> logm(const Eigen::Matrix<T, 3, 3>& val) {
  T acosinput = (val(0, 0) + val(1, 1) + val(2, 2) - 1.0) / 2.0;
  T h = 1e-10;
  if (acosinput >= 1.0) {
    return Eigen::Matrix<T, 3, 1>(0, 0, 0);
  } else if (acosinput <= -1.0) {
    Eigen::Matrix<T, 3, 1> omg;
    omg.setZero();
    if (fabs(1 + val(2,2)) > h) {
      omg = (1.0/sqrt(2.0*(1.0+val(2,2))))*Eigen::Matrix<T, 3, 1>(val(0,2), val(1,2), 1.0+val(2,2));
    } else if (fabs(1.0 + val(1,1)) > h) {
      omg = (1.0/sqrt(2.0*(1.0+val(1,1))))*Eigen::Matrix<T, 3, 1>(val(0,1), 1.0+val(1,1), val(2,1));
    } else {
      omg = (1.0/sqrt(2.0*(1.0+val(0,0))))*Eigen::Matrix<T, 3, 1>(1.0+val(0,0), val(1,0), val(2,0));
    }
    return M_PI*omg;
  } else {
    T theta = acos(acosinput);
    Eigen::Matrix<T, 3, 3> R = val - val.transpose();
    return theta * (1.0 / (2.0 * sin(theta))) * skewt(R);
  }
  return Eigen::Matrix<T, 3, 1>(0, 0, 0);
}

/*! @brief 左雅克比 w = jl*d_{so3}/d_{t}， 世界系角速度 = 左雅克比 * 李代数的导数 
 */
template<typename T>
Eigen::Matrix<T, 3, 3> jl(const Eigen::Matrix<T, 3, 1>& val) {
  T theta = sqrt(val(0)*val(0) + val(1)*val(1) + val(2)*val(2));
  T h = 1e-10;
  if (theta < h) {
    theta = h;
  }
  Eigen::Matrix<T, 3, 3> V = skew(val)/theta;
  Eigen::Matrix<T, 3, 3> I = Eigen::Matrix<T, 3, 3>::Identity();
  return (sin(theta)/theta)*I + (1.0-(sin(theta)/theta))*(V*V+I) + ((1.0-cos(theta))/(theta))*(V);
}

/*! @brief 左雅克比的逆
 */
template<typename T>
Eigen::Matrix<T, 3, 3> jlt(const Eigen::Matrix<T, 3, 1>& val) {
  T theta = sqrt(val(0)*val(0) + val(1)*val(1) + val(2)*val(2));
  T h = 1e-10;
  if (theta < h) {
    theta = h;
  }
  Eigen::Matrix<T, 3, 3> V = skew(val)/theta;
  Eigen::Matrix<T, 3, 3> I = Eigen::Matrix<T, 3, 3>::Identity();
  T cot = 1.0/tan(0.5*theta);
  return (0.5*cot*theta)*I + (1.0-(0.5*cot*theta))*(V*V+I) - (0.5*(theta))*(V);
}

/*! @brief 左雅克比的导数
 *  @param djlx d_jl_d_so3(1)
 *  @param djly d_jl_d_so3(2)
 *  @param djlz d_jl_d_so3(3)
 */
template<typename T>
void djl(const Eigen::Matrix<T, 3, 1>& val, Eigen::Matrix<T, 3, 3>* djlx, Eigen::Matrix<T, 3, 3>* djly, Eigen::Matrix<T, 3, 3>* djlz) {
  T theta = sqrt(val(0)*val(0) + val(1)*val(1) + val(2)*val(2));
  T h = 1e-10;
  if (theta < h) {
    theta = h;
  }

  Eigen::Matrix<T, 3, 3> V = skew(val)/theta;
  Eigen::Matrix<T, 3, 3> I = Eigen::Matrix<T, 3, 3>::Identity();
  T cot_ = 1.0/tan(0.5*theta);
  T sin_ = sin(theta);
  T cos_ = cos(theta);

  T theta2 = theta*theta;
  T theta3 = theta2*theta;
  // Eigen::Matrix<T, 3, 3> djl_Q = ((theta*cos_ - sin_)/(theta3))*I;

  // Eigen::Matrix<T, 3, 3> djl_W = 2.0*((theta - sin_)/(theta3))*I;

  // Eigen::Matrix<T, 3, 3> djl_E =  (-(2.0*theta - 3.0*sin_ + theta*cos_)/(theta3))*(V*V+I);

  // Eigen::Matrix<T, 3, 3> djl_R = ((theta*sin_ + (2.0*(cos_ - 1.0)))/(theta3))*V;
                                   
  Eigen::Matrix<T, 3, 3> djl_base = ((theta*sin_ + (2.0*(cos_ - 1.0)))/(theta3))*(V) + (-(2.0*theta - 3.0*sin_ + theta*cos_)/(theta3))*(V*V);

  Eigen::Matrix<T, 3, 3> skew_x = skew<T>(Eigen::Matrix<T, 3, 1>(1,0,0));
  Eigen::Matrix<T, 3, 3> skew_y = skew<T>(Eigen::Matrix<T, 3, 1>(0,1,0));
  Eigen::Matrix<T, 3, 3> skew_z = skew<T>(Eigen::Matrix<T, 3, 1>(0,0,1));

  Eigen::Matrix<T, 3, 3> djl_x = ((theta-sin_)/(theta3))*(skew_x*skew<T>(val)+skew<T>(val)*skew_x) + ((1.0-cos_)/(theta2))*skew_x;
  Eigen::Matrix<T, 3, 3> djl_y = ((theta-sin_)/(theta3))*(skew_y*skew<T>(val)+skew<T>(val)*skew_y) + ((1.0-cos_)/(theta2))*skew_y;
  Eigen::Matrix<T, 3, 3> djl_z = ((theta-sin_)/(theta3))*(skew_z*skew<T>(val)+skew<T>(val)*skew_z) + ((1.0-cos_)/(theta2))*skew_z;

  *(djlx) = djl_x + (djl_base)*val(0);
  *(djly) = djl_y + (djl_base)*val(1);
  *(djlz) = djl_z + (djl_base)*val(2);
  // return (0.5*cot*theta)*I + (1.0-(0.5*cot*theta))*(V*V+I) - (0.5*(theta))*(V);
}

/*! @brief 左雅克比的逆的导数
 *  @param djltx d_jlt_d_so3(1)
 *  @param djlty d_jlt_d_so3(2)
 *  @param djltz d_jlt_d_so3(3)
 */
template<typename T>
void djlt(const Eigen::Matrix<T, 3, 1>& val, Eigen::Matrix<T, 3, 3>* djltx, Eigen::Matrix<T, 3, 3>* djlty, Eigen::Matrix<T, 3, 3>* djltz) {

  Eigen::Matrix<T, 3, 3> jlt_ = jlt<T>(val);

  Eigen::Matrix<T, 3, 3> djl_x;
  Eigen::Matrix<T, 3, 3> djl_y;
  Eigen::Matrix<T, 3, 3> djl_z;

  djl<T>(val, &djl_x, &djl_y, &djl_z);

  *(djltx) = -jlt_*djl_x*jlt_;
  *(djlty) = -jlt_*djl_y*jlt_;
  *(djltz) = -jlt_*djl_z*jlt_;
  // return (0.5*cot*theta)*I + (1.0-(0.5*cot*theta))*(V*V+I) - (0.5*(theta))*(V);
}

#endif // FLOW_TOOL_H_