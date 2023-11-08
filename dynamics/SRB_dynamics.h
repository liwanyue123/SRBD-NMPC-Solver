#ifndef MPC_SRB_DYNAMICS_H_
#define MPC_SRB_DYNAMICS_H_

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Core>
#include <vector>
#include <string>
#include <iostream>
// #include "timer/timer.h"

/*! @brief SRBD 动力学模型
 */
class SRBDynamic
{
public:
  using scale_t = double;
  using Mat = Eigen::Matrix<scale_t, Eigen::Dynamic, Eigen::Dynamic>;
  using Vec = Eigen::Matrix<scale_t, Eigen::Dynamic, 1>;
  using Jal = Eigen::Matrix<scale_t, 3, 3>;
  using Val = Eigen::Matrix<scale_t, 3, 1>;
  using Ual = Eigen::Matrix<scale_t, 12, 1>;
  using Box = Eigen::Matrix<int, Eigen::Dynamic, 1>;
  using VecV = std::vector<Vec>;

  SRBDynamic();
  ~SRBDynamic();

  /*! @brief 设定两足端位置与姿态
   */
  void SetFoot(const Vec &pr, const Vec &pl, const Mat &R0, const Mat &R1);

  /*! @brief 设定质量
   */
  void SetMass(scale_t m);

  /*! @brief 设定时间间隔
   */
  void SetMPCdt(scale_t dt);

  /*! @brief 设定惯量的逆，
   *         即L = I^-1
   */
  void SetInertia(const Mat &L);

  /*! @brief 获取足端位置 0: right 1: left
   */
  Vec GetFoot(int N);

  /*! @brief 获取足端姿态 0: right 1: left
   */
  Mat GetFootR(int N);

  /*! @brief 计算动力学
   *  @param x 输入-当前状态 x = [世界系轴角,世界系角速度，世界系躯干位置，世界系躯干速度]
   *  @param x_next 输入-下一步状态
   *  @param u 输入-足底力 u = [世界系力R,世界系力矩R，世界系力L，世界系力矩L]
   *  @param A 输出-离散动力学状态矩阵，SQP迭代应满足 dx = A*dx + B*du + b
   *  @param B 输出-离散动力学输入矩阵
   *  @param b 输出-离散动力学误差矩阵
   */
  void GetDynamic(const Vec &x, const Vec &x_next, const Vec &u, Mat *pA, Mat *pB, Mat *pb, Mat *pf);

  /*! @brief 计算连续动力学
   *  @param x 输入-当前状态 x = [世界系轴角,世界系角速度，世界系躯干位置，世界系躯干速度]
   *  @param u 输入-足底力 u = [世界系力R,世界系力矩R，世界系力L，世界系力矩L]
   *  @param j_x 输入-dx_dx0
   *  @param j_u 输入-dx_du
   *  @param dx 输出-导数
   *  @param j_dx 输出-ddx_dx
   *  @param j_du 输出-ddx_du
   */
  void GetContinuousDynamic(const Vec &x, const Vec &u, const Mat &j_x, const Mat &j_u, Mat *dx, Mat *j_dx, Mat *j_du);

  /*! @brief 计算动力学,multi-shooting方法
   *  @param x 输入-当前状态 x = [世界系轴角,世界系角速度，世界系躯干位置，世界系躯干速度]
   *  @param x_next 输入-下一步状态
   *  @param u 输入-足底力 u = [世界系力R,世界系力矩R，世界系力L，世界系力矩L]
   *  @param A 输出-离散动力学状态矩阵，SQP迭代应满足 dx = A*dx + B*du + b
   *  @param B 输出-离散动力学输入矩阵
   *  @param b 输出-离散动力学误差矩阵
   */
  void GetShootingDynamic(const Vec &x, const Vec &x_next, const Vec &u, Mat *pA, Mat *pB, Mat *pb, Mat *pf);

  /*! @brief 计算摩擦锥约束
   *  @param u 输入-足底力 u = [世界系力R,世界系力矩R，世界系力L，世界系力矩L]
   *  @param Ac 输出-摩擦锥约束雅克比
   *  @param f 输出-摩擦锥约束， f > 0
   */
  void GetConstrain(const Vec &u, Mat &Ac, Mat &f);

  /*! @brief 标准障碍方程
   *  @param value 输入-约束
   *  @param mu 输入-系数mu
   *  @param theta 输入-系数theta
   *  @param b 输出-障碍方程
   *  @param db 输出-障碍方程一阶导数
   *  @param ddb 输出-障碍方程二阶导数
   */
  void Barrier(scale_t value, scale_t mu, scale_t theta, scale_t *b, scale_t *db, scale_t *ddb);

private:
  Vec pf_;       //! 足端位置
  Mat Rf_;       //! 足端姿态
  scale_t m_;    //! 质量
  scale_t dt_;   //! 时间间隔
  Mat Lbody_;    //! 惯量的逆
  scale_t mu_;   //! 摩擦系数
  scale_t Lfx_;  //! 即l_zmp
  scale_t Lfz_;  //! 即mu_angle
  scale_t fmax_; //! 最大足底力
  scale_t fmin_; //! 最小足底力

  int Nx_; // ！变量数目
  int Nu_; //! 输入量数目
  int Ng_; //! 约束数组
};

#endif // MPC_SRB_DYNAMICS_H_