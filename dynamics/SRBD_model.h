#ifndef MPC_SRBD_MODEL_H_
#define MPC_SRBD_MODEL_H_

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Core>
#include <vector>
#include <string>
#include <iostream>
// #include "timer/timer.h"

/*! @brief SRBD (Single Rigid Body Dynamics) model
 */
class SRBDModel
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

  SRBDModel();
  ~SRBDModel();

  /*! @brief Set positions and attitudes of both feet
   */
  void SetFoot(const Vec &pr, const Vec &pl, const Mat &R0, const Mat &R1);

  /*! @brief Set the mass
   */
  void SetMass(scale_t m);

  /*! @brief Set the time interval
   */
  void SetMPCdt(scale_t dt);

  /*! @brief Set the inverse of inertia,
   *         i.e., L = I^-1
   */
  void SetInertia(const Mat &L);

  /*! @brief Get foot position 0: right 1: left
   */
  Vec GetFoot(int N);

  /*! @brief Get foot attitude 0: right 1: left
   */
  Mat GetFootR(int N);

  /*! @brief Compute continuous dynamics
   *  @param x Input - current state x = [world axis angle, world angular velocity, world trunk position, world trunk velocity]
   *  @param u Input - foot forces u = [world force R, world torque R, world force L, world torque L]
   *  @param j_x Input - dx_dx0
   *  @param j_u Input - dx_du
   *  @param dx Output - derivative
   *  @param j_dx Output - ddx_dx
   *  @param j_du Output - ddx_du
   */
  void GetContinuousDynamic(const Vec &x, const Vec &u, const Mat &j_x, const Mat &j_u, Mat *dx, Mat *j_dx, Mat *j_du);

  /*! @brief Compute dynamics using multi-shooting method
   *  @param x Input - current state x = [world axis angle, world angular velocity, world trunk position, world trunk velocity]
   *  @param x_next Input - next state
   *  @param u Input - foot forces u = [world force R, world torque R, world force L, world torque L]
   *  @param A Output - discrete dynamics state matrix, SQP iteration should satisfy dx = A*dx + B*du + b
   *  @param B Output - discrete dynamics input matrix
   *  @param b Output - discrete dynamics error matrix
   */
  void GetShootingDynamic(const Vec &x, const Vec &x_next, const Vec &u, Mat *pA, Mat *pB, Mat *pb, Mat *pf);

  /*! @brief Compute friction cone constraints
   *  @param u Input - foot forces u = [world force R, world torque R, world force L, world torque L]
   *  @param Ac Output - Jacobian of friction cone constraints
   *  @param f Output - friction cone constraints, f > 0
   */
  void GetConstrain(const Vec &u, Mat &Ac, Mat &f);

  /*! @brief Standard barrier function
   *  @param value Input - constraint
   *  @param mu Input - coefficient mu
   *  @param theta Input - coefficient theta
   *  @param b Output - barrier function
   *  @param db Output - first-order derivative of the barrier function
   *  @param ddb Output - second-order derivative of the barrier function
   */
  void Barrier(scale_t value, scale_t mu, scale_t theta, scale_t *b, scale_t *db, scale_t *ddb);

private:
  Vec pf_;       //! Foot position
  Mat Rf_;       //! Foot attitude
  scale_t m_;    //! Mass
  scale_t dt_;   //! Time interval
  Mat Lbody_;    //! Inverse of inertia
  scale_t mu_;   //! Friction coefficient
  scale_t Lfx_;  //! Corresponding to l_zmp
  scale_t Lfz_;  //! Corresponding to mu_angle
  scale_t fmax_; //! Maximum foot force
  scale_t fmin_; //! Minimum foot force

  int Nx_; // ! Number of variables
  int Nu_; //! Number of inputs
  int Ng_; //! Constraint array
};

#endif // MPC_SRBD_MODEL_H_
