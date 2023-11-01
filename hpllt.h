#ifndef HPLLT_H_
#define HPLLT_H_

#include <eigen3/Eigen/Eigen>
#include "blasfeo.h"
/*! @brief 高性能cholesky分解，同样使用AVX指令集的情况下，求解速度约为Eigen的3倍 
 *  通过 solve(M, v) 求解 x = M/v
 */
class HpChol {
 public:
  HpChol() { 
    m_ = -1; 
  }
  ~HpChol() {
    if (m_ > 0) {
      blasfeo_free_dmat(&sA_);
      blasfeo_free_dmat(&sB_);
      blasfeo_free_dmat(&sC_);
    }
  }
  Eigen::VectorXd solve(Eigen::MatrixXd& M, Eigen::VectorXd& v) {
    int m = v.rows();
    if (m != m_) {
      if (m_ > 0) {
        blasfeo_free_dmat(&sA_);
        blasfeo_free_dmat(&sB_);
        blasfeo_free_dmat(&sC_);
      }
      blasfeo_allocate_dmat(m, m, &sA_);
      blasfeo_allocate_dmat(m, m, &sB_);
      blasfeo_allocate_dmat(m, 1, &sC_);
      m_ = m;
    }

    //
    blasfeo_pack_dmat(m,m,M.data(),m,&sA_,0,0);
    blasfeo_pack_dmat(m,1,v.data(),m,&sC_,0,0);
    blasfeo_dpotrf_l(m, &sA_, 0, 0, &sB_, 0, 0);
    blasfeo_dtrsm_llnn(m, 1, 1.0, &sB_, 0, 0, &sC_, 0, 0, &sC_, 0, 0);
    blasfeo_dtrsm_lltn(m, 1, 1.0, &sB_, 0, 0, &sC_, 0, 0, &sC_, 0, 0);
    Eigen::VectorXd solve;
    solve.resize(m);

    blasfeo_unpack_dmat(m,1,&sC_,0,0,solve.data(),m);
    return solve;
  }
 private:
  blasfeo_dmat sA_;
  blasfeo_dmat sB_;
  blasfeo_dmat sC_;
  int m_;
};

#endif  // HPLLT_H_