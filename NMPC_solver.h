#ifndef STANCE_NMPC_H
#define STANCE_NMPC_H

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <string>
#include <iostream>
#include "utils/timer.h"
#include "hpipm-cpp/include/hpipm-cpp/hpipm-cpp.hpp"
#include "dynamics/orientation_tool.h"
#include "dynamics/SRB_dynamics.h"

#include <yaml-cpp/yaml.h>
#include <thread>
#include <mutex>
#include "utils/cppTypes.h"

class NMPCSolver
{
public:
    explicit NMPCSolver(const std::string &config_file);
    ~NMPCSolver();

    void initialize();
    void controlLoop();

private:
    void setupDynamics();

    void setupReference();

    void prepareQpStructures(std::vector<hpipm::OcpQp> &qp);

    void solveQpProblems(std::vector<hpipm::OcpQp> &qp, std::vector<hpipm::OcpQpSolution> &solution);

    bool checkConvergence();

    bool linearSearch();

    void printOptimizationInfo(int sqp_max_loop_);

    bool readYaml(const std::string &config_file);

    // NMPC
    Eigen::Matrix<double, 12, 1> Q_read_;
    Eigen::Matrix<double, 12, 1> Qf_read_;
    double R_read_;

    Eigen::Matrix<double, 12, 12> Q_; // state
    Eigen::Matrix<double, 12, 12> S_;
    Eigen::Matrix<double, 12, 12> R_; // force
    Eigen::Matrix<double, 12, 12> Qf_;

    int N_test_rep_;   // 测试次数
    int sqp_max_loop_; // sqp 迭代最大次数
    int N_;            // 预测步长
    double dt_MPC_;

    Eigen::MatrixXd x_nmpc_;
    Eigen::MatrixXd u_nmpc_;

    Eigen::VectorXd x0_;
    Eigen::MatrixXd x_ref_;
    Eigen::Matrix<double, 12, 1> x_ref_k_;

    Eigen::MatrixXd x_solved_;
    Eigen::MatrixXd u_solved_;

    Eigen::VectorXd x_mpc_;
    Eigen::VectorXd x_mpc_next_;
    Eigen::VectorXd u_mpc_;

    hpipm::OcpQpIpmSolverSettings solver_settings_;

    SRBDynamic flow_dynamic_;
    Eigen::Matrix<double, 3, 1> L_;

    // barrier_function
    Eigen::MatrixXd A_constrain_, f_constrain_, f_constrain_all_, f_all_;
    Eigen::VectorXd b_constrain_, db_constrain_, ddb_constrain_;
    double mu_barrier_;
    double theta_barrier_;

    // linear search
    Eigen::MatrixXd Jphi_x_;
    Eigen::MatrixXd Jphi_u_;
    Eigen::VectorXd x_search_;
    Eigen::VectorXd x_search_next_;
    Eigen::VectorXd u_search_;
    Eigen::MatrixXd f_search_;
    Eigen::VectorXd f_search_v_;

    double theta_max_ = 1e-6;
    double theta_min_ = 5e-10;
    double eta_ = 1e-4;
    double byta_phi_ = 1e-6;
    double byta_theta_ = 1e-6;
    double byta_alpha_ = 0.5;
    double alpha_min_ = 1e-4;
    double alpha_ = 1.0;
    double theta_; // 等式约束
    double phi_;   // cost
    double dphi_;
};

#endif // STANCE_NMPC_H
