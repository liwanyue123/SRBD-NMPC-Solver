#ifndef STANCE_NMPC_H
#define STANCE_NMPC_H

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <string>
#include <iostream>
#include "timer/timer.h"
#include "hpipm-cpp/include/hpipm-cpp/hpipm-cpp.hpp"
#include "dynamic/flow_tool.h"
#include "dynamic/flow_dynamic.h"
#include "matplotlibcpp.h"
#include <yaml-cpp/yaml.h>
#include <thread>
#include <mutex>
template <typename T>
using Vec12 = Eigen::Matrix<T, 12, 1>;
template <typename T>
using Vec24 = Eigen::Matrix<T, 24, 1>;
class StanceNMPC
{
public:
    explicit StanceNMPC(const std::string &config_file);
    ~StanceNMPC();

    bool readYaml(const std::string &config_file);
    void initialize();
    void controlLoop();
    void drawPlot();

private:
    Eigen::Matrix<double, 12, 1> Q_read;
    Eigen::Matrix<double, 12, 1> QN_read;
    Eigen::Matrix<double, 3, 1> L_read;

    double Qu_read;
    double mu_barrier;
    double theta_barrier;
    double alpha_ocp_read;
    int loop_read;

    int Nrep; // 测试次数
    int N;    // 预测步长
    double T_step_read;
    bool show_flag_read;

    Eigen::Matrix<double, 12, 12> Q; // state
    Eigen::Matrix<double, 12, 12> S;
    Eigen::Matrix<double, 12, 12> R; // force
    Eigen::Matrix<double, 12, 12> QN;

    FlowDynamic flow_dynamic_;
    Eigen::MatrixXd x_nmpc;
    Eigen::MatrixXd u_nmpc;

    Eigen::VectorXd x0;
    Eigen::Matrix<double, 12, 1> x_ref_k;
    Eigen::MatrixXd x_ref;

    Eigen::MatrixXd get_solveX;

    Eigen::MatrixXd get_solveU;

    hpipm::OcpQpIpmSolverSettings solver_settings;

    void setupDynamics();

    void setupReference();

    void prepareQpStructures(std::vector<hpipm::OcpQp> &qp);

    void solveQpProblems(std::vector<hpipm::OcpQpSolution> &solution);

    void updateStateAndControl(std::vector<hpipm::OcpQpSolution> &solution);

    void printOptimizationInfo(int sqp_loop);
    bool checkConvergence();
    bool linearSearch();
    void prepareQpStep(int begin, int end);

    void initThread();
    void atomicPrint();
    // Helper methods for control calculations
    void solveNMPC();
    void applyControl(const Eigen::VectorXd &control);
    void updateState();

    Eigen::VectorXd x_mpc;
    Eigen::VectorXd x_mpc_next;
    Eigen::VectorXd u_mpc;

    Eigen::MatrixXd A_constrain, f_constrain;
    Eigen::MatrixXd f_constrain_all, f_all;
    // 障碍方程系数
    Eigen::VectorXd b_constrain, db_constrain, ddb_constrain; // 障碍方程
};

#endif // STANCE_NMPC_H
