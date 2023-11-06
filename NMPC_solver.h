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

class StanceNMPC {
public:
    explicit StanceNMPC(const std::string& config_file);
    ~StanceNMPC();

    bool readYaml(const std::string& config_file);
    void initialize();
    void controlLoop();
    void drawPlot();

private:
    Eigen::Matrix<double, 12, 1> Q_read;
    Eigen::Matrix<double, 12, 1> QN_read;
    Eigen::Matrix<double, 3, 1> L_read;

    double Qu_read;
    double mu_read;
    double theta_read;
    double alpha_ocp_read;
    int loop_read;

    int Nrep; // 测试次数
    int N;    // 预测步长
    double T_step_read;
    bool show_flag_read;

    Eigen::Matrix<double, 12, 12> Q;
    Eigen::Matrix<double, 12, 12> S;
    Eigen::Matrix<double, 12, 12> R;
    Eigen::Matrix<double, 12, 12> QN;

    FlowDynamic flow_dynamic_;
    Eigen::MatrixXd x_nmpc;
    Eigen::MatrixXd u_nmpc;

    Eigen::VectorXd x0;
    Eigen::Matrix<double, 12, 1> x_ref_k;
    Eigen::MatrixXd x_ref;

    hpipm::OcpQpIpmSolverSettings solver_settings;

    std::mutex control_mutex;

    // Helper methods for control calculations
    void solveNMPC();
    void applyControl(const Eigen::VectorXd& control);
    void updateState();
    // ... other helper methods as needed ...
};

#endif // STANCE_NMPC_H
