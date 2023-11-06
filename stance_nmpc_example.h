// #define EIGEN_USE_BLAS
// #include <ros/ros.h>

// #include <visualization_msgs/Marker.h>
// #include <visualization_msgs/MarkerArray.h>

#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <string>
#include <iostream>
#include "timer/timer.h"
// #include "hpllt.h"
#include "hpipm-cpp/include/hpipm-cpp/hpipm-cpp.hpp"
#include "dynamic/flow_tool.h"
#include "dynamic/flow_dynamic.h"
#include "matplotlibcpp.h"
#include <yaml-cpp/yaml.h>
// #include <cblas.h>
#include <thread>
#include <mutex>

Eigen::Matrix<double, 12, 1> Q_read;
Eigen::Matrix<double, 12, 1> QN_read;
Eigen::Matrix<double, 3, 1> L_read;

double Qu_read = 0.0;
double mu_read = 0.0;
double theta_read = 0.0;
double alpha_ocp_read = 1.0;
int loop_read = 0;

int Nrep; // 测试100次
int N;    // 预测步
double T_step_read = 0.01;
bool show_flag_read = false;

// control param
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

// hpipm求解参数
hpipm::OcpQpIpmSolverSettings solver_settings;