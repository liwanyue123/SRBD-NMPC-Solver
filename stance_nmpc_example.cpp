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

// void print(const Eigen::MatrixXf& H, const std::string& s = "") {
//   std::cout << s << std::endl << H << std::endl;
// }

int main()
{

  timer timer_test;
  timer_test.start();

  Eigen::Matrix<double, 12, 1> Q_read;
  Eigen::Matrix<double, 12, 1> QN_read;
  Eigen::Matrix<double, 3, 1> L_read;
  double Qu_read = 0.0;
  double mu_read = 0.0;
  double theta_read = 0.0;
  double alpha_ocp_read = 1.0;
  int loop_read = 0;
  bool show_flag_read = false;

  int N_rep_read = 0;
  int N_step_read = 0;
  double T_step_read = 0.01;

  YAML::Node config;
  try
  {
    config = YAML::LoadFile("../mpc_option.yaml");

    Q_read[0] = config["K_r"][0].as<double>();
    Q_read[1] = config["K_r"][1].as<double>();
    Q_read[2] = config["K_r"][2].as<double>();
    Q_read[3] = config["K_w"][0].as<double>();
    Q_read[4] = config["K_w"][1].as<double>();
    Q_read[5] = config["K_w"][2].as<double>();
    Q_read[6] = config["K_p"][0].as<double>();
    Q_read[7] = config["K_p"][1].as<double>();
    Q_read[8] = config["K_p"][2].as<double>();
    Q_read[9] = config["K_v"][0].as<double>();
    Q_read[10] = config["K_v"][1].as<double>();
    Q_read[11] = config["K_v"][2].as<double>();

    QN_read[0] = config["K_r_N"][0].as<double>();
    QN_read[1] = config["K_r_N"][1].as<double>();
    QN_read[2] = config["K_r_N"][2].as<double>();
    QN_read[3] = config["K_w_N"][0].as<double>();
    QN_read[4] = config["K_w_N"][1].as<double>();
    QN_read[5] = config["K_w_N"][2].as<double>();
    QN_read[6] = config["K_p_N"][0].as<double>();
    QN_read[7] = config["K_p_N"][1].as<double>();
    QN_read[8] = config["K_p_N"][2].as<double>();
    QN_read[9] = config["K_v_N"][0].as<double>();
    QN_read[10] = config["K_v_N"][1].as<double>();
    QN_read[11] = config["K_v_N"][2].as<double>();

    L_read[0] = config["Lbody"][0].as<double>();
    L_read[1] = config["Lbody"][1].as<double>();
    L_read[2] = config["Lbody"][2].as<double>();

    Qu_read = config["K_u"].as<double>();
    mu_read = config["mu_b"].as<double>();
    theta_read = config["theta_b"].as<double>();
    alpha_ocp_read = config["alpha_ocp"].as<double>();
    loop_read = config["sqp_loop_time"].as<int>();

    N_rep_read = config["N_rep"].as<int>();
    N_step_read = config["N_step"].as<int>();
    T_step_read = config["T_step"].as<double>();

    show_flag_read = config["show_flag"].as<bool>();
  }
  catch (YAML::BadFile &e)
  {
    std::cout << "read error!" << std::endl;
    return -1;
  }

  int Nrep = N_rep_read;     // 测试100次
  const int N = N_step_read; // 预测步

  FlowDynamic flow_dynamic_;

  Eigen::MatrixXd x_nmpc;
  Eigen::MatrixXd u_nmpc;
  x_nmpc.resize(12, N + 1);
  u_nmpc.resize(12, N);
  x_nmpc.setZero();
  u_nmpc.setZero();

  for (int nrep = 0; nrep < Nrep; ++nrep)
  {
    flow_dynamic_.SetM(15.0);                                                                                                                                             // 设置质量
    flow_dynamic_.SetT(T_step_read);                                                                                                                                      // 求解间隔
    flow_dynamic_.SetL(L_read.asDiagonal());                                                                                                                              // 惯量的逆
    flow_dynamic_.SetFoot(Eigen::Vector3d(0.0, -0.1, 0.0), Eigen::Vector3d(0.0, 0.1, 0.0), Eigen::Vector3d(1, 1, 1).asDiagonal(), Eigen::Vector3d(1, 1, 1).asDiagonal()); // 足端位置

    // double dt = 0.02;
    // double g = 9.8;
    // double z = 0.6;
    // double l_zmp = 0.05;

    Eigen::VectorXd x0;
    x0.resize(12, 1);
    x0 << 0, 0, 0 /*初始轴角*/, 0.0, 0, 0.0 /*初始角速度*/, -0.06, 0, 1.0 /*初始位置*/, 0, 0, 0 /*初始速度*/;

    // Eigen::MatrixXd x_nmpc;
    // Eigen::MatrixXd u_nmpc;
    // x_nmpc.resize(12, N+1);
    // u_nmpc.resize(12, N);
    x_nmpc.setZero();
    u_nmpc.setZero();
    u_nmpc.setOnes();
    u_nmpc *= 100.0;
    // Eigen::Matrix<double, 12, 12> S_u;
    // S_u.setZero();
    // S_u.diagonal() << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
    // u_nmpc = S_u*u_nmpc;

    Eigen::VectorXd x_mpc;
    Eigen::VectorXd x_mpc_next;
    Eigen::VectorXd u_mpc;
    x_mpc.resize(12);
    x_mpc_next.resize(12);
    u_mpc.resize(12);

    Eigen::Matrix<double, 12, 1> x_ref_k;
    x_ref_k.setZero();
    x_ref_k << 0.0, 0.0, 0.0 /*目标轴角*/, 0, 0, 0.0 /*目标角速度*/, 0.0, 0, 1.0 /*目标位置*/, 0, 0, 0 /*目标速度*/;
    Eigen::MatrixXd x_ref;
    x_ref.resize(12, N + 1);

    for (int k_ref = 0; k_ref < N + 1; ++k_ref)
    {
      x_ref.block<12, 1>(0, k_ref) = x_ref_k;
    }

    Eigen::Matrix<double, 12, 12> Q;
    Q.setZero();
    // Q.diagonal() << 10,10,10,0.001,0.001,0.001,500,500,500,1,1,1; // 设定状态权重
    Q = Q_read.asDiagonal();
    Eigen::Matrix<double, 12, 12> S;
    S.setZero();
    // Eigen::Matrix<double, 12, 12> R = 1e-4*Eigen::Matrix<double, 12, 12>::Identity(); // 力权重
    Eigen::Matrix<double, 12, 12> R = Qu_read * Eigen::Matrix<double, 12, 12>::Identity(); // 力权重

    Eigen::Matrix<double, 12, 12> QN;
    QN.setZero();
    // QN.diagonal() << 10*N,10*N,10*N,0.1*N,0.1*N,0.1*N,1000*N,1000*N,1000*N,10*N,10*N,10*N; // 设定状态权重
    // QN = 1000.0*Q;
    QN = static_cast<double>(N) * QN_read.asDiagonal();

    Eigen::MatrixXd A_dynamic, B_dynamic, b_dynamic, A_constrain, f_constrain;
    Eigen::MatrixXd f_constrain_all, f_all;

    f_constrain_all.resize(24, N);
    f_all.resize(12, N);

    double mu_barrier = mu_read;                              // 障碍方程系数
    double theta_barrier = theta_read;                        // 障碍方程系数
    Eigen::VectorXd b_constrain, db_constrain, ddb_constrain; // 障碍方程
    b_constrain.resize(24);
    db_constrain.resize(24);
    ddb_constrain.resize(24);

    // hpipm求解参数
    hpipm::OcpQpIpmSolverSettings solver_settings;
    solver_settings.mode = hpipm::HpipmMode::Speed;
    solver_settings.iter_max = 30;
    solver_settings.alpha_min = 1e-8;
    solver_settings.mu0 = 1e2;
    solver_settings.tol_stat = 1e-04;
    solver_settings.tol_eq = 1e-04;
    solver_settings.tol_ineq = 1e-04;
    solver_settings.tol_comp = 1e-04;
    solver_settings.reg_prim = 1e-12;
    solver_settings.warm_start = 0;
    solver_settings.pred_corr = 1;
    solver_settings.ric_alg = 0;
    solver_settings.split_step = 1;

    std::vector<hpipm::OcpQp> qp(N + 1);
    std::vector<hpipm::OcpQpSolution> solution(N + 1);

    for (int sqp_loop = 0; sqp_loop < loop_read; ++sqp_loop)
    {

      for (int i = 0; i < N; ++i)
      {
        x_mpc = x_nmpc.block<12, 1>(0, i);
        x_mpc_next = x_nmpc.block<12, 1>(0, i + 1);
        u_mpc = u_nmpc.block<12, 1>(0, i);
        flow_dynamic_.GetShootingDynamic(x_mpc, x_mpc_next, u_mpc, &A_dynamic, &B_dynamic, &b_dynamic, nullptr); // 获取离散动力学矩阵
        f_all.block<12, 1>(0, i) = -b_dynamic;
        flow_dynamic_.GetConstrain(u_mpc, A_constrain, f_constrain); // 获取约束方程
        f_constrain_all.block<24, 1>(0, i) = f_constrain;
        for (int k = 0; k < 24; ++k)
        {
          flow_dynamic_.Barrier(f_constrain(k), mu_barrier, theta_barrier, &b_constrain(k), &db_constrain(k), &ddb_constrain(k)); // 将约束方程转化为障碍方程
        }

        qp[i].A = A_dynamic;
        qp[i].B = B_dynamic;
        qp[i].b = b_dynamic;
        // qp[i].C = Eigen::Matrix<double, 20, 12>::Zero(); // 只有必须严格满足的不等式约束会启用此硬约束，开启后求解速度会大幅下降
        // qp[i].D = f_all;
        // qp[i].ug = ub;
        // qp[i].lg = lb;
        qp[i].Q = Q;
        qp[i].q = Q * (x_mpc - x_ref.block<12, 1>(0, i));
        qp[i].S = S;
        qp[i].R = R + A_constrain.transpose() * ddb_constrain.asDiagonal() * A_constrain;
        qp[i].r = R * u_mpc + A_constrain.transpose() * db_constrain; // eb = diag(db)*Ja
      }
      qp[N].Q = QN;
      qp[N].q = QN * (x_mpc_next - x_ref.block<12, 1>(0, N));

      hpipm::OcpQpIpmSolver solver(qp, solver_settings);
      const auto res = solver.solve(x0 - x_nmpc.block<12, 1>(0, 0), qp, solution); // 求解 qp问题

      Eigen::MatrixXd get_solveX;
      get_solveX.resize(12, N + 1);
      get_solveX.setZero();
      Eigen::MatrixXd get_solveU;
      get_solveU.resize(12, N);
      get_solveU.setZero();
      for (int i = 0; i <= N; ++i)
      {
        get_solveX.block<12, 1>(0, i) = solution[i].x;
      }

      for (int i = 0; i < N; ++i)
      {
        get_solveU.block<12, 1>(0, i) = solution[i].u;
      }

      // 线性搜索以提高收敛速度
      double theta = 0.0; // 等式约束
      double phi = 0.0;   // cost
      Eigen::MatrixXd Jphi_x = Eigen::MatrixXd::Zero(12, N + 1);
      Eigen::MatrixXd Jphi_u = Eigen::MatrixXd::Zero(12, N);

      Eigen::VectorXd x_search = Eigen::VectorXd::Zero(12);
      Eigen::VectorXd x_search_next = Eigen::VectorXd::Zero(12);
      Eigen::VectorXd u_search = Eigen::VectorXd::Zero(12);
      Eigen::MatrixXd f_search = Eigen::MatrixXd::Zero(12, 1);
      Eigen::VectorXd f_search_v = Eigen::VectorXd::Zero(12);

      for (int k_l = 0; k_l < N + 1; ++k_l)
      {
        x_search = x_nmpc.block<12, 1>(0, k_l);
        u_search = u_nmpc.block<12, 1>(0, k_l);
        if (k_l == N)
        {
          phi += 0.5 * (x_search - x_ref.block<12, 1>(0, k_l)).transpose() * QN * (x_search - x_ref.block<12, 1>(0, k_l));
          Jphi_x.block<12, 1>(0, k_l) = QN * (x_search - x_ref.block<12, 1>(0, k_l));
        }
        else
        {
          x_search_next = x_nmpc.block<12, 1>(0, k_l + 1);
          // 计算等式约束
          flow_dynamic_.GetShootingDynamic(x_search, x_search_next, u_search, nullptr, nullptr, nullptr, &f_search);
          f_search_v = f_search;
          theta += 0.5 * (f_search_v.transpose()) * f_search_v;
          // 计算x cost
          phi += 0.5 * (x_search - x_ref.block<12, 1>(0, k_l)).transpose() * Q * (x_search - x_ref.block<12, 1>(0, k_l));
          Jphi_x.block<12, 1>(0, k_l) = Q * (x_search - x_ref.block<12, 1>(0, k_l));
          // 计算 u cost
          flow_dynamic_.GetConstrain(u_search, A_constrain, f_constrain); // 获取约束方程
          for (int k_b = 0; k_b < 24; ++k_b)
          {
            flow_dynamic_.Barrier(f_constrain(k_b), mu_barrier, theta_barrier, &b_constrain(k_b), &db_constrain(k_b), &ddb_constrain(k_b)); // 将约束方程转化为障碍方程
          }
          phi += b_constrain.sum() + 0.5 * u_search.transpose() * R * u_search;
          Jphi_u.block<12, 1>(0, k_l) = A_constrain.transpose() * db_constrain + R * u_search;
        }
      }

      double theta_max = 1e-6;
      double theta_min = 5e-10;
      double eta = 1e-4;
      double byta_phi = 1e-6;
      double byta_theta = 1e-6;
      double byta_alpha = 0.5;
      double alpha_min = 1e-4;

      double alpha = 1.0;
      auto x_alpha = x_nmpc;
      auto u_alpha = u_nmpc;

      double dphi = 0.0; // 计算d_{cost}/d_{alpha}

      for (int k_l = 0; k_l < N + 1; ++k_l)
      {
        dphi += get_solveX.block<12, 1>(0, k_l).transpose() * Jphi_x.block<12, 1>(0, k_l);
        if (k_l < N)
        {
          dphi += get_solveU.block<12, 1>(0, k_l).transpose() * Jphi_u.block<12, 1>(0, k_l);
        }
      }

      while (alpha > alpha_min)
      {
        double theta_alpha = 0.0;
        double phi_alpha = 0.0;
        x_alpha = x_nmpc + alpha * get_solveX;
        u_alpha = u_nmpc + alpha * get_solveU;
        //
        for (int k_l = 0; k_l < N + 1; ++k_l)
        {
          x_search = x_alpha.block<12, 1>(0, k_l);
          u_search = u_alpha.block<12, 1>(0, k_l);
          if (k_l == N)
          {
            phi_alpha += 0.5 * (x_search - x_ref.block<12, 1>(0, k_l)).transpose() * QN * (x_search - x_ref.block<12, 1>(0, k_l));
          }
          else
          {
            x_search_next = x_alpha.block<12, 1>(0, k_l + 1);
            // 计算等式约束
            flow_dynamic_.GetShootingDynamic(x_search, x_search_next, u_search, nullptr, nullptr, nullptr, &f_search);
            f_search_v = f_search;
            theta_alpha += 0.5 * (f_search_v.transpose()) * f_search_v;
            // 计算x cost
            phi_alpha += 0.5 * (x_search - x_ref.block<12, 1>(0, k_l)).transpose() * Q * (x_search - x_ref.block<12, 1>(0, k_l));
            // 计算 u cost
            flow_dynamic_.GetConstrain(u_search, A_constrain, f_constrain); // 获取约束方程
            for (int k_b = 0; k_b < 24; ++k_b)
            {
              flow_dynamic_.Barrier(f_constrain(k_b), mu_barrier, theta_barrier, &b_constrain(k_b), &db_constrain(k_b), &ddb_constrain(k_b)); // 将约束方程转化为障碍方程
            }
            phi_alpha += b_constrain.sum() + 0.5 * u_search.transpose() * R * u_search;
          }
        }

        if (theta_alpha > theta_max)
        {
          if (theta_alpha < (1.0 - byta_theta) * theta)
          {
            x_nmpc = x_alpha;
            u_nmpc = u_alpha;
            break;
          }
        }
        else if ((std::max(theta_alpha, theta) < theta_min) and (dphi < 0.0))
        {
          if (phi_alpha < phi + eta * alpha * dphi)
          {
            x_nmpc = x_alpha;
            u_nmpc = u_alpha;
            break;
          }
        }
        else
        {
          if ((phi_alpha < phi - byta_phi * theta) or (theta_alpha < (1.0 - byta_theta) * theta))
          {
            x_nmpc = x_alpha;
            u_nmpc = u_alpha;
            break;
          }
        }
        alpha = byta_alpha * alpha;
        //
      }

      // TODO(2) : 轴角应使用SO3加法进行迭代，目前直接相加，有可能导致轴角圈数超过1
      // x_nmpc += 0.001*get_solveX;
      // u_nmpc += 0.001*get_solveU;

      std::cout << "i : " << sqp_loop << std::endl;
      std::cout << "phi = " << std::endl
                << phi << std::endl;
      std::cout << "dphi = " << std::endl
                << dphi << std::endl;
      std::cout << "theta = " << std::endl
                << theta << std::endl;
      std::cout << "alpha = " << std::endl
                << alpha << std::endl;
      // 优化完成
      if (dphi > -1e-3 and theta < 1e-6)
      {
        std::cout << "nmpc solve success!" << std::endl;
        break;
      }

      // std::cout << "x_nmpc = " << std::endl << x_nmpc << std::endl;
      // std::cout << "u_nmpc = " << std::endl << u_nmpc << std::endl;
      // std::cout << "f_constrain_all = " << std::endl << f_constrain_all << std::endl;
      // std::cout << "f_all = " << std::endl << f_all << std::endl;
    }
    std::cout << "迭代次数 : " << nrep << std::endl;
    std::cout << "求解状态 = " << std::endl
              << x_nmpc << std::endl; // 求解状态量
    std::cout << "求解力矩 = " << std::endl
              << u_nmpc << std::endl; // 求解力
    std::cout << "最大摩擦锥约束违反量（负数） = " << std::endl
              << f_constrain_all.minCoeff() << std::endl; //
    std::cout << "最大动力学方程违反量（正数） = " << std::endl
              << f_all.cwiseAbs().maxCoeff() << std::endl;
    std::cout << "动力学方程违反量 = " << std::endl
              << f_all << std::endl;
  }

  double box_x = 0.10;
  double box_y = 0.25;
  double box_z = 0.50;

  Eigen::Matrix<double, 3, 10> box_list_raw;
  Eigen::Matrix<double, 3, 6> box_list_;
  box_list_raw << Eigen::Matrix<double, 3, 1>(0.5 * box_x, 0.5 * box_y, 0.5 * box_z),
      Eigen::Matrix<double, 3, 1>(0.5 * box_x, 0.5 * box_y, -0.5 * box_z),
      Eigen::Matrix<double, 3, 1>(0.5 * box_x, -0.5 * box_y, -0.5 * box_z),
      Eigen::Matrix<double, 3, 1>(0.5 * box_x, -0.5 * box_y, 0.5 * box_z),
      Eigen::Matrix<double, 3, 1>(0.5 * box_x, 0.5 * box_y, 0.5 * box_z),
      Eigen::Matrix<double, 3, 1>(-0.5 * box_x, 0.5 * box_y, 0.5 * box_z),
      Eigen::Matrix<double, 3, 1>(-0.5 * box_x, 0.5 * box_y, -0.5 * box_z),
      Eigen::Matrix<double, 3, 1>(-0.5 * box_x, -0.5 * box_y, -0.5 * box_z),
      Eigen::Matrix<double, 3, 1>(-0.5 * box_x, -0.5 * box_y, 0.5 * box_z),
      Eigen::Matrix<double, 3, 1>(-0.5 * box_x, 0.5 * box_y, 0.5 * box_z);

  std::cout << "nmpc平均求解时间 = " << timer_test.get() / double(Nrep) << std::endl;

  if (show_flag_read)
  {
    for (int show_loop = 0; show_loop < N + 1; ++show_loop)
    {
      Eigen::Matrix<double, 3, 1> r_k = x_nmpc.block<3, 1>(0, show_loop);
      Eigen::Matrix<double, 3, 1> p_k = x_nmpc.block<3, 1>(6, show_loop);
      Eigen::Matrix<double, 3, 10> box_list = expm(r_k) * box_list_raw;
      box_list.colwise() += p_k;

      std::vector<double> x_box;
      std::vector<double> y_box;
      std::vector<double> z_box;

      x_box.resize(10);

      y_box.resize(10);

      z_box.resize(10);
      for (int k_box = 0; k_box < 10; ++k_box)
      {
        x_box[k_box] = box_list(0, k_box);
        y_box[k_box] = box_list(1, k_box);
        z_box[k_box] = box_list(2, k_box);
      }
      matplotlibcpp::figure(1);
      std::map<std::string, std::string> map_;
      matplotlibcpp::plot3(x_box, y_box, z_box, map_, 1);

      // d R*I*Rt*w = tau
      //  R*I*Rt*dw + [w]*R*I*Rt*w*dt = Tau*dt

      std::vector<double> x_sky = {0};
      std::vector<double> y_sky = {0};
      std::vector<double> z_sky = {1.0};

      std::vector<double> x_ground = {0};
      std::vector<double> y_ground = {0};
      std::vector<double> z_ground = {-0.05};

      matplotlibcpp::xlim(-0.5, 0.5);
      matplotlibcpp::ylim(-0.5, 0.5);

      // matplotlibcpp::text(0,0,2.0,"sky");

      matplotlibcpp::plot3(x_sky, y_sky, z_sky, map_, 1);
      matplotlibcpp::plot3(x_ground, y_ground, z_ground, map_, 1);
      matplotlibcpp::set_zlabel("sky");
      // matplotlibcpp::
      // matplotlibcpp::zlim(-1, 1);
      // matplotlibcpp::plot3({);
      // matplotlibcpp::show();

      matplotlibcpp::pause(0.1);
      matplotlibcpp::clf();
    }
  }
  // std::cout << QN_read << std::endl;

  Eigen::MatrixXd A_test, B_test, b_test, f_test;
  Eigen::VectorXd x_test;
  Eigen::VectorXd x_test_next;
  Eigen::VectorXd u_test;

  x_test.resize(12);
  x_test_next.resize(12);

  u_test.resize(12);

  x_test_next.setZero();

  double m = 15.0;
  double g = 9.8;
  x_test << 0.1, 0.15, 0.25, 0.1, 0.2, 0.3, 0.05, 0.06, 0.5, 0.5, 0.8, 0.1;
  u_test << 0, 0, 0.5 * m * g, 0.1, 0.2, 0.3, 0, 0, 0.5 * m * g, 0.4, 0.5, 0.6;
  flow_dynamic_.GetShootingDynamic(x_test, x_test_next, u_test, &A_test, &B_test, &b_test, &f_test); // 获取离散动力学矩阵

  // std::cout << "A_test = " << std::endl << A_test << std::endl; // 求解状态量
  // std::cout << "B_test = " << std::endl << B_test << std::endl; // 求解状态量

  // std::cout << "b_test = " << std::endl << b_test << std::endl; // 求解状态量
  // std::cout << "f_test = " << std::endl << f_test << std::endl; // 求解状态量

  /**/
  // Eigen::Matrix<double, 3, 1> val_test(0.1,0.2,0.3);

  // Eigen::Matrix<double, 3, 3> djl_x;
  // Eigen::Matrix<double, 3, 3> djl_y;
  // Eigen::Matrix<double, 3, 3> djl_z;
  // Eigen::Matrix<double, 3, 3> djlt_x;
  // Eigen::Matrix<double, 3, 3> djlt_y;
  // Eigen::Matrix<double, 3, 3> djlt_z;
  // djl<double>(val_test, &djl_x, &djl_y, &djl_z);
  // djlt<double>(val_test, &djlt_x, &djlt_y, &djlt_z);
  // /**/
  // std::cout << "djl_x = " << std::endl << djl_x << std::endl; // 求解状态量
  // std::cout << "djl_y = " << std::endl << djl_y << std::endl; // 求解状态量
  // std::cout << "djl_z = " << std::endl << djl_z << std::endl; // 求解状态量
  // std::cout << "djlt_x = " << std::endl << djlt_x << std::endl; // 求解状态量
  // std::cout << "djlt_y = " << std::endl << djlt_y << std::endl; // 求解状态量
  // std::cout << "djlt_z = " << std::endl << djlt_z << std::endl; // 求解状态量
  return 0;
}
