#include "NMPC_solver.h"

NMPCSolver::NMPCSolver(const std::string &config_file)
{
  if (!readYaml(config_file))
  {
    throw std::runtime_error("Failed to read configuration file.");
  }
  initialize();
}

NMPCSolver::~NMPCSolver()
{
  // Cleanup if required
}

// Read configuration from YAML file
bool NMPCSolver::readYaml(const std::string &config_file)
{
  YAML::Node config;
  try
  {
    config = YAML::LoadFile("../config/mpc_option.yaml");

    std::vector<double> Q_values = config["MPC"]["Q"].as<std::vector<double>>();
    Q_read_ = Eigen::Map<Vec12<double>>(Q_values.data());

    std::vector<double> Qf_values = config["MPC"]["Qf"].as<std::vector<double>>();
    Qf_read_ = Eigen::Map<Vec12<double>>(Qf_values.data());

    R_read_ = config["MPC"]["R"].as<double>();
    N_ = config["MPC"]["horizon_MPC"].as<int>();
    dt_MPC_ = config["MPC"]["dt_MPC"].as<double>();
    sqp_max_loop_ = config["MPC"]["sqp_max_loop"].as<int>();

    std::vector<double> L_values = config["Physical"]["Lbody"].as<std::vector<double>>();
    L_ = Eigen::Map<Vec3<double>>(L_values.data());

    mu_barrier_ = config["mu_b"].as<double>();
    theta_barrier_ = config["theta_b"].as<double>();

    N_test_rep_ = config["N_rep"].as<int>();
  }
  catch (YAML::BadFile &e)
  {
    std::cout << "read error!" << std::endl;
    return false;
  }
  return true;
}

// Initialize NMPC parameters and matrices
void NMPCSolver::initialize()
{
  S_.setZero();
  Q_ = Q_read_.asDiagonal();
  R_ = R_read_ * Eigen::Matrix<double, 12, 12>::Identity();
  Qf_ = static_cast<double>(N_) * Qf_read_.asDiagonal();

  x_nmpc_.resize(12, N_ + 1);
  x_nmpc_.setZero();
  u_nmpc_.resize(12, N_);
  u_nmpc_.setOnes();
  u_nmpc_ *= 100.0;

  x0_.resize(12, 1);
  x_ref_.resize(12, N_ + 1);

  // hpipm params
  solver_settings_.mode = hpipm::HpipmMode::Speed;
  solver_settings_.iter_max = 30;
  solver_settings_.alpha_min = 1e-8;
  solver_settings_.mu0 = 1e2;
  solver_settings_.tol_stat = 1e-04;
  solver_settings_.tol_eq = 1e-04;
  solver_settings_.tol_ineq = 1e-04;
  solver_settings_.tol_comp = 1e-04;
  solver_settings_.reg_prim = 1e-12;
  solver_settings_.warm_start = 0;
  solver_settings_.pred_corr = 1;
  solver_settings_.ric_alg = 0;
  solver_settings_.split_step = 1;

  x_solved_.resize(12, N_ + 1);
  x_solved_.setZero();
  u_solved_.resize(12, N_);
  u_solved_.setZero();

  x_mpc_.resize(12);
  x_mpc_next_.resize(12);
  u_mpc_.resize(12);
  f_constrain_all_.resize(24, N_);
  f_all_.resize(12, N_);
  b_constrain_.resize(24);
  db_constrain_.resize(24);
  ddb_constrain_.resize(24);

  x_solved_.resize(12, N_ + 1);
  x_solved_.setZero();
  u_solved_.resize(12, N_);
  u_solved_.setZero();

  Jphi_x_ = Eigen::MatrixXd::Zero(12, N_ + 1);
  Jphi_u_ = Eigen::MatrixXd::Zero(12, N_);

  x_search_ = Eigen::VectorXd::Zero(12);
  x_search_next_ = Eigen::VectorXd::Zero(12);
  u_search_ = Eigen::VectorXd::Zero(12);
  f_search_ = Eigen::MatrixXd::Zero(12, 1);
  f_search_v_ = Eigen::VectorXd::Zero(12);
}

void NMPCSolver::printOptimizationInfo()
{

  // std::cout << "Testing repetitions: " << N_test_rep_ << std::endl;
  // std::cout << "Solution state = " << std::endl
  //           << x_nmpc_ << std::endl; // Output the solved state variables
  // std::cout << "Solution torque = " << std::endl
  //           << u_nmpc_ << std::endl; // Output the solved torque
  // std::cout << "Maximum friction cone constraint violation (negative value) = " << std::endl
  //           << f_constrain_all_.minCoeff() << std::endl; //
  // std::cout << "Maximum dynamic equation violation (positive value) = " << std::endl
  //           << f_all_.cwiseAbs().maxCoeff() << std::endl;
  // std::cout << "Dynamic equation violation = " << std::endl
  //           << f_all_ << std::endl;

  // std::cout << "sqp_loop : " << sqp_loop << std::endl;
  // std::cout << "phi = " << std::endl
  //           << phi_ << std::endl;
  // std::cout << "dphi = " << std::endl
  //           << dphi_ << std::endl;
  // std::cout << "theta = " << std::endl
  //           << theta_ << std::endl;
  // std::cout << "alpha = " << std::endl
  //           << alpha_ << std::endl;
  std::cout << "-----------------------" << std::endl;
  std::cout << "Testing repetitions: " << N_test_rep_ << std::endl;
  std::cout << "NMPC horizon: " << N_ << std::endl;
  std::cout << "NMPC dt: " << dt_MPC_ << std::endl;
}

bool NMPCSolver::checkConvergence()
{
  return linearSearch();
}

// reference: Perceptive Locomotion through Nonlinear Model Predictive Control
bool NMPCSolver::linearSearch()
{

  theta_ = 0.0; // equality constraints
  phi_ = 0.0;   // cost

  auto x_alpha = x_nmpc_;
  auto u_alpha = u_nmpc_;

  dphi_ = 0.0; // cal d_{cost}/d_{alpha_}

  for (int k_l = 0; k_l < N_ + 1; ++k_l)
  {
    x_search_ = x_nmpc_.block<12, 1>(0, k_l);
    u_search_ = u_nmpc_.block<12, 1>(0, k_l);
    if (k_l == N_)
    {
      phi_ += 0.5 * (x_search_ - x_ref_.block<12, 1>(0, k_l)).transpose() * Qf_ * (x_search_ - x_ref_.block<12, 1>(0, k_l));
      Jphi_x_.block<12, 1>(0, k_l) = Qf_ * (x_search_ - x_ref_.block<12, 1>(0, k_l));
    }
    else
    {
      x_search_next_ = x_nmpc_.block<12, 1>(0, k_l + 1);
      // cal equality constraints
      flow_dynamic_.GetShootingDynamic(x_search_, x_search_next_, u_search_, nullptr, nullptr, nullptr, &f_search_);
      f_search_v_ = f_search_;
      theta_ += 0.5 * (f_search_v_.transpose()) * f_search_v_;
      // cal x cost
      phi_ += 0.5 * (x_search_ - x_ref_.block<12, 1>(0, k_l)).transpose() * Q_ * (x_search_ - x_ref_.block<12, 1>(0, k_l));
      Jphi_x_.block<12, 1>(0, k_l) = Q_ * (x_search_ - x_ref_.block<12, 1>(0, k_l));
      // cal u cost
      flow_dynamic_.GetConstrain(u_search_, A_constrain_, f_constrain_); // get constraint equations
      for (int k_b = 0; k_b < 24; ++k_b)
      {
        // Convert constraint equations into barrier equations
        flow_dynamic_.Barrier(f_constrain_(k_b), mu_barrier_, theta_barrier_, &b_constrain_(k_b), &db_constrain_(k_b), &ddb_constrain_(k_b));
      }
      phi_ += b_constrain_.sum() + 0.5 * u_search_.transpose() * R_ * u_search_;
      Jphi_u_.block<12, 1>(0, k_l) = A_constrain_.transpose() * db_constrain_ + R_ * u_search_;
    }
  }

  for (int k_l = 0; k_l < N_ + 1; ++k_l)
  {
    dphi_ += x_solved_.block<12, 1>(0, k_l).transpose() * Jphi_x_.block<12, 1>(0, k_l);
    if (k_l < N_)
    {
      dphi_ += u_solved_.block<12, 1>(0, k_l).transpose() * Jphi_u_.block<12, 1>(0, k_l);
    }
  }

  while (alpha_ > alpha_min_)
  {
    double theta_alpha = 0.0;
    double phi_alpha = 0.0;
    x_alpha = x_nmpc_ + alpha_ * x_solved_;
    u_alpha = u_nmpc_ + alpha_ * u_solved_;
    //
    for (int k_l = 0; k_l < N_ + 1; ++k_l)
    {
      x_search_ = x_alpha.block<12, 1>(0, k_l);
      u_search_ = u_alpha.block<12, 1>(0, k_l);
      if (k_l == N_)
      {
        phi_alpha += 0.5 * (x_search_ - x_ref_.block<12, 1>(0, k_l)).transpose() * Qf_ * (x_search_ - x_ref_.block<12, 1>(0, k_l));
      }
      else
      {
        x_search_next_ = x_alpha.block<12, 1>(0, k_l + 1);
        // cal equality constraints
        flow_dynamic_.GetShootingDynamic(x_search_, x_search_next_, u_search_, nullptr, nullptr, nullptr, &f_search_);
        f_search_v_ = f_search_;
        theta_alpha += 0.5 * (f_search_v_.transpose()) * f_search_v_;
        // cal x cost
        phi_alpha += 0.5 * (x_search_ - x_ref_.block<12, 1>(0, k_l)).transpose() * Q_ * (x_search_ - x_ref_.block<12, 1>(0, k_l));
        // cal u cost
        flow_dynamic_.GetConstrain(u_search_, A_constrain_, f_constrain_); // get constraint equations
        for (int k_b = 0; k_b < 24; ++k_b)
        {
          // Convert constraint equations into barrier equations
          flow_dynamic_.Barrier(f_constrain_(k_b), mu_barrier_, theta_barrier_, &b_constrain_(k_b), &db_constrain_(k_b), &ddb_constrain_(k_b));
        }
        phi_alpha += b_constrain_.sum() + 0.5 * u_search_.transpose() * R_ * u_search_;
      }
    }

    if (theta_alpha > theta_max_)
    {
      if (theta_alpha < (1.0 - byta_theta_) * theta_)
      {
        x_nmpc_ = x_alpha;
        u_nmpc_ = u_alpha;
        break;
      }
    }
    else if ((std::max(theta_alpha, theta_) < theta_min_) and (dphi_ < 0.0))
    {
      if (phi_alpha < phi_ + eta_ * alpha_ * dphi_)
      {
        x_nmpc_ = x_alpha;
        u_nmpc_ = u_alpha;
        break;
      }
    }
    else
    {
      if ((phi_alpha < phi_ - byta_phi_ * theta_) or (theta_alpha < (1.0 - byta_theta_) * theta_))
      {
        x_nmpc_ = x_alpha;
        u_nmpc_ = u_alpha;
        break;
      }
    }
    alpha_ = byta_alpha_ * alpha_;
    //
  }

  // 优化完成
  if (dphi_ > -1e-3 and theta_ < 1e-6)
  {
    std::cout << "nmpc solve success!" << std::endl;
    // break;
    return true;
  }
  return false;
}

void NMPCSolver::prepareQpStructures(std::vector<hpipm::OcpQp> &qp)
{
  for (int i = 0; i < N_; ++i)
  {
    Eigen::MatrixXd A_dynamic, B_dynamic, b_dynamic;

    x_mpc_ = x_nmpc_.block<12, 1>(0, i);
    x_mpc_next_ = x_nmpc_.block<12, 1>(0, i + 1);
    u_mpc_ = u_nmpc_.block<12, 1>(0, i);
    // Obtain discrete dynamic matrices for shooting method
    flow_dynamic_.GetShootingDynamic(x_mpc_, x_mpc_next_, u_mpc_, &A_dynamic, &B_dynamic, &b_dynamic, nullptr);
    f_all_.block<12, 1>(0, i) = -b_dynamic;
    flow_dynamic_.GetConstrain(u_mpc_, A_constrain_, f_constrain_);
    f_constrain_all_.block<24, 1>(0, i) = f_constrain_;
    for (int k = 0; k < 24; ++k)
    {
      // Convert constraint equations into barrier equations
      flow_dynamic_.Barrier(f_constrain_(k), mu_barrier_, theta_barrier_, &b_constrain_(k), &db_constrain_(k), &ddb_constrain_(k));
    }

    // Set QP structures for the current state
    qp[i].A = A_dynamic;
    qp[i].B = B_dynamic;
    qp[i].b = b_dynamic;
    // Hard constraints are only enabled for strict inequality constraints
    // qp[i].C = Eigen::Matrix<double, 20, 12>::Zero();
    // qp[i].D = f_all_;
    // qp[i].ug = ub;
    // qp[i].lg = lb;
    qp[i].Q = Q_;
    qp[i].q = Q_ * (x_mpc_ - x_ref_.block<12, 1>(0, i));
    qp[i].S = S_;
    qp[i].R = R_ + A_constrain_.transpose() * ddb_constrain_.asDiagonal() * A_constrain_;
    qp[i].r = R_ * u_mpc_ + A_constrain_.transpose() * db_constrain_; // eb = diag(db)*Ja
  }
  // Set terminal cost for QP
  qp[N_].Q = Qf_;
  qp[N_].q = Qf_ * (x_mpc_next_ - x_ref_.block<12, 1>(0, N_));
}

void NMPCSolver::solveQpProblems(std::vector<hpipm::OcpQp> &qp, std::vector<hpipm::OcpQpSolution> &solution)
{
  // Solve QP problems using the HPIPM solver
  hpipm::OcpQpIpmSolver solver(qp, solver_settings_);
  const auto res = solver.solve(x0_ - x_nmpc_.block<12, 1>(0, 0), qp, solution);

  for (int i = 0; i < N_; ++i)
  {
    // Extract the solved state and control from the solution
    x_solved_.block<12, 1>(0, i) = solution[i].x;
    u_solved_.block<12, 1>(0, i) = solution[i].u;
  }
  // Extract the solved terminal state
  x_solved_.block<12, 1>(0, N_) = solution[N_].x;
}

void NMPCSolver::setupDynamics()
{
  flow_dynamic_.SetMass(15.0);               // Set mass
  flow_dynamic_.SetMPCdt(dt_MPC_);           // Set MPC time step
  flow_dynamic_.SetInertia(L_.asDiagonal()); // Set inertia matrix
  flow_dynamic_.SetFoot(Eigen::Vector3d(0.0, -0.1, 0.0), Eigen::Vector3d(0.0, 0.1, 0.0),
                        Eigen::Vector3d(1, 1, 1).asDiagonal(), Eigen::Vector3d(1, 1, 1).asDiagonal());
}

void NMPCSolver::setupReference()
{
  // [axis angle, angular velocity, CoM position, CoM velocity]
  x0_ << 0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0;
  x_ref_k_ << 0, 0, 0.2, 0, 0, 0, 0.5, 0, 1.0, 0, 0, 0;
  for (int k_ref = 0; k_ref < N_ + 1; ++k_ref)
  {
    // Propagate the reference state across the horizon
    x_ref_.block<12, 1>(0, k_ref) = x_ref_k_;
  }
}
// Main control loop
void NMPCSolver::controlLoop()
{
  timer timer_test;
  timer_test.start();

  // Testing average speed
  for (int nrep = 0; nrep < N_test_rep_; ++nrep)
  {
    setupDynamics();
    setupReference();

    std::vector<hpipm::OcpQp> qp(N_ + 1);
    std::vector<hpipm::OcpQpSolution> solution(N_ + 1);
    // SQP loop for convergence
    for (int i = 0; i < sqp_max_loop_; ++i)
    {
      prepareQpStructures(qp);
      solveQpProblems(qp, solution); // Solve QP problems
      if (checkConvergence())
      {
        break;
      }
    }
  }
  printOptimizationInfo();
  std::cout << "Average NMPC solution time = " << timer_test.get() / double(N_test_rep_) << "ms" << std::endl;
  // return ;
}
int main()
{
  try
  {
    NMPCSolver nmpc("../mpc_option.yaml");
    nmpc.controlLoop();
  }
  catch (const std::exception &e)
  {
    std::cerr << "An exception occurred: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
