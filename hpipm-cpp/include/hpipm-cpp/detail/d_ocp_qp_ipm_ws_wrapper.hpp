#ifndef HPIPM_CPP_D_OCP_QP_IPM_WS_WRAPPER_HPP_
#define HPIPM_CPP_D_OCP_QP_IPM_WS_WRAPPER_HPP_

#include <memory>

extern "C" {
#include "hpipm_d_ocp_qp_ipm.h"
}

#include "hpipm-cpp/detail/d_ocp_qp_dim_wrapper.hpp"
#include "hpipm-cpp/detail/d_ocp_qp_ipm_arg_wrapper.hpp"


namespace hpipm {

///
/// @class d_ocp_qp_ipm_ws_wrapper
/// @brief A wrapper of d_ocp_qp_ipm_ws with a memory management.
///
class d_ocp_qp_ipm_ws_wrapper {
public:
  ///
  /// @brief Constructor. Allocates the hpipm resource.
  /// @param[in] dim Dimension.
  /// @param[in] ipm_arg Ipm solver argument.
  ///
  d_ocp_qp_ipm_ws_wrapper(const std::shared_ptr<d_ocp_qp_dim_wrapper>& dim, 
                          const std::shared_ptr<d_ocp_qp_ipm_arg_wrapper>& ipm_arg);

  ///
  /// @brief Default constructor. Does not allocate the hpipm resource.
  ///
  d_ocp_qp_ipm_ws_wrapper();

  ///
  /// @brief Destructor.
  ///
  ~d_ocp_qp_ipm_ws_wrapper();

  ///
  /// @brief Prohibit copy constructor.
  ///
  d_ocp_qp_ipm_ws_wrapper(const d_ocp_qp_ipm_ws_wrapper&) = delete;

  ///
  /// @brief Prohibit copy assign operator.
  ///
  d_ocp_qp_ipm_ws_wrapper& operator=(const d_ocp_qp_ipm_ws_wrapper&) = delete;

  ///
  /// @brief Custom move constructor.
  ///
  d_ocp_qp_ipm_ws_wrapper(d_ocp_qp_ipm_ws_wrapper&&) noexcept;

  ///
  /// @brief Custom move assign operator.
  ///
  d_ocp_qp_ipm_ws_wrapper& operator=(d_ocp_qp_ipm_ws_wrapper&&) noexcept;

  ///
  /// @brief Gets the pointer to the hpipm resource. Throw an exception if the 
  /// memory for the instance is not allocated.
  /// @return Pointer to the hpipm resource.
  ///
  d_ocp_qp_ipm_ws* get();

  ///
  /// @brief Gets the const pointer to the hpipm instance.
  /// @return const pointer to the hpipm resource.
  ///
  const d_ocp_qp_ipm_ws* get() const;

  ///
  /// @brief Resizes the hpipm resource.
  /// @param[in] dim Dimension.
  /// @param[in] ipm_arg Ipm solver argument.
  ///
  void resize(const std::shared_ptr<d_ocp_qp_dim_wrapper>& dim,
              const std::shared_ptr<d_ocp_qp_ipm_arg_wrapper>& ipm_arg);

private:
  std::shared_ptr<d_ocp_qp_dim_wrapper> dim_;
  std::shared_ptr<d_ocp_qp_ipm_arg_wrapper> ipm_arg_;
  d_ocp_qp_ipm_ws ocp_qp_ipm_ws_hpipm_;
  void *memory_ = nullptr;
  hpipm_size_t memsize_ = 0;

  void copy(const d_ocp_qp_ipm_ws_wrapper& other);
};

} // namespace hpipm

#endif // HPIPM_CPP_D_OCP_QP_IPM_WS_WRAPPER_HPP_