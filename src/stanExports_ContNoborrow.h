// Generated by rstantools.  Do not edit by hand.

/*
    psBayesborrow is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    psBayesborrow is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with psBayesborrow.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.26.1-4-gd72b68b7-dirty
#include <stan/model/model_header.hpp>
namespace model_ContNoborrow_namespace {
inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math;
using stan::math::pow; 
stan::math::profile_map profiles__;
static int current_statement__= 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in 'string', line 11, column 2 to column 13)",
                                                      " (in 'string', line 12, column 2 to column 15)",
                                                      " (in 'string', line 13, column 2 to column 17)",
                                                      " (in 'string', line 14, column 2 to column 22)",
                                                      " (in 'string', line 18, column 4 to column 53)",
                                                      " (in 'string', line 17, column 2 to line 18, column 53)",
                                                      " (in 'string', line 20, column 4 to column 53)",
                                                      " (in 'string', line 19, column 2 to line 20, column 53)",
                                                      " (in 'string', line 2, column 2 to column 19)",
                                                      " (in 'string', line 3, column 2 to column 19)",
                                                      " (in 'string', line 4, column 2 to column 17)",
                                                      " (in 'string', line 5, column 9 to column 12)",
                                                      " (in 'string', line 5, column 2 to column 18)",
                                                      " (in 'string', line 6, column 9 to column 12)",
                                                      " (in 'string', line 6, column 2 to column 18)",
                                                      " (in 'string', line 7, column 20 to column 23)",
                                                      " (in 'string', line 7, column 13 to column 14)",
                                                      " (in 'string', line 7, column 2 to column 25)",
                                                      " (in 'string', line 8, column 20 to column 23)",
                                                      " (in 'string', line 8, column 13 to column 14)",
                                                      " (in 'string', line 8, column 2 to column 25)",
                                                      " (in 'string', line 13, column 9 to column 10)"};
#include <stan_meta_header.hpp>
class model_ContNoborrow final : public model_base_crtp<model_ContNoborrow> {
private:
  int nCT;
  int nCC;
  int p;
  Eigen::Matrix<double, -1, 1> yCT;
  Eigen::Matrix<double, -1, 1> yCC;
  std::vector<Eigen::Matrix<double, 1, -1>> xCT;
  std::vector<Eigen::Matrix<double, 1, -1>> xCC;
 
public:
  ~model_ContNoborrow() { }
  
  inline std::string model_name() const final { return "model_ContNoborrow"; }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.26.1-4-gd72b68b7-dirty", "stancflags = "};
  }
  
  
  model_ContNoborrow(stan::io::var_context& context__,
                     unsigned int random_seed__ = 0,
                     std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "model_ContNoborrow_namespace::model_ContNoborrow";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 9;
      context__.validate_dims("data initialization","nCT","int",
          context__.to_vec());
      nCT = std::numeric_limits<int>::min();
      
      current_statement__ = 9;
      nCT = context__.vals_i("nCT")[(1 - 1)];
      current_statement__ = 9;
      current_statement__ = 9;
      check_greater_or_equal(function__, "nCT", nCT, 0);
      current_statement__ = 10;
      context__.validate_dims("data initialization","nCC","int",
          context__.to_vec());
      nCC = std::numeric_limits<int>::min();
      
      current_statement__ = 10;
      nCC = context__.vals_i("nCC")[(1 - 1)];
      current_statement__ = 10;
      current_statement__ = 10;
      check_greater_or_equal(function__, "nCC", nCC, 0);
      current_statement__ = 11;
      context__.validate_dims("data initialization","p","int",
          context__.to_vec());
      p = std::numeric_limits<int>::min();
      
      current_statement__ = 11;
      p = context__.vals_i("p")[(1 - 1)];
      current_statement__ = 11;
      current_statement__ = 11;
      check_greater_or_equal(function__, "p", p, 0);
      current_statement__ = 12;
      validate_non_negative_index("yCT", "nCT", nCT);
      current_statement__ = 13;
      context__.validate_dims("data initialization","yCT","double",
          context__.to_vec(nCT));
      yCT = Eigen::Matrix<double, -1, 1>(nCT);
      stan::math::fill(yCT, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> yCT_flat__;
        current_statement__ = 13;
        assign(yCT_flat__, nil_index_list(), context__.vals_r("yCT"),
          "assigning variable yCT_flat__");
        current_statement__ = 13;
        pos__ = 1;
        current_statement__ = 13;
        for (int sym1__ = 1; sym1__ <= nCT; ++sym1__) {
          current_statement__ = 13;
          assign(yCT, cons_list(index_uni(sym1__), nil_index_list()),
            yCT_flat__[(pos__ - 1)], "assigning variable yCT");
          current_statement__ = 13;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 14;
      validate_non_negative_index("yCC", "nCC", nCC);
      current_statement__ = 15;
      context__.validate_dims("data initialization","yCC","double",
          context__.to_vec(nCC));
      yCC = Eigen::Matrix<double, -1, 1>(nCC);
      stan::math::fill(yCC, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> yCC_flat__;
        current_statement__ = 15;
        assign(yCC_flat__, nil_index_list(), context__.vals_r("yCC"),
          "assigning variable yCC_flat__");
        current_statement__ = 15;
        pos__ = 1;
        current_statement__ = 15;
        for (int sym1__ = 1; sym1__ <= nCC; ++sym1__) {
          current_statement__ = 15;
          assign(yCC, cons_list(index_uni(sym1__), nil_index_list()),
            yCC_flat__[(pos__ - 1)], "assigning variable yCC");
          current_statement__ = 15;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 16;
      validate_non_negative_index("xCT", "nCT", nCT);
      current_statement__ = 17;
      validate_non_negative_index("xCT", "p", p);
      current_statement__ = 18;
      context__.validate_dims("data initialization","xCT","double",
          context__.to_vec(nCT, p));
      xCT = std::vector<Eigen::Matrix<double, 1, -1>>(nCT, Eigen::Matrix<double, 1, -1>(p));
      stan::math::fill(xCT, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> xCT_flat__;
        current_statement__ = 18;
        assign(xCT_flat__, nil_index_list(), context__.vals_r("xCT"),
          "assigning variable xCT_flat__");
        current_statement__ = 18;
        pos__ = 1;
        current_statement__ = 18;
        for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
          current_statement__ = 18;
          for (int sym2__ = 1; sym2__ <= nCT; ++sym2__) {
            current_statement__ = 18;
            assign(xCT,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              xCT_flat__[(pos__ - 1)], "assigning variable xCT");
            current_statement__ = 18;
            pos__ = (pos__ + 1);}}
      }
      current_statement__ = 19;
      validate_non_negative_index("xCC", "nCC", nCC);
      current_statement__ = 20;
      validate_non_negative_index("xCC", "p", p);
      current_statement__ = 21;
      context__.validate_dims("data initialization","xCC","double",
          context__.to_vec(nCC, p));
      xCC = std::vector<Eigen::Matrix<double, 1, -1>>(nCC, Eigen::Matrix<double, 1, -1>(p));
      stan::math::fill(xCC, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> xCC_flat__;
        current_statement__ = 21;
        assign(xCC_flat__, nil_index_list(), context__.vals_r("xCC"),
          "assigning variable xCC_flat__");
        current_statement__ = 21;
        pos__ = 1;
        current_statement__ = 21;
        for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
          current_statement__ = 21;
          for (int sym2__ = 1; sym2__ <= nCC; ++sym2__) {
            current_statement__ = 21;
            assign(xCC,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              xCC_flat__[(pos__ - 1)], "assigning variable xCC");
            current_statement__ = 21;
            pos__ = (pos__ + 1);}}
      }
      current_statement__ = 22;
      validate_non_negative_index("beta", "p", p);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += p;
      num_params_r__ += 1;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI, stan::require_vector_like_t<VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "model_ContNoborrow_namespace::log_prob";
(void) function__;  // suppress unused var warning
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      local_scalar_t__ theta;
      theta = DUMMY_VAR__;
      
      current_statement__ = 1;
      theta = in__.scalar();
      local_scalar_t__ gammaCC;
      gammaCC = DUMMY_VAR__;
      
      current_statement__ = 2;
      gammaCC = in__.scalar();
      Eigen::Matrix<local_scalar_t__, -1, 1> beta;
      beta = Eigen::Matrix<local_scalar_t__, -1, 1>(p);
      stan::math::fill(beta, DUMMY_VAR__);
      
      current_statement__ = 3;
      beta = in__.vector(p);
      local_scalar_t__ sigma;
      sigma = DUMMY_VAR__;
      
      current_statement__ = 4;
      sigma = in__.scalar();
      current_statement__ = 4;
      if (jacobian__) {
        current_statement__ = 4;
        sigma = stan::math::lb_constrain(sigma, 0, lp__);
      } else {
        current_statement__ = 4;
        sigma = stan::math::lb_constrain(sigma, 0);
      }
      {
        current_statement__ = 6;
        for (int i = 1; i <= nCT; ++i) {
          current_statement__ = 5;
          lp_accum__.add(
            normal_lpdf<propto__>(yCT[(i - 1)],
              ((gammaCC + theta) + multiply(xCT[(i - 1)], beta)), sigma));}
        current_statement__ = 8;
        for (int i = 1; i <= nCC; ++i) {
          current_statement__ = 7;
          lp_accum__.add(
            normal_lpdf<propto__>(yCC[(i - 1)],
              (gammaCC + multiply(xCC[(i - 1)], beta)), sigma));}
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr>
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "model_ContNoborrow_namespace::write_array";
(void) function__;  // suppress unused var warning
    (void) function__;  // suppress unused var warning
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      double theta;
      theta = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      theta = in__.scalar();
      double gammaCC;
      gammaCC = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      gammaCC = in__.scalar();
      Eigen::Matrix<double, -1, 1> beta;
      beta = Eigen::Matrix<double, -1, 1>(p);
      stan::math::fill(beta, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 3;
      beta = in__.vector(p);
      double sigma;
      sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      sigma = in__.scalar();
      current_statement__ = 4;
      sigma = stan::math::lb_constrain(sigma, 0);
      vars__.emplace_back(theta);
      vars__.emplace_back(gammaCC);
      for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
        vars__.emplace_back(beta[(sym1__ - 1)]);}
      vars__.emplace_back(sigma);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, stan::require_std_vector_t<VecVar>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void transform_inits_impl(const stan::io::var_context& context__,
                                   VecI& params_i__, VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.clear();
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      double theta;
      theta = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      theta = context__.vals_r("theta")[(1 - 1)];
      double gammaCC;
      gammaCC = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      gammaCC = context__.vals_r("gammaCC")[(1 - 1)];
      Eigen::Matrix<double, -1, 1> beta;
      beta = Eigen::Matrix<double, -1, 1>(p);
      stan::math::fill(beta, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> beta_flat__;
        current_statement__ = 3;
        assign(beta_flat__, nil_index_list(), context__.vals_r("beta"),
          "assigning variable beta_flat__");
        current_statement__ = 3;
        pos__ = 1;
        current_statement__ = 3;
        for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
          current_statement__ = 3;
          assign(beta, cons_list(index_uni(sym1__), nil_index_list()),
            beta_flat__[(pos__ - 1)], "assigning variable beta");
          current_statement__ = 3;
          pos__ = (pos__ + 1);}
      }
      double sigma;
      sigma = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      sigma = context__.vals_r("sigma")[(1 - 1)];
      double sigma_free__;
      sigma_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      sigma_free__ = stan::math::lb_free(sigma, 0);
      vars__.emplace_back(theta);
      vars__.emplace_back(gammaCC);
      for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
        vars__.emplace_back(beta[(sym1__ - 1)]);}
      vars__.emplace_back(sigma_free__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__.clear();
    names__.emplace_back("theta");
    names__.emplace_back("gammaCC");
    names__.emplace_back("beta");
    names__.emplace_back("sigma");
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.clear();
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(p)});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "theta");
    param_names__.emplace_back(std::string() + "gammaCC");
    for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "beta" + '.' + std::to_string(sym1__));
      }}
    param_names__.emplace_back(std::string() + "sigma");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "theta");
    param_names__.emplace_back(std::string() + "gammaCC");
    for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "beta" + '.' + std::to_string(sym1__));
      }}
    param_names__.emplace_back(std::string() + "sigma");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"theta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"gammaCC\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"beta\",\"type\":{\"name\":\"vector\",\"length\":" << p << "},\"block\":\"parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"theta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"gammaCC\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"beta\",\"type\":{\"name\":\"vector\",\"length\":" << p << "},\"block\":\"parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      std::vector<double> vars_vec(vars.size());
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i) {
        vars.coeffRef(i) = vars_vec[i];
      }
    }
    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }
    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
  
    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits_impl(context, params_i, params_r_vec, pstream);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i) {
        params_r.coeffRef(i) = params_r_vec[i];
      }
    }
    inline void transform_inits(const stan::io::var_context& context,
                                std::vector<int>& params_i,
                                std::vector<double>& vars,
                                std::ostream* pstream = nullptr) const final {
      transform_inits_impl(context, params_i, vars, pstream);
    }        
};
}
using stan_model = model_ContNoborrow_namespace::model_ContNoborrow;
#ifndef USING_R
// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_ContNoborrow_namespace::profiles__;
}
#endif
#endif
