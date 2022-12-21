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
// Code generated by stanc v2.26.1-1-g67504470
#include <stan/model/model_header.hpp>
namespace model_T2EFullborrowC0_namespace {
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
                                                      " (in 'string', line 20, column 2 to column 13)",
                                                      " (in 'string', line 21, column 2 to column 15)",
                                                      " (in 'string', line 22, column 2 to column 17)",
                                                      " (in 'string', line 23, column 2 to column 22)",
                                                      " (in 'string', line 26, column 2 to column 25)",
                                                      " (in 'string', line 28, column 4 to column 63)",
                                                      " (in 'string', line 27, column 2 to line 28, column 63)",
                                                      " (in 'string', line 30, column 4 to column 63)",
                                                      " (in 'string', line 29, column 2 to line 30, column 63)",
                                                      " (in 'string', line 32, column 4 to column 63)",
                                                      " (in 'string', line 31, column 2 to line 32, column 63)",
                                                      " (in 'string', line 34, column 4 to column 77)",
                                                      " (in 'string', line 33, column 2 to line 34, column 77)",
                                                      " (in 'string', line 36, column 4 to column 77)",
                                                      " (in 'string', line 35, column 2 to line 36, column 77)",
                                                      " (in 'string', line 2, column 2 to column 21)",
                                                      " (in 'string', line 3, column 2 to column 21)",
                                                      " (in 'string', line 4, column 2 to column 21)",
                                                      " (in 'string', line 5, column 2 to column 21)",
                                                      " (in 'string', line 6, column 2 to column 21)",
                                                      " (in 'string', line 7, column 2 to column 17)",
                                                      " (in 'string', line 8, column 9 to column 14)",
                                                      " (in 'string', line 8, column 2 to column 22)",
                                                      " (in 'string', line 9, column 9 to column 14)",
                                                      " (in 'string', line 9, column 2 to column 22)",
                                                      " (in 'string', line 10, column 9 to column 14)",
                                                      " (in 'string', line 10, column 2 to column 22)",
                                                      " (in 'string', line 11, column 9 to column 14)",
                                                      " (in 'string', line 11, column 2 to column 22)",
                                                      " (in 'string', line 12, column 9 to column 14)",
                                                      " (in 'string', line 12, column 2 to column 22)",
                                                      " (in 'string', line 13, column 22 to column 27)",
                                                      " (in 'string', line 13, column 13 to column 14)",
                                                      " (in 'string', line 13, column 2 to column 29)",
                                                      " (in 'string', line 14, column 22 to column 27)",
                                                      " (in 'string', line 14, column 13 to column 14)",
                                                      " (in 'string', line 14, column 2 to column 29)",
                                                      " (in 'string', line 15, column 22 to column 27)",
                                                      " (in 'string', line 15, column 13 to column 14)",
                                                      " (in 'string', line 15, column 2 to column 29)",
                                                      " (in 'string', line 16, column 22 to column 27)",
                                                      " (in 'string', line 16, column 13 to column 14)",
                                                      " (in 'string', line 16, column 2 to column 29)",
                                                      " (in 'string', line 17, column 22 to column 27)",
                                                      " (in 'string', line 17, column 13 to column 14)",
                                                      " (in 'string', line 17, column 2 to column 29)",
                                                      " (in 'string', line 22, column 9 to column 10)"};
#include <stan_meta_header.hpp>
class model_T2EFullborrowC0 final : public model_base_crtp<model_T2EFullborrowC0> {
private:
  int nCT_o;
  int nCT_c;
  int nCC_o;
  int nEC_o;
  int nEC_c;
  int p;
  Eigen::Matrix<double, -1, 1> yCT_o;
  Eigen::Matrix<double, -1, 1> yCT_c;
  Eigen::Matrix<double, -1, 1> yCC_o;
  Eigen::Matrix<double, -1, 1> yEC_o;
  Eigen::Matrix<double, -1, 1> yEC_c;
  std::vector<Eigen::Matrix<double, 1, -1>> xCT_o;
  std::vector<Eigen::Matrix<double, 1, -1>> xCT_c;
  std::vector<Eigen::Matrix<double, 1, -1>> xCC_o;
  std::vector<Eigen::Matrix<double, 1, -1>> xEC_o;
  std::vector<Eigen::Matrix<double, 1, -1>> xEC_c;
 
public:
  ~model_T2EFullborrowC0() { }
  
  inline std::string model_name() const final { return "model_T2EFullborrowC0"; }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.26.1-1-g67504470", "stancflags = "};
  }
  
  
  model_T2EFullborrowC0(stan::io::var_context& context__,
                        unsigned int random_seed__ = 0,
                        std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "model_T2EFullborrowC0_namespace::model_T2EFullborrowC0";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 16;
      context__.validate_dims("data initialization","nCT_o","int",
          context__.to_vec());
      nCT_o = std::numeric_limits<int>::min();
      
      current_statement__ = 16;
      nCT_o = context__.vals_i("nCT_o")[(1 - 1)];
      current_statement__ = 16;
      current_statement__ = 16;
      check_greater_or_equal(function__, "nCT_o", nCT_o, 1);
      current_statement__ = 17;
      context__.validate_dims("data initialization","nCT_c","int",
          context__.to_vec());
      nCT_c = std::numeric_limits<int>::min();
      
      current_statement__ = 17;
      nCT_c = context__.vals_i("nCT_c")[(1 - 1)];
      current_statement__ = 17;
      current_statement__ = 17;
      check_greater_or_equal(function__, "nCT_c", nCT_c, 1);
      current_statement__ = 18;
      context__.validate_dims("data initialization","nCC_o","int",
          context__.to_vec());
      nCC_o = std::numeric_limits<int>::min();
      
      current_statement__ = 18;
      nCC_o = context__.vals_i("nCC_o")[(1 - 1)];
      current_statement__ = 18;
      current_statement__ = 18;
      check_greater_or_equal(function__, "nCC_o", nCC_o, 1);
      current_statement__ = 19;
      context__.validate_dims("data initialization","nEC_o","int",
          context__.to_vec());
      nEC_o = std::numeric_limits<int>::min();
      
      current_statement__ = 19;
      nEC_o = context__.vals_i("nEC_o")[(1 - 1)];
      current_statement__ = 19;
      current_statement__ = 19;
      check_greater_or_equal(function__, "nEC_o", nEC_o, 1);
      current_statement__ = 20;
      context__.validate_dims("data initialization","nEC_c","int",
          context__.to_vec());
      nEC_c = std::numeric_limits<int>::min();
      
      current_statement__ = 20;
      nEC_c = context__.vals_i("nEC_c")[(1 - 1)];
      current_statement__ = 20;
      current_statement__ = 20;
      check_greater_or_equal(function__, "nEC_c", nEC_c, 1);
      current_statement__ = 21;
      context__.validate_dims("data initialization","p","int",
          context__.to_vec());
      p = std::numeric_limits<int>::min();
      
      current_statement__ = 21;
      p = context__.vals_i("p")[(1 - 1)];
      current_statement__ = 21;
      current_statement__ = 21;
      check_greater_or_equal(function__, "p", p, 1);
      current_statement__ = 22;
      validate_non_negative_index("yCT_o", "nCT_o", nCT_o);
      current_statement__ = 23;
      context__.validate_dims("data initialization","yCT_o","double",
          context__.to_vec(nCT_o));
      yCT_o = Eigen::Matrix<double, -1, 1>(nCT_o);
      stan::math::fill(yCT_o, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> yCT_o_flat__;
        current_statement__ = 23;
        assign(yCT_o_flat__, nil_index_list(), context__.vals_r("yCT_o"),
          "assigning variable yCT_o_flat__");
        current_statement__ = 23;
        pos__ = 1;
        current_statement__ = 23;
        for (int sym1__ = 1; sym1__ <= nCT_o; ++sym1__) {
          current_statement__ = 23;
          assign(yCT_o, cons_list(index_uni(sym1__), nil_index_list()),
            yCT_o_flat__[(pos__ - 1)], "assigning variable yCT_o");
          current_statement__ = 23;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 24;
      validate_non_negative_index("yCT_c", "nCT_c", nCT_c);
      current_statement__ = 25;
      context__.validate_dims("data initialization","yCT_c","double",
          context__.to_vec(nCT_c));
      yCT_c = Eigen::Matrix<double, -1, 1>(nCT_c);
      stan::math::fill(yCT_c, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> yCT_c_flat__;
        current_statement__ = 25;
        assign(yCT_c_flat__, nil_index_list(), context__.vals_r("yCT_c"),
          "assigning variable yCT_c_flat__");
        current_statement__ = 25;
        pos__ = 1;
        current_statement__ = 25;
        for (int sym1__ = 1; sym1__ <= nCT_c; ++sym1__) {
          current_statement__ = 25;
          assign(yCT_c, cons_list(index_uni(sym1__), nil_index_list()),
            yCT_c_flat__[(pos__ - 1)], "assigning variable yCT_c");
          current_statement__ = 25;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 26;
      validate_non_negative_index("yCC_o", "nCC_o", nCC_o);
      current_statement__ = 27;
      context__.validate_dims("data initialization","yCC_o","double",
          context__.to_vec(nCC_o));
      yCC_o = Eigen::Matrix<double, -1, 1>(nCC_o);
      stan::math::fill(yCC_o, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> yCC_o_flat__;
        current_statement__ = 27;
        assign(yCC_o_flat__, nil_index_list(), context__.vals_r("yCC_o"),
          "assigning variable yCC_o_flat__");
        current_statement__ = 27;
        pos__ = 1;
        current_statement__ = 27;
        for (int sym1__ = 1; sym1__ <= nCC_o; ++sym1__) {
          current_statement__ = 27;
          assign(yCC_o, cons_list(index_uni(sym1__), nil_index_list()),
            yCC_o_flat__[(pos__ - 1)], "assigning variable yCC_o");
          current_statement__ = 27;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 28;
      validate_non_negative_index("yEC_o", "nEC_o", nEC_o);
      current_statement__ = 29;
      context__.validate_dims("data initialization","yEC_o","double",
          context__.to_vec(nEC_o));
      yEC_o = Eigen::Matrix<double, -1, 1>(nEC_o);
      stan::math::fill(yEC_o, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> yEC_o_flat__;
        current_statement__ = 29;
        assign(yEC_o_flat__, nil_index_list(), context__.vals_r("yEC_o"),
          "assigning variable yEC_o_flat__");
        current_statement__ = 29;
        pos__ = 1;
        current_statement__ = 29;
        for (int sym1__ = 1; sym1__ <= nEC_o; ++sym1__) {
          current_statement__ = 29;
          assign(yEC_o, cons_list(index_uni(sym1__), nil_index_list()),
            yEC_o_flat__[(pos__ - 1)], "assigning variable yEC_o");
          current_statement__ = 29;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 30;
      validate_non_negative_index("yEC_c", "nEC_c", nEC_c);
      current_statement__ = 31;
      context__.validate_dims("data initialization","yEC_c","double",
          context__.to_vec(nEC_c));
      yEC_c = Eigen::Matrix<double, -1, 1>(nEC_c);
      stan::math::fill(yEC_c, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> yEC_c_flat__;
        current_statement__ = 31;
        assign(yEC_c_flat__, nil_index_list(), context__.vals_r("yEC_c"),
          "assigning variable yEC_c_flat__");
        current_statement__ = 31;
        pos__ = 1;
        current_statement__ = 31;
        for (int sym1__ = 1; sym1__ <= nEC_c; ++sym1__) {
          current_statement__ = 31;
          assign(yEC_c, cons_list(index_uni(sym1__), nil_index_list()),
            yEC_c_flat__[(pos__ - 1)], "assigning variable yEC_c");
          current_statement__ = 31;
          pos__ = (pos__ + 1);}
      }
      current_statement__ = 32;
      validate_non_negative_index("xCT_o", "nCT_o", nCT_o);
      current_statement__ = 33;
      validate_non_negative_index("xCT_o", "p", p);
      current_statement__ = 34;
      context__.validate_dims("data initialization","xCT_o","double",
          context__.to_vec(nCT_o, p));
      xCT_o = std::vector<Eigen::Matrix<double, 1, -1>>(nCT_o, Eigen::Matrix<double, 1, -1>(p));
      stan::math::fill(xCT_o, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> xCT_o_flat__;
        current_statement__ = 34;
        assign(xCT_o_flat__, nil_index_list(), context__.vals_r("xCT_o"),
          "assigning variable xCT_o_flat__");
        current_statement__ = 34;
        pos__ = 1;
        current_statement__ = 34;
        for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
          current_statement__ = 34;
          for (int sym2__ = 1; sym2__ <= nCT_o; ++sym2__) {
            current_statement__ = 34;
            assign(xCT_o,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              xCT_o_flat__[(pos__ - 1)], "assigning variable xCT_o");
            current_statement__ = 34;
            pos__ = (pos__ + 1);}}
      }
      current_statement__ = 35;
      validate_non_negative_index("xCT_c", "nCT_c", nCT_c);
      current_statement__ = 36;
      validate_non_negative_index("xCT_c", "p", p);
      current_statement__ = 37;
      context__.validate_dims("data initialization","xCT_c","double",
          context__.to_vec(nCT_c, p));
      xCT_c = std::vector<Eigen::Matrix<double, 1, -1>>(nCT_c, Eigen::Matrix<double, 1, -1>(p));
      stan::math::fill(xCT_c, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> xCT_c_flat__;
        current_statement__ = 37;
        assign(xCT_c_flat__, nil_index_list(), context__.vals_r("xCT_c"),
          "assigning variable xCT_c_flat__");
        current_statement__ = 37;
        pos__ = 1;
        current_statement__ = 37;
        for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
          current_statement__ = 37;
          for (int sym2__ = 1; sym2__ <= nCT_c; ++sym2__) {
            current_statement__ = 37;
            assign(xCT_c,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              xCT_c_flat__[(pos__ - 1)], "assigning variable xCT_c");
            current_statement__ = 37;
            pos__ = (pos__ + 1);}}
      }
      current_statement__ = 38;
      validate_non_negative_index("xCC_o", "nCC_o", nCC_o);
      current_statement__ = 39;
      validate_non_negative_index("xCC_o", "p", p);
      current_statement__ = 40;
      context__.validate_dims("data initialization","xCC_o","double",
          context__.to_vec(nCC_o, p));
      xCC_o = std::vector<Eigen::Matrix<double, 1, -1>>(nCC_o, Eigen::Matrix<double, 1, -1>(p));
      stan::math::fill(xCC_o, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> xCC_o_flat__;
        current_statement__ = 40;
        assign(xCC_o_flat__, nil_index_list(), context__.vals_r("xCC_o"),
          "assigning variable xCC_o_flat__");
        current_statement__ = 40;
        pos__ = 1;
        current_statement__ = 40;
        for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
          current_statement__ = 40;
          for (int sym2__ = 1; sym2__ <= nCC_o; ++sym2__) {
            current_statement__ = 40;
            assign(xCC_o,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              xCC_o_flat__[(pos__ - 1)], "assigning variable xCC_o");
            current_statement__ = 40;
            pos__ = (pos__ + 1);}}
      }
      current_statement__ = 41;
      validate_non_negative_index("xEC_o", "nEC_o", nEC_o);
      current_statement__ = 42;
      validate_non_negative_index("xEC_o", "p", p);
      current_statement__ = 43;
      context__.validate_dims("data initialization","xEC_o","double",
          context__.to_vec(nEC_o, p));
      xEC_o = std::vector<Eigen::Matrix<double, 1, -1>>(nEC_o, Eigen::Matrix<double, 1, -1>(p));
      stan::math::fill(xEC_o, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> xEC_o_flat__;
        current_statement__ = 43;
        assign(xEC_o_flat__, nil_index_list(), context__.vals_r("xEC_o"),
          "assigning variable xEC_o_flat__");
        current_statement__ = 43;
        pos__ = 1;
        current_statement__ = 43;
        for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
          current_statement__ = 43;
          for (int sym2__ = 1; sym2__ <= nEC_o; ++sym2__) {
            current_statement__ = 43;
            assign(xEC_o,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              xEC_o_flat__[(pos__ - 1)], "assigning variable xEC_o");
            current_statement__ = 43;
            pos__ = (pos__ + 1);}}
      }
      current_statement__ = 44;
      validate_non_negative_index("xEC_c", "nEC_c", nEC_c);
      current_statement__ = 45;
      validate_non_negative_index("xEC_c", "p", p);
      current_statement__ = 46;
      context__.validate_dims("data initialization","xEC_c","double",
          context__.to_vec(nEC_c, p));
      xEC_c = std::vector<Eigen::Matrix<double, 1, -1>>(nEC_c, Eigen::Matrix<double, 1, -1>(p));
      stan::math::fill(xEC_c, std::numeric_limits<double>::quiet_NaN());
      
      {
        std::vector<local_scalar_t__> xEC_c_flat__;
        current_statement__ = 46;
        assign(xEC_c_flat__, nil_index_list(), context__.vals_r("xEC_c"),
          "assigning variable xEC_c_flat__");
        current_statement__ = 46;
        pos__ = 1;
        current_statement__ = 46;
        for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
          current_statement__ = 46;
          for (int sym2__ = 1; sym2__ <= nEC_c; ++sym2__) {
            current_statement__ = 46;
            assign(xEC_c,
              cons_list(index_uni(sym2__),
                cons_list(index_uni(sym1__), nil_index_list())),
              xEC_c_flat__[(pos__ - 1)], "assigning variable xEC_c");
            current_statement__ = 46;
            pos__ = (pos__ + 1);}}
      }
      current_statement__ = 47;
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
    static const char* function__ = "model_T2EFullborrowC0_namespace::log_prob";
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
      local_scalar_t__ alpha;
      alpha = DUMMY_VAR__;
      
      current_statement__ = 4;
      alpha = in__.scalar();
      current_statement__ = 4;
      if (jacobian__) {
        current_statement__ = 4;
        alpha = stan::math::lb_constrain(alpha, 0, lp__);
      } else {
        current_statement__ = 4;
        alpha = stan::math::lb_constrain(alpha, 0);
      }
      {
        current_statement__ = 5;
        lp_accum__.add(exponential_lpdf<propto__>(alpha, 1));
        current_statement__ = 7;
        for (int i = 1; i <= nCT_o; ++i) {
          current_statement__ = 6;
          lp_accum__.add(
            weibull_lpdf<propto__>(yCT_o[(i - 1)], alpha,
              stan::math::exp(
                ((gammaCC + theta) + multiply(xCT_o[(i - 1)], beta)))));}
        current_statement__ = 9;
        for (int i = 1; i <= nCC_o; ++i) {
          current_statement__ = 8;
          lp_accum__.add(
            weibull_lpdf<propto__>(yCC_o[(i - 1)], alpha,
              stan::math::exp((gammaCC + multiply(xCC_o[(i - 1)], beta)))));}
        current_statement__ = 11;
        for (int i = 1; i <= nEC_o; ++i) {
          current_statement__ = 10;
          lp_accum__.add(
            weibull_lpdf<propto__>(yEC_o[(i - 1)], alpha,
              stan::math::exp((gammaCC + multiply(xEC_o[(i - 1)], beta)))));}
        current_statement__ = 13;
        for (int i = 1; i <= nCT_c; ++i) {
          current_statement__ = 12;
          lp_accum__.add(
            weibull_lccdf(yCT_c[(i - 1)], alpha,
              stan::math::exp(
                ((gammaCC + theta) + multiply(xCT_c[(i - 1)], beta)))));}
        current_statement__ = 15;
        for (int i = 1; i <= nEC_c; ++i) {
          current_statement__ = 14;
          lp_accum__.add(
            weibull_lccdf(yEC_c[(i - 1)], alpha,
              stan::math::exp((gammaCC + multiply(xEC_c[(i - 1)], beta)))));}
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
    static const char* function__ = "model_T2EFullborrowC0_namespace::write_array";
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
      double alpha;
      alpha = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      alpha = in__.scalar();
      current_statement__ = 4;
      alpha = stan::math::lb_constrain(alpha, 0);
      vars__.emplace_back(theta);
      vars__.emplace_back(gammaCC);
      for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
        vars__.emplace_back(beta[(sym1__ - 1)]);}
      vars__.emplace_back(alpha);
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
      double alpha;
      alpha = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      alpha = context__.vals_r("alpha")[(1 - 1)];
      double alpha_free__;
      alpha_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 4;
      alpha_free__ = stan::math::lb_free(alpha, 0);
      vars__.emplace_back(theta);
      vars__.emplace_back(gammaCC);
      for (int sym1__ = 1; sym1__ <= p; ++sym1__) {
        vars__.emplace_back(beta[(sym1__ - 1)]);}
      vars__.emplace_back(alpha_free__);
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
    names__.emplace_back("alpha");
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
    param_names__.emplace_back(std::string() + "alpha");
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
    param_names__.emplace_back(std::string() + "alpha");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"theta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"gammaCC\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"beta\",\"type\":{\"name\":\"vector\",\"length\":" << p << "},\"block\":\"parameters\"},{\"name\":\"alpha\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"theta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"gammaCC\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"beta\",\"type\":{\"name\":\"vector\",\"length\":" << p << "},\"block\":\"parameters\"},{\"name\":\"alpha\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
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
using stan_model = model_T2EFullborrowC0_namespace::model_T2EFullborrowC0;
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
  return model_T2EFullborrowC0_namespace::profiles__;
}
#endif
#endif
