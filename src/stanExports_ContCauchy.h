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
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_ContCauchy_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_ContCauchy");
    reader.add_event(34, 32, "end", "model_ContCauchy");
    return reader;
}
#include <stan_meta_header.hpp>
class model_ContCauchy
  : public stan::model::model_base_crtp<model_ContCauchy> {
private:
        int nCT;
        int nCC;
        int nEC;
        int p;
        vector_d yCT;
        vector_d yCC;
        vector_d yEC;
        std::vector<row_vector_d> xCT;
        std::vector<row_vector_d> xCC;
        std::vector<row_vector_d> xEC;
        double scale;
public:
    model_ContCauchy(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_ContCauchy(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_ContCauchy_namespace::model_ContCauchy";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "nCT", "int", context__.to_vec());
            nCT = int(0);
            vals_i__ = context__.vals_i("nCT");
            pos__ = 0;
            nCT = vals_i__[pos__++];
            check_greater_or_equal(function__, "nCT", nCT, 1);
            current_statement_begin__ = 3;
            context__.validate_dims("data initialization", "nCC", "int", context__.to_vec());
            nCC = int(0);
            vals_i__ = context__.vals_i("nCC");
            pos__ = 0;
            nCC = vals_i__[pos__++];
            check_greater_or_equal(function__, "nCC", nCC, 1);
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "nEC", "int", context__.to_vec());
            nEC = int(0);
            vals_i__ = context__.vals_i("nEC");
            pos__ = 0;
            nEC = vals_i__[pos__++];
            check_greater_or_equal(function__, "nEC", nEC, 1);
            current_statement_begin__ = 5;
            context__.validate_dims("data initialization", "p", "int", context__.to_vec());
            p = int(0);
            vals_i__ = context__.vals_i("p");
            pos__ = 0;
            p = vals_i__[pos__++];
            check_greater_or_equal(function__, "p", p, 1);
            current_statement_begin__ = 6;
            validate_non_negative_index("yCT", "nCT", nCT);
            context__.validate_dims("data initialization", "yCT", "vector_d", context__.to_vec(nCT));
            yCT = Eigen::Matrix<double, Eigen::Dynamic, 1>(nCT);
            vals_r__ = context__.vals_r("yCT");
            pos__ = 0;
            size_t yCT_j_1_max__ = nCT;
            for (size_t j_1__ = 0; j_1__ < yCT_j_1_max__; ++j_1__) {
                yCT(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 7;
            validate_non_negative_index("yCC", "nCC", nCC);
            context__.validate_dims("data initialization", "yCC", "vector_d", context__.to_vec(nCC));
            yCC = Eigen::Matrix<double, Eigen::Dynamic, 1>(nCC);
            vals_r__ = context__.vals_r("yCC");
            pos__ = 0;
            size_t yCC_j_1_max__ = nCC;
            for (size_t j_1__ = 0; j_1__ < yCC_j_1_max__; ++j_1__) {
                yCC(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 8;
            validate_non_negative_index("yEC", "nEC", nEC);
            context__.validate_dims("data initialization", "yEC", "vector_d", context__.to_vec(nEC));
            yEC = Eigen::Matrix<double, Eigen::Dynamic, 1>(nEC);
            vals_r__ = context__.vals_r("yEC");
            pos__ = 0;
            size_t yEC_j_1_max__ = nEC;
            for (size_t j_1__ = 0; j_1__ < yEC_j_1_max__; ++j_1__) {
                yEC(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 9;
            validate_non_negative_index("xCT", "p", p);
            validate_non_negative_index("xCT", "nCT", nCT);
            context__.validate_dims("data initialization", "xCT", "row_vector_d", context__.to_vec(nCT,p));
            xCT = std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic> >(nCT, Eigen::Matrix<double, 1, Eigen::Dynamic>(p));
            vals_r__ = context__.vals_r("xCT");
            pos__ = 0;
            size_t xCT_j_1_max__ = p;
            size_t xCT_k_0_max__ = nCT;
            for (size_t j_1__ = 0; j_1__ < xCT_j_1_max__; ++j_1__) {
                for (size_t k_0__ = 0; k_0__ < xCT_k_0_max__; ++k_0__) {
                    xCT[k_0__](j_1__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 10;
            validate_non_negative_index("xCC", "p", p);
            validate_non_negative_index("xCC", "nCC", nCC);
            context__.validate_dims("data initialization", "xCC", "row_vector_d", context__.to_vec(nCC,p));
            xCC = std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic> >(nCC, Eigen::Matrix<double, 1, Eigen::Dynamic>(p));
            vals_r__ = context__.vals_r("xCC");
            pos__ = 0;
            size_t xCC_j_1_max__ = p;
            size_t xCC_k_0_max__ = nCC;
            for (size_t j_1__ = 0; j_1__ < xCC_j_1_max__; ++j_1__) {
                for (size_t k_0__ = 0; k_0__ < xCC_k_0_max__; ++k_0__) {
                    xCC[k_0__](j_1__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 11;
            validate_non_negative_index("xEC", "p", p);
            validate_non_negative_index("xEC", "nEC", nEC);
            context__.validate_dims("data initialization", "xEC", "row_vector_d", context__.to_vec(nEC,p));
            xEC = std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic> >(nEC, Eigen::Matrix<double, 1, Eigen::Dynamic>(p));
            vals_r__ = context__.vals_r("xEC");
            pos__ = 0;
            size_t xEC_j_1_max__ = p;
            size_t xEC_k_0_max__ = nEC;
            for (size_t j_1__ = 0; j_1__ < xEC_j_1_max__; ++j_1__) {
                for (size_t k_0__ = 0; k_0__ < xEC_k_0_max__; ++k_0__) {
                    xEC[k_0__](j_1__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 12;
            context__.validate_dims("data initialization", "scale", "double", context__.to_vec());
            scale = double(0);
            vals_r__ = context__.vals_r("scale");
            pos__ = 0;
            scale = vals_r__[pos__++];
            check_greater_or_equal(function__, "scale", scale, 0);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 15;
            num_params_r__ += 1;
            current_statement_begin__ = 16;
            num_params_r__ += 1;
            current_statement_begin__ = 17;
            num_params_r__ += 1;
            current_statement_begin__ = 18;
            num_params_r__ += 1;
            current_statement_begin__ = 19;
            validate_non_negative_index("beta", "p", p);
            num_params_r__ += p;
            current_statement_begin__ = 20;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_ContCauchy() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 15;
        if (!(context__.contains_r("theta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable theta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("theta");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "theta", "double", context__.to_vec());
        double theta(0);
        theta = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(theta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable theta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 16;
        if (!(context__.contains_r("gammaCC")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable gammaCC missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("gammaCC");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "gammaCC", "double", context__.to_vec());
        double gammaCC(0);
        gammaCC = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(gammaCC);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable gammaCC: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 17;
        if (!(context__.contains_r("gammaEC")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable gammaEC missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("gammaEC");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "gammaEC", "double", context__.to_vec());
        double gammaEC(0);
        gammaEC = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(gammaEC);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable gammaEC: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 18;
        if (!(context__.contains_r("tau")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable tau missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("tau");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "tau", "double", context__.to_vec());
        double tau(0);
        tau = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, tau);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable tau: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 19;
        if (!(context__.contains_r("beta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "p", p);
        context__.validate_dims("parameter initialization", "beta", "vector_d", context__.to_vec(p));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta(p);
        size_t beta_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            beta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 20;
        if (!(context__.contains_r("sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "sigma", "double", context__.to_vec());
        double sigma(0);
        sigma = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 15;
            local_scalar_t__ theta;
            (void) theta;  // dummy to suppress unused var warning
            if (jacobian__)
                theta = in__.scalar_constrain(lp__);
            else
                theta = in__.scalar_constrain();
            current_statement_begin__ = 16;
            local_scalar_t__ gammaCC;
            (void) gammaCC;  // dummy to suppress unused var warning
            if (jacobian__)
                gammaCC = in__.scalar_constrain(lp__);
            else
                gammaCC = in__.scalar_constrain();
            current_statement_begin__ = 17;
            local_scalar_t__ gammaEC;
            (void) gammaEC;  // dummy to suppress unused var warning
            if (jacobian__)
                gammaEC = in__.scalar_constrain(lp__);
            else
                gammaEC = in__.scalar_constrain();
            current_statement_begin__ = 18;
            local_scalar_t__ tau;
            (void) tau;  // dummy to suppress unused var warning
            if (jacobian__)
                tau = in__.scalar_lb_constrain(0, lp__);
            else
                tau = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 19;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.vector_constrain(p, lp__);
            else
                beta = in__.vector_constrain(p);
            current_statement_begin__ = 20;
            local_scalar_t__ sigma;
            (void) sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma = in__.scalar_lb_constrain(0, lp__);
            else
                sigma = in__.scalar_lb_constrain(0);
            // model body
            current_statement_begin__ = 23;
            lp_accum__.add(cauchy_log<propto__>(tau, 0, scale));
            current_statement_begin__ = 24;
            lp_accum__.add(normal_log<propto__>(gammaCC, gammaEC, tau));
            current_statement_begin__ = 26;
            for (int i = 1; i <= nCT; ++i) {
                current_statement_begin__ = 27;
                lp_accum__.add(normal_log<propto__>(get_base1(yCT, i, "yCT", 1), ((gammaCC + theta) + multiply(get_base1(xCT, i, "xCT", 1), beta)), sigma));
            }
            current_statement_begin__ = 28;
            for (int i = 1; i <= nCC; ++i) {
                current_statement_begin__ = 29;
                lp_accum__.add(normal_log<propto__>(get_base1(yCC, i, "yCC", 1), (gammaCC + multiply(get_base1(xCC, i, "xCC", 1), beta)), sigma));
            }
            current_statement_begin__ = 30;
            for (int i = 1; i <= nEC; ++i) {
                current_statement_begin__ = 31;
                lp_accum__.add(normal_log<propto__>(get_base1(yEC, i, "yEC", 1), (gammaEC + multiply(get_base1(xEC, i, "xEC", 1), beta)), sigma));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("theta");
        names__.push_back("gammaCC");
        names__.push_back("gammaEC");
        names__.push_back("tau");
        names__.push_back("beta");
        names__.push_back("sigma");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(p);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_ContCauchy_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double theta = in__.scalar_constrain();
        vars__.push_back(theta);
        double gammaCC = in__.scalar_constrain();
        vars__.push_back(gammaCC);
        double gammaEC = in__.scalar_constrain();
        vars__.push_back(gammaEC);
        double tau = in__.scalar_lb_constrain(0);
        vars__.push_back(tau);
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta = in__.vector_constrain(p);
        size_t beta_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            vars__.push_back(beta(j_1__));
        }
        double sigma = in__.scalar_lb_constrain(0);
        vars__.push_back(sigma);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_ContCauchy";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "theta";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "gammaCC";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "gammaEC";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau";
        param_names__.push_back(param_name_stream__.str());
        size_t beta_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "theta";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "gammaCC";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "gammaEC";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "tau";
        param_names__.push_back(param_name_stream__.str());
        size_t beta_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_ContCauchy_namespace::model_ContCauchy stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
