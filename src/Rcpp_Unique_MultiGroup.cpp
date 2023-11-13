#include <cassert>
#include <ctime>
#include <iostream>
#include <random>
#include <vector>

// [[Rcpp::plugins(cpp17)]]

// Armadillo
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
//using namespace arma;

// UMG represents Unique Multi-group
void centerNumericMatrix_bis_UMG(Rcpp::NumericMatrix& X) {
  const unsigned int m = X.ncol();
  for (unsigned int j = 0; j < m; ++j) {
    X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
  }
}

void covRcpp_bis_UMG(Rcpp::NumericMatrix& Y,
                 Rcpp::NumericMatrix& cov,
                 double const& c_diag,
                 double const& c_prop, 
                 unsigned int const& r,
                 unsigned int const& K) {
  unsigned const int df = r - 101;
  
  // Centering the matrix!
  centerNumericMatrix_bis_UMG(Y);  // Defined in aux_functions
  
  // COV only updates the bottom right corner element!
  
  // Computing the covariance matrix
  for (unsigned int i = 0; i < K; ++i) {
    for (unsigned int j = 0; j <= i; ++j) { // I REMOVE THE FIRST 100 ROWS OF THE MCMC AND START FROM THE 101-st row.
      cov(i,j) = c_prop * Rcpp::sum(Y(Rcpp::_, i)*Y(Rcpp::_, j))/df;
      cov(j,i) = cov(i,j);
    }
    cov(i,i) += c_diag;
  }
}

void
  my_rmvnorm_final_UMG(Rcpp::NumericVector& alpha_prop,
                   arma::mat const& chol_R, 
                   Rcpp::NumericVector const& mean, 
                   unsigned int const& K){
    arma::vec sample(K);
    
    sample = Rcpp::as<arma::vec>(Rcpp::rnorm(K));
    sample = trans(sample.t() * chol_R);
    
    alpha_prop = Rcpp::wrap(sample);
    alpha_prop += mean;
  }

void 
  prior_informative_UMG(double& prior,
                    Rcpp::NumericVector const& x,
                    double const& mean, 
                    double const& sd,
                    unsigned int const& K){
    double sum_exp_x = 0.0;
    for (unsigned int k = 0; k < K; ++k) {
      sum_exp_x += exp(x[k]);
    }
    double log_disp = log(sum_exp_x);
    
    arma::mat mat(K, K);
    mat.zeros(); 
    // can I initialize a matrix like this with all 0's ???
    for (unsigned int k = 0; k < K; ++k) {
      mat(k,k) = 1;
      mat(K-1,k) = exp(x[k])/sum_exp_x;
    }
    
    // real part of the log_determinant
    //  std::complex<double> tmp = arma::log_det(mat);
    //  prior = tmp.real();
    // TO DO: USE R determinant
    
    double sign;
    arma::log_det(prior, sign, mat);
    
    prior += R::dnorm(log_disp, mean, sd, true);
    for (unsigned int k = 0; k < K-1; ++k) {
      prior += R::dnorm(x[k], mean - log(K), 10, true);
    }
  }


void 
  prior_non_informative_UMG(double& prior, 
                        Rcpp::NumericVector const& x,
                        unsigned int const& K){
    prior = 0.0;
    
    for (unsigned int k = 0; k < K; ++k) {
      prior += R::dnorm(x[k], 0, 10, true);
    }
  }

double
  ll_alpha_new_UMG(Rcpp::NumericMatrix const& pi, 
               Rcpp::NumericVector const& alpha,
               unsigned int const& N,
               unsigned int const& K)
  {
    double ll = 0.0;
    double sum_alpha = 0.0;
    
    for (unsigned int k = 0; k < K; ++k) {
      sum_alpha += alpha[k];
      
      ll -= N * R::lgammafn(alpha[k]);
      
      for (unsigned int i = 0; i < N; ++i) {
        ll += log(pi(i,k)) * ( alpha[k] - 1);
      }
    }
    ll += N * R::lgammafn(sum_alpha);
    
    if( ISNAN(ll) ){
      Rcout << "na ll is" << std::endl << ll << std::endl;
      ll = -std::numeric_limits<double>::infinity();
      Rcout << "ll becomes" << std::endl << ll << std::endl;
    }
    
    return ll;
  }
// I update the log-likelihood value directly (ll)
// check lgammafn function !!!

// [[Rcpp::export]]
Rcpp::List
  Rcpp_FULL_Unique_Multigroup( unsigned int const& K,
                    unsigned int& R,
                    unsigned int const& burn_in,
                    Rcpp::IntegerVector const& N,
                    unsigned int const& N_groups,
                    double const& mean_log_delta,
                    double const& sd_log_delta,
                    Rcpp::ListOf<Rcpp::NumericMatrix>& pi_new,
                    Rcpp::ListOf<Rcpp::NumericMatrix>& mcmc_alpha,
                    Rcpp::ListOf<Rcpp::NumericVector>& alpha_new,
                    Rcpp::ListOf<Rcpp::NumericMatrix>& chol,
                    Rcpp::NumericVector const& l, // double vector of len K
                    Rcpp::ListOf<Rcpp::IntegerMatrix> const& f,
                    Rcpp::IntegerMatrix const& exon_id, // integer matrix (K * J)
                    Rcpp::LogicalVector const& one_transcript) // 
  {
    // N_groups defines the number of groups to compare.
    
    // the proportionality constant for the ARW matrix.
    double c_prop = 0.6/K;

    // define the number of classes
    unsigned int J = exon_id.ncol();
    
    // initialize a few objects once only:
    bool cond;
    unsigned int cond_01;
    double pi_new_i_tot;
    double prob_tot;
    Rcpp::NumericMatrix Y_new(N_groups, K);
    Rcpp::NumericVector pi_new_i(K);

    // Initialize on top:
    Rcpp::NumericVector alpha_prop(K);
    double dir_sample_sum;
    double prior_alpha_prop = 0.0;
    Rcpp::NumericVector prior_alpha_new(N_groups);
    double ll_alpha_prop = 0.0;
    double ll_alpha_new_ = 0.0;
    double alpha;
    
    // store the log-posterior of the hyper-parameters | rest
    Rcpp::NumericVector ll(R);
    
    Rcpp::NumericVector dir_sample(K);
    Rcpp::NumericMatrix cv_alpha(K, K);
    Rcpp::NumericVector prob(K);
    Rcpp::IntegerVector n_multinom(K);
    
    Rcpp::NumericMatrix Y;
    
    // Initialize the prior:
    for (unsigned int n = 0; n < N_groups; ++n) { // n representing the group id
      if( mean_log_delta != 0){ // if prior has been specified
        prior_informative_UMG( prior_alpha_prop, alpha_new[n],  mean_log_delta, sd_log_delta, K  );
      }
      else{ // if prior has NOT been specified: TO DO: compute prior once only at the beginning!
        prior_non_informative_UMG( prior_alpha_prop, alpha_new[n], K  );
      }
      prior_alpha_new(n) = prior_alpha_prop;
    }
    
    // constant to add to the diagonal of the ARW matrix.
    double c_diag = 0.001;
    
    // MCMC iterator:
    for (unsigned int r = 0; r < R; ++r) {
      
      for (unsigned int n = 0; n < N_groups; ++n) { // n representing the group id
        for (int i = 0; i < N[n]; ++i) { // loop on all N[n] samples of group n
          // compute `pi_new_i` element-wise
          for (unsigned int k = 0; k < K; ++k) {
            pi_new_i[k] = pi_new[n](i,k) / l[k];
          }
          // compute sum of vector
          pi_new_i_tot = std::accumulate(pi_new_i.begin(), pi_new_i.end(), 0.0);

          // normalize elements of `pi_new_i`
          for (unsigned int k = 0; k < K; ++k) {
            pi_new_i[k] /= pi_new_i_tot;
          }
          
          for (unsigned int k = 0; k < K; ++k) {
            Y_new(n,k) = exp(alpha_new[n][k]);      
          }
          for (unsigned int j = 0; j < J; ++j) {
            // Group A:
            if (one_transcript[j] or f[n](j, i) == 0) {
              for (unsigned int k = 0; k < K; ++k) {
                Y_new(n,k) += f[n](j, i) * exon_id(k,j);
              }
            }
            else {
              // compute (elementwise): prob = pi_new_i * exon_id_col_j
              for (unsigned int k = 0; k < K; ++k) {
                prob[k] = pi_new_i[k] * exon_id(k,j);
              }
              // normalize elements of prob to sum to 1
              prob_tot = std::accumulate(prob.begin(), prob.end(), 0.0);
              for (unsigned int k = 0; k < K; ++k) {
                prob[k] /= prob_tot;
              }
              
              rmultinom(f[n](j, i), prob.begin(), K, n_multinom.begin());
              //          gsl_ran_multinomial(r, K, f(j, i), prob.begin(), (unsigned int *) n.begin());
              
              for (unsigned int k = 0; k < K; ++k) {
                //  prob.push_back(pi_new_i[k] * exon_id_col_j[k]);
                Y_new(n,k) += n_multinom[k];
              }
            }
          }
        
          // Dirichltet
          dir_sample_sum = 0;
  
          for (unsigned int k = 0; k < K; ++k) {
            dir_sample[k] = as<double>(Rcpp::rgamma(1, Y_new(n,k), 1));
            dir_sample_sum += dir_sample[k];
          }
          // normalize the gamma samples directly when assigning them to the matrix!
          // check that:  sum(is.na(pi)) == 0 AND all( pi > 10^{-100} )
          cond = true;
          for (unsigned int k = 0; k < K; ++k) {
            dir_sample[k] = dir_sample[k]/dir_sample_sum;
            
            if( ( Rcpp::NumericVector::is_na(dir_sample[k]) ) || ( dir_sample[k] < pow(10, -100) )  ) {
              cond = false;
            }
          }
          
          // Update pi_new:
          if(cond){
            for (unsigned int k = 0; k < K; ++k) {
              pi_new[n](i, k) = dir_sample[k];
            }
          }
        }
        
        if( r < 200 ){ // simple R.W. for the initial 200 iter.
          alpha_prop = Rcpp::rnorm(K, 0, 0.1);
          alpha_prop += alpha_new[n];
        }
        else{ // ARW for the following iterations:
          if( (r == 200) || (r == burn_in) ){ //round( r/100) == r/100 ){ // update the covariance matrix every 100 iterations
            // maybe I need to use arma:: mat to do element-wise summation
            Y = mcmc_alpha[n]( Range(100, r-1), Range(0, K-1) ); // SUB-SET of mcmc, mcmc[101:{r-2},1:K]
            covRcpp_bis_UMG(Y, cv_alpha, c_diag, c_prop, r, K);
            chol[n] = Rcpp::wrap( arma::chol(Rcpp::as<arma::mat>(cv_alpha)) );
  
            //Rcout << "chol_A is" << std::endl << chol_A << std::endl;
            //Rcout << "chol_B is" << std::endl << chol_B << std::endl;
          }
          // modify here to compute the cholesky decomposition above.
          my_rmvnorm_final_UMG(alpha_prop, Rcpp::as<arma::mat>(chol[n]), alpha_new[n], K);
        }
        
        if( mean_log_delta != 0){ // if prior has been specified
          prior_informative_UMG( prior_alpha_prop, alpha_prop,  mean_log_delta, sd_log_delta, K  );
        }
        else{ // if prior has NOT been specified: TO DO: compute prior once only at the beginning!
          prior_non_informative_UMG( prior_alpha_prop, alpha_prop, K  );
        }
        
        // WHEN ALL VALUES ARE REJECTED: error: chol(): decomposition failed
        
        // ADD conditions on A and B:
        cond = true;

        // check if there are NA's in the proposed values:
        for (unsigned int k = 0; k < K; ++k) {
          if( NumericVector::is_na(alpha_prop[k]) ){
            cond = false;
            Rcout << "1st control failed: " << std::endl << alpha_prop << std::endl;
          }
        }
        
        // only if there are no NA's:
        // check that:
        // 1) All exp(alpha_prop) are strictly > 0;
        // 2) the sum of exp(alpha_prop) is < 10^10
        if(cond){
          if( !is_true( all( exp(alpha_prop) > 0 ) ) || (sum(exp(alpha_prop)) >  pow(10, 10) ) ){
            cond = false;
            Rcout << "2nd control failed: " << std::endl << alpha_prop << std::endl;
          }
        }
        
        if(cond){
          ll_alpha_prop = ll_alpha_new_UMG( pi_new[n], exp(alpha_prop),    N[n], K);
          ll_alpha_new_ = ll_alpha_new_UMG( pi_new[n], exp(alpha_new[n]),  N[n], K);
          
          // ACCEPT/REJECT proposal for group A:
          alpha = std::min( ll_alpha_prop - ll_alpha_new_ + prior_alpha_prop - prior_alpha_new[n], 0.0);
          if( !NumericVector::is_na(alpha) ){ // check no NA's in alpha.
            cond_01 = Rcpp::rbinom(1, 1, exp(alpha))(0);
            
            if(cond_01 == 1){ // # Update parameter and prior
              alpha_new[n] = alpha_prop;
              prior_alpha_new[n] = prior_alpha_prop;
              ll_alpha_new_ = ll_alpha_prop;
            }
          }
        }
        
        // IMPORTANT:
        // LIKELIHOOD: put a control when res is nan and -inf, otherwise all values accepted!
        // DONE IT: it's in alpha (prob of accepting the new value).
        
        // update mcmc matrixes
        mcmc_alpha[n](r,_) = alpha_new[n];
        //      for (unsigned int k = 0; k < K; ++k) {
        //        mcmc_alpha_A(r,k) = alpha_new_A[k];
        //        mcmc_alpha_B(r,k) = alpha_new_B[k];
        //      }
        ll[r] += prior_alpha_new[n] + ll_alpha_new_;
      }
    }
    
    Rcpp::NumericVector row(K);
    Rcpp::NumericMatrix precision(R-burn_in, N_groups); //precision for group A
    
    for (unsigned int n = 0; n < N_groups; ++n) { // n representing the group id
      for (unsigned int r = burn_in; r < R; ++r) {
        row = Rcpp::exp(mcmc_alpha[n](r,_)); // exp to gain the alpha's
        precision(r-burn_in, n) = log(sum(row)); //precision for group A (before deviding row by l)
        
        // row = row/sum(row); // divide by their sum to get pi_bar
        row = row/l; // divide by the transcript effective length
        mcmc_alpha[n](r,_) = row/sum(row); // standardize to obtain proprotions again.
      }
    
      mcmc_alpha[n] = mcmc_alpha[n]( Range(burn_in, R-1), Range(0, K-1) );
    }
    
    // THIN HERE: return 10^4 values (1.2 * 10^4 for ll).
    return Rcpp::List::create(Rcpp::Named("mcmc") = mcmc_alpha,
                              Rcpp::Named("log-posterior") = ll[ Range(burn_in, R-1) ],
                              Rcpp::Named("log-precision") = precision);
  } // return the ll WITH the burn-in
// the original element pi_new is modified
