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

void centerNumericMatrix_bis(Rcpp::NumericMatrix& X) {
  const unsigned int m = X.ncol();
  for (unsigned int j = 0; j < m; ++j) {
    X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
  }
}

void covRcpp_bis(Rcpp::NumericMatrix& mcmc,
                 Rcpp::NumericMatrix& cov,
                 double const& c_diag,
                 double const& c_prop, 
                 unsigned int const& r,
                 unsigned int const& K) {
  unsigned const int df = r - 101;
  
  // do this outside, once only!
  Rcpp::NumericMatrix Y = mcmc( Range(100, r-1), Range(0, K-1) ); // SUB-SET of mcmc, mcmc[101:{r-2},1:K]
  
  // Centering the matrix!
  centerNumericMatrix_bis(Y);  // Defined in aux_functions
  
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
  my_rmvnorm_final(Rcpp::NumericVector& alpha_prop,
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
  prior_informative(double& prior,
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
  prior_non_informative(double& prior, 
                        Rcpp::NumericVector const& x,
                        unsigned int const& K){
    prior = 0.0;
    
    for (unsigned int k = 0; k < K; ++k) {
      prior += R::dnorm(x[k], 0, 10, true);
    }
  }

double
  ll_alpha_new(Rcpp::NumericMatrix const& pi, 
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
  Rcpp_Unique( unsigned int const& K,
               unsigned int& R,
               unsigned int const& burn_in,
               unsigned int const& N_1,
               unsigned int const& N_2,
               double const& mean_log_delta,
               double const& sd_log_delta,
               Rcpp::NumericVector const& l, // double vector of len K
               Rcpp::IntegerMatrix const& f, // integer matrix (J * 2*N)
               Rcpp::IntegerMatrix const& exon_id, // integer matrix (K * J)
               Rcpp::LogicalVector const& one_transcript) // 
  {
    // the proportionality constant for the ARW matrix.
    double c_prop = 0.6/K;
    
    // increment R by the burn-in
    R += burn_in;
    
    // define the MCMC matrixes
    Rcpp::NumericMatrix mcmc_alpha_A(R, K);
    Rcpp::NumericMatrix mcmc_alpha_B(R, K);
    
    // define the pi matrixes
    Rcpp::NumericMatrix pi_new_A(N_1, K);
    Rcpp::NumericMatrix pi_new_B(N_2, K);
    std::fill(pi_new_A.begin(), pi_new_A.end(), 1.0/K);
    std::fill(pi_new_B.begin(), pi_new_B.end(), 1.0/K);
    
    // define the number of classes
    unsigned int J = exon_id.ncol();
    
    // initialize a few objects once only:
    bool cond_A;
    bool cond_B;
    unsigned int cond_01;
    double pi_new_i_tot;
    double prob_tot;
    Rcpp::NumericVector Y_new(K);
    Rcpp::NumericVector pi_new_i(K);
    
    // Initialize on top:
    Rcpp::NumericVector alpha_new_A(K);
    Rcpp::NumericVector alpha_new_B(K);
    Rcpp::NumericVector alpha_prop_A(K);
    Rcpp::NumericVector alpha_prop_B(K);
    double dir_sample_sum;
    double prior_alpha_prop_A = 0.0;
    double prior_alpha_prop_B = 0.0;
    double prior_alpha_new_A = 0.0;
    double prior_alpha_new_B = 0.0;
    double ll_alpha_prop_A = 0.0;
    double ll_alpha_new_A = 0.0;
    double ll_alpha_prop_B = 0.0;
    double ll_alpha_new_B = 0.0;
    double alpha;
    
    // store the log-posterior of the hyper-parameters | rest
    Rcpp::NumericVector ll(R);
    
    Rcpp::NumericVector dir_sample(K);
    
    Rcpp::NumericMatrix cv_alpha_A(K, K);
    Rcpp::NumericMatrix cv_alpha_B(K, K);
    
    Rcpp::NumericVector prob(K);
    Rcpp::IntegerVector n(K);
    
    // constant to add to the diagonal of the ARW matrix.
    double c_diag = 0.001;
    
    arma::mat chol_A(K, K);
    arma::mat chol_B(K, K);
    // ARW and N and N:
    
    // Initialize the prior:
    if( mean_log_delta != 0){ // if prior has been specified
      prior_informative( prior_alpha_new_A, alpha_new_A,  mean_log_delta, sd_log_delta, K  );
      prior_informative( prior_alpha_new_B, alpha_new_B,  mean_log_delta, sd_log_delta, K  );
    }
    else{ // if prior has NOT been specified: TO DO: compute prior once only at the beginning!
      prior_non_informative( prior_alpha_new_A, alpha_new_A, K  );
      prior_non_informative( prior_alpha_new_B, alpha_new_B, K  );
    }
    
    // MCMC iterator:
    for (unsigned int r = 0; r < R; ++r) {
      
      // Group A, Multinomial sampling of X | pi
      for (unsigned int i = 0; i < N_1; ++i) {
        // compute `pi_new_i` element-wise
        for (unsigned int k = 0; k < K; ++k) {
          pi_new_i[k] = pi_new_A(i,k) / l[k];
        }
        // compute sum of vector
        pi_new_i_tot = std::accumulate(pi_new_i.begin(), pi_new_i.end(), 0.0);
        
        // normalize elements of `pi_new_i`
        for (unsigned int k = 0; k < K; ++k) {
          pi_new_i[k] /= pi_new_i_tot;
        }
        
        for (unsigned int k = 0; k < K; ++k) {
          Y_new[k] = exp(alpha_new_A[k]);      
        }
        
        for (unsigned int j = 0; j < J; ++j) {
          // Group A:
          if (one_transcript[j] or f(j, i) == 0) {
            for (unsigned int k = 0; k < K; ++k) {
              Y_new[k] += f(j, i) * exon_id(k,j);
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
            
            rmultinom(f(j, i), prob.begin(), K, n.begin());
            //          gsl_ran_multinomial(r, K, f(j, i), prob.begin(), (unsigned int *) n.begin());
            
            for (unsigned int k = 0; k < K; ++k) {
              Y_new[k] += n[k];
            }
          }
        }
        
        // Group A: Dirichlter sampling of pi | X, delta
        dir_sample_sum = 0;
        
        for (unsigned int k = 0; k < K; ++k) {
          dir_sample[k] = as<double>(Rcpp::rgamma(1, Y_new[k], 1));
          dir_sample_sum += dir_sample[k];
        }
        // normalize the gamma samples directly when assigning them to the matrix!
        // check that:  sum(is.na(pi)) == 0 AND all( pi > 10^{-100} )
        cond_A = true;
        for (unsigned int k = 0; k < K; ++k) {
          dir_sample[k] = dir_sample[k]/dir_sample_sum;
          
          if( ( Rcpp::NumericVector::is_na(dir_sample[k]) ) || ( dir_sample[k] < pow(10, -100) )  ) {
            cond_A = false;
          }
        }
        
        // Update pi in group A:
        if(cond_A){
          for (unsigned int k = 0; k < K; ++k) {
            pi_new_A(i, k) = dir_sample[k];
          }
        }
      }
      
      // Group B, Multinomial sampling of X | pi
      for (unsigned int i = 0; i < N_2; ++i) {
        // compute `pi_new_i` element-wise
        for (unsigned int k = 0; k < K; ++k) {
          pi_new_i[k] = pi_new_B(i,k) / l[k];
        }
        // compute sum of vector
        pi_new_i_tot = std::accumulate(pi_new_i.begin(), pi_new_i.end(), 0.0);
        
        // normalize elements of `pi_new_i`
        for (unsigned int k = 0; k < K; ++k) {
          pi_new_i[k] /= pi_new_i_tot;
        }
        
        for (unsigned int k = 0; k < K; ++k) {
          Y_new[k] = exp(alpha_new_B[k]);
        }
        for (unsigned int j = 0; j < J; ++j) {
          // Group B:
          if (one_transcript[j] or f(j, N_1+i) == 0) {
            for (unsigned int k = 0; k < K; ++k) {
              Y_new[k] += f(j, N_1+i) * exon_id(k,j);
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
            
            rmultinom(f(j, N_1+i), prob.begin(), K, n.begin());
            //          gsl_ran_multinomial(r, K, f(j, N+i), prob.begin(), (unsigned int *) n.begin());
            
            for (unsigned int k = 0; k < K; ++k) {
              Y_new[k] += n[k];
            }
          }
        }
        
        // Group B: Dirichlter sampling of pi | X, delta
        dir_sample_sum = 0;
        
        for (unsigned int k = 0; k < K; ++k) {
          dir_sample[k] = as<double>(Rcpp::rgamma(1, Y_new[k], 1));
          dir_sample_sum += dir_sample[k];
        }
        // normalize the gamma samples directly when assigning them to the matrix!
        // check that:  sum(is.na(pi)) == 0 AND all( pi > 10^{-100} )
        cond_B = true;
        for (unsigned int k = 0; k < K; ++k) {
          dir_sample[k] = dir_sample[k]/dir_sample_sum;
          
          if( ( Rcpp::NumericVector::is_na(dir_sample[k]) ) || ( dir_sample[k] < pow(10, -100) )  ) {
            cond_B = false;
          }
        }
        
        // Update pi in group B:
        if(cond_B){
          for (unsigned int k = 0; k < K; ++k) {
            pi_new_B(i, k) = dir_sample[k];
          }
        }
      }
      
      // Sample Dirichlet parameters | pi's
      if( r < 200 ){ // simple R.W. for the initial 200 iter.
        alpha_prop_A = Rcpp::rnorm(K, 0, 0.1);
        alpha_prop_A += alpha_new_A;
        alpha_prop_B = Rcpp::rnorm(K, 0, 0.1);
        alpha_prop_B += alpha_new_B;
      }
      else{ // ARW for the following iterations:
        if( (r == 200) || (r == burn_in) ){ //round( r/100) == r/100 ){ // update the covariance matrix every 100 iterations
          // maybe I need to use arma:: mat to do element-wise summation
          covRcpp_bis(mcmc_alpha_A, cv_alpha_A, c_diag, c_prop, r, K);
          covRcpp_bis(mcmc_alpha_B, cv_alpha_B, c_diag, c_prop, r, K);
          
          chol_A = arma::chol(Rcpp::as<arma::mat>(cv_alpha_A)); // cholesky decomposition of the covariance matrixes for the ARW
          chol_B = arma::chol(Rcpp::as<arma::mat>(cv_alpha_B));
        }
        // modify here to compute the cholesky decomposition above.
        my_rmvnorm_final(alpha_prop_A, chol_A, alpha_new_A, K);
        my_rmvnorm_final(alpha_prop_B, chol_B, alpha_new_B, K);
      }
      
      if( mean_log_delta != 0){ // if prior has been specified
        prior_informative( prior_alpha_prop_A, alpha_prop_A,  mean_log_delta, sd_log_delta, K  );
        prior_informative( prior_alpha_prop_B, alpha_prop_B,  mean_log_delta, sd_log_delta, K  );
      }
      else{ // if prior has NOT been specified: TO DO: compute prior once only at the beginning!
        prior_non_informative( prior_alpha_prop_A, alpha_prop_A, K  );
        prior_non_informative( prior_alpha_prop_B, alpha_prop_B, K  );
      }
      
      // WHEN ALL VALUES ARE REJECTED: error: chol(): decomposition failed
      
      // ADD conditions on A and B:
      cond_A = true;
      cond_B = true;
      
      // check if there are NA's in the proposed values:
      for (unsigned int k = 0; k < K; ++k) {
        if( NumericVector::is_na(alpha_prop_A[k]) ){
          cond_A = false;
          Rcout << "1st control failed for A: " << std::endl << alpha_prop_A << std::endl;
        }
        if( NumericVector::is_na(alpha_prop_B[k]) ){
          cond_B = false;
          Rcout << "1st control failed for B: " << std::endl << alpha_prop_B << std::endl;
        }
      }
      
      // only if there are no NA's:
      // check that:
      // 1) All exp(alpha_prop) are strictly > 0;
      // 2) the sum of exp(alpha_prop) is < 10^10
      if(cond_A){
        if( !is_true( all( exp(alpha_prop_A) > 0 ) ) || (sum(exp(alpha_prop_A)) >  pow(10, 10) ) ){
          cond_A = false;
          Rcout << "2nd control failed for A: " << std::endl << alpha_prop_A << std::endl;
        }
      }
      if(cond_B){
        if( !is_true( all( exp(alpha_prop_B) > 0 ) ) || (sum(exp(alpha_prop_B)) >  pow(10, 10) ) ){
          cond_B = false;
          Rcout << "2nd control failed for B: " << std::endl << alpha_prop_B << std::endl;
        }
      }
      
      if(cond_A){
        ll_alpha_prop_A = ll_alpha_new( pi_new_A, exp(alpha_prop_A), N_1, K);
        ll_alpha_new_A  = ll_alpha_new( pi_new_A, exp(alpha_new_A),  N_1, K);
        
        // ACCEPT/REJECT proposal for group A:
        alpha = std::min( ll_alpha_prop_A - ll_alpha_new_A + prior_alpha_prop_A - prior_alpha_new_A, 0.0);
        if( !NumericVector::is_na(alpha) ){ // check no NA's in alpha.
          cond_01 = Rcpp::rbinom(1, 1, exp(alpha))(0);
          
          if(cond_01 == 1){ // # Update parameter and prior
            alpha_new_A = alpha_prop_A;
            prior_alpha_new_A = prior_alpha_prop_A;
            ll_alpha_new_A = ll_alpha_prop_A;
          }
        }
      }
      
      if(cond_B){
        ll_alpha_prop_B = ll_alpha_new(pi_new_B, exp(alpha_prop_B), N_2, K);
        ll_alpha_new_B  = ll_alpha_new(pi_new_B, exp(alpha_new_B),  N_2, K);
        
        // ACCEPT/REJECT proposal for group B:
        alpha = std::min( ll_alpha_prop_B - ll_alpha_new_B + prior_alpha_prop_B - prior_alpha_new_B, 0.0);
        if( !NumericVector::is_na(alpha) ){ // check no NA's in alpha.
          cond_01 = Rcpp::rbinom(1, 1, exp(alpha))(0);
          
          if(cond_01 == 1){ // # Update parameter and prior
            alpha_new_B = alpha_prop_B;
            prior_alpha_new_B = prior_alpha_prop_B;
            ll_alpha_new_B = ll_alpha_prop_B;
          }
        }
      }
      
      // IMPORTANT:
      // LIKELIHOOD: put a control when res is nan and -inf, otherwise all values accepted!
      // DONE IT: it's in alpha (prob of accepting the new value).
      
      // update mcmc matrixes
      mcmc_alpha_A(r,_) = alpha_new_A;
      mcmc_alpha_B(r,_) = alpha_new_B;
      //      for (unsigned int k = 0; k < K; ++k) {
      //        mcmc_alpha_A(r,k) = alpha_new_A[k];
      //        mcmc_alpha_B(r,k) = alpha_new_B[k];
      //      }
      ll[r] = prior_alpha_new_A + prior_alpha_new_B + ll_alpha_new_A + ll_alpha_new_B;
    }
    
    Rcpp::NumericVector row(K);
    Rcpp::NumericVector precision_A(R-burn_in); //precision for group A
    Rcpp::NumericVector precision_B(R-burn_in); //precision for group B
    for (unsigned int r = burn_in; r < R; ++r) {
      row = Rcpp::exp(mcmc_alpha_A(r,_)); // exp to gain the alpha's
      precision_A(r-burn_in) = log(sum(row)); //precision for group A (before deviding row by l)
      row = row/l; // divide by the transcript effective length
      mcmc_alpha_A(r,_) = row/sum(row); // standardize to obtain proprotions again.
      
      row = Rcpp::exp(mcmc_alpha_B(r,_)); // exp to gain the alpha's
      precision_B(r-burn_in) = log(sum(row)); //precision for group B (before deviding row by l)
      row = row/l; // divide by the transcript effective length
      mcmc_alpha_B(r,_) = row/sum(row); // standardize to obtain proprotions again.
    }
    
    // THIN HERE: return 10^4 values (1.2 * 10^4 for ll).
    return Rcpp::List::create(Rcpp::Named("A") = mcmc_alpha_A( Range(burn_in, R-1), Range(0, K-1) ),
                              Rcpp::Named("B") = mcmc_alpha_B( Range(burn_in, R-1), Range(0, K-1) ),
                              Rcpp::Named("log-posterior") = ll[ Range(burn_in, R-1) ],
                              Rcpp::Named("prec-A") = precision_A,
                              Rcpp::Named("prec-B") = precision_B);
    
  } // return the ll WITH the burn-in
