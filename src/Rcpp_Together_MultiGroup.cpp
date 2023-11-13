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

// TMG represents Together Multi-group
void centerNumericMatrix_bisTMG(Rcpp::NumericMatrix& X) {
  const unsigned int m = X.ncol();
  for (unsigned int j = 0; j < m; ++j) {
    X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
  }
}

void covRcpp_bisTMG(Rcpp::NumericMatrix& Y,
                    Rcpp::NumericMatrix& cov,
                    double const& c_diag,
                    double const& c_prop, 
                    unsigned int const& r,
                    unsigned int const& K) {
  unsigned const int df = r - 101;
  
  // Centering the matrix!
  centerNumericMatrix_bisTMG(Y);  // Defined in aux_functions
  
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
  my_rmvnorm_finalTMG(Rcpp::NumericVector& alpha_prop,
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
  prior_informativeTMG(double& prior,
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
  prior_non_informativeTMG(double& prior, 
                           Rcpp::NumericVector const& x,
                           unsigned int const& K){
    prior = 0.0;
    
    for (unsigned int k = 0; k < K; ++k) {
      prior += R::dnorm(x[k], 0, 10, true);
    }
  }

double
  ll_alpha_newTMG(Rcpp::NumericMatrix const& pi, 
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
  Rcpp_FULL_Together_Multigroup(unsigned int& R,
                                unsigned int const& burn_in,
                                Rcpp::IntegerVector const& N,
                                unsigned int const& N_groups,
                                double const& mean_log_delta,
                                double const& sd_log_delta,
                                Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericMatrix>>& pi_new,
                                Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericMatrix>>& mcmc_alpha,
                                Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericVector>>& alpha_new,
                                Rcpp::ListOf<Rcpp::ListOf<Rcpp::NumericMatrix>>& chol,
                                Rcpp::ListOf<Rcpp::NumericMatrix> TOT_Y_new,
                                Rcpp::ListOf<Rcpp::NumericMatrix> precision,
                                Rcpp::IntegerVector const& K, // number of transcripts per gene
                                Rcpp::NumericVector const& l, // double vector of len K
                                Rcpp::ListOf<Rcpp::IntegerMatrix> const& f,
                                Rcpp::IntegerMatrix const& exon_id, // integer matrix (K * J)
                                Rcpp::LogicalVector const& One_transcript, // bool vector of len n_genes, it indicates if the gene has 1 transcript only
                                Rcpp::LogicalVector const& one_transcript // bool vector of len J, it indicates if the class has 1 transcript only
  ){
    // N_groups defines the number of groups to compare.
    
    // define the number of classes
    unsigned int J = exon_id.ncol();
    unsigned int K_tot = exon_id.nrow();
    unsigned int n_genes = K.length();
    unsigned int k_;
    
    Rcpp::NumericVector alpha;
    Rcpp::IntegerVector Y_new; // I initialize the matrix of counts.
    Rcpp::NumericVector dir_sample; // dirichlet sample.
    
    // the proportionality constant for the ARW matrix.
    Rcpp::NumericVector c_prop(n_genes);
    for (unsigned int g = 0; g < n_genes; ++g) {
      if(One_transcript[g] == false){ // if >1 transcript in the gene.  
        c_prop[g] = 0.6/K[g];
      }
    }
    
    
    // initialize a few objects once only:
    bool cond;
    unsigned int cond_01;
    double prob_tot;
    
    Rcpp::NumericVector pi_new_i(K_tot);
    
    // Initialize on top:
    Rcpp::NumericVector alpha_prop;
    double dir_sample_sum;
    double prior_alpha_prop = 0.0;
    Rcpp::NumericMatrix prior_alpha_new(N_groups, n_genes);
    double ll_alpha_prop = 0.0;
    double ll_alpha_new_ = 0.0;
    double acc_prob;
    
    // store the log-posterior of the hyper-parameters | rest
    Rcpp::NumericVector ll(R);
    Rcpp::IntegerMatrix X_new(J,K_tot);
    
    // Rcpp::NumericVector dir_sample(K);
    Rcpp::NumericMatrix cv_alpha;
    Rcpp::NumericVector prob(K_tot);
    Rcpp::IntegerVector n_multinom(K_tot);
    
    Rcpp::NumericMatrix Y;
    
    // Initialize the prior:
    for (unsigned int n = 0; n < N_groups; ++n) { // n representing the group id
      for (unsigned int g = 0; g < n_genes; ++g) {
        if(One_transcript[g] == false){ // if >1 transcript in the gene.  
          
          alpha_prop = alpha_new[n][g];
          if( mean_log_delta != 0){ // if prior has been specified
            prior_informativeTMG( prior_alpha_prop, alpha_prop,  mean_log_delta, sd_log_delta, K[g] );
          }
          else{ // if prior has NOT been specified: TO DO: compute prior once only at the beginning!
            prior_non_informativeTMG( prior_alpha_prop, alpha_prop, K[g] );
          }
          prior_alpha_new(n,g) = prior_alpha_prop;
        }
      }
    }
    
    // constant to add to the diagonal of the ARW matrix.
    double c_diag = 0.001;
    
    // MCMC iterator:
    for (unsigned int r = 0; r < R; ++r) {
      
      for (unsigned int n = 0; n < N_groups; ++n) { // n representing the group id
        for (int i = 0; i < N[n]; ++i) { // loop on all N[n] samples of group n
          k_ = 0;
          for (unsigned int g = 0; g < n_genes; ++g) { // loop on n_genes genes
            for (int k = 0; k < K[g]; ++k) {
              pi_new_i[k_] = TOT_Y_new[n](i,g) * pi_new[n][g](i,k) / l[k_];
              k_ += 1;
            }
          }
          
          for (unsigned int j = 0; j < J; ++j) {
            if (one_transcript[j] or f[n](j, i) == 0) {
              for (unsigned int k = 0; k < K_tot; ++k) {
                X_new(j,k) = f[n](j, i) * exon_id(k,j);
              }
            }
            else {
              // compute (elementwise): prob = pi_new_i * exon_id_col_j
              for (unsigned int k = 0; k < K_tot; ++k) {
                prob[k] = pi_new_i[k] * exon_id(k,j);
              }
              // normalize elements of prob to sum to 1
              prob_tot = std::accumulate(prob.begin(), prob.end(), 0.0);
              
              if( prob_tot > 0.0){
                for (unsigned int k = 0; k < K_tot; ++k) {
                  prob[k] /= prob_tot;
                }
                rmultinom(f[n](j, i), prob.begin(), K_tot, n_multinom.begin());
                //          gsl_ran_multinomial(r, K, f(j, i), prob.begin(), (unsigned int *) n.begin());
                
                for (unsigned int k = 0; k < K_tot; ++k) {
                  //  prob.push_back(pi_new_i[k] * exon_id_col_j[k]);
                  X_new(j,k) = n_multinom[k];
                }
              }
              else { // if prob == 0, then set the X counts to 0
                for (unsigned int k = 0; k < K_tot; ++k) {
                  X_new(j,k) = 0; // maybe not needed ? already 0 ? double-check!
                }
              }
            }
          }
          
          // Dirichlet sampling:
          
          k_=0;
          // add gene_id START!
          for (unsigned int g = 0; g < n_genes; ++g) {
            Y_new = Rcpp::IntegerVector(K[g]);
            
            TOT_Y_new[n](i,g) = 0;
            
            for (int k = 0; k < K[g]; ++k) { // g is needed to loop over the correct k
              for (unsigned int j = 0; j < J; ++j) { // sum all values over the corresponding equivalence classes.
                Y_new(k) += X_new(j,k_);
              }
              k_ += 1;
              
              TOT_Y_new[n](i,g) += Y_new(k); // total per sample and gene, added across all transcripts (k = 0 to K[g]) of all classes (j = 1 to J)
            }
            // add gene_id END!
            
            if(One_transcript[g] == false){ // if >1 transcript in the gene.
              alpha = Rcpp::NumericVector(K[g]);
              for (int k = 0; k < K[g]; ++k) { // g is needed to loop over the correct k
                alpha[k] = exp(alpha_new[n][g][k]) + Y_new(k);
              }
              
              // Dirichlter sampling
              // sample from a gamma (MAKE SURE RATE/SHAPE is correct!)
              dir_sample = Rcpp::NumericVector(K[g]); // dirichlet sample.
              
              dir_sample_sum = 0;
              
              for (int k = 0; k < K[g]; ++k) {
                dir_sample[k] = as<double>(Rcpp::rgamma(1, alpha[k], 1));
                dir_sample_sum += dir_sample[k];
              }
              
              // check that: all( pi's > 10^{-100} )
              cond = true;
              for (int k = 0; k < K[g]; ++k) {
                dir_sample[k] = dir_sample[k]/dir_sample_sum;
                
                if( ( Rcpp::NumericVector::is_na(dir_sample[k]) ) || ( dir_sample[k] < pow(10, -100) )  ) {
                  cond = false;
                }
              }
              
              if(cond){
                for (int k = 0; k < K[g]; ++k) {
                  pi_new[n][g](i, k) = dir_sample[k];
                }
              }
            }
          }
        } // end of pi & X Gibbs sampling.
        
        // Start of alpha Metropolis sampling:
        for (unsigned int g = 0; g < n_genes; ++g) {
          if(One_transcript[g] == false){ // if >1 transcript in the gene.
            
            if( r < 200 ){ // simple R.W. for the initial 200 iter.
              alpha_prop = Rcpp::rnorm(K[g], 0, 0.1);
              alpha_prop += alpha_new[n][g];
            }
            else{ // ARW for the following iterations:
              if( (r == 200) || (r == burn_in) ){ //round( r/100) == r/100 ){ // update the covariance matrix every 100 iterations
                cv_alpha = Rcpp::NumericMatrix(K[g], K[g]);
                
                // maybe I need to use arma:: mat to do element-wise summation
                Y = mcmc_alpha[g][n]( Range(100, r-1), Range(0, K[g]-1) ); // SUB-SET of mcmc, mcmc[101:{r-2},1:K]
                covRcpp_bisTMG(Y, cv_alpha, c_diag, c_prop[g], r, K[g]);
                chol[n][g] = Rcpp::wrap( arma::chol(Rcpp::as<arma::mat>(cv_alpha)) );
                
                //Rcout << "chol_A is" << std::endl << chol_A << std::endl;
                //Rcout << "chol_B is" << std::endl << chol_B << std::endl;
              }
              // modify here to compute the cholesky decomposition above.
              my_rmvnorm_finalTMG(alpha_prop, Rcpp::as<arma::mat>(chol[n][g]), alpha_new[n][g], K[g]);
            }
            
            if( mean_log_delta != 0){ // if prior has been specified
              prior_informativeTMG( prior_alpha_prop, alpha_prop,  mean_log_delta, sd_log_delta, K[g]  );
            }
            else{ // if prior has NOT been specified: TO DO: compute prior once only at the beginning!
              prior_non_informativeTMG( prior_alpha_prop, alpha_prop, K[g]  );
            }
            
            // WHEN ALL VALUES ARE REJECTED: error: chol(): decomposition failed
            
            // ADD conditions on A and B:
            cond = true;
            
            // check if there are NA's in the proposed values:
            for (int k = 0; k < K[g]; ++k) {
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
              ll_alpha_prop = ll_alpha_newTMG( pi_new[n][g], exp(alpha_prop), N[n], K[g]);
              ll_alpha_new_ = ll_alpha_newTMG( pi_new[n][g], exp(alpha_new[n][g]),  N[n], K[g]);
              
              // ACCEPT/REJECT proposal:
              acc_prob = std::min( ll_alpha_prop - ll_alpha_new_ + prior_alpha_prop - prior_alpha_new(n,g), 0.0);
              if( !NumericVector::is_na(acc_prob) ){ // check no NA's in alpha.
                cond_01 = Rcpp::rbinom(1, 1, exp(acc_prob))(0);
                
                if(cond_01 == 1){ // # Update parameter and prior
                  alpha_new[n][g] = alpha_prop;
                  prior_alpha_new(n,g) = prior_alpha_prop;
                  ll_alpha_new_ = ll_alpha_prop;
                }
              }
            }
            
            // IMPORTANT:
            // LIKELIHOOD: put a control when res is nan and -inf, otherwise all values accepted!
            // DONE IT: it's in alpha (prob of accepting the new value).
            
            // update mcmc matrixes
            mcmc_alpha[g][n](r,_) = alpha_new[n][g];
            //      for (unsigned int k = 0; k < K; ++k) {
            //        mcmc_alpha_A(r,k) = alpha_new_A[k];
            //        mcmc_alpha_B(r,k) = alpha_new_B[k];
            //      }
            ll[r] += prior_alpha_new(n,g) + ll_alpha_new_;
          }
        }
      }
    }
    
    Rcpp::IntegerVector cumsum_K = cumsum(K);
    Rcpp::NumericVector l_correct;
    Rcpp::NumericVector row;
    
    for (unsigned int g = 0; g < n_genes; ++g) {
      if(g == 0){
        l_correct = l[ Range(0, K[0]-1) ];
      }else{
        l_correct = l[ Range(cumsum_K[g-1]-1, cumsum_K[g]-1) ];
      }
      
      if(One_transcript[g] == false){ // if >1 transcript in the gene.
        row = Rcpp::NumericVector(K[g]);
        for (unsigned int n = 0; n < N_groups; ++n) { // n representing the group id
          
          for (unsigned int r = burn_in; r < R; ++r) {
            row = Rcpp::exp(mcmc_alpha[g][n](r,_)); // exp to gain the alpha's
            precision[g](r-burn_in, n) = log(sum(row)); //precision for group A (before deviding row by l)
            // row = row/sum(row); // divide by their sum to get pi_bar
            row = row/l_correct; // divide by the transcript effective length
            mcmc_alpha[g][n](r,_) = row/sum(row); // standardize to obtain proprotions again.
          }
          
          mcmc_alpha[g][n] = mcmc_alpha[g][n]( Range(burn_in, R-1), Range(0, K[g]-1) );
        }
      }
    }
    
    // THIN HERE: return 10^4 values (1.2 * 10^4 for ll).
    return Rcpp::List::create(Rcpp::Named("mcmc") = mcmc_alpha,
                              Rcpp::Named("log-posterior") = ll[ Range(burn_in, R-1) ],
                              Rcpp::Named("log-precision") = precision);
  } // return the ll WITH the burn-in
// the original element pi_new is modified
