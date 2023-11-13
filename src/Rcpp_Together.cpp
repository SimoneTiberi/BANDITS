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

void centerNumericMatrix_bis_Tog(Rcpp::NumericMatrix& X) {
  const unsigned int m = X.ncol();
  for (unsigned int j = 0; j < m; ++j) {
    X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
  }
}

void covRcpp_bis_Tog(Rcpp::NumericMatrix& Y,
                     Rcpp::NumericMatrix& cov,
                     double const& c_diag,
                     double const& c_prop, 
                     unsigned int const& r,
                     unsigned int const& K) {
  unsigned const int df = r - 101;
  
  // Centering the matrix!
  centerNumericMatrix_bis_Tog(Y);  // Defined in aux_functions
  
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

double
  ll_alpha_new_Tog(Rcpp::NumericMatrix const& pi, 
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

void
  my_rmvnorm_final_Tog(Rcpp::NumericVector& alpha_prop,
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
  prior_informative_Tog(double& prior,
                        Rcpp::NumericVector& x,
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
  prior_non_informative_Tog(double& prior, 
                            Rcpp::NumericVector& x,
                            unsigned int const& K){
    prior = 0.0;
    
    for (unsigned int k = 0; k < K; ++k) {
      prior += R::dnorm(x[k], 0, 10, true);
    }
    
  }

// [[Rcpp::export]]
Rcpp::List
  Rcpp_Together(  unsigned int& R,
                  unsigned int const& burn_in,
                  unsigned int const& N_1,
                  unsigned int const& N_2,
                  Rcpp::ListOf<Rcpp::NumericMatrix>& pi_new_A,
                  Rcpp::ListOf<Rcpp::NumericMatrix>& pi_new_B,
                  Rcpp::ListOf<Rcpp::NumericMatrix>& mcmc_alpha_A,
                  Rcpp::ListOf<Rcpp::NumericMatrix>& mcmc_alpha_B,
                  Rcpp::ListOf<Rcpp::NumericVector>& log_alpha_new_A,
                  Rcpp::ListOf<Rcpp::NumericVector>& log_alpha_new_B,
                  Rcpp::ListOf<Rcpp::NumericMatrix>& chol_A,
                  Rcpp::ListOf<Rcpp::NumericMatrix>& chol_B,
                  double const& mean_log_delta,
                  double const& sd_log_delta,
                  Rcpp::IntegerVector const& K, // number of transcripts per gene
                  Rcpp::NumericVector const& l, // list of n_genes elements, each with a double vector of len K[n]
                  Rcpp::IntegerMatrix const& f, // integer matrix (J * 2*N)
                  Rcpp::IntegerMatrix const& exon_id, // integer matrix (K * J)
                  Rcpp::LogicalVector const& One_transcript, // bool vector of len n_genes, it indicates if the gene has 1 transcript only
                  Rcpp::LogicalVector const& one_transcript // bool vector of len J, it indicates if the class has 1 transcript only
  ){
    // initialize a few objects once only:
    double prob_tot;
    unsigned int k_;
    
    unsigned int n_genes = K.length();
    
    R += burn_in;
    unsigned int J = exon_id.ncol();
    unsigned int K_tot = exon_id.nrow();
    
    Rcpp::NumericMatrix TOT_Y_new_A(N_1, n_genes);
    Rcpp::NumericMatrix TOT_Y_new_B(N_2, n_genes);
    std::fill(TOT_Y_new_A.begin(), TOT_Y_new_A.end(), 1.0);
    std::fill(TOT_Y_new_B.begin(), TOT_Y_new_B.end(), 1.0);
    
    Rcpp::NumericVector c_prop(n_genes);
    
    Rcpp::NumericVector prior_alpha_new_A(n_genes);
    Rcpp::NumericVector prior_alpha_new_B(n_genes);
    
    double prior_alpha_prop_A;
    double prior_alpha_prop_B;
    Rcpp::NumericVector alpha_prop_A;
    Rcpp::NumericVector alpha_prop_B;
    
    for (unsigned int g = 0; g < n_genes; ++g) {
      if(One_transcript[g] == false){ // if >1 transcript in the gene.  
        c_prop[g] = 0.6/K[g];
        
        alpha_prop_A = log_alpha_new_A[g];
        // prior initialization:
        if( mean_log_delta != 0){ // if prior has been specified
          prior_informative_Tog( prior_alpha_prop_A, alpha_prop_A, mean_log_delta, sd_log_delta, K[g]  );
        }
        else{ // if prior has NOT been specified:
          prior_non_informative_Tog( prior_alpha_prop_A, alpha_prop_A, K[g]  );
        }
        
        prior_alpha_new_A[g] = prior_alpha_prop_A; //same starting value.
        prior_alpha_new_B[g] = prior_alpha_prop_A;
      }
    }
    
    bool cond_A;
    bool cond_B;
    double dir_sample_sum;
    
    Rcpp::NumericVector pi_all(K_tot);
    Rcpp::NumericVector prob(K_tot);
    Rcpp::IntegerVector n(K_tot);
    // Matrix of J rows & K_tot columns (filled with 0)
    Rcpp::IntegerMatrix X_new(J,K_tot);
    
    // CAN these elements change size during the MCMC ??? (according to K[g])
    Rcpp::NumericMatrix cv_alpha;
    
    double ll_alpha_prop_A = 0.0;
    double ll_alpha_new_A = 0.0;
    double ll_alpha_prop_B = 0.0;
    double ll_alpha_new_B = 0.0;
    
    unsigned int cond_01;
    
    double alpha;
    
    Rcpp::NumericVector ll(R);
    
    Rcpp::NumericMatrix Y;
    
    // constant to add to the diagonal of the ARW matrix.
    double c_diag = 0.001;
    
    // MCMC iterator:
    for (unsigned int r = 0; r < R; ++r) {
      
      // Group A: X and pi sampling
      for (unsigned int i = 0; i < N_1; ++i) {
        // compute `pi_new_i` element-wise, un-normalized probs
        k_ = 0;
        for (unsigned int g = 0; g < n_genes; ++g) {
          for (int k = 0; k < K[g]; ++k) {
            pi_all[k_] = TOT_Y_new_A(i,g) * pi_new_A[g](i,k) / l[k_];
            k_ += 1;
          }
        }
        
        for (unsigned int j = 0; j < J; ++j) {
          // Group A:
          if (one_transcript[j] or f(j, i) == 0) { // if the class has a single transcript
            for (unsigned int k = 0; k < K_tot; ++k) {
              X_new(j,k) = f(j, i) * exon_id(k,j);
            }
          }
          else { // if the class has > 1 transcripts
            // compute (elementwise): prob[k] =  pi_all[k] * exon_id(k,j);
            for (unsigned int k = 0; k < K_tot; ++k) {
              prob[k] = pi_all[k] * exon_id(k,j);
            }
            // normalize elements of prob to sum to 1
            prob_tot = std::accumulate(prob.begin(), prob.end(), 0.0);
            // if prob > 0, then sample from the multinomial
            if( prob_tot > 0.0){
              for (unsigned int k = 0; k < K_tot; ++k) {
                prob[k] /= prob_tot;
              }
              
              // multinomial sampling:
              //            gsl_ran_multinomial(r, K_tot, f(j, i), prob.begin(), (unsigned int *) n.begin());
              rmultinom(f(j, i), prob.begin(), K_tot, n.begin());
              for (unsigned int k = 0; k < K_tot; ++k) {
                X_new(j,k) = n[k];
              }
            }
            else { // if prob == 0, then set the X counts to 0
              for (unsigned int k = 0; k < K_tot; ++k) {
                X_new(j,k) = 0; // maybe not needed ? already 0 ? double-check!
              }
            }
          }
        }
        
        k_=0;
        // add gene_id START!
        for (unsigned int g = 0; g < n_genes; ++g) {
          Rcpp::IntegerVector Y_new(K[g]); // I initialize the matrix of counts.
          TOT_Y_new_A(i,g) = 0;
          
          for (int k = 0; k < K[g]; ++k) { // g is needed to loop over the correct k
            for (unsigned int j = 0; j < J; ++j) { // sum all values over the corresponding equivalence classes.
              Y_new(k) += X_new(j,k_);
            }
            k_ += 1;
            
            TOT_Y_new_A(i,g) += Y_new(k); // total per sample and gene, added across all transcripts (k = 0 to K[g]) of all classes (j = 1 to J)
          }
          // add gene_id END!
          
          if(One_transcript[g] == false){ // if >1 transcript in the gene.
            Rcpp::NumericVector alpha(K[g]);
            for (int k = 0; k < K[g]; ++k) { // g is needed to loop over the correct k
              alpha[k] = exp(log_alpha_new_A[g][k]) + Y_new(k);
            }
            
            // Dirichlter sampling
            // sample from a gamma (MAKE SURE RATE/SHAPE is correct!):
            Rcpp::NumericVector dir_sample(K[g]); // dirichlet sample.
            
            // Dirichlter:
            dir_sample_sum = 0;
            
            for (int k = 0; k < K[g]; ++k) {
              dir_sample[k] = as<double>(Rcpp::rgamma(1, alpha[k], 1));
              dir_sample_sum += dir_sample[k];
            }
            
            // check that: all( pi's > 10^{-100} )
            cond_A = true;
            for (int k = 0; k < K[g]; ++k) {
              dir_sample[k] = dir_sample[k]/dir_sample_sum;
              
              if( ( Rcpp::NumericVector::is_na(dir_sample[k]) ) || ( dir_sample[k] < pow(10, -100) )  ) {
                cond_A = false;
              }
            }
            
            if(cond_A){
              for (int k = 0; k < K[g]; ++k) {
                pi_new_A[g](i, k) = dir_sample[k];
              }
            }
          }
        }
      } // end of pi & X Gibbs sampling.
      
      
      
      // Group B: X and pi sampling
      for (unsigned int i = 0; i < N_2; ++i) {
        // compute `pi_new_i` element-wise, un-normalized probs
        k_ = 0;
        for (unsigned int g = 0; g < n_genes; ++g) {
          for (int k = 0; k < K[g]; ++k) {
            pi_all[k_] = TOT_Y_new_B(i,g) * pi_new_B[g](i,k) / l[k_];
            k_ += 1;
          }
        }
        
        for (unsigned int j = 0; j < J; ++j) {
          // Group B:
          if (one_transcript[j] or f(j, i+N_1) == 0) {
            for (unsigned int k = 0; k < K_tot; ++k) {
              X_new(j,k) = f(j, i+N_1) * exon_id(k,j);
            }
          }
          else {
            // compute (elementwise): prob = pi_new_i * exon_id_col_j
            for (unsigned int k = 0; k < K_tot; ++k) {
              prob[k] = pi_all[k] * exon_id(k,j);
            }
            // normalize elements of prob to sum to 1
            prob_tot = std::accumulate(prob.begin(), prob.end(), 0.0);
            // if prob > 0, then sample from the multinomial
            if( prob_tot > 0.0){
              for (unsigned int k = 0; k < K_tot; ++k) {
                prob[k] /= prob_tot;
              }
              
              // multinomial sampling:
              //            gsl_ran_multinomial(r, K_tot, f(j, i+N_1), prob.begin(), (unsigned int *) n.begin());
              rmultinom(f(j, i+N_1), prob.begin(), K_tot, n.begin());
              for (unsigned int k = 0; k < K_tot; ++k) {
                X_new(j,k) = n[k];
              }
            }
            else{ // if prob_tot == 0, then set the X[j,] counts to 0
              for (unsigned int k = 0; k < K_tot; ++k) {
                X_new(j,k) = 0;
              }
            }
          }
        }
        
        k_=0;
        // add gene_id START!
        for (unsigned int g = 0; g < n_genes; ++g) {
          Rcpp::IntegerVector Y_new(K[g]); // I initialize the matrix of counts.
          TOT_Y_new_B(i,g) = 0;
          
          for (int k = 0; k < K[g]; ++k) { // g is needed to loop over the correct k
            for (unsigned int j = 0; j < J; ++j) { // sum all values over the corresponding equivalence classes.
              Y_new(k) += X_new(j,k_);
            }
            k_ += 1;
            
            TOT_Y_new_B(i,g) += Y_new(k); // total per sample and gene, added across all transcripts (k = 0 to K[g]) of all classes (j = 1 to J)
          }
          // add gene_id END!
          
          if(One_transcript[g] == false){ // if >1 transcript in the gene.
            Rcpp::NumericVector alpha(K[g]);
            for (int k = 0; k < K[g]; ++k) { // g is needed to loop over the correct k
              alpha[k] = exp(log_alpha_new_B[g][k]) + Y_new(k);
            }
            
            // Dirichlter sampling
            // sample from a gamma (MAKE SURE RATE/SHAPE is correct!)
            Rcpp::NumericVector dir_sample(K[g]); // dirichlet sample.
            
            // Dirichlter:
            dir_sample_sum = 0;
            
            for (int k = 0; k < K[g]; ++k) {
              dir_sample[k] = as<double>(Rcpp::rgamma(1, alpha[k], 1));
              dir_sample_sum += dir_sample[k];
            }
            
            // check that: all( pi's > 10^{-100} )
            cond_B = true;
            for (int k = 0; k < K[g]; ++k) {
              dir_sample[k] = dir_sample[k]/dir_sample_sum;
              
              if( ( Rcpp::NumericVector::is_na(dir_sample[k]) ) || ( dir_sample[k] < pow(10, -100) )  ) {
                cond_B = false;
              }
            }
            
            if(cond_B){
              for (int k = 0; k < K[g]; ++k) {
                pi_new_B[g](i, k) = dir_sample[k];
              }
            }
          }
        }
      } // end of pi & X Gibbs sampling.
      
      
      
      
      // Start of alpha Metropolis sampling:
      for (unsigned int g = 0; g < n_genes; ++g) {
        if(One_transcript[g] == false){ // if >1 transcript in the gene.
          if( r < 200 ){ // simple R.W. for the initial 200 iter.
            alpha_prop_A = Rcpp::rnorm(K[g], 0, 0.1);
            alpha_prop_A += log_alpha_new_A[g];
            alpha_prop_B = Rcpp::rnorm(K[g], 0, 0.1);
            alpha_prop_B += log_alpha_new_B[g];
          }
          else{ // ARW for the following iterations:
            if( (r == 200) || (r == burn_in) ){ //round( r/100) == r/100 ){ // update the covariance matrix every 100 iterations
              // maybe I need to use arma:: mat to do element-wise summation
              // do this outside, once only!
              cv_alpha = Rcpp::NumericMatrix(K[g], K[g]);
              
              Y = mcmc_alpha_A[g]( Range(100, r-1), Range(0, K[g]-1) ); // SUB-SET of mcmc, mcmc[101:{r-2},1:K]
              covRcpp_bis_Tog(Y, cv_alpha, c_diag, c_prop[g], r, K[g]);
              chol_A[g]  = Rcpp::wrap( arma::chol(Rcpp::as<arma::mat>(cv_alpha)) );
              
              Y = mcmc_alpha_B[g]( Range(100, r-1), Range(0, K[g]-1) ); // SUB-SET of mcmc, mcmc[101:{r-2},1:K]
              covRcpp_bis_Tog(Y, cv_alpha, c_diag, c_prop[g], r, K[g]);
              chol_B[g]  = Rcpp::wrap( arma::chol(Rcpp::as<arma::mat>(cv_alpha)) );
            }
            // modify here to compute the cholesky decomposition above.
            my_rmvnorm_final_Tog(alpha_prop_A, Rcpp::as<arma::mat>(chol_A[g]), log_alpha_new_A[g], K[g]);
            my_rmvnorm_final_Tog(alpha_prop_B, Rcpp::as<arma::mat>(chol_B[g]), log_alpha_new_B[g], K[g]);
          }
          
          if( mean_log_delta != 0){ // if prior has been specified
            prior_informative_Tog( prior_alpha_prop_A, alpha_prop_A,  mean_log_delta, sd_log_delta, K[g]  );
            prior_informative_Tog( prior_alpha_prop_B, alpha_prop_B,  mean_log_delta, sd_log_delta, K[g]  );
          }
          else{ // if prior has NOT been specified: TO DO: compute prior once only at the beginning!
            prior_non_informative_Tog( prior_alpha_prop_A, alpha_prop_A, K[g]  );
            prior_non_informative_Tog( prior_alpha_prop_B, alpha_prop_B, K[g]  );
          }
          
          // ADD conditions on A and B:
          cond_A = true;
          cond_B = true;
          
          // check if there are NA's in the proposed values:
          for (int k = 0; k < K[g]; ++k) {
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
            ll_alpha_prop_A = ll_alpha_new_Tog( pi_new_A[g], exp(alpha_prop_A), N_1, K[g]);
            ll_alpha_new_A  = ll_alpha_new_Tog( pi_new_A[g], exp(log_alpha_new_A[g]),  N_1, K[g]);
            
            // ACCEPT/REJECT proposal for group A:
            alpha = std::min( ll_alpha_prop_A - ll_alpha_new_A + prior_alpha_prop_A - prior_alpha_new_A[g], 0.0);
            if( !NumericVector::is_na(alpha) ){ // check no NA's in alpha.
              cond_01 = Rcpp::rbinom(1, 1, exp(alpha))(0);
              
              if(cond_01 == 1){ // # Update parameter and prior
                log_alpha_new_A[g] = alpha_prop_A;
                prior_alpha_new_A[g] = prior_alpha_prop_A;
                ll_alpha_new_A = ll_alpha_prop_A;
              }
            }
          }
          
          if(cond_B){
            ll_alpha_prop_B = ll_alpha_new_Tog( pi_new_B[g], exp(alpha_prop_B), N_2, K[g]);
            ll_alpha_new_B  = ll_alpha_new_Tog( pi_new_B[g], exp(log_alpha_new_B[g]),  N_2, K[g]);
            
            // ACCEPT/REJECT proposal for group A:
            alpha = std::min( ll_alpha_prop_B - ll_alpha_new_B + prior_alpha_prop_B - prior_alpha_new_B[g], 0.0);
            if( !NumericVector::is_na(alpha) ){ // check no NA's in alpha.
              cond_01 = Rcpp::rbinom(1, 1, exp(alpha))(0);
              
              if(cond_01 == 1){ // # Update parameter and prior
                log_alpha_new_B[g] = alpha_prop_B;
                prior_alpha_new_B[g] = prior_alpha_prop_B;
                ll_alpha_new_B = ll_alpha_prop_B;
              }
            }
          }
          
          // update mcmc matrixes
          mcmc_alpha_A[g](r,_) = log_alpha_new_A[g];
          mcmc_alpha_B[g](r,_) = log_alpha_new_B[g];
          ll[r] += (prior_alpha_new_A[g] + prior_alpha_new_B[g] + ll_alpha_new_A + ll_alpha_new_B);
        }
      } // end of alpha Metropolis sampling.
    } // end of MCMC loop
    
    // turn alpha's into pi_T's
    Rcpp::IntegerVector cumsum_K = cumsum(K);
    Rcpp::NumericVector l_correct;
    
    Rcpp::NumericMatrix precision_A(R-burn_in, n_genes); //precision for group A
    Rcpp::NumericMatrix precision_B(R-burn_in, n_genes); //precision for group B
    
    for (unsigned int g = 0; g < n_genes; ++g) {
      if(g == 0){
        l_correct = l[ Range(0, K[0]-1) ];
      }else{
        l_correct = l[ Range(cumsum_K[g-1]-1, cumsum_K[g]-1) ];
      }
      
      if(One_transcript[g] == false){ // if >1 transcript in the gene.
        Rcpp::NumericVector row(K[g]);
        
        for (unsigned int r = burn_in; r < R; ++r) {
          row = Rcpp::exp(mcmc_alpha_A[g](r,_)); // exp to gain the alpha's
          precision_A(r-burn_in, g) = log(sum(row)); //precision for group A (before deviding row by l)
          // row = row/sum(row); // divide by their sum to get pi_bar
          row = row/l_correct; // divide by the transcript effective length
          mcmc_alpha_A[g](r,_) = row/sum(row); // standardize to obtain proprotions again.
          
          row = Rcpp::exp(mcmc_alpha_B[g](r,_)); // exp to gain the alpha's
          precision_B(r-burn_in, g) = log(sum(row)); //precision for group B (before deviding row by l)
          // row = row/sum(row); // divide by their sum to get pi_bar
          row = row/l_correct; // divide by the transcript effective length
          mcmc_alpha_B[g](r,_) = row/sum(row); // standardize to obtain proprotions again.
        }

        // I remove the burn-in here:
        mcmc_alpha_A[g] = mcmc_alpha_A[g]( Range(burn_in, R-1), Range(0, K[g]-1) );
        mcmc_alpha_B[g] = mcmc_alpha_B[g]( Range(burn_in, R-1), Range(0, K[g]-1) );
      }
    }
    
    return Rcpp::List::create(Rcpp::Named("A") = mcmc_alpha_A,
                              Rcpp::Named("B") = mcmc_alpha_B,
                              Rcpp::Named("log-posterior") = ll[ Range(burn_in, R-1) ],
                              Rcpp::Named("prec-A") = precision_A,
                              Rcpp::Named("prec-B") = precision_B);
  }
