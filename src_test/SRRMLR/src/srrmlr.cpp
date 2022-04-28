#include "functions.h"
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


//' The base function of SRRMLR
//'
//' @param X the covariate matrix
//' @param Y the response matrix
//' @param s the sparsity
//' @param r the rank
//' @param e_in the stopping epsilon of inloop
//' @param e_out the stopping epsilon of outloop
//' @return C matrix aims to be solved
//' @return iters iterations taken when algorithm meet the stopping criterion 
//' @return Prob the probability matrix
//' @export 
//' @useDynLib SRRMLR
// [[Rcpp::export]]
Rcpp::List srrmlr_one(arma::dmat &X, arma::dmat &Y, 
                      arma::uword s, arma::uword r, 
                      double e_in=0.001, double e_out=0.001){
  if(r>s){
    return 0;
  }
  const arma::uword n = X.n_rows, p = X.n_cols, q = Y.n_cols;
  arma::dmat B(p, r, arma::fill::zeros);
  arma::dmat V(q-1, r, arma::fill::eye);
  arma::dmat pre_C(p, q-1), pr(n, q-1), M(n, q-1), Ga(p, r), obj(q-1, q-1), C(p, q-1);
  arma::uvec indices, ord, ordset, aset, iset, iset_;
  arma::dvec dB(p);
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::uword iter;
  double med;
  
  for(iter=0; iter<1000; iter++){
    pre_C = B*trans(V);
    for(arma::uword _=0; _<p; _++){dB(_) = arma::norm(B.row(_), 2);}
    iset_ = arma::find(dB==0);
    
    for(arma::uword k=0; k<n; k++){
      med = accu(exp(trans(pre_C) * trans(X.row(k))));
      for(arma::uword _=0; _<q-1; _++){
        pr(k,_) = accu(exp(trans(pre_C.col(_))*trans(X.row(k))))/(1+med);
      }
    }
    M = 2*(Y.cols(1,q-1)-pr)+X*pre_C;
    Ga.zeros();
    Ga.rows(iset_) = trans(X.cols(iset_))*(M*V-X*B) / n;
    
    for(arma::uword in_loop=0; in_loop<500; in_loop++){
      C = B*trans(V);
      for(arma::uword _=0; _<p; _++){dB(_) = arma::norm(B.row(_)+Ga.row(_), 2);}
      indices = arma::sort_index(dB, "descend");
      aset = arma::sort(indices.subvec(0, s-1));
      iset = arma::sort(indices.subvec(s, p-1));
      obj = trans(M)*X.cols(aset)*inv(trans(X.cols(aset))*X.cols(aset))*trans(X.cols(aset))*M;
      eig_gen(eigval, eigvec, obj);
      ord = arma::stable_sort_index(real(eigval), "descend");
      ordset = ord.subvec(0, r-1);
      V = real(eigvec.cols(ordset));
      B.zeros(); Ga.zeros();
      B.rows(aset) = inv(trans(X.cols(aset))*X.cols(aset))*trans(X.cols(aset))*M*V;
      Ga.rows(iset) = trans(X.cols(iset))*(M*V-X*B)/n;
      if(accu(square(B*trans(V)-C))<e_in) break;
    }
    if(iter>50){
      if(accu(square(pre_C-C))<e_out){
        break;
      }
    }
  }
  arma::dmat prob = arma::join_rows(arma::ones(n,1)-arma::sum(pr, 1), pr);
  C = B*trans(V);
  
  Rcpp::List out;
  out["C"] = C;
  out["iters"] = iter;
  out["Prob"] = prob;
  return out;
}



//' SRRMLR in simulation analysis
//'
//' @param train_X the covariate matrix of training data
//' @param train_Y the response matrix of training data
//' @param val_X the covariate matrix of validation data
//' @param val_Y the response matrix of validation data
//' @param test_X the covariate matrix of test data
//' @param test_Y the response matrix of test data
//' @param s_list the list of sparsity
//' @param r_list the list of rank
//' @param C_true the true matrix C
//' @param s_true the true sparsity
//' @param r_true the true rank
//' @return Pred the predictive accuracy with the negative log-likelihood
//' @return Est the estimation accuracy
//' @return r rank selected by SRRMLR, the estimated rank
//' @return s sparsity selected by SRRMLR, also the number of non-zero rows in C_hat
//' @return Sen the sensitivity
//' @return Spe the specificity
//' @return time running time of the algorithm with selected rank and sparsity
//' @export 
//' @useDynLib SRRMLR
// [[Rcpp::export]]
Rcpp::List srrmlr_simu(arma::dmat &train_X, arma::dmat &train_Y,
                       arma::dmat &val_X, arma::dmat &val_Y,
                       arma::dmat &test_X, arma::dmat &test_Y,
                       arma::uvec &s_list, arma::uvec &r_list,
                       arma::dmat &C_true,
                       arma::uword s_true, arma::uword r_true){
  const arma::uword rC = train_X.n_cols, cC = train_Y.n_cols;
  arma::uword now_s, now_r, best_s, best_r;
  
  double best_time, best_Pred=INFINITY, best_Est, best_Sen, best_Spe;
  arma::dmat best_C(rC, cC-1);

  for(arma::uword i=0; i<r_list.n_elem; i++){
    for(arma::uword j=0; j<s_list.n_elem; j++){
      now_s = s_list(j);
      now_r = r_list(i);
      if(now_r>now_s or now_r>=cC or now_s>=rC){
        continue;
      }
      SRR test = srrmlr_in(train_X, train_Y, val_X, val_Y, now_s, now_r);
      if(test.Pred<best_Pred){
        best_r = now_r;
        best_s = now_s;
        best_time = test.Stime;
        best_Pred = test.Pred;
        best_C = test.get_C;
      }
    }
  }
  arma::dvec temp_aset = arma::sum(best_C, 1);
  arma::uvec final_aset = arma::find(temp_aset!=0);
  arma::uvec final_iset = arma::find(temp_aset==0);
  
  temp_aset = arma::sum(C_true, 1);
  arma::uvec Pset = arma::find(temp_aset!=0);
  arma::uvec Nset = arma::find(temp_aset==0);

  best_Est = accu(square(best_C-C_true.cols(1,cC-1)))/(rC*(cC-1));
  best_Pred = compute_error(test_Y, test_X, best_C);
  arma::uvec Sen_temp = arma::intersect(final_aset, Pset);
  arma::uvec Spe_temp = arma::intersect(final_iset, Nset);
  double Sens = Sen_temp.n_elem, Spes = Spe_temp.n_elem;
  best_Sen = Sens/s_true;
  best_Spe = Spes/(rC-s_true);

  Rcpp::List out;
  out["Pred"] = best_Pred;
  out["Est"] = best_Est;
  out["r"] = best_r;
  out["s"] = best_s;
  out["Sen"] = best_Sen;
  out["Spe"] = best_Spe;
  out["time"] = best_time;
  return out;
}



//' SRRMLR in real data analysis
//'
//' @param train_X the covariate matrix of training data
//' @param train_Y the response matrix of training data
//' @param test_X the covariate matrix of test data
//' @param test_Y the response matrix of test data
//' @param s_list the list of sparsity
//' @param r_list the list of rank
//' @param cv the number of folds when conducting cross-validation, 5 in default
//' @return Pred the predictive accuracy with the negative log-likelihood
//' @return r rank selected by SRRMLR, the estimated rank
//' @return s sparsity selected by SRRMLR, also the number of non-zero rows in C_hat
//' @return A_set the index of features selected by SRRMLR.
//' @return time running time of the algorithm with selected rank and sparsity
//' @export 
//' @useDynLib SRRMLR
// [[Rcpp::export]]
Rcpp::List srrmlr_real(arma::dmat &train_X, arma::dmat &train_Y,
                       arma::dmat &test_X, arma::dmat &test_Y,
                       arma::uvec &s_list, arma::uvec &r_list,
                       arma::uword cv=5){
  const arma::uword rX = train_X.n_rows, rC = train_X.n_cols, cC = train_Y.n_cols;
  const arma::uword split_size=ceil(train_X.n_rows/cv);
  arma::uword temp, start, end, now_s, now_r, best_s, best_r;
  double best_time, best_Pred=INFINITY;
  arma::dmat best_C(rC, cC-1);
  SRR test;
  for(arma::uword i=0; i<r_list.n_elem; i++){
    for(arma::uword j=0; j<s_list.n_elem; j++){
      now_s = s_list(j);
      now_r = r_list(i);
      if(now_r>now_s or now_r>=cC or now_s>=rC){
        continue;
      }
      if(cv>1){
        double temp_sum=0.0;
        for(arma::uword spt=0; spt<cv; spt++){
          start = spt*split_size;
          temp = (spt+1)*split_size;
          end= (temp>=rX)*rX + (temp<rX)*temp -1;
          arma::dmat XX = train_X;
          arma::dmat YY = train_Y;
          XX.shed_rows(start, end);
          YY.shed_rows(start, end);
          test = srrmlr_in(XX, YY, 
                           train_X.rows(start, end),
                           train_Y.rows(start, end),
                           now_s, now_r);
          temp_sum += test.Pred;
        }
        if(temp_sum<best_Pred){
          best_r = now_r;
          best_s = now_s;
          best_Pred = temp_sum;
        }
      }
      else{
        SRR test = srrmlr_in(train_X, train_Y, train_X, train_Y, now_s, now_r);
        if(test.Pred<best_Pred){
          best_r = now_r;
          best_s = now_s;
          best_Pred = test.Pred;
        }
      }
    }
  }
  
  test = srrmlr_in(train_X, train_Y, train_X, train_Y, best_s, best_r);
  best_time = test.Stime;
  best_Pred = test.Pred;
  best_C = test.get_C;
  arma::dvec temp_aset = arma::sum(best_C, 1);
  arma::uvec final_aset = arma::find(temp_aset!=0);
  arma::uvec final_iset = arma::find(temp_aset==0);

  best_Pred = compute_error(test_Y, test_X, best_C);

  Rcpp::List out;
  out["Pred"] = best_Pred;
  out["r"] = best_r;
  out["s"] = best_s;
  out["A_set"] = final_aset;
  out["time"] = best_time;
  return out;
}