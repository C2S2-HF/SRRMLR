#ifndef SRC_FUNCTIONS_H
#define SRC_FUNCTIONS_H


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

inline double compute_error(const arma::dmat err_Y, const arma::dmat err_X, const arma::dmat err_C){
  arma::uword err_q = err_Y.n_cols;
  arma::dmat temp = err_X*err_C;
  arma::dvec temp1 = arma::sum(exp(temp), 1);
  temp1++;
  return accu(log(temp1))-arma::trace(temp*trans(err_Y.cols(1,err_q-1)));
}

typedef struct{
  arma::dmat get_C;
  double Pred;
  double Stime;
  arma::uword Ster;
}SRR;


SRR srrmlr_in(const arma::dmat XX, const arma::dmat YY,
              const arma::dmat val_XX, const arma::dmat val_YY,
              const arma::uword ss, const arma::uword rr,
              const double e_in=0.001, const double e_out=0.001){
  arma::wall_clock timer; timer.tic();
  const arma::uword n = XX.n_rows, p = XX.n_cols, q = YY.n_cols;
  
  arma::dmat B(p, rr, arma::fill::zeros);
  arma::dmat V(q-1, rr, arma::fill::eye);
  arma::dmat pre_C(p, q-1), pr(n, q-1), M(n, q-1), Ga(p, rr), obj(q-1, q-1), C(p, q-1);
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
      med = accu(exp(trans(pre_C) * trans(XX.row(k))));
      for(arma::uword _=0; _<q-1; _++){
        pr(k,_) = accu(exp(trans(pre_C.col(_))*trans(XX.row(k))))/(1+med);
      }
    }
    M = 2*(YY.cols(1,q-1)-pr)+XX*pre_C;
    Ga.zeros();
    Ga.rows(iset_) = trans(XX.cols(iset_))*(M*V-XX*B) / n;
    
    for(arma::uword in_loop=0; in_loop<500; in_loop++){
      C = B*trans(V);
      for(arma::uword _=0; _<p; _++){dB(_) = arma::norm(B.row(_)+Ga.row(_), 2);}
      indices = arma::sort_index(dB, "descend");
      aset = arma::sort(indices.subvec(0, ss-1));
      iset = arma::sort(indices.subvec(ss, p-1));
      obj = trans(M)*XX.cols(aset)*inv(trans(XX.cols(aset))*XX.cols(aset))*trans(XX.cols(aset))*M;
      eig_gen(eigval, eigvec, obj);
      ord = arma::stable_sort_index(real(eigval), "descend");
      ordset = ord.subvec(0, rr-1);
      V = real(eigvec.cols(ordset));
      B.zeros(); Ga.zeros();
      B.rows(aset) = inv(trans(XX.cols(aset))*XX.cols(aset))*trans(XX.cols(aset))*M*V;
      Ga.rows(iset) = trans(XX.cols(iset))*(M*V-XX*B)/n;
      if(accu(square(B*trans(V)-C))<e_in) break;
    }
    if(iter>50){
      if(accu(square(pre_C-C))<e_out){
        break;
      }
    }
  }
  
  SRR tmp;
  if(iter>=990){
    tmp.get_C = B*trans(V);
    tmp.Pred = INFINITY;
    tmp.Ster = iter;
    tmp.Stime = timer.toc();
  }
  else{
    tmp.get_C = B*trans(V);
    tmp.Pred = compute_error(val_YY, val_XX, tmp.get_C);
    tmp.Ster = iter;
    tmp.Stime = timer.toc();
  }
  return tmp;
}

#endif
