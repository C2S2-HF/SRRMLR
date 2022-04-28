#include <stdio.h>
#include <warm_start.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::List srrmlr_warm(arma::dmat &X, arma::dmat &Y, 
                      arma::uword s, arma::uword r){

  const arma::uword n = X.n_rows, p = X.n_cols, q = Y.n_cols;
  arma::dmat C(p, q-1);
  arma::umat A(s,10000,arma::fill::none);
  SRR temp;

  arma::uword now_r=1, now_s=1, pivot=1;
  temp = get_init(X, Y, s);
  while(true){
    
    if(now_r==r && now_s==s){
      temp = srrmlr_one(X,Y,temp.init_B, temp.init_V,now_s,now_r);
      break;
    }
    
    if(pivot%2==0){
      if(now_r!=r){
        temp = srrmlr_mode2(X,Y,temp.init_B,temp.init_V,now_s,now_r,s);
        now_r+=1;
        pivot+=1;
      }else{
        temp = srrmlr_mode1(X,Y,temp.init_B,temp.init_V,now_s,now_r,s);
        now_s+=1;
      }
    }else{
      if(now_s!=s){
        temp = srrmlr_mode1(X,Y,temp.init_B,temp.init_V,now_s,now_r,s);
        now_s+=1;
        pivot+=1;
      }else if(now_s==s){
        pivot+=1;
      }
    }

  }
  
  Rcpp::List out;
  out["C"] = temp.init_B*trans(temp.init_V);
  out["iters"] = temp.one_iter;
  out["Aset"] = temp.one_aset;
  return out;
}
    

// [[Rcpp::export]]
Rcpp::List srrmlr_once(arma::dmat &X, arma::dmat &Y,
                      arma::uword s, arma::uword r,
                      double e_in=0.001, double e_out=0.001){
  if(r>s){
    printf("r must not be bigger than s! \n");
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
  arma::umat A(s, 10000);
  arma::uword cnt=0;
  double med;

  for(iter=0; iter<1000; iter++){
    pre_C = B*trans(V);
    for(arma::uword _=0; _<p; _++){dB(_) = arma::norm(B.row(_), 2);}
    iset_.reset();
    arma::uvec iset_ = arma::find(dB==0);

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

      A.col(cnt) = aset;
      cnt+=1;
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
  out["iters"] = cnt;
  out["Prob"] = prob;
  out["Aset"] = A;
  return out;
}

