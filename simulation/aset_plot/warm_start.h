#include <stdio.h>
#include <armadillo>

typedef struct{
  arma::dmat init_B;
  arma::dmat init_V;
  arma::uword one_iter;
  arma::umat one_aset;
}SRR;

SRR srrmlr_mode1(arma::dmat &X, arma::dmat &Y, arma::dmat B, arma::dmat V,
                 arma::uword s, arma::uword r, arma::uword s_true,
                 double e_in=0.001, double e_out=0.001){
  if(r>s){
    printf("r must not be bigger than s! \n");
    std::exit(0);
  }
  const arma::uword n = X.n_rows, p = X.n_cols, q = Y.n_cols;
  arma::dmat pre_C(p, q-1), pr(n, q-1), M(n, q-1), Ga(p, r), obj(q-1, q-1);
  arma::uvec indices, ord, ordset, aset_, iset_;
  arma::dvec dC(p);
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::uword iter;
  arma::umat A(s_true, 10000, arma::fill::value(INFINITY));
  arma::uword cnt=0;
  double med;
  arma::dmat C = B*trans(V);
  
  for(iter=0; iter<1000; iter++){
    for(arma::uword _=0; _<p; _++){dC(_) = arma::norm(C.row(_), 2);}
    aset_ = arma::find(dC!=0);
    iset_ = arma::find(dC==0);

    for(arma::uword k=0; k<n; k++){
      med = accu(exp(trans(C) * trans(X.row(k))));
      for(arma::uword _=0; _<q-1; _++){
        pr(k,_) = accu(exp(trans(C.col(_))*trans(X.row(k))))/(1+med);
      }
    }
    M = 2*(Y.cols(1,q-1)-pr)+X*C;
    for(arma::uword in_loop=0; in_loop<500; in_loop++){
      pre_C = C;
      obj = trans(M)*X.cols(aset_)*inv(trans(X.cols(aset_))*X.cols(aset_))*trans(X.cols(aset_))*M;
      eig_gen(eigval, eigvec, obj);
      ord = arma::stable_sort_index(real(eigval), "descend");
      ordset = ord.subvec(0, r-1);
      V = real(eigvec.cols(ordset));
      B.zeros(); Ga.zeros();
      B.rows(aset_) = inv(trans(X.cols(aset_))*X.cols(aset_))*trans(X.cols(aset_))*M*V;
      Ga.rows(iset_) = trans(X.cols(iset_))*(M*V-X*B)/n;
      
      for(arma::uword _=0; _<p; _++){dC(_) = arma::norm(B.row(_)+Ga.row(_), 2);}
      indices = arma::sort_index(dC, "descend");
      aset_ = arma::sort(indices.subvec(0, s-1));
      iset_ = arma::sort(indices.subvec(s, p-1));
      C = B*trans(V);
      A.col(cnt).head(s) = aset_;
      cnt+=1;
      if(accu(square(B*trans(V)-C))<e_in) break;
    }
    if(accu(square(pre_C-C))<e_out){
      for(arma::uword _=0; _<p; _++){dC(_) = arma::norm(B.row(_)+Ga.row(_), 2);}
      indices = arma::sort_index(dC, "descend");
      arma::uvec aset1 = arma::sort(indices.subvec(0, s));
      arma::uvec iset1 = arma::sort(indices.subvec(s+1, p-1));
      arma::dmat obj1 = trans(M)*X.cols(aset1)*inv(trans(X.cols(aset1))*X.cols(aset1))*trans(X.cols(aset1))*M;
      eig_gen(eigval, eigvec, obj1);
      ord = arma::stable_sort_index(real(eigval), "descend");
      ordset = ord.subvec(0, r-1);
      V = real(eigvec.cols(ordset));
      B.zeros();
      B.rows(aset1) = inv(trans(X.cols(aset1))*X.cols(aset1))*trans(X.cols(aset1))*M*V;
      break;
    }
  }

  SRR temp;
  temp.init_B = B;
  temp.init_V = V;
  temp.one_aset = A;
  temp.one_iter = cnt;
  return temp;
}


SRR srrmlr_mode2(arma::dmat &X, arma::dmat &Y, arma::dmat B, arma::dmat V,
                 arma::uword s, arma::uword r, arma::uword s_true,
                 double e_in=0.001, double e_out=0.001){
  if(r>s){
    printf("r must not be bigger than s! \n");
    std::exit(0);
  }
  const arma::uword n = X.n_rows, p = X.n_cols, q = Y.n_cols;
  arma::dmat pre_C(p, q-1), pr(n, q-1), M(n, q-1), Ga(p, r), obj(q-1, q-1), B1(p,r+1), V1(q-1,r+1);
  arma::uvec indices, ord, ordset, aset_, iset_;
  arma::dvec dC(p);
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::uword iter;
  arma::umat A(s_true, 10000, arma::fill::value(INFINITY));
  arma::uword cnt=0;
  double med;
  
  arma::dmat C = B*trans(V);
  
  for(iter=0; iter<1000; iter++){
    for(arma::uword _=0; _<p; _++){dC(_) = arma::norm(C.row(_), 2);}
    aset_ = arma::find(dC!=0);
    iset_ = arma::find(dC==0);
    
    for(arma::uword k=0; k<n; k++){
      med = accu(exp(trans(C) * trans(X.row(k))));
      for(arma::uword _=0; _<q-1; _++){
        pr(k,_) = accu(exp(trans(C.col(_))*trans(X.row(k))))/(1+med);
      }
    }
    M = 2*(Y.cols(1,q-1)-pr)+X*C;
    
    for(arma::uword in_loop=0; in_loop<500; in_loop++){
      pre_C = C;
      obj = trans(M)*X.cols(aset_)*inv(trans(X.cols(aset_))*X.cols(aset_))*trans(X.cols(aset_))*M;
      eig_gen(eigval, eigvec, obj);
      ord = arma::stable_sort_index(real(eigval), "descend");
      ordset = ord.subvec(0, r-1);
      V = real(eigvec.cols(ordset));
      B.zeros(); Ga.zeros();
      B.rows(aset_) = inv(trans(X.cols(aset_))*X.cols(aset_))*trans(X.cols(aset_))*M*V;
      Ga.rows(iset_) = trans(X.cols(iset_))*(M*V-X*B)/n;
      
      for(arma::uword _=0; _<p; _++){dC(_) = arma::norm(B.row(_)+Ga.row(_), 2);}
      indices = arma::sort_index(dC, "descend");
      aset_ = arma::sort(indices.subvec(0, s-1));
      iset_ = arma::sort(indices.subvec(s, p-1));
      C = B*trans(V);
      A.col(cnt).head(s) = aset_;
      cnt+=1;
      if(accu(square(B*trans(V)-C))<e_in) break;
    }
    if(accu(square(pre_C-C))<e_out){
      obj = trans(M)*X.cols(aset_)*inv(trans(X.cols(aset_))*X.cols(aset_))*trans(X.cols(aset_))*M;
      eig_gen(eigval, eigvec, obj);
      ord = arma::stable_sort_index(real(eigval), "descend");
      arma::uvec ordset1 = ord.subvec(0, r);
      V1 = real(eigvec.cols(ordset1));
      B1.zeros();
      B1.rows(aset_) = inv(trans(X.cols(aset_))*X.cols(aset_))*trans(X.cols(aset_))*M*V1;
      break;
    }
  }
  
  SRR temp;
  temp.init_B = B1;
  temp.init_V = V1;
  temp.one_aset = A;
  temp.one_iter = cnt;
  return temp;
}


SRR srrmlr_one(arma::dmat &X, arma::dmat &Y,
               arma::dmat B, arma::dmat V,
               arma::uword s, arma::uword r,
               double e_in=0.001, double e_out=0.001){
  if(r>s){
    printf("r must not be bigger than s! \n");
    std::exit(0);
  }
  const arma::uword n = X.n_rows, p = X.n_cols, q = Y.n_cols;
  arma::dmat pre_C(p, q-1), pr(n, q-1), M(n, q-1), Ga(p, r), obj(q-1, q-1), B1(p,r+1), V1(q-1,r+1);
  arma::uvec indices, ord, ordset, aset_, iset_;
  arma::dvec dC(p);
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::uword iter;
  arma::umat A(s, 10000, arma::fill::value(INFINITY));
  arma::uword cnt=0;
  double med;
  arma::dmat C = B*trans(V);
  
  for(iter=0; iter<1000; iter++){
    for(arma::uword _=0; _<p; _++){dC(_) = arma::norm(C.row(_), 2);}
    aset_ = arma::find(dC!=0);
    iset_ = arma::find(dC==0);
    
    for(arma::uword k=0; k<n; k++){
      med = accu(exp(trans(C) * trans(X.row(k))));
      for(arma::uword _=0; _<q-1; _++){
        pr(k,_) = accu(exp(trans(C.col(_))*trans(X.row(k))))/(1+med);
      }
    }
    M = 2*(Y.cols(1,q-1)-pr)+X*C;
    
    for(arma::uword in_loop=0; in_loop<500; in_loop++){
      pre_C = C;
      obj = trans(M)*X.cols(aset_)*inv(trans(X.cols(aset_))*X.cols(aset_))*trans(X.cols(aset_))*M;
      eig_gen(eigval, eigvec, obj);
      ord = arma::stable_sort_index(real(eigval), "descend");
      ordset = ord.subvec(0, r-1);
      V = real(eigvec.cols(ordset));
      B.zeros(); Ga.zeros();
      B.rows(aset_) = inv(trans(X.cols(aset_))*X.cols(aset_))*trans(X.cols(aset_))*M*V;
      Ga.rows(iset_) = trans(X.cols(iset_))*(M*V-X*B)/n;
      
      for(arma::uword _=0; _<p; _++){dC(_) = arma::norm(B.row(_)+Ga.row(_), 2);}
      indices = arma::sort_index(dC, "descend");
      aset_ = arma::sort(indices.subvec(0, s-1));
      iset_ = arma::sort(indices.subvec(s, p-1));
      C = B*trans(V);
      A.col(cnt) = aset_;
      cnt+=1;
      if(accu(square(B*trans(V)-C))<e_in) break;
    }
    if(accu(square(pre_C-C))<e_out) break;
  }

  SRR temp;
  temp.init_B = B;
  temp.init_V = V;
  temp.one_aset = A;
  temp.one_iter = cnt;
  return temp;
}


SRR get_init(arma::dmat &X, arma::dmat &Y, arma::uword s_true){
  const arma::uword n = X.n_rows, p = X.n_cols, q = Y.n_cols;
  arma::dmat pr(n, q-1), M(n, q-1), C(p, q-1, arma::fill::zeros), obj(q-1, q-1);
  arma::uvec indices, ord;
  arma::dvec dB(p);
  arma::dvec eigval;
  arma::cx_vec eigval1;
  arma::cx_mat eigvec1;
  arma::dmat U1,V1;
  arma::dvec Sigma;
  arma::umat A(s_true, 1, arma::fill::value(INFINITY));
  arma::uword cnt, ordset;
  double med;

  for(arma::uword k=0; k<n; k++){
    med = accu(exp(trans(C) * trans(X.row(k))));
    for(arma::uword _=0; _<q-1; _++){
      pr(k,_) = accu(exp(trans(C.col(_))*trans(X.row(k))))/(1+med);
    }
  }
  M = 2*(Y.cols(1,q-1)-pr)+X*C;
  svd(U1, Sigma, V1, M);
  arma::dmat V = V1.col(0);
  dB = trans(X)*M*V;
  indices = arma::stable_sort_index(dB,"descend");
  arma::uword aset = indices(0);
  A(0,0) = aset;
  arma::dmat B(p, 1, arma::fill::zeros);
  B.row(aset) = inv(trans(X.col(aset))*X.col(aset)) * trans(X.col(aset)) *M*V;
  
  obj = trans(M)*X.col(aset)*inv(trans(X.col(aset))*X.col(aset))*trans(X.col(aset))*M;
  eig_gen(eigval1, eigvec1, obj);
  ord = arma::stable_sort_index(real(eigval1), "descend");
  ordset = ord(0);
  V = real(eigvec1.col(ordset));
  
  SRR temp;
  temp.init_B = B;
  temp.init_V = V;
  temp.one_aset = A;
  temp.one_iter = cnt;
  return temp;
}