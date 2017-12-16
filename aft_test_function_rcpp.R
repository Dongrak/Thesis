order=base::order
cppFunction('
            List W_omni(arma::vec b, arma::vec Time, arma::vec Delta, arma::vec Covari){
            int n = Time.size();
            int p = b.size();
            arma::vec e_i_beta = log(Time)+Covari*b;
            arma::vec order_resid = as<NumericVector>(wrap(arma::sort_index( e_i_beta )+1)) ;
            
            
            return List::create(order_resid);
            }',depends="RcppArmadillo")

W_omni(beta_hat_ln_aft,X_ln_cox,D_ln_cox,Z_ln_cox)


e_i_beta=as.vector(log(X_ln_cox)+Z_ln_cox%*%beta_hat_ln_aft)
order(e_i_beta)