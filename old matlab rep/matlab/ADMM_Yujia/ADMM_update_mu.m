function mu_out = ADMM_update_mu(mu_in,mu_tol,mu_inc,mu_dec,primal_res,dual_res)
    if (primal_res/dual_res) > mu_tol
        mu_out = mu_in*mu_inc;
    elseif (primal_res/dual_res) < 1/mu_tol
        mu_out = mu_in/mu_dec;
    else
        mu_out = mu_in;
    end
end