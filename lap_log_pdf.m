function y  = lap_log_pdf(x, mu, sigma)
    y = -abs(x-mu)./sigma - log(2.*sigma);