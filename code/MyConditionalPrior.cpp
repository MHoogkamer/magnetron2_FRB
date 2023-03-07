#include "MyConditionalPrior.h"
#include "DNest4/code/Utils.h"
#include "Data.h"
#include <cmath>
#include <gsl/gsl_cdf.h>
#include <iostream>
using namespace DNest4;

MyConditionalPrior::MyConditionalPrior(double x_min, double x_max, 
					double mu_min, double mu_max)
:x_min(x_min) // start observation 
,x_max(x_max) // end observation 
,mu_min(mu_min)
,mu_max(mu_max)
,min_width(0.1*0.983)
,max_width(10.0)
{

}

void MyConditionalPrior::from_prior(RNG& rng)
{
//	mu = tan(M_PI*(0.97*rng.rand() - 0.485)); // mean amplitude 
//	mu = exp(mu);
//	mu_widths = exp(log(0.33*0.983) + (log(5*0.983) - log(0.33*0.983))*rng.rand()); // mean risetime hardcoded with CHIME time resolution 983 us (CHIME data is in ms)

//	sig = 2.*rng.rand(); // width amplitude
//	sig_widths = 2.*rng.rand(); // width risetime

        const DNest4::Cauchy cauchy_mu(1.0, 2.0);
        do
        {
            mu = cauchy_mu.generate(rng);
        }while(mu >= 5.5 || mu <= - 2.3);
        mu = exp(mu);

        const DNest4::Cauchy cauchy_width(0.0, 1.0);
        do
        {
             mu_widths = cauchy_width.generate(rng);
        }while(mu_widths >= 2.3 || mu_widths <= -1.1);
        mu_widths = exp(mu_widths);


        sig = exp(-2 + 2.5*rng.rand());
        sig_widths = exp(-2 + 2.*rng.rand());

	a = -5. + 10.*rng.rand(); // mean skewness
	b = 2.5*rng.rand(); // width skewness 
}

double MyConditionalPrior::perturb_hyperparameters(RNG& rng)
{
	double logH = 0.;

	int which = rng.rand_int(6);

	if(which == 0)
	{
//		mu = log(mu);
//		mu = (atan(mu)/M_PI + 0.485)/0.97;
//		mu += pow(10., 1.5 - 6.*rng.rand())*rng.randn();
//		mu = mod(mu, 1.);
//		mu = tan(M_PI*(0.97*mu - 0.485));
//		mu = exp(mu);
               const DNest4::Cauchy cauchy_mu(1.0, 2.0);

               mu = log(mu);
               logH += cauchy_mu.perturb(mu, rng);
               if(mu >= 5.5 || mu <= - 2.3)
               {
                   mu=2.0;
                   return -1E300;
               }
               mu = exp(mu);
	}
	if(which == 1)
	{  
//                mu_widths = log(mu_widths);
//                mu_widths += (log(5*0.983) - log(0.33*0.983))*rng.randh();
//                mu_widths = mod(mu_widths + log(0.33*0.983), log(5*0.983) - log(0.33*0.983)) - log(0.33*0.983);
//                mu_widths = exp(mu_widths); 
               const DNest4::Cauchy cauchy_width(0.0, 1.0);

               mu_widths = log(mu_widths);
               logH += cauchy_width.perturb(mu_widths, rng);
               if(mu_widths >= 2.3 || mu_widths <= -1.1)
               {
                   mu_widths=1.0;
                   return -1E300;
               }
               mu_widths = exp(mu_widths);
	}
	if(which == 2)
	{
//		sig += 2.*rng.randh();
//		sig = mod(sig, 2 .);
                sig = log(sig);
                sig += 2.5*rng.randh();
                sig = mod(sig + 2., 2.5) - 2.;
                sig = exp(sig);
	}
	
	if(which == 3)
	{ 
//		sig_widths += 2.*rng.randh();
//		sig_widths = mod(sig_widths, 2.);
                sig_widths = log(sig_widths);
                sig_widths += 2.*rng.randh();
                sig_widths = mod(sig_widths + 2, 2.0) - 2.0;
                sig_widths = exp(sig_widths);
	}
        if(which == 4)
	{
		a += 10.*rng.randh();
		a = mod(a + 5., 10.) - 5.;
	}
	if(which == 5)
	{
		b += 2.5*rng.randh();
		b = mod(b, 2.5);
	}

	return logH;
}

double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
        //std::cout << min_width << " " << max_width << " " << mu_min << " " << mu_max << std::endl;
        //std::cout << "amp: " << vec[1] << ", rise: " << vec[2] << std::endl;
	if(vec[0] < x_min || vec[0] > x_max || vec[1] < 0.0 || vec[2] < min_width
                || vec[1] < mu_min || vec[1] > mu_max || vec[2] > max_width
		|| log(vec[3]) < (a-b) || log(vec[3]) > (a + b))
                //std::cout << "Rejecting component" << std::endl;
		return -1E300;

        //std::cout << "Not rejecting component" << std::endl;
	return	- log(vec[1]*sig) - 0.5*pow((log(vec[1]) - log(mu))/sig, 2)
		- log(vec[2]*sig_widths) - 0.5*pow((log(vec[2]) - log(mu_widths))/sig_widths, 2)
		- log(2.*b*vec[3]);
}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const // inverse CDF 
{
	vec[0] = x_min + (x_max - x_min)*vec[0]; // peak time
	vec[1] = exp(log(mu) + sig*gsl_cdf_ugaussian_Pinv(vec[1])); //-mu * log(1. - vec[1]); // left lognormal, rigth exp amplitude 
	vec[2] = exp(log(mu_widths) + sig_widths*gsl_cdf_ugaussian_Pinv(vec[2])); // risetime
	vec[3] = exp(a - b + 2.*b*vec[3]); // skewness
}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const // CDF
{
	vec[0] = (vec[0] - x_min)/(x_max - x_min);
	vec[1] = gsl_cdf_ugaussian_P((log(vec[1]) - log(mu))/sig); // 1. -exp(-vec[1] / mu); // left lognormal, rigth exp amplitude
	vec[2] = gsl_cdf_ugaussian_P((log(vec[2]) - log(mu_widths))/sig_widths);
	vec[3] = (log(vec[3]) + b - a)/(2.*b);
}

void MyConditionalPrior::print(std::ostream& out) const
{
//	out<<mu<<' '<<mu_widths<<' '<<a<<' '<<b<<' ';
	out<<mu<<' '<<sig<<' '<<mu_widths<<' '<<sig_widths<<' '<<a<<' '<<b<<' ';

}

