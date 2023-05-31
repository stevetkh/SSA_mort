#include <TMB.hpp>     

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  DATA_VECTOR(dm);
  DATA_VECTOR(Em);
  DATA_VECTOR(df);
  DATA_VECTOR(Ef);

  DATA_VECTOR(tp_common_mean);

  DATA_SPARSE_MATRIX(Dm);
  DATA_SPARSE_MATRIX(Df);
  DATA_SPARSE_MATRIX(Dm_age);
  DATA_SPARSE_MATRIX(Df_age);
  DATA_SPARSE_MATRIX(Dm_time);
  DATA_SPARSE_MATRIX(Df_time);

  DATA_SPARSE_MATRIX(tp_mat_m);
  DATA_SPARSE_MATRIX(tp_mat_f);
  DATA_SPARSE_MATRIX(tp_mat_m_common);
  DATA_SPARSE_MATRIX(tp_mat_f_common);

  DATA_SPARSE_MATRIX(intercept_mat_m);
  DATA_SPARSE_MATRIX(intercept_mat_f);

  DATA_SPARSE_MATRIX(penal_age);
  DATA_SPARSE_MATRIX(penal_age_cross);
  DATA_SPARSE_MATRIX(null_penal_age);

  DATA_SPARSE_MATRIX(penal_time1_m);
  DATA_SPARSE_MATRIX(penal_time2_m);
  DATA_SPARSE_MATRIX(penal_time1_f);
  DATA_SPARSE_MATRIX(penal_time2_f);
  DATA_SPARSE_MATRIX(penal_time_cross);
  DATA_SPARSE_MATRIX(null_penal_time_all);

  DATA_SPARSE_MATRIX(penal_age_2d);
  DATA_SPARSE_MATRIX(penal_time_2d);
  DATA_SPARSE_MATRIX(null_penal);
 
  DATA_SPARSE_MATRIX(penal_tp_common);
  DATA_SPARSE_MATRIX(null_penal_tp_common);
  DATA_SPARSE_MATRIX(penal_tp_0);
  DATA_SPARSE_MATRIX(penal_tp);
  DATA_SPARSE_MATRIX(null_penal_tp);

  DATA_SPARSE_MATRIX(rwanda_mat_m);
  DATA_SPARSE_MATRIX(rwanda_mat_f);

  DATA_SCALAR(sp_mean);
  DATA_SCALAR(sp_sd);

  PARAMETER(log_lambda_age_m);
  PARAMETER(log_lambda_age_f);
  PARAMETER(log_lambda_age_m_cross);
  PARAMETER(log_lambda_age_f_cross);

  PARAMETER(log_lambda_time_m1);
  PARAMETER(log_lambda_time_f1);
  PARAMETER(log_lambda_time_m2);
  PARAMETER(log_lambda_time_f2);
  PARAMETER(log_lambda_time_cross);

  PARAMETER(log_lambda_age_m_2d);
  PARAMETER(log_lambda_age_f_2d);
  PARAMETER(log_lambda_time_m_2d);
  PARAMETER(log_lambda_time_f_2d);

  PARAMETER(log_lambda_avg_m);
  PARAMETER(log_lambda_avg_f);

  PARAMETER(log_lambda_tp);
  PARAMETER(log_lambda_tp_common);
  PARAMETER(log_lambda_tp_common_0_inflated_sd);

  PARAMETER(log_dispersion_m);
  PARAMETER(log_dispersion_f);

  PARAMETER_VECTOR(spline_params_m_age);
  PARAMETER_VECTOR(spline_params_f_age);
  PARAMETER_VECTOR(spline_params_m_time);
  PARAMETER_VECTOR(spline_params_f_time);
  PARAMETER_VECTOR(spline_params_m_2d);
  PARAMETER_VECTOR(spline_params_f_2d);

  PARAMETER(tp_common_slope);
  PARAMETER(tp_common_5);
  PARAMETER(tp_common_10);

  PARAMETER(avg_common_m);
  PARAMETER(avg_common_f);
  PARAMETER_VECTOR(avg_m);
  PARAMETER_VECTOR(avg_f);
  
  PARAMETER_VECTOR(tips_params);
  PARAMETER_VECTOR(tips_params_common);
  
  PARAMETER_VECTOR(rwanda_geno_m);
  PARAMETER_VECTOR(rwanda_geno_f);

  SparseMatrix<Type> QQm, QQf, QQm_age, QQf_age, QQ_time, QQ_tp_common, QQ_tp;
  vector<Type> time_params_all(spline_params_m_time.size() + spline_params_f_time.size());
  time_params_all << spline_params_m_time, spline_params_f_time;

  Type nll;
  QQm = exp(log_lambda_age_m_2d) * penal_age_2d + exp(log_lambda_time_m_2d) * penal_time_2d + null_penal;
  QQf = exp(log_lambda_age_f_2d) * penal_age_2d + exp(log_lambda_time_f_2d) * penal_time_2d + null_penal;

  QQm_age = exp(log_lambda_age_m) * penal_age + exp(log_lambda_age_m_cross) * penal_age_cross + null_penal_age;
  QQf_age = exp(log_lambda_age_f) * penal_age + exp(log_lambda_age_f_cross) * penal_age_cross + null_penal_age;
  QQ_time = exp(log_lambda_time_m1)*penal_time1_m + exp(log_lambda_time_m2) * penal_time2_m + 
	    exp(log_lambda_time_f1) * penal_time1_f + exp(log_lambda_time_f2) * penal_time2_f + 
	    exp(log_lambda_time_cross) * penal_time_cross + null_penal_time_all;

  QQ_tp_common = exp(log_lambda_tp_common) * penal_tp_common + exp(-2 * log_lambda_tp_common_0_inflated_sd) * penal_tp_0 + null_penal_tp_common;
  QQ_tp = exp(log_lambda_tp) * penal_tp + null_penal_tp;

  nll += GMRF(QQm)(spline_params_m_2d) + GMRF(QQf)(spline_params_f_2d);
  nll += GMRF(QQm_age)(spline_params_m_age) + GMRF(QQf_age)(spline_params_f_age);
  nll += GMRF(QQ_time)(time_params_all);
  nll += GMRF(QQ_tp_common)(tips_params_common-tp_common_slope * tp_common_mean);
  nll += GMRF(QQ_tp)(tips_params);

  vector<Type> mum, muf, varm, varf, tips_params_common_after;
  tips_params_common_after = tips_params_common;
  tips_params_common_after(5) += tp_common_5;
  tips_params_common_after(10) += tp_common_10;

  mum = exp(avg_common_m + intercept_mat_m * avg_m + rwanda_mat_m * rwanda_geno_m +
	    Dm_age*spline_params_m_age + Dm_time * spline_params_m_time + Dm * spline_params_m_2d + 
	    tp_mat_m * tips_params + tp_mat_m_common * tips_params_common_after) * Em;
  muf = exp(avg_common_f + intercept_mat_f * avg_f + rwanda_mat_f * rwanda_geno_f + 
	    Df_age * spline_params_f_age + Df_time * spline_params_f_time + Df * spline_params_f_2d + 
	    tp_mat_f * tips_params + tp_mat_f_common * tips_params_common_after) * Ef;
  varm = mum * (1 + mum/exp(log_dispersion_m));
  varf = muf * (1 + muf/exp(log_dispersion_f));

  nll -= dnbinom2(dm, mum, varm, 1).sum() + dnbinom2(df, muf, varf, 1).sum();

  nll -= dnorm(log_lambda_age_m, sp_mean, sp_sd, 1) + dnorm(log_lambda_age_f, sp_mean, sp_sd, 1);
  nll -= dnorm(log_lambda_age_m_cross, Type(3.0), Type(5.0), 1) + dnorm(log_lambda_age_f_cross, Type(3.0), Type(5.0), 1);
  
  nll -= dnorm(log_lambda_time_m1, sp_mean, sp_sd, 1) + dnorm(log_lambda_time_f1, sp_mean, sp_sd, 1);
  nll -= dnorm(log_lambda_time_m2, sp_mean, sp_sd, 1) + dnorm(log_lambda_time_f2, sp_mean, sp_sd, 1);
  nll -= dnorm(log_lambda_time_cross, Type(3.0), Type(5.0), 1);

  nll -= dnorm(log_lambda_age_m_2d, sp_mean, sp_sd, 1) + dnorm(log_lambda_age_f_2d, sp_mean, sp_sd, 1);
  nll -= dnorm(log_lambda_time_m_2d, sp_mean, sp_sd, 1) + dnorm(log_lambda_time_f_2d, sp_mean, sp_sd, 1);

  nll -= dnorm(log_dispersion_m, Type(3.0), Type(10.0), 1);
  nll -= dnorm(log_dispersion_f, Type(3.0), Type(10.0), 1);

  nll -= dnorm(log_lambda_tp, Type(0.0), Type(5.0), 1);
  nll -= dnorm(log_lambda_tp_common, Type(0.0), Type(5.0), 1);
  nll -= dnorm(log_lambda_tp_common_0_inflated_sd, Type(0.0), Type(5.0), 1);
  nll -= dnorm(tp_common_slope, Type(0.0), Type(10.0), 1) + dnorm(tp_common_5, Type(0.0), Type(10.0), 1) + dnorm(tp_common_10, Type(0.0), Type(10.0), 1);

  nll -= dnorm(avg_common_m, Type(0.0), Type(10.0), 1) + dnorm(avg_common_f, Type(0.0), Type(10.0), 1);
  nll -= dnorm(avg_m, Type(0.0), exp(log_lambda_avg_m), 1).sum() + dnorm(avg_f, Type(0.0), exp(log_lambda_avg_f), 1).sum();

  nll -= dnorm(rwanda_geno_m, Type(0.0), Type(10.0), 1).sum() + dnorm(rwanda_geno_f, Type(0.0), Type(10.0), 1).sum();

  return nll;
}
