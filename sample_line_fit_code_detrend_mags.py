#sample python code for fitting a line to lightcurve data and then subtracting it from the data

x_data = dates_mjd_utc #xdata, y_data and y_data_unc are generic variable names to help with input to scipy.optimize.curve_fit function
y_data = mags
y_data_unc = mags_unc

function_lambda = lambda x, a, b: a*x + b
guess_params_numpy_array = np.array([1,1])
result = scipy.optimize.curve_fit(function_lambda, x_data, y_data, guess_params_numpy_array, y_data_unc, full_output=True)
pars, corr = result[0], result[1]

mags_detrended = (mags - function_lambda(x_data,*pars)) + np.median(mags)

