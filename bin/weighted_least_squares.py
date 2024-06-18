import pandas as pd

#create DataFrame
df = pd.DataFrame({'hours': [1, 1, 2, 2, 2, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 8],
                   'score': [48, 78, 72, 70, 66, 92, 93, 75, 75, 80, 95, 97,
                             90, 96, 99, 99]})

#view first five rows of DataFrame
print(df.head())

import statsmodels.api as sm

#define predictor and response variables
y = df['score']
x = df['hours']

#add constant to predictor variables
x = sm.add_constant(x)

#fit linear regression model
fit = sm.OLS(y, x).fit()

#view model summary
print(fit.summary())

#define weights to use
wt = 1 / smf.ols('fit.resid.abs() ~ fit.fittedvalues', data=df).fit().fittedvalues**2

#fit weighted least squares regression model
fit_wls = sm.WLS(y, X, weights=wt).fit()

#view summary of weighted least squares regression model
print(fit_wls.summary())

import numpy as np
import matplotlib.pyplot as plt
res1 = np.array([2,3,-4,1,6,-2,7])
# res2 = np.array([2,2,-4,1,-3,-9,7])
x = np.arange(7)
sin = np.sin(x)
cos = np.cos(x)
syn1 = np.arange(7) + 3*sin + 2*cos + 2
syn2 = np.arange(7) + 1*sin + 3*cos - 1
y1 = syn1 + res1
y2 = syn2 + res1
# y1[2] = np.nan
# y2[3] = np.nan
# y2[4] = np.nan

X = np.vstack([x, sin, cos, np.ones_like(x)]).transpose()
Xo = (np.vstack([x, sin, cos, np.ones_like(x)])/res1).transpose()

Y = np.vstack([y1, y2]).transpose()
Yo = (np.vstack([y1, y2])/res1).transpose()
# weights = np.vstack([1/res1**2, 1/res2**2])
# weights = np.vstack([1/res1**2, 1/res2**2])

olsfit = sm.OLS(Yo, Xo).fit()
wlsfit = sm.WLS(Y, X, weights=1/res1**2).fit()
# model = np.zeros_like(Y)*np.nan
# model[~np.isnan(Y)] = wlsfit.fittedvalues
plt.scatter(x, y1, s=100/res1**2)
plt.scatter(x, y2, s=100/res1**2)
plt.plot(x, syn1, label='syn1')
plt.plot(x, syn2, label='syn2')
# plt.plot(x, olsfit.fittedvalues, label='ols')
plt.plot(x, olsfit.fittedvalues[:,0]*res1, label='ols1')
plt.plot(x, olsfit.fittedvalues[:,1]*res1, label='ols2')

plt.plot(x, wlsfit.fittedvalues, label='wls')

plt.legend()
plt.show()

import statsmodels.api as sm
sin=np.sin(2*np.pi*dt_cum)
cos=np.cos(2*np.pi*dt_cum)
# X = np.vstack([dt_cum, sin, cos, np.ones_like(dt_cum)]).transpose()
X = np.vstack([dt_cum, np.ones_like(dt_cum)]).transpose()

Y = np.vstack([range_coefs, azi_coefs]).transpose()
flat_std = np.array([np.nanstd(flat_cum[i, :, :]) for i in np.arange(n_im)])
flat_std[0] = np.std(flat_std)
weights = 1 / flat_std**2
plt.scatter(epochs, range_coefs, s=200 * weights, label="range_coef")
plt.scatter(epochs, azi_coefs, s=200 * weights, label="azi_coef")
plt.plot(epochs, range_coefs, label="range_coef")
plt.plot(epochs, azi_coefs, label="azi_coef")
wlsfit = sm.WLS(Y, X, weights=weights).fit()

range_slope = []
azi_slope = []
for t in np.arange(4):
    window = np.logical_and(dt_cum>t, dt_cum<t+5)
    wlsfit = sm.WLS(Y[window], X[window], weights=weights[window]).fit()
    range_slope.append(wlsfit.params[0,0])
    plt.plot(np.array(epochs)[window], wlsfit.fittedvalues, label='wls')
    azi_slope.append(wlsfit.params[1,0])

    # olsfit = sm.OLS(Y[window], X[window]).fit()
    # range_slope.append(olsfit.params[0, 0])
    # azi_slope.append(olsfit.params[1, 0])
    # plt.plot(np.array(epochs)[window], olsfit.fittedvalues, label='ols')

range_slope_std = np.std(np.array(range_slope))
azi_slope_std = np.std(np.array(azi_slope))
print(range_slope_std, azi_slope_std)
plt.xlabel("Epoch")
plt.ylabel("ramp rate unit/pixel")
plt.title(args.cumfile + ', range_slope_std={:.2f}, azi_slope_std={:.2f}'.format(range_slope_std, azi_slope_std))
plt.legend()
plt.savefig(args.cumfile + "_ramp_coefs_fit.png")
plt.close()

plt.plot(epochs, wlsfit.resid, label='wls.resid')
plt.xlabel("Epoch")
plt.ylabel("ramp rate unit/pixel")
plt.title(args.cumfile)
plt.legend()
plt.savefig(args.cumfile + "_ramp_coefs_resid.png")


Y = np.vstack([range_coefs, azi_coefs]).transpose()
Y = np.vstack([range_coefs]).transpose()

sig = flat_std[:, np.newaxis]
G = X
d = Y
model = np.linalg.lstsq(G / sig, d / sig, rcond=None)[0]
cov_d = np.diag(np.square(sig).transpose()[0])
cov_m = np.linalg.inv(np.dot(np.dot(G.transpose(), np.linalg.inv(cov_d)), G))
wlsfit.resid = resid = d - np.dot(G, model)
chi_square = np.sum(np.square(resid / sig))
reduced_chi_square = chi_square / (len(d) - len(model))

reduced_chi_square = 0.0028
np.sqrt(reduced_chi_square) = 0.05
wlsfit.mse_model = 0.2140
wlsfit.mse_resid = 0.0028 = reduced_chi_square
wlsfit.mse_total = 0.003957672855930073
chi_square = 0.549
wlsfit.bse = np.array([0.01083044, 0.05236671])
wlsfit.HC0_se = array([0.01602831, 0.06917684])

cov_m = array([[ 0.0409631 , -0.17970594],
               [-0.17970594,  0.95765946]])

std_m1 =  0.20239343007258423
std_m2 =  0.9786007680860049

wlsfit.bse = np.sqrt(reduced_chi_square * cov_m)


from scipy.fft import fft, fftfreq
import numpy as np
# Number of sample points
N = 600
# sample spacing
T = 1.0 / 800.0
x = np.linspace(0.0, N*T, N, endpoint=False)
y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
yf = fft(y)
xf = fftfreq(N, T)[:N//2]
import matplotlib.pyplot as plt
plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
plt.grid()
plt.show()

print("Pylenin")
print("loves")
print("Python")
print("Pylenin", end="\r")
1+3
print("loves", end="\r")
print("Python")