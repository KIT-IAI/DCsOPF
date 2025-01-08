#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame as DF
import matplotlib.image as img
from ipywidgets import interact, interactive, fixed, interact_manual, Dropdown, IntSlider, widgets
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
import math
import scipy
import scipy.io

# gpflow
import gpflow
import tensorflow as tf
from gpflow.utilities import print_summary

# gpytorch
import torch
import gpytorch


PATH = './data'
PATH_SAVE = './graphics/GPR'


# ========== HELPER FUNCTIONS ============================


# preprocessing

def drop_nans(data):
    if isinstance(data, dict):
        return {el : drop_nans(data[el]) for el in data}
    elif isinstance(data, pd.DataFrame):
        return data.dropna(axis=1, how='all')
    
def get_numeric_columns(df):
    numeric_columns = list(df.select_dtypes(include=[np.number]).columns.values)
    return numeric_columns

def get_nonnumeric_columns(df):
    nonnumeric_columns = list(df.select_dtypes(exclude=[np.number]).columns.values)
    return nonnumeric_columns
    
def set_index_to_nonnumeric_cols_df(df):
    nonnumeric_columns = get_nonnumeric_columns(df)
    numeric_columns = get_numeric_columns(df)
    if nonnumeric_columns != []:
        if numeric_columns != []:
            df = df.set_index(nonnumeric_columns)
            df.index.name = '_'.join(nonnumeric_columns)
        else:
            df = df.set_index(nonnumeric_columns[0])
    else:
        df.index.name = 'hours'
    if is_multi_index(df.index):
        df = remove_multiindex_df(df, axis=0)
    return df

def set_index_to_nonnumeric_cols(data):
    if isinstance(data, dict):
        return {el : set_index_to_nonnumeric_cols(data[el]) for el in data}
    elif isinstance(data, pd.DataFrame):
        return set_index_to_nonnumeric_cols_df(data)

def is_multi_index(lst):
    return isinstance(lst[0], tuple)

def remove_multiindex_df(df, axis=0):
    if axis==0:
        if is_multi_index(df.index):
            df.index = ['_'.join(el) for el in list(df.index)]
    elif axis==1:
        if is_multi_index(df.columns):
            df.columns = ['_'.join(el) for el in list(df.columns)]
    return df

def transpose(data):
    if isinstance(data, dict):
        return {el : transpose(data[el]) for el in data}
    elif isinstance(data, pd.DataFrame):
        return remove_multiindex_df(data.transpose(), axis=1)

#-------------------------------------------------
# visualisation
    
def simple_plot(df, x_col, y_col):
    fig = plt.figure(figsize=(20,5))
    plt.plot(df[y_col])
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.xticks(rotation=90)
    plt.title(sheet)
    plt.show()

def visualise_data(sheet, y_col, x_ticks=False, sort=False):
    
    df = frame[sheet]
    
    if sort:
        df = df.sort_values(by=[y_col], ascending=False)#.reset_index(drop=True)
    
    # no numeric values in data frame
    if not (get_numeric_columns(df) == []):
        # plot
        fig = plt.figure(figsize=(20,5))
        for y in [y_col]:
            plt.scatter(df.index, df[y])
        plt.xlabel(df.index.name, fontsize=15)
        plt.ylabel(y_col, fontsize=15)
        if x_ticks:
            plt.xticks(rotation=90)
        else:
            plt.xticks([])
        plt.title(sheet, fontsize=20)
        plt.show()
    else:
        print('Nothing to plot here!')

#-------------------------------------------------
# sample data
        
def scale_data(data):
    # add x-axis of hours
    data['hours'] = data.index.astype('float64')
    data.columns = ['demand', 'hours']

    # scale to [0,1]
    scaler = MinMaxScaler()
    data['hours'] = DF(scaler.fit_transform(DF(data.hours)))
    data['demand'] = DF(scaler.fit_transform(DF(data['demand'])))
    return data

def sample_data_df(df, n=10):
    df_sampled = df.sample(n, axis=0).sort_values(by='hours')
    Y = np.array(df_sampled['demand']).reshape(-1,1)
    X = np.array(df_sampled['hours']).reshape(-1,1)
    return X, Y

def sample_data_dict(dic, n=10):
    samples = { county : sample_data_df(df, n) for county, df in dic.items() }
    return samples

def plot_sample(X, Y, SAVE_FIGS=False):
    plt.plot(X,Y, 'kx', mew=2)
    plt.xlabel('hours (scaled)')
    plt.ylabel('MWh (scaled)')
    plt.title('Sample data')
    if SAVE_FIGS:
        plt.savefig(str(PATH_SAVE + '/test_data_' + str(n) + '.png'))
    plt.show()

#------------------------------------------------------

# plot image of Turkish transmission system
def plot_Grid(PATH):
    fig = plt.figure(figsize=(10,5))
    plt.imshow(img.imread(PATH + '/Grid.png'), aspect='auto')
    plt.show()
    
def load_data(PATH, sheet=None):
    excel_file = pd.read_excel(PATH + '/Counties Hourly Demand.xlsx', header=[0], sheet_name=sheet)

    # delete NAN columns & turn multi-index into single-index
    data = set_index_to_nonnumeric_cols(drop_nans(excel_file))
    
    if sheet:
        data = data[sheet[0]]
        
    return data

def create_output_file(L, mean):
    out = dict()
    out['__header__'] = b'MATLAB 5.0 MAT-file, Platform: PCWIN64'
    out['__version__'] = '1.0'
    out['__globals__'] = []
    out['Lpost'] = L
    out['mu_post'] = mean
    return out

def calc_cholesky(cov, lower=True):
    try:
        L = scipy.linalg.cholesky(cov, lower=True)
        return L
    except:
        print(county + ": Covariance matrix is singular!")


# ========== GPFLOW ============================


def GPFLOW_set_params(model, kernel, parameters):
    likelihood_variance = parameters[0]
    kernel_lengthscales = parameters[1]
    kernel_variance = parameters[2]
    
    model.likelihood.variance.assign(likelihood_variance)
    if kernel.name not in ['sum', 'product']: # single kernel
        model.kernel.variance.assign(kernel_variance)
        if kernel.name not in ['linear', 'polynomial', 'constant', 'product']:
            model.kernel.lengthscales.assign(kernel_lengthscales)
    else: # multiple kernels
        for i in range(len(model.kernel.kernels)):
            model.kernel.kernels[i].variance.assign(kernel_lengthscales)
            if model.kernel.kernels[i].name not in ['linear', 'polynomial', 'constant', 'product']:
                model.kernel.kernels[i].lengthscales.assign(kernel_variance)
                
    return model

def GPFLOW_model(X, Y, gpflow_kernel, gpflow_meanf, params=[0.01, 0.1, 0.1]):
    m = gpflow.models.GPR(data=(X, Y), kernel=gpflow_kernel, mean_function=meanf)
    m = GPFLOW_set_params(m, gpflow_kernel, params)
    return m

def GPFLOW_optimize(m):
    opt = gpflow.optimizers.Scipy()
    opt_logs = opt.minimize(m.training_loss, m.trainable_variables, options=dict(maxiter=100))
    return opt_logs
        
def GPFLOW_predict(m, T):
    xx = np.linspace(-0.1, 1.1, T).reshape(T, 1)
    mean, var = m.predict_f(xx)
    _, cov_full = m.predict_f(xx, full_cov=True)
    cov = cov_full.numpy()[0,:,:]
    return mean.numpy(), var, cov

def GPFLOW_sample(m, xx):
    tf.random.set_seed(1)
    samples = m.predict_f_samples(xx, 10)  # (10, 100, 1)
    return samples

def GPFLOW_plot(m, X, Y, mean, var, cov, T):
    print_summary(m)
    
    xx = np.linspace(-0.1, 1.1, T).reshape(T, 1)
    samples = GPFLOW_sample(m, xx)
    
    fig = plt.figure(figsize=(12, 6))
    plt.plot(X, Y, "kx", mew=2)
    plt.plot(xx, mean, "C0", lw=2)
    plt.fill_between(xx[:, 0], mean[:, 0] - 1.96 * np.sqrt(var[:, 0]), 
                     mean[:, 0] + 1.96 * np.sqrt(var[:, 0]), color="C0", alpha=0.2)
    plt.plot(xx, samples[:, :, 0].numpy().T, "C0", linewidth=0.5)
    plt.ylim(-0.1, 1.1)
    plt.title([m.kernel.kernels[i].name for i in range(len(m.kernel.kernels))])
    plt.show()
    
    ax = sns.heatmap(cov)
    plt.show()
    
def GPFLOW_run(X, Y, kernel, meanf, params, T):
    m = GPFLOW_model(X, Y, kernel, meanf, params)
    opt_logs = GPFLOW_optimize(m)
    mean, var, cov = GPFLOW_predict(m, T)
    return m, mean, var, cov


# ========== GPYTORCH ============================


class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, kernel, meanf):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        lengthscale_prior = gpytorch.priors.GammaPrior(1, 5) # first = l (for smoothness), second = sigma (variance)
        self.mean_module = meanf
        self.covar_module = gpytorch.kernels.ScaleKernel(kernel)
    
    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
    
def GPYTORCH_prepare_data_for_model(X, Y):
    train_x, train_y = torch.squeeze(torch.tensor(X)).float(), torch.squeeze(torch.tensor(Y)).float()
    return train_x, train_y

def GPYTORCH_optimize(model, likelihood, train_x, train_y, itr=50, lr=0.1):
    likelihood.train()
    model.train()
    
    smoke_test = ('CI' in os.environ)
    training_iter = 2 if smoke_test else itr
    optimizer = torch.optim.Adam([{'params': model.parameters()}], lr) # Includes GaussianLikelihood parameters
    
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)
    for i in range(training_iter):
        optimizer.zero_grad()
        output = model(train_x)
        loss = -mll(output, train_y)
        loss.backward()
        optimizer.step()
        
    return model, likelihood

def GPYTORCH_sample(observed_pred):
    px = torch.distributions.Normal(observed_pred.mean, observed_pred.stddev)
    sample_ = px.sample()
    return sample_

def GPYTORCH_evaluate(likelihood, model, train_x, train_y, T):
    model.eval()
    likelihood.eval()

    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        test_x = torch.linspace(0, 1, T)
        out = model(test_x)
        observed_pred = likelihood(out)
        mean = observed_pred.mean
        var = observed_pred.variance
        cov = observed_pred.covariance_matrix
                
    return DF(torch.squeeze(mean).numpy()).values, var.numpy(), cov.numpy()

def GPYTORCH_plot_posterior(train_x, train_y, test_x, observed_pred, sample):
    with torch.no_grad():
        f, ax = plt.subplots(1, 1, figsize=(8, 6))
        lower, upper = observed_pred.confidence_region() # confidence bounds
        ax.plot(train_x.numpy(), train_y.numpy(), 'k*') # train data
        ax.plot(test_x.numpy(), observed_pred.mean.numpy(), 'b') # posterior mean
        #ax.plot(test_x.numpy(), sample.numpy(), 'r*') # posterior sample
        ax.fill_between(test_x.numpy(), lower.numpy(), upper.numpy(), color='y', alpha=0.3) # confidence interval
        ax.legend(['Observed Data', 'Mean', 'Sample Posterior', 'Confidence'])
        plt.title('GPyTorch mean and variance')
        plt.show()
        
def GPYTORCH_plot_cov(observed_pred):
    cov = observed_pred.covariance_matrix
    ax = sns.heatmap(cov)

def GPYTORCH_plot(model, train_x, train_y):
    T = len(train_x)
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        test_x = torch.linspace(0, 1, T)
        out = model(test_x)
        observed_pred = likelihood(out)
        sample = GPYTORCH_sample(observed_pred)
        
        GPYTORCH_plot_posterior(train_x, train_y, test_x, observed_pred, sample)
        GPYTORCH_plot_cov(observed_pred)
        
def GPYTORCH_run(train_x, train_y, kernel, meanf, likelihood):
    model                 = ExactGPModel(train_x, train_y, likelihood, kernel, meanf)
    model, likelihood     = GPYTORCH_optimize(model, likelihood, train_x, train_y)
    mean, var, cov = GPYTORCH_evaluate(likelihood, model, train_x, train_y, T)
    return model, mean, var, cov


# ========== RUN ============================


# outputs
PLOT = False
SAVE = True

PATH_output = './data/CovarianceMatrices'
plt.rcParams.update({'font.size': 18})

T = 200 # time horizon = size of covariance matrix
PACKAGE = 'GPFLOW' #['GPFLOW', 'GPYTORCH'] # GPYTORCH, GPFLOW
PACKAGE = 'GPYTORCH'

# -------------------------------------------

# GPFLOW
# kernels: RBF, Matern12, Matern32, Matern52, Exponential, Polynomial, Cosine, Linear, Constant, RationalQuadratic
#GPFLOW_kernel = gpflow.kernels.Sum([gpflow.kernels.RBF(), gpflow.kernels.Exponential(), gpflow.kernels.Cosine()])
GPFLOW_kernel = gpflow.kernels.Product([gpflow.kernels.Matern12(), gpflow.kernels.Matern32()])
GPFLOW_meanf = gpflow.mean_functions.Constant() # None, Linear, Constant
GPFLOW_params = [0.01, 0.1, 0.1] # [likelihood variance, kernel lengthscale, kernel variance]

# GPYTORCH
# kernels: MaternKernel(nu=2.5 or 1.5 or 0.5), LinearKernel, PeriodicKernel, CosineKernel, PolynomialKernel(power=20), RQKernel, SpectralMixtureKernel(num_mixtures=3)
GPYTORCH_lengthscale_prior = gpytorch.priors.GammaPrior(1, 5) # prior (lengthscale and variance)
GPYTORCH_kernel = gpytorch.kernels.CosineKernel(lengthscale_prior=GPYTORCH_lengthscale_prior) 
GPYTORCH_meanf = gpytorch.means.ConstantMean()
GPYTORCH_likelihood = gpytorch.likelihoods.GaussianLikelihood()

# -------------------------------------------


if __name__ == '__main__':
    
    # data
    #hourly_demand = load_data(PATH)['hourly demand conties (MWh)']
    substations = pd.read_excel('data/data_selected.xlsx', header=[0], sheet_name=['Substation Capacity Coordinates'])['Substation Capacity Coordinates']['Substation Nodes']
    gpr_data = {city : scale_data(DF(hourly_demand[city]).reset_index(drop=True)) for city in substations}
    sample_indices = gpr_data[next(iter(gpr_data_samples))].sample(T).index
    
    # create mean/cov file for each city in network
    for i, county in enumerate(list(gpr_data.keys())):
        
        df = gpr_data[county].loc[sample_indices]
        X_ = np.array(np.array([df['hours']]).reshape((T, 1)))
        Y_ = np.array(np.array([df['demand']]).reshape((T, 1)))
        # run
        if PACKAGE == 'GPFLOW':
            model, mean, var, cov = GPFLOW_run(X_, Y_, GPFLOW_kernel, GPFLOW_meanf, GPFLOW_params, T)
            if PLOT: GPFLOW_plot(model, X, Y, mean, var, cov, T)
        elif PACKAGE == 'GPYTORCH':
            train_x, train_y = GPYTORCH_prepare_data_for_model(X_, Y_)
            model, mean, var, cov = GPYTORCH_run(train_x, train_y, GPYTORCH_kernel, GPYTORCH_meanf, GPYTORCH_likelihood)
            if PLOT: GPYTORCH_plot(model, train_x, train_y)
            
        # save Cholesky
        L = calc_cholesky(cov)
        out = create_output_file(L, mean)
        scipy.io.savemat(str(PATH_output + '/' + PACKAGE + '_' + county + '.mat'), out)
        
        print('Round: ', str(i+1), ' of ', len(gpr_data), ', \t', county, 'done!')

