import numpy as np
import matplotlib.pyplot as plt
import sys    
import os
import xlwt
import xlrd
import pandas as pd
import sklearn.linear_model
regr = sklearn.linear_model.LinearRegression()
from scipy.optimize import curve_fit
import operator

#from unpack_and_split_data import extract_xl_BLF, yr_unwrap, split_years_unwrap, consider_lag

#-------------------------------------------------------------------------
#MODEL FUNCTIONS:
def avg_of_funcs(x,funcs):
    return sum([fun(x) for fun in funcs])/len(funcs)

#exponential and logistic
#def func_exp(x,a,b,c):
    #return a * np.exp(b * x) + c
#def func_logistic(x,c,a,b,d):
    #return c/(1+a*np.exp(b*x)) + d
def func_exp(x,a,b):
    return a * np.exp(b * x) 
def func_logistic(x,c,a,b,d):
    return c/(1+a*np.exp(b*x)) + d

import numpy, scipy.optimize
#polynomial approximations
def func_poly(x, *args):
    n, poly = 0,0
    for a in args:
        poly += a*(x**n)
        n += 1
    return poly
def func_linear(x,a0,a1):
    return func_poly(x,a0,a1)
def func_poly2(x,a0,a1,a2):
    return func_poly(x,a0,a1,a2)
def func_poly3(x,a0,a1,a2,a3):
    return func_poly(x,a0,a1,a2,a3)

def func_poly2_fixed_endpoints(xy1,xy2):
    x1,y1 = xy1[0],xy1[1]
    x2,y2 = xy2[0],xy2[1]
    def polyA(x,a):
        b = (y1-a*x1**2-y2+a*x2**2)/(x1-x2)
        c = y1-a*x1**2-b*x1
        return c+b*x+a**2
    return polyA

#-----------------------------------------------------------------------------
#PLOTTING FUNCTIONS
def model_type(funct):
    classA = [func_linear, func_poly,func_poly2,func_poly3]
    if funct in classA:
        return 'poly'
    else:
        return funct

def pre_plot(tupl,fit_type='poly'):#,y_axis=None):
    if fit_type=='poly':
        strg = 'f(x) = %5.3f + %5.3fx'
        for i in range(2, len(tupl)):
            strg += ' %5.3fx^'+str(i)
        return strg % tuple(tupl)
    elif fit_type==func_exp:
        return 'f(x) = %5.3fexp(%5.3fx)' % tuple(tupl)
    elif fit_type==func_logistic:
        return ' %5.3f/{1 + %5.3fexp(-%5.3fx)}+%5.3f' % tuple(tupl)
        
# Create linear regression object
#regr = sklearn.linear_model.LinearRegression()

def fit_2sets(X_series,Y_series, fit_func=func_linear, mask=None):
            
    X,Y = X_series, Y_series
    if type(mask)!=type(None):
        Xm = np.ma.masked_array(X,mask=mask)
        Ym = np.ma.masked_array(Y,mask=mask)
        X,Y = Xm.compressed(), Ym.compressed()
    popt, pcov = curve_fit(fit_func, X, Y)#,sigma=sigma)
    def newf(x): return fit_func(x,*popt)
    labe = pre_plot(tuple(popt),model_type(fit_func))
    dic1={}
    dic1.update({'function':newf})
    dic1.update({'parameters':popt})
    dic1.update({'print function':labe})
    return dic1

def array_span(Xob, function,dense=0,specify_points=0):
    if dense==0 and specify_points==0:
        Xspan = sorted(Xob)
        Yexp = [function(x) for x in Xob]
        return Xspan, Yexp
    else:
        pts = len(Xob) if specify_points==0 else specify_points
        x0,xN = min(Xob),max(Xob)
        Xspan = np.linspace(x0,xN,pts)
        Yexp = [function(x) for x in Xspan]
        return Xspan, Yexp        

m = [1,0,0,1,1,1,0,0,1,1,0]
x = [1,2,5,8,10,2,3,4,4,5,3]
y = [3,4,-1,9,0,1,1,2,3,4,5]

## with mask m
#dicc = fit_2sets(x,y,fit_func=func_linear,mask=m)
#fun = dicc['function']
#xs, yexp = array_span(x,fun,specify_points=20)
#Xm = np.ma.masked_array(x,mask=m)
#Ym = np.ma.masked_array(y,mask=m)
#plt.plot(Xm,Ym,'bo')
#plt.plot(xs,yexp,'r',label=dicc['print function'])

## with mask of all 1's (no values masked)
#dicc = fit_2sets(x,y,fit_func=func_linear,mask=[0,0,0,0,0,0,0,0,0,0,0])
#fun = dicc['function']
#Xm = np.ma.masked_array(x,mask=[0,0,0,0,0,0,0,0,0,0,0])
#Ym = np.ma.masked_array(y,mask=[0,0,0,0,0,0,0,0,0,0,0])
#xs, yexp = array_span(x,fun,specify_points=20)
#plt.plot(Xm,Ym,'g.')
#plt.plot(xs,yexp,'y--',label=dicc['print function'])

## with mask=None
#dicc = fit_2sets(x,y,fit_func=func_linear)
#fun = dicc['function']
#xs, yexp = array_span(x,fun,specify_points=20)
#plt.plot(x,y,'g.')
#plt.plot(xs,yexp,'y.',label=dicc['print function'])
##plot_fit(x,y,func_linear)
#plt.legend()
#plt.show()
