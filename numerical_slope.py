import numpy as np
import matplotlib.pyplot as plt
import sys    
import os
import xlwt
import xlrd
import pandas as pd
from xlrd.sheet import ctype_text 
import pickle
import sklearn.linear_model
regr = sklearn.linear_model.LinearRegression()
from scipy.optimize import curve_fit, fsolve
import operator
import scipy.stats
#from scipy.stats import norm
from numpy import linspace
from pylab import plot,show,hist,figure,title
import matplotlib
import matplotlib.cm as cm
import math
from scipy.signal import argrelextrema, savgol_filter

import sys
sys.path.insert(0, '/home/rachellonchar/Dropbox/python_work/00py_projects/peatland_project/marcell_data/py_code')
sys.path.insert(0, '/home/rachellonchar/Dropbox/python_work/00py_projects/peatland_project/marcell_data/py_code/obj')
parent_fold = '/home/rachellonchar/Dropbox/python_work/00py_projects/peatland_project/marcell_data/py_code/obj/'
cwd = os.getcwd()
def load():
    params = load_obj('new_parameters',parent_fold)
    names = load_obj('new_naming',parent_fold)
    return params, names
from helpful_functions import extract_xl, urlify,save_obj,load_obj,pic_namer
import fitting_funcs as ftf
#-----------------------------------------------------------------------
variables,naming = load()
all_years = naming['years']
from analysis_funcs_newdat import dict_call
def updater(*params,normalized_to_all_years='no',stats='no'):
    for param in params:
        dict_call(param,variables, naming,normalized_to_all_years,stats)
    return variables, naming
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
#-----------------------------------------------------------------------
def slope_2pts(xy1,xy2):
    x1,y1 = xy1[0],xy1[1]
    x2,y2 = xy2[0],xy2[1]
    m = (y2-y1)/(x2-x1)
    return m
def line_2pts(xy1,xy2,mb='no'):
    m = slope_2pts(xy1,xy2)
    b = xy1[1]-m*xy1[0]
    def linmxb(x): return m*x+b
    if mb=='yes':
        return m,b
    else:
        return linmxb
#-----------------------------------------------------------------------

def splitter(*others,main_param,years=all_years,x_axis='DoY',normalized_to_all_years='yes',
    extrema_tolerance=3,smooth_window=[19,5]):
    
    param=main_param
    years=years if type(years)==list else [years]
    v,n = updater(*others,param,normalized_to_all_years=normalized_to_all_years,
        stats='yes')
    maxpY = v['reference measures'][param]['max']
    minpY = v['reference measures'][param]['min']
    dic = {}
    slopes=[]
    split_slope={}
    for y in years:
        s,f = v['year marks'][y][0],v['year marks'][y][1]
        X,Y=v[x_axis][s:f],v[param][s:f]
        #smooth function:
        Ysm = savgol_filter(Y, smooth_window[0], smooth_window[1])
        max_id,min_id = np.argmax(Y),np.argmin(Y)
        maxp,minp=max(Ysm),min(Ysm)
        j=0
        loc_max = argrelextrema(np.array(Ysm), np.greater_equal,order=extrema_tolerance)[0]
        loc_min = argrelextrema(np.array(Ysm), np.less_equal,order=extrema_tolerance)[0]
        comb=sorted(np.append(loc_max,loc_min))
        tru_es = comb
        if 0 not in comb:
            comb=np.append(comb,[0])
        if len(Y)-1 not in comb:
            comb=np.append(comb,[len(Y)-1])
        exr = [e for e in comb if (e>7 and e<len(Y)-1-7) or (e==0) or (e==len(Y)-1)]
        extremas = sorted(exr)
        for idx in range(0,len(extremas[1:])):
            i1,i2 = extremas[idx-1],extremas[idx]
            xy1=[X[i1],Y[i1]]
            xy2=[X[i2],Y[i2]]
            slopes.append(slope_2pts(xy1,xy2))
        #----------------------------------------------------
        temp={}
        temp.update({'extremas':tru_es})
        temp.update({'adjusted extremas':extremas})
        temp.update({'local maximums':loc_max})
        temp.update({'local minimums':loc_min})
        split_slope.update({y:temp})
        #dic.update({'linear split':split_slope})
    dic.update({'linear split':split_slope})
    smx,smn = max(slopes),min(slopes)
    dic.update({'max slope':smx})
    dic.update({'min slope':smn})
    dic.update({'unique slopes all years':slopes})
    return dic,v,n
    
#------------------------------------------------------------------
#year dependent and max/min slope dependent returns 
#e.g. depend on smx, smn

def splitter_expanded(*others,split_by_param,years=all_years,x_axis='DoY',
    normalized_to_all_years='yes',extrema_tolerance=3,
    smooth_window=[11,6]):
    
    param=split_by_param
    model_type=ftf.func_linear
    years=years if type(years)==list else [years]
    dic,v,n = splitter(*others,main_param=param,years=years,x_axis=x_axis,
        normalized_to_all_years=normalized_to_all_years,extrema_tolerance=extrema_tolerance,
        smooth_window=smooth_window)
    for y in years:
        s,f = v['year marks'][y][0],v['year marks'][y][1]
        X,Y=v[x_axis][s:f],v[param][s:f]
        #X,Y=v[y][x_axis],v[y][param]
        extremas=dic['linear split'][y]['adjusted extremas']
        tru_extremas=dic['linear split'][y]['extremas']
        #if y==2010:
            #print(extremas,tru_extremas)
        slopesY,lines,line_funcs = [52 for y in Y],[52 for y in Y],[]
        j=1
        iv1=extremas[0]
        info, info_idx = {},[]
        split_srs=[]
        for idx in extremas[1:]:
            iv2=idx
            piecewise={}
            piecewise.update({'index range':[iv1,iv2]})
            piecewise.update({'x range (actual domain)':[X[iv1],X[iv2]]})
            #endpoints:
            funl = line_2pts([X[iv1],Y[iv1]],[X[iv2],Y[iv2]])
            mcon,bb = line_2pts([X[iv1],Y[iv1]],[X[iv2],Y[iv2]],mb='yes')
            fpr = ' %5.3f + %5.3fx'%tuple([bb,mcon])
            #def m(x): return mcon
            split_srs.append([X[iv1:iv2+1],Y[iv1:iv2+1]])
            funn={}
            funn.update({'function':funl})
            funn.update({'pretty print':param+'(t) ='+fpr})
            funn.update({'slope':mcon})
            funn.update({'y intercept':bb})
            piecewise.update({'model':funn})
            piecewise.update({'split':[X[iv1:iv2+1],Y[iv1:iv2+1]]})
            def splitfun(other_param):
                Z = v[other_param][s:f]
                return Z[iv1:iv2+1]
            for ot in others:
                piecewise.update({ot:splitfun(ot)})
            info.update({j:piecewise})
            info_idx.append(j)
            iv1=idx
            j+=1
        dic['linear split'][y].update({'index':info})
        dic['linear split'][y].update({'list indices':info_idx})
        #dic[y].update({'xy split series':split_srs})
    return dic,v,n

#from analysis_funcs_newdat import *
def fitter(param,split_by_param='NTs10',years=all_years,x_axis='DoY',split_x_axis='DoY',
    normalized_to_all_years='yes',extrema_tolerance=12,color_scheme=None,
    smooth_window=[19,5],id_by='val',subs='yes',marker_set=None,
    pic_name=None, pic_folder_name=None,save_or_show='show',cwd=cwd):
        
        
    #--------
    #zero_tol=0.04
    #----------------------------------------------------
    if color_scheme==None:
        if split_by_param=='WT' or split_by_param=='NWT':
            color_scheme='BrBG'
        elif split_by_param[:3]=='NTs' or split_by_param[:2]=='Ts':
            color_scheme='coolwarm'
        else:
            color_scheme='Greys'
    param_id=split_by_param
    num_plots=len(years)
    cols = math.floor(np.sqrt(num_plots))
    rows = math.ceil(num_plots/cols)
    if cols<rows:
        f1,f2=11.7,8.3
    else:
        f1,f2=8.3,11.7
    #print(rows,cols)
    if subs=='yes':
        fig, ax = plt.subplots(ncols=cols,nrows=rows,  sharex=True, sharey=True, figsize=(11.7,8.3))
    else:
        fig, ax = plt.subplots(figsize=(11.7,8.3))
    #----------------------------------------------------
    years=years if type(years)==list else [years]
    d,v,n = splitter_expanded(param,x_axis,split_by_param=param_id,years=years,x_axis=split_x_axis,
        normalized_to_all_years=normalized_to_all_years,extrema_tolerance=extrema_tolerance,
        smooth_window=smooth_window)
    xlow,xhigh = min(v[x_axis]),max(v[x_axis])
    ylow,yhigh = min(v[param]),max(v[param])
    if id_by=='slope':
        bigg = max([abs(d['max slope']),abs(d['min slope'])])
        smx,smn = bigg, -1*bigg
    else:
        smx0=v['reference measures'][param_id]['max']
        smn0=v['reference measures'][param_id]['min']
        #
        if split_by_param=='WT' or split_by_param=='NWT':
            if marker_set==None:
                bigg = max([abs(smx0),abs(smn0)])
                smx,smn = bigg, -1*bigg
                #print(v[split_by_param+' func'](10))
            else:
                if split_by_param=='NWT':
                    offs = v[split_by_param+' func'](marker_set/100)
                else:
                    offs = marker_set/100
                bigg = max([abs(smx0-abs(offs)),abs(smn0-abs(offs))])
                smx,smn = offs+bigg, offs-1*bigg
                #v[split_by_param+' func']
        else:
            smx,smn=smx0,smn0
    norm = matplotlib.colors.Normalize(vmin=smn, vmax=smx, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=plt.cm.get_cmap(color_scheme))
    yid = 0
    for yr in years:
        if subs=='yes':
            plt.subplot(cols,rows,yid+1)
        for j in d['linear split'][yr]['list indices']:
            Z = d['linear split'][yr]['index'][j][param]
            Xax = d['linear split'][yr]['index'][j][x_axis]
            #ar = d['linear split'][yr]['index'][j][x_axis]
            ar = d['linear split'][yr]['index'][j]['split']
            ar2 = d['linear split'][yr]['index'][j]['model']['function']
            ar3 = d['linear split'][yr]['index'][j]['model']['slope']
            slop_ar = [ar3 for i in ar]
            val_ar = [ar2(aa) for aa in ar[0]]
            if id_by=='slope':
                bill = [mapper.to_rgba(v) for v in slop_ar]
                plt.scatter(Xax,Z,c=mapper.to_rgba(slop_ar),marker='o',s=50,linewidth=0.5,
                    edgecolor='none',vmin=smn,vmax=smx)#,s=[(y+1)*10 for y in Y])
            else:
                bill = [mapper.to_rgba(v) for v in val_ar]
                plt.scatter(Xax,Z,c=bill,cmap=mapper,marker='o',linewidth=0.5,
                    edgecolor='none',vmin=smn,vmax=smx)#,s=[(y+1)*10 for y in Y])
            if param==param_id and id_by=='slope' and x_axis==split_x_axis:
                if ar3<=0:
                    coL='y'
                else:
                    coL='g'
                plt.plot([ar[0][0],ar[0][-1]],[Z[0],Z[-1]],coL,alpha=.5)
        ##----------------------------------------------------
        if subs=='yes' or yr==years[-1]:
            plt.xlim(xmin=xlow,xmax=xhigh)
            plt.ylim(ymin=ylow,ymax=yhigh)
            if subs=='yes':
                plt.title(years[yid])
            mapper.set_array([])
            plt.colorbar(mapper,label='slope between extrema')
        yid+=1
    fig.text(0.5, 0.04,n['full names'][x_axis], ha='center',fontdict=font)
    fig.text(0.04, 0.5, n['full names'][param], va='center', rotation='vertical',fontdict=font)
    plt.suptitle(n['full names'][param]+' vs. '+n['full names'][x_axis],fontsize=16,fontdict=font)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9,bottom=.1,left=.13)
    plt.grid()
    #-------
    if pic_name==None:
        if save_or_show=='save':
            foldID = 'slope_split' if id_by=='slope' else 'val_split'
            if subs=='yes':
                if marker_set==None or (split_by_param!='WT' and split_by_param!='NWT'):
                    plt.savefig(cwd+'/g/'+foldID+'/subs_split('+split_by_param+')_'+x_axis+'_'+param, bbox_inches='tight')
                else:
                    plt.savefig(cwd+'/g/'+foldID+'/subs_split('+split_by_param+str(marker_set)+')_'+x_axis+'_'+param, bbox_inches='tight')
            else:
                if marker_set==None or (split_by_param!='WT' and split_by_param!='NWT'):
                    plt.savefig(cwd+'/g/'+foldID+'/all_split('+split_by_param+')_'+x_axis+'_'+param, bbox_inches='tight')
                else:
                    plt.savefig(cwd+'/g/'+foldID+'/all_split('+split_by_param+str(marker_set)+')_'+x_axis+'_'+param, bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show()
    else:
        pnam = pic_namer(pic_name, pic_folder_name)
        plt.savefig(pnam, bbox_inches='tight')
    
#fitter('NCH4',x_axis='NWT',split_by_param='NTs10',marker_set=0,years=all_years,
    #subs='no',id_by='val',extrema_tolerance=10,
    #smooth_window=[19,5],save_or_show='save')


#splitter('NTs10')
      
