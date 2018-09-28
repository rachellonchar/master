from prelim_funcs import *

def plot_basic(X,Y,deviations_predictor=None,dev_fit_type=btf.func_exp,
    fit_type=btf.func_linear,num_of_outliers=30):
    
    fig, ax = plt.subplots(ncols=1,nrows=1)
    stable_fdic,stable_mask,outs_removed = find_stable_slope(X=X,Y=Y,deviations_predictor=deviations_predictor,
        dev_fit_type=dev_fit_type, fit_type=fit_type,num_of_outliers=num_of_outliers)
    Xsil,Ysil = notation_fix(X),notation_fix(Y)
    xSTAB, yexpSTAB = btf.array_span(Xsil,stable_fdic['function'],specify_points=20)
    if deviations_predictor!=None:
        dp = deviations_predictor
        typee = 'outlier removal, '+list_to_name(dp)+' vs. '+list_to_name(Y)
        Y='devs'
        fd = v['mask'][typee]
    if deviations_predictor!=None:
        Y='devs'
    Xv = np.ma.masked_array(Xsil,mask=stable_mask)
    Yv = np.ma.masked_array(Ysil,mask=stable_mask)
    func_dic, toss1, toss2 = general_fit_pre(X,Y,fit_type=fit_type)
    xs, yexp = btf.array_span(Xsil,func_dic['function'],specify_points=20)
    plt.scatter(Xsil,Ysil,c='r',edgecolor='r')
    plt.scatter(Xv,Yv,c='y')
    plt.plot(xSTAB,yexpSTAB,'y--',label=stable_fdic['print function']+'\nstable at '+str(outs_removed)+' outliers removed')
    plt.plot(xs,yexp,'r--',label=func_dic['print function'])
    plt.legend(loc=2,ncol=1, fancybox=True,prop={'size':10})
    plt.tight_layout()
    plt.show()


#plot_basic(X='NWT',Y=p2,deviations_predictor='NTs10')
##plot_basic(X='NWT',Y=p2)
#plt.show()

def plot_inundation_periods(X,Y,deviations_predictor=None,dev_fit_type=btf.func_exp,
    fit_type=btf.func_linear,num_of_outliers=60,
    mask_aerated=0,mask_inundated=0,mask_threshold=0):
    
    #cols,rows = 1,2
    #fig, ax = plt.subplots(ncols=cols,nrows=rows,  sharex=True, sharey=True, figsize=(15,8.3))
    
    #stable_fdic,stable_mask,outs_removed = find_stable_slope(X=X,Y=Y,deviations_predictor=deviations_predictor,
        #dev_fit_type=dev_fit_type, fit_type=fit_type,num_of_outliers=num_of_outliers)
    ##if deviations_predictor!=None:
        ##Y='devs'
    
    Xsil,Ysil = notation_fix(X),notation_fix(Y)
    if mask_aerated==1 or mask_inundated==1:
        a_or_i = 'mask aerated events' if mask_aerated==1 else 'mask inundated events'
        mask = notation_fix(['mask',a_or_i,mask_threshold])
        Xsil = np.ma.masked_array(Xsil,mask=mask)
        Ysil = np.ma.masked_array(Ysil,mask=mask)
    #print(Xsil)
    plt.plot(Xsil,Ysil,'ro')
    plt.show()
    #xSTAB, yexpSTAB = btf.array_span(Xsil,stable_fdic['function'],specify_points=20)
    #if deviations_predictor!=None:
        #dp = deviations_predictor
        #typee = 'outlier removal, '+list_to_name(dp)+' vs. '+list_to_name(Y)
        #out_X, out_Y = dp,Y
        #Y='devs'
    #else:
        #outlier_loop(X=X,Y=Y,fit_type=fit_type,num_of_outliers=num_of_outliers)
        #typee = 'outlier removal, '+list_to_name(X)+' vs. '+list_to_name(Y)
        #out_X, out_Y = X,Y
    #fd = v['mask'][typee]
    #Ysil = notation_fix(Y)
    
    ##plot graphs with model and stable model:
    #plt.subplot(cols,rows,1)
    #Xv = np.ma.masked_array(Xsil,mask=stable_mask)
    #Yv = np.ma.masked_array(Ysil,mask=stable_mask)
    #func_dic, toss1, toss2 = general_fit_pre(X,Y,fit_type=fit_type)
    ##if fit_type==btf.func_exp: #that is if the fit for WT vs devs is exponential then we want to see that
        ##def plotfun(xx): return np.exp(func_dic['function'](xx))
    ##else:
        ##plotfun=func_dic['function']
    #plotfun=func_dic['function']
    #xs, yexp = btf.array_span(Xsil,plotfun,specify_points=20)
    #plt.scatter(Xsil,Ysil,c='r',edgecolor='r')
    #plt.scatter(Xv,Yv,c='y')
    #plt.plot(xSTAB,yexpSTAB,'y--',label=stable_fdic['print function']+'\nstable at '+str(outs_removed)+' outliers removed')
    #plt.plot(xs,yexp,'r--',label=func_dic['print function'])
    #plt.legend(loc=2,ncol=1, fancybox=True,prop={'size':10})
    #plt.xlabel(X)
    #if deviations_predictor!=None:
        #plt.ylabel('deviations of '+list_to_name(out_Y)+'(observed-expected)'+
            #'\nbased off '+' of '+list_to_name(out_X)+' vs. '+list_to_name(out_Y))
    #else:
        #plt.ylabel(list_to_name(Y))
    
    ##plot outliers vs slope:
    #plt.subplot(cols,rows,2)
    #outs,slps = [],[]
    #for ii in range(1,num_of_outliers+2):
        #mask_ = fd[ii-1]['mask vector']
        #sl_op = general_fit_pre(X,Y,fit_type=fit_type,mask=mask_,just_slope=1)
        #outs.append(ii-1)
        #slps.append(sl_op)
    #plt.xlabel('outliers removed')
    #plt.ylabel('slope')
    #plt.plot(outs,slps)
    ##plt.show()

def demoo(Xvar='CH4_S1',Yvar='WT',color='NTs10',col_map='coolwarm',threshold=7,deviations_predictor='NTs10',dev_fit_type=btf.func_exp):
    
    cols,rows = 1,2
    fig, ax = plt.subplots(ncols=cols,nrows=rows,  sharex=True, sharey=True, figsize=(18,9))
    
    #just show WT time series, colored with temp and print a straight line for threshold 
    plt.subplot(cols,rows,1)
    mask = notation_fix(['mask','mask inundated events',threshold])
    #PoI = notation_fix(['period of aeration',threshold])
    if type(deviations_predictor) != type(None):
        Y = deviations_from_fit(deviations_predictor, Yvar, fit_type=dev_fit_type,mask=None)
    else:
        Y = notation_fix(Yvar)
    X = notation_fix(Xvar)
    #Y = np.ma.masked_array(notation_fix(Yvar),mask=mask)
    ColorV = np.ma.masked_array(notation_fix(color),mask=mask)
    maxD, minD = np.ma.MaskedArray.max(ColorV), np.ma.MaskedArray.min(ColorV)
    norm = matplotlib.colors.Normalize(vmin=minD, vmax=maxD, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=plt.cm.get_cmap(col_map))
    bill = [mapper.to_rgba(val) for val in notation_fix(color)]
    plt.plot(notation_fix('TotDays'),notation_fix('WT'))
    plt.scatter(notation_fix('TotDays'),notation_fix('WT'),s=30, c=bill,cmap=mapper,
        edgecolor=bill,vmin=minD,vmax=maxD)
    plt.hlines(threshold,0,notation_fix('TotDays')[-1])
    
    plt.subplot(cols,rows,2)
    holdX,holdY = [],[]
    all_holdsX,all_holdsY = [],[]
    color_sets,setC = [],[]
    for idx in range(0,len(mask)):
        if mask[idx]==True:
            if len(holdX)!=0:
                all_holdsX.append(holdX)
                all_holdsY.append(holdY)
                color_sets.append(setC)
                holdX,holdY,setC = [],[],[]
        else:
            holdX.append(X[idx])
            holdY.append(Y[idx])
            setC.append(bill[idx])
    slopes, periods = [],[]
    color_starts, color_ends = [],[]
    for period in range(0,len(all_holdsX)):
        slop = btf.lin_fit(all_holdsX[period],all_holdsY[period])
        slopes.append(slop)
        periods.append(len(all_holdsX[period]))
        #color_starts.append(color_sets[period][0])
        color_starts.append(color_sets[period][0])
        color_ends.append(color_sets[period][-1])
    plt.scatter(periods,slopes,marker='s',s=90,c=color_ends,cmap=mapper,
        edgecolor=color_ends,vmin=minD,vmax=maxD,label='end temp')
    plt.scatter(periods,slopes,marker='>',s=60,c=color_starts,cmap=mapper,
        edgecolor=color_starts,vmin=minD,vmax=maxD,label='start temp')
    #plt.plot(periods,slopes,'ro')
    mapper.set_array([])
    plt.ylabel('slope of WT vs methane residuals\nin this inundated period')
    plt.xlabel('days inundated as determinded by the threshold '+str(threshold))
    plt.colorbar(mapper,label='temp')
    plt.legend()
    plt.show()
        
        
    #Xsil = np.ma.masked_array(Xsil,mask=mask)
    #Ysil = np.ma.masked_array(Ysil,mask=mask)
    
demoo(threshold=5)
#plot_outliers(X='NWT',Y='devs',deviations_predictor='NTs10')

p2 = 'CH4_S1'
#plot_inundation_periods(X=['period of inundation',3],Y=p2,deviations_predictor='NTs10',mask_aerated=0)
#plot_inundation_periods(X='WT',Y=p2,deviations_predictor='NTs10',mask_inundated=1)
##plot_outliers_vslope(X='WT',Y=p2,deviations_predictor='Ts10')
##plot_outliers_vslope(X='NWT',Y=p2)
#plt.show()

