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

def plot_outliers_vslope(X,Y,deviations_predictor=None,dev_fit_type=btf.func_exp,
    fit_type=btf.func_linear,num_of_outliers=60):
    
    cols,rows = 1,2
    fig, ax = plt.subplots(ncols=cols,nrows=rows,  sharex=True, sharey=True, figsize=(15,8.3))
    
    stable_fdic,stable_mask,outs_removed = find_stable_slope(X=X,Y=Y,deviations_predictor=deviations_predictor,
        dev_fit_type=dev_fit_type, fit_type=fit_type,num_of_outliers=num_of_outliers)
    #if deviations_predictor!=None:
        #Y='devs'
    Xsil,Ysil = notation_fix(X),notation_fix(Y)
    xSTAB, yexpSTAB = btf.array_span(Xsil,stable_fdic['function'],specify_points=20)
    if deviations_predictor!=None:
        dp = deviations_predictor
        typee = 'outlier removal, '+list_to_name(dp)+' vs. '+list_to_name(Y)
        out_X, out_Y = dp,Y
        Y='devs'
    else:
        outlier_loop(X=X,Y=Y,fit_type=fit_type,num_of_outliers=num_of_outliers)
        typee = 'outlier removal, '+list_to_name(X)+' vs. '+list_to_name(Y)
        out_X, out_Y = X,Y
    fd = v['mask'][typee]
    Ysil = notation_fix(Y)
    
    #plot graphs with model and stable model:
    plt.subplot(cols,rows,1)
    Xv = np.ma.masked_array(Xsil,mask=stable_mask)
    Yv = np.ma.masked_array(Ysil,mask=stable_mask)
    func_dic, toss1, toss2 = general_fit_pre(X,Y,fit_type=fit_type)
    xs, yexp = btf.array_span(Xsil,func_dic['function'],specify_points=20)
    plt.scatter(Xsil,Ysil,c='r',edgecolor='r')
    plt.scatter(Xv,Yv,c='y')
    plt.plot(xSTAB,yexpSTAB,'y--',label=stable_fdic['print function']+'\nstable at '+str(outs_removed)+' outliers removed')
    plt.plot(xs,yexp,'r--',label=func_dic['print function'])
    plt.legend(loc=2,ncol=1, fancybox=True,prop={'size':10})
    plt.xlabel(X)
    if deviations_predictor!=None:
        plt.ylabel('deviations of '+list_to_name(out_Y)+'(observed-expected)'+
            '\nbased off '+' of '+list_to_name(out_X)+' vs. '+list_to_name(out_Y))
    else:
        plt.ylabel(list_to_name(Y))
    
    #plot outliers vs slope:
    plt.subplot(cols,rows,2)
    outs,slps = [],[]
    for ii in range(1,num_of_outliers+2):
        mask_ = fd[ii-1]['mask vector']
        sl_op = general_fit_pre(X,Y,fit_type=fit_type,mask=mask_,just_slope=1)
        outs.append(ii-1)
        slps.append(sl_op)
    plt.xlabel('outliers removed')
    plt.ylabel('slope')
    plt.plot(outs,slps)
    #plt.show()
        
#plot_outliers(X='NWT',Y='devs',deviations_predictor='NTs10')

p2 = 'CH4_S1'
#plot_outliers_vslope(X=['period of inundation',0],Y=p2,deviations_predictor='NTs10')
plot_outliers_vslope(X='WT',Y=p2,deviations_predictor='NTs10')
#plot_outliers_vslope(X='NWT',Y=p2)
plt.show()

