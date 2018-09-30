from prelim_funcs import *

def slope_at_thresholds(thresholds_array=[0,5,7,8,10],Xvar='period of aeration',Yvar='CH4_S1',color='NTs10',col_map='coolwarm',
    deviations_predictor='NTs10',dev_fit_type=btf.func_exp,mask_events='inundated',show_masked_pts=0):
        
    #thresholds_array = np.append(-5,thresholds_array)
    if Xvar[:6]!='period':
        thresholds_array = np.append(max(v['WT'])+1,thresholds_array) if mask_events=='inundated' else np.append(min(v['WT'])-1,thresholds_array)
    if Xvar=='period of aeration':
        mask_events = 'inundated'
    if Xvar=='period of inundation':
        mask_events = 'aerated'
    t_e = 'aerated' if mask_events=='inundated' else 'inundated'
    num_plots=len(thresholds_array)
    cols = math.floor(np.sqrt(num_plots))
    rows = math.ceil(num_plots/cols)
    if cols<rows:
        f1,f2=11.7,8.3
    else:
        f1,f2=8.3,11.7
    #print(rows,cols)
    ColorV = notation_fix(color)
    maxD, minD = max(ColorV), min(ColorV)
    norm = matplotlib.colors.Normalize(vmin=minD, vmax=maxD, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=plt.cm.get_cmap(col_map))
    mapper.set_array([])
    if type(deviations_predictor) != type(None):
        Y = deviations_from_fit(deviations_predictor, Yvar, fit_type=dev_fit_type,mask=None)
    else:
        Y = notation_fix(Yvar)
    minY,maxY = min(Y),max(Y)
    fig, ax = plt.subplots(ncols=cols,nrows=rows,  sharex=True, sharey=True, figsize=(f1,f2))
    ct = 1
    for threshold in thresholds_array:
        plt.subplot(cols,rows,ct)
        if mask_events=='inundated' or mask_events=='aerated':
            mask,dic,ndic = notation_fix(['mask','mask '+mask_events+' events',threshold],ret_dics=1)
        else:
            mask = [0 for idx in range(0,len(Y))] #masks nothing
            dic,ndic = v,n
        def nf(string): return notation_fix(string,dic=dic,naming_dic=ndic)
        if Xvar=='period of aeration':
            Xvar = ['period of aeration',threshold]
        elif Xvar=='period of inundation':
            Xvar = ['period of inundation',threshold]
        fun, print_fun = btf.lin_fit(nf(Xvar),Y,mask=mask,type_return='function and print')
        Xexp,Yexp = btf.array_span(nf(Xvar), fun,dense=1,specify_points=20)
        if show_masked_pts==1:
            plt.scatter(nf(Xvar),Y,c='white',edgecolor='darkgrey')
        
        X = np.ma.masked_array(nf(Xvar),mask=mask)
        Ymsk = np.ma.masked_array(Y,mask=mask)
        plt.plot(Xexp,Yexp,'r',label=print_fun)
        col_msk = np.ma.masked_array(notation_fix(color),mask=mask)
        bill = [mapper.to_rgba(val) for val in col_msk.compressed()]
        plt.scatter(X.compressed(),Ymsk.compressed(),s=30,c=bill,cmap=mapper,
            edgecolor=bill,vmin=minD,vmax=maxD,label='all '+t_e+' events')
        plt.legend(loc=2,ncol=1, fancybox=True,prop={'size':10})
        #if ct>1:
        plt.title('threshold='+str(threshold))
        if type(Xvar)==str and (Xvar=='WT' or Xvar=='NWT'):
            plt.vlines(threshold,minY,maxY)
        #else:
            #plt.title('no masking')
        plt.ylim(ymin=minY,ymax=maxY)
        if type(Xvar)==list:
            if Xvar[0]=='period of aeration':
                plt.xlim(xmin=0,xmax=41)
        else:
            plt.xlim(xmin=min(Xexp),xmax=max(Xexp))
        plt.colorbar(mapper,label='soil temp')
        ct+=1
    fig.text(0.5, 0.04,list_to_name(Xvar,base=1)+' of '+t_e+' events', ha='center',fontdict=font)
    fig.text(0.04, 0.5, 'CH4 residuals (based on soil temp at -10 cm)', va='center', rotation='vertical',fontdict=font)
    plt.suptitle('Sensitivity to water table at different threshold definitions of '+t_e+' events',fontsize=16,fontdict=font)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9,bottom=.1,left=.13)
    plt.grid()
    plt.show()
        
    #Xsil = np.ma.masked_array(Xsil,mask=mask)
    #Ysil = np.ma.masked_array(Ysil,mask=mask)
    
slope_at_thresholds(Xvar='WT')

#slope_per_period(threshold=15)
#plot_outliers(X='NWT',Y='devs',deviations_predictor='NTs10')

#p2 = 'CH4_S1'
#plot_inundation_periods(X=['period of inundation',3],Y=p2,deviations_predictor='NTs10',mask_aerated=0)
#plot_inundation_periods(X='WT',Y=p2,deviations_predictor='NTs10',mask_inundated=1)
##plot_outliers_vslope(X='WT',Y=p2,deviations_predictor='Ts10')
##plot_outliers_vslope(X='NWT',Y=p2)
#plt.show()

