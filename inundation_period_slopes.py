from prelim_funcs import *

def plot_inundation_periods(X,Y,deviations_predictor=None,dev_fit_type=btf.func_exp,
    fit_type=btf.func_linear,num_of_outliers=60,
    mask_aerated=0,mask_inundated=0,mask_threshold=0):
    Xsil,Ysil = notation_fix(X),notation_fix(Y)
    if mask_aerated==1 or mask_inundated==1:
        a_or_i = 'mask aerated events' if mask_aerated==1 else 'mask inundated events'
        mask = notation_fix(['mask',a_or_i,mask_threshold])
        Xsil = np.ma.masked_array(Xsil,mask=mask)
        Ysil = np.ma.masked_array(Ysil,mask=mask)
    #print(Xsil)
    plt.plot(Xsil,Ysil,'ro')
    plt.show()

def slope_per_period(Xvar='CH4_S1',Yvar='WT',color='NTs10',col_map='coolwarm',threshold=7,
    deviations_predictor='NTs10',dev_fit_type=btf.func_exp):
    
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

def slope_at_thresholds(thresholds_array=[0,3,5,7,8,10],Xvar='period of aeration',Yvar='CH4_S1',color='NTs10',col_map='coolwarm',
    deviations_predictor='NTs10',dev_fit_type=btf.func_exp,mask_events='inundated'):
        
    #thresholds_array = np.append(-5,thresholds_array)
    t_e = 'aerated' if mask_events=='inundated' else 'inundated'
    num_plots=len(thresholds_array)
    cols = math.floor(np.sqrt(num_plots))
    rows = math.ceil(num_plots/cols)
    if cols<rows:
        f1,f2=11.7,8.3
    else:
        f1,f2=8.3,11.7
    #print(rows,cols)
    fig, ax = plt.subplots(ncols=cols,nrows=rows,  sharex=True, sharey=True, figsize=(f1,f2))
    ct = 1
    for threshold in thresholds_array:
        plt.subplot(cols,rows,ct)
        #mask = notation_fix(['mask','mask '+mask_events+' events',threshold]) #if ct>1 else None #np.zeros_like(notation_fix(Xvar))
        
    #PoI = notation_fix(['period of aeration',threshold])
        if type(deviations_predictor) != type(None):
            Y = deviations_from_fit(deviations_predictor, Yvar, fit_type=dev_fit_type,mask=None)
        else:
            Y = notation_fix(Yvar)
        
        #added in meeting:
        mask = [0 if notation_fix(['period of aeration',threshold])[idx]>10 else 1 for idx in range(0,len(Y))]
        newp = np.ma.masked_array(notation_fix(['period of aeration',threshold]),mask=mask)
        newy = np.ma.masked_array(Y,mask=mask)
        plt.plot(notation_fix(['period of aeration',threshold]),Y,'ro')
        plt.plot(newp,newy,'bo')
        fun, print_fun = btf.lin_fit(notation_fix(['period of aeration',threshold]),Y,mask=mask,type_return='function and print')
        Xexp,Yexp = btf.array_span(notation_fix(['period of aeration',threshold]), fun,dense=1,specify_points=20)
        plt.plot(Xexp,Yexp,'r',label=print_fun)
        plt.legend(loc=2,ncol=1, fancybox=True,prop={'size':10})
        
        ##print(btf.lin_fit(notation_fix(Xvar),Y,mask=mask,type_return='function and print'))
        #Xvar = ['period of aeration',threshold]
        #fun, print_fun = btf.lin_fit(notation_fix(Xvar),Y,mask=mask,type_return='function and print')
        ##print(mask)
        #Xexp,Yexp = btf.array_span(notation_fix(Xvar), fun,dense=1,specify_points=20)
        #X = np.ma.masked_array(notation_fix(Xvar),mask=mask)
        #Y = np.ma.masked_array(Y,mask=mask)
        ##plt.plot(X,Y,'ro')
        #plt.plot(Xexp,Yexp,'r',label=print_fun)
        ##ColorV = np.ma.masked_array(notation_fix(color),mask=mask)
        #ColorV = notation_fix(color)
        ##maxD, minD = np.ma.MaskedArray.max(ColorV), np.ma.MaskedArray.min(ColorV)
        #maxD, minD = max(ColorV), min(ColorV)
        #norm = matplotlib.colors.Normalize(vmin=minD, vmax=maxD, clip=True)
        #mapper = cm.ScalarMappable(norm=norm, cmap=plt.cm.get_cmap(col_map))
        #bill = [mapper.to_rgba(val) for val in notation_fix(color)]
        #plt.scatter(X,Y,s=30,c=bill,cmap=mapper,
            #edgecolor=bill,vmin=minD,vmax=maxD,label='all '+t_e+' events')
        #plt.legend(loc=2,ncol=1, fancybox=True,prop={'size':10})
        #mapper.set_array([])
        #if ct>1:
        plt.title('threshold='+str(threshold))
        plt.vlines(threshold,-1,1)
        #else:
            #plt.title('no masking')
        ##plt.ylabel('CH4 residuals (based on soil temp at -10 cm)')
        ##plt.xlabel('Water table of '+t_e+' events')
        ##plt.ylim(ymin=-0.5,ymax=3.5)
        ##plt.xlim(xmin=-.02,xmax=.16)
        #plt.colorbar(mapper,label='soil temp')
        ct+=1
    #fig.text(0.5, 0.04,'Water table of '+t_e+' events', ha='center',fontdict=font)
    #fig.text(0.04, 0.5, 'CH4 residuals (based on soil temp at -10 cm)', va='center', rotation='vertical',fontdict=font)
    #plt.suptitle('Sensitivity to water table at different threshold definitions of '+t_e+' events',fontsize=16,fontdict=font)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9,bottom=.1,left=.13)
    plt.grid()
    plt.show()
        
    #Xsil = np.ma.masked_array(Xsil,mask=mask)
    #Ysil = np.ma.masked_array(Ysil,mask=mask)
    
#slope_at_thresholds()
slope_per_period(threshold=15)
#plot_outliers(X='NWT',Y='devs',deviations_predictor='NTs10')

p2 = 'CH4_S1'
#plot_inundation_periods(X=['period of inundation',3],Y=p2,deviations_predictor='NTs10',mask_aerated=0)
#plot_inundation_periods(X='WT',Y=p2,deviations_predictor='NTs10',mask_inundated=1)
##plot_outliers_vslope(X='WT',Y=p2,deviations_predictor='Ts10')
##plot_outliers_vslope(X='NWT',Y=p2)
#plt.show()

