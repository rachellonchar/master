
#from numerical_slope import *
#import fitting_funcs as ftf
#cwd = os.getcwd()

#v,n = variables,naming
#from analysis_funcs_newdat import dict_call
#def updater(*params,normalized_to_all_years='no',stats='no'):
    #for param in params:
        #dict_call(param,variables, naming,normalized_to_all_years,stats)
    #return variables, naming

##plt.plot(v['TotDays'],v['CH4_S2'],'b.')
##plt.ylim(ymax=.068)
##plt.show()

#def threshold_mask(threshold=0,param='WT'):
    #mask_above = np.ones_like(v[param])
    #mask_below = np.ones_like(v[param])
    #above,below = 0,0
    #days_above, days_below = np.zeros_like(v[param]),np.zeros_like(v[param])
    #yr=1
    #below, above = 0,0
    #for idx in range(0, len(v[param])):
        #pt = v[param][idx]
        ##winter condition:
        #if v['Ts10'][idx]<1:
            ##winter...ground too frozen - always mask this
            #below, above = 0,0
            #if not_winter:
                #not_winter=False
        #else:
            #if v['Ts10'][idx-1]>1:
                #not_winter=True
        ##threshold condition for parameter of choice:
        #if pt>threshold:
            ##if v[param][idx-1]<=threshold:
                ##above = 0
            #mask_above[idx] = 0
            #above += 1
            #days_above[idx] = above
            #days_below[idx] = below
            ##below = 0
        #else:
            #if v[param][idx-1]>threshold:
                #below = 0
            #mask_below[idx] = 0
            #below += 1
            #days_above[idx] = above
            #days_below[idx] = below
            #above = 0
    #return mask_above, mask_below, days_above, days_below
from prelim_funcs import *

#def pot_CH(independent='NTs10',dependent='NCH4_S2',log_trans='yes'):
    #v,n = updater(dependent, independent)
    #newf,newfD,labe,popt = ftf.fit_series(v[independent],v[dependent], 
        #fit_func=ftf.func_exp,
        #period='none',ampl='none',plot_pts=50,return_series='no',
        #add_w=[],add_w2=[],specified_weights=[],
        #fixed_endpoints='no',phas='none',f=r'$CH_4$',log_trans=log_trans)
    #X = v[independent]
    #Y = v[dependent] if log_trans=='no' else [np.log(de) for de in v[dependent]]
    #model = [newf(t) for t in X]
    #all_devs = [Y[i]-model[i] for i in range(0,len(X))]
    #nolog_model = [np.exp(newf(t)) for t in X]
    #return all_devs,[X,Y,model,labe,nolog_model]

#devv = [v['CH4'][i]-Ch_T(v['Ts10'][i]) for i in range(0, len(v['Ts10']))]


#plt.plot(v['NTs10'],v['NCH4'],'ro')
#plt.plot(v['NTs10'],devv,'bo')
#plt.show()

#def plot_outliers_vslope(X,Y,deviations_predictor=None,dev_fit_type=btf.func_exp,
    #fit_type=btf.func_linear,num_of_outliers=60):
    
    #cols,rows = 1,2
    #fig, ax = plt.subplots(ncols=cols,nrows=rows,  sharex=True, sharey=True, figsize=(15,8.3))
    
    #stable_fdic,stable_mask,outs_removed = find_stable_slope(X=X,Y=Y,deviations_predictor=deviations_predictor,
        #dev_fit_type=dev_fit_type, fit_type=fit_type,num_of_outliers=num_of_outliers)
    ##if deviations_predictor!=None:
        ##Y='devs'
    #Xsil,Ysil = notation_fix(X),notation_fix(Y)
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

def subplots(*xymask_pairs,thres=10,thres_param='WT',deviations_predictor=None,def_cmap='Greys',
    CH4_use='NCH4_S1',
    pic_name=None, pic_folder_name=None,save_or_show='show',cwd=cwd):
        
    #for xym in xymask_pairs:

    #ma, mb, days_above, days_below = threshold_mask(threshold=thres,param=thres_param)
    #devv,array = pot_CH(independent='NTs10',dependent=CH4_use)
    #v.update({'days above':days_above})
    #v.update({'days below':days_below})
    #v.update({'mask none':np.zeros_like(v[thres_param])})
    #v.update({'mask above':ma})
    #v.update({'mask below':mb})
    #v.update({CH4_use+' deviations':devv})
    ##print(CH4_use+' deviations')
    #v.update({'log '+CH4_use: array[1]})
    #v.update({'predicted log '+CH4_use: array[2]})
    #v.update({'predicted '+CH4_use: array[4]})
    #labe = array[3]
    num_plots = 0
    for xx in xymask_pairs:
        if type(xx)==str:
            num_plots -= 1
        else:
            num_plots += 1
    cols = math.floor(np.sqrt(num_plots))
    rows = math.ceil(num_plots/cols)
    if cols<rows:
        f1,f2=15,8.3
    else:
        f1,f2=8.3,15
    fig, ax = plt.subplots(ncols=cols,nrows=rows,  sharex=True, sharey=True, figsize=(f1,f2))
    ni,ni_last,brg=1,1,0
    
    for i in range(0,len(xymask_pairs)):
        vec = xymask_pairs[i] if type(xymask_pairs[i-1])!= type(btf.func_linear) else new_vec
        brgy = ['b','r','g','y']
subplots(['WT','CH4'])
        ##print(vec)
        #if type(vec)==str:
            #ni -= 1
            #if vec=='hold':
                #None
            #elif vec=='linear fit':
                #newf,newfD,labe,popt = ftf.fit_series(prevX,prevY,
                    #fit_func=np.ma.polyfit,
                    #period='none',ampl='none',plot_pts=50,return_series='no',
                    #add_w=[],add_w2=[],specified_weights=[],
                    #fixed_endpoints='no',phas='none',f=xymask_pairs[i-1][1],log_trans='no',
                    #poly_order=1)
                #lin_model = [newf(x) for x in prevX]
                #vsb.update({'lin model':lin_model})
                #new_vec=[xymask_pairs[i-1][0],'lin model','old mask dummy','solid line']
            #elif vec=='hor thres':
                #plt.hlines(thres,0,v['TotDays'][-1])
                ##ni += 1
        #elif vec[0]=='skip':
            #ni += 1
        #else:
            #plt.subplot(cols,rows,ni)
            #if ni_last!=ni:
                #brg=0
            #vsb,nsb = updater(*vec)
            #X0,Y0 = vsb[vec[0]],vsb[vec[1]]
            #if len(vec)>=3:
                #if xymask_pairs[i-1]=='linear fit':
                    #None
                #else:
                    #mask = v['mask '+vec[2]]
            #else:
                #mask = v['mask none']
            #X = np.ma.masked_array(vsb[vec[0]],mask=mask)
            #Y = np.ma.masked_array(vsb[vec[1]],mask=mask)
            #if len(vec)>=4: #colorbar given
                #if vec[3]!='solid line':
                    #daysA = np.ma.masked_array(vsb[vec[3]],mask=mask)
                    #maxD, minD = np.ma.MaskedArray.max(daysA), np.ma.MaskedArray.min(daysA)
                    #norm = matplotlib.colors.Normalize(vmin=minD, vmax=maxD, clip=True)
                    #if len(vec)==5:
                        #col_map = vec[4]
                    #else:
                        #col_map = def_cmap
                    #mapper = cm.ScalarMappable(norm=norm, cmap=plt.cm.get_cmap(col_map))
                    #bill = [mapper.to_rgba(val) for val in vsb[vec[3]]]
                    #plt.scatter(X,Y,c=bill,cmap=mapper,edgecolor=bill,vmin=minD,vmax=maxD)
                    #mapper.set_array([])
                    #plt.colorbar(mapper,label=vec[3])
                #else:
                    #plt.plot(X,Y,brgy[brg],label=labe)
                    #labe=array[3]
                    #plt.legend(loc=2,ncol=1, fancybox=True,prop={'size':10})
                    #brg+=1
                    #if brg>3:
                        #brg=0
            #else:
                #if len(vec)>=3:
                    #mask_lab = vec[2]+' '+str(thres)+' cm'
                #else:
                    #mask_lab = ''
                #plt.scatter(X,Y,c=brgy[brg],edgecolor=brgy[brg],label=vec[1]+' ' +mask_lab)
                #plt.legend(loc=2,ncol=1, fancybox=True,prop={'size':10})
                #brg+=1
                #if brg>3:
                    #brg=0
            #plt.xlabel(vec[0])
            #if xymask_pairs[i-1]!='linear fit':
                #plt.ylabel(vec[1])
            #else:
                #plt.legend(loc=2,ncol=1, fancybox=True,prop={'size':10})
            #plt.grid()
            #ni += 1
            #prevX,prevY = X,Y
        #ni_last=ni
    #plt.suptitle('Threshold = '+str(thres)+' cm',fontsize=16,fontdict=font)
    ##plt.suptitle(CH4_use,fontsize=16,fontdict=font)
    #plt.tight_layout()
    ##plt.show()
    ##-------
    #if pic_name==None:
        #if save_or_show=='save':
            
            #plt.savefig(cwd+'/g/'+CH4_use, bbox_inches='tight')
            #plt.close(fig)
        #else:
            #plt.show()
    #else:
        #pnam = pic_namer(pic_name, pic_folder_name)
        #plt.savefig(pnam, bbox_inches='tight')

#C_type = ''
#C_type = '_S1'
##C_type = '_S2'

##subplots(['TotDays','NCH4'+C_type,'above'],'hold',['TotDays','predicted NCH4'+C_type,'above'],
    ##['NWT','NCH4'+C_type+' deviations','above','Ts10','viridis'],'linear fit',[],
    ###['above
    ##CH4_use='NCH4'+C_type,save_or_show='show')
#abe = 'above'
#abe2 = 'below'
#subplots(
    #['NTs10','NCH4'+C_type,abe],'hold',
    #['NTs10','NCH4'+C_type,abe2],'hold',['NTs10','predicted NCH4'+C_type,'none','solid line'],
    
    #['NTs10','log NCH4'+C_type,abe],'hold',
    #['NTs10','log NCH4'+C_type,abe2],'hold',['NTs10','predicted log NCH4'+C_type,'none','solid line'],
    
    #['NWT','NCH4'+C_type+' deviations',abe,'Ts10','viridis'],'linear fit',[],
    #thres=5,CH4_use='NCH4'+C_type,save_or_show='show')

#subplots(['TotDays','CH4'])

#plt.scatter(v['CH4_S1'],v['CH4_S2'])
#plt.show()

