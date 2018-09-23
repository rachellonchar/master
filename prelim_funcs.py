
from numerical_slope import *
import basic_fitting_funcs as btf
import imageio
cwd = os.getcwd()

v,n = variables,naming
from analysis_funcs_newdat import dict_call
def updater(*params,normalized_to_all_years='no',stats='no'):
    for param in params:
        dict_call(param,variables, naming,normalized_to_all_years,stats)
    return variables, naming

#adding a default mask of masking nothing
#naming
mask_typesV = ['outlier removal']
n.update({'mask types':mask_typesV})
defN = {}
defN.update({'outlier removal':[0]})
n.update({'mask subtypes':defN})

#mask info:
mask_typesD = {}
mask0,mask1 = {},{}
mask1.update({'rationale':'removing outliers in a model should lead to a steady slope at some number of removed outliers'})
mask0.update({'mask vector':np.zeros_like(v['Ts10'])})
mask1.update({0:mask0})
mask_typesD({'outlier removal':mask1})
v.update({'mask':mask_typesD})
#to get the mask:  
   # v['mask'][*mask type*][*ID of mask within mask type*]['mask vector']
#to get rationale behind mask type:
   # v['mask'][*mask type*]['rationale']

#to get ALL mask type names:
   # n['mask types']
#to get ALL mask SUBtype names:
   # n['mask subtypes'][*mask type*]


def threshold_mask(threshold=0,param='WT',dic=v,naming_dic=n):
    #all points start off NOT masked:
    dic, naming_dic = updater(param)
    mask_above = np.zeros_like(dic[param])
    mask_below = np.zeros_like(dic[param])
    above,below = 0,0
    days_above, days_below = np.zeros_like(dic[param]),np.zeros_like(dic[param])
    yr=1
    below, above = 0,0
    dic2 = {}
    for idx in range(0, len(dic[param])):
        pt = dic[param][idx]
        #winter condition:
        if dic['Ts10'][idx]<1:
            #winter...ground too frozen - always mask this
            below, above = 0,0
            if not_winter:
                not_winter=False
        else:
            if dic['Ts10'][idx-1]>1:
                not_winter=True
        #threshold condition for parameter of choice:
        if pt>threshold: # inundated
            if dic[param][idx-1]<=threshold:
                above = 0
            mask_above[idx] = 1 # masks the inundated pt
            above += 1
            days_above[idx] = above
            days_below[idx] = 0#below
            #below = 0
        else: #aerated
            if dic[param][idx-1]>threshold:
                below = 0
            mask_below[idx] = 1 # masks the aerated pt
            below += 1
            days_above[idx] = 0#above
            days_below[idx] = below
            #above = 0
    #days inundated/aerated:
    try:
        dic['period of inundation'].update({threshold:days_above})
        dic['period of aeration'].update({threshold:days_above})
        vec_thresholds = naming_dic['threholds']
    except:
        ta,tb = {},{}
        ta.update({threshold:days_above})
        tb.update({threshold:days_below})
        dic.update({'period of inundation':ta})
        dic.update({'period of aeration':tb})
        vec_thresholds = []
    vec_thresholds.append(threshold)
    naming_dic.update({'threholds':vec_thresholds})
    #masks:
    try:
        dum1,dum2 = {},{}
        dum1.update({'mask vector':mask_above})
        dum2.update({'mask vector':mask_below})
        dic['mask']['mask inundated events'].update({threshold:dum1})
        dic['mask']['mask aerated events'].update({threshold:dum2})
    except:
        hold1,hold2 = {},{}
        hold11,hold22 = {},{}
        hold11.update({'mask vector':mask_above})
        hold22.update({'mask vector':mask_below})
        hold1.update({threshold:hold11})
        hold2.update({threshold:hold22})
        hold1.update({'rationale':'WT may only matter when the soil is not saturated'})
        hold2.update({'rationale':'WT may only matter when the soil is not saturated'})
        dic['mask'].update({'mask inundated events':hold1})
        dic['mask'].update({'mask aerated events':hold2})
        vec_types = n['mask types']
        vec_types.append('mask inundated events')
        vec_types.append('mask aerated events')
        naming_dic.update({'mask types':vec_types})
    naming_dic['mask subtypes'].update({'mask inundated events':vec_thresholds})
    naming_dic['mask subtypes'].update({'mask aerated events':vec_thresholds})
    return dic, naming_dic
v,n = threshold_mask()

#YOU GOT TO FIX FROM HERE DOWN --------------------------------------------------------
#outlier_analysis_funcs = {}
def expected_CH4(Tsoil='NTs10',CH4='NCH4_S1',fit_type=btf.func_exp,dic=v,naming_dic=n,mask=None):
    dic,naming_dic = updater(Tsoil,CH4)
    fic = btf.fit_2sets(dic[Tsoil],dic[CH4], fit_func=fit_type, mask=mask)
    exp_fun = fic['function']
    exp_labe = fic['print function']
    popt = fic['parameters']
    dic_Ts,dic_Ts2 = fic,{}
    dicT = {}
    dicT.update({'f':dic_Ts})
    fun = exp_fun
    if fit_type == btf.func_exp:
        lin_labe = btf.pre_plot(popt) #default model is polynomial
        def lin_fun(xx): return np.log(popt[0])+popt[1]*xx
        lin_popt = tuple((np.log(popt[0]),popt[1]))
        dic_Ts2.update({'function':lin_fun})
        dic_Ts2.update({'print function':lin_labe})
        dic_Ts2.update({'parameters':lin_popt})
        dicT.update({'logf':dic_Ts2})
        fun = lin_fun
    #find furthest
    res_squared = [(dic[CH4][idx]-fun(dic[Tsoil][idx]))**2 for idx in range(0,len(dic[Tsoil]))]
    mark_outlier = np.argmax(res_squared)
    dicT.update({'outlier index':mark_outlier})
    if mask==None:
        new_mask = np.zeros_like(dic[Tsoil])
    else:
        new_mask = mask
    new_mask[mark_outlier] = 1
    outs_removed = 1-sum(new_mask)
    return dicT, outs_removed


def carbon_predictions(Tsoil='NTs10',CH4='NCH4_S1',dictt=v,naming_dic=n,mask=None):
    
    dic = {}
    dicT, outs_removed = expected_CH4(Tsoil=Tsoil,CH4=CH4,fit_type=btf.func_exp,
        dic=dic,naming_dic=naming_dic,mask=mask)
    dic_Ts, dic_Ts2 = dicT['f'], dicT['logf']
    try:
        dic['f'][CH4].update({Tsoil:dic_Ts})
        dic['logf'][CH4].update({Tsoil:dic_Ts2})
    except:
        dic_CH,dic_CH2 = {},{}
        dic_CH.update({Tsoil:dic_Ts})
        dic_CH2.update({Tsoil:dic_Ts2})
        try:
            dic['f'].update({CH4:dic_CH})
            dic['logf'].update({CH4:dic_CH2})
        except:
            dic_f,dic_logf = {},{}
            dic_f.update({CH4:dic_CH})
            dic_logf.update({CH4:dic_CH2})
            dic.update({'f':dic_f})
            dic.update({'logf':dic_logf})
            try:
                f_types = naming_dic['function types']
            except:
                f_types = []
            f_types.append('f')
            f_types.append('logf')
            naming_dic.update({'function types':f_types})
    return dic,naming_dic

#no mask for Ts10 vs. CH4 curve fit:
#
v,n = carbon_predictions(Tsoil='NTs10',CH4='NCH4_S1',dic=v,naming_dic=n)

def general_fit_pre(X,Y,fit_type=btf.func_linear,mask=None):
    v,n = updater(X,Y)
    func_fit_dic, outs_removed = expected_CH4(Tsoil=X,CH4=Y,fit_type=fit_type,
        dic=v,naming_dic=n,mask=mask)
    mark = func_fit_dic['outlier index']
    if mask==None:
        new_mask = np.zeros_like(dic[Tsoil])
    else:
        new_mask = mask
    new_mask[mark] = 1
    popt = func_fit_dic['logf']['parameters'] if (fit_type == btf.func_exp) else func_fit_dic['f']['parameters']
    if fit_type==btf.func_linear or fit_type==btf.func_exp:
        def slopef(xx): return popt[1]
    elif fit_type==btf.func_poly2:
        def slopef(xx): return popt[1]+2*popt[2]*xx
    elif fit_type==btf.func_poly3:
        def slopef(xx): return popt[1]+2*popt[2]*xx+3*popt[3]*xx**2
    return outs_removed, slopef, new_mask

def linear_deviations(X,Y,dic=v,naming_dic=n,fit_type=btf.func_linear,mask=None):
    #check if X,Y linear (or exponential) model already exists 
    if mask==None:
        try:
            
    
    #if fit_type==btf.func_linear or fit_type==btf.func_exp:
    dic,naming_dic = updater(X,Y)
    dicT, outs_removed = expected_CH4(Tsoil=X,CH4=Y,fit_type=fit_type,
        dic=dic,naming_dic=naming_dic,mask=mask)

def general_fit(X,Y,fit_type=btf.func_linear,mask=None):
        

#def out_slope_plot(X,Y,outliers_removed=10):
    #Xn,Yn = X,Y
    #slopes = []
    #for ou in range(0, outliers_removed):
        #Xx,Yy,slop = pot_CH(Xn,Yn)
        #slopes.append(slop)
        #Xn,Yn=Xx,Yy
    #plt.plot(np.arange(1,outliers_removed+1),slopes,'bo')
    #plt.show()
   
#v,n = updater('NWT')
#expect = [fun['NCH4_S1']['NTs10'](x) for x in v['NTs10']]
#devs = [v['NCH4_S1'][idx]-expect[idx] for idx in range(0, len(expect))] #observed - model
##out_slope_plot(v['WT'],devs)
#plt.plot(v['NWT'],devs,'ro')
#newf,newfD,labe,popt = ftf.fit_series(v['WT'],devs, 
    #fit_func=ftf.func_linear, return_series='no',log_trans='no')
#plt.plot(v['NWT'],[newf(wt) for wt in v['WT']],'g',label=labe)
#plt.legend()

#plt.plot(v['TotDays'],v['NCH4_S1'],'r.')
#plt.plot(v['TotDays'],expect,'g.')
#plt.show()


#plt.plot(v['Ts10'],v['CH4_S1'],'y.')
#plt.plot(X,Y,'g.')
#plt.show()
#def pickout_outlier1(X,Y):
    #logF, regF = pot_CH(X,Y)

#fun = pot_CH()
#plt.plot(v['NTs10'],[fun['log']['NCH4_S2']['NTs10'](x) for x in v['NTs10']],'b')
#plt.plot(v['NTs10'],[fun['NCH4_S2']['NTs10'](x) for x in v['NTs10']],'r')
#plt.show()

    
    
    #function_dic.update({
    #X = v[independent]
    #Y = v[dependent] if log_trans=='no' else [np.log(de) for de in v[dependent]]
    #model = [newf(t) for t in X]
    #all_devs = [Y[i]-model[i] for i in range(0,len(X))]
    #nolog_model = [np.exp(newf(t)) for t in X]
    #return all_devs,[X,Y,model,labe,nolog_model]

