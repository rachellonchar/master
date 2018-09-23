
from numerical_slope import *
import fitting_funcs as ftf
import imageio
cwd = os.getcwd()

v,n = variables,naming
from analysis_funcs_newdat import dict_call
def updater(*params,normalized_to_all_years='no',stats='no'):
    for param in params:
        dict_call(param,variables, naming,normalized_to_all_years,stats)
    return variables, naming

def threshold_mask(threshold=0,param='WT'):
    #all points start off NOT masked:
    mask_above = np.zeros_like(v[param])
    mask_below = np.zeros_like(v[param])
    above,below = 0,0
    days_above, days_below = np.zeros_like(v[param]),np.zeros_like(v[param])
    yr=1
    below, above = 0,0
    dic2 = {}
    for idx in range(0, len(v[param])):
        pt = v[param][idx]
        #winter condition:
        if v['Ts10'][idx]<1:
            #winter...ground too frozen - always mask this
            below, above = 0,0
            if not_winter:
                not_winter=False
        else:
            if v['Ts10'][idx-1]>1:
                not_winter=True
        #threshold condition for parameter of choice:
        if pt>threshold: # inundated
            if v[param][idx-1]<=threshold:
                above = 0
            mask_above[idx] = 1 # masks the inundated pt
            above += 1
            days_above[idx] = above
            days_below[idx] = 0#below
            #below = 0
            #print(pt)
        else: #aerated
            if v[param][idx-1]>threshold:
                below = 0
            mask_below[idx] = 1 # masks the aerated pt
            below += 1
            days_above[idx] = 0#above
            days_below[idx] = below
            #above = 0
    dic2.update({'period of inundation':days_above})
    dic2.update({'period of aeration':days_below})
    dic2.update({'mask aerated events':mask_below})
    dic2.update({'mask inundated events':mask_above})
    return dic2

empty_fun_dic = {}
empty_fun_dic.update({'IDVs':[]})
empty_fun_dic.update({'DVs':[]})
empty_fun_dic.update({'log':{}})
def pot_CH(funD=empty_fun_dic,independent='NTs10',dependent='NCH4_S2'):
    v,n = updater(dependent, independent)
    newf,newfD,labe,popt = ftf.fit_series(v[independent],v[dependent], 
        fit_func=ftf.func_exp, return_series='no',log_trans='yes')
    def regf(Ts): return np.exp(newf(Ts))
    IDVs = funD['IDVs']
    DVs = funD['DVs']
    IDVs.append(independent)
    DVs.append(dependent)
    funD.update({'IDVs':IDVs})
    funD.update({'DVs':DVs})
    hold, hold2 = {},{}
    hold.update({independent:regf})
    hold2.update({independent:newf})
    funD.update({dependent:hold})
    funD['log'].update({dependent:hold2})
    return funD

fun = pot_CH()
plt.plot(v['NTs10'],[fun['log']['NCH4_S2']['NTs10'](x) for x in v['NTs10']],'b')
plt.plot(v['NTs10'],[fun['NCH4_S2']['NTs10'](x) for x in v['NTs10']],'r')
plt.show()

    
    
    #function_dic.update({
    #X = v[independent]
    #Y = v[dependent] if log_trans=='no' else [np.log(de) for de in v[dependent]]
    #model = [newf(t) for t in X]
    #all_devs = [Y[i]-model[i] for i in range(0,len(X))]
    #nolog_model = [np.exp(newf(t)) for t in X]
    #return all_devs,[X,Y,model,labe,nolog_model]

