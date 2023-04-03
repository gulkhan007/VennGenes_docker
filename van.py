#!/usr/bin/python
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_unweighted
from matplotlib import pyplot as plt
#%matplotlib inline
import os
import pandas as pd
def Get_VanData():
    dko =pd.read_csv('./data/DE_G0CA_WT.csv', sep=',', index_col=0)
    print(dko.head())
    hdac = pd.read_csv('./data/DE_G2CA_WT.csv', sep=',', index_col=0)
    print(hdac.head())
    gps2 = pd.read_csv('./data/GoCA_G2CA_LFC.csv', sep=',', index_col=0)
    print(gps2.head())
    
    glist_dict = {'GOCAWT':list(dko.index),'GoCAG2CA':list(gps2.index),'G2CAWT':list(hdac.index)}
    
    gdata_dict = {'GOCAWT':dko,'GoCAG2CA':gps2,'G2CAWT':hdac}
    
    return glist_dict,gdata_dict
def display_dict_pandas(data_p):
    for item in data_p:
        print(f'{item} : {(data_p[item].shape)}')
def display_dict_list(g_list):
    for item in g_list:
        print(f'{item} : {len(g_list[item])}')
    #set(gdata[item])
    
def Van_Prepration(glist, gdata, title="",wflag=False):
    i=0
    tmp=""
    for item in glist:
        #print(item)

        if i==0:
            tmp=set(glist[item])
            i=1
        else:
            tmp = set.intersection(tmp,set(glist[item]))
    glist['common'] = list(tmp)
    
    glist[f'{list(glist.keys())[0]}_unique'] = list((set(glist[list(glist.keys())[0]])-set(glist[list(glist.keys())[1]]))-set(glist[list(glist.keys())[2]]))
    glist[f'{list(glist.keys())[1]}_unique'] = list((set(glist[list(glist.keys())[1]])-set(glist[list(glist.keys())[0]]))-set(glist[list(glist.keys())[2]]))
    glist[f'{list(glist.keys())[0]}:{list(glist.keys())[1]}'] = list(set.intersection(set(glist[list(glist.keys())[0]]),set(glist[list(glist.keys())[1]]))-set(glist['common']))
    glist[f'{list(glist.keys())[2]}_unique'] = list((set(glist[list(glist.keys())[2]])-set(glist[list(glist.keys())[1]]))-set(glist[list(glist.keys())[0]]))
    #----------------------------------------------------------------------------------------------------------------------------------------------------------
    
    glist[f'{list(glist.keys())[0]}:{list(glist.keys())[2]}'] = list(set.intersection(set(glist[list(glist.keys())[0]]),set(glist[list(glist.keys())[2]]))-set(glist['common']))
    glist[f'{list(glist.keys())[2]}:{list(glist.keys())[1]}'] = list(set.intersection(set(glist[list(glist.keys())[1]]),set(glist[list(glist.keys())[2]]))-set(glist['common']))
    del glist['common']
    glist['common'] = list(tmp)
    #----------------------------------------- End of Van-diagram with varification---------------------------------------------------------------------------------------------------
    
    for item in gdata:
      #  print(item)
        tmp = gdata[item].rename(columns={'log2FoldChange':f'{item}.LFC','pvalue':f'{item}.pvalue','padj':f'{item}.padj'})
        gdata[item]=tmp
    j=0
    t=list(glist.keys())[0:3]
    gdata['common']=pd.merge(gdata[t[0]],gdata[t[1]],left_index=True, right_index=True)
    gdata['common']=pd.merge(gdata['common'],gdata[t[2]],left_index=True, right_index=True) 
    
   # print(gdata['common'].head())
    
    for item in glist:
        if ':' in item:
            [x,n] = item.split(':')
            tmp_df = pd.merge(gdata[x],gdata[n],left_index=True, right_index=True)
            gdata[item] = tmp_df[~(tmp_df.index.isin(gdata['common'].index))]
            print(f'they key in {item} : {gdata[item].shape}')
        elif '_unique' in item:
            t=list(glist.keys())[0:3]        
            [x,n] = item.split('_')
            
            
            t.remove(x)
            
            cindex = list(list(gdata[t[0]].index)+list(gdata[t[1]].index)+list(gdata['common'].index))
            
            
            gdata[item] = gdata[x][~gdata[x].index.isin(cindex)]
            print(f'the value of key is {gdata[x].shape} and after removal of index are {gdata[item].shape}')
            
            
            #print(f'they key in {item} : {gdata[item].shape}')
 
    
        
    
    names =list(glist.keys())[-7:]
    venn3_unweighted(subsets = (len(glist[names[0]]), len(glist[names[1]]), len(glist[names[2]]), len(glist[names[3]]), len(glist[names[4]]), len(glist[names[5]]),len(glist[names[6]])), set_labels = list(glist.keys())[0:3], alpha = 0.5);
    plt.title('Venn Diagram of all DE Genes');
   
    
    
    if wflag:
        path = "_".join(list(gdata.keys())[0:3])
        if not os.path.exists('data/Venn_output'):
            os.mkdir('data/Venn_output')
        Write_dict_excel(lsel_dictionary=gdata,path=f'data/Venn_output/{path}.xlsx')
        plt.savefig('data/Venn_output/destination_path.eps', format='eps')
    
    return glist,gdata

def Write_dict_excel(lsel_dictionary,path='data/tmp.xlsx'):
    writer = pd.ExcelWriter(path, engine='xlsxwriter')
    for item in lsel_dictionary:
        tmp = lsel_dictionary[item]
        tmp.to_excel(writer,item.replace(':','_'))
        #workbook.write(tmp)
        print(tmp.shape)
    writer.save()
    
    
glistx,gdatax = Get_VanData()
glistz,gdataz = Van_Prepration(glist=glistx,gdata=gdatax,wflag=True)

