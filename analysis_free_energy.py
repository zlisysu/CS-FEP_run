from threading import Thread
import alchemlyb 
from pymbar import BAR as BAR_
from glob import glob
import os
from os.path import join
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

class optParser():
    def __init__(self, fakeArgs):
        parser = OptionParser()
        parser.add_option('-d', '--dir', dest='file_directory', help='Directory in which datafiles are stored. Default: current directory', default='.')
        parser.add_option('-f', '--fraction', dest='fraction', help='The fraction of the data number that will be used for calculation.', default='0.75')
        parser.add_option('-o', '--output', dest='output_csv_filename', help='The output csv filename of the result. Default: free_ene.csv', default='free_ene_all.csv')
        parser.add_option('-u', '--energy_unit', dest='energy_unit', help='The unit of the output energy. Default: kcal/mol', default='kcal/mol')
        if fakeArgs:
            self.option, self.args = parser.parse_args(fakeArgs)
        else:
            self.option, self.args = parser.parse_args()


class PLOTTING():
    
    def __init__(self,):
        pass
    
    def unpack_esti_dir(self, esti_df):
        frame_ratio = np.around(np.array(esti_df.index), decimals=2)
        fe = np.array(esti_df.iloc[:,0])
        std = np.array(esti_df.iloc[:,1])
        fe_up = fe+std
        fe_down = fe-std
        #print(frame_ratio, fe, fe_up, fe_down)
        return frame_ratio, fe, fe_up, fe_down

    def plot_dU_distribution(self, d_u,png_file=None, ifplt=True, bins_=50):
        plt.clf()
        plt.figure(num=1,figsize=(5.4,3.6))#创建画图，序号为1，图片大小为2.7*1.8
        plt.rcParams['axes.unicode_minus'] = False#使用上标小标小一字号
        # plt.rcParams['font.sans-serif']=['Times New Roman']
        count,bin_edges = np.histogram(d_u,bins=bins_,density=True)
        xu = []
        for i in range(len(bin_edges)-1):
            xu.append((bin_edges[i]+bin_edges[i+1])/2)
        plt.plot(xu,count,"o", color = "green", markersize=4)
        if ifplt:
            if png_file:
                plt.savefig(png_file, format="png",dpi=600, transparent=True)
                plt.show()
            else:
                plt.show()
        return count, xu

    def plot_weighted_dU_distribution(self, d_u, weights_array, png_file=None, ifplt=False, bins_=50):
        plt.clf()
        plt.figure(num=1,figsize=(5.4,3.6))#创建画图，序号为1，图片大小为2.7*1.8
        plt.rcParams['axes.unicode_minus'] = False#使用上标小标小一字号
        # plt.rcParams['font.sans-serif']=['Times New Roman']
        count,bin_edges = np.histogram(d_u,bins=bins_,density=True, weights = weights_array)
        xu = []
        for i in range(len(bin_edges)-1):
            xu.append((bin_edges[i]+bin_edges[i+1])/2)
        plt.plot(xu,count,"o", color = "green", markersize=4)
        if ifplt:
            if png_file:
                plt.savefig(png_file, format="png",dpi=600, transparent=True)
                plt.show()
            else:
                plt.show()
        else:
            plt.clf()
        return count, xu

    def plot_resample_dU_distribution(self, ori_d_u, resample_dU, png_file=None, ifplt=True,bins=100):
        plt.clf()
        pre_y,pre_bin = self.plot_dU_distribution(ori_d_u,png_file=None,ifplt=False,bins_=bins)
        resample_count, xu = self.plot_dU_distribution(resample_dU,png_file=None, ifplt=None,bins_=bins)
        plt.figure(figsize=(8,6))
        ax=plt.gca()
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(0)
        ax.spines['top'].set_linewidth(0)
        plt.plot(pre_bin, pre_y,"o",label='original target dU distrib',color="green",markersize=5)
        plt.plot(xu,resample_count,"x",label='resampled target dU distrib',color="blue",markersize=5)
        plt.legend(loc="best",fontsize=9)
        # ori_d_u_mean, ori_d_u_std = np.mean(ori_d_u), np.std(ori_d_u)
        # resample_d_u_mean, resample_d_u_std = np.mean(resample_dU), np.std(resample_dU)
        # plt.text(union_x[-1],0.5*max(prob_y),'mean_pre_d_u : {:.3f}\nstd_pre_d_u : {:.5f}\nmean_resample_d_u : {:.3f}\nstd_resample_d_u : {:.5f}'.format\(ori_d_u_mean,ori_d_u_std,resample_d_u_mean,resample_d_u_std),\
        #          horizontalalignment='left',verticalalignment='center')#文本注释mean值和std值
        if ifplt:
            if png_file:
                plt.title(f'{png_file}', fontsize=15)#设置图片显示标题
                plt.savefig(png_file, dpi=600, format='png', transparent=True, bbox_inches='tight')
            else:
                plt.show()
        else:
            plt.clf()

    def plot_fe_time_serial(self, png_file_name, **fe_he_std_dir,):
        figsize=8,6
        figure, ax = plt.subplots(figsize=figsize)
        plt.rcParams['axes.unicode_minus'] = False#使用上标小标小一字号
        plt.rcParams['font.sans-serif']=['Times New Roman']#设置全局字体，可选择需要的字体替换掉‘Times New Roman’
        #plt.rcParams['font.sans-serif']=['SimHei']#使用黑体'SimHei'作为全局字体，可以显示中文
        font1={'family': 'Times New Roman', 'weight': 'bold', 'size': 14}#设置字体模板，
        font2={'family': 'Times New Roman', 'weight': 'bold', 'size': 20}#wight为字体的粗细，可选 ‘normal\bold\light’等
        font3={'family': 'Times New Roman', 'weight': 'light', 'size': 12}#size为字体大小
        plt.minorticks_on()#开启小坐标
        plt.tick_params(which='major',width=3, length=6)#设置大刻度的大小
        plt.tick_params(which='minor',width=2, length=4)#设置小刻度的大小
        ax=plt.gca();#获得坐标轴的句柄
        ax.spines['bottom'].set_linewidth(4);#设置底部坐标轴的粗细
        ax.spines['left'].set_linewidth(4);#设置左边坐标轴的粗细
        ax.spines['right'].set_linewidth(4);#设置右边坐标轴的粗细
        ax.spines['top'].set_linewidth(4);#设置上部坐标轴的粗细
        ######
        if len(fe_he_std_dir) == 3:
            dir_move = fe_he_std_dir['moving']
            dir_forw = fe_he_std_dir['forward']
            dir_reve = fe_he_std_dir['reverse']
            frame_ratio_move, fe_move, fe_up_move, fe_down_move = self.unpack_esti_dir(dir_move)
            frame_ratio_forw, fe_forw, fe_up_forw, fe_down_forw = self.unpack_esti_dir(dir_forw)
            frame_ratio_reve, fe_reve, fe_up_reve, fe_down_reve = self.unpack_esti_dir(dir_reve)
            y_min = fe_up_forw[-1]-2
            y_max = fe_up_forw[-1]+2
            ###moving estimate plot
            plt.plot(frame_ratio_move, fe_move, '-', lw=2, color='#75b84f', label='moving estimate', alpha=1)#数据主体
#             plt.plot(frame_ratio_move, fe_up_move, '-', lw=0.001, color='#ffffff')
#             plt.plot(frame_ratio_move, fe_down_move, '-', lw=0.001, color='#ffffff')
            plt.fill_between(frame_ratio_move, fe_down_move, fe_up_move, where=fe_down_move <= fe_up_move,
                     facecolor='#a9f971', interpolate=True,alpha=0.5)
            ###forward estimate plot
            plt.plot(frame_ratio_forw, fe_forw, '-', lw=2, color='#C11B17', label='forward estimate', alpha=1)#数据主体
#             plt.plot(frame_ratio_forw, fe_up_forw, '-', lw=0.001, color='#ffffff')
#             plt.plot(frame_ratio_forw, fe_down_forw, '-', lw=0.001, color='#ffffff')
            plt.fill_between(frame_ratio_forw, fe_down_forw, fe_up_forw, where=fe_down_forw <= fe_up_forw,
                     facecolor='#ff9a8a', interpolate=True, alpha=0.5)
            ###reverse estimate plot
            plt.plot(frame_ratio_reve, fe_reve, '-', lw=2, color='#736AFF', label='reverse estimate', alpha=1)#数据主体
#             plt.plot(frame_ratio_reve, fe_up_reve, '-', lw=0.001, color='#ffffff')
#             plt.plot(frame_ratio_reve, fe_down_reve, '-', lw=0.001, color='#ffffff')
            plt.fill_between(frame_ratio_reve, fe_down_reve, fe_up_reve, where=fe_down_reve <= fe_up_reve,
                     facecolor='#a2cffe', interpolate=True, alpha=0.5)
        elif len(fe_he_std_dir) == 2:
            dir_forw = fe_he_std_dir['forward']
            dir_reve = fe_he_std_dir['reverse']
            frame_ratio_forw, fe_forw, fe_up_forw, fe_down_forw = self.unpack_esti_dir(dir_forw)
            frame_ratio_reve, fe_reve, fe_up_reve, fe_down_reve = self.unpack_esti_dir(dir_reve)
            y_min = fe_up_forw[-1]-2
            y_max = fe_up_forw[-1]+2
            fe_x_plot = np.linspace(0,1,10000)
            estimate_fe = [np.array(fe_reve)[-1] for i in range(0, len(fe_x_plot))]
            estimate_std_range_up = [np.array(fe_up_reve)[-1] for i in range(0, len(fe_x_plot))]
            estimate_std_range_down = [np.array(fe_down_reve)[-1] for i in range(0, len(fe_x_plot))]
            ###estimate_fe_horizontal_line
            plt.plot(fe_x_plot, estimate_fe, '-', lw=2, color='#FFE11A', label='BAR estimate result', alpha=1)
            plt.fill_between(fe_x_plot, estimate_std_range_down, estimate_std_range_up, where=estimate_std_range_down <= estimate_std_range_up,
                     facecolor='#FFE11A', interpolate=True, alpha=0.5)
            ###forward estimate plot
            plt.plot(frame_ratio_forw, fe_forw, '-', lw=2, color='#C11B17', label='forward estimate', alpha=1)#数据主体
            plt.fill_between(frame_ratio_forw, fe_down_forw, fe_up_forw, where=fe_down_forw <= fe_up_forw,
                     facecolor='#ff9a8a', interpolate=True, alpha=0.5)
            ###reverse estimate plot
            plt.plot(frame_ratio_reve, fe_reve, '-', lw=2, color='#736AFF', label='reverse estimate', alpha=1)#数据主体
            plt.fill_between(frame_ratio_reve, fe_down_reve, fe_up_reve, where=fe_down_reve <= fe_up_reve,
                     facecolor='#a2cffe', interpolate=True, alpha=0.5)

            print(f'fe:{np.array(fe_reve)[-1]} kcal/mol')

        elif len(fe_he_std_dir) == 1:
            aly_stra = list(fe_he_std_dir.keys())[0]
            dir_ = fe_he_std_dir[aly_stra]
            frame_ratio_, fe_, fe_up_, fe_down_ = self.unpack_esti_dir(dir_)
            y_min = fe_up_[-1]-2
            y_max = fe_up_[-1]+2
            ### estimate plot
            if aly_stra == 'forward':
                plt.plot(frame_ratio_, fe_, '-', lw=2, color='#C11B17', label='forward estimate', alpha=1)#数据主体
                plt.fill_between(frame_ratio_, fe_down_, fe_up_, where=fe_down_ <= fe_up_,
                        facecolor='#ff9a8a', interpolate=True, alpha=0.5)
            elif aly_stra == 'reverse':
                plt.plot(frame_ratio_, fe_, '-', lw=2, color='#736AFF', label='reverse estimate', alpha=1)#数据主体
                plt.fill_between(frame_ratio_, fe_down_, fe_up_, where=fe_down_ <= fe_up_,
                        facecolor='#a2cffe', interpolate=True, alpha=0.5)
            elif aly_stra == 'moving':
                plt.plot(frame_ratio_, fe_, '-', lw=2, color='#75b84f', label='moving estimate', alpha=1)#数据主体
                plt.fill_between(frame_ratio_, fe_down_, fe_up_, where=fe_down_ <= fe_up_,
                        facecolor='#a9f971', interpolate=True,alpha=0.5)
        ######
        #plt.xlim(-15,-5)#限制x轴的大小
        # print('___________ymin:',y_min,"\n____________",y_max)
        plt.ylim(y_min,y_max)#限制y轴的大小
        #plt.title("Statistical analysis on the restraint strategies used in the test system ",fontdict=font2)#标题
        #plt.xlabel(r'$\Delta G_{exp} $ (kcal/mol)', fontdict=font2)#x轴标签
        #plt.ylabel(r'$\Delta G_{cal}^{MM-GBSA} $ (kcal/mol)',fontdict=font2)#y轴标签
        plt.legend(loc="best",scatterpoints=1,prop=font2,shadow=True,frameon=False)#添加图例，\
        # # loc控制图例位置，“best”为最佳位置，“bottom”,"top"，“topringt"等，\
        # # shadow为图例边框阴影，frameon控制是否有边框
        plt.tick_params(\
            axis='x',#设置x轴
            direction='in',# 小坐标方向，in、out
            which='both',      # 主标尺和小标尺一起显示，major、minor、both
            bottom=True,      #底部标尺打开
            top=False,         #上部标尺关闭
            labelbottom=True, #x轴标签打开
            labelsize=20) #x轴标签大小
        plt.tick_params(\
            axis='y',
            direction='in',
            which='both',
            left=True,
            right=False,
            labelbottom=True,
            labelsize=20)
        plt.yticks(fontproperties='Times New Roman', size=20,weight='bold')#设置x，y坐标轴字体大小及加粗
        plt.xticks(fontproperties='Times New Roman', size=20,weight='bold')#设置x，y坐标轴字体大小及加粗
        plt.ticklabel_format(axis='both',style='sci')#sci文章的风格
        plt.tight_layout(rect=(0,0,1,1))#rect=[left,bottom,right,top]
        plt.tight_layout()#防止因为ticklabel或者title等过长、或者过大导致的显示不全
        plt.savefig(png_file_name,format="png",dpi=600,transparent=True)#存储png
        #plt.show()

def tuple_to_str(tuple_):
    if type(tuple_) == tuple:        
        str_tmp = str(list(tuple_)).replace('[','(',1)
        str_tmp = str_tmp.replace(']',')',1)
        return str_tmp
    else:
        return str(tuple_)

class FEP():
    def __init__(self, d_u):
        self.u_l=np.array(d_u)
        self.u_std=np.std(self.u_l)
        exponential=np.exp(-self.u_l)
        expave=exponential.mean()
        self.ene=-np.log(expave)

class TIME_SERIAL_DATA():
    def __init__(self, time_serial_type, lambda_seq=[], columns=[]):
        self.time_serial_type = time_serial_type
        self.data_dict ={}
        self.lambda_seq = lambda_seq
        if lambda_seq==[]:
            pass
        else:
            print('This is my lambda_seq:{}'.format(str(self.lambda_seq)))
        self.columns_list=columns
        self.columns_list.insert(0, 'time_ratio')
        for i in lambda_seq:
            self.data_dict[i]=[]
        self.df_dict={}
        self.df_dict_with_std = {}
        
    
    def update_data(self, time_ratio, data_):
        '''Updata the free energy data, note that please send me the free energy data with unit of kcal*mol^(-1)
        
        '''
        lambda_seq = data_.index
        for i in lambda_seq:
            data_lst = list(data_.loc[i,:])
            data_lst.insert(0, time_ratio)
            self.data_dict[i].append(data_lst)
    
    def generate_df(self,):
        for i in self.data_dict.keys():
            onetime_data_df = pd.DataFrame(self.data_dict[i])
            onetime_data_df.columns=self.columns_list
            onetime_data_df.index = onetime_data_df.iloc[:,0]
#             print(onetime_data_df)
            self.df_dict[i] = onetime_data_df
    
    def cal_rolling_std(self, std_step_ratio, iffep=True):        
        for key in self.df_dict.keys():
            onetime_data_df = self.df_dict[key]
            std_step = int(np.floor(len(onetime_data_df)*std_step_ratio))
            if iffep:
                to_cal_df = onetime_data_df.loc[:, 'FEP_forward_bythislambda(kcal/mol)':'FEP_reverse_bythislambda(kcal/mol)']
                std_columns = ['FEP_forward_Rolling_STD(kcal/mol)', 'FEP_reverse_Rolling_STD(kcal/mol)']
            else:
                to_cal_df = pd.DataFrame(onetime_data_df.loc[:, 'free_energy(kcal/mol)'])
                std_columns = ['Rolling_STD(kcal/mol)',]
            std_df = to_cal_df.rolling(std_step).std()
#             print(to_cal_df)
#             print(std_df)
            for i in range(std_step):
                if iffep:
                    std_df.iloc[i,0]=std_df.iloc[std_step, 0]
                    std_df.iloc[i,1]=std_df.iloc[std_step, 1]
                else:
                    std_df.iloc[i,0]=std_df.iloc[std_step, 0]
                
            std_df.columns = std_columns
            onetime_data_df = pd.concat([onetime_data_df,std_df], axis = 1)
            self.df_dict_with_std[key] = onetime_data_df
            
    def output_csv(self, csv_prefix, ):
        for key in self.df_dict_with_std.keys():
            onetime_data_df = self.df_dict_with_std[key]
            csv_name = csv_prefix+'_'+str(key)+'_'+self.time_serial_type+'_esti'+'.csv'
            onetime_data_df.to_csv(csv_name)

    def output_half_time_dg(self, iffep=False):
        if iffep:
            print('Current not supported for FEP time series!')
        else:
            all_process_key = self.lambda_seq[-1]
            all_process_time_data_df = self.df_dict_with_std[all_process_key]
            half_line_num = int(np.around(all_process_time_data_df.shape[0]*0.5,0))
            self.half_time_dg = all_process_time_data_df.iloc[half_line_num, 1]
            return self.half_time_dg


class CAL_FREE_ENERGY():
    
    def __init__(self, u_nk_pd, wanted_win_lst=False, scale_f=0.75,  ):
        '''Calculate the free energy according to the multiple lambda internal energy under different force field.
        Parameters
        ----------
        u_nk_pd: pd.DataFrame
            The dataframe obtained by alchemlyb extraction from amber prod.out
        wanted_win_lst: List
            A list assigns simulation windows' data used in the following calculation
        scale_f: float
            To tail the percentage of the every dataframe data used for the calculation
     
        Key properties
        ----------
        self.mbar_lambda_dict: dict
            The keys of it are the sequential integers and its values are the calculated lambda for every simlation
        self.lambda_list: list
            The list used to group the u_nk_pd
        self.simulation_lambda: dict
            The keys of it are the sequential integers and its values are the actual lambda of the simulation
        self.lambda_range: int
            The number of the actual simulations
        self.all_data: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        self.all_data_unk: pd.DataFrame
            The dataframe obtained by concating the list(self.all_data.values())
        '''
        self.u_nk_pd = u_nk_pd
#         self.temperature = temperature
        self.scale_f = scale_f
#         self.ene_unit = ene_unit
        #get mbar index
        self.mbar_lambda_dict = {}
        idx = 0 
        for i in self.u_nk_pd.columns:
            self.mbar_lambda_dict[idx] = i
            idx+=1
        
        ori_index_list = list(self.u_nk_pd.index.names)
        ori_index_list.pop(0)
        self.lambda_list = ori_index_list
        
        a = self.u_nk_pd.groupby(self.lambda_list)
        lambda_dict = {}
        K=0
        if wanted_win_lst is False:
            for i,j in a:
                every_dataframe = pd.DataFrame(j)
                lambda_dict[K]=every_dataframe.tail(n=math.floor(every_dataframe.shape[0]*self.scale_f))
                K+=1
        else:
            for i,j in a:
                if wanted_win_lst[K]==i:
                    every_dataframe = pd.DataFrame(j)
                    lambda_dict[K]=every_dataframe.tail(n=math.floor(every_dataframe.shape[0]*self.scale_f))
                    K+=1
                    if K == len(wanted_win_lst):
                        break
                else:
                    pass
            
        self.all_data = lambda_dict
        self.all_data_unk = pd.concat(list(lambda_dict.values()),sort=True)
        
        lamb_value_dict = {}
        for key in self.all_data.keys():
            index_list = list(self.all_data[key].index[0])
            index_list.pop(0)
            if len(index_list) == 1:
                lamb_value_dict[key] = index_list[0]
            else:
                
                lamb_value_dict[key] = tuple(index_list)
        
        self.simulation_lambda = lamb_value_dict
        self.lambda_range = len(self.all_data)
        

        
    def get_deltaU_in_lambda_by_tuplekey(self, lambda_idx, delta_U_whominus_tuple_key):
        U_key_1 = delta_U_whominus_tuple_key[0]
        U_key_2 = delta_U_whominus_tuple_key[1]
        U_at_lambda_key1 = np.array(self.all_data[lambda_idx][U_key_1])
        U_at_lambda_key2 = np.array(self.all_data[lambda_idx][U_key_2])
        dU_at_lambda = U_at_lambda_key1 - U_at_lambda_key2
        return dU_at_lambda
    
        
    def get_deltaU_in_lambda(self, lambda_idx, delta_U_whominus, filter_=True):
        '''
        Get the speicific d_U in the specified lambda, which was only included in the self.simulation_lambda.
        Parameters
        ----------
        lambda_idx: int, to specify a lambda simulation.
        delta_U_whominus: tuple, int, a two elements tuple that the first element is the lambda index of reduced U.
        filter_: bool, if it equals to True, do the filteration of d_U based on its mean and std.  
        Return
        ----------
        dU_at_lambda: np.array, float, shape=(N,)
                    N is the flames including in each lambda window.
                    Each value in this np.array is the speicific d_U in the specified lambda. 
                    If filter == True, only remain the d_U in the range (d_U_mean-2*std, d_U_mean+2*std).
        '''
        if filter_ == True:
            result_dict = {}
            key_1 = delta_U_whominus[0]
            key_2 = delta_U_whominus[1]
            U_key_1 = self.simulation_lambda[key_1]
            U_key_2 = self.simulation_lambda[key_2]
            U_at_lambda_key1 = np.array(self.all_data[lambda_idx][U_key_1])
            U_at_lambda_key2 = np.array(self.all_data[lambda_idx][U_key_2])
            dU_at_lambda = U_at_lambda_key1 - U_at_lambda_key2
            d_U_mean = dU_at_lambda.mean()
            d_U_std = dU_at_lambda.std()
            up_edge = d_U_mean+2*d_U_std
            down_edge = d_U_mean-2*d_U_std
            bool_index_less = dU_at_lambda<up_edge
            bool_index_more = dU_at_lambda>down_edge
            bool_index_all = bool_index_less*bool_index_more
            dU_at_lambda = dU_at_lambda[bool_index_all]
            result_dict['dU_at_lambda'] = dU_at_lambda
            result_dict['bool_index'] = bool_index_all 
            return result_dict
        else:
            key_1 = delta_U_whominus[0]
            key_2 = delta_U_whominus[1]
            U_key_1 = self.simulation_lambda[key_1]
            U_key_2 = self.simulation_lambda[key_2]
            U_at_lambda_key1 = np.array(self.all_data[lambda_idx][U_key_1])
            U_at_lambda_key2 = np.array(self.all_data[lambda_idx][U_key_2])
            dU_at_lambda = U_at_lambda_key1 - U_at_lambda_key2
            return dU_at_lambda
    
    def cal_FE(self, filename="free_ene.csv",unit='kbT'):
        '''
        ene_unit: str
            To determine the final output csv energy unit, 'kbT' or 'kcal/mol', default is 'kbT'
        '''
        init_data_df = pd.DataFrame()
        last_lambda_key = self.mbar_lambda_dict[len(self.mbar_lambda_dict)-1]
        first_lambda_key = self.mbar_lambda_dict[0]
        for i in range(self.lambda_range):
            if i == self.lambda_range-1:
                if last_lambda_key in self.simulation_lambda.values():
                    single_df = pd.DataFrame()
                else:
                    esti_ene, esti_std = self.cal_FE_last_window_bar()
                    index_info=tuple_to_str(self.simulation_lambda[i])+' to '+ tuple_to_str(last_lambda_key)
                    single_df = pd.DataFrame(columns=['delta_A_what_to_what', 'free_energy(kbT)', 'estimated std'], index =[0,] )
                    single_df['delta_A_what_to_what'] = index_info
                    single_df['free_energy(kbT)'] = esti_ene
                    single_df[ 'estimated std'] = esti_std
            elif i == 0:
                if first_lambda_key in self.simulation_lambda.values():
                    esti_ene, esti_std= self.cal_FE_middle_window_bar(i, )
                    index_info=tuple_to_str(self.simulation_lambda[i])+' to '+ tuple_to_str(self.simulation_lambda[i+1])
                    single_df = pd.DataFrame(columns=['delta_A_what_to_what', 'free_energy(kbT)', 'estimated std'], index =[0,] )
                    single_df['delta_A_what_to_what'] = index_info
                    single_df['free_energy(kbT)'] = esti_ene
                    single_df[ 'estimated std'] = esti_std                    
                else:
                    esti_ene_0, esti_std_0 = self.cal_FE_first_window_bar()
                    index_info_0=tuple_to_str(first_lambda_key)+' to '+ tuple_to_str(self.simulation_lambda[i])
                    esti_ene_1, esti_std_1 = self.cal_FE_middle_window_bar(i, )
                    index_info_1=tuple_to_str(self.simulation_lambda[i])+' to '+ tuple_to_str(self.simulation_lambda[i+1])
                    single_df = pd.DataFrame(columns=['delta_A_what_to_what', 'free_energy(kbT)', 'estimated std'], index =[0,1] )
                    single_df.iloc[0,0]=index_info_0
                    single_df.iloc[0,1]=esti_ene_0
                    single_df.iloc[0,2]=esti_std_0
                    single_df.iloc[1,0]=index_info_1
                    single_df.iloc[1,1]=esti_ene_1
                    single_df.iloc[1,2]=esti_std_1
            else:
                esti_ene, esti_std= self.cal_FE_middle_window_bar(i, )
                index_info=tuple_to_str(self.simulation_lambda[i])+' to '+ tuple_to_str(self.simulation_lambda[i+1])
                single_df = pd.DataFrame(columns=['delta_A_what_to_what', 'free_energy(kbT)', 'estimated std'], index =[0,] )
                single_df['delta_A_what_to_what'] = index_info
                single_df['free_energy(kbT)'] = esti_ene
                single_df[ 'estimated std'] = esti_std
            init_data_df = pd.concat([init_data_df, single_df], axis = 0 )
        sum_df = pd.DataFrame(columns=['delta_A_what_to_what', 'free_energy(kbT)', 'estimated std'], index =[0,] )
        sum_df['delta_A_what_to_what']=tuple_to_str(first_lambda_key)+' to '+tuple_to_str(last_lambda_key)
        sum_df['free_energy(kbT)']=init_data_df.iloc[:,1].sum()
        sum_df['estimated std']=((init_data_df.iloc[:,2]**2).sum())**0.5
        init_data_df=pd.concat([init_data_df,sum_df],axis=0)
        init_data_df.index = init_data_df['delta_A_what_to_what']
        init_data_df=init_data_df.drop(columns=['delta_A_what_to_what',])
        #init_data_df['free_energy(kbT)'].astype(float)
        #init_data_df['estimated std'].astype(float)
        if unit=='kbT':
            pass
            #print('aaaa')
#             init_data_df.to_csv(filename)
        elif unit=='kcal/mol':
            #print('bbbb')
            init_data_df.columns=['free_energy(kcal/mol)', 'estimated std']
            init_data_df['free_energy(kcal/mol)']=init_data_df['free_energy(kcal/mol)'].astype(float)
            init_data_df['estimated std']=init_data_df['estimated std'].astype(float)
            #print(init_data_df.dtypes)
            init_data_df[init_data_df.select_dtypes(include=['float64']).columns] *= 0.5922
            #print(init_data_df)
        if filename:
            init_data_df.to_csv(filename)
        
        return init_data_df

    def cal_FE_FEP(self, filename='free_ene_fep.csv', unit='kbT'):
        '''calculate the free energy difference between the simulaion lambda by FEP
        Parameters
        ----------
        filename: str
            The output csv filename of the result dataframe, default is 'free_ene_fep.csv'
        ene_unit: str
            To determine the final output csv energy unit, 'kbT' or 'kcal/mol', default is 'kbT'
        
        Return    
        ----------
        init_data_df: pd.DataFrame
            The final output of calculated free energy
        '''
        init_data_df = pd.DataFrame()
        last_lambda_key = self.mbar_lambda_dict[len(self.mbar_lambda_dict)-1]
        first_lambda_key = self.mbar_lambda_dict[0]
        for i in range(self.lambda_range):
            lambda_value = self.simulation_lambda[i]
            if i == 0:
                if first_lambda_key in self.simulation_lambda.values():
                    fe_forward = self.cal_FE_FEP_forward(i)
                    fe_reverse = pd.NA
                else:
                    fe_forward = self.cal_FE_FEP_forward(i)
                    fe_reverse = -self.cal_FE_first_window_bar()[0]
            elif i == self.lambda_range-1:
                if last_lambda_key in self.simulation_lambda.values():
                    fe_forward = pd.NA
                    fe_reverse = -self.cal_FE_FEP_reverse(i)
                else:
                    fe_forward = self.cal_FE_last_window_bar()[0]
                    fe_reverse = -self.cal_FE_FEP_reverse(i)
            else:
                fe_forward = self.cal_FE_FEP_forward(i)
                fe_reverse = -self.cal_FE_FEP_reverse(i)
            single_df = pd.DataFrame(columns=['lambda_value', 'FEP_forward_bythislambda(kbT)', 'FEP_reverse_bythislambda(kbT)'], index=[0,])
            single_df['lambda_value'] = lambda_value
            single_df['FEP_forward_bythislambda(kbT)'] = fe_forward
            single_df['FEP_reverse_bythislambda(kbT)'] = fe_reverse
            init_data_df = pd.concat([init_data_df,single_df], axis = 0)
        sum_df = pd.DataFrame(columns=['lambda_value', 'FEP_forward_bythislambda(kbT)', 'FEP_reverse_bythislambda(kbT)'], index=[0,])
        sum_df['lambda_value'] = 'sum_of_all'
        sum_df['FEP_forward_bythislambda(kbT)'] = init_data_df.iloc[:,1].sum()
        sum_df['FEP_reverse_bythislambda(kbT)'] = init_data_df.iloc[:,2].sum()
        
        init_data_df=pd.concat([init_data_df,sum_df],axis=0)
        init_data_df.index = init_data_df['lambda_value']
        init_data_df=init_data_df.drop(columns=['lambda_value',])
        
        if unit=='kbT':
            pass
        elif unit=='kcal/mol':
            init_data_df.columns=['FEP_forward_bythislambda(kcal/mol)','FEP_reverse_bythislambda(kcal/mol)']
            init_data_df['FEP_forward_bythislambda(kcal/mol)']=init_data_df['FEP_forward_bythislambda(kcal/mol)'].astype(float)
            init_data_df['FEP_reverse_bythislambda(kcal/mol)']=init_data_df['FEP_reverse_bythislambda(kcal/mol)'].astype(float)
            init_data_df[init_data_df.select_dtypes(include=['float64']).columns] *= 0.5922
        if filename:
            init_data_df.to_csv(filename)
            
        return init_data_df
        
    def cal_FE_FEP_forward(self, lambda_idx):
        i=lambda_idx
        ori_du_f = self.get_deltaU_in_lambda(i,(i+1,i),False)
        obj_f = FEP(ori_du_f)
        fe_fep = obj_f.ene
        return fe_fep
    
    def cal_FE_FEP_reverse(self, lambda_idx):
        i=lambda_idx
        ori_du_b = self.get_deltaU_in_lambda(i, (i-1,i),False)
        obj_b = FEP(ori_du_b)
        fe_fep = obj_b.ene
        return fe_fep

    def cal_FE_first_window_bar(self, ):
        i=0
        first_key = self.simulation_lambda[0]
        second_key = self.mbar_lambda_dict[0]
        delta_U_tuplekey = (first_key, second_key)
        ori_du_b = self.get_deltaU_in_lambda_by_tuplekey(0,delta_U_tuplekey)
        #fep_exp
        fep_obj_b = FEP(ori_du_b)
        fep_exp_fe_b = fep_obj_b.ene
        return fep_exp_fe_b, 0.0

    def cal_FE_last_window_bar(self,):
        i=self.lambda_range-1
        first_key = self.mbar_lambda_dict[len(self.mbar_lambda_dict)-1]
        second_key = self.simulation_lambda[i]
        delta_U_tuplekey = (first_key, second_key)
        ori_du_f = self.get_deltaU_in_lambda_by_tuplekey(i,delta_U_tuplekey)
        #fep_exp
        fep_obj_f = FEP(ori_du_f)
        fep_exp_fe_f = fep_obj_f.ene
        return fep_exp_fe_f, 0.0
    
    def cal_FE_middle_window_bar(self,lambda_idx):
        i=lambda_idx
        ori_du_f = self.get_deltaU_in_lambda(i,(i+1,i),False)
#         obj_f=FEP(ori_du_f)
#         std_f=obj_f.u_std
        ori_du_b = self.get_deltaU_in_lambda(i+1,(i,i+1),False)
#         obj_b = FEP(ori_du_b)
#         std_b=obj_b.u_std
        fe_bar, dfe_bar = BAR_(ori_du_f, ori_du_b, method='self-consistent-iteration', maximum_iterations=1000, verbose=False)
        return fe_bar, dfe_bar
        
    
    def plot_overlap_matrix(self, png_file=None):
        from alchemlyb.estimators import MBAR
        '''
        Plot the overlap matrix using MBAR weight values.
        '''
        if png_file == None: 
            mbar__ = MBAR()
            mbar__.fit(self.all_data_unk)
            #overlap_matrix
            overlap_matx=mbar__.overlap_matrix
        else:
            mbar__ = MBAR()
            mbar__.fit(self.all_data_unk)
            #overlap_matrix
            overlap_matx=mbar__.overlap_matrix
            ax1=alchemlyb.visualisation.plot_mbar_overlap_matrix(matrix=overlap_matx)
            ax1.figure.savefig(png_file, dpi=600, format='png', transparent=True)
            plt.show()
            
        return overlap_matx

    def DC_MBAR_overlap_values(self, lambda_0, lambda_1,):
        '''
        Use the relationship 
        \int \frac{\rho_{0} \times \rho_{1}}{\rho_{1}+\rho_{0}} d q^{N}=overlap_0=overlap_1 to calculate the overlap_0 and overlap_1 for checking the overlap degree between two windows,
        which: 
        overlap_0=\left\langle\frac{1}{1+\exp \left(\Delta U_{1-0}-\Delta G_{1-0}\right)}\right\rangle_{0}
        overlap_1=\left\langle\frac{1}{1+\exp \left(-\Delta U_{1-0}+\Delta G_{1-0}\right)}\right\rangle_{1}
        Both of two overlap values should be between 0 and 0.5. 
        The closer the value is to 0.5, the better the overlap of the corresponding two distributions. 
        The closer the value is to 0, the two distributions corresponding to the surface have almost no overlap.
        Note that two lambda windows sample point must be same!
        
        Parameters
        ----------
        lambda_0: float(in amber rbfe), the key in the simulation_lambda(dict) to assign the first window
        lambda_1: float(in amber rbfe), the key in the simulation_lambda(dict) to assign the second window
        
        Generated key variables
        ----------
        delta_U_1to0_in_lambda_0: np.array, float, shape=(N,). \Delta U_{back-front} in front lambda window.
        delta_U_1to0_in_lambda_1: np.array, float, shape=(N,). \Delta U_{back-front} in back lambda window.
            N is the number of the flames in each lambda window.
            
        Return
        ----------
        overlap_0: float, the overlap value calculated from lambda 0
        overlap_1: float, the overlap value calculated from lambda 1
        dff:       float, the bar std calculated bwteen lambda 0 and lambda 1
        '''
        U_key_0 = self.simulation_lambda[lambda_0]
        U_key_1 = self.simulation_lambda[lambda_1]
        U_0_in_lambda_0 = np.array(self.all_data[lambda_0][U_key_0])
        U_1_in_lambda_0 = np.array(self.all_data[lambda_0][U_key_1])
        U_0_in_lambda_1 = np.array(self.all_data[lambda_1][U_key_0])
        U_1_in_lambda_1 = np.array(self.all_data[lambda_1][U_key_1])
        delta_U_1to0_in_lambda_0 = U_1_in_lambda_0 - U_0_in_lambda_0
        delta_U_1to0_in_lambda_1 = U_1_in_lambda_1 - U_0_in_lambda_1
        delta_U_0to1_in_lambda_1 = -delta_U_1to0_in_lambda_1
        df, dff = BAR_(delta_U_1to0_in_lambda_0, delta_U_0to1_in_lambda_1, method='self-consistent-iteration', maximum_iterations=1000, verbose=False)
        
        delta_U01_minus_df_0 = delta_U_1to0_in_lambda_0-df
        df_minus_delta_U01_1 = df-delta_U_1to0_in_lambda_1
        fermi_0 = 1/(1+np.exp(delta_U01_minus_df_0))
        fermi_1 = 1/(1+np.exp(df_minus_delta_U01_1))
        overlap_0 = np.mean(fermi_0)
        overlap_1 = np.mean(fermi_1)
        return overlap_0, overlap_1, dff
    
    def check_DC_MBAR_overlap_overall(self, ):
        '''
        Calculate the overlap value between two adjacent windows through all the windows.
        '''
        for i in range(0, self.lambda_range-1):            
            lambda_0 = i
            lambda_1 = i+1
            overlap_value = self.DC_MBAR_overlap_values(lambda_0, lambda_1)
            # print('This window lambda value is {}, its index is {}'.format(self.simulation_lambda[lambda_0], lambda_0))
            print('The overlap value between {} and {} is {}'.format(self.simulation_lambda[lambda_0], self.simulation_lambda[lambda_1], overlap_value))
            print('') 

class ANA_FEP_TOOLS():
    
    def __init__(self, all_df, ):
        '''The Class used to analyze the one tail process's free energy convergence.
        Parameters
        ----------
        all_df: pd.DataFrame
            The dataframe obtained by alchemlyb extraction from amber prod.out
               
        Key properties
        ----------
        self.forward_esti_obj: TIME_SERIAL_DATA obj
            The TIME_SERIAL_DATA object that contains the three demension data(fe&std, lambda, time_ratio) of the forward estimate
        self.reverse_esti_obj: TIME_SERIAL_DATA obj
            The TIME_SERIAL_DATA object that contains the three demension data(fe&std, lambda, time_ratio) of the reverse estimate
        self.moving_esti_obj: TIME_SERIAL_DATA obj
            The TIME_SERIAL_DATA object that contains the three demension data(fe&std, lambda, time_ratio) of the moving estimate
        '''
        self.all_df = all_df
        ori_index_list = list(self.all_df.index.names)
        ori_index_list.pop(0)
        self.lambda_list = ori_index_list
        self.forward_esti_obj = TIME_SERIAL_DATA('forward')
        self.reverse_esti_obj = TIME_SERIAL_DATA('reverse')
        self.moving_esti_obj = TIME_SERIAL_DATA('moving')
  

    def generate_data_dict(self, wanted_win_lst, scale_f):
        '''Apply the specific window list and scaling factor to get the data dict that used for time_serial analysis or free energy calculation
        Parameters
        ----------
        wanted_win_lst: list
            A list assigns simulation windows' data used in the following calculation
        scale_f: float
            To tail the percentage of the every dataframe data used for the calculation
        
        Return
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        '''
        ###apply scaling###
        tail_df_list = [j.tail(n=math.floor(j.shape[0]*scale_f)) for i,j in self.all_df.groupby(self.lambda_list)]
        tail_df = pd.concat(tail_df_list,axis=0)
        groupobj = tail_df.groupby(self.lambda_list)
        ###
        ###apply partial lambda###
        lambda_dict = {}
        K=0
        if wanted_win_lst is False:
            for i,j in groupobj:
                every_dataframe = pd.DataFrame(j)
                lambda_dict[K] = every_dataframe
                K+=1
        else:
            for i,j in groupobj:
                if wanted_win_lst[K]==i:
                    every_dataframe = pd.DataFrame(j)
                    lambda_dict[K] = every_dataframe
                    K+=1
                    if K == len(wanted_win_lst):
                        break
                else:
                    pass
        ###
        lambda_dict = lambda_dict
        return lambda_dict

    
    def plot_time_serial(self, png_prefix, lambda_, get_which_two_columns, plot_plan, ):
        '''Using the PLOTTING object to plot the time serial analysis for specifc lambda_ (lambda_ may be 'sum_of_all' or '0.0 to 1.0')
        Parameters
        ----------
        png_prefix: str
            The string that given for the output png's name, usually contain the information about the target, run-1 or run-2, which pair(38-60), which tail(cM2A)
        lambda_: str or float
            To specify which window for analysis
            For EXP(FEP) estimator, lambda_ can be one of [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 'sum_of_all']
            For BAR estimator, lambda_ can be one of ['0.0 to 0.05', '0.05 to 0.1', '0.1 to 0.2', '0.2 to 0.3', '0.3 to 0.4', '0.4 to 0.5', '0.5 to 0.6', 
            '0.6 to 0.7', '0.7 to 0.8', '0.8 to 0.9', '0.9 to 0.95', '0.95 to 1.0', '0.0 to 1.0']
        get_which_two_columns: [str, str]
            To specify which free energy and std columns for analysis
            For EXP(FEP) estimator, get_which_two_columns can be ['FEP_forward_bythislambda(kcal/mol)', 'FEP_forward_Rolling_STD(kcal/mol)'] 
            or ['FEP_reverse_bythislambda(kcal/mol)', 'FEP_reverse_Rolling_STD(kcal/mol)']
            For BAR estimator, get_which_two_columns will be ['free_energy(kcal/mol)', 'Rolling_STD(kcal/mol)']
        '''
        plot_obj = PLOTTING()
        fe_columns_name = get_which_two_columns[0]
        std_columns_name = get_which_two_columns[1]
        if plot_plan == 2:
            for_df = self.forward_esti_obj.df_dict_with_std[lambda_]
            rev_df = self.reverse_esti_obj.df_dict_with_std[lambda_]
            forward_two_col_df = pd.concat([for_df[fe_columns_name], for_df[std_columns_name]],axis=1)
            reverse_two_col_df = pd.concat([rev_df[fe_columns_name], rev_df[std_columns_name]],axis=1)
            fe_he_std_dir = {'forward':forward_two_col_df, 
                             'reverse':reverse_two_col_df,}
            plot_obj.plot_fe_time_serial('{}_{}_2.png'.format(png_prefix, lambda_), **fe_he_std_dir)
        elif plot_plan == 3:
            mov_df = self.moving_esti_obj.df_dict_with_std[lambda_]
            for_df = self.forward_esti_obj.df_dict_with_std[lambda_]
            rev_df = self.reverse_esti_obj.df_dict_with_std[lambda_]
            forward_two_col_df = pd.concat([for_df[fe_columns_name], for_df[std_columns_name]],axis=1)
            reverse_two_col_df = pd.concat([rev_df[fe_columns_name], rev_df[std_columns_name]],axis=1)
            moving_two_col_df = pd.concat([mov_df[fe_columns_name], mov_df[std_columns_name]], axis=1)
            fe_he_std_dir = {'forward':forward_two_col_df, 
                             'reverse':reverse_two_col_df,
                             'moving':moving_two_col_df}
            plot_obj.plot_fe_time_serial('{}_{}_3.png'.format(png_prefix, lambda_), **fe_he_std_dir)
        elif plot_plan == 1:
            mov_df = self.moving_esti_obj.df_dict_with_std[lambda_]
            moving_two_col_df = pd.concat([mov_df[fe_columns_name], mov_df[std_columns_name]], axis=1)
            fe_he_std_dir = {'moving':moving_two_col_df, }
            plot_obj.plot_fe_time_serial('{}_{}_1.png'.format(png_prefix, lambda_), **fe_he_std_dir)
        else:
            print('Error! No such plot plan!')
        
    
    def lambda_dict_to_df(self, lambda_dict):
        '''Convert the lambda_dict to lambda_df, which is the key input for the CAL_FREE_ENERGY object.
        Parameters
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        
        Return 
        ----------
        processed_df: pd.DataFrame
            The dataframe generated by concating every windows data 
        '''
        
        df_list = []
        for single_df in lambda_dict.values():
            df_list.append(single_df)
        processed_df = pd.concat(df_list, sort=True)
        return processed_df
    
    def single_interval_estimate(self, lambda_dict, start_frame, end_frame, iffep=True, ):
        '''Calculate free energy according to the given lambda_dict, time_interval of [start_frame, end_frame]
        Parameters
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        start_frame: int
            The integer (actually the index) to represent the start frame of the trajactory
        end_frame: int
            The integer (actually the index) to represent the end frame of the trajactory
        iffep: bool
            To determine if use EXP(FEP) estimator to calculate, default: True
            
        Return
        ----------
        result_fe: pd.DataFrame
            For iffep==True: result_fe.index may be like [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 'sum_of_all'] with index name of 'lambda_value' 
                             result_fe.columns be like ['FEP_forward_bythislambda(kcal/mol)', 'FEP_reverse_bythislambda(kcal/mol)']
            For iffep==False:result_fe.index may be like ['0.0 to 0.05', '0.05 to 0.1', '0.1 to 0.2', '0.2 to 0.3', '0.3 to 0.4', '0.4 to 0.5', '0.5 to 0.6', 
            '0.6 to 0.7', '0.7 to 0.8', '0.8 to 0.9', '0.9 to 0.95', '0.95 to 1.0', '0.0 to 1.0'] with index name of 'delta_A_what_to_what'
                             result_fe.columns be like ['free_energy(kcal/mol)', 'estimated std']
        '''
        df_list = []
        for single_df in lambda_dict.values():
            single_df = single_df
            new_single_df = single_df.iloc[start_frame: end_frame,:]
            df_list.append(new_single_df)
        processed_df = pd.concat(df_list, sort=True)
#         return processed_df
#         processed_df = self.lambda_dict_to_df(lambda_dict)

        fep_obj = CAL_FREE_ENERGY(processed_df, False, 1,)
        if iffep:
            result_fe = fep_obj.cal_FE_FEP(None, 'kcal/mol')
        else:
            result_fe = fep_obj.cal_FE(None, 'kcal/mol')
        return result_fe
        
   
    def moving_estimate(self, lambda_dict, std_step_ratio, divided_ratio, width_ratio = 0.2, iffep=False):
        '''Check the time serial convergence by moving estimate using the time interval with fixed width to scan the time serial data forward. 
        Parameters
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        std_step_ratio: float
            To assign the width that is used in the standard deviation of every single time_ratio point, which uses the pd.df.rolling(std_step_ratio).std()
            to calculate std (e.g. std_step_ratio: 0.1, The data forward 10% of each point will be used to calculate std)
        divided_ratio: float
            To determine the number of points in the time series will be used to show the convergence of the time series
            (e.g. divided_ratio: 0.01, Will generate 100 points like 0.01, 0.02, 0.03, 0.04, ..., 1.00)
        width_ratio: float
            To determine the width of the fixed width of the moving time interval
        iffep: bool
            To determine if use EXP(FEP) estimator to calculate, default: False
        
        Update
        ----------
        self.moving_esti_obj: TIME_SERIAL_DATA obj
            The TIME_SERIAL_DATA object that contains the three demension data(fe&std, lambda, time_ratio) of the moving estimate
        '''
        ###get frame number###
        count = 0
        for j in lambda_dict.values():
            count+=1
            if count == 1:
                frame_num = len(j)
                break
        ######
        label_width = int(np.floor(frame_num*divided_ratio))
        width = int(np.floor(frame_num*width_ratio))
        end_point_list = [i for i in range(0, frame_num, label_width)]
        end_point_list.pop(0)
        end_point_list.append(frame_num)
        count_=0
        for end_frame in end_point_list:
            start_frame = end_frame-width
            if start_frame < 0:
                start_frame = 0
            time_ratio = end_frame/frame_num
            data_ = self.single_interval_estimate(lambda_dict, start_frame,end_frame, iffep)
            if count_ == 0:
                lambda_seq = list(data_.index)
                columns = list(data_.columns)
                time_serial_data_obj = TIME_SERIAL_DATA('moving', lambda_seq, columns)
            else:
                pass
            time_serial_data_obj.update_data(time_ratio, data_)
            count_+=1
        time_serial_data_obj.generate_df()
        time_serial_data_obj.cal_rolling_std(std_step_ratio, iffep)        
        
        self.moving_esti_obj = time_serial_data_obj
        

    def forward_estimate(self, lambda_dict, std_step_ratio, divided_ratio=0.1, iffep=False):
        '''Check the time serial convergence by accumulating data points from the first frame forward with a fixed length
        Parameters
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        std_step_ratio: float
            To assign the width that is used in the standard deviation of every single time_ratio point, which uses the pd.df.rolling(std_step_ratio).std()
            to calculate std (e.g. std_step_ratio: 0.1, The data forward 10% of each point will be used to calculate std)
        divided_ratio: float
            To determine the number of points in the time series will be used to show the convergence of the time series
            (e.g. divided_ratio: 0.01, Will generate 100 points like 0.01, 0.02, 0.03, 0.04, ..., 1.00)
        iffep: bool
            To determine if use EXP(FEP) estimator to calculate, default: False
        
        Update
        ----------
        self.forward_esti_obj: TIME_SERIAL_DATA obj
            The TIME_SERIAL_DATA object that contains the three demension data(fe&std, lambda, time_ratio) of the forward estimate
        '''
        ###get frame number###
        count = 0
        for j in lambda_dict.values():
            count+=1
            if count == 1:
                frame_num = len(j)
                break
        ######
        end_point_inter_width = int(np.floor(frame_num*divided_ratio))
        end_point_list = [i for i in range(0, frame_num, end_point_inter_width)]
        end_point_list.pop(0)
        end_point_list.append(frame_num)
        start_frame=0
        count_=0
        for end_frame in end_point_list:
            time_ratio = end_frame/frame_num
            data_ = self.single_interval_estimate(lambda_dict, start_frame,end_frame, iffep)
            
            if count_ == 0:
                lambda_seq = list(data_.index)
                columns = list(data_.columns)
                time_serial_data_obj = TIME_SERIAL_DATA('forward', lambda_seq, columns)
            else:
                pass
            time_serial_data_obj.update_data(time_ratio, data_)
            count_+=1
        time_serial_data_obj.generate_df()
        time_serial_data_obj.cal_rolling_std(std_step_ratio, iffep)        
        
        self.forward_esti_obj = time_serial_data_obj
    
    def reverse_estimate(self, lambda_dict, std_step_ratio, divided_ratio=0.1, iffep=False):
        '''Check the time serial convergence by accumulating data points from the last frame backward with a fixed length
        Parameters
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        std_step_ratio: float
            To assign the width that is used in the standard deviation of every single time_ratio point, which uses the pd.df.rolling(std_step_ratio).std()
            to calculate std (e.g. std_step_ratio: 0.1, The data forward 10% of each point will be used to calculate std)
        divided_ratio: float
            To determine the number of points in the time series will be used to show the convergence of the time series
            (e.g. divided_ratio: 0.01, Will generate 100 points like 0.01, 0.02, 0.03, 0.04, ..., 1.00)
        iffep: bool
            To determine if use EXP(FEP) estimator to calculate, default: False
            
        Update
        ----------
        self.reverse_esti_obj: TIME_SERIAL_DATA obj
            The TIME_SERIAL_DATA object that contains the three demension data(fe&std, lambda, time_ratio) of the reverse estimate
        '''
        ###get frame number###
        count = 0
        for j in lambda_dict.values():
            count+=1
            if count == 1:
                frame_num = len(j)
                break
        ######
        start_point_inter_width = int(np.floor(frame_num*divided_ratio))
        start_point_list = [i for i in range(0, frame_num, start_point_inter_width)]
        start_point_list = start_point_list[::-1]
        end_frame =frame_num        
        count_ = 0
        for start_frame in start_point_list:
            time_ratio = 1 - start_frame/frame_num
            data_ = self.single_interval_estimate(lambda_dict, start_frame,end_frame, iffep)
            lambda_seq = list(data_.index)
            columns = list(data_.columns)
            if count_ == 0:
                lambda_seq = list(data_.index)
                columns = list(data_.columns)
                time_serial_data_obj = TIME_SERIAL_DATA('reverse', lambda_seq, columns)
            else:
                pass
            time_serial_data_obj.update_data(time_ratio, data_)
            count_+=1
        time_serial_data_obj.generate_df()
        time_serial_data_obj.cal_rolling_std(std_step_ratio, iffep)        
        
        self.reverse_esti_obj = time_serial_data_obj

        
    def use_fep_check_time_serial(self, plot_plan, png_prefix, use_forward=True,):
        '''Use fep calculated free energy value to check time serial convergence, which will output the figure and update self.forward_esti_obj, self.reverse_esti_obj, self.moving_esti_obj
        '''
        lambda_dict = self.generate_data_dict(False, 1)
        self.moving_estimate(lambda_dict, 0.1, 0.01, 0.2, True)
        self.forward_estimate(lambda_dict, 0.1, 0.01, True)
        self.reverse_estimate(lambda_dict, 0.1, 0.01, True)
        lambda_lst = self.forward_esti_obj.df_dict_with_std.keys()
        if use_forward:
            get_which_two_columns=['FEP_forward_bythislambda(kcal/mol)', 'FEP_forward_Rolling_STD(kcal/mol)']
            png_prefix = png_prefix+'_useForward'
            for lambda_ in lambda_lst:
                if lambda_ == 1.0:
                    pass
                else:
                    self.plot_time_serial(png_prefix, lambda_, get_which_two_columns, plot_plan, )
        else:
            get_which_two_columns=['FEP_reverse_bythislambda(kcal/mol)', 'FEP_reverse_Rolling_STD(kcal/mol)']
            png_prefix = png_prefix+'_useReverse'
            for lambda_ in lambda_lst:
                if lambda_ == 0.0:
                    pass
                else:
                    self.plot_time_serial(png_prefix, lambda_, get_which_two_columns, plot_plan, )
    
    
    def use_bar_check_time_serial(self, fraction=1):
        '''Use fep calculated free energy value to check time serial convergence, which will update self.forward_esti_obj, self.reverse_esti_obj, self.moving_esti_obj
        '''
        lambda_dict = self.generate_data_dict(False, fraction)
        self.moving_estimate(lambda_dict, 0.1, 0.01, 0.2, False)
        self.forward_estimate(lambda_dict, 0.1, 0.01, False)
        self.reverse_estimate(lambda_dict, 0.1, 0.01, False)
        
                            
    def output_fe_and_std(self, lambda_dict, scale_f, fe_mode, std_mode,):
        '''Output the calculated free energy and standard deviation according to the given data dict
        Parameters
        ----------
        lambda_dict: dict
            The keys of it are the sequential integers and its values are the every dataframe of the single lambda simulation
        scale_f: float
            Need to assign as the same scale_f value used in generate lambda_dict, when the std_mode is 'time_serial' 
        fe_mode: str
            'BAR' or 'FEP';
            'BAR': Use Bennett Acceptance Ratio estimator to calculate free energy
            'FEP': Use EXPonential averaging (EXP) estimator to calculate free energy
        std_mode: str
            'time_serial' or 'bar_std'
            'time_serial': Output the standard deviation calculated by the last 100*scale_f% moving estimated free energies' std
            'bar_std': Output the standard deviation calculated by BAR estimator. Note: the bar_std may underesitmate the std
         '''
        lambda_unk_df = pd.concat(list(lambda_dict.values()),sort=True)
        cal_fe_obj = CAL_FREE_ENERGY(lambda_unk_df, False, 1)
#         print(cal_fe_obj.cal_FE(None, 'kcal/mol'))
#         print(cal_fe_obj.cal_FE_FEP(None, 'kcal/mol'))
        if fe_mode == 'BAR':
            data_df = cal_fe_obj.cal_FE(None, 'kcal/mol')
            fe = data_df.iloc[:,0][-1]
#             print(data_df)
            if std_mode == 'time_serial':
                self.moving_estimate(lambda_dict, 0.1, 0.01, 0.2, False)
                std = self.moving_esti_obj.df_dict_with_std['0.0 to 1.0'].loc[1-scale_f:1, 'free_energy(kcal/mol)'].std()
            elif std_mode == 'bar_std':
                std = data_df.iloc[:,1][-1]
        elif fe_mode == 'FEP': 
            print('Note that the FEP so far calculation mode may not calculate the whole process free energy difference, only for convergence test!')
            fe = cal_fe_obj.cal_FE_FEP(None, 'kcal/mol').iloc[:,0][-1]
            self.moving_estimate(lambda_dict, 0.1, 0.01, 0.2, True)
#             print(self.moving_esti_obj.df_dict_with_std['sum_of_all'].loc[:scale_f, 'FEP_forward_bythislambda(kcal/mol)'])
            std = self.moving_esti_obj.df_dict_with_std['sum_of_all'].loc[1-scale_f:1, 'FEP_forward_bythislambda(kcal/mol)'].std()
        return fe, std            

    
    def check_part_lambda(self, part_lambda_list, scale_f, fe_mode, std_mode, filename):
        ###generate lambda_dict
        lambda_dict_list = []
        for part_lambda in part_lambda_list:
            lambda_dict_list.append(self.generate_data_dict(part_lambda, scale_f))
        fe_data_df_list = []
        for lambda_dict_ in lambda_dict_list:
            fe, std = self.output_fe_and_std(lambda_dict_, scale_f, fe_mode, std_mode,)
            fe_data_df_list.append([fe, std])
        fe_data_df = pd.DataFrame(fe_data_df_list, )
        fe_data_df.columns = ['fe', 'std']
        fe_data_df.index = [str(i) for i in part_lambda_list]
        fe_data_df.index.names = ['lambda_scheme']
        if filename:
            fe_data_df.to_csv(filename)
            
        return fe_data_df
            
    
#     def overlap_matrix_each_pair(self, lambda_dict, scale_f):
        
#         for i in range(self.all_df_fe_obj.lambda_range-1):
#             win_lst = [self.all_df_fe_obj.simulation_lambda[i], self.all_df_fe_obj.simulation_lambda[i+1]]
#             self.overlap_matrix_for_specific_pair(scale_f, win_lst)
    
#     def overlap_matrix_for_specific_pair(self, scale_f, win_lst):    
#         print(win_lst)
#         adjacent_fe_obj = CAL_FREE_ENERGY(self.all_df, win_lst, scale_f)
#         overlap_matrix = adjacent_fe_obj.plot_overlap_matrix()
#         adj_overlap_matrix = np.array([[overlap_matrix.loc[win_lst[0],win_lst[0]], overlap_matrix.loc[win_lst[1],win_lst[0]]], 
#                                        [overlap_matrix.loc[win_lst[0],win_lst[1]], overlap_matrix.loc[win_lst[1],win_lst[1]]]])
#         print(adj_overlap_matrix)

class ANA_ALL_TAIL():
    
    def __init__(self, pair_path,):
        self.pair_path = pair_path
        self.cM2A_path_pattern = os.path.join(pair_path, 'complex', 'cM2A', '?.??', 'prod.out')
        self.cM2B_path_pattern = os.path.join(pair_path, 'complex', 'cM2B', '?.??', 'prod.out')
        self.lM2A_path_pattern = os.path.join(pair_path, 'ligands', 'lM2A', '?.??', 'prod.out')
        self.lM2B_path_pattern = os.path.join(pair_path, 'ligands', 'lM2B', '?.??', 'prod.out')
        
        self.cM2A_u_nk_pd = READ_PROD_OUT(self.cM2A_path_pattern).extract_data()
        self.cM2B_u_nk_pd = READ_PROD_OUT(self.cM2B_path_pattern).extract_data()
        self.lM2A_u_nk_pd = READ_PROD_OUT(self.lM2A_path_pattern).extract_data()
        self.lM2B_u_nk_pd = READ_PROD_OUT(self.lM2B_path_pattern).extract_data()
    
        self.cM2A_ANA_FEP_TOOLS = ANA_FEP_TOOLS(self.cM2A_u_nk_pd, )
        self.cM2B_ANA_FEP_TOOLS = ANA_FEP_TOOLS(self.cM2B_u_nk_pd, )
        self.lM2A_ANA_FEP_TOOLS = ANA_FEP_TOOLS(self.lM2A_u_nk_pd, )
        self.lM2B_ANA_FEP_TOOLS = ANA_FEP_TOOLS(self.lM2B_u_nk_pd, )
    
    
    def get_data_after_time_serial_anal(self, ana_fep_tools, lambda_, get_which_two_columns):
        fe_columns_name = get_which_two_columns[0]
        std_columns_name = get_which_two_columns[1]
        mov_df = ana_fep_tools.moving_esti_obj.df_dict_with_std[lambda_]
        for_df = ana_fep_tools.forward_esti_obj.df_dict_with_std[lambda_]
        rev_df = ana_fep_tools.reverse_esti_obj.df_dict_with_std[lambda_]
        forward_two_col_df = pd.concat([for_df[fe_columns_name], for_df[std_columns_name]],axis=1)
        reverse_two_col_df = pd.concat([rev_df[fe_columns_name], rev_df[std_columns_name]],axis=1)
        moving_two_col_df = pd.concat([mov_df[fe_columns_name], mov_df[std_columns_name]], axis=1)
        fe_he_std_dir = {'forward':forward_two_col_df, 
                         'reverse':reverse_two_col_df,
                         'moving':moving_two_col_df}
        return fe_he_std_dir
    
    def get_festd_df_cal_all_fe(self, df_cM2A, df_cM2B, df_lM2A, df_lM2B,):
        fe_df = df_cM2B.iloc[:,0]-df_cM2A.iloc[:,0]-df_lM2B.iloc[:,0]+df_lM2A.iloc[:,0]
        std_df = (df_cM2A.iloc[:,1]**2+df_cM2B.iloc[:,1]**2+df_lM2A.iloc[:,1]**2+df_lM2B.iloc[:,1]**2)**0.5
        all_estimated_FE_df = pd.concat([fe_df, std_df], axis = 1)
        return all_estimated_FE_df
        
    def time_serial_analysis(self, png_prefix, plot_plan = 3, ifall=True, ifflog_csv=True, fraction=1):
        get_which_two_columns=['free_energy(kcal/mol)', 'Rolling_STD(kcal/mol)']
        cM2A_png_prefix = png_prefix+'_cM2A'
        cM2B_png_prefix = png_prefix+'_cM2B'
        lM2A_png_prefix = png_prefix+'_lM2A'
        lM2B_png_prefix = png_prefix+'_lM2B'
        ###check the time convergence of the BAR calculated all processe free energy  first
        self.cM2A_ANA_FEP_TOOLS.use_bar_check_time_serial(fraction)
        self.cM2B_ANA_FEP_TOOLS.use_bar_check_time_serial(fraction)
        self.lM2A_ANA_FEP_TOOLS.use_bar_check_time_serial(fraction)
        self.lM2B_ANA_FEP_TOOLS.use_bar_check_time_serial(fraction)
        #output half time dg
        half_time_dg_cM2A_forw = self.cM2A_ANA_FEP_TOOLS.forward_esti_obj.output_half_time_dg(False)
        half_time_dg_cM2A_back = self.cM2A_ANA_FEP_TOOLS.reverse_esti_obj.output_half_time_dg(False)
        half_time_dg_cM2B_forw = self.cM2B_ANA_FEP_TOOLS.forward_esti_obj.output_half_time_dg(False)
        half_time_dg_cM2B_back = self.cM2B_ANA_FEP_TOOLS.reverse_esti_obj.output_half_time_dg(False)
        half_time_dg_lM2A_forw = self.lM2A_ANA_FEP_TOOLS.forward_esti_obj.output_half_time_dg(False)
        half_time_dg_lM2A_back = self.lM2A_ANA_FEP_TOOLS.reverse_esti_obj.output_half_time_dg(False)
        half_time_dg_lM2B_forw = self.lM2B_ANA_FEP_TOOLS.forward_esti_obj.output_half_time_dg(False)
        half_time_dg_lM2B_back = self.lM2B_ANA_FEP_TOOLS.reverse_esti_obj.output_half_time_dg(False)        
        half_time_dg_forw_back_minus = abs((half_time_dg_cM2B_forw-half_time_dg_cM2A_forw-half_time_dg_lM2B_forw+half_time_dg_lM2A_forw)-(half_time_dg_cM2B_back-half_time_dg_cM2A_back-half_time_dg_lM2B_back+half_time_dg_lM2A_back))
        with open(f'{png_prefix}_time_conv.txt', 'w') as f:
            print(half_time_dg_forw_back_minus, file=f)

        ###log_csv
        if ifflog_csv:
            self.cM2A_ANA_FEP_TOOLS.forward_esti_obj.output_csv(cM2A_png_prefix+'_bar_',)
            self.cM2A_ANA_FEP_TOOLS.moving_esti_obj.output_csv(cM2A_png_prefix+'_bar_', )
            self.cM2A_ANA_FEP_TOOLS.reverse_esti_obj.output_csv(cM2A_png_prefix+'_bar_',)
            self.cM2A_ANA_FEP_TOOLS.forward_esti_obj.output_csv(cM2B_png_prefix+'_bar_',)
            self.cM2B_ANA_FEP_TOOLS.moving_esti_obj.output_csv(cM2B_png_prefix+'_bar_', )
            self.cM2B_ANA_FEP_TOOLS.reverse_esti_obj.output_csv(cM2B_png_prefix+'_bar_',)
            self.cM2A_ANA_FEP_TOOLS.forward_esti_obj.output_csv(lM2A_png_prefix+'_bar_',)
            self.lM2A_ANA_FEP_TOOLS.moving_esti_obj.output_csv(lM2A_png_prefix+'_bar_', )
            self.lM2A_ANA_FEP_TOOLS.reverse_esti_obj.output_csv(lM2A_png_prefix+'_bar_',)
            self.cM2A_ANA_FEP_TOOLS.forward_esti_obj.output_csv(lM2B_png_prefix+'_bar_',)
            self.lM2B_ANA_FEP_TOOLS.moving_esti_obj.output_csv(lM2B_png_prefix+'_bar_', )
            self.lM2B_ANA_FEP_TOOLS.reverse_esti_obj.output_csv(lM2B_png_prefix+'_bar_',)
        ###cM2A
        cM2A_fe_he_std_dir = self.get_data_after_time_serial_anal(self.cM2A_ANA_FEP_TOOLS, '0.0 to 1.0', get_which_two_columns)
        forward_df_cM2A = cM2A_fe_he_std_dir['forward']
        reverse_df_cM2A = cM2A_fe_he_std_dir['reverse']
        moving_df_cM2A  = cM2A_fe_he_std_dir['moving']
        ###cM2B
        cM2B_fe_he_std_dir = self.get_data_after_time_serial_anal(self.cM2B_ANA_FEP_TOOLS, '0.0 to 1.0', get_which_two_columns)
        forward_df_cM2B = cM2B_fe_he_std_dir['forward']
        reverse_df_cM2B = cM2B_fe_he_std_dir['reverse']
        moving_df_cM2B  = cM2B_fe_he_std_dir['moving']
        ###lM2A
        lM2A_fe_he_std_dir = self.get_data_after_time_serial_anal(self.lM2A_ANA_FEP_TOOLS, '0.0 to 1.0', get_which_two_columns)
        forward_df_lM2A = lM2A_fe_he_std_dir['forward']
        reverse_df_lM2A = lM2A_fe_he_std_dir['reverse']
        moving_df_lM2A  = lM2A_fe_he_std_dir['moving']
        ###lM2B
        lM2B_fe_he_std_dir = self.get_data_after_time_serial_anal(self.lM2B_ANA_FEP_TOOLS, '0.0 to 1.0', get_which_two_columns)
        forward_df_lM2B = lM2B_fe_he_std_dir['forward']
        reverse_df_lM2B = lM2B_fe_he_std_dir['reverse']
        moving_df_lM2B  = lM2B_fe_he_std_dir['moving']
        ###PLOTTING SUM FREE ENERGY
        plot_obj = PLOTTING()
        if plot_plan == 2:
            fe_he_std_dir = {'forward':self.get_festd_df_cal_all_fe(forward_df_cM2A, forward_df_cM2B, forward_df_lM2A, forward_df_lM2B),
                             'reverse':self.get_festd_df_cal_all_fe(reverse_df_cM2A, reverse_df_cM2B, reverse_df_lM2A, reverse_df_lM2B),}
            plot_obj.plot_fe_time_serial('{}_bar_sum_all_2.png'.format(png_prefix), **fe_he_std_dir)
        elif plot_plan == 3:
            fe_he_std_dir = {'forward':self.get_festd_df_cal_all_fe(forward_df_cM2A, forward_df_cM2B, forward_df_lM2A, forward_df_lM2B), 
                             'reverse':self.get_festd_df_cal_all_fe(reverse_df_cM2A, reverse_df_cM2B, reverse_df_lM2A, reverse_df_lM2B),
                             'moving' :self.get_festd_df_cal_all_fe(moving_df_cM2A, moving_df_cM2B, moving_df_lM2A, moving_df_lM2B)}
            plot_obj.plot_fe_time_serial('{}_bar_sum_all_3.png'.format(png_prefix), **fe_he_std_dir)
        elif plot_plan == 1:
            fe_he_std_dir_moving = {'moving' :self.get_festd_df_cal_all_fe(moving_df_cM2A, moving_df_cM2B, moving_df_lM2A, moving_df_lM2B),}
            plot_obj.plot_fe_time_serial('{}_bar_sum_all_1_moving.png'.format(png_prefix), **fe_he_std_dir_moving)
            fe_he_std_dir_forward = {'forward':self.get_festd_df_cal_all_fe(forward_df_cM2A, forward_df_cM2B, forward_df_lM2A, forward_df_lM2B), }
            plot_obj.plot_fe_time_serial('{}_bar_sum_all_1_forward.png'.format(png_prefix), **fe_he_std_dir_forward)
            fe_he_std_dir_reverse = {'reverse':self.get_festd_df_cal_all_fe(reverse_df_cM2A, reverse_df_cM2B, reverse_df_lM2A, reverse_df_lM2B),}
            plot_obj.plot_fe_time_serial('{}_bar_sum_all_1_reverse.png'.format(png_prefix), **fe_he_std_dir_reverse)
        else:
            print('Error! No such plot plan!')
        ###check if output every simulation windows convergence    
        if ifall:
            self.cM2A_ANA_FEP_TOOLS.use_fep_check_time_serial(3, cM2A_png_prefix, use_forward=True,)
            self.cM2B_ANA_FEP_TOOLS.use_fep_check_time_serial(3, cM2B_png_prefix, use_forward=True,)
            self.lM2A_ANA_FEP_TOOLS.use_fep_check_time_serial(3, lM2A_png_prefix, use_forward=True,)
            self.lM2B_ANA_FEP_TOOLS.use_fep_check_time_serial(3, lM2B_png_prefix, use_forward=True,)
            if ifflog_csv:
                self.cM2A_ANA_FEP_TOOLS.moving_esti_obj.output_csv(cM2A_png_prefix+'_fep_', )
                self.cM2A_ANA_FEP_TOOLS.reverse_esti_obj.output_csv(cM2A_png_prefix+'_fep_',)
                self.cM2B_ANA_FEP_TOOLS.moving_esti_obj.output_csv(cM2B_png_prefix+'_fep_', )
                self.cM2B_ANA_FEP_TOOLS.reverse_esti_obj.output_csv(cM2B_png_prefix+'_fep_',)
                self.lM2A_ANA_FEP_TOOLS.moving_esti_obj.output_csv(lM2A_png_prefix+'_fep_', )
                self.lM2A_ANA_FEP_TOOLS.reverse_esti_obj.output_csv(lM2A_png_prefix+'_fep_',)
                self.lM2B_ANA_FEP_TOOLS.moving_esti_obj.output_csv(lM2B_png_prefix+'_fep_', )
                self.lM2B_ANA_FEP_TOOLS.reverse_esti_obj.output_csv(lM2B_png_prefix+'_fep_',)

    

    
    def cal_final_fe_and_std(self, part_lambda_list, scale_f, fe_mode, std_mode, filename='part_lambda.csv'):
        # scale_f = 0.75
        # fe_mode = 'BAR'
        # std_mode = 'time_serial'
        cM2A_part_lambda_df = self.cM2A_ANA_FEP_TOOLS.check_part_lambda(part_lambda_list, scale_f, fe_mode, std_mode, None)
        cM2B_part_lambda_df = self.cM2B_ANA_FEP_TOOLS.check_part_lambda(part_lambda_list, scale_f, fe_mode, std_mode, None)
        lM2A_part_lambda_df = self.lM2A_ANA_FEP_TOOLS.check_part_lambda(part_lambda_list, scale_f, fe_mode, std_mode, None)
        lM2B_part_lambda_df = self.lM2B_ANA_FEP_TOOLS.check_part_lambda(part_lambda_list, scale_f, fe_mode, std_mode, None)
        all_estimated_FE_df = cM2B_part_lambda_df.iloc[:,0]-cM2A_part_lambda_df.iloc[:,0]-lM2B_part_lambda_df.iloc[:,0]+lM2A_part_lambda_df.iloc[:,0]
        all_std_df = (cM2B_part_lambda_df.iloc[:,1]**2+cM2A_part_lambda_df.iloc[:,1]**2+lM2B_part_lambda_df.iloc[:,1]**2+lM2A_part_lambda_df.iloc[:,1]**2)**0.5
        all_estimated_FE_df = pd.concat([all_estimated_FE_df, all_std_df], axis = 1)
        if filename:
            all_estimated_FE_df.to_csv(filename)
        return all_estimated_FE_df


if  __name__ == '__main__':
    from READ_PROD_OUT import *
    ##only keep this for jupyter notebook testing
    opts = optParser('')
    # fakeArgs = '-d F:\\rbfe_test\\rbfe_test\\mcl1\\run-1\\38-60\\FEP -f 0.75 -o 38-60_part_lambda_fe.csv'
    # opts = optParser(fakeArgs.strip().split())
    file_str = str(opts.option.file_directory)
    all_ts_obj = ANA_ALL_TAIL(file_str)
    part_lambda_list = [[0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.00], ]
    scale_f = float(opts.option.fraction)
    fe_mode = 'BAR'
    std_mode = 'time_serial'
    estimate_df_filename = str(opts.option.output_csv_filename)
    pair_name = estimate_df_filename.split('_')[0]
    all_ts_obj.time_serial_analysis(pair_name, 1, True, True, scale_f)
    all_ts_obj.cal_final_fe_and_std(part_lambda_list, scale_f, fe_mode, std_mode, estimate_df_filename)

