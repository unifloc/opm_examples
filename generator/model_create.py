from generator_commands import DataParser
import os
from ecl.summary import EclSum
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import rips
import subprocess
import time
import glob
from PIL import Image
from IPython.display import display
from plotly.offline import plot, iplot
import shutil
import math as m
"""
@alexeyvodopyan
23.09
@sevrn
"""


def clear_folders():
    subprocess.call("rm -f model_folder/*", shell=True)
    subprocess.call("rm -f csv_folder/*", shell=True)


class ModelGenerator:
    """
    
    """

    def __init__(self, init_file_name='RIENM1_INIT.DATA', start_date="1 'SEP' 2020", mounths = 12, days = 30, nx=100, ny=100, nz=5, dx=500, dy=500, dz=20, por=0.3, permx=100,
                 permy=100, permz=10, prod_names=None, prod_xs=None, prod_ys=None, prod_z1s=None, prod_z2s=None, prod_q_oil=None,
                 inj_names=None, inj_xs=None, inj_ys=None, inj_z1s=None, inj_z2s=None, inj_bhp=None, skin=None, density=None,
                 p_depth = 2510, p_init = 320, o_w_contact = 2600, pc_woc = 0, g_o_contact = 2450, pc_goc = 0, tops_depth = 2500, rezim = 'ORAT', prod_bhp = None, horizontal = None, y_stop = None, only_prod = False,
                 lgr = False, lx = None, ly = None, cells_cy = None, cells_v = None, cells_cx=None, indicator=None, upr_rezim_water=False, upr_rezim_gas=False, neodn_plast=False,
                 neodn_plast_tamp=None, rw=None):
        # паарметры расчета
        self.start_date = start_date
        self.mounths = mounths
        self.days = days

        # размеры модели
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.tops_depth = tops_depth
        
        # физические свойства
        self.por = por
        self.permx = permx
        self.permy = permy
        self.permz = permz

        # EQUILIBRIUM DATA
        self.p_depth = p_depth
        self.p_init = p_init
        self.o_w_contact = o_w_contact
        self.pc_woc = pc_woc
        self.g_o_contact = g_o_contact
        self.pc_goc = pc_goc

        # свойства нефти
        self.density = density

        # параметры скважин
        # добывающие
        self.prod_names = prod_names
        self.rezim = rezim
        self.prod_xs = prod_xs
        self.prod_ys = prod_ys
        self.prod_z1s = prod_z1s
        self.prod_z2s = prod_z2s
        self.q_oil = prod_q_oil
        self.prod_bhp = prod_bhp

        # нагнетательные 
        self.inj_names = inj_names
        self.inj_ys = inj_ys
        self.inj_xs = inj_xs
        self.inj_z1s = inj_z1s
        self.inj_z2s = inj_z2s
        self.inj_bhp = inj_bhp

        self.skin = skin # скин
        self.rw = rw # радиус скважины
        
        # горизонтальные (стоит взглянуть на compdat)
        self.horizontal = horizontal
        self.y_stop = y_stop # начало задается prod_ys

        # LGR (в разработке)
        self.lgr = lgr 
        self.lx = lx
        self.ly = ly
        self.cells_cy = cells_cy
        self.cells_cx = cells_cx
        self.cells_v = cells_v

        # индикаторы режимов (в разработке)
        self.only_prod = only_prod # оставляет только добывающие скважины
        self.upr_rezim_water = upr_rezim_water # умножает на 1000 пористость последнего пропласта (предварительно необходимо установить ВНК)
        self.upr_rezim_gas = upr_rezim_gas # умножает на 100 пористость первого пропластка (предварителньо необходимо установить ГНК)
        # моделирование тампонирования
        self.neodn_plast_tamp = neodn_plast_tamp # значение проницаемости в ПЗ высокопрониц. слоя доб. скв. 
        self.neodn_plast = neodn_plast # два несвязанных пропласта с разной проницаемостью

        # переменные для расчета
        self.indicator = indicator
        self.result_df = None
        self.fig = None
        self.fig_npv = None
        self.dir = None
        self.init_file_name = init_file_name
        self.filter_initial_data()
        self.parser = None
        self.initialize_parser(self.init_file_name, self.start_date, self.mounths, self.days, self.nx, self.ny, self.nz, self.dx, self.dy, self.dz, self.por,
                               self.permx, self.permy, self.permz, self.prod_names, self.prod_xs, self.prod_ys,
                               self.prod_z1s, self.prod_z2s, self.q_oil, self.inj_names, self.inj_xs, self.inj_ys, self.inj_z1s,
                               self.inj_z2s, self.inj_bhp, self.skin, self.density, self.p_depth, self.p_init, self.o_w_contact, self.pc_woc, 
                               self.g_o_contact, self.pc_goc, self.tops_depth, self.rezim, self.prod_bhp, self.horizontal, self.y_stop, self.only_prod,
                               self.lgr, self.lx, self.ly, self.cells_cy, self.cells_v, self.cells_cx, self.upr_rezim_water, self.upr_rezim_gas, self.neodn_plast,
                               self.neodn_plast_tamp, self.rw)

     
    def initialize_parser(self, init_file_name, start_date, mounths, days, nx, ny, nz, dx, dy, dz, por, permx, permy, permz, prod_names, prod_xs,
                          prod_ys, prod_z1s, prod_z2s, q_oil, inj_names, inj_xs, inj_ys, inj_z1s, inj_z2s, inj_bhp, skin, density, p_depth, 
                          p_init, o_w_contact, pc_woc, g_o_contact, pc_goc, tops_depth, rezim, prod_bhp, horizontal, y_stop, only_prod,
                          lgr, lx, ly, cells_cy, cells_v, cells_cx, upr_rezim_water, upr_rezim_gas, neodn_plast, neodn_plast_tamp, rw):
        init_file = open(init_file_name)
        self.parser = DataParser(init_file, start_date, mounths, days, nx, ny, nz, dx, dy, dz, por, permx, permy, permz, prod_names, prod_xs,
                                 prod_ys, prod_z1s, prod_z2s, q_oil, inj_names, inj_xs, inj_ys, inj_z1s, inj_z2s, inj_bhp, skin, density, p_depth, 
                                 p_init, o_w_contact, pc_woc, g_o_contact, pc_goc, tops_depth, rezim, prod_bhp, horizontal, y_stop, only_prod,
                                 lgr, lx, ly, cells_cy, cells_v, cells_cx, upr_rezim_water, upr_rezim_gas, neodn_plast, neodn_plast_tamp, rw)

    def filter_initial_data(self):
        if max(self.prod_xs) > self.nx:
            print('Х-координата скважин больше Х-размерности модели. Проверьте свои данные')
        if max(self.prod_ys) > self.ny:
            print('Y-координата скважин больше Y-размерности модели. Проверьте свои данные')
        if max(self.prod_z2s) > self.nz:
            print('Z2-координата скважин больше Z-размерности модели. Проверьте свои данные')

    def create_model(self, name, result_name, keys):
        self.parser.parse_file('DIMENS')
        self.parser.parse_file('START')
        self.parser.parse_file('DX')
        self.parser.parse_file('DY')
        self.parser.parse_file('DZ')
        self.parser.parse_file('TOPS')
        self.parser.parse_file('PORO')
        self.parser.parse_file('PERMX')
        self.parser.parse_file('PERMY')
        self.parser.parse_file('PERMZ')
        self.parser.parse_file('DENSITY')
        self.parser.parse_file('EQUIL')
        self.parser.parse_file('WELSPECS')
        self.parser.parse_file('COMPDAT')
        self.parser.parse_file('WCONPROD')
        if not self.only_prod:
            self.parser.parse_file('WCONINJE')
        self.parser.parse_file('TSTEP')
        self.save_file(name=name)
        self.calculate_file(name)
        self.create_result(name=name, keys=keys, indicator=self.indicator)
        self.read_result(name=result_name)
        self.export_snapshots(name=name)   


    def save_file(self, name):
        with open('model_folder/'+ name + '.DATA', "w") as file:
            for line in self.parser.content:
                print(line, file=file)

    @staticmethod
    def calculate_file(name):
        os.system("flow model_folder/%s.DATA" % name)

    @staticmethod
    def create_result(name, keys, indicator):
        summary = EclSum('model_folder/%s.DATA' % name)
        dates = summary.dates
        results = []
        all_keys = []

        if keys is None:
            keys = ["WOPR:*"]

        for key in keys:
            key_all_wells = summary.keys(key)
            all_keys = all_keys + list(key_all_wells)

        for key in all_keys:
            results.append(list(summary.numpy_vector(key)))

        if len(results) == 0:
            return print('Результаты из модели не загрузились. Файл с результатами не был создан')

        result_df = pd.DataFrame(data=np.array(results).T, index=dates, columns=all_keys)
        result_df.index.name = 'Time'
        if indicator is not None:
            if os.path.exists(f'csv_folder/{indicator}') is False:
                os.mkdir(f'csv_folder/{indicator}')
            result_df.to_csv(f'csv_folder/{indicator}/%s.csv' % name)
        else:
            result_df.to_csv('csv_folder/%s.csv' % name)
        print('%s_RESULT.csv is created' % name)


    def read_result(self, name):
        if self.indicator is not None:
            self.result_df = pd.read_csv(f'csv_folder/{self.indicator}/%s.csv' % name, parse_dates=[0], index_col=[0])
        else:
            self.result_df = pd.read_csv('csv_folder/%s.csv' % name, parse_dates=[0], index_col=[0])
        print('%s.csv is read' % name)
        

# метод построения графика по заданному параметру
    def summ_plot(self, parameters=None, mode='lines', x_axis=None, y_axis=None, title=None, name=None, npv_ind=False, l_list=None):
        directory = "csv_folder/"
        npv_list = []
        model_list = []
        self.fig = go.Figure()
        files = [f for f in os.listdir(directory)]
        files.sort(key=lambda x:int(x.split('.')[0][-1]))
        for file in files:
            df = pd.read_csv('csv_folder/%s' % file, parse_dates=[0], index_col=[0])
            i = int(file.split('.')[0][-1])
            if npv_ind == True:
                npv_ = self.npv_method(df, l_list[i])
                npv_list.append(npv_)
                model_list.append(f'Модель: {name[i]}')
            self.fig.add_trace(go.Scatter(
                    x=df.index,
                    y=df[parameters[0]],
                    mode=mode,
                    name='Модель:' + name[i]))
        if x_axis is None:
            x_axis = 'Дата'
        if y_axis is None:
            y_axis = ''
        if title is None:
            title = 'Динамика показателей по моделям'
        self.fig.update_xaxes(tickformat='%d.%m.%y')
        self.fig.update_layout(title=go.layout.Title(text=title),
                               xaxis_title=x_axis,
                               yaxis_title=y_axis)
        colors = ['lightslategray',] * 10
        colors[6] = 'crimson'
        data = [go.Bar(
            x = model_list,
            y = npv_list,
            marker_color=colors)]
        self.fig_npv = go.Figure(data=data)
        self.fig_npv.update_layout(title='NPV по моделям')
        self.fig_npv.update_yaxes(type="log")


# метод для создания скриншота с сеткой
    def export_snapshots(self, name):
        process = subprocess.Popen('exec ResInsight --case "model_folder/%s.EGRID"' % name, shell=True)
        time.sleep(5)
        resinsight = rips.Instance.find()
        case = resinsight.project.cases()[0]
        resinsight.set_main_window_size(width=400, height=150)
        property_list = ['PRESSURE', 'SOIL']
        case_path = case.file_path
        folder_name = os.path.dirname(case_path)

        # create a folder to hold the snapshots
        dirname = os.path.join(folder_name, f"snapshots/")
        self.dir = dirname
        # if os.path.exists(dirname) is False:
        #     os.mkdir(f"model_folder/snapshots/{name}")
        # shutil.rmtree(dirname)

        print("Exporting to folder: " + dirname)
        resinsight.set_export_folder(export_type='SNAPSHOTS', path=dirname)

        view = case.views()[0]
        time_steps = case.time_steps()
        l = len(time_steps) - 1
        for property in property_list:
            view.apply_cell_result(result_type='DYNAMIC_NATIVE', result_variable=property)
            view.set_time_step(time_step=l)
            view.export_snapshot()

        process.kill()

    def iplot_fig(self, npv=False):
        if npv == True:
            iplot(self.fig_npv)
        elif npv == False:
            iplot(self.fig)

    def plot_fig(self):
        plot(self.fig)

# метод для отображения сетки 
    def display_grids(self):
        images = []
        image_paths = glob.glob(self.dir + '*')

        for path in image_paths:
            images.append(Image.open(path))

        for grid in images:
            display(grid)

# методы расчета NPV для исследования скважин
    def npv_method(self, df, l):
        ci = 170*10**6 # руб, капитальные затраты на строительство скважины c поверхностным обустройством;
        cap_l = 40000 # РУБ, стоимость 1 метра горизонтального ствола;
        p = 4000 # руб/м3, net-baсk цена нефти за вычетом НПДИ и подготовку нефти; 
        opex = 10**6 # руб/год, операционные затраты на скважину;
        r = 0.12 # ставка дисконтирования;
        to = df.index[0]
        i = 0
        npv = -ci - l*cap_l
        j = 0
        q = 0
        for t in df.index:
            if abs((t - to).days) >= 365:
                i += 1
                to = t
                q = df['FOPT'][j] - q
                dcf = (q*p - opex)/(1 + r)**i
                npv += dcf 
            j += 1

        return round(npv, 0)

    def npv_plotn_method(self, df, l, A):
        ci = 170*10**6 # руб, капитальные затраты на строительство скважины c поверхностным обустройством;
        cap_l = 40000 # РУБ, стоимость 1 метра горизонтального ствола;
        p = 4000 # руб/м3, net-baсk цена нефти за вычетом НПДИ и подготовку нефти; 
        opex = 10**6 # руб/год, операционные затраты на скважину;
        r = 0.12 # ставка дисконтирования;
        to = df.index[0]
        i = 0
        npv = (-ci - l*cap_l)/A
        j = 0
        q = 0
        for t in df.index:
            if abs((t - to).days) >= 365:
                i += 1
                to = t
                q = df['FOPT'][j] - q
                dcf = (q*p - opex)/A/(1 + r)**i
                npv += dcf 
            j += 1

        return round(npv, 0)

# метод для отображения графика с оптимальной плотностью сетки
    def summ_plot_plotn(self, parameters=None, mode='lines', x_axis=None, y_axis=None, title=None, A=None, name=None, npv_ind=False, l=None):
        directory = "csv_folder/"
        npv_list = []
        model_list = []
        self.fig = go.Figure()
        files = [f for f in os.listdir(directory)]
        files.sort(key=lambda x:int(x.split('.')[1]))
        i = 0
        for file in files:
            df = pd.read_csv('csv_folder/%s' % file, parse_dates=[0], index_col=[0])
            j = int(file.split('.')[1])
            y = df[parameters[0]]/A[i]
            if npv_ind == True:
                npv_ = self.npv_plotn_method(df, l, A[i])
                npv_list.append(npv_)
                model_list.append(f'Модель: {name[i]}')
            self.fig.add_trace(go.Scatter(
                    x=df.index,
                    y=y,
                    mode=mode,
                    name='Модель:' + name[j]))
            i += 1
        if x_axis is None:
            x_axis = 'Дата'
        if y_axis is None:
            y_axis = ''
        if title is None:
            title = 'Динамика показателей по моделям'
        self.fig.update_xaxes(tickformat='%d.%m.%y')
        self.fig.update_layout(title=go.layout.Title(text=title),
                               xaxis_title=x_axis,
                               yaxis_title=y_axis)
        colors = ['lightslategray',] * 6
        colors[2] = 'crimson'
        data = [go.Bar(
            x = model_list,
            y = npv_list,
            marker_color=colors)]
        self.fig_npv = go.Figure(data=data)
        self.fig_npv.update_layout(title='NPV по моделям')

# метод для построения графика с оптимальным соотношнием сторон
    def sootn_plot(self, ls, par, k, h, mu, A):
        self.fig = go.Figure()
        x_opt = []
        y_opt = []
        for i in range(0, 6):
            P = []
            directory = f"csv_folder/{i}"
            files = [f for f in os.listdir(directory)]
            files.sort(key=lambda x:int(x.split('.')[2]))
            
            for file in files:
                df = pd.read_csv(f'csv_folder/{i}/%s' % file, parse_dates=[0], index_col=[0])
    
                val = k*h*(df['FPR'][-1]-df['WBHP:P1'][-1])/(df['WOPR:P1'][-1]*mu)*10**5*86400
                val = val + h*(10)**0.5/2/3.14/400*(m.log(h*(10)**0.5/(2*3.14*0.5*0.146*(1+(10)**0.5)*m.sin(3.14/2)))+ 0) # 0 в конце это скин
                P.append(val)
            self.fig.add_trace(go.Scatter(
                x=ls[i],
                y=P,
                mode='lines',
                name='Параметр - A/L^2 = ' + str(par[i])))
            x_opt.append(A[i]*10000/(A[i]*10000+400**2))
            y_opt.append(m.log(1+A[i]*10000/(400**2))*4) # убрать *2
        x_opt.append(384*10000/(384*10000+400**2))
        y_opt.append(m.log(1+384*10000/(400**2))*4) # убрать *2
        self.fig.add_trace(go.Scatter(
                x=x_opt,
                y=y_opt,
                mode='lines',
                name='Оптимальное соотношение геометрических размеров'
            ))
        x_axis = 'Соотношение геометрических размеров области дренирования (ширина/длина)'
        y_axis = 'Безразмерный перепад давлений P*'
        title = 'Безразмерный перепад давления как функция соотношения геометрических размеров области дренирования'
        self.fig.update_layout(title=go.layout.Title(text=title),
                               xaxis_title=x_axis,
                               yaxis_title=y_axis)