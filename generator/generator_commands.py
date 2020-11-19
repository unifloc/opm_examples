"""
@alexeyvodopyan
23.09

"""


class DataParser:
    """

    """
    def __init__(self, init_file, start_date, mounths, days, nx, ny, nz, dx, dy, dz, por, permx, permy, permz, prod_names, prod_xs,
                 prod_ys, prod_z1s, prod_z2s, q_oil, inj_names, inj_xs, inj_ys, inj_z1s, inj_z2s, inj_bhp, skin, density, p_depth,
                 p_init, o_w_contact, pc_woc, g_o_contact, pc_goc, tops_depth, rezim, prod_bhp, horizontal, y_stop, only_prod,
                 lgr, lx, ly, cells_c):
        
        self.content = init_file.readlines()
        'Удалим переносы в конце строк'
        self.content = [line.rstrip('\n') for line in self.content]
        self.content = [line.strip() for line in self.content]
        self.start = None
        self.TSTEP = None
        self.dimens = None
        self.dx_dim = None
        self.dy_dim = None
        self.dz_dim = None
        self.tops_dim = None
        self.por_dim = None
        self.permx_dim = None
        self.permy_dim = None
        self.permz_dim = None
        self.den = None
        self.equil = None
        self.welspecs = None
        self.compdat = None
        self.wconprod = None
        self.wconinje = None
        self.start_date = start_date
        self.mounths = mounths
        self.days = days
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.tops_depth = tops_depth
        self.por = por
        self.permx = permx
        self.permy = permy
        self.permz = permz
        self.density = density
        self.p_depth = p_depth
        self.p_init = p_init
        self.o_w_contact = o_w_contact
        self.pc_woc = pc_woc
        self.g_o_contact = g_o_contact
        self.pc_goc = pc_goc
        self.prod_names = prod_names
        self.rezim = rezim
        self.prod_xs = prod_xs
        self.prod_ys = prod_ys
        self.prod_z1s = prod_z1s
        self.prod_z2s = prod_z2s
        self.q_oil = q_oil
        self.prod_bhp = prod_bhp
        self.inj_names = inj_names
        self.inj_xs = inj_xs
        self.inj_ys = inj_ys
        self.inj_z1s = inj_z1s
        self.inj_z2s = inj_z2s
        self.inj_bhp = inj_bhp
        if not only_prod:
            self.all_well_names = self.prod_names + self.inj_names
            self.all_well_xs = self.prod_xs + self.inj_xs
            self.all_well_ys = self.prod_ys + self.inj_ys
            self.all_well_z1s = self.prod_z1s + self.inj_z1s
            self.all_well_z2s = self.prod_z2s + self.inj_z2s
        else:
            self.all_well_names = self.prod_names
            self.all_well_xs = self.prod_xs
            self.all_well_ys = self.prod_ys 
            self.all_well_z1s = self.prod_z1s
            self.all_well_z2s = self.prod_z2s
        self.all_well_fluid = ['OIL' for _ in range(len(self.prod_names))] + ['WAT' for _ in range(len(self.inj_names))]
        self.skin = skin
        self.horizontal = horizontal
        self.y_stop = y_stop
        self.lgr = lgr 
        self.lx = lx
        self.ly = ly
        self.cells_c = cells_c

    def parse_file(self, keyword):
        keyword_start_flag = False
        i = 0
        for line in self.content:
            if line[:2] == '--':  # comment in file detected
                i += 1
                continue
            elif keyword in line:
                keyword_start_flag = True
                print('%s detected' % keyword)
            elif keyword_start_flag and keyword == 'DIMENS':
                self.create_dimensions()
                self.content[i] = self.dimens
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'START':
                self.create_start_date()
                self.content[i] = self.start
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'DX':
                self.create_dx_dim()
                self.content[i] = self.dx_dim 
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'DY':
                self.create_dy_dim()
                self.content[i] = self.dy_dim
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'DZ':
                self.create_dz_dim()
                self.content[i] = self.dz_dim
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'TOPS':
                self.create_tops()
                self.content[i] = self.tops_dim
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'PORO':
                self.create_porosity()
                self.content[i] = self.por_dim
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'PERMX':
                self.create_permx()
                self.content[i] = self.permx_dim
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'PERMY':
                self.create_permy()
                self.content[i] = self.permy_dim
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'PERMZ':
                self.create_permz()
                self.content[i] = self.permz_dim
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'DENSITY':
                self.create_density()
                self.content[i] = self.den
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'EQUIL':
                self.create_equil()
                self.content[i] = self.equil
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'WELSPECS':
                for prod, x, y, fluid in zip(self.all_well_names, self.all_well_xs,
                                             self.all_well_ys, self.all_well_fluid):
                    self.create_welspecs(prod, x, y, fluid)
                    self.content.insert(i, self.welspecs)
                    i += 1
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'COMPDAT':
                for well, x, y, z1, z2, skin, horizontal, y_stop in zip(self.all_well_names, self.all_well_xs, self.all_well_ys,
                                              self.all_well_z1s, self.all_well_z2s, self.skin, self.horizontal, self.y_stop):
                    self.create_compdat(well, x, y, z1, z2, skin, horizontal, y_stop)
                    self.content.insert(i, self.compdat)
                    i += 1
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'WCONPROD':
                for prod, rezim, q_oil, prod_bhp in zip(self.prod_names, self.rezim, self.q_oil, self.prod_bhp):
                    self.create_wconprod(prod, rezim, q_oil, prod_bhp)
                    self.content.insert(i, self.wconprod)
                    i += 1
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'WCONINJE':
                for inj, inj_bhp in zip(self.inj_names, self.inj_bhp):
                    self.create_wconinje(inj, inj_bhp)
                    self.content.insert(i, self.wconinje)
                    i += 1
                keyword_start_flag = self.keyword_read(keyword)
            elif keyword_start_flag and keyword == 'TSTEP':
                self.create_TSTEP()
                self.content[i] = self.TSTEP
                keyword_start_flag = self.keyword_read(keyword)
            i += 1

    @staticmethod
    def keyword_read(keyword):
        print('%s written' % keyword)
        keyword_start_flag = False
        return keyword_start_flag
    
    @staticmethod
    def setca(nx, lx, s):
        k = 1
        l = 0
        while abs(l-lx) > 10:
            k += 0.001
            n = round((nx-s)/2)
            l = 8 + 2*(8*k*(k**n-1)/(k-1))
            if abs(l-lx) < 10:
                x = []
                for i in range(1, n+1):
                    x.append(round(8*k**i))
                x = x[::-1] + [8]*s + x
                return x

    def create_tops(self):
        self.tops_dim = str(self.nx*self.ny) + '*'+ str(self.tops_depth) + ' /'  

    def create_start_date(self):
        self.start = self.start_date + ' /'

    def create_dimensions(self):
        self.dimens = str(self.nx) + ' ' + str(self.ny) + ' ' + str(self.nz) + ' /'

    def create_dx_dim(self):
        if self.lgr:
            dx_lgr = self.setca(self.nx, self.lx, self.cells_c)
            self.dx_dim = str(dx_lgr[0]) + ' \n'
            for i in range(2, self.nx + 1):    
                self.dx_dim += str(dx_lgr[i-1]) + ' \n'
            self.dx_dim = self.dx_dim*self.nx
            self.dx_dim += '/\n'
        else:
            self.dx_dim = str(self.nx*self.ny*self.nz) + '*' + str(self.dx) + ' /'

    def create_dy_dim(self):
        if self.lgr:
            dy_lgr = self.setca(self.ny, self.ly, self.cells_c)
            self.dy_dim = ''
            for i in range(1, self.ny + 1):    
                dim = str(dy_lgr[i-1]) + ' \n'
                self.dy_dim += dim*self.ny
            self.dy_dim += '/\n'
        else:
            self.dy_dim = str(self.nx*self.ny*self.nz) + '*' + str(self.dy) + ' /'

    def create_dz_dim(self):
            self.dz_dim = str(self.nx*self.ny*self.nz) + '*' + str(self.dz) + ' /'

    def create_porosity(self):
        self.por_dim = str(self.nx*self.ny*self.nz) + '*' + str(self.por) + ' /'

    def create_permx(self):
        self.permx_dim = str(self.nx*self.ny*self.nz) + '*' + str(self.permx) + ' /'

    def create_permy(self):
        self.permy_dim = str(self.nx*self.ny*self.nz) + '*' + str(self.permy) + ' /'

    def create_permz(self):
        self.permz_dim = str(self.nx*self.ny*self.nz) + '*' + str(self.permz) + ' /'

    def create_density(self):
        self.den = str(self.density[0]) + ' ' + str(self.density[1]) + ' ' + str(self.density[2]) + ' /'
        
    def create_equil(self):
        self.equil = str(self.p_depth) + ' ' + str(self.p_init) + ' ' + str(self.o_w_contact) + ' ' + str(self.pc_woc) + ' ' + str(self.g_o_contact) + ' ' + str(self.pc_goc) + ' ' +' 1 1* 1*/'

    def create_welspecs(self, name, x, y, fluid):
        self.welspecs = name + ' G1 ' + str(x) + ' ' + str(y) + ' * ' + fluid + ' /'

    def create_compdat(self, name, x, y, z1, z2, skin, horizontal, y_stop):
        if horizontal:
            #self.compdat = name + ' ' + '2*' + ' ' + str(z1) + ' ' + str(z2) + ' OPEN	1*	1*	0.5  1* ' + str(skin) + ' /\n'
            self.compdat = name + ' ' + str(x) + ' ' + str(y) + ' ' + str(z2) + ' ' + str(z2) + ' OPEN	1*	1*	0.5  1* ' + str(skin) + ' /\n'
            for i in range(y+1, y_stop+1):
                  self.compdat += name + ' ' + str(x) + ' ' + str(i) + ' ' + str(z2) + ' ' + str(z2) + ' OPEN	1*	1*	0.5  1* ' + str(skin) + ' /\n'
        else:
            self.compdat = name + ' ' + str(x) + ' ' + str(y) + ' ' + str(z1) + ' ' + str(z2) + ' OPEN	1*	1*	0.5  1* ' + str(skin) + ' /'

    def create_wconprod(self, name, rezim, q_oil, prod_bhp):
        self.wconprod = name + ' OPEN ' + rezim + ' ' + str(q_oil) + ' 4* ' + str(prod_bhp) + ' /'

    def create_wconinje(self, name, inj_bhp):
        self.wconinje = name + ' WAT OPEN BHP ' + str(inj_bhp) + ' 1* /'

    def create_TSTEP(self):
        self.TSTEP = str(self.mounths) + f'*{self.days} /'
