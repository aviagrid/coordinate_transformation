from math import pi, sin, sqrt, cos, asin, trunc
import numpy as np


class CoordTransform:
    """ Класс для выполнения преобразования координат"""
    def __init__(self, lat_grad: list, long_grad: list, alt_m: float):
        """Входные параметры для координат

        Args:
            lat_grad (list): 
                Широта точки в градусах может быть в формате (DD, MM, SS) или в (DD.DDDD)
            long_grad (list): 
                Долгота точки в градусах может быть в формате (DD, MM, SS) или в (DD.DDDD)
            alt_m (float): 
                Высота точки в метрах      
        """
        # ----- Входные параметры ----- #
        self.lat_grad = lat_grad
        self.long_grad = long_grad
        self.alt_m = alt_m
        
        # ----- Параметры элипсоида ПЗ-90.11 ----- #
        self.a_ellips_pz = (6378136, )
        self.a_compr_pz = (1/298.25784, )
        self.square_eccent_pz = ((2*self.a_compr_pz[0]) - (self.a_compr_pz[0] ** 2), )
        self.eccent_pz = (sqrt(self.square_eccent_pz[0]), )
        
        # ----- Параметры элипсоида WGS-84 ----- #
        self.a_ellips_wgs = (6378137, )
        self.a_compr_wgs = (1/298.257223563, )
        self.square_eccent_wgs = ((2*self.a_compr_wgs[0]) - (self.a_compr_wgs[0] ** 2), )
        self.eccent_wgs = (sqrt(self.square_eccent_wgs[0]), )
        
        # ----- Параметры элипсоида СК-42 ----- #
        self.a_ellips_sk = (6378245, )
        self.a_compr_sk = (1/298.3, )
        self.square_eccent_sk = ((2*self.a_compr_sk[0]) - (self.a_compr_sk[0] ** 2), )
        self.eccent_sk = (sqrt(self.square_eccent_sk[0]), )
        
        # ----- Угловые коэффициенты СК-42, ПЗ-90.11 ---- #
        omega_x_sk_pz = 1.115071 * (10 ** -8)
        omega_y_sk_pz = 1.679685 * (10 ** -6)
        omega_z_sk_pz = 3.850439 * (10 ** -6)
        
        # ----- Матрицы СК-42 в ПЗ-90.11 ----- #
        self.koef_sk_pz = 1 + ((-0.228) * (10 ** -6))
        
        self.matrix_trans_sk_pz = np.array([[1, - omega_z_sk_pz, omega_y_sk_pz],
                                            [omega_z_sk_pz, 1, - omega_x_sk_pz],
                                            [- omega_y_sk_pz, omega_x_sk_pz, 1]], dtype=np.float64)
        
        self.matrix_delta_sk_pz = np.array([[23.557],
                                            [-140.844],
                                            [-79.778]], dtype=np.float64)
        
        # ----- Матрицы ПЗ-90.11 в СК-42 ----- #
        self.koef_pz_sk = 1 - ((-0.228) * (10 ** -6))
        
        self.matrix_trans_pz_sk = np.array([[1, omega_z_sk_pz, - omega_y_sk_pz],
                                            [- omega_z_sk_pz, 1, omega_x_sk_pz],
                                            [omega_y_sk_pz, - omega_x_sk_pz, 1]], dtype=np.float64)
        
        self.matrix_delta_pz_sk = - np.array([[23.557],
                                              [-140.844],
                                              [-79.778]], dtype=np.float64)
        
        # ----- Угловые коэффициенты WGS-84, ПЗ-90.11 ----- #
        omega_x_wgs_pz = 1.115071 * (10 ** -8)
        omega_y_wgs_pz = 1.716240 * (10 ** -8)
        omega_z_wgs_pz = 2.041066 * (10 ** -8)
        
        # ----- Матрицы WGS-84 в ПЗ-90.11 ----- #
        self.koef_wgs_pz = 1 + ((-0.008) * (10 ** -6))
        
        self.matrix_trans_wgs_pz = np.array([[1, - omega_z_wgs_pz, - omega_y_wgs_pz],
                                             [omega_z_wgs_pz, 1, - omega_x_wgs_pz],
                                             [omega_y_wgs_pz, omega_x_wgs_pz, 1]], dtype=np.float64)
        
        self.matrix_delta_wgs_pz = np.array([[-0.013],
                                             [0.106],
                                             [0.022]], dtype=np.float64)
        
        # ----- Матрицы ПЗ-90.11 в WGS-84 ----- #
        self.koef_pz_wgs = 1 - ((-0.008) * (10 ** -6))
        
        self.matrix_trans_pz_wgs = np.array([[1, omega_z_wgs_pz, omega_y_wgs_pz],
                                             [- omega_z_wgs_pz, 1, omega_x_wgs_pz],
                                             [- omega_y_wgs_pz, - omega_x_wgs_pz, 1]], dtype=np.float64)
        
        self.matrix_delta_pz_wgs = - np.array([[-0.013],
                                               [0.106],
                                               [0.022]], dtype=np.float64)
    
    
    
    @staticmethod
    def trans_to_dd_ddd(lat_grad, long_grad):
        """
        Функция преобразования формата (DD, MM, SS) в (DD.DDDD)
        в радианах
        """
        if len(lat_grad) == 3:
            if lat_grad[0] > 0:
                lat_rad = round(((lat_grad[2] / 3600) + (lat_grad[1] / 60) + lat_grad[0]) * (pi / 180), 8)
            else:
                lat_rad = round((- (lat_grad[2] / 3600) - (lat_grad[1] / 60) + lat_grad[0]) * (pi / 180), 8)
        else:
            lat_rad = round(lat_grad[0] * (pi/180), 8)
            
        if len(long_grad) == 3:
            if long_grad[0] > 0:
                long_rad = round(((long_grad[2] / 3600) + (long_grad[1] / 60) + long_grad[0]) * (pi / 180), 8)
            else:
                long_rad = round((- (long_grad[2] / 3600) - (long_grad[1] / 60) + long_grad[0]) * (pi / 180), 8)
        else:
            long_rad = round(long_grad[0] * (pi/180), 8)
            
        return lat_rad, long_rad
        
    
    @staticmethod
    def trans_to_ddmmss(lat_rad, long_rad):
        """[summary]

        Функция преобразования формата (DD.DDDD) в (DD, MM, SS)
        в градусах
        """
        lat_rad = round(lat_rad, 8)
        long_rad = round(long_rad, 8)
        
        lat_grad = lat_rad * (180 / pi)
        long_grad = long_rad * (180 / pi)
        
        lat_dd = trunc(lat_grad)
        lat_mm = trunc((lat_grad - lat_dd) * 60)
        lat_ss = round(((lat_grad - lat_dd) * 60 - lat_mm) * 60, 2)
        
        long_dd = trunc(long_grad)
        long_mm = trunc((long_grad - long_dd) * 60)
        long_ss = round(((long_grad - long_dd) * 60 - long_mm) * 60, 2)
        
        lat_out = [lat_dd, lat_mm, lat_ss]
        long_out = [long_dd, long_mm, long_ss]
        
        return lat_out, long_out
        
        
    @staticmethod
    def __data_checking(lat_grad: list, long_grad: list, alt_m: float):
        """Функция проверки введенных данных

        Args:
            lat_grad (list): 
                Широта точки в градусах может быть в формате (DD, MM, SS) или в (DD.DDDD)
            long_grad (list): 
                Долгота точки в градусах может быть в формате (DD, MM, SS) или в (DD.DDDD)
            alt_m (float): 
                Высота точки в метрах
        """
        if len(lat_grad) == 3:
            lat_grad_dd = lat_grad[0] + lat_grad[1]/60 + lat_grad[2]/3600
            if lat_grad_dd > 90 or lat_grad_dd < -90:
                raise ValueError('Максимальный диапазон широты +/- 90 градусов')
        else:
            if lat_grad[0] > 90 or lat_grad[0] < -90:
                raise ValueError('Максимальный диапазон широты +/- 90 градусов')
            
        if len(long_grad) == 3:
            long_grad_dd = long_grad[0] + long_grad[1]/60 + long_grad[2]/3600
            if long_grad_dd > 180 or long_grad_dd < -180:
                raise ValueError('Максимальный диапазон долготы +/- 180 градусов')
        else:
            if long_grad[0] > 180 or long_grad[0] < -180:
                raise ValueError('Максимальный диапазон долготы +/- 180 градусов')
            
        if alt_m > 10000:
            raise ValueError('Максимальная высота точки 10"000 метров')
    
    
    def georaph_to_rsc(self, sk_type: str) -> tuple:
        """Функция преобразования географических координат в
        прямоугольные пространственные координаты

        Args:
            sk_type (str): 
                'PZ' - с/к ПЗ-90.11
                'WGS' - с/к WGS-84
                'SK' - с/к СК-42

        Returns:
            tuple: (x_rsc, y_rsc, z_rsc)
        """
        self.__data_checking(self.lat_grad, self.long_grad, self.alt_m)
        
        lat_geograph, long_geograph = self.trans_to_dd_ddd(self.lat_grad, self.long_grad)
        alt_geograph = self.alt_m
        
        if sk_type == 'PZ':
            a_ellips = self.a_ellips_pz[0]
            eccent_ellips = self.square_eccent_pz[0]
            
            n_ellips = a_ellips / sqrt(1 - eccent_ellips * (sin(lat_geograph) ** 2))
        
        elif sk_type == 'WGS':
            a_ellips = self.a_ellips_wgs[0]
            eccent_ellips = self.square_eccent_wgs[0]
            
            n_ellips = a_ellips / sqrt(1 - eccent_ellips * (sin(lat_geograph) ** 2))
        
        elif sk_type == 'SK':
            a_ellips = self.a_ellips_sk[0]
            eccent_ellips = self.square_eccent_sk[0]
            
            n_ellips = a_ellips / sqrt(1 - eccent_ellips * (sin(lat_geograph) ** 2))
        
        x_rsc = (n_ellips + self.alt_m) * cos(lat_geograph) * cos(long_geograph)
        y_rsc = (n_ellips + self.alt_m) * cos(lat_geograph) * sin(long_geograph)
        z_rsc = (((1 - eccent_ellips) * n_ellips) + alt_geograph) * sin(lat_geograph)
        
        return x_rsc, y_rsc, z_rsc
    
    
    def rsc_to_geograph(self, x_rsc: float, y_rsc: float, z_rsc: float, sk_type: str, stop_val: float = 0.000001) -> tuple:
        """Функция преобразования прямоугольных пространственных координат в
        географиеские координаты

        Args:
            x_rsc (float): положение по Х, рад.
            y_rsc (float): положение по Y, рад.
            z_rsc (float): положение по Z, рад.
            sk_type (str): 
                'PZ' - с/к ПЗ-90.11
                'WGS' - с/к WGS-84
                'SK' - с/к СК-42

        Returns:
            tuple: (lat рад., long рад., alt м.)
        """
        d_help = sqrt((x_rsc ** 2) + (y_rsc ** 2))
        l_a_help = abs(asin(y_rsc / d_help))
        
        if sk_type == 'PZ':
            a_ellips = self.a_ellips_pz[0]
            eccent_ellips = self.square_eccent_pz[0]
            
        elif sk_type == 'WGS':
            a_ellips = self.a_ellips_wgs[0]
            eccent_ellips = self.square_eccent_wgs[0]
            
        elif sk_type == 'SK':
            a_ellips = self.a_ellips_sk[0]
            eccent_ellips = self.square_eccent_sk[0]
            
        if d_help == 0:
            lat_geograph = (pi/2) * (z_rsc/abs(z_rsc))
            long_geograph = 0
            alt_geograph = z_rsc * sin(lat_geograph) - a_ellips * sqrt(1 - eccent_ellips * (sin(lat_geograph) ** 2))
        
        else:
           if y_rsc < 0 and x_rsc > 0:
               long_geograph = (2 * pi) - l_a_help  
            
           elif y_rsc < 0 and x_rsc < 0:
               long_geograph = pi + l_a_help
               
           elif y_rsc > 0 and x_rsc < 0:
               long_geograph = pi - l_a_help
               
           elif y_rsc > 0 and x_rsc > 0:
               long_geograph = l_a_help
               
           elif y_rsc == 0 and x_rsc > 0:
               long_geograph = 0
               
           elif y_rsc == 0 and x_rsc < 0:
               long_geograph = pi
               
           if z_rsc == 0:
               lat_geograph = 0
               alt_geograph = d_help - a_ellips
           else:
               r_help = sqrt((x_rsc ** 2) + (y_rsc ** 2) + (z_rsc ** 2))
               c_help = asin(z_rsc / r_help)
               p_help = (eccent_ellips * a_ellips) / (2 * r_help)
               s1_help = 0
            
               while True:
                   b_help = c_help + s1_help
                   s2_help = asin((p_help * sin(2 * b_help)) / sqrt(1 - eccent_ellips * (sin(b_help) ** 2)))
                   d_stop = abs(s2_help - s1_help)
                   if d_stop <= stop_val:
                       break
                   else:
                       s1_help = s2_help
                   
               lat_geograph = b_help
               alt_geograph = round(d_help * cos(lat_geograph) + z_rsc * sin(lat_geograph) - a_ellips * sqrt(1 - eccent_ellips * (sin(lat_geograph) ** 2)), 0)
        
               return lat_geograph, long_geograph, alt_geograph
        
        
    def pz_to_wgs(self) -> tuple:
        """Функция преобразования ПЗ-90.11 в WGS-84

        Returns:
            tuple: (lat_geograph град., long_geograph град., alt_geograph м.)
        """
        x_rsc, y_rsc, z_rsc = self.georaph_to_rsc(sk_type='PZ')
    
        coord_pz = np.array([[x_rsc],
                             [y_rsc],
                             [z_rsc]], dtype=np.float64)
            
        coord_wgs = np.dot(np.dot(self.koef_pz_wgs, self.matrix_trans_pz_wgs), coord_pz) + self.matrix_delta_pz_wgs
        
        x_rsc_out = coord_wgs[0][0]
        y_rsc_out = coord_wgs[1][0]
        z_rsc_out = coord_wgs[2][0]
    
        lat_geograph, long_geograph, alt_geograph = self.rsc_to_geograph(x_rsc_out, y_rsc_out, z_rsc_out, sk_type='WGS')
        lat_grad, long_grad = self.trans_to_ddmmss(lat_geograph, long_geograph)
        
        return lat_grad, long_grad, alt_geograph
    
    
    def wgs_to_pz(self) -> tuple:
        """Функция преобразования WGS-84 в ПЗ-90.11

        Returns:
            tuple: (lat_geograph град., long_geograph град., alt_geograph м.)
        """
        x_rsc, y_rsc, z_rsc = self.georaph_to_rsc(sk_type='WGS')
    
        coord_wgs = np.array([[x_rsc],
                              [y_rsc],
                              [z_rsc]], dtype=np.float64)
            
        coord_pz = np.dot(np.dot(self.koef_wgs_pz, self.matrix_trans_wgs_pz), coord_wgs) + self.matrix_delta_wgs_pz
        
        x_rsc_out = coord_pz[0][0]
        y_rsc_out = coord_pz[1][0]
        z_rsc_out = coord_pz[2][0]
    
        lat_geograph, long_geograph, alt_geograph = self.rsc_to_geograph(x_rsc_out, y_rsc_out, z_rsc_out, sk_type='PZ')
        lat_grad, long_grad = self.trans_to_ddmmss(lat_geograph, long_geograph)
        
        return lat_grad, long_grad, alt_geograph
        
    
    def pz_to_sk(self) -> tuple:
        """Функция преобразования ПЗ-90.11 в СК-42

        Returns:
            tuple: (lat_geograph град., long_geograph град., alt_geograph м.)
        """  
        x_rsc, y_rsc, z_rsc = self.georaph_to_rsc(sk_type='PZ')
    
        coord_pz = np.array([[x_rsc],
                             [y_rsc],
                             [z_rsc]], dtype=np.float64)
            
        coord_sk = np.dot(np.dot(self.koef_pz_sk, self.matrix_trans_pz_sk), coord_pz) + self.matrix_delta_pz_sk
        
        x_rsc_out = coord_sk[0][0]
        y_rsc_out = coord_sk[1][0]
        z_rsc_out = coord_sk[2][0]
    
        lat_geograph, long_geograph, alt_geograph = self.rsc_to_geograph(x_rsc_out, y_rsc_out, z_rsc_out, sk_type='SK')
        lat_grad, long_grad = self.trans_to_ddmmss(lat_geograph, long_geograph)
        
        return lat_grad, long_grad, alt_geograph
    
    
    def sk_to_pz(self) -> tuple:
        """Функция преобразования СК-42 в ПЗ-90.11

        Returns:
            tuple: (lat_geograph град., long_geograph град., alt_geograph м.)
        """  
        x_rsc, y_rsc, z_rsc = self.georaph_to_rsc(sk_type='SK')
    
        coord_sk = np.array([[x_rsc],
                             [y_rsc],
                             [z_rsc]], dtype=np.float64)
            
        coord_pz = np.dot(np.dot(self.koef_sk_pz, self.matrix_trans_sk_pz), coord_sk) + self.matrix_delta_sk_pz
        
        x_rsc_out = coord_pz[0][0]
        y_rsc_out = coord_pz[1][0]
        z_rsc_out = coord_pz[2][0]
    
        lat_geograph, long_geograph, alt_geograph = self.rsc_to_geograph(x_rsc_out, y_rsc_out, z_rsc_out, sk_type='PZ')
        lat_grad, long_grad = self.trans_to_ddmmss(lat_geograph, long_geograph)
        
        return lat_grad, long_grad, alt_geograph
    
    
    def wgs_to_sk(self) -> tuple:
        """Функция преобразования WGS-84 в СК-42

        Returns:
            tuple: (lat_geograph град., long_geograph град., alt_geograph м.)
        """ 
        self.lat_grad, self.long_grad, self.alt_m = self.wgs_to_pz()
        lat_grad, long_grad, alt_geo = self.pz_to_sk()
        
        return lat_grad, long_grad, alt_geo
    
    
    def sk_to_wgs(self) -> tuple:
        """Функция преобразования СК-42 в WGS-84

        Returns:
            tuple: (lat_geograph град., long_geograph град., alt_geograph м.)
        """ 
        self.lat_grad, self.long_grad, self.alt_m = self.sk_to_pz()
        lat_grad, long_grad, alt_geo = self.pz_to_wgs()
        
        return lat_grad, long_grad, alt_geo
    
   
