from math import pi, sin, sqrt, cos

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
        self.lat_grad = lat_grad
        self.long_grad = long_grad
        self.alt_m = alt_m
        
    @staticmethod    
    def __trans_to_dd_ddd(lat_grad: list, long_grad: list) -> tuple:
        """Функция преобразования формата (DD, MM, SS) в (DD.DDDD)

        Args:
            lat_grad (list): 
                Широта точки в градусах может быть в формате (DD, MM, SS) или в (DD.DDDD)
            long_grad (list): 
                Долгота точки в градусах может быть в формате (DD, MM, SS) или в (DD.DDDD)

        Returns:
            tuple: (Широта в радианах в формате DD.DDDD, Долгота в радианах в формате DD.DDD)
        """
        if len(lat_grad) == 3:
            lat_grad = lat_grad[0] + lat_grad[1]/60 + lat_grad[2]/3600
            lat_grad = lat_grad * (pi/180)
        else:
            lat_grad = lat_grad[0] * (pi/180)
            
        if len(long_grad) == 3:
            long_grad = long_grad[0] + long_grad[1]/60 + long_grad[2]/3600
            long_grad = long_grad * (pi/180)
        else:
            long_grad = long_grad[0] * (pi/180)
            
        return lat_grad, long_grad
    
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
        
        
    def __geograph_to_rsc(self, a_ellips: int, a_compr: float) -> tuple:
        """Пересчет координат из географических координат в прямоугольные 
        пространственные координаты

        Args:
            a_ellips (int): Большая полуось эллипсоида в выбранной с/к, метры
            a_compr (float): Коэффициент сжатия эллипсоида

        Returns:
            tuple: (x, y, z)
        """
        eccent_ellips = (2 * a_compr) - (a_compr ** 2)
        r_curv_first_vert = a_ellips / (sqrt(1 - (eccent_ellips * ((sin(self.lat_grad)) ** 2))))
        
        x_rsc = (r_curv_first_vert + self.alt_m) * cos(self.lat_grad) * cos(self.long_grad)
        y_rsc = (r_curv_first_vert + self.alt_m) * cos(self.lat_grad) * sin(self.long_grad)
        z_rsc = ((1 - eccent_ellips) * r_curv_first_vert + self.alt_m )* sin(self.lat_grad)
        
        return x_rsc, y_rsc, z_rsc
        
