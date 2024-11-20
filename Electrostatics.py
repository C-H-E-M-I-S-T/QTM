import math
from random import randint
from Macroscopic import Thermodynamics
class Electrostatics(Thermodynamics):
    def __init__(self, P, T, Cv, Cp, V, m):
#Постоянные
#постоянная больцмана
        self.k = float(1.380649e-23)
#постоянная авогадро с перерасчетом на допустимое число операций
        self.operations_limit = 536870912
        
        self.Na = float(6.022045e23)*(self.operations_limit/6.022045e23)
        self.Na_ = float(6.022045e23)
#Интенсивные параметры термодинамической системы
        self.P,self.T,self.Cv,self.Cp = (1.0, 1.0, 1.0, 1.0)
#Экстенсивные параметры т.с. (объм выражен в кубических метрах)
        self.m,self.V = (m,V)
#Функции термодинамических процессов
        #self.Q, self.A = (1.0, 1.0)
        
#Блок-Физикохимическое моделирование взаимодействий между ионами
    def Niu_coulomb_interaction_Li_F(self,dist=float, FORCE=bool):
        power = float(0)
        #УЧТЕНА ДИЭЛЕКТРИЧЕСКАЯ ПРОНИЦАЕМОСТЬ ДИСТИЛЛИРОВАННОЙ ВОДЫ
        diel=80.2
        e_z=8.85*10e-12
        #Формула взята не из физ химии Герасимова, но из источника по электростатике
        print(f"{-259.714} - кдж/моль по справочным данным")
        print(f"{((-1)*(1.602176634e-19)**2/(diel*dist))*(self.Na_)} - кдж/моль расчет без K")
        print(f"{((-1)*(1.602176634e-19)**2/(diel*dist*4*math.pi*e_z))*(self.Na_)} - кдж/моль расчет с K")
        print(f"сила составила {9*10e9*(1.602176634e-19)**2/((dist**2)*diel)}")
        if (FORCE==True):
            return (9*10e9*(1.602176634e-19)**2/(dist**2*diel))
        else: return((-1)*(1.602176634e-19)**2/(diel*dist*4*math.pi*e_z))
    
    #Чтобы соответствовать этим понятиям, недостаточно лишь электростатических расчетов, 
    #но образование кристаллической решетки ионного строения
    def energy_release_from_ions_bonding(mole=float):
        #Термохимически рассчитанная величина энергии спаривания Li+ и F-
        return [-259.714*mole, -584.1*mole,35.9*mole]
    def part_energy_release_from_ions_bonding(self):
        #Термохимически рассчитанная величина энергии спаривания 
        #буквально двух частиц, Li+ и F-
        return [-259.714*(1/(self.Na)), -584.1*(1/(self.Na)),35.9*(1/(self.Na))]
    
    def modulate_ionic_bonding(self,volume=float, particles_radius=float):
        #согласно расчетам протяженность вдоль координатных осей около 1,1;
        #принимаю равной 1.
        Li = [randint(0,1) for i in range(3)]
        F = [randint(0,1) for i in range(3)]
        distance = self.above_than((((((Li[0]-F[0])**2+(Li[1]-F[1]))**2+(Li[2]-F[2])**2)**0.5)*10e-10),2*93.75e-12)
        diffusion_force = (8.9*10e-4)*10**-9
        if self.Niu_coulomb_interaction_Li_F(distance, True)>diffusion_force:
            _dQ=self.Niu_coulomb_interaction_Li_F(distance, False)
            _dS=_dQ/self.T
            return [self.dH(_dQ, self.P, self.P, self.V, self.V), self.dG_isobaric_isotermic(_dQ,self.T,_dS), _dS]
        else: return[0.0,0.0,0.0]
    def modulate_dissociaition(self,calculation=False):
            return [259.714*(1/self.Na), 584.1*(1/self.Na),-35.9*(1/self.Na)]
    def modulate_equal_charges_repulsion(self,volume=float, particles_radius=float):
        #согласно расчетам протяженность вдоль координатных осей около 1,1;
        #принимаю равной 1.
        Li = [randint(0,1) for i in range(3)]
        F = [randint(0,1) for i in range(3)]
        distance = self.above_than((((((Li[0]-F[0])**2+(Li[1]-F[1]))**2+(Li[2]-F[2])**2)**0.5)*10e-10),2*93.75e-12)
        if self.Niu_coulomb_interaction_Li_F(distance,True) > self.Niu_Diffusion:
            _dQ=float(abs(self.Niu_coulomb_interaction_Li_F(distance, False)))            
            _dS=_dQ/self.T
            #При постоянстве давления и объёма макроскопическая 
            # работа по расширению системы равна нулю, след-но 
            # изменение св.вн.эн. равно изменению теплоты
            return [self.dH(_dQ, self.P, self.P, self.V, self.V), self.dG(_dQ, self.T,self.T,_dS+self.S,self.S,self.P,self.P,self.V,self.V), _dS]
        
        else: return [0.0,0.0,0.0]

   
    
#Блок-термодинамическое моделирование взаимодействия между 
# ионами и молекулами

    def Niu_coulomb_interaction_LiF_ion(self,dist=float, mol_center = [],li__or_f = [], FORCE=bool):
        power = float(0)
        #!!! - УЧТЕНА ДИЭЛЕКТРИЧЕСКАЯ ПРОНИЦАЕМОСТЬ ДИСТИЛЛИРОВАННОЙ ВОДЫ
        LiF_bond_dipol = 2.110593454703e-29 #Кл*М
        diel=80.2
        e_z=8.85*10e-12
        if (((mol_center[0]**2+mol_center[1]**2+mol_center[2]**2)**0.5)*((li__or_f[0]**2+li__or_f[1]**2+li__or_f[2]**2)**0.5)) > 0:
            cos_alpha = (mol_center[0]*li__or_f[0]+mol_center[1]*li__or_f[1]+mol_center[2]*li__or_f[2])/(((mol_center[0]**2+mol_center[1]**2+mol_center[2]**2)**0.5)*((li__or_f[0]**2+li__or_f[1]**2+li__or_f[2]**2)**0.5))
        else: cos_alpha = 0
        #Формула взята не из физ химии Герасимова, но из источника по электростатике
        if (FORCE==True):
            return (1/(4*math.pi*8.85418782e-12))*((1.602176634e-19)*LiF_bond_dipol)/(diel*(dist**3))*((3*cos_alpha)+1)**0.5
        else: return(1/(4*math.pi*8.85418782e-12))*((1.602176634e-19)*LiF_bond_dipol)/(diel*(dist**2))*cos_alpha
    def modulate_ions_to_molecule_adsorbtion(self):
        Li_or_F = [randint(0,1) for i in range(3)]
        LiF_center = [randint(0,1) for i in range(3)]
        distance = self.above_than((((((LiF_center[0]-Li_or_F[0])**2+(LiF_center[1]-Li_or_F[1]))**2+(LiF_center[2]-Li_or_F[2])**2)**0.5)*10e-10),93.75e-12)
        print(f"LiF coordinates {LiF_center} ion coordinates {Li_or_F}")
        print(f"Расстояние между молекулой и ионом {distance} метра")
        print(f"сила диффузии {self.Niu_Diffusion}")
        print(f"кулоновская сила притяжения иона к молекуле {self.Niu_coulomb_interaction_LiF_ion(distance, LiF_center, Li_or_F, True)}")
        if self.Niu_coulomb_interaction_LiF_ion(distance, LiF_center, Li_or_F, True)>self.Niu_Diffusion:
            _dQ=self.Niu_coulomb_interaction_LiF_ion(distance, LiF_center, Li_or_F, False)
            _dS=_dQ/self.T
            #При постоянстве давления и объёма макроскопическая 
            # работа по расширению системы равна нулю, след-но 
            # изменение св.вн.эн. равно изменению теплоты
            return [self.dH(_dQ, self.P, self.P, self.V, self.V), self.dG_isobaric_isotermic(_dQ,self.T,_dS), _dS]
        else: return [0.0,0.0,0.0]


    
