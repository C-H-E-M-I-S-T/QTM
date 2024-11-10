from Electrostatics import Electrostatics
from Combiner import execute
from Processes import IonsInteraction
from Processes import Molecules_Interaction
class Kinetics(Electrostatics):
    def __init__(self, P, T, Cv, Cp, V, m):
#Постоянные
#постоянная больцмана
        self.k = float(1.380649e-23)
#постоянная авогадро с перерасчетом на допустимое число операций
       # self.operations_limit = 536870912
        self.operations_limit = 10000
        self.Na = float(6.022045e23)*(self.operations_limit/6.022045e23)
        self.Na_ = float(6.022045e23)
#Интенсивные параметры термодинамической системы
        self.P,self.T,self.Cv,self.Cp = (1.0, 1.0, 1.0, 1.0)
#Экстенсивные параметры т.с. (объм выражен в кубических метрах)
        self.m,self.V = (m,V)
#Функции термодинамических процессов
        #self.Q, self.A = (1.0, 1.0)
    #Блок по химической кинетики для
    # дальнйшей сборки массива из взаимодействий

    #Общие выражения для реакций второго порядка
    t_ = 0 #секунд
    t_1 = -1 # время начала образования осадка (данные эксперимента)
    t_2 = -1 # время достижения равновесия

    c_0 = 0.5 # концентрация ионов лития, фтора, численно равные конечной 
    # концентрации фторида лития

    k = 1 #данные эксперимента, константа скорости
    av_speed = 0.01 #средняя скорость до начала образования осадка
    def concentrations(self,t):
        #Функция возвращает кортеж: концентрации ионов лития и фтора ; концентрация фторида лития
        if t_1==-1 or t_2==-1:
            t_1=((0.0038)**0.5)/(self.k*(self.c_0-self.av_speed*t)**2)
            t_2=(self.c_0-(0.0038)**0.5)/(self.k*(self.c_0-self.av_speed*t)**2)
        if (t<t_1 and t_1>0) or (t_1 == -1 and t < ((0.0038)**0.5))/(self.k*(self.c_0-self.av_speed*t)**2):
            return ((self.c_0-t*(self.k*(self.c_0-self.av_speed*t)**2)),t*(self.k*(self.c_0-self.av_speed*t)**2))
        if (t_1<t<t_2 and t_2>0) or (t_2==-1 and t_1<t<(self.c_0-(0.0038)**0.5)/(self.k*(self.c_0-self.av_speed*t)**2)):
            return ((self.c_0-t*(self.k*(self.c_0-self.av_speed*t)**2)),(0.0038)**0.5)
        if (t>t_2 and t_2>0) or (t_2==-1 and t>(self.c_0-(0.0038)**0.5)/(self.k*(self.c_0-self.av_speed*t)**2)):
            return ((0.0038)**0.5, (0.0038)**0.5)   
    def reaction_speed(self,t):
        return t*(self.k*(self.c_0-self.av_speed*t)**2)
                
    top_electrostatic_energy = 2.8780165475731916883e-20
    def pl_ions_interaction(system,calculation=False)->float:
        if (calculation==True):
            return system.modulate_ionic_bonding(10e-30,93.75e-12)    
        else: 
            _dQ=float(-1*system.top_electrostatic_energy)
            _dS=_dQ/system.T
            return [system.dH(_dQ, system.P, system.P, system.V, system.V), system.dG(_dQ, system.T,system.T,_dS+system.S,system.S,system.P,system.P,system.V,system.V), _dS]
    def ds_ions_interaction(system,calculation=False)->float:
        if (calculation==True):
            return system.modulate_equal_charges_repulsion(10e-30,93.75e-12)
        else: 
            _dQ=float(system.top_electrostatic_energy)
            _dS=_dQ/system.T
            return [system.dH(_dQ, system.P, system.P, system.V, system.V), system.dG(_dQ, system.T,system.T,_dS+system.S,system.S,system.P,system.P,system.V,system.V), _dS]
    
    def modulate(self,t0,t1):
        integral_values = [0.0,0.0,0.0]
        for t in range(t0,t1,1):
            print(f"Модулирую {t} секунду")
            #limit = self.above_than(self.reaction_speed(t)*self.Na,self.operations_limit)
            limit = self.operations_limit
            processes = [(self.pl_ions_interaction,self.ds_ions_interaction) for i in range(int((limit/2)))]
            result = execute(self,processes)
            integral_values += [result[0]*(self.Na_/limit), result[1]*(self.Na_/limit),result[2]*(self.Na_/limit)]
        print(f"В результате реакции с {t0} по {t1} секунд состояние изменилось следующим образом: {integral_values}")


test_system = Kinetics(101325, 298.15, 1,1,10e-6,100)
test_system.modulate(1,5)
            