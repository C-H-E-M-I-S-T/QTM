import math
from random import randint


class Thermodynamics:
    def __init__(self, P, T, Cv, Cp, V, m):
#Постоянные
#постоянная больцмана
        self.k = float(1.380649e-23)
#постоянная авогадро с перерасчетом на допустимое число операций
        self.operations_limit = 536870912
        self.Na = float(6.022045e23)*(self.operations_limit/6.022045e23)
#Интенсивные параметры термодинамической системы
        self.P,self.T,self.Cv,self.Cp = (1.0, 1.0, 1.0, 1.0)
#Экстенсивные параметры т.с. (объм выражен в кубических метрах)
        self.m,self.V = (m,V)
#Функции термодинамических процессов
        #self.Q, self.A = (1.0, 1.0)

#Макроскопические функции состояния системы

#классическое определение работы системы как работы
#по расширению газа
    def dA(self,P, V2, V1)->float:
        return P*(V2-V1)

#определение свободной внутренней энергии системы
    def U(self,Q, A)->float:
        return float(Q-A)
    def dU(self,dQ,dA)->float:
        return float(dQ-dA)

# определение изменения теплоты через изохорную
# теплоемкость и изотермическую теплоту расширения
# тела (при постоянном давлении)
    def dQ(self,l, V2, V1, Cv, T2, T1)->float:
        return float (l*(V2-V1)+Cv*(T2-T1))

# определение изменения свободной внутренней энергии
# через изотермическую теплоту расширения тела и 
# изохорную теплоёмкость
    def dU(self,l, P, V2, V1, Cv, T2, T1)->float:
        return float((l-P)*(V2-V1)+Cv*(T2-T1))
# (аналогичного определения для изобарической 
# теплоемкости нет, так как при непостоянстве объёма 
# имеет место символьно неопределенная работа по 
# расширению тела)

# определение изменения теплоты через изобарную 
# теплоемкость и изотермическую теплоту возрастания 
# давления
    def dQ(self,h, P2, P1, Cp, T2, T1)->float:
        return float(h*(P2-P1)+Cp*(T2-T1))

#определение энтальпии
    def H(self,U,P,V)->float:
        return float(U+P*V)
    def dH(self,dU,P1,P2,V1,V2)->float:
        return float(dU+(P2-P1)*V2+(V2-V1)*P2)

#определение 
    def F(self,U,T,S)->float:
        return float(U-T*S)
    def dF(dU, T2, T1, S2, S1):
        return float(dU-(T2-T1)*S2-(S2-S1)*T2)

#определение изобарно-изотермического потенциала или
#функции Гиббса
    def G(self,U,T,S,P,V)->float:
        return float(U-T*S+P*V)
    def dG(self,dU, T2, T1, S2, S1, P2, P1, V2, V1):
        return float(dU-(T2-T1)*S2-(S2-S1)*T2+
                 (P2-P1)*V2+(V2-V1)*P2)
    def dG_isobaric_isotermic(self,dH, T, dS):
        return float(dH-T*dS)
    
#классическое определение энтропии как приведенной теплоты
    S = float(0)
    def dS(self,dQ,T)->float:
        if (T>0):
            return dQ/T
        else: return 0

#Определения химического потенциала
#(принимаются равными при равновесных условиях)
    def nu_i(self,component, U, Q, G, F, H, n):
        if n>0:
            return (component, U/n, Q/n, G/n, F/n, H/n)
    
# минимальный объём термодинамической системы определен 
# наименьшим допустимым числом микрочастиц (10^18) и равен 10^-9 л, 
# т.е. одной миллионной миллилитра; при этом минимальное число 
# различимых микросостояний определяется вероятностью самопроизвольного 
# теплоового возбуждения молекул
    def minimal_microstates_amount(self,T, dE):
        return 10**18*math.e**(-(dE/(k*T)))


#Нужно по идеи создать функцию, которая 
# в качестве аргуемента принимала бы кортеж микроскопических 
# процессов как функций и обрабатывала бы суммарное изменение 
# системы в результате их совокупного протекания

# В свою очередь, что должны возвращать функции микросостояния? 
# Могут ли они сами возвращать кортеж из основных 
# термодинамических функций, например изменения энтропии, 
# энтальпии и изобарно-изотермического потенциала?

#Примеры различных функций
    def preset_microprocess_ex(self,calculate=False):
        if calculate == True:
            print("calculating preset microinteraction...")
            return [1.4, -0.6, -0.9]  #calculated result
        else: return[2.0, -1.0, -1.1] #appoximate result
    def preset_disbalancing_microprocess_ex(self,calculate=False):
        if calculate == True:
            print("calculating disbalancing microinteraction...")
            return [1.4, 0.4, 0.5]
        else: return[2.0, 0.96, 0.94]
    def preset_disbalancing_microprocess_ex_2(self,calculate=False):
        if calculate == True:
            print("calculating disbalancing microinteraction...")
            return [1.4, 0.4, 0.5]
        else: return[2.0, 1.96, 1.94]
    def above_than(self,*args):
        print(args)
        v,b=args
        if (v>b):
            return v
        else: return b
    #Сила диффузии молекул в воде при н.у.
    Niu_Diffusion = (8.9*10e-4)*10**-9
    