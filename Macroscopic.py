import math
from random import randint

#Постоянные

#постоянная больцмана
k = float(1.380649e-23)
#постоянная авогадро
Na = float(6.022045e23)
#Интенсивные параметры термодинамической системы
P,T,Cv,Cp = (1.0, 1.0, 1.0, 1.0)

#Экстенсивные параметры т.с. (объм выражен в кубических метрах)
m,V = (float(1),float(10-4))

#Функции термодинамических процессов
Q, A = (1.0, 1.0)

#Макроскопические функции состояния системы

#классическое определение работы системы как работы
#по расширению газа
def dA(P, V2, V1)->float:
    return P*(V2-V1)

#определение свободной внутренней энергии системы
def U(Q, A)->float:
    return float(Q-A)
def dU(dQ,dA)->float:
    return float(dQ-dA)

# определение изменения теплоты через изохорную
# теплоемкость и изотермическую теплоту расширения
# тела (при постоянном давлении)
def dQ(l, V2, V1, Cv, T2, T1)->float:
    return float (l*(V2-V1)+Cv*(T2-T1))

# определение изменения свободной внутренней энергии
# через изотермическую теплоту расширения тела и 
# изохорную теплоёмкость
def dU(l, P, V2, V1, Cv, T2, T1)->float:
    return float((l-P)*(V2-V1)+Cv*(T2-T1))
# (аналогичного определения для изобарической 
# теплоемкости нет, так как при непостоянстве объёма 
# имеет место символьно неопределенная работа по 
# расширению тела)

# определение изменения теплоты через изобарную 
# теплоемкость и изотермическую теплоту возрастания 
# давления
def dQ(h, P2, P1, Cp, T2, T1)->float:
    return float(h*(P2-P1)+Cp*(T2-T1))

#определение энтальпии
def H(U,P,V)->float:
    return float(U+P*V)
def dH(dU,P1,P2,V1,V2)->float:
    return float(dU+(P2-P1)*V2+(V2-V1)*P2)

#определение 
def F(U,T,S)->float:
    return float(U-T*S)
def dF(dU, T2, T1, S2, S1):
    return float(dU-(T2-T1)*S2-(S2-S1)*T2)

#определение изобарно-изотермического потенциала или
#функции Гиббса
def G(U,T,S,P,V)->float:
    return float(U-T*S+P*V)
def dG(dU, T2, T1, S2, S1, P2, P1, V2, V1):
    return float(dU-(T2-T1)*S2-(S2-S1)*T2+
                 (P2-P1)*V2+(V2-V1)*P2)
def dG_isobaric_isotermic(dH, T, dS):
    return float(dH-T*dS)
    
#классическое определение энтропии как приведенной теплоты
S = float(0)
def dS(dQ,T)->float:
    if (T>0):
        return dQ/T
    else: return 0

#Определения химического потенциала
#(принимаются равными при равновесных условиях)
def nu_i(component, U, Q, G, F, H, n):
    if n>0:
        return (component, U/n, Q/n, G/n, F/n, H/n)
    
# минимальный объём термодинамической системы определен 
# наименьшим допустимым числом микрочастиц (10^18) и равен 10^-9 л, 
# т.е. одной миллионной миллилитра; при этом минимальное число 
# различимых микросостояний определяется вероятностью самопроизвольного 
# теплоового возбуждения молекул
def minimal_microstates_amount(T, dE):
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
def preset_microprocess_ex(calculate=False):
    if calculate == True:
        print("calculating preset microinteraction...")
        return [1.4, -0.6, -0.9]  #calculated result
    else: return[2.0, -1.0, -1.1] #appoximate result
def preset_disbalancing_microprocess_ex(calculate=False):
    if calculate == True:
        print("calculating disbalancing microinteraction...")
        return [1.4, 0.4, 0.5]
    else: return[2.0, 0.96, 0.94]
def preset_disbalancing_microprocess_ex_2(calculate=False):
    if calculate == True:
        print("calculating disbalancing microinteraction...")
        return [1.4, 0.4, 0.5]
    else: return[2.0, 1.96, 1.94]
def above_than(v=float(0),b=float(0)):
    if (v>b):
        return v
    else: return b
#Сила диффузии молекул в воде при н.у.
Niu_Diffusion = (8.9*10e-4)*10**-9
#Блок-Физикохимическое моделирование взаимодействий между ионами
def Niu_coulomb_interaction_Li_F(dist=float):
    power = float(0)
    #!!! - УЧТЕНА ДИЭЛЕКТРИЧЕСКАЯ ПРОНИЦАЕМОСТЬ ДИСТИЛЛИРОВАННОЙ ВОДЫ
    diel=80.2
    #Формула взята не из физ химии Герасимова, но из источника по электростатике
    return (((1.602176634e-19)**2)/dist**2)*(1/(4*math.pi*8.85418782e-12*diel))
def energy_release_from_ions_bonding(mole=float):
    #Термохимически рассчитанная величина энергии спаривания Li+ и F-
    return [-259.714*mole, -584.1*mole,35.9*mole]
def part_energy_release_from_ions_bonding():
    #Термохимически рассчитанная величина энергии спаривания 
    #буквально двух частиц, Li+ и F-
    return [-259.714*(1/Na), -584.1*(1/Na),35.9*(1/Na)]
def modulate_ionic_bonding(volume=float, particles_radius=float):
    #согласно расчетам протяженность вдоль координатных осей около 1,1;
    #принимаю равной 1.
    Li = [randint(0,1) for i in range(3)]
    F = [randint(0,1) for i in range(3)]
    distance = above_than((((((Li[0]-F[0])**2+(Li[1]-F[1]))**2+(Li[2]-F[2])**2)**0.5)*10e-10),2*93.75e-12)
    S = (volume**(1/3))**2
    if Niu_coulomb_interaction_Li_F(distance) > Niu_Diffusion:
        return part_energy_release_from_ions_bonding()
    else: return [0.0,0.0,0.0]
def modulate_dissociaition(calculation=False):
        return [259.714*(1/Na), 584.1*(1/Na),-35.9*(1/Na)]
def modulate_equal_charges_repulsion(volume=float, particles_radius=float):
    #согласно расчетам протяженность вдоль координатных осей около 1,1;
    #принимаю равной 1.
    Li = [randint(0,1) for i in range(3)]
    F = [randint(0,1) for i in range(3)]
    distance = above_than((((((Li[0]-F[0])**2+(Li[1]-F[1]))**2+(Li[2]-F[2])**2)**0.5)*10e-10),2*93.75e-12)
    if Niu_coulomb_interaction_Li_F(distance) > Niu_Diffusion:
        #Логика расчета такая: протяженность фазовой ячейки, вдоль 
        # которой будут растолканы одноименные заряды, есть длина (М), 
        # умножение которой на силу кулоновского взаимодействия (Н)
        # получают Дж
        _dQ=float(distance*Niu_coulomb_interaction_Li_F(distance))
        #При постоянстве давления и объёма макроскопическая 
        # работа по расширению системы равна нулю, след-но 
        # изменение св.вн.эн. равно изменению теплоты
        _dS=_dQ/T
        return [dH(_dQ, P, P, V, V), dG(_dQ, T,T,_dS+S,S,P,P,V,V), _dS]
    else: return [0.0,0.0,0.0]

#Рассчитана подставлением в качестве расстояния в выражение закона Кулона 
# двух радиусов 10е-10
niu_qulon_at_maximum = 1.4383276499573694846e-22
def pl_ions_interaction(calculation=False):
    if (calculation==True):
        return modulate_ionic_bonding(10e-30,93.75e-12)    
    else: 
        _dQ=-float(10e-10*niu_qulon_at_maximum)
        _dS=_dQ/T
        return [dH(_dQ, P, P, V, V), dG(_dQ, T,T,_dS+S,S,P,P,V,V), _dS]
def ds_ions_interaction(calculation=False):
    if (calculation==True):
        return modulate_equal_charges_repulsion(10e-30,93.75e-12)
    else: 
        _dQ=float(10e-10*niu_qulon_at_maximum)
        _dS=_dQ/T
        return [dH(_dQ, P, P, V, V), dG(_dQ, T,T,_dS+S,S,P,P,V,V), _dS]
    
#Блок-термодинамическое моделирование взаимодействия между 
# ионами и молекулами

def Niu_coulomb_interaction_LiF_ion(dist=float):
    power = float(0)
    #!!! - УЧТЕНА ДИЭЛЕКТРИЧЕСКАЯ ПРОНИЦАЕМОСТЬ ДИСТИЛЛИРОВАННОЙ ВОДЫ
    diel=80.2
    LiF_bond_dipol = 2.110593454703e-29 #Кл*М
    #LiF_bond_dipol = 6.3274 #Дебая
    #return ((2*((1.602176634e-19)**2)*LiF_bond_dipol)/dist**3)*(1/(4*math.pi*8.85418782e-12*diel))
    #return (2*LiF_bond_dipol*(1.602176634e-19)*dist)/(dist**4)
    return (LiF_bond_dipol*(1.602176634e-19))/(dist**2)*(1/(4*math.pi*8.85418782e-12*diel))
def modulate_ions_to_molecule_adsorbtion():
    Li_or_F = [randint(0,1) for i in range(3)]
    LiF = [randint(0,1) for i in range(3)]
    #!!! Не учитывается возможная ориентация диполя относительно заряда.
    distance = above_than((((((LiF[0]-Li_or_F[0])**2+(LiF[1]-Li_or_F[1]))**2+(LiF[2]-Li_or_F[2])**2)**0.5)*10e-10)-(93.75e-12),93.75e-12)
    if Niu_coulomb_interaction_LiF_ion(distance)>Niu_Diffusion:
        #Логика расчета такая: протяженность фазовой ячейки, вдоль 
        # которой будет действовать сила электрост. притяжения
        # , есть длина (М), 
        # умножение которой на силу кулоновского взаимодействия (Н)
        # получают Дж
        _dQ=(-1)*float(distance*Niu_coulomb_interaction_Li_F(distance))
        #При постоянстве давления и объёма макроскопическая 
        # работа по расширению системы равна нулю, след-но 
        # изменение св.вн.эн. равно изменению теплоты
        _dS=_dQ/T
        return [dH(_dQ, P, P, V, V), dG(_dQ, T,T,_dS+S,S,P,P,V,V), _dS]
    else: return [0.0,0.0,0.0]




#Рассчитана подставлением в качестве расстояния в выражение закона Кулона 
# двух радиусов 10е-10
niu_qulon_at_maximum = 1.4383276499573694846e-22
def pl_ions_to_molecule_adsorbtion(calculation=False):
    if (calculation==True):
        return modulate_ionic_bonding(10e-30,93.75e-12)    
    else: 
        _dQ=-float(10e-10*niu_qulon_at_maximum)
        _dS=_dQ/T
        return [dH(_dQ, P, P, V, V), dG(_dQ, T,T,_dS+S,S,P,P,V,V), _dS]
def ds_ions_to_molecule_adsorbtion(calculation=False):
    if (calculation==True):
        return modulate_equal_charges_repulsion(10e-30,93.75e-12)
    else: 
        _dQ=float(10e-10*niu_qulon_at_maximum)
        _dS=_dQ/T
        return [dH(_dQ, P, P, V, V), dG(_dQ, T,T,_dS+S,S,P,P,V,V), _dS]



#Большой вопрос - должна ли энергия взаимодействия превысить силу 
#диффузии или же они могут складываться из-за своей направленности?

def test():
 Li = [randint(0,1) for i in range(3)]
 F = [randint(0,1) for i in range(3)]
 distance = above_than((((((Li[0]-F[0])**2+(Li[1]-F[1]))**2+(Li[2]-F[2])**2)**0.5)*10e-10)-2*93.75e-12,2*93.75e-12)
 print(f"Li coordinates {Li} F coordinates {F}")
 print(f"distance equals {distance} метра")
 print(f"сила диффузии {(8.9*10e-4)*10**-9}")
 print(f"сила кулоновского взаимодействия {Niu_coulomb_interaction_Li_F(distance)}")
 print(f"На минимально допустимой дистанции ")
 
def test2():
    Li_or_F = [randint(0,1) for i in range(3)]
    LiF = [randint(0,1) for i in range(3)]
    distance = above_than((((((LiF[0]-Li_or_F[0])**2+(LiF[1]-Li_or_F[1]))**2+(LiF[2]-Li_or_F[2])**2)**0.5)*10e-10),93.75e-12)
    print(f"LiF coordinates {LiF} ion coordinates {Li_or_F}")
    print(f"Расстояние между молекулой и ионом {distance} метра")
    print(f"сила диффузии {Niu_Diffusion}")
    print(f"кулоновская сила притяжения иона к молекуле {Niu_coulomb_interaction_LiF_ion(distance)}")
    if Niu_coulomb_interaction_LiF_ion(distance)>Niu_Diffusion:
        #Логика расчета такая: протяженность фазовой ячейки, вдоль 
        # которой будет действовать сила электрост. притяжения
        # , есть длина (М), 
        # умножение которой на силу кулоновского взаимодействия (Н)
        # получают Дж
        _dQ=(-1)*float(distance*Niu_coulomb_interaction_Li_F(distance))
        #При постоянстве давления и объёма макроскопическая 
        # работа по расширению системы равна нулю, след-но 
        # изменение св.вн.эн. равно изменению теплоты
        _dS=_dQ/T
        #return [dH(_dQ, P, P, V, V), dG(_dQ, T,T,_dS+S,S,P,P,V,V), _dS]
        return [dH(_dQ, P, P, V, V), dG_isobaric_isotermic(_dQ,T,_dS), _dS]
    else: return [0.0,0.0,0.0]

print(f"Параметры изменились следующим образом: {test2()}")





#Блок по химической кинетики для
# дальнйшей сборки массива из взаимодействий

