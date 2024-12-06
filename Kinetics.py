import math
import matplotlib.pyplot as plt
from Electrostatics import Electrostatics
from Sampler import KineticsInterpreter 
from Sampler import CalculatedInteraction
from Combiner import execute
from Processes import IonsInteraction
from Processes import Molecules_Interaction
class Kinetics(Electrostatics):
            
    top_electrostatic_energy = 2.8780165475731916883e-20
    def pl_ions_interaction(system,calculation=False, args=(0.0,0.0))->float:
        Li,F = args
        if (calculation==True):
            return system.modulate_ionic_bonding_(args)    
        else: 
            _dQ=float(-1*system.top_electrostatic_energy)
            _dS=_dQ/system.T
            return [system.dH(_dQ, system.P, system.P, system.V, system.V), system.dG(_dQ, system.T,system.T,_dS+system.S,system.S,system.P,system.P,system.V,system.V), _dS]
    def ds_ions_interaction(system,calculation=False, args=(0.0,0.0))->float:
        if (calculation==True):
            return system.modulate_equal_charges_repulsion_(args)
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
    D_0=1.33*10**-9 #коэффициент диффузии раствора ионов, до образования коллоидных частиц
    
    #(10e10*10e-14)/0.0616#4*math.pi*D_0*9.375e-11 #данные эксперимента, константа скорости
    r_LiF=(145+50)**-12
    r_LiF_Col0=9.641e-8
    r_LiF_Col=r_LiF_Col0
    D=D_0
    av_speed = 0.01 #средняя скорость до начала образования осадк
    
    K = 0.6667 #данные эксперимента, константа скорости
    def reaction_speed(self, k,C_A, C_B, n_A, n_B):
        #n_A частный порядки реакции по фтору
        #n_B частный порядок реакции по литию
        return k*(C_A**n_A)*(C_B**n_B)
    
    def interpret_reaction_speed(self, interpreter, speed):
        return_values = [0.0,0.0,0.0]
        particles_created = self.Na_*(speed*self.V) # частиц/моль*моль/сек=частиц(продукта)/сек;
            # убыль концентрации в молях умножаю на постоянную Авогадро и 
            # получаю убыль частиц каждого из реагентов
            # так как каждое взаимодействие включает две частицы реагентов, 
            # число взаимодействий принимаю за половину взаимодействующих частиц
        interactions_amount=3*particles_created #пусть на каждое образование 
        #пар столкновения приходится отталкивание одноименных зарядов
            
            #Буду исходить из предположения, что все взаимодействия равновероятны
        overall_amount=0
        for calcul_inter in interpreter.interactions:
            overall_amount+=len(calcul_inter.impacts)
        each_interaction_repetitions = interactions_amount/overall_amount
        for calcul_inter in interpreter.interactions:
                #Пусть все исходы взаимодействий также равновероятны    
            each_hypothetical_interaction_output_repetition = each_interaction_repetitions
            print(f"REPET {each_hypothetical_interaction_output_repetition}")
            for output in calcul_inter.impacts:
                    #Пусть все варианты взаимодействий на квантовом уровне у каждого исхода взаимодействий равновероятны
                macro = output[0]
                quantum = output[1]
                each_quantum_interaction_repetition = each_hypothetical_interaction_output_repetition/len(quantum)
                return_values[0]+=macro[0]*each_hypothetical_interaction_output_repetition
                return_values[1]+=macro[1]*each_hypothetical_interaction_output_repetition
                return_values[2]+=macro[2]*each_hypothetical_interaction_output_repetition
                for quantum_variant in quantum:
                    return_values[0]+=0
        return return_values
        

    def run(self,t):  
        print(f"Просчитываются элементарные взаимодействия...")
        inter_Li_F = KineticsInterpreter()
        inter_LiF = KineticsInterpreter()      
        quantum_int=[0.0]
        macroscopic_effects_Li_F_positive=[(self.part_energy_release_from_ions_bonding(),[0.0])]
        macroscopic_effects_Li_F_negative=[([0.0,0.0,0.0],[0.0])]
        macroscopic_effects_LiF=[([0.0,0.0,0.0],[0.0])]
        CnREG=[0.5]
        CnPR=[0.0]
        dC_col=[0.0]
        C_col = [0.0]
        C_col_ = 0.0
        r_col=[0.0]
        d_H=[0.0]
        d_G=[0.0]
        d_S=[0.0]
        for p_Li_X in range(2):
            for p_Li_Y in range(2):
                for p_Li_Z in range(2):
                    for p_F_X in range(2):
                        for p_F_Y in range(2):
                            for p_F_Z in range(2):
                                eff=self.pl_ions_interaction(True, ([p_Li_X,p_Li_Y,p_Li_Z],[p_F_X,p_F_Y,p_F_Z]))
                               # if eff[0]<0:
                                macroscopic_effects_Li_F_positive.append((eff, quantum_int))
                                #else:                                    
                                    #macroscopic_effects_Li_F_negative.append((eff, quantum_int))
                                macroscopic_effects_Li_F_negative.append((self.ds_ions_interaction(True, ([p_Li_X,p_Li_Y,p_Li_Z],[p_F_X,p_F_Y,p_F_Z])), quantum_int))
                                macroscopic_effects_Li_F_negative.append((self.ds_ions_interaction(True, ([p_Li_X,p_Li_Y,p_Li_Z],[p_F_X,p_F_Y,p_F_Z])), quantum_int))
                                macroscopic_effects_LiF.append((self.modulate_ions_to_molecule_adsorbtion_(([p_Li_X,p_Li_Y,p_Li_Z],[p_F_X,p_F_Y,p_F_Z])), quantum_int))
        inter_Li_F.interactions.append(CalculatedInteraction(macroscopic_effects_Li_F_positive))
        inter_Li_F.interactions.append(CalculatedInteraction(macroscopic_effects_Li_F_negative))
        #inter_Li_F.interactions.append(CalculatedInteraction(macroscopic_effects_LiF))
        # V=L*C(Li)*C(F-)
        C_li0=0.5
        C_F0=0.5 #моль/л противоионов на момент начала реакции t=0, численно равна концентрации фторида лития, которая установилась бы в растворе, если он был бы неограниченно растворим
        C_li=C_li0
        C_saved = 0.0
        C_f=C_F0
        C_LiF=0
        d_parmateres=[0.0,0.0,0.0]
        print("==================================================")
        print(len(inter_Li_F.interactions[0].impacts))
        print(inter_Li_F.interactions[0].impacts)
        print("==================================================")
        print(len(inter_Li_F.interactions[1].impacts))
        print(inter_Li_F.interactions[1].impacts)
        print("==================================================")
        #print(len(inter_Li_F.interactions[0].impacts))
        #print(inter_Li_F.interactions[0].impacts)
        check=input()
        dc=self.reaction_speed(self.K, C_li,C_f,1,1)
        dCn = [dc]
        Dn = [self.D]
        for tt in range(0,t):
            print("================================================================") 
            print(f"T={tt} сек======================================================")
            #Расчет скорости при с учетом радиусов ионов и коэффициента
            # диффузии ионов в идеальном растворе
            if (C_li-dc>0.043 and C_f-dc>0.043):
                dc=self.reaction_speed(self.K, C_li,C_f,1,1)
            else:
                dc=0
            #if (C_LiF>0.043):
                #dc=self.reaction_speed(self.K, C_li,C_f,1,1)-(1/self.K)*C_LiF**2
            print(f"Скорость образования пар столкновения ионов = {dc}")
            d_par=self.interpret_reaction_speed(inter_Li_F,dc)
            d_parmateres[0]+=d_par[0]
            d_parmateres[1]+=d_par[1]
            d_parmateres[2]+=d_par[2]
            C_LiF+=dc
            C_li-=dc
            C_f-=dc
            nu_0=8.9*10**-4
            nu_=nu_0
            CnREG.append(C_li)
            print(f"Концентрации: {C_LiF}(+{dc})")
            
            #Прирост среднего радиуса коллоидных частиц в секунду
            dr_dt=2.4619e-8*(C_LiF-0.052)
            self.r_LiF_Col+=dr_dt
            
            #dc_coagulation=2.1587*C_LiF**2
            dc_coagulation=0
            
            
            if (C_LiF>0.052):
                #Преодоление концентрации насыщенного раствора и начало образования 
                #коллоидных частиц            
                
                
                
                
                fi=((4/3)*math.pi*(self.r_LiF_Col)**3)*C_col_*0.1*self.Na_/0.001
                #self.K=4*math.pi*((2*8.314*self.T)/(3*math.pi*nu_*9.375e-11))*9.375e-11
                print(f"Коэффициент диффузии = {self.D}")
                nu_=nu_0*(1+2.5*fi+14.1*(fi**2))
                self.D=(2*8.314*self.T)/(3*math.pi*nu_*self.r_LiF_Col)
                #dc_coagulation=(5.16*(0.1*C_LiF)*2.64*25.94)/(72.8*(10**-3)*3*self.r_LiF_Col*0.1)
                
                #Убыль концентрации фторида лития в связи с пересыщением раствора
                dc_coagulation=(23.603*(0.1*C_LiF)*2.64*self.r_LiF_Col)/((3*72.86*10e-3*25.94)*0.1*10e-6)
                C_LiF-=dc_coagulation
                C_col_+=dc_coagulation
                
                if (dc>0 and C_LiF>C_saved):
                    d_parmateres[1]+=21.603*(dc)*0.1
                    #d_parmateres[2]+=25.95816*(dc)*0.1/self.T
                    C_saved=C_LiF
                d_parmateres[2]-=25.95816*(dc_coagulation)*0.1/self.T
                #self.r_LiF_Col+=self.dr_dc*(dc)
                #C_col_+=dc_coagulation
                C_col.append(C_col_)                                
                dC_col.append(dc_coagulation)
                Dn.append(self.D)
                print(f"Скорость образования коллоидных частиц = {dc_coagulation}")
                #self.r_LiF+=self.dr_dc*dc_coagulation
                
                #Скорость образования и роста центров кристаллизации
            r_col.append(self.r_LiF_Col)
            dCn.append(dc)
            CnPR.append(C_LiF)
            d_H.append(d_parmateres[0])
            d_G.append(d_parmateres[1])
            d_S.append(d_parmateres[2])
        x = [i for i in range(0,len(dCn))]
        fig, axs = plt.subplots(2, 5)
        axs[0, 0].plot(x, dCn)
        axs[0, 0].set_title('Скорость образования пар столкновения LiF')
        axs[0, 1].plot(x, CnREG, 'tab:orange')
        axs[0, 1].set_title('Концентрация Li, F')
        axs[0, 2].plot(x, CnPR, 'tab:green')
        axs[0, 2].set_title('Концентрация LiF')
        x = [i for i in range(0,len(dC_col))]
        axs[0, 3].plot(x, dC_col, 'tab:blue')
        axs[0, 3].set_title('Скорость образования колл. частиц')
        x = [i for i in range(0,len(d_H))]
        axs[0, 4].plot(x, d_H, 'tab:blue')
        axs[0, 4].set_title('Энтальпия')
        axs[1, 3].plot(x, d_G, 'tab:blue')
        axs[1, 3].set_title('Энергия Гиббса')
        axs[1, 4].plot(x, d_S, 'tab:blue')
        axs[1, 4].set_title('Энтропия')
        #x = [i for i in range(0,len(Dn))]
        #axs[1, 0].plot(x, Dn, 'tab:red')f
        #axs[1, 0].set_title('Коэффициент диффузии')
        x = [i for i in range(0,len(C_col))]
        axs[1, 1].plot(x, C_col, 'tab:red')
        axs[1, 1].set_title('Концентрация  колл. частиц')
        x = [i for i in range(0,len(r_col))]
        axs[1, 2].plot(x, r_col, 'tab:red')
        axs[1, 2].set_title('Средний радиус  колл. частиц')

        for ax in axs.flat:
            ax.set(xlabel='время', ylabel='величина')
        print (f"Результат {d_parmateres}")
        plt.show()
        #plt.ylabel("Скорость, моль/л*с")
        #plt.title("Скорость образования пар столкновения LiF")
        #plt.ylabel("Концентрация, моль/л")
        #plt.title("Концентрация Li, F")
        #plt.plot(x,CnREG)
        #plt.ylabel("Концентрация, моль/л")
       # plt.title("Концентрация LiF")
        #plt.plot(x,CnPR)
       # plt.ylabel("D")
        #plt.title("Коэффициент диффузии")
        #plt.plot(x,Dn)
        #plt.ylabel("Скорость, моль/л*с")
        #plt.title("Скорость образования коллоидных частиц")
        #plt.plot(x,dC_col)
        #plt.ylabel("Концентрация, моль/л")
        #plt.title("Концентрация  коллоидных частиц")
        #plt.plot(x,C_col)
        #plt.ylabel("Радиус, м")
        #plt.title("Средний радиус  коллоидных частиц")
        #plt.plot(x,r_col)
           
        
        return d_parmateres
    
test_system = Kinetics(101325, 298.15, 1,10e-6,0.1,100)
print (f"Результат {test_system.run(90)}")