from Kinetics import Kinetics
class CalculatedInteraction:
    # Предполагается, что для микропроцесса просчитаны 
    # все возможные исходы
    # impacts = {'описание результата взаимодействия', (dH_dG_dS_electrost, quantum_impact[])}
    def __init__(self, **impacts):
        self.impacts = impacts
        self.pull_size = len(impacts.items())
    # Возвращает массив вида dH,dG,dS
    def summorize(self,variant, quantum_variant):
        try:
            if (isinstance(variant,int)==True):
                for macro,quantum in self.impacts[variant]:
                    return [macro_impact[0]+quantum_impact[quantum_variant], macro_impact[1], macro_impact[2]]
            if (isinstance(variant,str)==True):
                for name, effects in self.impacts.items():
                    if variant==name:
                        macro_impact, quantum_impact = effects
                        return [macro_impact[0]+quantum_impact[quantum_variant], macro_impact[1], macro_impact[2]]                
        except: 
            print(f"Exception occured while summorizing effect of an interaction {self.impacts[variant]}")
            return 0
class KineticsInterpreter:
    # interactions = {'описание взаимодействия', CalculatedInteraction}
    def __init__(self, **interactions):
        self.interactions = interactions
    
    #Самоконструктор, не требующий ввода стандартных значений
    #предполагается, что он самостоятельно будет собирать 
    #матрицы возможных взаимодействий
    def __init__(self):
        self.interactions = {}
    
    def interpret_reaction_speed(self, t1,t2 , sys=Kinetics):
        reaction_speed_gradient = [sys.reaction_speed(sys, t) for t in range(t1,t2)]
        return_values = [0,0,0]
        for sec in range(0,len(reaction_speed_gradient)):
            interacting_particles = sys.Na_*reaction_speed_gradient[sec] # частиц/моль*моль/сек=частиц/сек;
            # убыль концентрации в молях умножаю на постоянную Авогадро и 
            # получаю убыль частиц каждого из реагентов
            # так как каждое взаимодействие включает две частицы реагентов, 
            # число взаимодействий принимаю за половину взаимодействующих частиц
            interactions_amount=0.5*interacting_particles
            
            #Буду исходить из предположения, что все взаимодействия равновероятны
            each_interaction_repetitions = interactions_amount/len(self.interactions.items())
            
            for name,calcul_inter in self.interactions.items():
                #Пусть все исходы взаимодействий также равновероятны                
                each_hypothetical_interaction_output_repetition = each_interaction_repetitions/len(calcul_inter.impacts.items())
                for nname, output in calcul_inter.impacts.items():
                    #Пусть все варианты взаимодействий на квантовом уровне у каждого исхода взаимодействий равновероятны
                    macro,quantum = output
                    each_quantum_interaction_repetition = each_hypothetical_interaction_output_repetition/len(quantum)
                    for quantum_variant in quantum:
                        return_values+=[macro[0]+quantum_variant,macro[1],macro[2]]*each_quantum_interaction_repetition
        return return_values
                    
                
            