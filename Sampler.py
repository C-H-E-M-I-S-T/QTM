
class CalculatedInteraction:
    # Предполагается, что для микропроцесса просчитаны 
    # все возможные исходы
    # impacts = [(dH_dG_dS_electrost, quantum_impact[])]
    def __init__(self, impacts=[([0.0,0.0,0.0],[0.0])]):
        self.impacts = impacts
    # Возвращает массив вида dH,dG,dS
    def summorize(self,variant, quantum_variant):
        try:
            if (isinstance(variant,int)==True):
                for macro,quantum in self.impacts[variant]:
                    return [macro[0]+quantum[quantum_variant], macro[1], macro[2]]
        except: 
            print(f"Exception occured while summorizing effect of an interaction {self.impacts[variant]}")
            return 0
class KineticsInterpreter:
    # interactions = [CalculatedInteraction]
    def __init__(self, *interactions):
        self.interactions = interactions
    
    #Самоконструктор, не требующий ввода стандартных значений
    #предполагается, что он самостоятельно будет собирать 
    #матрицы возможных взаимодействий
    def __init__(self):
        self.interactions = []
    
