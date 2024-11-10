from Electrostatics import Electrostatics
from random import randint
from Combiner import execute
class IonsInteraction:    
    #Рассчитана подставлением в качестве расстояния в выражение закона Кулона 
    # двух радиусов 10е-10
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

class Molecules_Interaction:
    #Рассчитана подставлением в качестве расстояния в выражение закона Кулона 
    # среднего радиуса 10е-10
    top_el_en_molecule_to_ion = 3.7912941330972668921e-20
    def pl_ions_to_molecule_adsorbtion(self,calculation=False)->float:
        if (calculation==True):
            return self.modulate_ionic_bonding(10e-30,93.75e-12)    
        else: 
            _dQ=-float(10e-10*self.top_electrostatic_energy)
            _dS=_dQ/self.T
            return [self.dH(_dQ, self.P, self.P, self.V, self.V), self.dG(_dQ, self.T,self.T,_dS+self.S,self.S,self.P,self.P,self.V,self.V), _dS]
    def ds_ions_to_molecule_adsorbtion(self,calculation=False)->float:
        if (calculation==True):
            return self.modulate_equal_charges_repulsion(10e-30,93.75e-12)
        else: 
            _dQ=float(10e-10*self.top_electrostatic_energy)
            _dS=_dQ/self.T
            return [self.dH(_dQ, self.P, self.P, self.V, self.V), self.dG(_dQ, self.T,self.T,_dS+self.S,self.S,self.P,self.P,self.V,self.V), _dS]



    #Большой вопрос - должна ли энергия взаимодействия превысить силу 
    #диффузии или же они могут складываться из-за своей направленности?

    def test(self):
        print("===========================================")
        Li = [randint(0,1) for i in range(3)]
        F = [randint(0,1) for i in range(3)]
        dis=(((((Li[0]-F[0])**2+(Li[1]-F[1]))**2+(Li[2]-F[2])**2)**0.5)*10e-10)
        distance = self.above_than(dis,(2*93.75e-12))
        print(f"Li coordinates {Li} F coordinates {F}")
        print(f"distance equals {distance} метра")
        print(f"сила диффузии {(8.9*10e-4)*10**-9}")
        print(f"сила кулоновского взаимодействия {self.Niu_coulomb_interaction_Li_F(distance, True)}")
        _dQ=self.Niu_coulomb_interaction_Li_F(distance, False)
        _dS=_dQ/self.T
        return [self.dH(_dQ, self.P, self.P, self.V, self.V), self.dG_isobaric_isotermic(_dQ,self.T,_dS), _dS]
    def test2(self):
        print("===========================================")
        Li_or_F = [randint(0,1) for i in range(3)]
        LiF_center = [randint(0,1) for i in range(3)]
        distance = self.above_than(((((LiF_center[0]-Li_or_F[0])**2+(LiF_center[1]-Li_or_F[1]))**2+(LiF_center[2]-Li_or_F[2])**2)**0.5)*10e-10,93.75e-12)
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

