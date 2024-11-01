import Macroscopic as m

def execute(micro_processes):
    d_result = [0.0,0.0,0.0] # dH, dG, dS
    print(micro_processes)
    cut_=False
    for cut in range(int((len(micro_processes)/2)*0.8)):
     print(f"Cut parameter {cut, int((len(micro_processes)/2)*0.8)}")
     if cut_==True:
        micro_processes = micro_processes[cut:(len(micro_processes)-cut)]
        print(f"Cutted to {micro_processes}")
        d_result = [0.0,0.0,0.0]
     for microfunc in micro_processes:
       try:
            for i in range (0,3):
                d_result[i] += microfunc[i]
       except TypeError:    
            solution = microfunc()
            for i in range (0,3):
                d_result[i] += solution[i]
     print(d_result)
     print(f"dG_macro minus integral dG equals {m.dG_isobaric_isotermic(d_result[0],m.T,d_result[2])-d_result[1]}")
     if (m.dG_isobaric_isotermic(d_result[0],m.T,d_result[2]) < 0 or d_result[1] < 0) and d_result[2] > 0 and d_result[0] < 0:
        print(f"Combination of processes is allowed to occur")
        d_result = [0.0,0.0,0.0]
        for microfunc in micro_processes:
            p = microfunc(True)            
            for i in range (0,3):
                d_result[i] += p[i]
        print(d_result)
        break
     else: 
        print(f"Cutting as returned {d_result[0]}<0 or {d_result[1]}>0 or {d_result[2]}>0")
        cut_=True 
            
l = [m.preset_microprocess_ex for x in range(10)]    
l[::2]=[m.preset_disbalancing_microprocess_ex for x in range(5)]
l+=[m.preset_disbalancing_microprocess_ex_2 for x in range(2)]
execute(l)



#Сборку массива процессов нужно осуществлять на каждый момент времени,
# исходя из кинетических выражений для реакций второго порядка 
# и обладающего наибольшей термодинамической вероятностью равномерного 
# распределения
def combine():