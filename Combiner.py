import Macroscopic as m

def _h(val=float,mid = float):
    if (val<mid):
        return val
    else: return val-mid;
def _h2(val=float,mid=float, low=float,peak=float):
    if (val<mid): 
        return low
    else: return peak
def execute(self,micro_processes):
    d_result = [0.0,0.0,0.0] # dH, dG, dS
    #print(micro_processes)
    #print(micro_processes[2])
    cut_=False
    for cut in range(int((len(micro_processes)/2)*0.8)):
     print(f"Cut parameter {cut, int((len(micro_processes)/2)*0.8)}")
     if cut_==True:
        micro_processes = micro_processes[cut:(len(micro_processes)-cut)]
        d_result = [0.0,0.0,0.0]
        cut_=False
        if (len(micro_processes)==0):
            print("Failed to run sequence")
            break
     for microfunc in micro_processes:
       try:
            for i in range (0,6):
                d_result[_h(i,3)] += microfunc[_h2(i,3,0,1)][_h(i,3)]
       except TypeError:    
            solution = (microfunc[0](),microfunc[1]())
            for i in range (0,6):
                d_result[_h(i,3)] += solution[_h2(i,3,0,1)][_h(i,3)]
     print(d_result)
     print(f"dG_macro minus integral dG equals {self.dG_isobaric_isotermic(d_result[0],self.T,d_result[2])-d_result[1]}")
     if (self.dG_isobaric_isotermic(d_result[0],self.T,d_result[2]) < 0 or d_result[1] < 0) and d_result[2] > 0 and d_result[0] < 0:
        print(f"Combination of processes is allowed to occur")
        d_result = [0.0,0.0,0.0]
        for microfunc in micro_processes:
            p = microfunc(True)            
            for i in range (0,3):
                d_result[i] += p[i]
        print(d_result)
        break
     else: 
        print(f"Cutting as returned {d_result[0]}>0 or {d_result[1]}>0 or {d_result[2]}<0")
        cut_=True 
            
#l = [m.preset_microprocess_ex for x in range(10)]    
#l[::2]=[m.preset_disbalancing_microprocess_ex for x in range(5)]
#l+=[m.preset_disbalancing_microprocess_ex_2 for x in range(2)]
#execute(l)

