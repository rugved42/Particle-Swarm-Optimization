#!/usr/bin/env python
import math
import numpy as np
from numpy import inf
import random
import matplotlib.pyplot as plt

class pso(): 
    def __init__(self,N, n, domain): 
        self.N = N  #Total number of particles
        self.n = n # Total design variable
        self.alpha = 0.5 #Inertia Component [0,1]
        self.beta1 = 1.49 
        self.beta2 = 1.49
        self.domain = domain
        self.pbest = [inf]*self.N
        self.gbest = inf
        self.gbest_iter = []
        self.gbest_iter_pos = []
        self.mean_obj = []
        self.vel = [[0]*self.n for i in range(self.N)]
        self.pos = self.particles()
        self.pbest_pos = [[0]*self.n for i in range(self.N)]
        self.gbest_pos = []
        self.f1 = [0]*self.N
        
    def particles(self): 
        pop = [list(random.uniform(self.domain[0], self.domain[1]) for i in range(self.n)) for t in range(self.N)]
        return pop

    def obj_fun(self): 
        tmp = self.pos    
#         self.f1 = [sum(np.array(i)**2) for i in tmp] # "De-jong's Function"
        self.f1 = [(sum(100*(i[x+1] - i[x]**2)**2 + (1-i[x])**2 for x in range(len(i) - 1))) for i in tmp] #Rosenbrock Function
#         self.f1 = [(10*self.n + sum((i[x]**2 - 10*np.cos(2 * np.pi * i[x])) for x in range(len(i)))) for i in tmp] #Rastrigin Function
#         self.f1 = [(sum(-1*i[x] * np.sin(np.sqrt(abs(i[x]))) for x in range(len(i)))) for i in tmp] #Schwefel Function
#         self.f1 = [(-1 *20 * np.exp( -0.2 * np.sqrt(np.abs(sum(i[x]**2 for x in range(len(i))) * 1/self.n))) -1 * np.exp(sum(np.cos(np.pi*2*i[x]) for x in range(len(i))) * 1/self.n) + 20 + np.exp(1))for i in tmp] #Ackley Function
#         self.f1 = [(-1 *sum(np.sin(i[x]) * (np.sin((x * (i[x])**2) * 1/np.pi))**20 for x in range(len(i)))) for i in tmp]
#         self.f1 = self.griewangk_function()
    
    def griewangk_function(self):
        tmp = self.pos
        f = []
        for x in tmp:
            m,s = 0,0
            for i,j in enumerate(x):
                m *= np.cos(j/np.sqrt(i+1))
                s += (j**2)/4000
            total = s - m 
            f.append(total)
        return f
            
                
    
        
    def fitness_ros(self):
        tmp = self.f1
        tmp_1 = self.pos
        for i in range(len(tmp)):
            if tmp[i] < self.pbest[i]:
                self.pbest[i] = tmp[i]
                for j in range(self.n):
                    self.pbest_pos[i][j] = tmp_1[i][j]
            else:
                continue
#         print("P_best array {}".format(self.pbest))
        self.gbest = min(self.pbest)
        t = np.where(np.array(self.pbest) == self.gbest)[0][0]
        self.gbest_pos = self.pbest_pos[t]
        self.gbest_iter.append(self.gbest)
        self.gbest_iter_pos.append(tuple(self.gbest_pos))
        print("position of g_best_ {}" .format(self.gbest_pos))
        
    def calc_ros(self): 
        tmp1 = self.pos
        pos_new = self.pos
        for i in range(len(tmp1)):
            # print(i)
            v_new = (self.alpha*np.array(self.vel[i])) + ((self.beta1*random.random()) * (np.array(self.pbest_pos[i]) - np.array(tmp1[i]))) + ((random.random()*self.beta2) * (np.array(self.gbest_pos) - np.array(tmp1[i])))
            # print("vnew",v_new)
            self.vel[i] = v_new
        # print("Velocity", vel)
        for i in range(len(tmp1)):
            pos_new[i] = tmp1[i] + self.vel[i]
            for c in range(len(pos_new[i])):
                if pos_new[i][c] > self.domain[1]:
                    pos_new[i][c] = self.domain[1]
                if pos_new[i][c] < self.domain[0]:
                    pos_new[i][c] = self.domain[0]
        self.pos = pos_new
        # print(pos_new)
    
    def run_pso(self):
        i = 0
        while i < 25:
#             print("Particles..:", self.pos)
            self.obj_fun()
            self.fitness_ros()
            self.calc_ros()
            i += 1
        plt.plot(self.gbest_iter)
        plt.show()
if __name__ == "__main__":
    pso_run = pso(100, 3, (-2.048,2.04))
#     f = open("pso_ackley.txt", "w")
#     f.write(" Function: {}\n Number of Particles: {}\n Number of Design Variables: {}\n Domain:{}".format("Ackley_Function",pso_run.N,pso_run.n,pso_run.domain))
#     f.write("\n")
    run = 0
    while run <= 10:
        pso_run.run_pso()
#         f = open("pso_ackley.txt", "a")
#         f.write("Run:{}/{} Global_Best:{} Global_Best_Position:{}".format(run,10,pso_run.gbest,pso_run.gbest_pos))
#         f.write("\n")
#         plt.boxplot(pso_run.gbest_iter,showfliers=False)
        run += 1