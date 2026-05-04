import numpy as np
from fractions import Fraction as fr
import math
from math import ceil, floor
import time, datetime
from multiprocessing import Pool

def in_ZZ(Fraction):
    return Fraction.denominator==1

def untwist(r,c,d,e,beta):
  R = r
  C = c + (beta)*R
  D = - fr(beta**2,2)*R + beta*C + d
  E = fr(beta**3,6)*R - fr(beta**2,2)*C + beta*D + e
  return [R,C,D,E]

def twist(R,C,D,E,beta):
  r = R
  c = C - (beta)*R
  d = fr(beta**2,2)*R - beta*C + D
  e = -fr(beta**3,6)*R + fr(beta**2,2)*C - beta*D + E
  return [r,c,d,e]

def verify(tuple,beta):
    r = tuple[0]
    c = tuple[1]
    d = tuple[2]
    e = tuple[3]
    a = d - fr(c**2, 2) + (beta)*(1-r)*(fr(beta, 2)*r + c)
    b = e - fr(c, 6) + (beta)*(d - fr(r, 6)) + (fr(beta**2, 2))*c + (fr(beta**3, 6))*r
    #c = 2*e - c*d + fr(c**3, 6) + (beta)*(d*(2-r) + (c**2)*(3*r-1)) + fr(beta**2, 2)*c*(2 + r*(r-3)) + fr(beta**3, 6)*r*(r-1)*(r-2)
    c1 = 2*e - c*d + fr(c**3, 6) + (beta)*(d*(2-r) + (c**2)*(fr(r-2,2)) ) + fr(beta**2, 2)*c*(2 + r*(r-3)) + fr(beta**3, 6)*r*(r-1)*(r-2)
    return in_ZZ(a) and in_ZZ(b) and in_ZZ(c1)

def grandverify(list,beta):
    a = len(list)
    granddefinitive = []
    #walls = []
    for i in range(a):
        r1 = math.ceil(list[i][0])
        r2 = math.floor(list[i][1])
        definitive = ([], list[i][2],list[i][3],list[i][4])
        wl =  math.sqrt(fr(6*list[i][4], list[i][2])) 
        # This if comes from a conjecture. Change to 1 or 0 if needed. 
        # More than alpha_0^2 plus epsilon
        if wl > (fr(1,2)):
            for j in range(r1,r2+1):
                c=(j,list[i][2],list[i][3],list[i][4])
                if verify(c,beta):
                    definitive[0].append(j)
        if len(definitive[0])>0:
            #walls.append(wl)
            # untwist stop: comment next line if this does now work
            #for j in range(len(definitive[0])):
                #unt = untwist(definitive[0][j], list[i][2], list[i][3], list[i][4], beta)
            # untwisted line
            granddefinitive.append(definitive)
            # granddefinitive.append(unt)
    return granddefinitive

class Sheaf:
    def __init__(self, R, D, k):
        self.R = -R
        self.D = D
        self.chern = (R, 0, D, 0)
        self.k = k
        if k==1:
            self.beta = 0
        else:
            self.beta = fr(1, k)
        self.walls = set()
        self.dest = []
        self.maxwall = 0
        self.maxobject = []
    
    def possible_c6e(self):
        solutions = []
        for d in range(1, 2*((self.k)**2)*self.D):
            trued = fr(d, 2*((self.k)**2))
            limit1 = 4*(trued**2)
            limit2 = 4*((self.D - trued)**2)
            min_limit = math.floor( ((self.k)**4)*min(limit1,limit2) )

            for c in range(1, min_limit + 1):
                for e in range(1, min_limit + 1):
                    if 0 < c * e <= min_limit:
                        truec = fr(c, self.k)
                        truee = fr(e, 6*((self.k)**3))
                        # Ranks
                        rank1 = fr(-truec*(2*self.D - 2*trued), 6*truee) - self.R
                        rank2 = fr(truec*2*trued, 6*truee)
                        solutions.append((rank1, rank2, truec, trued, truee))
        return solutions
    
    def num_dest(self, number_division):
        list = self.possible_c6e()
        a = len(list)
        lenlist = math.floor( fr(a, number_division) )
        listblocks = [list[x:x+lenlist] for x in range(0, a, lenlist)]
        b = len(listblocks)

        if __name__ == "__main__":
            with Pool(processes=8) as pool:
                results = [pool.apply_async(grandverify, (listblocks[i], self.beta)) for i in range(b)]
                output = [res.get() for res in results]
                # Non-optimized but usually don't take a lot of time
                output2 = [x for x in output if x!=[]]
                self.dest = [output2[i][j] for i in range(len(output2)) for j in range(len(output2[i]))]
        return self.dest
    
    def walls_dest(self):
        for i in range(len(self.dest)):
            wl = math.sqrt(fr(6*self.dest[i][3], self.dest[i][1]))
            self.walls.add(wl)
            if wl > self.maxwall:
                self.maxwall = wl
                self.maxobject = self.dest[i]
        return sorted(self.walls, reverse=True)

t1 = time.time()
a = Sheaf(0,4,2)
print(a.num_dest(100))
print(a.walls_dest())
t2 = time.time()
total = t2-t1
print(str(datetime.timedelta(seconds=total)))