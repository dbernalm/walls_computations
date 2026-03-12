import numpy as np
from fractions import Fraction as fr
import math
from math import ceil, floor
import time, datetime
from multiprocessing import Pool

def in_ZZ(Fraction):
    return Fraction.denominator==1 

def verify2(tuple, beta, amb_sp):
    r = tuple[0]
    c = tuple[1]
    d = tuple[2]
    e = tuple[3]
    # Expressions
    b = 2*beta*c*(r-1) + (c**2-d)  
    c2 = (beta**3)*r*(r-1)*(r-2) + 3*(beta**2)*c*(r-1)*(r-2) + 3*beta*(r-2)*(c**2 - d) + 2*e + c**3 - 3*c*d
    # Fractions of the expressions above
    b2 = fr(b,2)
    c3 = fr(c2,6)
    if amb_sp=="P3":
        a = (beta**3)*r + 3*(beta**2)*(c + 2*r) + 3*beta*(d + 4*c + fr(11,3)*r) + e + 6*d + 11*c 
    elif amb_sp=="Ab3":
        a = 6*e + 18*beta*d + 18*(beta**2)*c + 6*r*(beta**3)
    elif amb_sp=="Q3":
        a = 2*(beta**3)*r + 6*(beta**2)*(c + fr(3,2)*r) + 6*beta*(d + 3*c + fr(13,6)*r) + 2*e + 9*d + 13*c
    elif amb_sp=="V5":
        a = 5*(beta**3)*r + 15*(beta**2)*(c + r) + 15*beta*(d + 2*c + fr(2,3)*r + fr(2,5)*r) + 5*e + 15*d + 16*c
    elif amb_sp=="V22":
        a = 22*(beta**3)*r + 66*(beta**2)*(c + fr(r,2)) + 66*beta*(d + c + fr(r,6) + fr(r,11)) + 22*e + 33*d + 23*c
    a2 = fr(a, 6)
    return in_ZZ(a2) and in_ZZ(b2) and in_ZZ(c3)

#def verify(tuple,beta,amb_sp):
    if amb_sp=="P3":
        r = tuple[0]
        c = tuple[1]
        d = tuple[2]
        e = tuple[3]
        a = d - fr(c**2, 2) + (beta)*c*(1-r) + (fr(beta**2, 2))*r*(1-r)
        b = e - fr(c, 6) + (beta)*(d - fr(r, 6)) + (fr(beta**2, 2))*c + (fr(beta**3, 6))*r
        c_result = 2*e - c*d + fr(c**3, 6) + (beta)*(d*(2-r) + (c**2)*(3*r-1)) + fr(beta**2, 2)*c*(2 + r*(r-3)) + fr(beta**3, 6)*r*(r-1)*(r-2)
    elif amb_sp=="Ab3":
        r = tuple[0]
        c = tuple[1]
        d = tuple[2]
        e = tuple[3]
        # a = r*(beta**3) + 3*(beta**2)*c + 3*beta*d + e
        b = r*(beta**2) + 2*beta*c + fr(c**2-d,2) 
        c_result = r*beta + fr(c**3 + 2*e - 3*c*d, 6)
    return in_ZZ(b) and in_ZZ(c_result) 

def grandverify(list,beta, amb_sp):
    a = len(list)
    granddefinitive = []
    #walls = []
    for i in range(a):
        r1 = math.ceil(list[i][0])
        r2 = math.floor(list[i][1])
        definitive = ([], list[i][2],list[i][3],list[i][4])
        wl = fr(list[i][4], list[i][2]) 
        for j in range(r1, r2+1):
            c=(j,list[i][2],list[i][3],list[i][4])
            if verify2(c, beta, amb_sp):
                definitive[0].append(j)
        if len(definitive[0])>0:
            #walls.append(wl)
            granddefinitive.append(definitive)
    return granddefinitive

class Sheaf:
    def __init__(self, R, D, k, amb_sp):
        self.R = -R
        self.D = D
        self.chern = (R, 0, D, 0)
        self.k = k
        self.amb_sp = amb_sp
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
        for d in range(1, ((self.k)**2)*self.D+1):
            trued = fr(d, ((self.k)**2))
            limit1 = trued**2
            limit2 = (self.D - trued)**2
            min_limit = math.floor( ((self.k)**4)*min(limit1,limit2) )
            for c in range(1, min_limit + 1):
                for e in range(1, min_limit + 1):
                    if 0 < c * e <= min_limit:
                        truec = fr(c, self.k)
                        truee = fr(e, (self.k)**3)
                        # Ranks
                        rank1 = fr(-truec*(self.D - trued), truee) - self.R
                        rank2 = fr(truec*trued,truee)
                        solutions.append((rank1,rank2,truec,trued,truee))
        return solutions

    def num_dest(self, number_division):
        list = self.possible_c6e()
        a = len(list)
        lenlist = math.floor( fr(a, number_division) )
        listblocks = [list[x:x+lenlist] for x in range(0, a, lenlist)]
        b = len(listblocks)

        if __name__ == "__main__":
            with Pool(processes=8) as pool:
                results = [pool.apply_async(grandverify, (listblocks[i], self.beta, self.amb_sp)) for i in range(b)]
                output = [res.get() for res in results]
                # Non-optimized but usually don't take a lot of time
                output2 = [x for x in output if x!=[]]
                self.dest = [output2[i][j] for i in range(len(output2)) for j in range(len(output2[i]))]
        return self.dest

t1 = time.time()
# Change the ones that have beta non zero if beta=0 works
a = Sheaf(0,2,-2, "Q3")
print(a.num_dest(1))
t2 = time.time()
total = t2-t1
print(str(datetime.timedelta(seconds=total)))