from itertools import chain, combinations

with open("raysS5.dat" ,'r') as rF:
    rFLines = rF.readlines()
    raysS5 = list(map(lambda x: vector(map(lambda y: int(y), x.split(" "))), rFLines))
# rays of TS5 from p. 21 of the paper.

with open("linealityS5.dat" ,'r') as lF:
    lFLines = lF.readlines()
    linealityS5 = list(map(lambda x: vector(map(lambda y: int(y), x.split(" "))), lFLines))
# lineality of TS5 from p. 21 of the paper.

def coneMatrix(rays, lineality, cone):
    return matrix(QQ, [rays[i] for i in cone]+lineality)
# matrix whose rows are the rays and lineality of the cone, where the cone is given by a list of ray indices.

E5 = list(map(lambda x: tuple(x), [[0], [1], [2], [3], [4], 
[0,1,2], [0,1,3], [0,2,3], [1,2,3], [0,1,4], [0,2,4], [1,2,4],
[0,3,4], [1,3,4], [2,3,4], [0,1,2,3,4] ]))
# odd subsets of [5].

def S5ToSE5(sigma):
    sigmaE5 = list(map(lambda x: tuple(sorted(list(map(lambda y: sigma[y] , x )))), E5   ))
    return tuple(map(lambda x: E5.index(x), sigmaE5))
# given a permutation of S5, returns a permutation of SE5 = S16.  

def sigmaOnVector(sigma, v):
    sigmaInv = [sigma.index(i) for i in range(len(sigma))]
    sigmaE5 = S5ToSE5(sigmaInv)
    return vector(map(lambda i: v[sigmaE5[i]], range(len(v)) ))
# given a permutation of S5 and vector in R^ES5 = R^16, returns the vector sigma.v by permuting the coordinates.


def twistBymu(mu):
    twistmu = list(map(lambda x: tuple(sorted(list(set(x).symmetric_difference(set(mu))))), E5))
    return tuple(map(lambda x: E5.index(x), twistmu))
# given an even subset mu of [5] of S5, returns a permutation of SE5 = S16 induced by twisting the odd subsets of [5].  

def twistOnVector(mu, v):
    twistmu = twistBymu(mu)
    return vector(map(lambda i: v[twistmu[i]], range(len(v)) ))
# given an even subset mu [5] and vector v in R^ES5 = R^16, returns the vector obtained by permuting the coordinates induced by the twist of mu on the odd subsets of [5].

def sigmaOnCone(tau, sigma):
    rs =  list(map(lambda x: vector(x), tau.rays()))
    ls =  list(map(lambda x: vector(x), tau.lines()))
    sigmars = list(map(lambda v: sigmaOnVector(sigma, v), rs))   
    return Polyhedron(rays = sigmars, lines = ls, base_ring=QQ)
# given cone tau and permutation sigma of [5] returns new cone induced by the permutation of the coordinates by sigma.

def twistmuCone(tau, mu):
    rs =  list(map(lambda x: vector(x), tau.rays()))
    ls =  list(map(lambda x: vector(x), tau.lines()))
    twistrs = list(map(lambda v: twistOnVector(mu, v), rs))   
    return Polyhedron(rays = twistrs, lines = ls, base_ring=QQ)
# given cone tau and even subset mu of [5] returns new cone induced by the permutation of the coordinates by twisting by mu.

    
def sigmaTwistCone(taui, sigma, mu):
    return sigmaOnCone(twistmuCone(indexToCone(taui), mu), sigma)
# given cone as a tuple of ray indices taui, permutation sigma of [5], and even subset mu, returns the cone obtained by twisting by mu, then permuting by sigma. 

def indexToCone(taui):
    return Polyhedron(rays = [raysS5[i] for i in taui], lines = linealityS5, base_ring=QQ)
# given a tuple of indices of rays, returns the cone with these rays.
 
 

# Atau1 is a dictionary of the data from Table 7.1 in the paper. 
Atau1 = {}
Atau1[tuple([])] = [[(0,1,2,3,4), set([2,3]), (4,7,8,9,26)], [(0,3,2,1,4), set([0,3]), (4,7,8,9,26)]]
Atau1[tuple([26])] = [[(0,2,1,3,4), set([]), (4,7,8,9,26)], [(0,1,2,4,3), set([]), (4,7,8,9,26)], [(0,2,1,4,3), set([]), (4,7,8,9,26)]]
Atau1[tuple([4])] = [[(0,4,1,3,2), set([1,4]), (4,7,8,9,26)], [(2,1,4,3,0), set([0,3]), (4,7,8,9,26)]]
Atau1[(3,4)] = [[(2,0,1,3,4), set([3,4]), (3,4,5,7,9)], [(2,1,0,4,3), set([]), (3,4,5,7,9)]]
Atau1[(4,5)]  = [[(3,4,0,2,1), set([2,3]), (4,7,8,9,26)], [(2,3,0,1,4), set([]), (4,7,8,9,26)]]
Atau1[(4,26)] = [[(0,1,2,3,4), set([]), (4,7,8,9,26)], [(1,2,0,3,4), set([]), (4,7,8,9,26)], [(0,2,1,3,4), set([]), (4,7,8,9,26)]]
Atau1[(4,5,9)] = [[(0,1,2,3,4), set([]), (3,4,5,7,9)], [(0,1,4,2,3), set([2,3]), (3,4,5,6,25)]]
Atau1[(4,9,26)] = [[(0,1,2,3,4), set([]), (4,7,8,9,26)], [(3,4,1,0,2), set([0,2]), (3,4,5,6,25)], [(3,4,1,0,2), set([0,3]), (3,4,5,6,25)]]
Atau1[(4,6,26)] = [[(1,2,0,3,4), set([]), (4,7,8,9,26)], [(0,2,1,3,4), set([]), (4,7,8,9,26)]]
Atau1[(3,4,5)] = [[(0,1,2,3,4), set([]), (3,4,5,7,9)], [(0,1,2,3,4), set([]), (3,4,5,6,25)]]
Atau1[(3,4,26)] = [[(4,3,0,2,1), set([0,3]), (3,4,5,6,25)], [(3,4,0,2,1), set([0,3]), (3,4,5,6,25)]]
Atau1[(4,7,8,9)] = [[(0,1,2,3,4), set([]), (4,7,8,9,26)], [(0,1,3,2,4), set([]), (4,7,8,9,26)]]
Atau1[(3,4,5,25)] = [[(0,1,2,3,4), set([]), (3,4,5,6,25)], [(0,1,2,4,3), set([]), (3,4,5,6,25)]]
Atau1[(3,4,5,6)] = [[(0,2,1,3,4), set([]), (3,4,5,7,9)], [(1,0,2,4,3), set([]), (3,4,5,7,9)]]
Atau1[(3,4,6,26)] = [[(3,4,2,0,1), set([0,3]), (3,4,5,6,25)], [(3,4,2,1,0), set([0,3]), (3,4,5,6,25)]]

# Atau2 is a similar dictionary, just that the permutation and twist are carried out by hand. 
Atau2 = {}
Atau2[tuple([])] = [(0,1,14,15,16), (2,7,9,13,18)]
Atau2[tuple([26])] = [(4, 6, 8, 10, 26), (3, 6, 10, 11, 26), (3, 7, 9, 11, 26)]
Atau2[tuple([4])] = [(0, 4, 8, 15, 20), (4, 6, 7, 11, 22)]
Atau2[(3,4)] = [(3, 4, 7, 11, 15), (3, 4, 5, 8, 10)]
Atau2[(4,5)]  = [(3, 4, 5, 15, 23), (4, 5, 6, 14, 28)]
Atau2[(4,26)] = [(4, 7, 8, 9, 26), (4, 6, 7, 11, 26), (4, 6, 8, 10, 26)]
Atau2[(4,5,9)] = [(3, 4, 5, 7, 9), (4, 5, 9, 14, 25)]
Atau2[(4,9,26)] = [(4, 7, 8, 9, 26), (4, 6, 8, 9, 26), (3, 4, 7, 9, 26)]
Atau2[(4,6,26)] = [(4, 6, 7, 11, 26), (4, 6, 8, 10, 26)]
Atau2[(3,4,5)] = [(3, 4, 5, 7, 9), (3, 4, 5, 6, 25)]
Atau2[(3,4,26)] = [(3, 4, 7, 11, 26), (3, 4, 8, 10, 26)]
Atau2[(4,7,8,9)] = [(4, 7, 8, 9, 26), (4, 7, 8, 9, 27)]
Atau2[(3,4,5,25)] = [(3, 4, 5, 6, 25), (3, 4, 5, 9, 25)]
Atau2[(3,4,5,6)] = [(3, 4, 5, 6, 10), (3, 4, 5, 6, 11)]
Atau2[(3,4,6,26)] = [(3, 4, 6, 10, 26), (3, 4, 6, 11, 26)]	


# verification that the dictionaries Atau and Atau2 give the same cones. 
for taui in Atau1.keys():
    C1s = list(map(lambda x: sigmaTwistCone(x[2], x[0], x[1]), Atau1[taui]))
    C2s = list(map(lambda y: indexToCone(y), Atau2[taui]))
    [C1s[i] == C2s[i] for i in range(len(C1s))]
    
# Test to see if Lemma 7.2 is satisfied for tau and Atau.
def testAtau(tau,AtauV):
    matricesR = map( lambda s: coneMatrix(raysS5, linealityS5, s), AtauV )
    matricesS = map( lambda R: kernel(transpose(R)).matrix(), matricesR  )
    S = block_matrix( [[ transpose(Si) for Si in matricesS ]] , subdivide = False  )
    tauM = coneMatrix(raysS5, linealityS5, tau)
    return 16-rank(tauM) == rank(S)

# verification of Lemma 7.8.
list(map(lambda x: testAtau(x,Atau2[x]), Atau2.keys()))
