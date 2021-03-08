import itertools as itt
import time


with open("raysS5gfanv2.dat" ,'r') as rF:
    rFLines = rF.readlines()
    raysS5 = list(map(lambda x: vector(map(lambda y: int(y), x.split(" "))), rFLines))

    
with open("maxConesS5gfanv2.dat" ,'r') as mcF:
    mcFLines = mcF.readlines()
    maxConesS5 = list(map(lambda x: tuple(map(lambda y: int(y), x.split(" "))), mcFLines))

with open("coneOrbitsS5gfanv2.dat" ,'r') as coF:
    coFLines = coF.readlines()
    coneOrbitsS5 = list(map(lambda x: tuple(map(lambda y: int(y), x.split(" "))), coFLines))


with open("raysS5.dat" ,'r') as rF:
    rFLines = rF.readlines()
    raysS5s = list(map(lambda x: vector(map(lambda y: int(y), x.split(" "))), rFLines))

with open("linealityS5.dat" ,'r') as lF:
    lFLines = lF.readlines()
    linealityS5 = list(map(lambda x: vector(map(lambda y: int(y), x.split(" "))), lFLines))

with open("coneOrbitsS5.dat" ,'r') as coF:
    coFLines = coF.readlines()
    coneOrbitsS5s = list(map(lambda x: tuple(map(lambda y: int(y), x.split(" "))), coFLines))


with open("symE5v2.dat",'r') as sF:
    sFLines = sF.readlines()
    symE5 = list(map(lambda x: list(map(lambda y: int(y), x.split(" "))), sFLines))


oldToNewRayIndex = {
0:16, 1:17, 2:20, 3:35, 4:18, 5:21, 6:34, 7:23,
8:32, 9:29, 10:19, 11:22, 12:33, 13:24, 14:31, 15:28,
16:25, 17:30, 18:27, 19:26, 20:0, 21:15, 22:14, 23:13,
24:11, 25:8, 26:1, 27:12, 28:10, 29:7, 30:2, 31:9,
32:6, 33:3, 34:5, 35:4
}

newToOldRayIndex = {v: k for k, v in oldToNewRayIndex.items()}


def oldToNewConeIndex(oldc):
    return tuple(sorted([oldToNewRayIndex[i] for i in oldc]))

def newToOldConeIndex(newc):
    return tuple(sorted([newToOldRayIndex[i] for i in newc]))



maxConesNew = list(map( lambda c: oldToNewConeIndex(c) , maxConesS5))

starCone = {};
for c in coneOrbitsS5s:
    starCone[c] = [mc for mc in maxConesNew if set(c).issubset(set(mc))]

def makeStar(c):
    [mc for mc in maxConesNew if set(c).issubset(set(mc))]


def coneSpan(rays, lineality, cone):
    return span(QQ, [rays[i] for i in cone]+lineality)

def coneMatrix(rays, lineality, cone):
    return matrix(QQ, [rays[i] for i in cone]+lineality)


def testAtau(tau,Atau):
    matricesR = map( lambda s: coneMatrix(raysS5s, linealityS5, s), Atau  )
    matricesS = map( lambda R: kernel(transpose(R)).matrix(), matricesR  )
    S = block_matrix( [[ transpose(Si) for Si in matricesS ]] , subdivide = False  )
    tauM = coneMatrix(raysS5s, linealityS5, tau)
    return 16-rank(tauM) == rank(S)



maxCone2OrbitRepOld = {}
for i in range(0,120):
    maxCone2OrbitRepOld[maxConesS5[i]] = maxConesS5[119]
for i in range(120,312):
    maxCone2OrbitRepOld[maxConesS5[i]] = maxConesS5[311]
for i in range(312,720):
    maxCone2OrbitRepOld[maxConesS5[i]] = maxConesS5[719]
for i in range(720,912):
    maxCone2OrbitRepOld[maxConesS5[i]] = maxConesS5[911]

maxCone2OrbitRepNew = {}
for mc in maxConesNew:
    oldmc = newToOldConeIndex(mc); 
    oldmcOrb = maxCone2OrbitRepOld[oldmc];
    maxCone2OrbitRepNew[mc] = oldToNewConeIndex(oldmcOrb)
    
Atau = {}
Atau[tuple([])] = [(0,1,14,15,16), (2,7,9,13,18)]
Atau[tuple([26])] = [(3, 6, 10, 11, 26), (3, 7, 9, 11, 26), (4, 6, 8, 10, 26) ]
Atau[tuple([4])] = [(0, 4, 8, 15, 20), (4, 6, 7, 11, 22) ]
Atau[(3,4)] = [(3, 4, 7, 11, 15), (3, 4, 5, 8, 10) ]
Atau[(4,5)]  = [(3, 4, 5, 15, 23), (4, 5, 6, 14, 29)]
Atau[(4,26)] = [(4, 6, 7, 11, 26), (4, 6, 8, 10, 26), (4, 7, 8, 9, 26) ]
Atau[(4,5,9)] = [(4, 5, 8, 9, 14), (4, 5, 9, 14, 25)]
Atau[(4,9,26)] = [(4, 7, 8, 9, 26), (4, 6, 8, 9, 26), (4, 6, 7, 9, 26) ]
Atau[(4,6,26)] = [(4, 6, 7, 11, 26), (4, 6, 8, 10, 26)  ]
Atau[(3,4,5)] = [(3, 4, 5, 6, 25), (3, 4, 5, 7, 9)]
Atau[(3,4,26)] = [(3, 4, 6, 9, 26), (3, 4, 7, 11, 26)]
Atau[(4,7,8,9)] = [(4, 7, 8, 9, 20), (4, 7, 8, 9, 26)]
Atau[(3,4,5,25)] = [(3, 4, 5, 6, 25), (3, 4, 5, 9, 25)]
Atau[(3,4,5,6)] = [(3, 4, 5, 6, 10), (3, 4, 5, 6, 11)]
Atau[(3,4,6,26)] = [(3, 4, 6, 9, 26), (3, 4, 6, 10, 26)]	





list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[tuple([])] ))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[tuple([26])] ))
list(map( lambda mc: maxCone2OrbitRepNew[mc] ,Atau[tuple([4])]  ))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[(3,4)] ))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[(4,5)] ))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[(4,26)] ))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[(4,5,9)] ))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[(4,9,26)]))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[(4,6,26)]))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[(3,4,5)]))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[(3,4,26)]))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[(4,7,8,9)] ))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[(3,4,5,25)] ))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[(3,4,5,6)] ))
list(map( lambda mc: maxCone2OrbitRepNew[mc] , Atau[(3,4,6,26)] ))






def indexToCone(c, rs, lin):
    return Polyhedron(rays = [vector(rs[i]) for i in c], lines = lin, base_ring=QQ)



def checkRays(i):
    oldRay = i
    newRay = oldToNewRayIndex[oldRay]
    C1 = indexToCone([oldRay], raysS5, linealityS5 )
    C2 = indexToCone([newRay], raysS5s, linealityS5 )
    return -C1==C2


indexToCone([0], raysS5, linealityS5 ).rays()
indexToCone([16], raysS5s, linealityS5 ).rays()








[b for b in a for a in *list(Atau.values())]
*list(Atau.values())

for k in Atau.keys():
    if not testAtau(k,Atau[k]):
        print(k)



list(set.union(*[set(s) for s in Atau.values()]))



















'''
rayToIndex = {}
for i in range(len(raysS5s)):
    rayToIndex[tuple(raysS5s[i])] = i


def sOnVec(s,v):
    return vector([v[s[i]] for i in range(len(s))])


def sOnRayIndex(s,ri):
    r = raysS5s[ri]
    sr = sOnVec(s,r)
    return rayToIndex[tuple(sr)]

def sOnConeIndex(s,c):
    return tuple(sorted(list( map(lambda y: sOnRayIndex(s,y), list(c) ))))



maxCones = {}
for i in [14,15,16,17]:
    new = set(map( lambda s: sOnConeIndex(s,coneOrbitsS5s[i]), symE5  ))
    print(new)
    maxCones = maxCones + new
    

sOnConeIndex(s,c)

'''























































def intersectLinearSpansRelDim(starDict, rays, lineality, sigma):
    sigmaSpan = coneSpan(rays,lineality,sigma)
    sigmaSpanDim = sigmaSpan.dimension()
    maxCones = starDict[sigma]
    intersectConesSpan = coneSpan(rays,lineality,maxCones[0])
    i=0
    while intersectConesSpan.dimension() > sigmaSpanDim and i<len(maxCones):
        newConeSpan = coneSpan(rays,lineality,maxCones[i])
        intersectConesSpan = intersectConesSpan.intersection(newConeSpan)
        i+=1
    return intersectConesSpan.dimension() - sigmaSpanDim


starCone = {};
for c in coneOrbitsS5:
    starCone[c] = [mc for mc in maxConesS5 if set(c).issubset(set(mc))]



def indexToCone(c, rs, lin):
    return Polyhedron(rays = [vector(rs[i]) for i in c], lines = lin, base_ring=QQ)


Cone1 = indexToCone((11,12,14,15,32), raysS5, linealityS5)
Cone2 = indexToCone((11,12,14,15,35), raysS5, linealityS5)
Cone3 = indexToCone((32,35), raysS5, linealityS5)

C13 = Cone1.intersection(Cone3)
    
Cone1s = indexToCone((3,4,6,9,25), raysS5s, linealityS5)
Cone2s = indexToCone((3,4,6,9,26), raysS5s, linealityS5)


allTests = []
for j in [1,2,3,4,5]: 
    cones_j = secondaryRepsDim[j] 
    for sigma in cones_j:
        allTests.append(intersectLinearSpansRelDim(starsDimensionDict[j],rays37,lineality37,sigma))
print(allTests)



def groupFromGenerators(n,gens):
# Input: a positive integer "n", and a set "gens" of tuples, each tuple is a permutation of (0,1,...,n-1). 
# Output: a list of tuples, each tuple is a permutation of (0,1,...,n-1). These are the elements of S_n generated by the generators "gens" 
    def one_iteration(n,gens,to_act,group,exhausted):
        if to_act == set():
            return list(group)
        else:
            new=set.union(*[set([tuple([g[s[i]] for i in range(n)]) for s in to_act]) for g in gens])
            group.update(new)
            exhausted.update(to_act)
            to_act = new - exhausted
            return one_iteration(n,gens,to_act,group,exhausted)

    return one_iteration(n,gens,gens.copy(),gens.copy(),gens.copy())

def g_inv(g):
# Input: a tuple "g" which is a permutation of (0,1,...,n-1), for n=len(g).
# Output: a tuple that is the inverse of the permutation g.
    inv={g[i]:i for i in range(len(g))}
    return tuple([inv[gi] for gi in range(len(g))])


def s_on_ijk(s,triple):
# Input: a permutation "s"  and a tuple "triple" (i,j,k) representing a basis of a rank 3 matroid.
# Output: a new tuple permuted by s, namely (s[i],s[j],s[k]) (reordered in increasing order).
    return tuple(sorted([s[triple[0]],s[triple[1]],s[triple[2]]]))

def s_on_ray(s,ray):
# Input: a ray "ray" and a tuple "s" which is a permutation on (0,1,...,len(ray)-1).
# Output: new ray with coordinates permuted by s.
    g=g_inv(s)
    return tuple([ray[g[i]] for i in range(len(ray))])

def s_to_NBases(n,Bases,s):
# Input: a positive integer n (the size of the base set of the matroid), the list of "Bases" of a rank 3 matroid on [n], and a tuple "s" representing a permutation 
#  of (0,1,...,n-1).
# Output: a tuple which is a permutation of (0,1,...,len(Bases)-1) given by the induced permutation of s on the list Bases. More specifically, "Bases" is originally ordered in revlex, so this corresponds to the permutation (0,1,...,len(Bases)-1).  Then we use "s_on_ijk" to let s act on each basis, and record its original position in revlex. This gives a permutation of (0,1,...,len(Bases)-1).
    Bases_index={Bases[i]:i for i in range(len(Bases))}

    def s_on_NBases(s,i):
    # Input: and a tuple "s" representing a permutation of (0,1,...,len(Bases)-1).
    # Output: the position of (s[i],s[j],s[k]) in "Bases".
        p_ijk = Bases[i]
        p_sijk= s_on_ijk(s,p_ijk)
        return Bases_index[p_sijk]

    return tuple([s_on_NBases(s,i) for i in range(len(Bases))])


def s_to_Nrays(n,Bases,s,RAYS):
# Input: a positive integer n (the size of the base set of the matroid), the list of "Bases" of a rank 3 matroid on [n], a tuple "s" representing a permutation 
#  of (0,1,...,n-1), and a list "RAYS" of rays, each ray is a tuple.
# Output: a tuple which is a permutation of (0,1,...,len(RAYS)-1) given by the induced permutation of s on the RAYS. More specifically, "RAYS" is a given fixed list, and his corresponds to the permutation (0,1,...,len(Bases)-1).  Then we use "s_on_ray" to let s act on each ray, and record its position in RAYS. This gives a permutation of (0,1,...,len(RAYS)-1).

    sBases=s_to_NBases(n,Bases,s)
    rays_index = {RAYS[i]:i for i in range(len(RAYS))}

    def g_on_ray_index(g,i):
    # Input: and a tuple "g" representing a permutation of (0,1,...,len(RAYS)-1).
    # Output: the position of s_on_ray(g,ray) in "RAYS".
        ray=RAYS[i]
        gray=tuple(s_on_ray(g,ray))
        return rays_index[gray]

    return tuple([g_on_ray_index(sBases,i) for i in range(len(rays))])

def G_in_NRays(n,Bases,gensG,RAYS):
# Input: a positive integer n (the size of the base set of the matroid), the list of "Bases" of a rank 3 matroid on [n], a set "gensG" tuples "s" which are permutations of (0,1,...,n-1) (This is a set of generators of a group G), and a list "RAYS" of rays, each ray is a tuple.

# Output: a list of tuples giving a group "GROUP." This is a subgroup of the permutations on (0,1,...,len(RAYS)-1) generated by the elements of "gensG" where each such generator is converted to a permutation on (0,1,...,len(RAYS)-1).

    G=groupFromGenerators(n,gensG)

    Bases_index={Bases[i]:i for i in range(len(Bases))}
    rays_index={RAYS[i]:i for i in range(len(RAYS))}

    def s_on_NBases(s,i):
    # Input: and a tuple "s" representing a permutation of (0,1,...,len(Bases)-1).
    # Output: the position of (s[i],s[j],s[k]) in (0,1,...,len(Bases)-1) with respect to revLex.
        p_ijk = Bases[i]
        p_sijk= s_on_ijk(s,p_ijk)
        return Bases_index[p_sijk]

    NBases=len(Bases)
    G_in_NBases=[tuple([s_on_NBases(s,i) for i in range(NBases)]) for s in G]

    def g_on_ray_index(g,i):
    # Input: and a tuple "g" representing a permutation of (0,1,...,len(RAYS)-1).
    # Output: the position of s_on_ray(g,ray) in "RAYS".
        ray=RAYS[i]
        gray=tuple(s_on_ray(g,ray))
        return rays_index[gray]

    GROUP=[tuple([g_on_ray_index(g,i) for i in range(len(RAYS))]) for g in G_in_NBases]

    return GROUP

def g_on_tuple(g,x):
# Input: a tuple "g" that is a permutation of (0,1,...,N-1), and a tuple "x" consisting of numbers from {0,1,...,N-1}, sorted in increasing order.
# Output: a new tuple derived from "x" given by the permutation "g", sorted in increasing order.
    return tuple(sorted([g[i] for i in x]))

def G_orbits(G,X):
# Input: a group "G" as a list of permutations of (0,1,...,N-1) (giving a group under composition), and "X" a set of tuples of range(N), each tuple is ordered in increasing order
# Output: a dictionary, keys = representatives, and values = the orbit of the representative.
    remaining = X.copy()
    orbits = {}
    while len(remaining)>0:
        x=remaining.pop()
        orbits[x]=set([x])
        for g in G:
            gx=g_on_tuple(g,x)
            orbits[x].add(gx)
            remaining.discard(gx)
    return orbits

def act_by_G(G,Xreps):
# Input: a group "G" as a list of permutations of (0,1,...,N-1) (giving a group under composition), and "Xreps" a set of tuples of range(N), each tuple is ordered in increasing order. 
# Output: a dictionary, keys = x in Xreps, and values = the orbit of x by the action of G.
    orbits={}
    for x in Xreps:
        orbits[x]=set([g_on_tuple(g,x) for g in G])
    return orbits
