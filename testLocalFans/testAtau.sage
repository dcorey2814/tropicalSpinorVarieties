with open("raysS5.dat" ,'r') as rF:
    rFLines = rF.readlines()
    raysS5s = list(map(lambda x: vector(map(lambda y: int(y), x.split(" "))), rFLines))

with open("linealityS5.dat" ,'r') as lF:
    lFLines = lF.readlines()
    linealityS5 = list(map(lambda x: vector(map(lambda y: int(y), x.split(" "))), lFLines))

with open("coneOrbitsS5.dat" ,'r') as coF:
    coFLines = coF.readlines()
    coneOrbitsS5s = list(map(lambda x: tuple(map(lambda y: int(y), x.split(" "))), coFLines))


def coneSpan(rays, lineality, cone):
    return span(QQ, [rays[i] for i in cone]+lineality)

def coneMatrix(rays, lineality, cone):
    return matrix(QQ, [rays[i] for i in cone]+lineality)


def testAtau(tau,AtauV):
    matricesR = map( lambda s: coneMatrix(raysS5s, linealityS5, s), AtauV )
    matricesS = map( lambda R: kernel(transpose(R)).matrix(), matricesR  )
    S = block_matrix( [[ transpose(Si) for Si in matricesS ]] , subdivide = False  )
    tauM = coneMatrix(raysS5s, linealityS5, tau)
    return 16-rank(tauM) == rank(S)

    
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


list(map(lambda x: testAtau(x,Atau[x]), Atau.keys()))


