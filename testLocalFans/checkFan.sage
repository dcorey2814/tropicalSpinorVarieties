import itertools as itt
import time


with open("linealityS5.dat" ,'r') as lF:
    lFLines = lF.readlines()
    linealityS5 = list(map(lambda x: vector(map(lambda y: int(y), x.split(" "))), lFLines))

with open("raysS5gfanv2.dat" ,'r') as rF:
    rFLines = rF.readlines()
    raysS5 = list(map(lambda x: vector(map(lambda y: int(y), x.split(" "))), rFLines))

with open("raysS5.dat" ,'r') as rF:
    rFLines = rF.readlines()
    raysS5s = list(map(lambda x: vector(map(lambda y: int(y), x.split(" "))), rFLines))


oldToNewRayIndex = {
0:16, 1:17, 2:20, 3:35, 4:18, 5:21, 6:34, 7:23,
8:32, 9:29, 10:19, 11:22, 12:33, 13:24, 14:31, 15:28,
16:25, 17:30, 18:27, 19:26, 20:0, 21:15, 22:14, 23:13,
24:11, 25:8, 26:1, 27:12, 28:10, 29:7, 30:2, 31:9,
32:6, 33:3, 34:5, 35:4
}



def indexToCone(c, rs, lin):
    return Polyhedron(rays = [vector(rs[i]) for i in c], lines = lin, base_ring=QQ)

def checkRays(i):
    oldRay = i
    newRay = oldToNewRayIndex[oldRay]
    C1 = indexToCone([oldRay], raysS5, linealityS5 )
    C2 = indexToCone([newRay], raysS5s, linealityS5 )
    return -C1==C2



list(map(checkRays, range(36)))
