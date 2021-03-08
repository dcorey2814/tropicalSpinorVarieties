needsPackage "Polyhedra"
needsPackage "Matroids"
needsPackage "Graphs"
needsPackage "SpechtModule"


stepwiseSaturate = (I,g) -> (
    for i from 0 to #g-1 do (
        I=saturate(I,sub(g#i, ring I), Strategy=>Bayer);
	);
    return I
    )
-- I is an ideal of a polynomial ring.
-- g is a list of homogeneous elements of ring I.
-- returns the saturation of I with respect to the elements in g.


stepwiseSaturateNoBayer = (I,g) -> (
    for i from 0 to #g-1 do (
        I=saturate(I,sub(g#i, ring I));
        );
    return I
    )
-- I is an ideal of a polynomial ring.
-- g is a list of homogeneous elements of ring I.
-- returns the saturation of I with respect to the elements in g.



inw=(w,I)->(
    R := ring I;
    Rw := newRing(R, Weights=>w, Global => false);
    inwI := sub(ideal(leadTerm(1, sub(I,Rw))),R);
    return stepwiseSaturate(inwI, gens R)
    )
--w is the vector in TV(I), in the min convention, revLex ordered.
--I is the ideal of the thin Schubert cell.



-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
-- The following functions use Pluecker coordinates to: create the ideals of thin schubert cells and limits of thin Schubert cells
-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------

evenSubsets = n-> (
    k := n//2;
    return flatten( for i from 0 to k list subsets(n,2*i))
    )

oddSubsets = n-> (
    k := n//2;
    return flatten( for i from 0 to k list subsets(n,2*i+1))
    )

ringSpin = (n,F) -> (
    eS := {};
    if n%2==0 then eS=evenSubsets(n) else eS=oddSubsets(n) ;
    F[apply(eS, s-> value("p"|concatenate(apply(s, j->toString(j)))|"e"))]
    )


sublistToCoord = (s,R) -> (
    if s == {} then return sub(value("pe"),R) else if #s != #(set(s)) then return 0_R else  (
        sortS := sort(s);
	zeroIndex := hashTable for i from 0 to #sortS-1 list {sortS#i,i};
	zeroIndexs := apply(s, i->zeroIndex#i); 
	sign := permutationSign(zeroIndexs);
    	return sign*sub(value("p"|concatenate(apply(sortS, j->toString(j)))|"e"),R)
	)
    )



Pmunu = (mu, nu, n, F) -> (
    rS := ringSpin(n,F); 
    s1 := sum(for i from 0 to #nu-1 list ((-1)^i)*sublistToCoord(join({nu#i}, mu),rS)*sublistToCoord(drop(nu,{i,i}),rS));
    s2 := sum(for j from 0 to #mu-1 list ((-1)^j)*sublistToCoord(drop(mu,{j,j}),rS)*sublistToCoord(join({mu#j}, nu),rS));
    return sub(s1+s2, rS)
    )


spin = (n,R) -> (
    oddn = {};
    if n%2==0 then oddn = oddSubsets(n) else oddn = evenSubsets(n);
    F=coefficientRing(R);
    ideal(mingens ideal( flatten apply(oddn, mu -> apply(oddn, nu -> sub(Pmunu(mu,nu,n,F), R)))))
    )



sOnList = (s,l) -> sort apply(l, i->s#i)

sOnListSign = (s,l) -> (
    sl := apply(l, i->s#i);
    sortS := sort(sl);
    zeroIndex := hashTable for i from 0 to #sortS-1 list {sortS#i,i};
    zeroIndexs := apply(sl, i->zeroIndex#i); 
    return permutationSign(zeroIndexs)
    )


sOnCoords = (s,n) -> (
    eS := {};
    if n%2==0 then eS= evenSubsets(n) else eS=oddSubsets(n);
    return apply(eS, l -> position(eS, i->i==sOnList(s,l)))
    )

sOnCoordsSigns = (s,n) -> (
    eS := {};
    if n%2==0 then eS= evenSubsets(n) else eS=oddSubsets(n);
    return apply(eS, l -> sOnListSign(s,l))
    )





--fix
TSC = (Bs,n,R) -> (
    F=coefficientRing(R);
    rS := ringSpin(n,F);
    spinn := spin(n,rS); 
    nBs := toList(set(evenSubsets(n)) + set({{}}) - set(Bs)); 
    output := spinn; 
    if #nBs != 0 then (
	pBases := apply(Bs, b-> sublistToCoord(b,rS));
	pNonBases := apply(nBs, nb -> sublistToCoord(nb,rS));
    	SM := eliminate(spinn + ideal(pNonBases), pNonBases);
    	output = sub(stepwiseSaturate(SM,pBases), F[pBases]); 
	);
    return output
    )
-- Input: a matroid S
-- output: the ideal of the thin Schubert cell. Note that we saturate with respect to the 
-- pijk where ijk is a basis of S




--fix
TSCpreSaturate = (Bs,n,R) -> (
    F=coefficientRing(R);
    rS := ringSpin(n,F);
    spinn := spin(n,rS);
    nBs := toList(set(evenSubsets(n)) + set({{}}) - set(Bs));
    output := spinn; 
    if #nBs != 0 then (
	pBases := apply(Bs, b-> sublistToCoord(b,rS));
	pNonBases := apply(nBs, nb -> sublistToCoord(nb,rS));
    	SM := eliminate(spinn + ideal(pNonBases), pNonBases);
    	output = sub(SM, F[pBases]); 
	);
    return output
    )

-- Input: a matroid S
-- output: the ideal of the thin Schubert cell. Note that we do not saturate with respect to the 
-- pijk where ijk is a basis of S



--fix
limitTSC = (SD, n, R) -> (
    idealsPreSat := SD / (Si -> TSCpreSaturate(Si,n,R));
    idealsPreSat = idealsPreSat / (I -> sub(I,R));
    idealPreSat := sum idealsPreSat;
    return stepwiseSaturate(idealPreSat, gens R)
    )
-- input: ringS - the polynomial ring where Gr_S lives and SD - a list of matroids representing a 
-- subdivision of Delta_S
-- output: the ideal of ringS generated by the ideals of the thin Schubert cells from the matroids in SD
-- this is the limit of the thin Schubert cells over the subdivision.



-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
-- The following functions use gfan to compute the tropicalizations of thin Schubert cells. These functions 
-- use gfan0.5 for all tropical computations
-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------


makeInitialFormsInFile = (R,M,w) -> (  
    foutName := temporaryFileName () | ".txt";
    fout := openOut(foutName);
    fout << substring(toString(R), 1, #(toString(R)) -1) << endl;
    fout << toString(M) << endl; 
    if #w > 1 then (
	fout << toString(toSequence(w)) << endl; 
	) else (
	fout << "("|toString(w#0)|")";
	);
    close fout;
    return foutName
    )
-- Input: R a ring, M a list of elements of R, w a list of length (gens R).
-- Output: makes temporary file to compute the initial ideal generated by M wrt w using gfan_initialforms. Uses max convention. Returns the name of the file.



makeTropicalStartingConeFile = (R,M) -> (
    foutName := temporaryFileName () | ".txt";
    fout := openOut(foutName); 
    fout << "Q[";
    fout << substring(toString(gens R), 1, #(toString(gens R)) -2);
    fout << "]" << endl;
    fout << toString(M) << endl; 
    close fout;
    return foutName
    )
-- Input: R a ring, M a list of elements of R generating  a prime ideal.
-- Output: makes temporary file to compute starting cone using gfan_tropicalstartingcone. Returns name of file. 



    
makeTropicalTraverseFile = (startFile, sym, signs) -> (
    fout := openOutAppend(startFile); 
    fout << toString(sym) << endl;
    fout << toString(signs) << endl;
    close fout;
    return startFile
    )
-- Input: startFile a file to use for gfan_tropicaltraverse, sym a list of symmetries of ideal, signs a list of signs of the symmetries, see gfan documentation. 
-- Output: makes temporary file to compute tropicalization using gfan_tropicaltraverse. Returns name of file. 


makeGroebnerConeFile = (R,M,w) -> (
    inFormsFile := makeInitialFormsInFile(R,M,w); 
    outInFormsFile :=  temporaryFileName () | ".txt";
    grobConeFile := temporaryFileName () | ".txt"; 
    run ("gfan_initialforms --ideal --pair < "|inFormsFile|" > "|outInFormsFile);
    run ("gfan_groebnercone --pair --asfan < "|outInFormsFile|" > "|grobConeFile);
    return grobConeFile
    )



gfanRay = (r,lr) -> (
    tr := separate("#", r);
    tr2 := separate(" ", tr#0);
    for i in (0..lr-1) list value(tr2#i)
    )
-- Input: r a string coming from gfan_tropicaltraverse output representing a RAY, lr an integer representing the ambient dimension of r. 
-- Output: a list of integers representing the ray r.



gfanRaysToM2 = (rs, lr) -> for i from 0 to #rs-1 list gfanRay(rs#i, lr)
-- Input: rs a list of strings r as in gfanRay, lr an integer representing the ambient dimension of each r.
-- Output: a list of lists of integers, representing the rays of the tropicalization. 

gfanLinealityToM2 = (ls, d) -> (
    tls := for i from 0 to #ls-1 list separate(" ", ls#i);
    for i from 0 to #ls-1 list( for j from 0 to d-1 list value(tls#i#j)  )
    )


gfanConeLineToList = l -> (
    Cstr := (separate("#", l))#0;
    rightBracket := position(characters Cstr, i->i== "}");
    CstrNoBrackets := substring(Cstr, 1, rightBracket-1);
    if #CstrNoBrackets == 0 then return {} else return apply(separate(" ", CstrNoBrackets), i->value(i))
    )

parseTropFile = f -> (
     lf := lines get f;
     T := new MutableHashTable;
     ambDimLoc := position(lf, i-> i == "AMBIENT_DIM"); 
     dimLoc := position(lf, i-> i == "DIM"); 
     nRaysLoc := position(lf, i-> i == "N_RAYS"); 
     linDimLoc := position(lf, i-> i == "LINEALITY_DIM");
     
     T#"ambDim" = value(lf#(ambDimLoc+1));
     T#"dim" = value(lf#(dimLoc+1)); 
     T#"linealityDim" = value(lf#(linDimLoc+1)); 
     T#"nRays" = value(lf#(nRaysLoc+1));
     if T#"nRays"==0 then T#"rays" = {} else (
	 
         raysLoc := position(lf, i-> i== "RAYS"); 
 	 conesLoc := position(lf, i-> i=="CONES_ORBITS"); 
     	 maxConesLoc := position(lf, i-> i=="MAXIMAL_CONES_ORBITS");
    	 multLoc := position(lf, i-> i == "MULTIPLICITIES_ORBITS");
	 conesLines := for i in (conesLoc+1..maxConesLoc-2) list lf#i;
         maxConesLines := for i in (maxConesLoc+1..multLoc-2) list lf#i; 


         T#"rays" = gfanRaysToM2(for i from raysLoc+1 to raysLoc+T#"nRays" list lf#i, T#"ambDim"); 
	 T#"conesOrbits" = for i in (0..#conesLines-1) list gfanConeLineToList(conesLines#i); 
         T#"relativeInteriorVectors" = for i in (1..#(T#"conesOrbits")-1) list sum(for j in T#"conesOrbits"#i list T#"rays"#j);    
	 T#"maximalConesOrbits" = for i in (0..#maxConesLines-1) list gfanConeLineToList(maxConesLines#i); 
	 ); 
     
     if T#"linealityDim" == 0 then T#"lineality" = {} else (
	 linLoc := position(lf, i-> i=="LINEALITY_SPACE"); 
	 T#"lineality" = gfanLinealityToM2(for i from linLoc+1 to linLoc+T#"linealityDim" list lf#i, T#"ambDim");
	 ); 
     
     return T
    )

-- Input: file f, the output of gfan_tropicaltraverse
-- Output: a MutableHashTable T:
--    T#"ambDim" = integer, the ambient dimension of the tropicalization
--    T#"dim" = integer, the dimension of the tropicalization
--    T#"linealityDim" = integer, dimension of the lineality space
--    T#"nRays" = integer, number of rays
--    T#"rays" = list of lists of integers, the rays of the tropicalization
--    T#"conesOrbits" = list of lists of integers,  the cones of the tropicalization (as lists of indices of rays)
--    T#"relativeInteriorVectors" = list of lists of integers, the sum of all rays of a given cone, this is a vector in the relative interior of the cone. 
--    T#"lineality" = list of lists of integers, span of these vectors form the lineality space. 
 


parseGroebnerConeFile = f -> (
     lf := lines get f;
     GC := new MutableHashTable;
     
     ambDimLoc := position(lf, i-> i == "AMBIENT_DIM");
     dimLoc := position(lf, i-> i == "DIM");
     nRaysLoc := position(lf, i-> i== "N_RAYS");
     linDimLoc := position(lf, i-> i=="LINEALITY_DIM");
     
     GC#"ambDim" = value(lf#(ambDimLoc+1));
     GC#"dim" = value(lf#(dimLoc+1));
     GC#"linealityDim" = value(lf#(linDimLoc+1));
     GC#"nRays" = value(lf#(nRaysLoc+1));
    
     
     if GC#"nRays"==0 then GC#"rays" = {} else (
	 raysLoc := position(lf, i-> i== "RAYS");
         GC#"rays" = gfanRaysToM2(for i from raysLoc+1 to raysLoc+GC#"nRays" list lf#i, GC#"ambDim");
         );
     if GC#"linealityDim" == 0 then GC#"lineality" = {} else (
	 linLoc := position(lf, i-> i=="LINEALITY_SPACE"); 
	 GC#"lineality" = gfanLinealityToM2(for i from linLoc+1 to linLoc+GC#"linealityDim" list lf#i, GC#"ambDim");
	 ); 
     
     return GC
    )


groebnerCone = (R, M, w) -> (
    f:=makeGroebnerConeFile(R,M,w); 
    d:=parseGroebnerConeFile(f);
    outp := false;
    if (d#"nRays" != 0 and d#"linealityDim" != 0) then (
	outp = coneFromVData(transpose matrix d#"rays", transpose matrix d#"lineality");
	) else if (d#"nRays" != 0 and d#"linealityDim" == 0 ) then (
	outp = coneFromVData(transpose matrix d#"rays");
	) else if (d#"nRays" == 0 and d#"linealityDim" != 0 ) then (
	z := for i from 0 to d#"ambDim"-1 list 0; 
	outp = coneFromVData(transpose matrix {z}, transpose matrix d#"lineality");
	) else outp = false;
    return outp
    )



groebnerConeData = (R, M, w) -> (
    f:=makeGroebnerConeFile(R,M,w); 
    return parseGroebnerConeFile(f)
    )



tropicalize = (I, sym, signs) -> (
    inFileName := makeTropicalStartingConeFile(ring I, I_*); 
    outFileName := temporaryFileName () | ".txt"; 
    tropFileName :=  temporaryFileName () | ".txt"; 
    run ("gfan_tropicalstartingcone < "|inFileName|" > "|outFileName);
    makeTropicalTraverseFile(outFileName, sym, signs);
    run ("gfan_tropicaltraverse --symmetry --symsigns --nocones < "|outFileName|" > "|tropFileName);
    T := parseTropFile(tropFileName);
    return T
    )
-- Input: GrS an ideal of a thin schubert cell (should be prime), sym list of sequences, signs list of sequences. These are for --symmetry and --symsigns for gfan_tropicaltraverse.
-- Output: a MutableHashTable from parseTropFile

tropicalizeWithVector = (I, w, sym, signs) -> (
    inFileName := makeInitialFormsInFile(ring I, I_*, w);
    outFileName := temporaryFileName () | ".txt";
    tropFileName :=  temporaryFileName () | ".txt";
    run ("gfan_initialforms --ideal --pair < "|inFileName|" > "|outFileName);
    makeTropicalTraverseFile(outFileName, sym, signs);
    run ("gfan_tropicaltraverse --symmetry --symsigns --nocones< "|outFileName|" > "|tropFileName);
    T := parseTropFile(tropFileName);
    return T
    )
-- Input: GrS an ideal of a thin schubert cell (should be prime), w a list of integers (should be in the relative interior of a maximal cone), sym list of sequences, signs list of sequences. These are for --symmetry and --symsigns for gfan_tropicaltraverse.
-- Output: a MutableHashTable from parseTropFile





Ry=QQ[y0,y1,y2,y3,z0,z1,z2,z3]
IF=ideal(-y0-y1-y2+z3, -y0-y1+y3+z2, -y0+y2+y3+z1, y1+y2+y3+z0, 
         y0-z1+z2-z3, y1+z0-z2+z3, y2-z0+z1-z3, y3+z0-z1+z2)

symF = {{4,5,6,7,0,1,2,3},{1,2,3,0,5,6,7,4}, {1,0,2,3,5,4,6,7}}
signsF = {{-1,1,-1,1,-1,1,-1,1}, {1,1,1,-1,1,1,1,-1}, {-1,-1,1,1,1,1,1,1}}

symF = {{4,5,6,7,0,1,2,3},{1,2,3,0,5,6,7,4}}
signsF = {{-1,1,-1,1,-1,1,-1,1}, {1,1,1,-1,1,1,1,-1}}
TF = tropicalize(IF, symF, signsF)

TF#"rays"


tropicalize(IF2,{{0,1,2,3,4,5,6,7}},{{1,1,1,1,1,1,1,1}})

dim(IF)
ring IF
IF_*

inw(TF#"relativeInteriorVectors"#0,IF)
inw(TF#"relativeInteriorVectors"#1,IF)
inw(TF#"relativeInteriorVectors"#2,IF)

Mt = matrix{{1,0,0,0,0,1,1,1},{0,1,0,0,-1,0,1,1},{0,0,1,0,-1,-1,0,1},{0,0,0,1,-1,-1,-1,0}}
MatF = matroid(Mt)
circuits(MatF)

inw(  , IF)
