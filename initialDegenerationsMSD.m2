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
    fout << substring(toString(R), 1, #(toString(R)) -1) << endl;
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
    tr := (separate("\t#", r))#0;
    tr2 := separate(" ", tr); 
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







-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
-- The following functions compute matroid subdivisions of matroid polytopes 
-- by interfacing with polymake 3.2. 
-------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------


basisToPolymakeVector = (b,n) -> (
    str := toString({1}|for i in (0..n-1) list (if any(toList b, x-> x ==  i) then 1 else 0));
    return ("["|substring(str, 1, #str-2)|"]")
    )
    
makePolyFile = (Bs,n,w) -> (
    polyFileName := temporaryFileName () | ".poly";
    f := openOut(polyFileName);
    f << "use application 'polytope';" << endl; 
    f << ("my $w = -new Vector<Rational>(["|substring(toString(w), 1, #(toString(w))-2) |"]);") << endl; 
    f << ("my $p = new Polytope(POINTS => ["|demark(", ", for i in (0..#Bs-1) list basisToPolymakeVector(Bs#i,n))|"]);") << endl;
    f << "my $msd = new fan::SubdivisionOfPoints(POINTS=>$p->VERTICES, WEIGHTS => $w);" << endl; 
    f << "print $msd -> N_MAXIMAL_CELLS;";
    f << "print '\n';";
    f << "print $msd -> MAXIMAL_CELLS;";
    f << "print $msd -> POLYHEDRAL_COMPLEX -> DUAL_GRAPH -> EDGES";
    close f;
    return polyFileName
    )


matroidSubdivision =  (Bs,n,w) -> (
    f := makePolyFile(Bs,n,w);    
    o := temporaryFileName () | ".txt";
    run ("polymake --script "|f|">"|o);
    msdStrs := lines get o; 
    NMaxCones := value(msdStrs#0);
    	
    msdLines := for i in (1..#msdStrs-1) list ( apply(separate(" ", substring(msdStrs#i, 1, #(msdStrs#i)-2)), x->value(x)));
    msdIndices := for i in (0..NMaxCones-1) list msdLines#i;
    msdEdges := for i in (NMaxCones..#msdLines-1) list msdLines#i;
           
    basesS :=  Bs / ( x-> sort toList x); 
    msdBases := for i in (0..#msdIndices-1) list (for j in msdIndices#i list basesS#j);
    
    return {msdBases, graph(msdEdges)}
    )







-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
-- The following functions are used to determine dual graphs to subdivisions and which ones have no leaves
-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------

intersectionSets = l -> (
    x := l#0; 
    (0..#l-1) / (i -> x = x * l#i);
    return x
    )
-- Input: l - a list of sets.
-- Output: the intersection of these sets.



symDiff = (A,B) -> (A - B) + (B-A)

noLeaves = G -> not any(vertexSet(G), v->degree(G,v)==1)
-- Input: a connected graph G
-- Output: true if there are no one-valent vertices, false otherwise.



varToList = (p,R) -> (
    sp := toString(p);
    trimsp := substring((1,#sp-2),sp);
    lp := for i from 0 to #trimsp-1 list value(trimsp#i);
    return lp
    )


twist = (J,p,R) -> (
    lp := varToList(p,R); 
    --clp := set(toList(0..n-1)) - set lp;
    s1 := (-1_R)^(#(set(J)*set(lp))); 
    Jmp := sort toList(set J - set lp); 
    pmJ := sort toList(set lp - set J);
    return s1*sublistToCoord(Jmp|pmJ, R)
    )


twistNoSign = (J,p,R) -> (
    lp := varToList(p,R); 
    Jmp := sort toList(set J - set lp); 
    pmJ := sort toList(set lp - set J);
    return value (toList(factor(sublistToCoord(Jmp|pmJ, R))))#0
    )

twistSign = (J,p,R) -> (
    lp := varToList(p,R); 
    s1 := (-1_R)^(#(set(J)*set(lp)));
    Jmp := sort toList(set J - set lp); 
    pmJ := sort toList(set lp - set J);
    factored := toList(factor(s1*sublistToCoord(Jmp|pmJ, R)));
    --if #factored == 1 then return 1 else return value factored#1
    return 1
    )


twistToPerm = (J,R) -> (
    twistedGensNoSigns := apply(gens R, p -> twistNoSign(J, p, R));
    apply(twistedGensNoSigns, p -> position(gens R, r->r==p))
    )


twistToSigns = (J,R) -> (
    return apply(gens R, p -> twistSign(J, p, R))
    )


S5OnE5 = g -> (
    E5 := oddSubsets(5);
    gE5 := apply(E5, e-> sort apply(e, i->g#i));
    apply(gE5, e-> position(E5,i->i==e) )
    )

S5OnRay = (g,v) -> (
    gE5 := S5OnE5(g);
    return apply(toList(0..#v-1), i->v#(gE5#i))
    )

twistRay = (J,R,v) -> (
    JE5 := twistToPerm(J,R);
    return apply(toList(0..#v-1), i->v#(JE5#i))
    )



end


-- secondary fan






uniformS4 = {{}} | evenSubsets(4)
sp4 = spin(4,ringSpin(4,QQ));

sym4={(0,1,3,2,5,4,6,7),(0,3,5,6,1,2,4,7)} | 
apply(evenSubsets(4), J-> twistToPerm(J,ring sp4));

sgn4={(1,-1,1,1,1,1,1,-1),(1,1,1,1,-1,-1,-1,-1)} | 
apply(evenSubsets(4), J-> twistToSigns(J,ring sp4));

TS4 = tropicalize(spin(4,ringSpin(4,QQ)), sym4, sgn4 );
msd4 = apply(TS4#"relativeInteriorVectors", w -> matroidSubdivision(uniformS4,4,w));

msd4#0
msd4#1
peek TS4





















uniformS5 = oddSubsets(5)

sp5 = spin(5,ringSpin(5,QQ))

sym5={(1, 2, 3, 4, 0, 8, 11, 13, 14, 5, 6, 7, 9, 10, 12, 15),
(1, 0, 2, 3, 4, 5, 6, 8, 7, 9, 11, 10, 13, 12, 14, 15)} | 
apply(evenSubsets(5), J-> twistToPerm(J,ring sp5));


sgn5={(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
      (1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, 1, 1, 1, 1, -1)} | 
apply(evenSubsets(5), J-> twistToSigns(J,ring sp5));



TS5 = tropicalize(spin(5,ringSpin(5,QQ)), sym5, sgn5 );
msd5 = apply(TS5#"relativeInteriorVectors", w -> matroidSubdivision(uniformS5,5,w));




TS5#"conesOrbits"#12
TS5#"conesOrbits"#19

inw(TS5#"relativeInteriorVectors"#11, sp5) == inw(TS5#"relativeInteriorVectors"#18, sp5)

sp5

listToString = l -> (
    s:="";
    for i in (0..#l-1) do (s = s|toString(l#i));
    return s
    )
    
matroidToStrs = bl -> apply(bl, l->listToString(l))



L={
    {1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1},
    {0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1},
    {0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1},
    {0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1},
    {0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},   
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2}
    }


r25 = {0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}
r26 = {-1, -1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
wt={0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0}

r25+r26 -(1/2)*(-L#0-L#1+L#2+L#3+L#4-{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})

(L#0+L#1+L#2+L#3+L#4)-2*L#5 == {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}

r-L#3+{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}
+L#5



(TS5#"rays"#0-L#0-L#1-L#2+L#3+L#4+2*L#5-2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#1-L#0-L#1+L#2-L#3+L#4+2*L#5-2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#2-L#0-L#1+L#2+L#3-L#4+2*L#5-2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#3-L#0-L#1+L#2+L#3+L#4-2*L#5+2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#4-L#0+L#1-L#2-L#3+L#4+2*L#5-2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#5-L#0+L#1-L#2+L#3-L#4+2*L#5-2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#6-L#0+L#1-L#2+L#3+L#4-2*L#5+2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#7-L#0+L#1+L#2-L#3-L#4+2*L#5-2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#8-L#0+L#1+L#2-L#3+L#4-2*L#5+2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#9-L#0+L#1+L#2+L#3-L#4-2*L#5+2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#10+L#0-L#1-L#2-L#3+L#4+2*L#5-2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#11+L#0-L#1-L#2+L#3-L#4+2*L#5-2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#12+L#0-L#1-L#2+L#3+L#4-2*L#5+2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#13+L#0-L#1+L#2-L#3-L#4+2*L#5-2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#14+L#0-L#1+L#2-L#3+L#4-2*L#5+2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#15+L#0-L#1+L#2+L#3-L#4-2*L#5+2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#16+L#0+L#1-L#2-L#3-L#4+2*L#5-2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#17+L#0+L#1-L#2-L#3+L#4-2*L#5+2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#18+L#0+L#1-L#2+L#3-L#4-2*L#5+2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4
(TS5#"rays"#19+L#0+L#1+L#2-L#3-L#4-2*L#5+2*{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1})/4

(TS5#"rays"#20-3*L#0-L#1-L#2-L#3-L#4+4*L#5)/8
(TS5#"rays"#21+L#0+L#1+L#2+L#3+L#4-4*L#5)/8
(TS5#"rays"#22+L#0+L#1-L#2-L#3-L#4)/8
(TS5#"rays"#23+L#0-L#1+L#2-L#3-L#4)/8
(TS5#"rays"#24+L#0-L#1-L#2+L#3-L#4)/8
(TS5#"rays"#25+L#0-L#1-L#2-L#3+L#4)/8
(TS5#"rays"#26-L#0-3*L#1-L#2-L#3-L#4+4*L#5)/8
(TS5#"rays"#27-L#0+L#1+L#2-L#3-L#4)/8
(TS5#"rays"#28-L#0+L#1-L#2+L#3-L#4)/8
(TS5#"rays"#29-L#0+L#1-L#2-L#3+L#4)/8
(TS5#"rays"#30-L#0-L#1-3*L#2-L#3-L#4+4*L#5)/8
(TS5#"rays"#31-L#0-L#1+L#2+L#3-L#4)/8
(TS5#"rays"#32-L#0-L#1+L#2-L#3+L#4)/8
(TS5#"rays"#33-L#0-L#1-L#2-3*L#3-L#4+4*L#5)/8
(TS5#"rays"#34-L#0-L#1-L#2+L#3+L#4)/8
(TS5#"rays"#35-L#0-L#1-L#2-L#3-3*L#4+4*L#5)/8







qs = {"f_{0}","f_{1}","f_{2}","f_{3}","f_{4}",
    "f_{012}","f_{013}","f_{023}","f_{123}","f_{014}","f_{024}",
    "f_{124}","f_{034}","f_{134}","f_{234}","f_{01234}"}
    


toLatex = l ->(
    presum := apply(0..#l-1,i-> if l#i==-1 then "-"|qs#i 
	else if l#i==1 then "+"|qs#i
	);
    -- sumqs := for i from 0 to #l-1 list (
    -- 	if l#i==0 then continue;
    -- 	if l#i==-1 then "-"|qs#i;
    -- 	if l#i==1 then "+"|qs#i;
    -- 	);
    -- return sum(presum)   				    
    return concatenate(presum) 
    )




sp5_9
inw(-{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},sp5)

wt = -{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
Rwt = newRing(ring sp5, Weights=>wt, Global => false)
leadTerm(1, sub(sp5_9,Rwt))

    Rw := newRing(R, Weights=>w, Global => false);
    inwI := sub(ideal(leadTerm(1, sub(I,Rw))),R);
    return stepwiseSaturate(inwI, gens R)
    )




basesToIndex = (Bs,n) -> (
    En := {};
    if n%2==0 then En = evenSubsets(n) else En = oddSubsets(n);
    return apply(Bs, B -> position(En, j->j==B)) 
    )


invertPerm = p -> apply(toList(0..#p-1), j-> position(p,i->i==j))

pOnBases = (p,Bsi,n) -> (
    pinv = invertPerm(p);
    return sort apply(Bsi, j->pinv#j)
    )
    
    
pOnMSDi = (p, msdi, n) ->  sort apply(msdi, Bsi-> pOnBases(p,Bsi,n))

pOnCone = (p,raysIndex,n) ->(
    raysTS5 := TS5#"rays"; 
    rays := apply(raysIndex, i-> raysTS5#i);
    prays := apply(rays, r -> for i in (0..#r-1) list r#(p#i));
    return sort apply(prays, pr-> position(raysTS5, j->pr==j ))
    )


symE5 = apply(lines get "symE5v2.dat", l -> apply(separate(" ", l), i->value(i)));


coneToMSDi = new MutableHashTable
for i from 1 to #(TS5#"conesOrbits")-1 do (
    ci := TS5#"conesOrbits"#i; 
    msd := matroidSubdivision(uniformS5,5,TS5#"relativeInteriorVectors"#(i-1)); 
    msdi := sort apply(msd#0, Bs -> basesToIndex(Bs,5)); 
    for j from 0 to #symE5-1 do (
	coneToMSDi#(pOnCone(symE5#j,ci,5)) = pOnMSDi(symE5#j, msdi,5); 
	)
    )

--msd5 = apply(TS5#"relativeInteriorVectors", w -> matroidSubdivision(uniformS5,5,w));

coneToMSDDirect = c -> (
    r := sum apply(c, i-> TS5#"rays"#i);
    msd := matroidSubdivision(uniformS5,5,r);
    return msd#1
    --return sort apply(msd#0, Bs -> basesToIndex(Bs,5))
    )



A = coneToMSDDirect({16,31,32,33,35})
B = coneToMSDDirect({19,31,32,33,35})
A#4 == B#4


MSDs = values coneToMSDi;
conesS5 = keys coneToMSDi;
msdiToCone = new MutableHashTable
Dr2Gr = new MutableHashTable
for i in (0..#conesS5-1) do (
    ci := conesS5#i; 
    msdi := coneToMSDi#ci; 
    if msdiToCone#?(msdi) then msdiToCone#(msdi) = sort toList(set msdiToCone#(msdi) + set ci) else msdiToCone#(msdi) = ci;
    if Dr2Gr#?(msdi) then Dr2Gr#msdi = append(Dr2Gr#msdi, ci) else Dr2Gr#msdi = {ci}; 
    )


msdis = keys msdiToCone;
Dr2GrV2 = hashTable for i in (0..#msdis-1) list {msdiToCone#(msdis#i),Dr2Gr#(msdis#i)};

DrCones = keys Dr2GrV2;
DrConesOrbits = new MutableHashTable;
seen = set();
for i in (0..#DrCones-1) do (
    c := DrCones#i; 
    if not member(c,seen) then (
	pcs := set apply(symE5,p->pOnCone(p,c,5) ); 
        DrConesOrbits#c = toList(pcs);
	seen = seen + pcs;
	)
    else continue; 
    )

keys DrConesOrbits

conesS5


s = {11,12,14,15,32,35}
for i in (0..#conesS5-1) list if isSubset(conesS5#i, s) and #conesS5#i>4 then conesS5#i else continue
















    

keys DrConesOrbits
Dr2GrV2






unique apply(evenSubsets(5), e->S5OnRay({4,1,0,3,0},  twistRay(e, ring sp5, {1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, 1, 1}))  )

S5OnRay({1,4,0,3,2},  twistRay({1,3}, ring sp5, {1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, 1, 1}))

 



niceRays = hashTable{
   0=> {0,1,1,1,1,1,1,0,0,0,0,1,0,0,0,0},
   1=> {0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,0},
   2=> {0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,1},
   3=> {0,0,0,0,0,0,0,1,1,1,1,0,1,0,1,1},
   4=> {0,0,0,0,0,0,0,1,1,1,1,0,0,1,1,1},
   5=> {1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0},
   6=> {0,0,0,1,0,1,0,0,1,0,0,1,1,1,0,1},
   7=> {0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,1},
   8=> {0,0,0,0,0,0,0,1,1,0,1,0,1,1,1,1},
   9=> {0,0,0,0,1,1,0,0,0,0,1,1,0,1,1,1},
   10=> {0,0,0,1,0,0,1,0,0,1,0,1,1,0,1,1},
   11=> {0,0,0,0,0,0,0,1,0,1,1,0,1,1,1,1},
   12=> {0,0,0,0,1,0,1,0,0,0,1,1,0,1,1,1},
   13=> {0,0,1,0,0,0,1,0,0,1,0,1,1,0,1,1},
   14=> {0,0,0,0,0,1,1,0,0,0,1,1,0,1,1,1},
   15=> {0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1},
   32=> {0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0},
   33=> {0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1},
   34=> {0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,1},
   35=> {0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1}
    }


inw(niceRays#11+niceRays#12+niceRays#14+niceRays#15+niceRays#32,sp5) == 
inw(niceRays#11+niceRays#12+niceRays#14+niceRays#15+niceRays#35,sp5)



inw(TS5#"relativeInteriorVectors"#0,sp5) == inw(niceRays#15,sp5)
inw(TS5#"relativeInteriorVectors"#1,sp5) == inw(niceRays#35,sp5)
inw(TS5#"relativeInteriorVectors"#2,sp5) == inw(niceRays#13+niceRays#15,sp5)
inw(TS5#"relativeInteriorVectors"#3,sp5) == inw(niceRays#14+niceRays#15,sp5)
inw(TS5#"relativeInteriorVectors"#4,sp5) == inw(niceRays#15+niceRays#35,sp5)

inw(TS5#"relativeInteriorVectors"#5,sp5) == inw(niceRays#12+niceRays#14+niceRays#15,sp5)
inw(TS5#"relativeInteriorVectors"#6,sp5) == inw(niceRays#11+niceRays#15+niceRays#35,sp5)
inw(TS5#"relativeInteriorVectors"#7,sp5) == inw(niceRays#12+niceRays#15+niceRays#35,sp5)
inw(TS5#"relativeInteriorVectors"#8,sp5) == inw(niceRays#13+niceRays#14+niceRays#15,sp5)
inw(TS5#"relativeInteriorVectors"#9,sp5) == inw(niceRays#14+niceRays#15+niceRays#35,sp5)

inw(TS5#"relativeInteriorVectors"#10,sp5) == inw(niceRays#5+niceRays#9+niceRays#11+niceRays#15,sp5)
inw(TS5#"relativeInteriorVectors"#11,sp5) == inw(niceRays#11+niceRays#12+niceRays#14+niceRays#15,sp5)
inw(TS5#"relativeInteriorVectors"#12,sp5) == inw(niceRays#13+niceRays#14+niceRays#15+niceRays#32,sp5)
inw(TS5#"relativeInteriorVectors"#13,sp5) == inw(niceRays#12+niceRays#13+niceRays#14+niceRays#15,sp5)
inw(TS5#"relativeInteriorVectors"#14,sp5) == inw(niceRays#12+niceRays#14+niceRays#15+niceRays#35,sp5)

inw(TS5#"relativeInteriorVectors"#15,sp5) == inw(niceRays#5+niceRays#9+niceRays#11+niceRays#15+niceRays#35,sp5)
inw(TS5#"relativeInteriorVectors"#16,sp5) == inw(niceRays#9+niceRays#11+niceRays#13+niceRays#14+niceRays#15,sp5)
inw(TS5#"relativeInteriorVectors"#17,sp5) == inw(niceRays#12+niceRays#13+niceRays#14+niceRays#15+niceRays#32,sp5)
inw(TS5#"relativeInteriorVectors"#18,sp5) == inw(niceRays#11+niceRays#12+niceRays#14+niceRays#15+niceRays#35,sp5)









A=(TS5#"lineality")|{TS5#"rays"#35}
matrix A

ringSpin(5,QQ)
r0={0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1}
r1={0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1}

TS5#"conesOrbits"
A = TS5#"rays"
#A
permuted = (g,l)->apply(g,i->l#i)

apply(sym5, g-> permuted(g,r0))

sym5#0
r0

















sp6 = spin(6,ringSpin(6,QQ));
rays6 = lines get "raysn6.txt";
indexToRays = rs -> apply(for j in rs list rays6#j, l -> for i in (0..31) list value((separate(" ",l))#i));
cones6 = lines get "conesn6.txt";
uniformS6 = {{}} | evenSubsets(6)

indexToCone = i -> apply(separate(" ", cones6#i), j -> value(j))

commonBasis = CONE -> (
    w:= sum indexToRays(CONE);
    msd := (matroidSubdivision(uniformS6,6,w))#0;
    return intersectionSets(apply(msd, s-> set s))
    )

wt = sum indexToRays(indexToCone(9))
matroidSubdivision(uniformS6,6,wt)

cBs = for i in (1..50) list commonBasis(indexToCone(i));

cBs

w= sum indexToRays({439, 4887, 7096, 27257, 65295, 117937, 147050, 147057, 184497, 219050, 219057, 269937, 278577, 316017, 327537, 339057})

toString w
w2 = sum indexToRays({179})

uniformS6 = {{}} | evenSubsets(6)
msd6 =  matroidSubdivision(uniformS6,6,w)
msd62 =  matroidSubdivision(uniformS6,6,w2)



intersectionSets(for i in (1,2,3,4,6,7,8,9,10,12,13,14,15,17,18,21,22,23,25,26,27,28,29,30) list set msd6#0#i)




for i in (0..18) do print msd5#i
msd5#1#1

R=QQ[t, Weights => {-1}, Global=> false]
S = R[flatten apply(0..4,i->apply(i+1..5, j-> value("x"|i|j)))]


tM = sub(genericSkewMatrix(S,x01,6), {
	x01 => t,
	x02 => t^2,
	x03 => t^3,
	x04 => t^4,
	x05 => t^9,
	x12 => t^7,
	x13 => t^5,
	x14 => t^2,
	x15 => t,
	x23 => t^11,
	x24 => t^2,
	x25 => t^3,
	x34 => t,
	x35 => t^7,
	x45 => t^13})


trop = M -> (
    n := numgens source M;
    eS := evenSubsets(n);
    return {0}|apply(eS, s -> (degree(leadTerm det submatrix(M,s,s)))#0)/2
    ) 




sym6={sOnCoords((1,2,3,4,5,0),6), sOnCoords((1,0,2,3,4,5),6)} | 
apply(evenSubsets(6), J-> twistToPerm(J,ring sp6))


sgn6={sOnCoordsSigns((1,2,3,4,5,0),6), sOnCoordsSigns((1,0,2,3,4,5),6)} | 
apply(evenSubsets(6), J-> twistToSign(J,ring sp6))



tropicalizeWithVector(sp6, ww, sym6, sgns6)

TS6 = parseTropFile("n6trop.txt")





-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
-- the following functions use affine coordinates to compute ideals of thin schubert cells, limits of thin Schubert cells 
-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------



affineTSC = (S, B, k) -> (
    oB := sort toList B;
    D := new MutableHashTable;
    d := rank(S);
    n := #(S.groundSet);
    xij := for j in (0..n-4) list (for i in (0..d-1) list x_(i,j) );
    R := k[flatten xij]; 
    I := entries id_(R^d); 
    xij = for j in (0..n-4) list (for i in (0..d-1) list sub(x_(i,j), R) ); 
    scan(toList (0..d-1), c-> xij = insert(oB#c, I#c, xij) ); 
    X := transpose matrix xij; 
    D#"ideal" = ideal(0_R);
    basesXij := bases(S) / (b -> det submatrix(X,, toList b));
    nonBasesXij := {}; 
        
    if #(nonbases(S)) != 0 then (
	variablesNonbases := flatten apply (for i in (0..#(nonbases(S))-1) list ( if #(((nonbases(S))#i) * B) == 2 then (nonbases(S))#i else continue), nb -> support det submatrix(X,, toList nb));
	basesXij = unique(apply(basesXij, f  ->  sub(f, variablesNonbases / (g -> g=>0))));
	nonBasesXij = (nonbases(S))/ (nb -> det submatrix(X,, toList nb));
	nonBasesXij = unique(apply(nonBasesXij, f  ->  sub(f, variablesNonbases/ (g -> g=>0 ))) );
	IMpreS := ideal nonBasesXij;
    	D#"ideal" = stepwiseSaturateNoBayer(IMpreS, basesXij);  
        );
    
    D#"ring" = R;
    D#"matrix" = X;
    D#"basesX" = basesXij;
    D#"nonBasesX" = nonBasesXij;
    return D
    )
-- Input: a matroid S, basis B of S, and field k
-- Output: a MutableHashTable D. Say S is a rank d matroid on n elements.  
--    D#"matrix" is a d-by-n matrix whose columns from B form the d-by-d identity matrix, and remaining entries filled in with x_(i,j).
--    D#"ring" = k[x_(i,j)]
--    D#"nonBasesX" = minors of D#"matrix" corresponding to nonbases of S, after setting x_(i,j) corresponding to nonbases = 0.
--    D#"basesX" = minors of D#"matrix" corresponding to bases of S, after setting x_(i,j) corresponding to nonbases = 0.
--    D#"ideal" = ideal generated by D#"nonBasesX" and saturated wrt D#"basesX"
 



 
affineLimitTSC = (S, SD, B, k) -> (
    D := new MutableHashTable;
    tscS := affineTSC(S,B,k);
    tscSD := SD/(Si ->  affineTSC(Si, B, k));
    ideals  := for i in (0..#SD-1) list sub(tscSD#i#"ideal", tscS#"ring");
    limitIdealPreSat := sum(for i in (0..#SD-1) list sub(tscSD#i#"ideal", tscS#"ring"));
    basesSD := unique apply( flatten for i in (0..#SD-1) list tscSD#i#"basesX", b -> sub(b,tscS#"ring"));
    
    D#"matrix" = tscS#"matrix";
    D#"ring" = tscS#"ring";
    D#"idealPreSat" = limitIdealPreSat;
    D#"ideal" = stepwiseSaturateNoBayer(limitIdealPreSat, basesSD);
    D#"basesX" = basesSD;
    D#"tsc" = tscSD;
    D#"graph" = graph(adjacencyMatrixSubDivision(SD));
    D#"ideals" = ideals;
    return D
    )
-- Input: a matroid S, a list of matroids SD (representing the maximal cells of a subdivision of S), basis B of S common to all matroids in SD, and field k
-- Output: a MutableHashTable D. Say S is a rank d matroid on n elements.  
--    D#"matrix" is a d-by-n matrix whose columns from B form the d-by-d identity matrix, and remaining entries filled in with x_(i,j).
--    D#"ring" = k[x_(i,j)]
--    D#"tsc" is a list of MutableHashTables affineTSC(Si, B ,k) for Si in SD. 
--    D#"basesX" = list of "basesX" from the elements of D#"tsc".
--    D#"idealPreSat" = ideal of D#"ring", this is the sum of the ideals from D#"tsc".
--    D#"ideal" = D#"idealPreSat" saturated with respect to D#"basesX".
--    D#"graph" = dual graph to the subdivision SD.
--    D#"ideals" = list of the ideals D#"tsc" as ideals of D#"ring".


