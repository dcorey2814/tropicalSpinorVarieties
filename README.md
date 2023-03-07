# Initial degenerations of Spinor varieties

Here, we record all of the code used in the paper <a href="https://arxiv.org/abs/2104.03442">Initial degenerations of spinor varieties</a> by <a href="https://www.danieljcorey.com/">Daniel Corey</a>. 

We compute the tropicalizations of <b>S</b><sub>4</sub>&deg; and <b>S</b><sub>5</sub>&deg; using `gfan` version 0.6.2. 


## Computing the tropicalization of <b>S</b><sub>4</sub>&deg; using  gfan.
First, compute the starting cone by running
```
gfan_tropicalstartingcone < startingConeInput4.txt > outputFile.txt
```
To record the symmetries, append the lines to the outputFile.txt: 
```
{{0,1,4,5,2,3,6,7}, {0,4,5,1,6,2,3,7}, {1,0,4,5,2,3,7,6}}
{{1,-1,1,1,1,1,1,-1}, {1,1,1,-1,1,-1,-1,-1}, {1,1,1,1,1,1,1,1}}
```
The result is stored in the file tropicalTraverseInput4.txt. To get the tropicalization of <b>S</b><sub>4</sub>&deg;, run
```
gfan_tropicaltraverse  --symmetry --symsigns --nocones < tropicalTraverseInput4.txt
```
The result is contained in the file TS4.txt.

## Computing the tropicalization of <b>S</b><sub>5</sub>&deg; using  gfan.
The computation is similar to the n=4 case above. First, compute the starting cone by running
```
gfan_tropicalstartingcone < startingConeInput5.txt > outputFile.txt
```
To record the symmetries, append the lines to the outputFile.txt: 
```
{{1, 2, 3, 4, 0, 8, 11, 13, 14, 5, 6, 7, 9, 10, 12, 15}, {1, 0, 2, 3, 4, 5, 6, 8, 7, 9, 11, 10, 13, 12, 14, 15}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, {1, 0, 5, 6, 9, 2, 3, 8, 7, 4, 11, 10, 13, 12, 15, 14}, {2, 5, 0, 7, 10, 1, 8, 3, 6, 11, 4, 9, 14, 15, 12, 13}, {5, 2, 1, 8, 11, 0, 7, 6, 3, 10, 9, 4, 15, 14, 13, 12}, {3, 6, 7, 0, 12, 8, 1, 2, 5, 13, 14, 15, 4, 9, 10, 11}, {6, 3, 8, 1, 13, 7, 0, 5, 2, 12, 15, 14, 9, 4, 11, 10}, {7, 8, 3, 2, 14, 6, 5, 0, 1, 15, 12, 13, 10, 11, 4, 9}, {4, 9, 10, 12, 0, 11, 13, 14, 15, 1, 2, 5, 3, 6, 7, 8}, {9, 4, 11, 13, 1, 10, 12, 15, 14, 0, 5, 2, 6, 3, 8, 7}, {10, 11, 4, 14, 2, 9, 15, 12, 13, 5, 0, 1, 7, 8, 3, 6}, {12, 13, 14, 4, 3, 15, 9, 10, 11, 6, 7, 8, 0, 1, 2, 5}, {8, 7, 6, 5, 15, 3, 2, 1, 0, 14, 13, 12, 11, 10, 9, 4}, {11, 10, 9, 15, 5, 4, 14, 13, 12, 2, 1, 0, 8, 7, 6, 3}, {13, 12, 15, 9, 6, 14, 4, 11, 10, 3, 8, 7, 1, 0, 5, 2}, {14, 15, 12, 10, 7, 13, 11, 4, 9, 8, 3, 6, 2, 5, 0, 1}, {15, 14, 13, 11, 8, 12, 10, 9, 4, 7, 6, 3, 5, 2, 1, 0}}
{{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}
```
The result is stored in the file tropicalTraverseInput5.txt. To get the tropicalization of <b>S</b><sub>5</sub>&deg;, run
```
gfan_tropicaltraverse  --symmetry --symsigns --nocones < tropicalTraverseInput5.txt
```
The result is contained in the file TS5.txt.

## From gfan computation to Trop <b>S</b><sub>5</sub>&deg; in the paper. 
The rays of Trop <b>S</b><sub>5</sub>&deg; computed by gfan are complicated expressions, so in the paper we use different primitive generators modulo the lineality space. To verify that these are the same, see the Oscar notebook `check_fan_oscar.ipynb`. 

## Proof of Lemma 7.10.
To prove Lemma 7.10, we apply Lemma 7.3 to each nonmaximal cone &tau; of &Sigma;<sub>5</sub>'. The collections of cones A<sub>&tau;</sub> are listed in Table 7.1. To verify Equation 7.1 for these cones, see the Oscar notebook `lem-7-10.ipynb` (this works with Oscar 0.11.3). 

## Matroid subdivisions of <b>S</b><sub>5</sub>&deg;.
This is done in polymake, and works with version 4.0. See the notebook `subdivisionsS5.ipynb`. 

## References to software

We use the following programs. 

<a href="https://users-math.au.dk/jensen/software/gfan/gfan.html">gfan</a>

```
@Misc{gfan,
	author = {Jensen, Anders N.},
	title = {Gfan, a software system for {G}r{\"o}bner fans and tropical varieties},
	howpublished = {{\tt http://home.imf.au.dk/jensen/software/gfan/gfan.html}}
} 
```

<a href="https://www.oscar-system.org/">Oscar</a>

```
@misc{OSCAR,
  key          = {OSCAR},
  organization = {The OSCAR Team},
  title        = {{OSCAR} -- Open Source Computer Algebra Research system,
                  Version  0.10.1},
  year         = {2022},
  url          = {https://oscar.computeralgebra.de},
  }


@Book{OSCAR-book,
  editor = {Eder, Christian and Decker, Wolfram and Fieker, Claus and Horn, Max and Joswig, Michael},
  title = {The {OSCAR} book},
  year = {2024},
}
```

<a href="https://www.polymake.org/doku.php">polymake</a>

```
@incollection{polymake:2000,
    AUTHOR = {Gawrilow, Ewgenij and Joswig, Michael},
     TITLE = {{\tt polymake}: a framework for analyzing convex polytopes},
 BOOKTITLE = {Polytopes---combinatorics and computation ({O}berwolfach, 1997)},
    SERIES = {DMV Sem.},
    VOLUME = {29},
     PAGES = {43--73},
 PUBLISHER = {Birkh\"auser, Basel},
      YEAR = {2000},
   MRCLASS = {52B55 (68U05)},
  MRNUMBER = {1785292},
}  

@InProceedings{polymakeJL,
    author={Kaluba, Marek and Lorenz, Benjamin and Timme, Sascha},
    editor={Bigatti, Anna Maria and Carette, Jacques and Davenport, James H. and Joswig, Michael and de Wolff, Timo},
    title={Polymake.jl: A New Interface to polymake},
    booktitle={Mathematical Software -- ICMS 2020},
    year={2020},
    publisher={Springer International Publishing},
    address={Cham},
    pages={377--385},
    abstract={We present the Julia interfaceto polymake, a software for research in polyhedral geometry. We describe the technical design and how the integration into Julia makes it possible to combine polymake with state-of-the-art numerical software.},
    isbn={978-3-030-52200-1}
}

```  

In previous versions of this project, we used sage.

<a href="https://www.sagemath.org/">sage</a>

```
@manual{sagemath,
  Key          = {SageMath},
  Author       = {{The Sage Developers}},
  Title        = {{S}ageMath, the {S}age {M}athematics {S}oftware {S}ystem ({V}ersion 9.0)},
  note         = {{\tt https://www.sagemath.org}},
  Year         = {2020},
}
```
