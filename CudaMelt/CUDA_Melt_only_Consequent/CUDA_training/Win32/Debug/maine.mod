V24 maine
123 C:\Programming\CUDA fortran\CUDA_Melt_DIFFERENT_CONDITIONS\DOUBLE_PECISION\CUDA_Melt_only_Consequent\CUDA_training\Cons.f90 S582 0
06/09/2015  11:05:40
use user public 0 direct
use var public 0 direct
enduse
D 147 21 9 3 19 21 0 0 0 0 0
 0 13 3 3 13 13
 0 13 13 3 13 13
 0 13 20 3 13 13
S 582 24 0 0 0 6 1 0 4659 10005 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 294 0 0 0 0 0 0 maine
S 601 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 99 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 603 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 9801 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 604 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 970299 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 605 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 9901 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 699 23 5 0 0 0 700 582 5124 0 2000000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 setup1
S 700 14 5 0 0 0 1 699 5124 0 2400000 0 0 0 10 0 0 0 0 0 0 0 0 0 0 0 0 0 301 0 582 0 0 0 0 setup1
F 700 0
S 701 23 5 0 0 0 702 582 5131 0 2000000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 setup2
S 702 14 5 0 0 0 1 701 5131 0 2400000 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 374 0 582 0 0 0 0 setup2
F 702 0
S 703 23 5 0 0 0 709 582 5138 0 2000000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 solve_1
S 704 1 3 0 0 6 1 703 5146 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ntimes
S 705 1 3 0 0 6 1 703 4806 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ist
S 706 1 3 0 0 6 1 703 4810 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 jst
S 707 1 3 0 0 6 1 703 4814 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 kst
S 708 7 3 0 0 147 1 703 4972 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 f
S 709 14 5 0 0 0 1 703 5138 0 2400000 0 0 0 12 5 0 0 0 0 0 0 0 0 0 0 0 0 442 0 582 0 0 0 0 solve_1
F 709 5 705 706 707 704 708
S 710 23 5 0 0 0 711 582 5153 0 2000000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 solve
S 711 14 5 0 0 0 1 710 5153 0 2400000 0 0 0 18 0 0 0 0 0 0 0 0 0 0 0 0 0 576 0 582 0 0 0 0 solve
F 711 0
S 712 23 5 0 0 0 713 582 5159 0 2000000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 reset
S 713 14 5 0 0 0 1 712 5159 0 2400000 0 0 0 19 0 0 0 0 0 0 0 0 0 0 0 0 0 599 0 582 0 0 0 0 reset
F 713 0
A 13 2 0 0 0 6 601 0 0 0 13 0 0 0 0 0 0 0 0 0
A 19 2 0 0 0 6 605 0 0 0 19 0 0 0 0 0 0 0 0 0
A 20 2 0 0 0 6 603 0 0 0 20 0 0 0 0 0 0 0 0 0
A 21 2 0 0 0 6 604 0 0 0 21 0 0 0 0 0 0 0 0 0
Z
Z