## 1. The input format
The input graph is asked to be in DIMACS ascii format for undirected graphs.

A small example brock400_4.clq is shown in the supplement material:

## 2. statement

For classic instances, since their scale is relatively small, to accelerate the operations in Nukplex, both the adjacency matrix and adjacency list are adopted to represent a graph. The code is Nukplex_classic.cpp. 

As for the large-sparse instances, their scale is very large, so we can only adopt the adjacency list to represent a graph. The code is Nukplex_massive.cpp.

## 3. compile

For classic instances: g++ Nukplex_classic.cpp -O3 -o Nukplex
For classic instances: g++ Nukplex_massive.cpp -O3 -o Nukplex

## 4. Usage

The command to run Nukplex is:
```
./Nukplex -f <filename> -s <parameter k> -t <cutoff time (s)> -o<objective value> -c<random seed>
```
Nukplex terminates when reaching the cutoff time.

## 3. Output

 Instance: brock400_4.clq
 k: 2
 seed: 1
 best solution size: 31
 best solution time: 0
 Best solution: 288 117 313 105 251 315 383 189 56 384 179 76 10 14 0 343 52 377 299 309 165 162 67 362 353 186 253 88 1 119 95