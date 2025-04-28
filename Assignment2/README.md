## Assignment 2

### Densest Subgraph Discovery

For Assignment 2, we have implemented the following densest subgraph discovery algorithms:

- [Exact Algorithm (Algorithm 1)] (https://arxiv.org/pdf/1906.00341)
- [Core Exact Algorithm (Algorithm 4)] (https://arxiv.org/pdf/1906.00341)

🔗 Link to website: https://advsinha17.github.io/DAAWebpage

### Dataset Preparation

We have implemented our algorithms on the following three datasets:

- [Yeast Network] (https://networkrepository.com/bio-yeast.php)
- [CAIDA AS Relationships] (https://snap.stanford.edu/data/as-Caida.html)
- [Autonomous Systems AS-733 (large)] (https://snap.stanford.edu/data/as-733.html)
- [Autonomous Systems AS-733 (small)] (https://snap.stanford.edu/data/as-733.html)
- [High Energy Physics - Theory collaboration network] (https://snap.stanford.edu/data/ca-HepTh.html)
- [Netscience] (https://networks.skewed.de/net/netscience)

Our implementations require the h value to be provided in the first line of the input text file, the second line must contain the number of vertices (n) and number of edges (m) separated by a space,
followed by m lines representing the edges of the undirected graph with the two vertices separated by a space.

The data can be found here: https://drive.google.com/drive/folders/1y2CIn9R2hgs0DlkdJKVn74Ost6avxIiH?usp=sharing

The files are labelled as (dataset name)\_h.txt, with h as the size of the h-Cliques used by the algorithm. This value is the same as the first line in the input text file and can be changed by editing the input file.

### Execution Instructions

Clone the repository, or download the `.cpp` file for the algorithm you wish to run

```
git clone https://github.com/advsinha17/DAA-Assignments.git
```

#### Exact Algorithm: `exact.cpp`

1. Move to the directory containing the file (Run the command below if you cloned the repository)

```
cd DAA-Assignments/Assignment2
```

2. Compile the program

```
g++ --std=c++17 -O3 exact.cpp -o exact.out
```

3. Run the executable, providing the **absolute** path of the dataset file as the command line argument

```
./exact.out <path-to-dataset-txt>
```

#### Core Exact Algorithm: `coreex.cpp`

1. Move to the directory containing the file (Run the command below if you cloned the repository)

```
cd DAA-Assignments/Assignment2
```

2. Compile the program

```
g++ --std=c++17 -O3 coreex.cpp -o coreex.out
```

3. Run the executable, providing the **absolute** path of the dataset file as the command line argument

```
./coreex.out <path-to-dataset-txt>
```

### Group Members

1. Srija Guda (2022A7PS2003H)
2. Khushi Kamtam (2022A7PS2010H)
3. Pranathi Voora (2022A7PS2012H)
4. Advik Sinha (2022A7PS2004H)

### Contributions

- Exact Algorithm (Algorithm 1): Srija Guda (f20222003@hyderabad.bits-pilani.ac.in | 2022A7PS2003H) and Pranathi Voora (f20222012@hyderabad.bits-pilani.ac.in | 2022A7PS2012H)
- Core Exact Algorithm (Algorithm 4): Advik Sinha (f20222004@hyderabad.bits-pilani.ac.in | 2022A7PS2004H) and Khushi Kamtam (f20222010@hyderabad.bits-pilani.ac.in | 2022A7PS2010H)
- Webpage: Advik Sinha (f20222004@hyderabad.bits-pilani.ac.in | 2022A7PS2004H) and Pranathi Voora (f20222012@hyderabad.bits-pilani.ac.in | 2022A7PS2012H)
- Report: Srija Guda (f20222003@hyderabad.bits-pilani.ac.in | 2022A7PS2003H) and Khushi Kamtam (f20222010@hyderabad.bits-pilani.ac.in | 2022A7PS2010H)
