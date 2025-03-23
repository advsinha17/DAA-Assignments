## Assignment 1

### Enumeration Of Maximal Cliques

For Assignment 1, we have implemented the following maximal clique finding algorithms:

- [Bron-Kerbosch Algorithm] (https://arxiv.org/abs/1006.5440)
- [CLIQUE Algorithm (Tomita)] (https://www.sciencedirect.com/science/article/pii/S0304397506003586)
- [Chiba Algorithm] (https://www.cs.cornell.edu/courses/cs6241/2019sp/readings/Chiba-1985-arboricity.pdf)

🔗 Link to website: https://advsinha17.github.io/DAAPage

### Dataset Preparation

We have implemented our algorithms on the following three datasets:

- [Email Enron Network] (https://snap.stanford.edu/data/email-Enron.html)
- [Autonomous Systems by Skitter] (https://snap.stanford.edu/data/as-Skitter.html)
- [Wikipedia Vote Network] (https://snap.stanford.edu/data/wiki-Vote.html)

Our implementations expect the format of input provided by the above 3 datasets (in particular, one or more lines beginning with a `#`, with one of them containin the number of nodes followed by the number of edges in the graph. Following these should be lines containing two integers, separated by a space). Hence, our algorithms do not require any specific dataset preparation for the provided datasets.

The data can be found here: https://drive.google.com/drive/folders/1K-wc67ojHsqgaeSR74S-1lfFU6qpxdoJ?usp=drive_link

### Execution Instructions

Clone the repository, or download the `.cpp` file for the algorithm you wish to run

```
git clone https://github.com/advsinha17/DAA-Assignments.git
```

#### Bron-Kerbosch Algorithm: `els.cpp`

1. Move to the directory containing the file (Run the command below if you cloned the repository)

```
cd DAA-Assignments/Assignment1
```

2. Compile the program

```
g++ --std=c++17 -O3 els.cpp -o els.out
```

3. Run the executable, providing the **absolute** path of the dataset file as the command line argument

```
./els.out <path-to-dataset-txt>
```

#### CLIQUE Algorithm (Tomita): `tomita.cpp`

1. Move to the directory containing the file (Run the command below if you cloned the repository)

```
cd DAA-Assignments/Assignment1
```

2. Compile the program

```
g++ --std=c++17 -O3 tomita.cpp -o tomita.out
```

3. Run the executable, providing the **absolute** path of the dataset file as the command line argument

```
./tomita.out <path-to-dataset-txt>
```

#### Chiba Algorithm: `chiba.cpp`

1. Move to the directory containing the file (Run the command below if you cloned the repository)

```
cd DAA-Assignments/Assignment1
```

2. Compile the program

```
g++ --std=c++17 -O3 chiba.cpp -o chiba.out
```

3. Run the executable, providing the **absolute** path of the dataset file as the command line argument (NOTE: Chiba must be run as sudo as it requires allocation of extra stack space)

```
sudo ./chiba.out <path-to-dataset-txt>
```

### Contributions

- Listing All Maximal Cliques in Sparse Graphs in Near-optimal Time: Advik Sinha (f20222004@hyderabad.bits-pilani.ac.in | 2022A7PS2004H)
- The worst-case time complexity for generating all maximal cliques and computational experiments: Khushi Kamtam (f20222010@hyderabad.bits-pilani.ac.in | 2022A7PS2010H) and Pranathi Voora (f20222012@hyderabad.bits-pilani.ac.in | 2022A7PS2012H)
- Arboricity and Subgraph Listing Algorithms: Srija Guda (f20222003@hyderabad.bits-pilani.ac.in | 2022A7PS2003H)
- Webpage: Advik Sinha (f20222004@hyderabad.bits-pilani.ac.in | 2022A7PS2004H) and Pranathi Voora (f20222012@hyderabad.bits-pilani.ac.in | 2022A7PS2012H)
- Report: Srija Guda (f20222003@hyderabad.bits-pilani.ac.in | 2022A7PS2003H) and Khushi Kamtam (f20222010@hyderabad.bits-pilani.ac.in | 2022A7PS2010H)
