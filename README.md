# A-HyperCut
This is implementation of Hypergraph Local Clusteringg algorithm **A-HyperCut** based on sweep cuts of **Averaging-based Personalized Page Rank for Hypergraphs (APPRH)**. The code can be run on Linux operating system also. But we recommend to use Visual Studio on Windows operating system. The following are required before you can download and run this project on your system.
## Requirements
**Visual Studio 2022**  
**.NET 6.0**  
**MathNet** library  
**Linq** Library  
# How to launch the project in Visual Studio
Assuming that your system satisfies the requirements listed above. And, you have extracted this project to the folder **A-HyperCut**. Go to the folder **A-HyperCut**, and double click on the file **HyperGraphClustering.sln**. This will automatically launch the project in Visual Studio. You need to build the solution by clicking on *Build Solution* in *Build* menu of Visual Studio before you run the project. Finally, click on *Star Without Debugging* in the *Debug* menu of Visual Studio. The output of the code is saved in *.txt file in A-HyperCut/HyperGraphClustering/results/
# Datasets
The datasets are saved as *.txt files in the folder A-HyperCut/HyperGraphClustering/instance/  
Each line of any dataset file represents a hyperedge. A line of any dataset file consists of the vertices contained in the corresponding hyperedge. Any two vertices in the same line are separated by a space. The last number of every line is weight of the corresponding hyperedge. The first 10 lines of the dataset file **graphprod_LCC.txt** are shown below.  
  
21 39 40 41 42 1  
21 8 1  
21 8 9 1  
15 2  
15 5 1  
15 29 30 1  
15 31 32 1  
15 10 3  
15 10 28 1  
15 33 1  

