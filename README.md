# HDPF
High dimensional peaks finding via topological data analysis.

# Install
```
  git clone https://github.com/fightinggggg/HDPF.git
  cd HDPF
  g++ PersistentHomology.cpp -std=c++11 -o ph
```
# Input 
The input file should be a plain text file in tsv format with n columns. With first n-1 columns denotes the coordinates while the nth column denotes the cooresponding probability landscape. One thing to remember that the coordiantes must be integers while the probability can be float numbers.
```
  cat test.txt
```

# Run
```
  ./ph test.txt
```

# Output
The output pf **HDPF** is a plain text in tsv format with 3 columns. 
Each column denotes the birth, death and lenght of persistent barcode.


# Acknowledgement
The **HDPF** program relies on **Eigen** library to perform matrix calculation. If you find **HDPF** helpful in your research, please do not forget to cite **Eigen** at the same time.

**Eigen**: http://eigen.tuxfamily.org/
