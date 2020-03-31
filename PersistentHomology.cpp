#include <iostream>
#include <fstream>
#include <sstream> 
#include <vector>
#include <map>
#include "Eigen/Dense"
using namespace Eigen;
using namespace std;

MatrixXd load_csv (const string& path);
MatrixXd rowSelection(MatrixXd mat, VectorXi select_rows);
MatrixXd colSelection(MatrixXd mat, VectorXi select_cols);
VectorXi IndexSelection(VectorXi select_index);
VectorXd sort_decending(VectorXd vec);
MatrixXd Neighbors(MatrixXd Points);
MatrixXd addRow(MatrixXd mat, RowVectorXd vec);
VectorXi addElement(VectorXi vec, int element);
MatrixXd rmRow(MatrixXd mat, int row);
MatrixXd rmDupRow(MatrixXd mat);
MatrixXd rmDupRows(MatrixXd mat);
bool MatrixContainRowVector(MatrixXd mat, RowVectorXd vec);
bool VectorContainElement(VectorXi vec, int element);
VectorXi extract_keys(map<int,VectorXi> mymap);
VectorXi extract_keys_double(map<int,double> mymap);
template<typename T> void vecDupRM(vector<T> vec);
VectorXi Intersect(VectorXi x, VectorXi y);
VectorXi Union(VectorXi x, VectorXi y);
MatrixXd UnionMatrix(MatrixXd A, MatrixXd B);
bool ManhattanDist(MatrixXd PointSet, MatrixXd n_Pos);
void PersistentHomologyAnalysis(MatrixXd mat);


// Main Function
int main(int argc, char *argv[])
{
	MatrixXd A = load_csv(argv[1]);
    //MatrixXd B = MatrixXd::Random(4,7);
    //TestMatrix(A,B);
    PersistentHomologyAnalysis(A);
    return 0;
}


// Functions
MatrixXd load_csv (const string& path) 
{
    ifstream indata;
    indata.open(path);
    string line;
    vector<double> values;
    uint rows = 0;
    while (getline(indata, line)) 
    {
        stringstream lineStream(line);
        string cell;
        while (getline(lineStream, cell, '\t')) 
        {
            values.push_back(stod(cell));
        }
        ++rows;
    }

    MatrixXd Mat = Map<Matrix<double,Dynamic,Dynamic,RowMajor> >(values.data(), rows, values.size()/rows);
    indata.close();
    return Mat;
}



MatrixXd rowSelection(MatrixXd mat, VectorXi select_rows)
{
    MatrixXd mat_sel(select_rows.sum(), mat.cols());
    int rownew = 0;
    for (int i = 0; i < mat.rows(); ++i) 
    {
        if (select_rows[i]) 
        {       
            mat_sel.row(rownew) = mat.row(i);
            rownew++;
        }
    }
    return mat_sel;
}



MatrixXd colSelection(MatrixXd mat, VectorXi select_cols)
{
    MatrixXd mat_sel(mat.rows(),select_cols.sum());
    int colnew = 0;
    for (int i = 0; i < mat.cols(); ++i) 
    {
        if (select_cols[i]) 
        {       
            mat_sel.col(colnew) = mat.col(i);
            colnew++;
        }
    }
    return mat_sel;
}



VectorXi IndexSelection(VectorXi select_index)
{
    VectorXi idx_sel(select_index.sum());
    int item = 0;
    for (int i = 0; i < select_index.size(); ++i) 
    {
        if (select_index[i]) 
        {       
            idx_sel(item) = i;
            item++;
        }
    }
    return idx_sel;
}



VectorXd sort_decending(VectorXd vec)
{
    ArrayXd sorted = vec.array();
    std::sort(sorted.begin(), sorted.end(), std::greater<double>());
    VectorXd result = sorted;
    return result;
}



MatrixXd Neighbors(MatrixXd Points)
{
    MatrixXd NeighborPoints;
    RowVectorXd Point;
    RowVectorXd Neighbor;
    for(int i = 0; i < Points.rows(); ++i)
    {   
        Point = Points.row(i);
        for(int j = 0; j < Point.size(); ++j)
        {
            Neighbor = Point;
            Neighbor[j] = Point[j] + 1;
            NeighborPoints.conservativeResize(NeighborPoints.rows() + 1, Points.cols());
            NeighborPoints.row(NeighborPoints.rows()-1) = Neighbor;


            Neighbor = Point;
            Neighbor[j] = Neighbor[j] - 1;
            NeighborPoints.conservativeResize(NeighborPoints.rows() + 1, Points.cols());
            NeighborPoints.row(NeighborPoints.rows()-1) = Neighbor;
        }
    }
    MatrixXd Result = rmDupRows(NeighborPoints);
    return Result;
}



MatrixXd addRow(MatrixXd mat, RowVectorXd vec) 
{
    mat.conservativeResize(mat.rows() + 1, mat.cols());
    mat.row(mat.rows()-1) = vec;
    return mat;
}



VectorXi addElement(VectorXi vec, int element) 
{
    vector<int> tmp(&vec[0], vec.data() + vec.size());
    tmp.push_back(element);
    Map<VectorXi> result(tmp.data(), tmp.size());

    return result;
}



MatrixXd rmRow(MatrixXd mat, int row) 
{
    int numRows = mat.rows() - 1;
    int numCols = mat.cols();
    
    if( row < numRows ) 
    {
        mat.block(row,0,numRows-row,numCols) = mat.block(row+1,0,numRows-row,numCols);
    }
    mat.conservativeResize(numRows,numCols);
    return mat;
}



MatrixXd rmDupRow(MatrixXd mat)
{
    for(int i = 0; i < mat.rows(); ++i)
    {
        for(int j = i + 1; j < mat.rows(); ++j)
        {
            if(mat.row(i) == mat.row(j))
            {
                mat = rmRow(mat,j);
                break;
            }

        }
    }
    return mat;
}



MatrixXd rmDupRows(MatrixXd mat)
{
    MatrixXd Update = mat;
    while(true)
    {
        int Before = Update.rows();
        Update = rmDupRow(Update);
        int After = Update.rows();
        if(Before == After)
        {
            break;
        }
    }
    return Update;
}



bool MatrixContainRowVector(MatrixXd mat, RowVectorXd vec)
{
    int contain = false;
    for(int i = 0; i < mat.rows(); ++i)
    {
        if(mat.row(i) == vec)
        {
            contain = true;
        }
    }
    return contain;
}

bool VectorContainElement(VectorXi vec, int element)
{
    int contain = false;
    for(int i = 0; i < vec.size(); ++i)
    {
        if(vec[i] == element)
        {
            contain = true;
        }
    }
    return contain;
}


VectorXi extract_keys(map<int,VectorXi> mymap)
{
    vector<int> v;
    for (map<int,VectorXi>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
    {
        v.push_back(it->first);
    }
    Map<VectorXi> result(v.data(), v.size());
    return result;
}

VectorXi extract_keys_double(map<int,double> mymap)
{
    vector<int> v;
    for (map<int,double>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
    {
        v.push_back(it->first);
    }
    Map<VectorXi> result(v.data(), v.size());
    return result;
}


vector<int> vecDupRM(vector<int> vec)
{
    std::sort(vec.begin(),vec.end());
    vec.erase(std::unique(vec.begin(),vec.end()),vec.end());
    return(vec);
}



VectorXi Intersect(VectorXi x, VectorXi y)
{
    vector<int> result;
    for(int i = 0; i<x.size(); ++i)
    {
        for(int j = 0; j<y.size(); ++j)
        {
            if(x[i] == y[j])
            {
                result.push_back(x[i]);
            }
        }
    }
    result = vecDupRM(result);
    Map<VectorXi> Intersection(result.data(), result.size());
    return Intersection;
}



VectorXi Union(VectorXi x, VectorXi y)
{
    vector<int> result;
    for(int i = 0; i<x.size(); ++i)
    {
        for(int j = 0; j<y.size(); ++j)
        {
            result.push_back(x[i]);
            result.push_back(y[j]);
        }
    }
    result = vecDupRM(result);
    Map<VectorXi> UnionSet(result.data(), result.size());
    return UnionSet;
}



MatrixXd UnionMatrix(MatrixXd A, MatrixXd B)
{
    MatrixXd C(A.rows()+B.rows(), A.cols()); 
    C << A, B;
    MatrixXd D = rmDupRows(C);
    return D;
}



bool ManhattanDist(MatrixXd PointSet, MatrixXd n_Pos)
{
    for(int i=0; i<=PointSet.rows(); ++i)
    {
        PointSet.row(i) = PointSet.row(i) - n_Pos;
    }
    MatrixXd ABS = PointSet.cwiseAbs();
    double Minimum = ABS.rowwise().sum().minCoeff();

    if(Minimum <= 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}


void PersistentHomologyAnalysis(MatrixXd mat)
{
    // Read density landscape file
    MatrixXd Dimension = mat.block(0,0,mat.rows(),mat.cols()-1);
    VectorXd Density = mat.block(0,mat.cols()-1,mat.rows(),1);
    int dim = Dimension.cols();
    cout << "Data dimension is = " << dim << endl;

    double Maximum = Density.maxCoeff();
    double Minimum = Density.minCoeff();
    cout << "Maximum Density = " << Maximum << endl;
    cout << "Minimum Density = " << Minimum << endl;

    map<int,VectorXi> PointIndexCluster;
    map<int,double> BarcodeBirth;
    map<int,double> BarcodeDeath;
    map<int,VectorXi> Peaks;
    VectorXi IndexNotes;
    VectorXi BeforeUpdate;
    VectorXi AfterUpdate;
    MatrixXd Point;
    VectorXi Keys;
    VectorXi I;
    VectorXi J;
    VectorXi Intersection;
    VectorXi BarcodeBirthKey;
    VectorXi BarcodeDeathKey;
    VectorXd NewRow;

    VectorXi select_rows_1st = (mat.col(mat.cols()-1).array() == Maximum).cast<int>();
    MatrixXd HighestLevelSet = rowSelection(Dimension, select_rows_1st);
    VectorXi Index = IndexSelection(select_rows_1st);

    cout << "HighestLevelSet: " << HighestLevelSet << endl;
    cout << "highest Index: " << Index << endl << endl;

//  Find the highest point in the density landscape

    for(int i = 0; i <  Index.size(); ++i)
    {
        Point = HighestLevelSet(i,all);
        PointIndexCluster[i] = Index;
        Peaks[i] = Index;
        BarcodeBirth[i] = Maximum;
    }
//  Superlevel set filtration 
    VectorXd Den = sort_decending(Density);

    for(int i = 1; i < Den.size(); ++i)

    {
        VectorXi select_rows_2nd = (mat.col(mat.cols()-1).array() == Den[i]).cast<int>();
        MatrixXd LevelSet = rowSelection(Dimension, select_rows_2nd);
        VectorXi Index = IndexSelection(select_rows_2nd);
        
        cout << "Density: " << Den[i] << endl;

        for(int n = 0; n < LevelSet.rows(); ++n)
        {
            MatrixXd n_pos = LevelSet.row(n);
            VectorXi n_index(1);
            n_index << Index(n);
            IndexNotes = extract_keys(PointIndexCluster);
            int knowncluster = IndexNotes(IndexNotes.size()-1);

   
            for(int key = 0; key < IndexNotes.size(); ++key)
            {
                VectorXi idx = PointIndexCluster[key];
                MatrixXd PointPosition = Dimension(idx,all);
                if( MatrixContainRowVector(Neighbors(PointPosition),n_pos) )
                {
                    PointIndexCluster[key] = addElement(PointIndexCluster[key],Index[n]);
                }
                
                else
                {
                    PointIndexCluster[knowncluster + 1] = n_index;
                    Peaks[knowncluster + 1] = n_index;
                    BarcodeBirth[knowncluster + 1] = Den[i];
                }         
            }
        }


        BeforeUpdate = extract_keys(PointIndexCluster);
        Keys = BeforeUpdate;
        if(Keys.size() > 1)
        {

            for(int x = 0; x < Keys.size(); ++x)
            {
                for(int y = x + 1; y < Keys.size(); ++y)
                {
<<<<<<< HEAD
                    I = PointIndexCluster[x];
                    J = PointIndexCluster[y];
                    IP = PointPositionCluster[x];
                    JP = PointPositionCluster[y];
=======
                    I = PointIndexCluster[i];
                    J = PointIndexCluster[j];
>>>>>>> c3e79e816509e5f2e6d3ff85ecd6e81900d3a7ad
                    Intersection = Intersect(I,J);

                    if(Intersection.size() > 0)
                    {

<<<<<<< HEAD
                        BarcodeDeath[y] = Den[i];
                        PointIndexCluster[x] = Union(I,J);
                        PointPositionCluster[x] = UnionMatrix(IP,JP);                        
                        PointIndexCluster.erase(y);
                        PointPositionCluster.erase(y);
=======
                        BarcodeDeath[j] = Den[i];
                        PointIndexCluster[i] = Union(I,J);                    
                        PointIndexCluster.erase(j);
>>>>>>> c3e79e816509e5f2e6d3ff85ecd6e81900d3a7ad

                    }
                }
            }
        }
    }

    BarcodeBirthKey = extract_keys_double(BarcodeBirth);
    BarcodeDeathKey = extract_keys_double(BarcodeDeath);
    cout << BarcodeBirthKey << "\nend" << endl;
    cout << BarcodeDeathKey << "\nend" << endl << endl;


    for(int i = 0; i<BarcodeBirthKey.size(); ++i)
    {
        int key = BarcodeBirthKey[i];
        cout << VectorContainElement(BarcodeDeathKey,key) << endl;
        if(not VectorContainElement(BarcodeDeathKey,key) )
        {
            BarcodeDeath[key] = Minimum;
            cout << BarcodeDeath[key] << endl;
        }
    }


    ofstream output;
    output.open ("PersistentDiagram.txt");
    output << "BarcodeBirth\tBarcodeDeath\tBarcodeLength" << endl;
    for(int i = 0; i<BarcodeBirthKey.size(); ++i)
    {
        int key = BarcodeBirthKey[i];
        output << BarcodeBirth[key] << "\t" << BarcodeDeath[key] << "\t" << BarcodeBirth[key] - BarcodeDeath[key]<<endl;
    }
    output.close();
    cout << "Persistent Diagram Calculating Finished." << endl;
}


