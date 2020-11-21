#ifndef BB_TREE_H
#define BB_TREE_H

#include "KNode.h"
#include "BubbleBucket.h"
#include <vector>

class BB_Tree{
    vector<int> kary;
    int k; //size of inner search nodes
    int m; //dimensionality of data
public:
    //constructor

    //insert 
    /*
    1. nagivate through inner search nodes by dimension designated by level
    2. if super bubble bucket -> reorganize
       else If bubble bucket full -> create super bubble bucket
       else stick it in a bubble bucket 

    */

    //delete

    //search & range search?
};

#endif