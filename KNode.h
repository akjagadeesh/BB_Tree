#ifndef KNODE_H
#define KNODE_H
#include <vector>

class Node{
public:
    int k;
    int level;
    Node(int _k): k(_k), level(-1) {}
};

class KNode: public Node{
    vector<int> delims;
    vector<Node*> children;
public:
    KNode(int _k, int _level): Node(_k), level(_level){
        delims.resize(k-1);
        children.resize(k);
    }
};

#endif