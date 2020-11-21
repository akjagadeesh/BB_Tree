#ifndef BUBBLE_BUCKET_H
#define BUBBLE_BUCKET_H

#include "KNode.h"

struct Point{
    int k;
    vector<int> coords;
    Point(int _k, vector<int> _coords): k(_k), coords(_coords) {}
};

class BubbleBucket: public Node{
    int bmax;
    vector<Point> data;
public:
    BubbleBucket(int _k, int _level, int _bmax): Node(_k), level(_level), bmax(_bmax) {}
};

#endif