#ifndef BB_TREE_H
#define BB_TREE_H

#include "BubbleBucket.h"
#include <vector>
#include <map>
#include <algorithm>
#include <random>

class BB_Tree{
    vector<int> kary;
    map<int,Bucket*> buckets;
    int size = 0; // actual size of the inner search tree
    int count = 0; // how many data points are stored in the tree
    int k; // size of inner search nodes
    int m; // dimensionality of data
    int bmax; // maxmimum allowed in a BubbleBucket
    double Rsample = 0.1;
    int tid = 0;

    int parent(int i){ return (k-1)*floor((i-k+1)/k); }

    int childi(int p, int i){ return k*p+(i+1)*(k-1); }

    bool is_internal(int i){
        bool flag = false;
        for(int child = 0; child < k; child++){
            flag ||= true ? childi(i,child) < size && kary[childi(i,child)] != -1: false;
        }
        return flag;
    }

    bool is_leaf(int i){
        return !is_internal(i);
    }
public:
    BB_Tree(int _k, int _m, int _bmax, double _Rsample): k(_k), m(_m), bmax(_bmax), Rsample(_Rsample) {
        buckets[0] = new BubbleBucket(k,m,bmax,0);
    }

    int insert(Point p){
        int level = 0;
        int node = 0;
        bool changed = false;
        tid++;
        while(!is_leaf(node)){
            for(int i = 0; i < k-1; i++){
                if(kary[node+i] <= p.coords[level]){
                    node = childi(node,i);
                    changed = true;
                }
            }
            if(!changed) node = childi(node,k-1);
            level++;
        }
        int flag = buckets[node]->insert(p,tid);
        if(flag == -1){
            if(buckets[node]->superId()){
                //reorganize inner search tree
            
                //1. Calculate new # of BB and level for IST (n/.1*bmax)
                int newBBs = count/(0.1*bmax);
                int height = ceil(log(newBBs)/log(k));
                kary.resize((pow(k,height-1)+1)*(k-1));

                //2. Randomly sample and sort samples
                int samples = Rsample*n;
                vector<Point> data;
                data.push_back(p); // don't forget the new point
                for(auto bucket: buckets){
                    vector<Point> temp = bucket->dump();
                    data.insert(data.end(), temp.begin(), temp.end());
                    delete *bucket;
                    *bucket = NULL;
                }
                std::random_device rd;
                std::mt19937 g(rd());
                std::shuffle(data.begin(), data.end(), g);
                vector<Point> sample;
                std::sample(data.begin(), data.end(), std::back_inserter(sample), samples, g);

                //3. Determine knode values based on the range of the sample data
                int kmin, kmax;
                for(int l = 0; l < height; l++){ //for each level
                    kmin = 32767;
                    kmax = -1;
                    for(int idx = 0; idx < sample.size(); idx++){ //determine min/max of dim every
                        if(sample[idx].coords[l % m] < kmin) kmin = sample[idx].coords[l % m];
                        if(sample[idx].coords[l % m] > kmax) kmax = sample[idx].coords[l % m];
                    }
                    for(int i = 0; i < (k-1)*pow(k,level); i++){
                        kary.push_back(kmin+(i+1)*(kmax-kmin)/(k-1));
                    }
                }

                //4. Insert into the new BBs
                buckets.erase(buckets.begin(), buckets.end());
                for(int i = kary.size(); i < kary.size()+newBBs; i++){
                    buckets[i] = new BubbleBucket(k,m,bmax,i);
                }
                for(auto point: data){
                    insert(point);
                }
            }
            else{
                Bucket* holder = buckets[node];
                Bucekt* newBB = new SuperBubbleBucket(*holder,p,tid);
                delete holder;
                holder = NULL;
                buckets[node] = newBB;
            }
        }
        count++;
        return tid;
    }

    //delete
    Point remove(Point p, int _tid){
        int level = 0;
        int node = 0;
        bool changed = false;
        while(!is_leaf(node)){
            for(int i = 0; i < k-1; i++){
                if(kary[node+i] <= p.coords[level]){
                    node = childi(node,i);
                    changed = true;
                }
            }
            if(!changed) node = childi(node,k-1);
            level++;
        }
        Point ret = buckets[node]->remove(_tid);
        count--;
        return ret;
    }

    //Point remove(Point p) {}

    //search & range search?
};

#endif