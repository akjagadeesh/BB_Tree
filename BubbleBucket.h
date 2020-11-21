#ifndef BUBBLE_BUCKET_H
#define BUBBLE_BUCKET_H

#include <vector>
#include <exception>
#include <string>
#include <sstream>
#include <algorithm>
using namespace std;

class Bucket{
public:
    int bmax; // maximum points stored in a bucket (bmax for BubbleBucket, k*bmax for SuperBubbleBucket)
    int k; // size of inner search nodes (used mostly in SuperBB)
    int m; // dimensionality of data points
    int id; // key for indexing from BBTree
    int level; // used for insert/search decisions

    Bucket(int _k, int _m, int _bmax, int _id): k(_k), m(_m), bmax(_bmax), id(_id), level(int(log(id)/log(k))%m) {}

    virtual int insert(Point p, int tid) = 0;
    virtual Point remove(int tid) = 0;
    virtual Point find(int tid) = 0;
};

struct Point{
    int m;
    vector<int> coords;

    Point(): m(-1){}

    Point(int _m, vector<int> _coords): m(_m) {
        if(_coords.size() != m){
            stringstream ss;
            ss << "Error: Point coords: {";
            for(auto i: _coords){
                ss << i << " ";
            }
            ss << "} not matching m: " << m;
            throw ss.str();
        }
        coords = _coords;
    }
};

class BubbleBucket: public Bucket{
    vector<Point> data;
    vector<int> tids;
    int fill;

public:
    BubbleBucket(int _k, int _m, int _bmax, int _id): Bucket(_k,_m,_bmax,_id), fill(0) {
        data.resize(bmax);
        tids.resize(bmax);
    }

    int insert(Point p, int tid){
        if(fill == bmax){ return -1; } //overflow -> create SuperBB
        else{
            data[fill++] = p;
            tids[fill++] = tid;
        }
        return 1; //-1 for overflow and 1 for success
    }

    Point remove(int tid){
        if(empty()) throw "Error: Remove of empty BubbleBucket";
        for(int i = 0; i < fill; i++){
            if(tid == tids[i]){
                Point victim = data[i];
                for(int j = i; j < fill-1; j++){
                    data[j] = data[j+1];
                    tids[j] = tids[j+1];
                }
                fill--;
                return victim;
            }
        }
        stringstream ss;
        ss << "Error: BubbleBucket(" << id << ") remove didn\'t find tid: " << tid;
        throw ss.str();
    }

    Point find(int tid){
        for(int i = 0; i < fill; i++){
            if(tid = tids[i]){
                return data[i];
            }
        }
        stringstream ss;
        ss << "Error: BubbleBucket(" << id << ") find didn\'t find tid: " << tid;
        throw ss.str();
    }

    bool empty(){ return fill == 0; }

    pair<int,Point> pop(){ //used for morphing to SuperBB
        if(empty()) throw "Error: Pop of empty BubbleBucket";
        return make_pair(tids[fill-1],data[--fill]);
    }
};

class SuperBubbleBucket: public Bucket{
    vector<int> knode;
    vector<BubbleBucket> buckets;

public:
    SuperBubbleBucket(BubbleBucket old, Point p, int tid): Bucket(old.k,old.m,old.bmax,old.id) {
        vector<pair<int,Point>> temp;
        temp.push_back(make_pair(tid,p));

        int kmin = p.coords[level], kmax = p.coords[level];
        while(!old.empty()){
            temp.push_back(old.pop());
            if(temp[temp.size()-1].second.coords[level] < kmin) kmin = temp[temp.size()-1].second.coords[level];
            if(temp[temp.size()-1].second.coords[level] > kmax) kmax = temp[temp.size()-1].second.coords[level];
        }
        for(int i = 0; i < k-1; i++){
            knode.push_back(kmin+(i+1)*(kmax-kmin)/(k-1));
        }

        for(int i = 0; i < k; i++){
            buckets.push_back(BubbleBucket(old.k,old.m,old.bmax,old.id));
        }
        for(auto p: temp){
            insert(p.second,p.first);
        }
    }

    int insert(Point p, int tid){
        int delim = p.coords[level];
        int check = 0;
        for(int i = 0; i < knode.size(); i++){
            if(delim <= knode[i]){
                check = buckets[i].insert(p,tid);
            }
        }
        if(check == 0) check = buckets[buckets.size()-1].insert(p,tid);
        return check; //-1 for overflow and 1 for success
    }

    Point remove(Point p, int tid){
        int delim = p.coords[level];
        for(int i = 0; i < knode.size(); i++){
            if(delim <= knode[i]){
                try{
                    return buckets[i].remove(tid);
                }
                catch(string s){}
            }
        }
        try{
            return buckets[buckets.size()-1].remove(tid);
        }
        catch(string s){
            stringstream ss;
            ss << "Error: SuperBubbleBucket(" << id << ") remove didn\'t find tid: " << tid;
            throw ss.str();
        }
        catch(...){
            stringstream ss;
            ss << "Unknown Error in SuperBubbleBucket(" << id << ")";
            throw ss.str();
        }
    }

    Point remove(int tid){
        for(auto bucket: buckets){
            try{
                return bucket.remove(tid);
            }
            catch(string s){}
        }
        stringstream ss;
        ss << "Error: SuperBubbleBucket(" << id << ") remove didn\'t find tid: " << tid;
        throw ss.str();
    }

    Point find(int tid){
        for(auto bucket: buckets){
            try{
                return bucket.find(tid);
            }
            catch(...){}
        }
        stringstream ss;
        ss << "Error: SuperBubbleBucket(" << id << ") find didn\'t find tid: " << tid;
        throw ss.str();
    }
};

#endif