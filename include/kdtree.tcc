/**
 * @file kdtree.tcc
 * Implementation of kd_tree class.
 */

#include "kdtree.h"
#define square(x) ((x)*(x))
template <int Dim>
bool kd_tree<Dim>::smaller_in_dimension(const point<Dim>& first,
                                       const point<Dim>& second,
                                       int curDim) const
{
   
    if (second[curDim] < first[curDim])
    {
        return false;
    }
    else
    {
        return true;
    }
}

template <int Dim>
bool kd_tree<Dim>::should_replace(const point<Dim>& target,
                                 const point<Dim>& current_best,
                                 const point<Dim>& potential) const
{   
   
if(euclideanDistance(target, potential) < euclideanDistance(target, current_best)){
        return true;
    }
    else if(euclideanDistance(target, potential) == euclideanDistance(target, current_best)){
        if(potential < current_best){
            return  true;
        }
        else{
            return false;
        }
    }
    else{
        return false;
    }
}

template <int Dim> 
int kd_tree<Dim>::euclideanDistance(const point<Dim> &first, const point<Dim> &second)const{
    int distance = 0;
    for(int i = 0; i < Dim; i++){
        distance = square(first[i] - second[i]);
    }
    return distance;
}

template <int Dim>
void kd_tree<Dim>::build_kd_tree(std::vector< point<Dim> > & newpoints, int left, int right, int curDim){
    if(left > right){
        return;
    }
    int mid = (left + right) / 2;
    points[mid] = select(newpoints, left, right, mid, curDim);

    build_kd_tree(newpoints, left, mid - 1, (curDim + 1)%Dim);
    build_kd_tree(newpoints, mid + 1, right, (curDim + 1)%Dim);
}


template < int Dim>
point<Dim> kd_tree<Dim>::select(std::vector< point<Dim> > &list, int left, int right, int k, int curDim){
    while(left < right){

        int pivotIndex = partition(list, left, right, k, curDim);

        if(pivotIndex == k){
            return list[pivotIndex];
        }
        else if(k < pivotIndex){
            right = pivotIndex - 1;
        }
        else{
            left = pivotIndex + 1;
        }
    }
    return list[left];
}

template < int Dim>
int kd_tree<Dim>::partition(std::vector< point<Dim> > & list, int left, int right, int pivotIndex, int curDim){
    point< Dim> pivotValue = list[pivotIndex];
    std::swap(list[pivotIndex], list[right]);
    int storeIndex = left;
    int i = left;

    while(i < right){

        if(smaller_in_dimension(list[i], pivotValue, curDim)){

            std::swap(list[i], list[storeIndex]);
            storeIndex++;
        }

        i++;
    }
    std::swap(list[right], list[storeIndex]);

    return storeIndex;
}

template <int Dim>
kd_tree<Dim>::kd_tree(const std::vector<point<Dim>>& newpoints)
{
    points = newpoints;

    if(points.size() == 0){
        return;
    }
    build_kd_tree(points, 0, (int) (points.size() - 1), 0);
   
}


template <int Dim>
point<Dim> kd_tree<Dim>::find_nearest_neighbor(const point<Dim>& query) const
{
   
}
