#pragma once

#include <type_traits>
#include <RcppArmadillo.h>

//template <class K, class V> 
//using map = std::unordered_map<K, V>;

//template <class T>
//using set = std::unordered_set<T>;

typedef int ent_id;
typedef int val_id;
typedef int rec_id;
typedef int file_id;
typedef int attr_id;

// Default empty values
// https://stackoverflow.com/questions/40070169/how-to-return-an-reference-to-an-empty-object
template <typename T>
struct emptyval {
  static const typename std::decay<T>::type value;
};

template <typename T>
const typename std::decay<T>::type emptyval<T>::value = T{};

typedef arma::ivec attr_vec;