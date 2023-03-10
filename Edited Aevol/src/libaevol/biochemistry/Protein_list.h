// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

#ifndef AEVOL_PROTEINLIST_H_
#define AEVOL_PROTEINLIST_H_

// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cstdlib>

#include <vector>
#include <memory>
#include <string>
#include <cinttypes>
#include <getopt.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <zlib.h>
#include <err.h>
#include <cerrno>
#include <sys/stat.h>
#include <unistd.h>
#include <list>
#include <iostream>
#include <fstream>

#include "aevol.h"

using namespace aevol;

// =================================================================
//                          Class declarations
// =================================================================
class protein_par{
public:
  protein_par(){}
  int32_t start_pos_;
  int32_t prot_id_;
  int32_t prot_parent_id_;
  int32_t length_;
  double basal_level_;
  int32_t hamming_dist_;
  int32_t dist_next_prot_;
  long double m, h, w, c;
  std::string functional;
  std::string strand;
  
  bool operator < (const protein_par& str) const
  {
    return (start_pos_ < str.start_pos_);
  }
  
  bool operator > (const protein_par& str) const
  {
    return (start_pos_ > str.start_pos_);
  }
  
};

class Protein_List {
public:
  Protein_List(int nb_prot, int dna_length) {
    nb_prot_ = nb_prot;
    /*start_pos_.resize(nb_prot_);
     length_.resize(nb_prot_);
     basal_level_.resize(nb_prot_);
     hamming_dist_.resize(nb_prot_);
     dist_next_prot_.resize(nb_prot_);
     prot_id_.resize(nb_prot_);*/
    pro_list_.resize(nb_prot_);
    dna_length_ = dna_length;
  }
  
  /*void hello(int32_t start_pos, int32_t length, double basal_level, int32_t hamming_dist, int32_t dist_next_prot, std::string prot_id){
   start_pos_[cpt_] = start_pos;
   length_[cpt_] = length;
   basal_level_[cpt_] = basal_level;
   hamming_dist_[cpt_] = hamming_dist;
   dist_next_prot_[cpt_] = dist_next_prot;
   prot_id_[cpt_] = prot_id;
   cpt_++;
   std::cout << "hello world";
  }*/
  void sort_protein_list(){
    std::sort(pro_list_.begin(), pro_list_.end());
  }
  
  void initialize_protein(int32_t start_pos, int32_t length, double basal_level, int32_t hamming_dist, int32_t dist_next_prot, int32_t prot_id, int32_t prot_parent_id) {
    pro_list_[cpt_].start_pos_ = start_pos;
    pro_list_[cpt_].prot_id_ = prot_id;
    pro_list_[cpt_].prot_parent_id_ = prot_parent_id;
    pro_list_[cpt_].prot_parent_id_ = 0;
    pro_list_[cpt_].length_ = length;
    pro_list_[cpt_].basal_level_ = basal_level;
    pro_list_[cpt_].hamming_dist_ = hamming_dist;
    pro_list_[cpt_].dist_next_prot_ = dist_next_prot;
    cpt_++;
  }
  
  void insert_protein(int32_t start_pos, int32_t length, double basal_level, int32_t hamming_dist, int32_t dist_next_prot, int32_t prot_id, int32_t prot_parent_id) {
    //resize vectors
    nb_prot_++;
    /*start_pos_.resize(nb_prot_);
     length_.resize(nb_prot_);
     basal_level_.resize(nb_prot_);
     hamming_dist_.resize(nb_prot_);
     dist_next_prot_.resize(nb_prot_);*/
    pro_list_.resize(nb_prot_);
    //insert values
    /*start_pos_[cpt_] = start_pos;
     prot_id_[cpt_] = (prot_id);
     length_[cpt_] = length;
     basal_level_[cpt_] = basal_level;
     hamming_dist_[cpt_] = hamming_dist;
     dist_next_prot_[cpt_] = dist_next_prot;*/
    pro_list_[cpt_].start_pos_ = start_pos;
    pro_list_[cpt_].prot_id_ = prot_id;
    pro_list_[cpt_].prot_parent_id_ = prot_parent_id;
    pro_list_[cpt_].length_ = length;
    pro_list_[cpt_].basal_level_ = basal_level;
    pro_list_[cpt_].hamming_dist_ = hamming_dist;
    pro_list_[cpt_].dist_next_prot_ = dist_next_prot;
    cpt_++;
  }
  
  void modify_protein(int32_t start_pos, int32_t length, double basal_level, int32_t hamming_dist, int32_t dist_next_prot, int32_t prot_id, int32_t prot_parent_id) {
    int index;
    for(int i = 0; i < pro_list_.size(); i++){
      if (pro_list_[i].prot_id_ == prot_id)
        index = i;
    }
    /*start_pos_[index] = start_pos;
     length_[index] = length;
     basal_level_[index] = basal_level;
     hamming_dist_[index] = hamming_dist;
     dist_next_prot_[index] = dist_next_prot;*/
    pro_list_[index].start_pos_ = start_pos;
    pro_list_[index].prot_parent_id_ = prot_parent_id;
    pro_list_[index].length_ = length;
    pro_list_[index].basal_level_ = basal_level;
    pro_list_[index].hamming_dist_ = hamming_dist;
    pro_list_[index].dist_next_prot_ = dist_next_prot;
  }
  
  void remove_protein(int32_t start_pos, int32_t length) {
    int index;
    for(int i = 0; i < pro_list_.size(); i++){
      if (pro_list_[i].start_pos_ == start_pos)
        index = i;
    }
    pro_list_.erase(pro_list_.begin() + index);
    /*start_pos_.erase(start_pos_.begin() + index);
     prot_id_.erase(prot_id_.begin() + index);
     length_.erase(length_.begin() + index);
     basal_level_.erase(basal_level_.begin() + index);
     hamming_dist_.erase(hamming_dist_.begin() + index);
     dist_next_prot_.erase(dist_next_prot_.begin() + index);*/
    cpt_--;
    nb_prot_--;
  }
  
  std::vector<protein_par> pro_list_;
  /*std::vector<int32_t> start_pos_;
   std::vector<std::string> prot_id_;
   std::vector<int32_t> length_;
   std::vector<double> basal_level_;
   std::vector<int32_t> hamming_dist_;
   std::vector<int32_t> dist_next_prot_;*/
  int nb_prot_;
  int cpt_ = 0;
  int dna_length_;
};



#endif // AEVOL_PROTEINLIST_H_
