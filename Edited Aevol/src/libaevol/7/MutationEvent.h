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


#ifndef AEVOL_MUTATIONEVENT_H
#define AEVOL_MUTATIONEVENT_H


#include <cstdint>

namespace aevol {

enum MutationEventType {
    DO_SWITCH         = 0,
    SMALL_INSERTION   = 1,
    SMALL_DELETION    = 2,
    DUPLICATION       = 3,
    DELETION          = 4,
    TRANSLOCATION     = 5,
    INVERSION         = 6
};

class MutationEvent {

 public:
    MutationEvent() = default;
    ~MutationEvent();

    void switch_pos(int32_t pos
    #ifdef BASE_4
    , char new_base
    #endif
    );

    void small_insertion(int32_t pos, int32_t number, char* seq);
    void small_deletion(int32_t pos, int32_t number);

    void duplication(int32_t pos1, int32_t pos2, int32_t pos3);

    void translocation(int32_t pos1, int32_t pos2, int32_t pos3,
                       int32_t pos4, bool invert);

    void inversion(int32_t pos1, int32_t pos2);

    void deletion(int32_t pos1, int32_t pos2);

    int32_t type() { return type_; };

    int32_t pos_1() { return pos_1_; }
    int32_t pos_2() { return pos_2_; }
    int32_t pos_3() { return pos_3_; }
    int32_t pos_4() { return pos_4_; }

    int32_t number() { return number_; }

    int32_t invert() { return invert_; }

    char* seq() { return seq_; }

    #ifdef BASE_4
    char base() { return base_; }
    #endif

 private:
    int32_t type_;

    int32_t pos_1_,pos_2_,pos_3_,pos_4_;

    int32_t number_; // insertion or deletion

    bool invert_;

    char* seq_ = nullptr;

    #ifdef BASE_4
    char base_; // switch
    #endif
};
}


#endif //RAEVOL_CUDA_MUTATIONEVENT_H
