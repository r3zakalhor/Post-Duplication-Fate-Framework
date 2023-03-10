//
// Created by arrouan on 19/11/15.
//

#ifndef AEVOL_REGUL_HABITATFACTORY_H
#define AEVOL_REGUL_HABITATFACTORY_H

#ifndef __REGUL
#include "Habitat.h"
#else
#include "raevol/Habitat_R.h"
#endif

#if __cplusplus == 201103L
#include "make_unique.h"
#endif

namespace aevol {

class HabitatFactory {
 public:
#ifndef __REGUL
  static std::unique_ptr<Habitat> create_unique_habitat(Habitat& habitat,
                                                        bool share_phenotypic_target);
  static std::unique_ptr<Habitat>create_unique_habitat(gzFile backup_file,
                            std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler);
#else
  static std::unique_ptr<Habitat_R> create_unique_habitat(Habitat_R& habitat,
														  bool share_phenotypic_target);
  static std::unique_ptr<Habitat_R> create_unique_habitat(gzFile backup_file,
                 					  std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler);
#endif

};

}
#endif //AEVOL_REGUL_HABITATFACTORY_H
