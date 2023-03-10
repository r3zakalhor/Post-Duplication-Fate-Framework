//
// Created by duazel on 16/04/2020.
//

#include "IOJson.h"

#include "ExpManager.h"
#include "Individual.h"
#include "JumpingMT.h"
#include "MutationParams.h"
#include "ParamLoader.h"
#include "json.hpp"

#include <algorithm>
#include <ctime>
#include <fstream>
#include <io/ParamLoader.h>
#include <iostream>
#include <population/IndividualFactory.h>

using namespace aevol;

static const int8_t STRAIN_NAME_DEFAULT_SIZE  = 20;
static const int8_t STRAIN_NAME_LOGIN_SIZE    = 10;


void IOJson::init() {
  // ----------------------------------------- PseudoRandom Number Generators
  seed_           = 0;
  mut_seed_       = 0;
  stoch_seed_     = 0;
  env_var_seed_   = 0;
  env_noise_seed_ = 0;

  // ------------------------------------------------------------ Constraints
  min_genome_length_  = 1;
  max_genome_length_  = 10000000;
  max_triangle_width_ = 0.033333333;

  // ----------------------------------------------------- Initial conditions
  chromosome_initial_length_  = 5000;
  init_method_            = ONE_GOOD_GENE | CLONE;
  init_pop_size_          = 1024;

  // ------------------------------------------------------------- Strain name
  char* login_name = new char[LOGIN_NAME_MAX+1];
  // Try get user login. If fail, replace by default value
  if(getlogin_r(login_name, LOGIN_NAME_MAX) != 0)
    strcpy(login_name, "anon");

  // Copy login into strain name with at most STRAIN_NAME_LOGIN_SIZE characters
  strain_name_ = login_name;
  //strncpy(strain_name_, login_name, STRAIN_NAME_LOGIN_SIZE);
  delete [] login_name;

  // Null-terminate the c-string if the max number of characters were copied
  if (strain_name_[STRAIN_NAME_LOGIN_SIZE] != 0)
    strain_name_[STRAIN_NAME_LOGIN_SIZE + 1] = 0;

  // -------------------------------------------------------- Phenotypic target
  env_sampling_ = 300;

  // ------------------------------------ Phenotypic target x-axis segmentation
  env_axis_nb_segments_         = 1;
  env_axis_segment_boundaries_  = NULL;
  env_axis_features_            = NULL;
  env_axis_separate_segments_   = false;

  // ---------------------------------------------- Phenotypic target variation
  env_var_method_ = NO_VAR;
  env_var_sigma_  = 0;
  env_var_tau_    = 0;

  // ----------------------------------------------- Phenotypic Stochasticity
  with_stochasticity_ = false;

  // -------------------------------------------------- Phenotypic target noise
  env_noise_method_       = NO_NOISE;
  env_noise_alpha_        = 0;
  env_noise_sigma_        = 0;
  env_noise_prob_         = 0;
  env_noise_sampling_log_ = 0;

  // --------------------------------------------------------- Mutation rates
  point_mutation_rate_  = 1e-5;
  small_insertion_rate_ = 1e-5;
  small_deletion_rate_  = 1e-5;
  max_indel_size_       = 6;
  with_4pts_trans_            = true;

  // -------------------------------------------- Rearrangements and Transfer
  with_HT_                    = false;
  repl_HT_with_close_points_  = false;
  HT_ins_rate_                = 0.0;
  HT_repl_rate_               = 0.0;
  repl_HT_detach_rate_        = 0.0;

  // ------------------------------ Rearrangement rates (without alignements)
  duplication_rate_   = 1e-5;
  deletion_rate_      = 1e-5;
  translocation_rate_ = 1e-5;
  inversion_rate_     = 1e-5;

  // --------------------------------- Rearrangement rates (with alignements)
  with_alignments_          = false;
  neighbourhood_rate_       = 5e-5;
  duplication_proportion_   = 0.3;
  deletion_proportion_      = 0.3;
  translocation_proportion_ = 0.3;
  inversion_proportion_     = 0.3;

  // ------------------------------------------------------------ Alignements
  align_fun_shape_    = SIGMOID;
  align_sigm_lambda_  = 4;
  align_sigm_mean_    = 50;
  align_lin_min_      = 0;
  align_lin_max_      = 100;

  align_max_shift_      = 20;
  align_w_zone_h_len_   = 50;
  align_match_bonus_    = 1;
  align_mismatch_cost_  = 2;

  // -------------------------------------------------------------- Selection
  selection_scheme_   = RANK_EXPONENTIAL;
  selection_pressure_ = 0.998;

  selection_scope_   = SCOPE_LOCAL;
  selection_scope_x_ = 3;
  selection_scope_y_ = 3;

  fitness_function_ = FITNESS_EXP;
  fitness_function_x_ = 3;
  fitness_function_y_ = 3;

  // -------------------------------------------------------------- Secretion
  with_secretion_               = false;
  secretion_contrib_to_fitness_ = 0;
  secretion_diffusion_prop_     = 0;
  secretion_degradation_prop_   = 0;
  secretion_cost_               = 0;
  secretion_init_               = 0;

  // --------------------------------------------------------------- Plasmids
  allow_plasmids_             = false;
  plasmid_initial_length_     = -1;
  plasmid_initial_gene_       = 0;
  plasmid_minimal_length_     = -1;
  plasmid_maximal_length_     = -1;
  chromosome_minimal_length_  = -1;
  chromosome_maximal_length_  = -1;
  prob_plasmid_HT_            = 0;
  tune_donor_ability_         = 0;
  tune_recipient_ability_     = 0;
  donor_cost_                 = 0;
  recipient_cost_             = 0;
  compute_phen_contrib_by_GU_ = false;
  swap_GUs_         = false;

  // ------------------------------------------------------- Translation cost
  translation_cost_ = 0;

  // ---------------------------------------------------------------- Outputs
  stats_            = 0;
  delete_old_stats_ = false;

  // Backups
  backup_step_      = 500;
  big_backup_step_  = 10000;

  // Tree
  record_tree_  = false;
  tree_step_    = 100;

  //LightTree
  record_light_tree_ = false;

  // Dumps
  make_dumps_ = false;
  dump_step_  = 1000;

  // Logs
  logs_ = 0;

  // Other
  more_stats_ = false;

  fuzzy_flavor_ = 0;

  world_width_  = 32;
  world_heigth_ = 32;
  well_mixed = false;
  partial_mix_nb_permutations_ = 0;
}

IOJson::IOJson(const std::string & filename) {
  std::ifstream file(filename);
  json input_file_;
  file >> input_file_;
  init();

  setStrainName(input_file_["param_in"].value("strain_name", "default"));
  setSeed(input_file_["param_in"]["rng"]["seed"]);
  setInitPopSize(input_file_["param_in"]["initial_conditions"]["init_pop_size"]);
  setWorldHeigth(input_file_["param_in"]["world_size"][0]);
  setWorldWidth(input_file_["param_in"]["world_size"][1]);

  //setChromosomeInitialLength(input_file_["param_in"]["initial_conditions"]["chromosome_inital_length"]);
  setMinGenomeLength(input_file_["param_in"]["constraints"]["min_genome_length"]);
  setFuzzyFlavor(input_file_["param_in"]["fuzzy_flavor"]);

  /*json::array_t json_init_method = input_file_["param_in"]["initial_conditions"]["init_method"];
  int init_method_tmp = 0;
  for (const auto & item : json_init_method) {
    std::string item_str = item;
    std::transform(item_str.begin(), item_str.end(), item_str.begin(),::toupper);
    if (strcmp(item_str.c_str(), "ONE_GOOD_GENE") == 0)
    {
      init_method_tmp |= ONE_GOOD_GENE;
    }
    else if (strcmp(item_str.c_str(), "CLONE") == 0)
    {
      init_method_tmp |= CLONE;
    }
    else if (strcmp(item_str.c_str(), "WITH_INS_SEQ") == 0)
    {
      init_method_tmp |= WITH_INS_SEQ;
    }
  }
  setInitMethod(init_method_tmp);*/


  setPointMutationRate(input_file_["param_in"]["mutations"].value("point_mutation_rate", 1e-7));
  setSmallDeletionRate(input_file_["param_in"]["mutations"].value("small_deletion_rate", 1e-7));
  setSmallInsertionRate(input_file_["param_in"]["mutations"].value("small_insertion_rate", 1e-7));
  setDeletionRate(input_file_["param_in"]["mutations"].value("deletion_rate", 1e-7));
  setDuplicationRate(input_file_["param_in"]["mutations"].value("duplication_rate", 1e-7));
  setTranslocationRate(input_file_["param_in"]["mutations"].value("translocation_rate", 1e-7));
  setInversionRate(input_file_["param_in"]["mutations"].value("inversion_rate", 1e-7));
  setWithAlignments(input_file_["param_in"]["mutations"].value("with_aligmnents", false));

  setEnvSampling(input_file_["param_in"]["env"].value("env_sampling",300));
  setEnvNoiseAlpha(input_file_["param_in"]["env"]["noise"].value("env_noise_alpha",0));
  setEnvNoiseProb(input_file_["param_in"]["env"]["noise"].value("env_noise_prob",0));
  setEnvNoiseSamplingLog(input_file_["param_in"]["env"]["noise"].value("env_noise_sampling_log",0));
  setEnvNoiseSigma(input_file_["param_in"]["env"]["noise"].value("env_noise_sigma",0));
  //setEnvVarMethod(input_file_["param_in"]["env"]["variation"].value("env_var_method",NO_VAR));
  setEnvVarSeed(input_file_["param_in"]["env"]["variation"].value("env_var_seed",0));
  setEnvVarSigma(input_file_["param_in"]["env"]["variation"].value("env_var_sigma",0));
  setEnvVarTau(input_file_["param_in"]["env"]["variation"].value("env_var_tau",0));
  setEnvNoiseSeed(input_file_["param_in"]["env"]["noise"].value("env_noise_seed",0));


  std::string selection_scheme_type_str = input_file_["param_in"]["selection"].value("selection_scheme_type", "default");
  SelectionScheme selection_scheme_tmp;
  if (selection_scheme_type_str == "RANK_EXPONENTIAL") {
    selection_scheme_tmp = SelectionScheme::RANK_EXPONENTIAL;
  } else if (selection_scheme_type_str == "RANK_LINEAR") {
    selection_scheme_tmp = SelectionScheme::RANK_LINEAR;
  } else if (selection_scheme_type_str == "FITNESS_PROPORTIONATE") {
    selection_scheme_tmp = SelectionScheme::FITNESS_PROPORTIONATE;
  } else if (selection_scheme_type_str == "FITTEST") {
    selection_scheme_tmp = SelectionScheme::FITTEST;
  }
  setSelectionScheme(selection_scheme_tmp);

  setSelectionPressure(input_file_["param_in"]["selection"].value("selection_pressure",0.998));
  setMaxIndelSize(input_file_["param_in"]["mutations"]["max_indel_size"]);

  setEnvSampling(input_file_["param_in"]["env"]["env_sampling"]);
  json::array_t json_env_add_gaussian = input_file_["param_in"]["env"]["env_add_gaussian"];
  std::list<Gaussian> env_gaussian_tmp;
  for (const auto & item : json_env_add_gaussian) {
    env_gaussian_tmp.emplace_back(item[0], item[1], item[2]);
  }
  setEnvAddGaussian(env_gaussian_tmp);
  //setMaxTriangleWidth(input_file_["param_in"]["max_triangle_width"]);
  setMaxTriangleWidth(input_file_["param_in"]["constraints"]["max_triangle_width"]);

  setBackupStep(input_file_["param_in"]["output"]["backup_step"]);
  setRecordTree(input_file_["param_in"]["output"].value("record_tree", false));
  //setMoreStats(input_file_["param_in"]["output"].value("more_stats", false));

  json individuals_json = input_file_["indivs"];

  for (auto & indiv : individuals_json) {
    uint32_t seed = input_file_["param_in"]["rng"].value("seed", seed_);
    uint32_t max_size = indiv.value("max_genome_length", max_genome_length_);
    bool allow_plasmids = indiv.value("allow_plasmids",allow_plasmids_);
    uint32_t id = indiv.value("id", 0);
    int minimum_size = indiv.value("min_genome_length", min_genome_length_);
    //double max_triangle_width = input_file_["param_in"]["constraints"]["max_triangle_width"];
    int age  = indiv.value("generation", 0);
    //int age  = indiv.value("generation", 0);
    std::string strain_name = indiv.value("strain_name", strain_name_);
    auto mut_prng = new JumpingMT(seed);
    auto stoch_prng = new JumpingMT(seed);
    std::shared_ptr<JumpingMT> mut_prng_ptr(mut_prng);
    std::shared_ptr<JumpingMT> stoch_prng_ptr(stoch_prng);
    auto *mut_param = new MutationParams();
    std::shared_ptr<MutationParams> mut_params_ptr(mut_param);
    mut_param->set_duplication_rate(indiv.value("duplication_rate",duplication_rate_));
    mut_param->set_deletion_rate(indiv.value("deletion_rate",deletion_rate_));
    mut_param->set_translocation_rate(indiv.value("translocation_rate",translocation_rate_));
    mut_param->set_small_insertion_rate(indiv.value("small_insertion_rate",small_insertion_rate_));
    mut_param->set_small_deletion_rate(indiv.value("small_deletion_rate",small_deletion_rate_));
    mut_param->set_inversion_rate(indiv.value("inversion_rate",inversion_rate_));
    mut_param->set_point_mutation_rate(indiv.value("point_mutation_rate",point_mutation_rate_));
    mut_param->set_max_indel_size(indiv.value("max_indel_size",max_indel_size_));
    Individual *individual = new Individual(new ExpManager(), mut_prng_ptr, stoch_prng_ptr, mut_params_ptr, max_triangle_width_,
                          minimum_size, max_size, allow_plasmids, id, strain_name.c_str(), age);
    auto *expSetup = new ExpSetup(nullptr);
    FuzzyFactory::fuzzyFactory = new FuzzyFactory(expSetup);


    int x = indiv.value("x_pos", 0);
    int y = indiv.value("y_pos", 0);

    std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler = std::make_shared<PhenotypicTargetHandler>();
    phenotypic_target_handler->set_sampling(input_file_["param_in"]["env"].value("env_sampling",300));
    phenotypic_target_handler->set_noise_alpha(input_file_["param_in"]["env"]["noise"].value("env_noise_alpha",0));
    phenotypic_target_handler->set_noise_prob(input_file_["param_in"]["env"]["noise"].value("env_noise_prob",0));
    phenotypic_target_handler->set_noise_sampling_log(input_file_["param_in"]["env"]["noise"].value("env_noise_sampling_log",0));
    uint32_t noise_seed = input_file_["param_in"]["env"]["noise"].value("env_noise_seed",0);
    std::shared_ptr<JumpingMT> noise_prng = std::make_shared<JumpingMT>(noise_seed);
    phenotypic_target_handler->set_noise_prng(noise_prng);
    phenotypic_target_handler->set_noise_sigma(input_file_["param_in"]["env"]["noise"].value("env_noise_sigma",0));
    uint32_t var_seed = input_file_["param_in"]["env"]["variation"].value("env_var_seed",0);
    std::shared_ptr<JumpingMT> var_prng = std::make_shared<JumpingMT>(var_seed);
    phenotypic_target_handler->set_var_prng(var_prng);
    phenotypic_target_handler->set_var_sigma(input_file_["param_in"]["env"]["variation"].value("env_var_sigma",0));
    phenotypic_target_handler->set_var_tau(input_file_["param_in"]["env"]["variation"].value("env_var_tau",0));

#ifdef __REGUL
    Habitat_R* habitat = new Habitat_R();
    habitat->set_phenotypic_target_handler(phenotypic_target_handler);
    std::unique_ptr<Habitat_R> habitat_ptr = std::unique_ptr<Habitat_R>(habitat);
    GridCell* grid_cell = new GridCell(x,y,std::move(habitat_ptr),individual,mut_prng_ptr,stoch_prng_ptr);

#else
    Habitat* habitat = new Habitat();
    habitat->set_phenotypic_target_handler(phenotypic_target_handler);
    std::unique_ptr<Habitat> habitat_ptr = std::unique_ptr<Habitat>(habitat);
    GridCell* grid_cell = new GridCell(x,y,std::move(habitat_ptr),individual,mut_prng_ptr,stoch_prng_ptr);
    grid_cell->set_individual(individual);
    individual->set_grid_cell(grid_cell);
    int16_t height = input_file_["param_in"]["world_size"][0];
    int16_t width = input_file_["param_in"]["world_size"][1];
    std::shared_ptr<JumpingMT> prng_ptr = std::make_shared<JumpingMT>(seed);
    individual->exp_m()->InitializeWorld(width,height,prng_ptr,mut_prng_ptr,stoch_prng_ptr,*habitat,true);

    json gu_list = indiv["GU"];
    for (auto & a_gu : gu_list) {
      std::string genome_str = a_gu.value("seq", "test");
      ;
      //a_gu.value("seq","0").c_str();
      char *genome = new char[genome_str.length() + 1];
      strcpy(genome, genome_str.c_str());
      individual->add_GU(genome, strlen(genome));
      individual->genetic_unit_nonconst(0).set_min_gu_length(0);
      individual->genetic_unit_nonconst(0).set_max_gu_length(1000000);
    }
    addIndividual(individual,gu_list);
#endif
  }
}

IOJson::IOJson(const std::string &param_in, const std::string &chromosome) {
  ParamLoader paramLoader(param_in.c_str());
  setStrainName(paramLoader.strain_name_);
  setSeed(paramLoader.seed_);
  setMutSeed(paramLoader.mut_seed_);
  setStochSeed(paramLoader.stoch_seed_);
  setEnvNoiseSeed(paramLoader.env_noise_seed_);
  setEnvVarSeed(paramLoader.env_var_seed_);
  // ------------------------------------------------------------ Constraints
  setMinGenomeLength(paramLoader.min_genome_length_);
  setMaxGenomeLength(paramLoader.max_genome_length_);
  setMaxTriangleWidth(paramLoader.w_max_);

  // ----------------------------------------------------- Initial conditions
  setChromosomeInitialLength(chromosome_initial_length_);
  setInitMethod(paramLoader.init_method_);
  setInitPopSize(paramLoader.init_pop_size_);
  setEnvSampling(paramLoader.env_sampling_);

  // ------------------------------------ Phenotypic target x-axis segmentation
  setEnvAxisNbSegments(paramLoader.env_axis_nb_segments_);
  setEnvAxisSegmentBoundaries(paramLoader.env_axis_segment_boundaries_);
  setEnvAxisSeparateSegments(paramLoader.env_axis_separate_segments_);
  setEnvAxisFeatures(paramLoader.env_axis_features_);

  // ---------------------------------------------- Phenotypic target variation
  setEnvVarMethod(paramLoader.env_var_method_);
  setEnvVarSigma(paramLoader.env_var_sigma_);
  setEnvVarTau(paramLoader.env_var_tau_);

  // ----------------------------------------------- Phenotypic Stochasticity
  setWithStochasticity(paramLoader.with_stochasticity_);

  // -------------------------------------------------- Phenotypic target noise
  setEnvNoiseMethod(paramLoader.env_noise_method_);
  setEnvNoiseAlpha(paramLoader.env_noise_alpha_);
  setEnvNoiseSigma(paramLoader.env_noise_sigma_);
  setEnvNoiseProb(paramLoader.env_noise_prob_);
  setEnvNoiseSamplingLog(paramLoader.env_noise_sampling_log_);

  // --------------------------------------------------------- Mutation rates
  setPointMutationRate(paramLoader.point_mutation_rate_);
  setSmallInsertionRate(paramLoader.small_insertion_rate_);
  setSmallDeletionRate(paramLoader.small_deletion_rate_);
  setMaxIndelSize(paramLoader.max_indel_size_);
  setWith4PtsTrans(paramLoader.with_4pts_trans_);

  // -------------------------------------------- Rearrangements and Transfer
  setWithHt(paramLoader.with_HT_);
  setReplHtWithClosePoints(paramLoader.repl_HT_with_close_points_);
  setHtInsRate(paramLoader.HT_ins_rate_);
  setHtReplRate(paramLoader.HT_repl_rate_);
  setReplHtDetachRate(paramLoader.repl_HT_detach_rate_);

  // ------------------------------ Rearrangement rates (without alignements)
  setDuplicationRate(paramLoader.duplication_rate_);
  setDeletionRate(paramLoader.deletion_rate_);
  setTranslocationRate(paramLoader.translocation_rate_);
  setInversionRate(paramLoader.inversion_rate_);

  // --------------------------------- Rearrangement rates (with alignements)
  setWithAlignments(paramLoader.with_alignments_);
  setNeighbourhoodRate(paramLoader.neighbourhood_rate_);
  setDuplicationProportion(paramLoader.duplication_proportion_);
  setDeletionProportion(paramLoader.deletion_proportion_);
  setTranslocationProportion(paramLoader.translocation_proportion_);
  setInversionProportion(paramLoader.inversion_proportion_);

  // ------------------------------------------------------------ Alignements
  setAlignFunShape(paramLoader.align_fun_shape_);
  setAlignSigmLambda(paramLoader.align_sigm_lambda_);
  setAlignSigmMean(paramLoader.align_sigm_mean_);
  setAlignLinMax(paramLoader.align_lin_max_);
  setAlignLinMin(paramLoader.align_lin_min_);

  setAlignMaxShift(paramLoader.align_max_shift_);
  setAlignWZoneHLen(paramLoader.align_w_zone_h_len_);
  setAlignMatchBonus(paramLoader.align_match_bonus_);
  setAlignMismatchCost(paramLoader.align_mismatch_cost_);

  // -------------------------------------------------------------- Selection
  setSelectionScheme(paramLoader.selection_scheme_);
  setSelectionPressure(paramLoader.selection_pressure_);
  setSelectionScope(paramLoader.selection_scope_);
  setSelectionScopeX(paramLoader.selection_scope_x_);
  setSelectionScopeY(paramLoader.selection_scope_y_);

  setFitnessFunction(paramLoader.fitness_function_);
  setFitnessFunctionX(paramLoader.fitness_function_x_);
  setFitnessFunctionY(paramLoader.fitness_function_y_);

  // -------------------------------------------------------------- Secretion
  setWithSecretion(paramLoader.with_secretion_);
  setSecretionContribToFitness(paramLoader.secretion_contrib_to_fitness_);
  setSecretionDiffusionProp(paramLoader.secretion_diffusion_prop_);
  setSecretionDegradationProp(paramLoader.secretion_degradation_prop_);
  setSecretionCost(paramLoader.secretion_cost_);
  setSecretionInit(paramLoader.secretion_init_);

  // --------------------------------------------------------------- Plasmids
  setAllowPlasmids(paramLoader.allow_plasmids_);
  setPlasmidInitialGene(paramLoader.plasmid_initial_gene_);
  setPlasmidInitialLength(paramLoader.plasmid_initial_length_);
  setPlasmidMaximalLength(paramLoader.plasmid_maximal_length_);
  setPlasmidMinimalLength(paramLoader.plasmid_minimal_length_);
  setChromosomeInitialLength(paramLoader.chromosome_initial_length_);
  setChromosomeMaxLength(paramLoader.chromosome_maximal_length_);
  setChromosomeMinLength(paramLoader.chromosome_minimal_length_);
  setProbPlasmidHt(paramLoader.prob_plasmid_HT_);
  setTuneDonorAbility(paramLoader.tune_donor_ability_);
  setTuneRecipientAbility(paramLoader.tune_recipient_ability_);
  setDonorCost(paramLoader.donor_cost_);
  setRecipientCost(paramLoader.recipient_cost_);
  setComputePhenContribByGu(paramLoader.compute_phen_contrib_by_GU_);
  setSwapGUs(paramLoader.swap_GUs_);

  // ------------------------------------------------------- Translation cost
  setTranslationCost(paramLoader.translation_cost_);

  // ---------------------------------------------------------------- Outputs
  setStats(paramLoader.stats_);
  setDeleteOldStats(paramLoader.delete_old_stats_);

  // Backups
  setBackupStep(paramLoader.backup_step_);
  setBigBackupStep(paramLoader.big_backup_step_);

  // Tree
  setRecordTree(paramLoader.record_tree_);
  setTreeStep(paramLoader.tree_step_);

  //LightTree
  setRecordLightTree(paramLoader.record_light_tree_);

  // Dumps
  setMakeDumps(paramLoader.make_dumps_);
  setDumpStep(paramLoader.dump_step_);

  // Logs
  setLogs(paramLoader.logs_);

  // Other
  setMoreStats(paramLoader.more_stats_);

  setFuzzyFlavor(paramLoader._fuzzy_flavor);

  setWorldWidth(paramLoader.grid_width_);
  setWorldHeigth(paramLoader.grid_height_);
  setWellMixed(paramLoader.well_mixed);
  setPartialMixNbPermutations(paramLoader.partial_mix_nb_permutations);

  std::ifstream sequence_file(chromosome);
  std::string sequence;
  sequence_file >> sequence;
  json indivs = json::array();
  json first_indiv;
  json gen_unit;
  gen_unit["seq"] = sequence;
  first_indiv["GU"] = gen_unit;
  indivs[0] =first_indiv;
  json_file_["indivs"] = indivs;

  std::list<Gaussian> gaussian_vector;
  for (auto & item : paramLoader.std_env_gaussians) {
    gaussian_vector.emplace_back(item.height(),item.mean(),item.width());
  }
  setEnvAddGaussian(gaussian_vector);
}

IOJson::IOJson(ExpManager * exp_m) {

  ExpSetup* exp_s = exp_m->exp_s();
  Selection* sel = exp_m->sel();
  OutputManager* output_m = exp_m->output_m();
  output_m->InitStats();
  World* world = exp_m->world();
  std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler = world->phenotypic_target_handler();
  std::shared_ptr<aevol::MutationParams> param_mut = exp_s->mut_params();

  init();

  // ---------------------------------------------------------------- Selection
  setSelectionScheme(sel->selection_scheme());
  setSelectionPressure(sel->selection_pressure());

  setSelectionScope(sel->selection_scope());
  setSelectionScopeX(sel->selection_scope_x());
  setSelectionScopeX(sel->selection_scope_y());

  setFitnessFunction(sel->fitness_func());
  setFitnessFunctionX(sel->fitness_function_scope_x());
  setFitnessFunctionY(sel->fitness_function_scope_y());
  // ----------------------------------------------------------------- Transfer
  setWithHt(exp_s->with_HT());
  setReplHtWithClosePoints(exp_s->repl_HT_with_close_points());
  setHtInsRate(exp_s->HT_ins_rate());
  setHtReplRate(exp_s->HT_repl_rate());
  setReplHtDetachRate(exp_s->repl_HT_detach_rate());

  // ----------------------------------------------------------------- Plasmids
  setAllowPlasmids(exp_s->with_plasmids());
  setProbPlasmidHt(exp_s->prob_plasmid_HT());
  setTuneDonorAbility(exp_s->tune_donor_ability());
  setTuneRecipientAbility(exp_s->tune_recipient_ability());
  setDonorCost(exp_s->donor_cost());
  setRecipientCost(exp_s->recipient_cost());
  setSwapGUs(exp_s->swap_GUs());
  setComputePhenContribByGu(output_m->compute_phen_contrib_by_GU());
  setPartialMixNbPermutations(world->partial_mix_nb_permutations());

  // ---------------------------------------------------------------- Secretion
  setWithSecretion(exp_s->with_secretion());
  setSecretionContribToFitness(exp_s->secretion_contrib_to_fitness());
  setSecretionCost(exp_s->secretion_cost());
  setSecretionDegradationProp(world->secretion_degradation_prop());
  setSecretionDiffusionProp(world->secretion_diffusion_prop());

  setFuzzyFlavor(exp_s->get_fuzzy_flavor());

  //------------------------------------------------------------------ Constraints
  setMinGenomeLength(exp_s->min_genome_length());
  setMaxGenomeLength(exp_s->max_genome_length());
  setMaxTriangleWidth(world->indiv_at(0, 0)->w_max());

  //------------------------------------------------------------------ Alignements
  setDeletionProportion(param_mut->deletion_proportion());
  setDuplicationProportion(param_mut->duplication_proportion());
  setInversionProportion(param_mut->inversion_proportion());
  setTranslocationProportion(param_mut->translocation_proportion());
  setNeighbourhoodRate(param_mut->neighbourhood_rate());
  setWithAlignments(param_mut->with_alignments());

  //------------------------------------------------------------------ Mutations
  setDeletionRate(param_mut->deletion_rate());
  setDuplicationRate(param_mut->duplication_rate());
  setInversionRate(param_mut->inversion_rate());
  setMaxIndelSize(param_mut->max_indel_size());
  setPointMutationRate(param_mut->point_mutation_rate());
  setSmallDeletionRate(param_mut->small_deletion_rate());
  setSmallInsertionRate(param_mut->small_insertion_rate());
  setTranslocationRate(param_mut->translocation_rate());
  setWith4PtsTrans(param_mut->with_4pts_trans());

  //------------------------------------------------------------------ Env
  setEnvSampling(phenotypic_target_handler->sampling());
  setEnvNoiseAlpha(phenotypic_target_handler->noise_alpha());
  setEnvNoiseProb(phenotypic_target_handler->noise_prob());
  setEnvNoiseSamplingLog(phenotypic_target_handler->noise_sampling_log());
  setEnvNoiseSeed(phenotypic_target_handler->noise_prng().use_count());
  setEnvNoiseSigma(phenotypic_target_handler->noise_sigma());

  //------------------------------------------------------------------ Variation
  setEnvVarMethod(phenotypic_target_handler->var_method());
  setEnvVarSeed(phenotypic_target_handler->var_prng().use_count());
  setEnvVarSigma(phenotypic_target_handler->var_sigma());
  setEnvVarTau(phenotypic_target_handler->var_tau());



  setRecordTree(output_m->record_tree());
  setRecordLightTree(output_m->record_light_tree());
  //setTreeStep((int) output_m->tree_step());

  setEnvAddGaussian(phenotypic_target_handler->current_gaussians());
  setWorldWidth(world->width());
  setWorldHeigth(world->height());
  setInitPopSize(world->nb_indivs());

  setBackupStep(output_m->backup_step());
  setBigBackupStep(output_m->big_backup_step());

  setStrainName(exp_m->best_indiv()->strain_name());
  setAllowPlasmids(exp_m->best_indiv()->allow_plasmids());

  setSeed(world->prng()->initial_seed());

  json_file_["indivs"];

}

void IOJson::load(ExpManager * exp_m, bool verbose,
                  char* chromosome, int32_t lchromosome,
                  char* plasmid, int32_t lplasmid) {
  // Check consistency of min, max and initial length of chromosome and plasmid
  // Default for by GU minimal or maximal size is -1.
  // If equal to -1, maximal sizes of each GU will be replaced by total maximal size for the whole genome
//  CheckConsistency();

  // Initialize prng_
  // This one will be used to create the initial genome(s) and to generate seeds for other prng
  prng_ = std::make_shared<JumpingMT>(seed_);

  // Initialize mut_prng, stoch_prng, world_prng :
  // if mut_seed (respectively stoch_seed) not given in param.in, choose it at random
  if (mut_seed_ == 0) {
    mut_seed_ = prng_->random(1000000);
  }
  if (stoch_seed_ == 0) {
    stoch_seed_ = prng_->random(1000000);
  }
  auto mut_prng   = std::make_shared<JumpingMT>(mut_seed_);
  auto stoch_prng = std::make_shared<JumpingMT>(stoch_seed_);
  auto world_prng = std::make_shared<JumpingMT>(prng_->random(1000000));

  // Create aliases
  ExpSetup* exp_s = exp_m->exp_s();
  Selection* sel = exp_m->sel();
  OutputManager* output_m = exp_m->output_m();
  output_m->InitStats();


  // 1) ------------------------------------- Initialize the experimental setup


  // ---------------------------------------------------------------- Selection
  sel->set_selection_scheme(selection_scheme_);
  sel->set_selection_pressure(selection_pressure_);

  sel->set_selection_scope(selection_scope_);
  sel->set_selection_scope_x(selection_scope_x_);
  sel->set_selection_scope_y(selection_scope_y_);

  sel->set_fitness_function(fitness_function_);
  sel->set_fitness_function_scope_x(fitness_function_x_);
  sel->set_fitness_function_scope_y(fitness_function_y_);
  // ----------------------------------------------------------------- Transfer
  exp_s->set_with_HT(with_HT_);
  exp_s->set_repl_HT_with_close_points(repl_HT_with_close_points_);
  exp_s->set_HT_ins_rate(HT_ins_rate_);
  exp_s->set_HT_repl_rate(HT_repl_rate_);
  exp_s->set_repl_HT_detach_rate(repl_HT_detach_rate_);

  // ----------------------------------------------------------------- Plasmids
  exp_s->set_with_plasmids(allow_plasmids_);
  exp_s->set_prob_plasmid_HT(prob_plasmid_HT_);
  exp_s->set_tune_donor_ability(tune_donor_ability_);
  exp_s->set_tune_recipient_ability(tune_recipient_ability_);
  exp_s->set_donor_cost(donor_cost_);
  exp_s->set_recipient_cost(recipient_cost_);
  exp_s->set_swap_GUs(swap_GUs_);
  output_m->set_compute_phen_contrib_by_GU(compute_phen_contrib_by_GU_);

  // ---------------------------------------------------------------- Secretion
  exp_s->set_with_secretion(with_secretion_);
  exp_s->set_secretion_contrib_to_fitness(secretion_contrib_to_fitness_);
  exp_s->set_secretion_cost(secretion_cost_);

  exp_s->set_fuzzy_flavor(fuzzy_flavor_);

  //------------------------------------------------------------------ Parameter for SIMD
  exp_s->set_min_genome_length(min_genome_length_);
  exp_s->set_max_genome_length(max_genome_length_);

#ifdef __REGUL
  printf("IOJSON Not working yet with RAevol\n");
  exit(-11);
#endif

  if (FuzzyFactory::fuzzyFactory == NULL)
    FuzzyFactory::fuzzyFactory = new FuzzyFactory(exp_s);

  // 2) --------------------------------------------- Create and init a Habitat
#ifndef __REGUL
  Habitat habitat;
#else
  Habitat_R habitat;
#endif

  // Shorthand for phenotypic target handler
#ifndef __REGUL
  PhenotypicTargetHandler& phenotypic_target_handler =
      habitat.phenotypic_target_handler_nonconst();
#else
  PhenotypicTargetHandler_R& phenotypic_target_handler =
      habitat.phenotypic_target_handler_nonconst();
#endif

  // Move the gaussian list from the parameters to the phen target handler
#ifndef __REGUL
  phenotypic_target_handler.set_gaussians(env_add_gaussian_);
#else
//  phenotypic_target_handler.set_gaussians(_env_gaussians_list);
//  phenotypic_target_handler.set_signals_models(_signals_models);
//  phenotypic_target_handler.set_signals(_env_signals_list);
#endif

  // Copy the sampling
  phenotypic_target_handler.set_sampling(env_sampling_);

  // Set phenotypic target segmentation

  if((env_axis_features_ != NULL) && (env_axis_segment_boundaries_ != NULL)) {
    // if param.in contained a line starting with ENV_AXIS_FEATURES,
    // we use the values indicated on this line
    phenotypic_target_handler.set_segmentation(env_axis_nb_segments_,
                                               env_axis_segment_boundaries_,
                                               env_axis_features_,
                                               env_axis_separate_segments_);
  }
  // else we leave the segmentation as it is by default
  // (one "metabolic" segment from X_MIN to X_MAX)


  // Set phenotypic target variation
  if (env_var_method_ != NO_VAR)
  {
    phenotypic_target_handler.set_var_method(env_var_method_);
    phenotypic_target_handler.set_var_prng(std::make_shared<JumpingMT>(env_var_seed_));
    phenotypic_target_handler.set_var_sigma_tau(env_var_sigma_, env_var_tau_);
#ifdef __REGUL
//    phenotypic_target_handler.set_switch_probability(_env_switch_probability);
#endif
  }

  // Set phenotypic target noise
  if (env_noise_method_ != NO_NOISE)
  {
    phenotypic_target_handler.set_noise_method(env_noise_method_);
    phenotypic_target_handler.set_noise_sampling_log(env_noise_sampling_log_);
    phenotypic_target_handler.set_noise_prng(std::make_shared<JumpingMT>(env_noise_seed_));
    phenotypic_target_handler.set_noise_alpha(env_noise_alpha_);
    phenotypic_target_handler.set_noise_sigma(env_noise_sigma_);
    phenotypic_target_handler.set_noise_prob(env_noise_prob_);
  }

  // Build the phenotypic target
#ifndef __REGUL
  phenotypic_target_handler.BuildPhenotypicTarget();
#else
//  printf("Init phenotypic target with %d\n",_nb_indiv_age);
//  phenotypic_target_handler.InitPhenotypicTargetsAndModels( _nb_indiv_age );
#endif

  if (verbose) {
#ifndef __REGUL
    printf("Entire geometric area of the phenotypic target : %f\n",
           phenotypic_target_handler.get_geometric_area());
#else
//    phenotypic_target_handler.print_geometric_areas();
#endif
  }


  // 3) --------------------------------------------- Create the new population
  list<Individual *> indivs;
  // Generate a model ae_mut_param object
  auto param_mut = std::make_shared<MutationParams>();
  param_mut->set_point_mutation_rate(point_mutation_rate_);
  param_mut->set_small_insertion_rate(small_insertion_rate_);
  param_mut->set_small_deletion_rate(small_deletion_rate_);
  param_mut->set_max_indel_size(max_indel_size_);
  param_mut->set_with_4pts_trans(with_4pts_trans_);
  param_mut->set_with_alignments(with_alignments_);
  param_mut->set_with_HT(with_HT_);
  param_mut->set_repl_HT_with_close_points(repl_HT_with_close_points_);
  param_mut->set_HT_ins_rate(HT_ins_rate_);
  param_mut->set_HT_repl_rate(HT_repl_rate_);
  param_mut->set_repl_HT_detach_rate(repl_HT_detach_rate_);
  param_mut->set_duplication_rate(duplication_rate_);
  param_mut->set_deletion_rate(deletion_rate_);
  param_mut->set_translocation_rate(translocation_rate_);
  param_mut->set_inversion_rate(inversion_rate_);
  param_mut->set_neighbourhood_rate(neighbourhood_rate_);
  param_mut->set_duplication_proportion(duplication_proportion_);
  param_mut->set_deletion_proportion(deletion_proportion_);
  param_mut->set_translocation_proportion(translocation_proportion_);
  param_mut->set_inversion_proportion(inversion_proportion_);

  exp_s->set_mutation_parameters(param_mut);

  Individual * indiv = NULL;
  int32_t id_new_indiv = 0;

  if (chromosome != NULL)
  {
    printf("Option -c is used: chromosome will be loaded from a text file\n");
#ifndef __REGUL
    printf("construction de indiv :\n");
    Individual * indiv = new Individual(exp_m,
                                        mut_prng,
                                        stoch_prng,
                                        param_mut,
                                        max_triangle_width_,
                                        min_genome_length_,
                                        max_genome_length_,
                                        allow_plasmids_,
                                        id_new_indiv++,
                                        strain_name_.data(),
                                        0);
#else
    double w_max_=0.333;
    Individual_R * indiv = new Individual_R(exp_m,
                                         mut_prng,
                                         stoch_prng,
                                         param_mut,
                                         w_max_,
                                         min_genome_length_,
                                         max_genome_length_,
                                         allow_plasmids_,
                                         id_new_indiv++,
                                         strain_name_.data(),
                                         0);

#endif

    indiv->add_GU(chromosome, lchromosome);
    indiv->genetic_unit_nonconst(0).set_min_gu_length(chromosome_minimal_length_);
    indiv->genetic_unit_nonconst(0).set_max_gu_length(chromosome_maximal_length_);
    if (plasmid != NULL)
    {
      printf("Option -p is used: plasmid will be loaded from a text file\n");
      if (! allow_plasmids_)
      {
        printf("ERROR: option -p requires ALLOW_PLASMIDS set to true\n");
        exit(EXIT_FAILURE);
      }
      indiv->add_GU(plasmid, lplasmid);
      indiv->genetic_unit_nonconst(1).set_min_gu_length(plasmid_minimal_length_);
      indiv->genetic_unit_nonconst(1).set_max_gu_length(plasmid_maximal_length_);
    }
    else if (allow_plasmids_)
    {
      printf("ERROR: if you use option -c and ALLOW_PLASMIDS is set to true, you must also use option -p. \n For now loading a genetic unit from text file and generating the other is not supported.\n");
      exit(EXIT_FAILURE);
    }
    indiv->set_with_stochasticity(with_stochasticity_);
    indiv->compute_statistical_data();
    indiv->EvaluateInContext(habitat);
    printf("Starting with a clonal population of individual with metabolic error %f and secretion error %f \n",indiv->dist_to_target_by_feature(METABOLISM),indiv->dist_to_target_by_feature(SECRETION));
    indivs.push_back(indiv);

    // Make the clones and add them to the list of individuals
    for (int32_t i = 1 ; i < init_pop_size_ ; i++)
    {
#ifndef __REGUL
      Individual * clone = Individual::CreateClone(indiv, id_new_indiv++);
#else
      Individual_R * clone = Individual_R::CreateClone(indiv, id_new_indiv++);
#endif
      clone->EvaluateInContext(habitat);
      indivs.push_back(clone);
    }
  }
  else if (plasmid != NULL)
  {
    printf("ERROR: option -p can only be used in combination with option -c for now\n");
    exit(EXIT_FAILURE);
  }
  else if (init_method_ & ONE_GOOD_GENE)
  {
    if (init_method_ & CLONE)
    {
      // Create an individual with a "good" gene (in fact, make an indiv whose
      // fitness is better than that corresponding to a flat phenotype)
      // and set its id
      indiv = IndividualFactory::create_random_individual(
          exp_m,
          id_new_indiv++,
          param_mut,
          mut_prng,
          stoch_prng,
          habitat,
          max_triangle_width_,
          min_genome_length_,
          max_genome_length_,
          chromosome_initial_length_,
          allow_plasmids_,
          plasmid_initial_gene_,
          plasmid_initial_length_,
          const_cast<char*>(strain_name_.data()),
          prng_,
          true);
      indiv->genetic_unit_nonconst(0).set_min_gu_length(chromosome_minimal_length_);
      indiv->genetic_unit_nonconst(0).set_max_gu_length(chromosome_maximal_length_);

      if (allow_plasmids_)
      {
        indiv->genetic_unit_nonconst(1).set_min_gu_length(plasmid_minimal_length_);
        indiv->genetic_unit_nonconst(1).set_max_gu_length(plasmid_maximal_length_);
      }

      indiv->set_with_stochasticity(with_stochasticity_);

      // Add it to the list
      indivs.push_back(indiv);

      // Make the clones and add them to the list of individuals
      for (int32_t i = 1 ; i < init_pop_size_ ; i++)
      {
        // Add new clone to the list
#ifndef __REGUL
        Individual * clone = Individual::CreateClone(indiv, id_new_indiv++);
#else
        Individual_R * clone = Individual_R::CreateClone(dynamic_cast<Individual_R*>(indiv), id_new_indiv++);
#endif
        clone->EvaluateInContext(habitat);
        indivs.push_back(clone);
      }
    }
    else // if (! CLONE)
    {
      for (int32_t i = 0 ; i < init_pop_size_ ; i++)
      {
        // Create an individual and set its id
        indiv = IndividualFactory::create_random_individual(
            exp_m,
            id_new_indiv++,
            param_mut,
            mut_prng,
            stoch_prng,
            habitat,
            max_triangle_width_,
            min_genome_length_,
            max_genome_length_,
            chromosome_initial_length_,
            allow_plasmids_,
            plasmid_initial_gene_,
            plasmid_initial_length_,
            const_cast<char*>(strain_name_.data()),
            prng_,
            true);
        indiv->genetic_unit_nonconst(0).set_min_gu_length(chromosome_minimal_length_);
        indiv->genetic_unit_nonconst(0).set_max_gu_length(chromosome_maximal_length_);
        if (allow_plasmids_)
        {
          indiv->genetic_unit_nonconst(1).set_min_gu_length(plasmid_minimal_length_);
          indiv->genetic_unit_nonconst(1).set_max_gu_length(plasmid_maximal_length_);
        }

        // Add it to the list
        indivs.push_back(indiv);
      }
    }
  }
  else // if (! ONE_GOOD_GENE)
  {
    if (init_method_ & CLONE)
    {
      // Create a random individual and set its id
      indiv = IndividualFactory::create_random_individual(
          exp_m,
          id_new_indiv++,
          param_mut,
          mut_prng,
          stoch_prng,
          habitat,
          max_triangle_width_,
          min_genome_length_,
          max_genome_length_,
          chromosome_initial_length_,
          allow_plasmids_,
          plasmid_initial_gene_,
          plasmid_initial_length_,
          const_cast<char*>(strain_name_.data()),
          prng_,
          false);
      indiv->genetic_unit_nonconst(0).set_min_gu_length(chromosome_minimal_length_);
      indiv->genetic_unit_nonconst(0).set_max_gu_length(chromosome_maximal_length_);
      if (allow_plasmids_)
      {
        indiv->genetic_unit_nonconst(1).set_min_gu_length(plasmid_minimal_length_);
        indiv->genetic_unit_nonconst(1).set_max_gu_length(plasmid_maximal_length_);
      }

      // Add it to the list
      indivs.push_back(indiv);

      // Make the clones and add them to the list of individuals
      for (int32_t i = 1 ; i < init_pop_size_ ; i++)
      {
        // Add clone to the list
#ifndef __REGUL
        Individual * clone = Individual::CreateClone(indiv, id_new_indiv++);
#else
        Individual_R * clone = Individual_R::CreateClone(dynamic_cast<Individual_R*>(indiv), id_new_indiv++);
#endif
        clone->EvaluateInContext(habitat);
        indivs.push_back(clone);
      }
    }
    else // if (! CLONE)
    {
      for (int32_t i = 0 ; i < init_pop_size_ ; i++)
      {
        // Create a random individual and set its id
        indiv = IndividualFactory::create_random_individual(
            exp_m,
            id_new_indiv++,
            param_mut,
            mut_prng,
            stoch_prng,
            habitat,
            max_triangle_width_,
            min_genome_length_,
            max_genome_length_,
            chromosome_initial_length_,
            allow_plasmids_,
            plasmid_initial_gene_,
            plasmid_initial_length_,
            const_cast<char*>(strain_name_.data()),
            prng_,
            false);
        indiv->genetic_unit_nonconst(0).set_min_gu_length(chromosome_minimal_length_);
        indiv->genetic_unit_nonconst(0).set_max_gu_length(chromosome_maximal_length_);
        if (allow_plasmids_)
        {
          indiv->genetic_unit_nonconst(1).set_min_gu_length(plasmid_minimal_length_);
          indiv->genetic_unit_nonconst(1).set_max_gu_length(plasmid_maximal_length_);
        }

        // Add it to the list
        indivs.push_back(indiv);
      }
    }
  }

  // -------------------------------------------------------- Spatial structure
  exp_m->InitializeWorld(world_width_, world_heigth_,
                         world_prng,mut_prng,stoch_prng,
                         habitat,
                         true);
  World* world = exp_m->world();

  for (int16_t x = 0; x < exp_m->grid_width(); x++) {
    for (int16_t y = 0; y < exp_m->grid_height(); y++) {
      int32_t seed = prng_->random(1000000);
#if __cplusplus == 201103L
      exp_m->world()->grid(x,y)->set_reprod_prng(make_unique<JumpingMT>(seed));
      exp_m->world()->grid(x,y)->set_reprod_prng_simd(make_unique<JumpingMT>(seed));
#else
      exp_m->world()->grid(x,y)->set_reprod_prng(std::make_unique<JumpingMT>(seed));
      exp_m->world()->grid(x,y)->set_reprod_prng_simd(std::make_unique<JumpingMT>(seed));
#endif
    }
  }

  world->set_secretion_degradation_prop(secretion_degradation_prop_);
  world->set_secretion_diffusion_prop(secretion_diffusion_prop_);
  world->set_is_well_mixed(well_mixed);
  world->set_partial_mix_nb_permutations(partial_mix_nb_permutations_);

 //Set each individual's position on the grid
  int16_t x, y;
  int16_t x_max = exp_m->grid_width();
  int16_t y_max = exp_m->grid_height();

  for (const auto& indiv: indivs) {
    do {
      x = exp_m->world()->prng()->random(x_max);
      y = exp_m->world()->prng()->random(y_max);
    } while (world->indiv_at(x, y) != NULL);

    world->PlaceIndiv(indiv, x, y, true);
  }

  world->set_best(0, 0);



  // 4) ------------------------------------------ Set the recording parameters
  output_m->set_backup_step(backup_step_);
  output_m->set_big_backup_step(big_backup_step_);

  if (record_tree_)
  {
    output_m->init_tree(exp_m, tree_step_);
  }

  output_m->init_light_tree(record_light_tree_,exp_m, tree_step_);

  if (make_dumps_)
  {
    output_m->set_dump_step(dump_step_);
  }
  output_m->set_logs(logs_);
}

uint32_t IOJson::getSeed() const { return seed_; }
void IOJson::setSeed(uint32_t seed) {
  seed_ = seed;
  json_file_["param_in"]["rng"]["seed"] = seed;
}

int IOJson::getInitPopSize() const { return init_pop_size_; }
void IOJson::setInitPopSize(int initPopSize) {
  init_pop_size_ = initPopSize;
  json_file_["param_in"]["initial_conditions"]["init_pop_size"] = initPopSize;
}

int IOJson::getWorldHeigth() const {return world_heigth_; }
void IOJson::setWorldHeigth(int worldHeigth) {
  world_heigth_ = worldHeigth;
  json_file_["param_in"]["world_size"][0] = worldHeigth;
}

int IOJson::getWorldWidth() const { return world_width_; }
void IOJson::setWorldWidth(int worldWidth) {
  world_width_ = worldWidth;
  json_file_["param_in"]["world_size"][1] = worldWidth;
}

int IOJson::getInitMethod() const {
  return init_method_;
}
void IOJson::setInitMethod(int initMethod) {
  init_method_ = initMethod;
  json::array_t json_array;
  if (initMethod & GenomeInitializationMethod::CLONE) {
    json_array.emplace_back("CLONE");
  }
  if (initMethod & GenomeInitializationMethod::ONE_GOOD_GENE) {
    json_array.emplace_back("ONE_GOOD_GENE");
  }
  if (initMethod & GenomeInitializationMethod::WITH_INS_SEQ) {
    json_array.emplace_back("WITH_INS_SEQ");
  }
  json_file_["param_in"]["initial_conditions"]["init_method"] = json_array;
}

int IOJson::getChromosomeInitialLength() const {
  return chromosome_initial_length_;
}
void IOJson::setChromosomeInitialLength(int chromosomeInitialLength) {
  chromosome_initial_length_ = chromosomeInitialLength;
  json_file_["param_in"]["initial_conditions"]["chromosome_inital_length"] = chromosomeInitialLength;
}

int IOJson::getChromosomeMinLength() const { return min_genome_length_; }
void IOJson::setChromosomeMinLength(int chromosome_min_length) {
  min_genome_length_ = chromosome_min_length;
  json_file_["param_in"]["plasmids"]["chromosome_min_length"] = chromosome_min_length;
}

int IOJson::getFuzzyFlavor() const { return fuzzy_flavor_; }
void IOJson::setFuzzyFlavor(int fuzzyFlavor) {
  fuzzy_flavor_ = fuzzyFlavor;
  json_file_["param_in"]["fuzzy_flavor"] = fuzzy_flavor_;
}

const SelectionScheme IOJson::getSelectionScheme() const {
  return selection_scheme_;
}

void IOJson::setSelectionScheme(const SelectionScheme selectionScheme) {
  selection_scheme_ = selectionScheme;

  switch (selectionScheme) {
  case SelectionScheme::RANK_EXPONENTIAL:
    json_file_["param_in"]["selection"]["selection_scheme_type"] = "RANK_EXPONENTIAL";
    break;
  case SelectionScheme::RANK_LINEAR:
    json_file_["param_in"]["selection"]["selection_scheme_type"] = "RANK_LINEAR";
    break;
  case SelectionScheme::FITNESS_PROPORTIONATE:
    json_file_["param_in"]["selection"]["selection_scheme_type"] = "FITNESS_PROPORTIONATE";
    break;
  case SelectionScheme::FITTEST:
    json_file_["param_in"]["selection"]["selection_scheme_type"] = "FITTEST";
    break;
  default:
    json_file_["param_in"]["selection"]["selection_scheme_type"] = "UNKNOWN";
    break;
  }
}

int IOJson::getEnvSampling() const { return env_sampling_; }
void IOJson::setEnvSampling(int envSampling) {
  env_sampling_ = envSampling;
  json_file_["param_in"]["env"]["env_sampling"] = envSampling;
}

const list<Gaussian> &IOJson::getEnvAddGaussian() const {
  return env_add_gaussian_;
}
void IOJson::setEnvAddGaussian(const list<Gaussian> &envAddGaussian) {
  env_add_gaussian_.clear();
  for (const auto & item: envAddGaussian) {
    env_add_gaussian_.emplace_back(item);
  }
  json::array_t json_array;
  for(const auto & item : envAddGaussian) {
    json::array_t gaussian_array;
    gaussian_array.emplace_back(item.height());
    gaussian_array.emplace_back(item.mean());
    gaussian_array.emplace_back(item.width());
    json_array.emplace_back(gaussian_array);
  }
  json_file_["param_in"]["env"]["env_add_gaussian"] = json_array;
}

double IOJson::getMaxTriangleWidth() const { return max_triangle_width_; }
void IOJson::setMaxTriangleWidth(double maxTriangleWidth) {
  max_triangle_width_                                         = maxTriangleWidth;
  json_file_["param_in"]["constraints"]["max_triangle_width"] = maxTriangleWidth;
}

int IOJson::getBackupStep() const { return backup_step_; }
void IOJson::setBackupStep(int backupStep) {
  backup_step_ = backupStep;
  json_file_["param_in"]["output"]["backup_step"] = backupStep;
}

bool IOJson::isRecordTree() const { return record_tree_; }
void IOJson::setRecordTree(bool recordTree) {
  record_tree_ = recordTree;
  json_file_["param_in"]["output"]["record_tree"] = record_tree_;
}

bool IOJson::isMoreStats() const { return more_stats_; }
void IOJson::setMoreStats(bool moreStats) {
  more_stats_ = moreStats;
  json_file_["param_in"]["output"]["more_stats"] = more_stats_;
}

double IOJson::getPointMutationRate() const { return point_mutation_rate_; }
void IOJson::setPointMutationRate(double pointMutationRate) {
  point_mutation_rate_ = pointMutationRate;
  json_file_["param_in"]["mutations"]["point_mutation_rate"] = point_mutation_rate_;
}
double IOJson::getSmallInsertionRate() const { return small_insertion_rate_; }
void IOJson::setSmallInsertionRate(double smallInsertionRate) {
  small_insertion_rate_ = smallInsertionRate;
  json_file_["param_in"]["mutations"]["small_insertion_rate"] = small_insertion_rate_;
}
double IOJson::getSmallDeletionRate() const { return small_deletion_rate_; }
void IOJson::setSmallDeletionRate(double smallDeletionRate) {
  small_deletion_rate_ = smallDeletionRate;
  json_file_["param_in"]["mutations"]["small_deletion_rate"] = small_deletion_rate_;
}

int IOJson::getMaxIndelSize() const { return max_indel_size_; }
void IOJson::setMaxIndelSize(int maxIndelSize) {
  max_indel_size_ = maxIndelSize;
  json_file_["param_in"]["mutations"]["max_indel_size"] = max_indel_size_;
}

double IOJson::getDuplicationRate() const { return duplication_rate_; }
void IOJson::setDuplicationRate(double duplicationRate) {
  duplication_rate_ = duplicationRate;
  json_file_["param_in"]["mutations"]["duplication_rate"] = duplication_rate_;
}

double IOJson::getDeletionRate() const { return deletion_rate_; }
void IOJson::setDeletionRate(double deletionRate) {
  deletion_rate_ = deletionRate;
  json_file_["param_in"]["mutations"]["deletion_rate"] = deletion_rate_;
}

double IOJson::getTranslocationRate() const { return translocation_rate_; }
void IOJson::setTranslocationRate(double translocationRate) {
  translocation_rate_ = translocationRate;
  json_file_["param_in"]["mutations"]["translocation_rate"] = translocation_rate_;
}

double IOJson::getInversionRate() const { return inversion_rate_; }
void IOJson::setInversionRate(double inversionRate) {
  inversion_rate_ = inversionRate;
  json_file_["param_in"]["mutations"]["inversion_rate"] = inversion_rate_;
}

bool IOJson::isWithAlignments() const { return with_alignments_; }
void IOJson::setWithAlignments(bool withAlignments) {
  with_alignments_ = withAlignments;
  json_file_["param_in"]["mutations"]["alignments"]["with_alignments"] = with_alignments_;
}

int IOJson::getChromosomeMaxLength() const { return max_genome_length_; }
void IOJson::setChromosomeMaxLength(int chromosomeMaxLength) {
  max_genome_length_ = chromosomeMaxLength;
  json_file_["param_in"]["plasmids"]["chromosome_max_length"] = chromosomeMaxLength;
}

bool IOJson::isAllowPlasmids() const { return allow_plasmids_; }
void IOJson::setAllowPlasmids(bool allowPlasmids) {
  allow_plasmids_ = allowPlasmids;
  json_file_["param_in"]["plasmids"]["allow_plasmids"] = allowPlasmids;
}

AlignmentFunctionShape IOJson::getAlignFunShape() const {
  return align_fun_shape_;
}
void IOJson::setAlignFunShape(AlignmentFunctionShape alignFunShape) {
  align_fun_shape_ = alignFunShape;
  json_file_["param_in"]["mutations"]["alignments"]["align_fun_shape"] = alignFunShape;
}

int IOJson::getAlignLinMax() const { return align_lin_max_; }
void IOJson::setAlignLinMax(int alignLinMax) {
  align_lin_max_ = alignLinMax;
  json_file_["param_in"]["mutations"]["alignments"]["align_lin_max_"] = alignLinMax;
}

int IOJson::getAlignLinMin() const { return align_lin_min_; }
void IOJson::setAlignLinMin(int alignLinMin) {
  align_lin_min_ = alignLinMin;
  json_file_["param_in"]["mutations"]["alignments"]["align_lin_min"] = alignLinMin;
}

int IOJson::getAlignMatchBonus() const { return align_match_bonus_; }
void IOJson::setAlignMatchBonus(int alignMatchBonus) {
  align_match_bonus_ = alignMatchBonus;
  json_file_["param_in"]["mutations"]["alignments"]["align_match_bonus"] = alignMatchBonus;
}

int IOJson::getAlignMaxShift() const { return align_max_shift_; }
void IOJson::setAlignMaxShift(int alignMaxShift) {
  align_max_shift_ = alignMaxShift;
  json_file_["param_in"]["mutations"]["alignments"]["align_max_shift"] = alignMaxShift;
}
int IOJson::getAlignMismatchCost() const { return align_mismatch_cost_; }
void IOJson::setAlignMismatchCost(int alignMismatchCost) {
  align_mismatch_cost_ = alignMismatchCost;
  json_file_["param_in"]["mutations"]["alignments"]["align_mismatch_cost"] = alignMismatchCost;
}
int IOJson::getAlignSigmLambda() const { return align_sigm_lambda_; }
void IOJson::setAlignSigmLambda(int alignSigmLambda) {
  align_sigm_lambda_ = alignSigmLambda;
  json_file_["param_in"]["mutations"]["alignments"]["align_sigm_lambda"] = alignSigmLambda;
}
int IOJson::getAlignSigmMean() const { return align_sigm_mean_; }
void IOJson::setAlignSigmMean(int alignSigmMean) {
  align_sigm_mean_ = alignSigmMean;
  json_file_["param_in"]["mutations"]["alignments"]["align_sigm_mean"] = alignSigmMean;
}
int IOJson::getAlignWZoneHLen() const { return align_w_zone_h_len_; }
void IOJson::setAlignWZoneHLen(int alignWZoneHLen) {
  align_w_zone_h_len_ = alignWZoneHLen;
  json_file_["param_in"]["mutations"]["alignments"]["align_w_zone_h_len"] = alignWZoneHLen;
}
int IOJson::getBigBackupStep() const { return big_backup_step_; }
void IOJson::setBigBackupStep(int bigBackupStep) {
  big_backup_step_ = bigBackupStep;
  json_file_["param_in"]["output"]["big_backup_step"] = bigBackupStep;
}
int IOJson::getChromosomeMinimalLength() const {
  return chromosome_minimal_length_;
}
void IOJson::setChromosomeMinimalLength(int chromosomeMinimalLength) {
  chromosome_minimal_length_ = chromosomeMinimalLength;
  json_file_["param_in"]["plasmids"]["chromosome_minimal_length"] = chromosomeMinimalLength;
}
int IOJson::getChromosomeMaximalLength() const {
  return chromosome_maximal_length_;
}
void IOJson::setChromosomeMaximalLength(int chromosomeMaximalLength) {
  chromosome_maximal_length_ = chromosomeMaximalLength;
  json_file_["param_in"]["plasmids"]["chromosome_maximal_length"] = chromosomeMaximalLength;
}
bool IOJson::isComputePhenContribByGu() const {
  return compute_phen_contrib_by_GU_;
}
void IOJson::setComputePhenContribByGu(bool computePhenContribByGu) {
  compute_phen_contrib_by_GU_ = computePhenContribByGu;
  json_file_["param_in"]["plasmids"]["compute_phen_contrib_by_GU"] = computePhenContribByGu;
}
bool IOJson::isDeleteOldStats() const { return delete_old_stats_; }
void IOJson::setDeleteOldStats(bool deleteOldStats) {
  delete_old_stats_ = deleteOldStats;
  json_file_["param_in"]["output"]["delete_old_stats"] = deleteOldStats;
}
double IOJson::getDeletionProportion() const { return deletion_proportion_; }
void IOJson::setDeletionProportion(double deletionProportion) {
  deletion_proportion_ = deletionProportion;
  json_file_["param_in"]["mutations"]["alignments"]["deletion_proportion"] = deletionProportion;
}
double IOJson::getDonorCost() const { return donor_cost_; }
void IOJson::setDonorCost(double donorCost) {
  donor_cost_ = donorCost;
  json_file_["param_in"]["mutations"]["alignments"]["donor_cost"] = donorCost;
}
int IOJson::getDumpStep() const { return dump_step_; }
void IOJson::setDumpStep(int dumpStep) {
  dump_step_ = dumpStep;
  json_file_["param_in"]["output"]["dump_step"] = dumpStep;
}
double IOJson::getDuplicationProportion() const {
  return duplication_proportion_;
}
void IOJson::setDuplicationProportion(double duplicationProportion) {
  duplication_proportion_ = duplicationProportion;
  json_file_["param_in"]["mutations"]["alignments"]["duplication_proportion"] = duplicationProportion;
}
int IOJson::getEnvAxisNbSegments() const { return env_axis_nb_segments_; }
void IOJson::setEnvAxisNbSegments(int envAxisNbSegments) {
  env_axis_nb_segments_ = envAxisNbSegments;
  json_file_["param_in"]["env"]["segmentation"]["env_axis_nb_segments"] = envAxisNbSegments;
}
PhenotypicFeature *IOJson::getEnvAxisFeatures() const {
  return env_axis_features_;
}
void IOJson::setEnvAxisFeatures(PhenotypicFeature *envAxisFeatures) {
  env_axis_features_ = new PhenotypicFeature[env_axis_nb_segments_];
  json::array_t json_array;
  if (envAxisFeatures != NULL) {
    for (int i = 0; i < env_axis_nb_segments_; ++i) {
      env_axis_features_[i] = envAxisFeatures[i];
      switch (envAxisFeatures[i]) {
      case PhenotypicFeature::DONOR:
        json_array.emplace_back("DONOR");
        break;
      case PhenotypicFeature::METABOLISM:
        json_array.emplace_back("METABOLISM");
        break;
      case PhenotypicFeature::NEUTRAL:
        json_array.emplace_back("NEUTRAL");
        break;
      case PhenotypicFeature::SECRETION:
        json_array.emplace_back("SECRETION");
        break;
      case PhenotypicFeature::RECIPIENT:
        json_array.emplace_back("RECIPIENT");
        break;
      }
    }
  }
  json_file_["param_in"]["env"]["segmentation"]["env_axis_features"] = json_array;
}
double *IOJson::getEnvAxisSegmentBoundaries() const {
  return env_axis_segment_boundaries_;
}
void IOJson::setEnvAxisSegmentBoundaries(double *envAxisSegmentBoundaries) {
  env_axis_segment_boundaries_ = envAxisSegmentBoundaries;
  /*if (envAxisSegmentBoundaries != NULL) {
    json::array_t json_array;
    for (int i = 0; i < env_axis_nb_segments_+1; ++i) {
      json_array.emplace_back(envAxisSegmentBoundaries[i]);
    }
    json_file_["param_in"]["env"]["segmentation"]["env_axis_boundaries"] = json_array;
  }*/
}
bool IOJson::isEnvAxisSeparateSegments() const {
  return env_axis_separate_segments_;
}
void IOJson::setEnvAxisSeparateSegments(bool envAxisSeparateSegments) {
  env_axis_separate_segments_ = envAxisSeparateSegments;
  json_file_["param_in"]["env"]["segmentation"]["env_axis_separate_segments"] = envAxisSeparateSegments;
}
double IOJson::getEnvNoiseAlpha() const { return env_noise_alpha_; }
void IOJson::setEnvNoiseAlpha(double envNoiseAlpha) {
  env_noise_alpha_ = envNoiseAlpha;
  json_file_["param_in"]["env"]["noise"]["env_noise_alpha"] = envNoiseAlpha;
}
int IOJson::getEnvNoiseSamplingLog() const { return env_noise_sampling_log_; }
void IOJson::setEnvNoiseSamplingLog(int envNoiseSamplingLog) {
  env_noise_sampling_log_ = envNoiseSamplingLog;
  json_file_["param_in"]["env"]["noise"]["env_noise_sampling_log"] = envNoiseSamplingLog;
}
PhenotypicTargetNoiseMethod IOJson::getEnvNoiseMethod() const {
  return env_noise_method_;
}
void IOJson::setEnvNoiseMethod(PhenotypicTargetNoiseMethod envNoiseMethod) {
  env_noise_method_ = envNoiseMethod;
  switch (envNoiseMethod) {
  case PhenotypicTargetNoiseMethod::NO_NOISE:
    json_file_["param_in"]["env"]["noise"]["env_noise_method"] = "NO_NOISE";
    break;
  case PhenotypicTargetNoiseMethod::FRACTAL:
    json_file_["param_in"]["env"]["noise"]["env_noise_method"] = "FRACTAL";
    break;
  }
}
double IOJson::getEnvNoiseProb() const { return env_noise_prob_; }
void IOJson::setEnvNoiseProb(double envNoiseProb) {
  env_noise_prob_ = envNoiseProb;
  json_file_["param_in"]["env"]["noise"]["env_noise_prob"] = envNoiseProb;
}
int IOJson::getEnvNoiseSeed() const { return env_noise_seed_; }
void IOJson::setEnvNoiseSeed(int envNoiseSeed) {
  env_noise_seed_ = envNoiseSeed;
  json_file_["param_in"]["env"]["noise"]["env_noise_seed"] = envNoiseSeed;
}
double IOJson::getEnvNoiseSigma() const { return env_noise_sigma_; }
void IOJson::setEnvNoiseSigma(double envNoiseSigma) {
  env_noise_sigma_ = envNoiseSigma;
  json_file_["param_in"]["env"]["noise"]["env_noise_sigma"] = envNoiseSigma;
}
PhenotypicTargetVariationMethod IOJson::getEnvVarMethod() const {
  return env_var_method_;
}
void IOJson::setEnvVarMethod(PhenotypicTargetVariationMethod envVarMethod) {
  env_var_method_ = envVarMethod;
  switch (envVarMethod) {
  case PhenotypicTargetVariationMethod::SWITCH_IN_A_LIST :
    json_file_["param_in"]["env"]["variation"]["env_var_method"] = "SWITCH_IN_A_LIST";
    break;
  case PhenotypicTargetVariationMethod::AUTOREGRESSIVE_HEIGHT_VAR :
    json_file_["param_in"]["env"]["variation"]["env_var_method"] = "AUTOREGRESSIVE_HEIGHT_VAR";
    break;
  case PhenotypicTargetVariationMethod::AUTOREGRESSIVE_MEAN_VAR :
    json_file_["param_in"]["env"]["variation"]["env_var_method"] = "AUTOREGRESSIVE_MEAN_VAR";
    break;
  case PhenotypicTargetVariationMethod::ONE_AFTER_ANOTHER :
    json_file_["param_in"]["env"]["variation"]["env_var_method"] = "ONE_AFTER_ANOTHER";
    break;
  case PhenotypicTargetVariationMethod::NO_VAR :
    json_file_["param_in"]["env"]["variation"]["env_var_method"] = "NO_VAR";
    break;
  case PhenotypicTargetVariationMethod::LOCAL_GAUSSIANS_VAR :
    json_file_["param_in"]["env"]["variation"]["env_var_method"] = "LOCAL_GAUSSIANS_VAR";
    break;
  }
}
double IOJson::getEnvVarSigma() const { return env_var_sigma_; }
void IOJson::setEnvVarSigma(double envVarSigma) {
  env_var_sigma_ = envVarSigma;
  json_file_["param_in"]["env"]["variation"]["env_var_sigma"] = envVarSigma;
}
int IOJson::getEnvVarSeed() const { return env_var_seed_; }
void IOJson::setEnvVarSeed(int envVarSeed) {
  env_var_seed_ = envVarSeed;
  json_file_["param_in"]["env"]["variation"]["env_var_seed"] = envVarSeed;
}
int IOJson::getEnvVarTau() const { return env_var_tau_; }
void IOJson::setEnvVarTau(int envVarTau) {
  env_var_tau_ = envVarTau;
  json_file_["param_in"]["env"]["variation"]["env_var_tau"] = envVarTau;
}
const std::string &IOJson::getFilename() const { return filename_; }
void IOJson::setFilename(const std::string &filename) {
  filename_ = filename;
}
FitnessFunction IOJson::getFitnessFunction() const { return fitness_function_; }
void IOJson::setFitnessFunction(FitnessFunction fitnessFunction) {
  fitness_function_ = fitnessFunction;
  switch (fitnessFunction) {
  case FitnessFunction::FITNESS_EXP:
    json_file_["param_in"]["selection"]["fitness"]["fitness_function"] = "FITNESS_EXP";
    break;
  case FitnessFunction::FITNESS_GLOBAL_SUM:
    json_file_["param_in"]["selection"]["fitness"]["fitness_function"] = "FITNESS_GLOBAL_SUM";
    break;
  case FitnessFunction::FITNESS_LOCAL_SUM:
    json_file_["param_in"]["selection"]["fitness"]["fitness_function"] = "FITNESS_LOCAL_SUM";
    break;
  }
}
int IOJson::getFitnessFunctionX() const { return fitness_function_x_; }
void IOJson::setFitnessFunctionX(int fitnessFunctionX) {
  fitness_function_x_ = fitnessFunctionX;
  json_file_["param_in"]["selection"]["fitness"]["fitness_function_x"] = fitnessFunctionX;
}
int IOJson::getFitnessFunctionY() const { return fitness_function_y_; }
void IOJson::setFitnessFunctionY(int fitnessFunctionY) {
  fitness_function_y_ = fitnessFunctionY;
  json_file_["param_in"]["selection"]["fitness"]["fitness_function_y"] = fitnessFunctionY;
}
double IOJson::getHtInsRate() const { return HT_ins_rate_; }
void IOJson::setHtInsRate(double htInsRate) {
  HT_ins_rate_ = htInsRate;
  json_file_["param_in"]["mutations"]["HT"]["HT_ins_rate"] = htInsRate;
}
double IOJson::getHtReplRate() const { return HT_repl_rate_; }
void IOJson::setHtReplRate(double htReplRate) {
  HT_repl_rate_ = htReplRate;
  json_file_["param_in"]["mutations"]["HT"]["HT_repl_rate"] = htReplRate;
}
double IOJson::getInversionProportion() const { return inversion_proportion_; }
void IOJson::setInversionProportion(double inversionProportion) {
  inversion_proportion_ = inversionProportion;
  json_file_["param_in"]["mutations"]["alignments"]["inversion_proportion"] = inversionProportion;
}
int IOJson::getLogs() const { return logs_; }
void IOJson::setLogs(int logs) {
  logs_ = logs;
  json_file_["param_in"]["output"]["logs"] = logs;
}
bool IOJson::isMakeDumps() const { return make_dumps_; }
void IOJson::setMakeDumps(bool makeDumps) {
  make_dumps_ = makeDumps;
  json_file_["param_in"]["output"]["make_dumps"] = makeDumps;
}
int IOJson::getMaxGenomeLength() const { return max_genome_length_; }
void IOJson::setMaxGenomeLength(int maxGenomeLength) {
  max_genome_length_ = maxGenomeLength;
  json_file_["param_in"]["constraints"]["max_genome_length"] = maxGenomeLength;
}
int IOJson::getMinGenomeLength() const { return min_genome_length_; }
void IOJson::setMinGenomeLength(int minGenomeLength) {
  min_genome_length_ = minGenomeLength;
  json_file_["param_in"]["constraints"]["min_genome_length"] = minGenomeLength;
}
int IOJson::getMutSeed() const { return mut_seed_; }
void IOJson::setMutSeed(int mutSeed) {
  mut_seed_ = mutSeed;
  json_file_["param_in"]["rng"]["mut_seed"] = mutSeed;
}
double IOJson::getNeighbourhoodRate() const { return neighbourhood_rate_; }
void IOJson::setNeighbourhoodRate(double neighbourhoodRate) {
  neighbourhood_rate_ = neighbourhoodRate;
  json_file_["param_in"]["mutations"]["alignments"]["neighbourhood_rate"] = neighbourhoodRate;
}
int32_t IOJson::getPartialMixNbPermutations() const {
  return partial_mix_nb_permutations_;
}
void IOJson::setPartialMixNbPermutations(int32_t partialMixNbPermutations) {
  partial_mix_nb_permutations_ = partialMixNbPermutations;
  json_file_["param_in"]["partial_mix_nb_permutations"] = partialMixNbPermutations;
}
int IOJson::getPlasmidInitialGene() const { return plasmid_initial_gene_; }
void IOJson::setPlasmidInitialGene(int plasmidInitialGene) {
  plasmid_initial_gene_ = plasmidInitialGene;
  json_file_["param_in"]["plasmids"]["plasmid_initial_gene"] = plasmidInitialGene;
}
int IOJson::getPlasmidInitialLength() const { return plasmid_initial_length_; }
void IOJson::setPlasmidInitialLength(int plasmidInitialLength) {
  plasmid_initial_length_ = plasmidInitialLength;
  json_file_["param_in"]["plasmids"]["plasmid_initial_length_"] = plasmidInitialLength;
}
int IOJson::getPlasmidMaximalLength() const { return plasmid_maximal_length_; }
void IOJson::setPlasmidMaximalLength(int plasmidMaximalLength) {
  plasmid_maximal_length_ = plasmidMaximalLength;
  json_file_["param_in"]["plasmids"]["plasmid_maximal_length"] = plasmidMaximalLength;
}
int IOJson::getPlasmidMinimalLength() const { return plasmid_minimal_length_; }
void IOJson::setPlasmidMinimalLength(int plasmidMinimalLength) {
  plasmid_minimal_length_ = plasmidMinimalLength;
  json_file_["param_in"]["plasmids"]["plasmid_minimal_length"] = plasmidMinimalLength;
}
double IOJson::getProbPlasmidHt() const { return prob_plasmid_HT_; }
void IOJson::setProbPlasmidHt(double probPlasmidHt) {
  prob_plasmid_HT_ = probPlasmidHt;
  json_file_["param_in"]["plasmids"]["prob_plasmid_HT"] = probPlasmidHt;
}
const std::shared_ptr<JumpingMT> &IOJson::getPrng() const { return prng_; }
void IOJson::setPrng(const std::shared_ptr<JumpingMT> &prng) {
  prng_ = prng;
}
double IOJson::getRecipientCost() const { return recipient_cost_; }
void IOJson::setRecipientCost(double recipientCost) {
  recipient_cost_ = recipientCost;
  json_file_["param_in"]["plasmids"]["recipient_cost"] = recipientCost;
}
bool IOJson::isRecordLightTree() const { return record_light_tree_; }
void IOJson::setRecordLightTree(bool recordLightTree) {
  record_light_tree_ = recordLightTree;
  json_file_["param_in"]["output"]["record_light_tree"] = recordLightTree;
}
double IOJson::getReplHtDetachRate() const { return repl_HT_detach_rate_; }
void IOJson::setReplHtDetachRate(double replHtDetachRate) {
  repl_HT_detach_rate_ = replHtDetachRate;
  json_file_["param_in"]["mutations"]["HT"]["repl_HT_detach_rate"] = replHtDetachRate;
}
bool IOJson::isReplHtWithClosePoints() const {
  return repl_HT_with_close_points_;
}
void IOJson::setReplHtWithClosePoints(bool replHtWithClosePoints) {
  repl_HT_with_close_points_ = replHtWithClosePoints;
  json_file_["param_in"]["mutations"]["HT"]["repl_HT_with_close_points"] = replHtWithClosePoints;
}
double IOJson::getSecretionContribToFitness() const {
  return secretion_contrib_to_fitness_;
}
void IOJson::setSecretionContribToFitness(double secretionContribToFitness) {
  secretion_contrib_to_fitness_ = secretionContribToFitness;
  json_file_["param_in"]["secretion"]["secretion_contrib_to_fitness"] = secretionContribToFitness;
}
double IOJson::getSecretionCost() const { return secretion_cost_; }
void IOJson::setSecretionCost(double secretionCost) {
  secretion_cost_ = secretionCost;
  json_file_["param_in"]["secretion"]["secretion_cost"] = secretionCost;
}
double IOJson::getSecretionDegradationProp() const {
  return secretion_degradation_prop_;
}
void IOJson::setSecretionDegradationProp(double secretionDegradationProp) {
  secretion_degradation_prop_ = secretionDegradationProp;
  json_file_["param_in"]["secretion"]["secretion_degradation_prop"] = secretionDegradationProp;
}
double IOJson::getSecretionDiffusionProp() const {
  return secretion_diffusion_prop_;
}
void IOJson::setSecretionDiffusionProp(double secretionDiffusionProp) {
  secretion_diffusion_prop_ = secretionDiffusionProp;
  json_file_["param_in"]["secretion"]["secretion_diffusion_prop"] = secretionDiffusionProp;
}
double IOJson::getSecretionInit() const { return secretion_init_; }
void IOJson::setSecretionInit(double secretionInit) {
  secretion_init_ = secretionInit;
  json_file_["param_in"]["secretion"]["secretion_init"] = secretionInit;
}
double IOJson::getSelectionPressure() const { return selection_pressure_; }
void IOJson::setSelectionPressure(double selectionPressure) {
  selection_pressure_ = selectionPressure;
  json_file_["param_in"]["selection"]["selection_pressure"] = selectionPressure;
}

SelectionScope IOJson::getSelectionScope() const { return selection_scope_; }
void IOJson::setSelectionScope(SelectionScope selectionScope) {
  selection_scope_ = selectionScope;
  json_file_["param_in"]["selection"]["selection_scope"] = selectionScope;
}
int32_t IOJson::getSelectionScopeX() const { return selection_scope_x_; }
void IOJson::setSelectionScopeX(int32_t selectionScopeX) {
  selection_scope_x_ = selectionScopeX;
  json_file_["param_in"]["selection"]["selection_scope_x"] = selectionScopeX;
}
int32_t IOJson::getSelectionScopeY() const { return selection_scope_y_; }
void IOJson::setSelectionScopeY(int32_t selectionScopeY) {
  selection_scope_y_ = selectionScopeY;
  json_file_["param_in"]["selection"]["selection_scope_y"] = selectionScopeY;
}
int IOJson::getStats() const { return stats_; }
void IOJson::setStats(int stats) {
  stats_ = stats;
  json_file_["param_in"]["output"]["stats"] = stats;
}
int IOJson::getStochSeed() const { return stoch_seed_; }
void IOJson::setStochSeed(int stochSeed) {
  stoch_seed_ = stochSeed;
  json_file_["param_in"]["rng"]["stoch_seed"] = stochSeed;
}
bool IOJson::isSwapGUs() const { return swap_GUs_; }
void IOJson::setSwapGUs(bool swapGUs) {
  swap_GUs_ = swapGUs;
  json_file_["param_in"]["plasmids"]["swap_GUs_"] = swapGUs;
}
double IOJson::getTranslationCost() const { return translation_cost_; }
void IOJson::setTranslationCost(double translationCost) {
  translation_cost_ = translationCost;
  json_file_["param_in"]["plasmids"]["translation_cost"] = translationCost;
}
double IOJson::getTranslocationProportion() const {
  return translocation_proportion_;
}
void IOJson::setTranslocationProportion(double translocationProportion) {
  translocation_proportion_ = translocationProportion;
  json_file_["param_in"]["mutations"]["alignments"]["translocation_proportion"] = translocationProportion;
}
int IOJson::getTreeStep() const { return tree_step_; }
void IOJson::setTreeStep(int treeStep) {
  tree_step_ = treeStep;
  json_file_["param_in"]["output"]["tree_step"] = treeStep;
}
double IOJson::getTuneDonorAbility() const { return tune_donor_ability_; }
void IOJson::setTuneDonorAbility(double tuneDonorAbility) {
  tune_donor_ability_ = tuneDonorAbility;
  json_file_["param_in"]["plasmids"]["tune_donor_ability"] = tuneDonorAbility;
}
double IOJson::getTuneRecipientAbility() const { return tune_recipient_ability_; }
void IOJson::setTuneRecipientAbility(double tuneRecipientAbility) {
  tune_recipient_ability_ = tuneRecipientAbility;
  json_file_["param_in"]["plasmids"]["tune_recipient_ability"] = tuneRecipientAbility;
}
bool IOJson::isWithHt() const { return with_HT_; }
void IOJson::setWithHt(bool withHt) {
  with_HT_ = withHt;
  json_file_["param_in"]["mutations"]["HT"]["with_HT"] = withHt;
}
bool IOJson::isWithStochasticity() const { return with_stochasticity_; }
void IOJson::setWithStochasticity(bool withStochasticity) {
  with_stochasticity_ = withStochasticity;
  json_file_["param_in"]["with_stochasticity"] = withStochasticity;
}
bool IOJson::isWithSecretion() const { return with_secretion_; }
void IOJson::setWithSecretion(bool withSecretion) {
  with_secretion_ = withSecretion;
  json_file_["param_in"]["secretion"]["with_secretion"] = withSecretion;
}
bool IOJson::isWith4PtsTrans() const { return with_4pts_trans_; }
void IOJson::setWith4PtsTrans(bool with4PtsTrans) {
  with_4pts_trans_ = with4PtsTrans;
  json_file_["param_in"]["mutations"]["with_4pts_trans"] = with4PtsTrans;
}
bool IOJson::isWellMixed() const { return well_mixed; }
void IOJson::setWellMixed(bool wellMixed) {
  well_mixed = wellMixed;
  json_file_["param_in"]["well_mixed"] = wellMixed;
}
void IOJson::setStrainName(const std::string &strainName) {
  strain_name_ = strainName;
  json_file_["param_in"]["strain_name"] = strainName;
}
const std::string &IOJson::getStrainName() const {
  return strain_name_;
}
void IOJson::write(const std::string &filename) const {
  std::ofstream file(filename);
  file << json_file_.dump(2) << std::endl;
}
void IOJson::addIndividual(Individual *indiv, json gu_list) {
    individuals_.push_back(indiv);
    json my_indiv;
    my_indiv["GU"] = gu_list;
    my_indiv["id"] = indiv->id();
    my_indiv["generation"] = indiv->age();
    if(indiv->grid_cell() != nullptr) {
      my_indiv["x_pos"] = indiv->grid_cell_->x();
      my_indiv["y_pos"] = indiv->grid_cell_->y();
    }
    json_file_["indivs"].push_back(my_indiv);

}
const vector<Individual*> IOJson::getIndividuals() const{
  return individuals_;
}

void IOJson::setIndividuals(const vector<Individual*> &individuals){
  individuals_ = individuals;
}
std::string IOJson::getIndividualSequence(int32_t index, int32_t gu) const{
  return json_file_["indivs"][index]["GU"][gu]["seq"];
}
void IOJson::setIndividualSequence(int32_t index, int32_t gu, const char* seq){
  json_file_["indivs"][index]["GU"][gu]["seq"] = seq;
}
int32_t IOJson::getNbrIndividuals() const {
  return json_file_["indivs"].size();
}
